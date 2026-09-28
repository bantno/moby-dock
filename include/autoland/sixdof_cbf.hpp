#pragma once
#include <Eigen/Dense>
#include <array>
#include <memory>
#include <string>
#include <vector>
#include "autoland/beaver_dynamics.hpp"
#include "autoland/config.hpp"
#include "autoland/lie_taylor.hpp"
#include "autoland/qp_solver.hpp"
#include "autoland/world_geometry.hpp"

// =============================================================================
// 6-DOF CBF-QP safety filter on the Beaver plant: SURFACES ONLY.
//
// Decision variables are the three surfaces U = [de, da, dr]; the throttle
// stays with the nominal and enters the barrier model as a KNOWN constant
// (it is computed before the QP), so no engine-lag state is needed and every
// row below has uniform relative degree 2 in U.
//
// Control-affine split at the current step (x0, u_prev, dT):
//   xdot(x, u) = f(x) + g(x) U + F(x) (U - u_prev)
//   f(x)  = plant(x, surfaces = u_prev, dT) - g(x) u_prev      (drift)
//   g(x)  = MOMENT columns of the surfaces: de -> qdot (Cm_de); da, dr ->
//           pdot, rdot through the Ixz inertia solve (Cl_da + Cl_daa alpha,
//           Cn_da, Cl_dr, Cn_dr). All affine in U -- so f + g U reproduces the
//           plant EXACTLY with the moments at U and the FORCES at u_prev.
//   F(x)  = the direct control-FORCE columns (Cz_de + Cz_deb2 beta^2, Cy_da,
//           Cy_dr + Cy_dra alpha, Cx_dr) that are left out of the control
//           map so the rows keep a clean relative degree.
// The model mismatch is therefore F (U - u_prev), bounded per step by the
// surface rate limit: |dU_k| <= rate dt. It enters the degree-2 condition
//   hddot + (c1 + c2) hdot + c1 c2 h >= 0
// at BOTH orders: hdot_true = L_f h + (L_F h) dU and hddot_true = L_f^2 h +
// L_g L_f h U + (L_F L_f h + L_f L_F h) dU (ZOH: dU is constant within the
// step). Each row therefore carries the ROBUSTNESS MARGIN
//   sum_k ( |D2_k| + (c1 + c2) |D1_k| ) * min(rate dt, range_k)
// with D1_k = L_{F_k} h and D2_k = L_{F_k} L_f h + L_f L_{F_k} h, both exact:
// they are the first-order coefficients in epsilon of the Lie stack of the
// perturbed drift f + epsilon F_k, obtained by seeding epsilon as an autodiff
// dual through the Taylor jet. An optional constant per row covers the
// sampled-data gap. The guarantee is for the modelled system with that
// margin; it is stated as such everywhere.
//
// Lie derivatives are exact Taylor jets of the actual plant polynomials
// (BeaverDynamics::xdotCoreT is scalar-generic); nothing is transcribed by
// hand. Rows are stacked as in lon_cbf_filter: hard rows without slack, soft
// rows with quadratic slack, best-effort minimum-violation fallback when the
// hard set is infeasible (never a silent relaxation), row normalization by
// max|a| (energy_reach_cbf lesson), input box = deflection AND rate limits
// about u_prev.
//
// WIND: the aerodynamics see the AIR-relative velocity. The model takes an
// earth-frame wind estimate W_e (north tailwind+, east+, updraft+ -- the
// plant's own convention) and rotates it into the body frame inside the jet,
// so alpha, beta, V and every force are air-relative exactly as in the
// plant, with the wind assumed constant over the step. In the sim W_e is the
// true plant wind: the idealization is a perfect wind estimate (ground
// velocity minus air velocity, which any GPS + air-data autopilot forms);
// without it a crosswind crab reads as sideslip and the rows fight the crab.
//
// Phase 1 rows (envelope): AoA ceiling alpha_max - alpha; airspeed floor and
// ceiling V - V_min, V_max - V; bank +/-phi_max -/+ phi; sideslip
// +/-beta_max -/+ beta -- all relative degree 2.
// Phase 2 row (terrain keep-out), BRAKING form (default): with the clearance
//   z = h - h_T(N, E) - margin   and the closure rate   zdot = hdot - grad h_T . v,
//   b = zdot + sqrt(v_safe^2 + 2 a_brk z_+) + k_neg z_- >= 0,
// where h_T is the CLIMB-CAPABILITY LOOKAHEAD terrain: the smooth maximum,
// over the aircraft's position and four points up to terrain_lookahead ahead
// along each of a FAN of directions about the ground course, of
// h_T(s) - s * terrain_climb_grad -- the height needed NOW to clear the
// terrain at s on a sustainable climb. The fan centre (the course) is FROZEN
// at the current step so the row depends on position only and keeps
// relative degree 2; the fan's spread keeps the barrier continuous under
// small course changes (a single frozen ray swept across a peak made the
// barrier jump by hundreds of metres between steps and the hard row flung
// the aircraft into a stall). (A local-gradient barrier is
// myopic: it let the aircraft settle at the margin on the flat and then meet
// a slope no climb can fly -- closure from a slope is not "arrested" by a
// pull-up. With the capability lookahead the row acts from far enough back,
// and turning changes the closure against the terrain ahead, which is what
// lets it steer around or through terrain),
//   z_+ = smooth max(z, 0),  z_- = z - z_+  (smooth min(z, 0)),
// i.e. the aircraft may close on the terrain no faster than a recovery at
// a_brk can arrest within the remaining clearance (the lon filter's descent
// envelope, extended with the terrain slope), and below the margin it must
// CLIMB at k_neg |z_-| -- without that term the lon form, built to ARRIVE at
// z = 0 at v_safe, keeps sinking through the margin at v_safe. Degree 2. A pure
// clearance row z >= 0 with a linear class-K chain (degree 3, kept as the
// `terrain_degree3` option) is far too permissive far from the ground -- at
// 57 m it tolerates a 57 m/s sink -- and its pull-up comes too late for a
// 1-2 g airplane; the braking form encodes the flyable safe set directly.
// The terrain field's exact Taylor jets (world_geometry.hpp) serve both.
// The margin for degree R generalizes to sum_k sum_{j=1..R} coeff_j |D_j,k| dU_k
// with coeff the class-K elementary-symmetric coefficients and D_j,k the
// epsilon-coefficient of the order-j Lie derivative under f + eps F_k.
// See documentation/corridor_landing_roadmap.md.
// =============================================================================
namespace autoland {

constexpr int NUS = 3;  // surfaces: DE, DA, DR
using SurfVec = Eigen::Vector3d;

struct SixDofCbfConfig {
  bool enabled{false};
  bool filter{true};          // false = monitor only (A/B twin)
  // Row limits. NOTE they must sit a margin ABOVE the nominal's own clamps
  // (guidance phi_max, AoA band) so the rows are rare interventions by
  // construction rather than always-binding duplicates of the clamps.
  bool row_alpha{true}, row_V{true}, row_phi{true}, row_beta{true};
  // Terrain keep-out row (needs `terrain`): h - h_T(N,E) - terrain_margin >= 0,
  // evaluated only while the clearance is below terrain_gate (cost).
  bool row_terrain{false};
  const TerrainField* terrain{nullptr};
  double terrain_margin{5.0};     // [m], see terrainRowMargin()
  double terrain_gate{200.0};     // [m] clearance below which the row is built
  bool terrain_degree3{false};    // pure clearance row (degree 3) instead of braking
  double terrain_a_brk{3.0};      // recovery acceleration the braking form assumes [m/s^2]
  double terrain_v_safe{0.1};     // residual closure rate allowed at zero clearance [m/s]
  double terrain_z_eps{1.0};      // smooth-max width for z_+ [m]
  double terrain_k_neg{1.0};      // recovery rate below the margin [1/s]
  double terrain_lookahead{0.0};  // [m] along the ground course (0 = local only)
  double terrain_smax_eps{2.0};   // smooth-max width for the lookahead heights [m]
  double terrain_climb_grad{0.05}; // sustainable climb gradient assumed ahead [-]
  double terrain_fan_half{20.0 * M_PI / 180.0};  // lookahead fan half-angle [rad]
  int terrain_fan_n{5};            // fan directions (odd), 1 = single ray
  std::array<double, 3> c_terrain{1.0, 1.0, 1.0};  // degree-2 form uses [0], [1]
  bool hard_terrain{true};
  double w_terrain{1e3};
  double alpha_max{12.0 * M_PI / 180.0};   // LR-556 trusted band (stall ~16)
  double V_min{33.0}, V_max{55.0};         // LR-556 30-55 m/s band, floor + 3
  double phi_max{30.0 * M_PI / 180.0};
  double beta_max{12.0 * M_PI / 180.0};
  // Linear class-K gains (degree 2 -> two per row) [1/s].
  std::array<double, 2> c_alpha{2.0, 2.0}, c_V{1.0, 1.0}, c_phi{3.0, 3.0},
      c_beta{3.0, 3.0};
  // Policy: hard = no slack; soft = quadratic slack at the weight below.
  bool hard_alpha{false}, hard_V{false}, hard_phi{false}, hard_beta{false};
  double w_alpha{1e3}, w_V{1e2}, w_phi{1e2}, w_beta{1e2};
  // Best-effort: hard-row slack penalty. Rows are normalized to |a| <= 1 and
  // the tracking weights are O(1), so 1e2 already dominates; 1e4 and above
  // make OSQP (eps 1e-8) hit its iteration limit on the violated-row
  // problem and the filter would silently fall back to the nominal.
  double best_effort_penalty{1e2};
  // Tracking weights on (U - U_nom), each divided by the surface range^2.
  double w_de{1.0}, w_da{1.0}, w_dr{1.0};
  // Robustness margins: force-term sensitivity (computed) and an extra
  // constant per row [row units/s^2] for the sampled-data gap.
  bool margin_force_terms{true};
  double margin_extra{0.0};
  SurfaceLimits limits;   // deflection + rate limits (from the nominal)
  bool rate_box{true};    // include the one-step rate limit in the input box
  double dt{0.01};
};

// Per-row diagnostics.
struct SixDofCbfRow {
  std::string name;
  double h{0}, Lfh{0}, Lf2h{0};   // barrier and its drift Lie stack
  SurfVec a{SurfVec::Zero()};     // a . U <= rhs (unnormalized)
  double rhs{0}, margin{0};
  SurfVec D1{SurfVec::Zero()};    // L_{F_k} h            (force-term leak, order 1)
  SurfVec D2{SurfVec::Zero()};    // L_F L_f h + L_f L_F h (order 2)
  SurfVec D3{SurfVec::Zero()};    // order 3 (degree-3 rows only)
  double Lf3h{0};                 // degree-3 rows only
  int degree{2};
  bool hard{false}, finite{true};
};
struct SixDofCbfDiag {
  bool active{false}, best_effort{false}, solved{true};
  int qp_status{0};         // solver status of the last solve attempt
  double min_h{1e30};       // min barrier value over the rows this step
  double margin_max{0};     // largest robustness margin applied this step
  double t_us{0};           // wall time of the Lie derivatives + QP [us]
  std::vector<SixDofCbfRow> rows;
  // Per-family barrier values for logging (min over the +/- pair).
  double h_alpha{0}, h_V{0}, h_phi{0}, h_beta{0}, h_terrain{0};
};

// The control-affine barrier model at one step.
class SixDofCbfModel {
 public:
  SixDofCbfModel(const BeaverDynamics& dyn, const SurfVec& u_prev, double dT,
                 const Eigen::Vector3d& W_earth = Eigen::Vector3d::Zero())
      : dyn_(dyn), up_(u_prev), dT_(dT), We_(W_earth) {}

  // Air-relative body velocity (ua, va, wa) at X: inertial minus the wind
  // rotated into the body frame (same DCM rows as BeaverDynamics::xdot).
  template <class T>
  std::array<T, 3> airRel(const std::array<T, NX>& X) const {
    using std::cos;
    using std::sin;
    const T ct = cos(X[THETA]), st = sin(X[THETA]);
    const T cp = cos(X[PHI]), sp = sin(X[PHI]);
    const T cy = cos(X[PSI]), sy = sin(X[PSI]);
    const double Wn = We_[0], We = We_[1], Wd = -We_[2];
    const T Wbx = ct * cy * Wn + ct * sy * We - st * Wd;
    const T Wby = (sp * st * cy - cp * sy) * Wn + (sp * st * sy + cp * cy) * We + sp * ct * Wd;
    const T Wbz = (cp * st * cy + sp * sy) * Wn + (cp * st * sy - sp * cy) * We + cp * ct * Wd;
    return {X[U] - Wbx, X[V] - Wby, X[W] - Wbz};
  }

  // Moment columns of the surfaces at X: rows P, Q, R of g (NX x 3).
  template <class T>
  std::array<std::array<T, 3>, NUS> momentColumns(const std::array<T, NX>& X) const {
    using std::atan;
    using std::sqrt;
    const BeaverAeroCoef& k = dyn_.config().aero;
    const auto ar = airRel<T>(X);
    const T Vt = sqrt(ar[0] * ar[0] + ar[1] * ar[1] + ar[2] * ar[2] + T(1e-12));
    const T alpha = atan(ar[2] / ar[0]);
    const T qS = (0.5 * dyn_.rho()) * Vt * Vt * BeaverGeom::S;
    const double Ix = BeaverGeom::Ix, Iy = BeaverGeom::Iy, Iz = BeaverGeom::Iz,
                 Ixz = BeaverGeom::Ixz, Gamma = Ix * Iz - Ixz * Ixz;
    const double b = BeaverGeom::b, c = BeaverGeom::c;
    auto solve = [&](const T& L, const T& N) {
      return std::array<T, 3>{(Iz * L + Ixz * N) / Gamma, T(0.0),
                              (Ix * N + Ixz * L) / Gamma};
    };
    std::array<std::array<T, 3>, NUS> g;
    g[0] = {T(0.0), qS * c * k.Cm_de / Iy, T(0.0)};                        // de
    g[1] = solve(qS * b * (k.Cl_da + k.Cl_daa * alpha), qS * b * k.Cn_da);  // da
    g[2] = solve(qS * b * k.Cl_dr, qS * b * k.Cn_dr);                       // dr
    return g;
  }

  // Drift f(X) = plant(X, u_prev, dT) - g(X) u_prev. Scalar-generic: T may be
  // double, autodiff::dual, or a Taylor jet.
  template <class T>
  std::array<T, NX> operator()(const std::array<T, NX>& X) const {
    std::array<T, NU> u;
    u[DE] = T(up_[0]); u[DA] = T(up_[1]); u[DR] = T(up_[2]); u[DT] = T(dT_);
    const auto ar = airRel<T>(X);
    std::array<T, NX> xd = dyn_.xdotCoreT<T>(X, ar[0], ar[1], ar[2], u);
    const auto g = momentColumns<T>(X);
    for (int k = 0; k < NUS; ++k) {
      xd[P] = xd[P] - g[k][0] * up_[k];
      xd[Q] = xd[Q] - g[k][1] * up_[k];
      xd[R] = xd[R] - g[k][2] * up_[k];
    }
    return xd;
  }

  // g columns as full NX vectors (double) at x0.
  std::array<std::array<double, NX>, NUS> gColumns(const std::array<double, NX>& x0) const;
  // Neglected direct-force columns F_k (per rad of surface), scalar-generic.
  template <class T>
  std::array<std::array<T, NX>, NUS> forceColumnsT(const std::array<T, NX>& X) const {
    using std::asin;
    using std::atan;
    using std::sqrt;
    const BeaverAeroCoef& k = dyn_.config().aero;
    const auto ar = airRel<T>(X);
    const T Vt = sqrt(ar[0] * ar[0] + ar[1] * ar[1] + ar[2] * ar[2] + T(1e-12));
    const T alpha = atan(ar[2] / ar[0]);
    const T beta = asin(ar[1] / Vt);
    const T qS = (0.5 * dyn_.rho()) * Vt * Vt * BeaverGeom::S;
    const double m = BeaverGeom::mass;
    std::array<std::array<T, NX>, NUS> F;
    for (auto& c : F) for (auto& v : c) v = T(0.0);
    F[0][W] = qS * (k.Cz_de + k.Cz_deb2 * beta * beta) / m;   // elevator lift
    F[1][V] = qS * k.Cy_da / m;                                // aileron side force
    F[2][U] = qS * k.Cx_dr / m;                                // rudder drag
    F[2][V] = qS * (k.Cy_dr + k.Cy_dra * alpha) / m;           // rudder side force
    return F;
  }
  std::array<std::array<double, NX>, NUS> forceColumns(const std::array<double, NX>& x0) const {
    return forceColumnsT<double>(x0);
  }

  const SurfVec& uPrev() const { return up_; }
  double dT() const { return dT_; }
  const Eigen::Vector3d& wind() const { return We_; }

 private:
  const BeaverDynamics& dyn_;
  SurfVec up_;
  double dT_;
  Eigen::Vector3d We_;
};

class SixDofCbfFilter {
 public:
  SixDofCbfFilter(const BeaverDynamics& dyn, const SixDofCbfConfig& cfg);
  ~SixDofCbfFilter();

  // Filter the nominal surfaces at state x, with the throttle dT the nominal
  // chose this step and the surfaces applied last step. Returns the safe
  // surfaces (== U_nom when nothing binds).
  SurfVec filter(const SurfVec& U_nom, const StateVec& x, double dT,
                 const SurfVec& u_prev, SixDofCbfDiag* diag = nullptr,
                 const Eigen::Vector3d& W_earth = Eigen::Vector3d::Zero()) const;

  // Build the rows only (no QP) -- for tests and the feasibility sweep.
  std::vector<SixDofCbfRow> rows(const StateVec& x, double dT,
                                 const SurfVec& u_prev,
                                 const Eigen::Vector3d& W_earth = Eigen::Vector3d::Zero()) const;
  // Is the HARD row set plus the input box feasible at this state?
  bool hardSetFeasible(const StateVec& x, double dT, const SurfVec& u_prev,
                       const Eigen::Vector3d& W_earth = Eigen::Vector3d::Zero()) const;

  const SixDofCbfConfig& config() const { return cfg_; }

 private:
  const BeaverDynamics& dyn_;
  SixDofCbfConfig cfg_;
  std::unique_ptr<QPSolver> solver_;
};

}  // namespace autoland
