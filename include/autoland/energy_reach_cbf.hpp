#pragma once
#include <array>
#include <cmath>
#include <memory>
#include <Eigen/Dense>
#include "autoland/beaver_dynamics.hpp"
#include "autoland/lie_taylor.hpp"
#include "autoland/qp_solver.hpp"

// =============================================================================
// Reach-the-target ENERGY CBF on a 3-DOF point-mass glide model (EnergyCBF
// slide deck, 2026-09). The barrier asks: does the vehicle hold enough total
// specific energy to glide to the target and arrive at V_min?
//
//   h_E = (1/2 V^2 + g z)  -  (1/2 V_min^2 + g d / (L/D))
//
// with BOTH improvements of the reformulation:
//   * FULL alpha-dependent glide ratio  L/D = C_La*at / (C_D0 + C_Da2*at^2)
//     (at = alpha - alpha_0L, alpha measured from the zero-lift line), so the
//     barrier knows how changing alpha changes the energy needed -- instead of
//     the old constant (L/D)_max that pushed toward the min-drag dive.
//   * HEADING-dependent circular-arc distance d = R*phi/sin(phi), with
//     phi = chi - atan2(dE, dN) the heading error to the target bearing and R
//     the straight-line NE distance (the arc through the current heading that
//     ends at the target). This brings chi into the barrier, which drops the
//     sideslip control to relative degree 1 via chidot = qbar S C_Yb b/(m V cg).
//
// MODEL (3-DOF point mass, unpowered/skid-to-turn -- the deck's reduced model):
//   State  X = [V, gamma, chi, N, E, z, alpha]        (NXER = 7)
//   Control U = [alpha_dot, beta]                     (NUER = 2)
//   Vdot     = -D/m - g sin(gamma)
//   gammadot =  L/(mV) - g cos(gamma)/V               (wings level)
//   chidot   =  qbar S C_Yb beta / (m V cos(gamma))   (side-force turn)
//   Ndot, Edot, zdot = kinematics;  alphadot = u_alpha
// Thrust is ABSENT: the drift models gliding flight. Fit the polar at the idle
// (windmilling) throttle so C_D0 absorbs the propeller drag; under power the
// true energy decays slower than modeled and the barrier is conservative.
//
// The whole CBF is uniform relative degree 1, enforced as the single QP row
//   L_f h + L_g h . U >= -ck * h
// with exact Lie derivatives from the flow Taylor-jet engine (lie_taylor.hpp),
// per house rules -- the slide-7 hand algebra is never transcribed.
//
// SMOOTH GUARDS (all C^1, autodiff-safe):
//   * at is floored at alpha_lo by a smooth-max (sqrt form): as at -> 0 the
//     energy needed to glide anywhere diverges (no lift), and at < 0 flips the
//     sign; the floor caps the barrier's pessimism while keeping d(h)/d(alpha)
//     exact above the floor.
//   * phi is saturated through phi_sat*tanh(phi/phi_sat) (phi_sat < pi), which
//     bounds the arc length at ~5R for a target behind the vehicle instead of
//     the sin(phi) -> 0 blowup, while keeping the gradient alive for the QP.
//   * phi/sin(phi) switches to its series 1 + phi^2/6 + 7 phi^4/360 below
//     |phi| = 0.05 (removable singularity; the seam mismatch is ~1e-11).
//   * within sqrt(dN^2+dE^2) < d_near of the target the arc collapses to the
//     straight-line distance (the bearing -- and atan2 -- is ill-defined at
//     the target; heading no longer matters there).
// =============================================================================
namespace autoland {

constexpr int NXER = 7;
enum ReachState { ERV = 0, ERGAM = 1, ERCHI = 2, ERN = 3, ERE = 4, ERZ = 5, ERAL = 6 };
constexpr int NUER = 2;
enum ReachCtrl { ERUA = 0, ERBE = 1 };

using ReachStateArr = std::array<double, NXER>;

// Physical constants, fitted glide polar, target, and the smoothing guards.
struct ReachParams {
  double m{2288.231}, S{23.23}, rho{1.225}, g{9.80665};
  // Quadratic glide polar in at = alpha - alpha0L (fitBeaverReachPolar).
  double CLa{5.2}, CD0{0.05}, CDa2{1.5}, alpha0L{0.0};
  double CYb{-0.7678};      // side-force slope [1/rad] (skid-turn authority)
  // Target (north/east, z datum at the target altitude) and arrival speed.
  double Nt{0.0}, Et{0.0};
  double V_min{30.0};
  // Guards.
  double alpha_lo{0.5 * M_PI / 180.0};      // smooth floor on at [rad]
  double alpha_smooth{0.5 * M_PI / 180.0};  // smooth-max width [rad]
  double phi_sat{2.6};   // tanh saturation on the arc half-angle [rad] (< pi)
  double R_eps{1.0};     // smooth floor inside the straight-line distance [m]
  double d_near{3.0};    // arc -> straight-line switch radius at the target [m]
  double V_eps{0.5};     // smooth airspeed floor [m/s] (never active in flight)
};

// --- Point-mass drift f(X) (U = 0), templated for the Lie jet engine. --------
struct ReachDrift {
  ReachParams p;
  explicit ReachDrift(const ReachParams& prm) : p(prm) {}

  template <class T>
  std::array<T, NXER> operator()(const std::array<T, NXER>& X) const {
    using std::sin;
    using std::cos;
    using std::sqrt;
    const T V = sqrt(X[ERV] * X[ERV] + p.V_eps * p.V_eps);
    const T gam = X[ERGAM], chi = X[ERCHI];
    const T at = X[ERAL] - p.alpha0L;  // RAW alpha-from-zero-lift: the drift
                                       // models the actual vehicle, floor-free
    const T qbar = (0.5 * p.rho) * V * V;
    const T Lift = qbar * p.S * (p.CLa * at);
    const T Drag = qbar * p.S * (p.CD0 + p.CDa2 * at * at);
    const T sg = sin(gam), cg = cos(gam);
    std::array<T, NXER> dX;
    dX[ERV] = (-1.0) * Drag * (1.0 / p.m) - p.g * sg;
    dX[ERGAM] = Lift / (p.m * V) - (p.g * cg) / V;
    dX[ERCHI] = T(0.0);
    dX[ERN] = V * cg * cos(chi);
    dX[ERE] = V * cg * sin(chi);
    dX[ERZ] = V * sg;
    dX[ERAL] = T(0.0);
    return dX;
  }
};

// Control matrix g(X): u_alpha is the alpha integrator; beta turns chi through
// the side force. Evaluated at the (double) state.
inline Eigen::Matrix<double, NXER, NUER> reachGMatrix(const ReachParams& p,
                                                      const ReachStateArr& X) {
  Eigen::Matrix<double, NXER, NUER> G = Eigen::Matrix<double, NXER, NUER>::Zero();
  const double V = std::sqrt(X[ERV] * X[ERV] + p.V_eps * p.V_eps);
  const double qbar = 0.5 * p.rho * V * V;
  G(ERAL, ERUA) = 1.0;
  G(ERCHI, ERBE) = qbar * p.S * p.CYb / (p.m * V * std::cos(X[ERGAM]));
  return G;
}

// --- The reach-energy barrier h_E(X), templated for the Lie jet engine. ------
struct ReachEnergyBarrier {
  ReachParams p;
  explicit ReachEnergyBarrier(const ReachParams& prm) : p(prm) {}

  // phi/sin(phi) with the series branch near the removable singularity.
  template <class T>
  T sincRatio(const T& phi) const {
    using std::sin;
    if (std::abs(scalarValue(phi)) < 0.05) {
      const T p2 = phi * phi;
      return 1.0 + p2 * (1.0 / 6.0) + (p2 * p2) * (7.0 / 360.0);
    }
    return phi / sin(phi);
  }

  // Circular-arc distance to the target through the current heading.
  template <class T>
  T arcDistance(const T& N, const T& E, const T& chi) const {
    using std::sin;
    using std::cos;
    using std::sqrt;
    using std::tanh;
    using std::atan2;
    const T dN = p.Nt - N, dE = p.Et - E;
    const T R = sqrt(dN * dN + dE * dE + p.R_eps * p.R_eps);
    if (scalarValue(dN) * scalarValue(dN) + scalarValue(dE) * scalarValue(dE) <
        p.d_near * p.d_near)
      return R;  // at the target the bearing (and the arc) is meaningless
    // phi = chi - atan2(dE, dN), wrapped to (-pi, pi] by construction.
    const T dpar = dN * cos(chi) + dE * sin(chi);
    const T dperp = dN * sin(chi) - dE * cos(chi);
    const T phi = atan2(dperp, dpar);
    const T phis = p.phi_sat * tanh(phi * (1.0 / p.phi_sat));
    return R * sincRatio(phis);
  }

  // Smooth-max(at, alpha_lo) with width alpha_smooth.
  template <class T>
  T alphaFloored(const T& alpha) const {
    using std::sqrt;
    const T s = (alpha - p.alpha0L) - p.alpha_lo;
    return p.alpha_lo +
           0.5 * (s + sqrt(s * s + p.alpha_smooth * p.alpha_smooth));
  }

  template <class T>
  T operator()(const std::array<T, NXER>& X) const {
    const T V = X[ERV];
    const T ats = alphaFloored(X[ERAL]);
    const T d = arcDistance(X[ERN], X[ERE], X[ERCHI]);
    const T DoverL = (p.CD0 + p.CDa2 * ats * ats) / (p.CLa * ats);
    return 0.5 * V * V + p.g * X[ERZ] - 0.5 * p.V_min * p.V_min -
           p.g * d * DoverL;
  }
};

// --- Exact RD-1 Lie bundle: {h, L_f h} and the control row L_g h. ------------
struct ReachLie {
  double h{0}, Lfh{0};
  double A_ua{0}, A_be{0};  // L_g h columns (alpha_dot, beta)
};
inline ReachLie reachLie(const ReachParams& p, const ReachStateArr& X) {
  const ReachDrift f(p);
  const ReachEnergyBarrier b(p);
  ReachLie out;
  const std::array<double, 2> Lf = lieDrift<1, NXER>(f, b, X);
  out.h = Lf[0];
  out.Lfh = Lf[1];
  const Eigen::Matrix<double, NXER, NUER> G = reachGMatrix(p, X);
  for (int c = 0; c < NUER; ++c) {
    std::array<double, NXER> dir;
    for (int i = 0; i < NXER; ++i) dir[i] = G(i, c);
    const double v = lieAlong<1, NXER>(f, b, X, dir);
    (c == ERUA ? out.A_ua : out.A_be) = v;
  }
  return out;
}

// --- CBF-QP filter -----------------------------------------------------------
// z = [u_alpha, beta]:
//   min 1/2 [ w_ua (u_a - u_a_nom)^2 + w_be (b - b_nom)^2 ]
//   s.t.  -A_ua u_a - A_be b <= L_f h + ck h              (energy reach, HARD)
//         u_a in [max(-ua_max, -c_alpha (alpha - alpha_cmd_min)),
//                 min( ua_max,  c_alpha (alpha_cmd_max - alpha))]  (AoA band,
//                            RD-1 rows folded into the variable bounds)
//         b   in [-beta_max, beta_max]
// The energy row is HARD; when it is infeasible against the input bounds (the
// set is genuinely not holdable -- e.g. a target abeam at low altitude) the
// filter falls back to the CLOSED-FORM minimum-violation action, railing each
// control toward its max-hdot bound (u_i = bound in the direction of A_i).
// A quadratic slack was tried and rejected: pricing the TRUE violation needs
// w*s^2 on the normalized row (s ~ 1e4), which wrecks OSQP's conditioning,
// while pricing the normalized violation makes the row soft as mush.
// Tracking weights are RANGE-SCALED (w = 1/u_max^2) so a full-scale deflection
// of each channel is priced equally (cf. lon_cbf_filter.hpp).
struct EnergyReachCfg {
  bool enabled{true};
  ReachParams prm;
  double ck{0.3};            // linear class-K gain [1/s]
  double c_alpha{2.0};       // class-K for the AoA band rows [1/s]
  double ua_max{5.0 * M_PI / 180.0};    // |alpha_dot| bound [rad/s]
  double beta_max{10.0 * M_PI / 180.0}; // |beta| bound [rad]
  // Polar-validity AoA band. NOTE the Beaver glides FAST at idle: the flaps-up
  // trim alpha is ~10.5 deg at 40 m/s and ~13.5 deg at 35 m/s, so the ceiling
  // must clear the whole approach envelope (attached flow holds to ~15 deg).
  double alpha_cmd_min{-4.0 * M_PI / 180.0};
  double alpha_cmd_max{14.0 * M_PI / 180.0};
};

struct EnergyReachDiag {
  double h{0}, Lfh{0}, A_ua{0}, A_be{0};
  double d{0}, phi{0}, LD{0};  // arc distance, heading error, instantaneous L/D
  double ua{0}, be{0};
  double viol{0};           // energy-row deficit at the returned u [m^2/s^3]
  bool active{false};       // the filter moved the command
  bool best_effort{false};  // row infeasible -> closed-form max-hdot action
  bool finite{true};        // false -> non-finite Lie data (pass-through)
};

class EnergyReachFilter {
 public:
  EnergyReachFilter() : solver_(std::make_unique<OsqpSolver>()) {}
  explicit EnergyReachFilter(EnergyReachCfg cfg,
                             std::unique_ptr<QPSolver> solver = nullptr)
      : cfg_(cfg),
        solver_(solver ? std::move(solver)
                       : std::unique_ptr<QPSolver>(std::make_unique<OsqpSolver>())) {}

  // Filter (u_alpha_nom, beta_nom) at point-mass state X. Returns [u_a, beta].
  Eigen::Vector2d filter(double ua_nom, double be_nom, const ReachStateArr& X,
                         EnergyReachDiag* diag = nullptr) const;

  const EnergyReachCfg& config() const { return cfg_; }

 private:
  EnergyReachCfg cfg_;
  std::unique_ptr<QPSolver> solver_;
};

// --- Idle-glide polar fit from the Beaver polynomials ------------------------
// Least-squares of the wind-axis lift/drag over an alpha sweep at zero rates
// and zero surface deflections, at the plant's flap setting and the dpt of
// throttle dT at V_ref (dT at idle -> windmilling-propeller drag folded into
// C_D0):  C_L ~ CLa*(alpha - alpha0L),  C_D ~ CD0 + CDa2*(alpha - alpha0L)^2.
struct ReachPolarFit {
  double CLa{0}, CD0{0}, CDa2{0}, alpha0L{0};
  double alpha_star{0};  // best-L/D alpha-from-zero-lift: sqrt(CD0/CDa2)
  double LD_max{0};      // CLa / (2 sqrt(CD0*CDa2))
  double rms_CL{0}, rms_CD{0};
};
ReachPolarFit fitBeaverReachPolar(const BeaverDynamics& dyn, double V_ref,
                                  double dT, double alpha_min = -3.0 * M_PI / 180.0,
                                  double alpha_max = 14.0 * M_PI / 180.0,
                                  int n = 61);

}  // namespace autoland
