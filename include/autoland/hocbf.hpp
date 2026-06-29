#pragma once
#include <array>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "autoland/lie_taylor.hpp"
#include "autoland/lon_augmented.hpp"

// =============================================================================
// High-Order Control Barrier Functions for the augmented longitudinal system.
//
// Two pieces, kept separate (per design):
//  1) Per-barrier EXACT Lie-derivative provider `barrierLie<R>()` -> the drift
//     stack {b, L_f b, ..., L_f^R b} and the control row L_g L_f^{R-1} b, via the
//     flow Taylor-jet engine in lie_taylor.hpp (no finite differences).
//  2) A standalone, generic HOCBF constraint builder `hocbfRow()` that turns a
//     Lie-derivative bundle + class-K functions (depth r) into the QP row.
//     Default class-K is linear, alpha_i(s) = c_i * s, with the c_i as inputs.
//
// QP row convention (OSQP A U <= b form):  a . U <= rhs.
// =============================================================================
namespace autoland {

// --- Barriers (templated on element type so the Lie engine can autodiff them) -
// Descent (soft-landing) barrier  b = V sin(gamma) + sqrt(v_safe^2 + 2 a_brk h),
// with the braking acceleration now SPEED / PATH-ANGLE dependent (lift + drag +
// gravity), not a constant:
//   a_brk(V,gamma) = (rho V^2 S / 2m)[CLmax cos(gamma) - CDmaxlift sin(gamma)] - g
// Thrust is deliberately omitted: putting the thrust state T / pitch theta into b
// would drop the barrier's relative degree below 3 and break the augmented HOCBF
// alignment (that theta-coupling is the section 3.3 mixed-degree problem). Lift
// and drag depend only on (V, gamma) -- already in b -- so degree 3 is preserved
// (h enters only the radicand, at relative degree 4, so it does not change that).
// Build via makeDescentBarrier() so the frozen-aero constants get filled in.
struct DescentBarrier {
  double v_safe2{0};
  double rho{0}, Sref{0}, mass{0}, g{0};  // environment / geometry
  double CLmax{0}, CDmaxlift{0};          // stall ceiling + its drag coefficient
  // Numerical regularizer [m^2/s^2] on the kinematic radicand -- NOT a tuning
  // knob. When the braking envelope momentarily has no real solution
  // (v_safe2 + 2*a_brk*h < 0 -- e.g. the touchdown sample where h dips below 0,
  // or a_brk <= 0 at low airspeed) the raw radicand goes negative and sqrt()
  // returns NaN. The filter's finiteness guard would then SILENTLY DROP this row
  // -- deleting the only hard sink-rate guarantee with no trace, and the NaN is
  // masked in the psi diagnostics (std::min ignores NaN). The smooth positive
  // floor below keeps b finite and C-infinity, so the exact-Lie Taylor jet and
  // the QP row stay well-posed and the barrier degrades gracefully (it stays in
  // the QP and drives the strongest available flare) instead of vanishing.
  double eps_r{0.01};
  template <class T>
  T operator()(const std::array<T, NXA>& X) const {
    using std::cos;
    using std::sin;
    using std::sqrt;
    const T V = X[LV], gam = X[LGAM];
    const T a_brk = (0.5 * rho * Sref / mass) * (V * V) *
                        (CLmax * cos(gam) - CDmaxlift * sin(gam)) - g;
    const T arg = v_safe2 + (2.0 * a_brk) * X[LH];
    // arg_pos = 0.5(arg + sqrt(arg^2 + 4 eps_r^2)) > 0 for all arg; -> arg for
    // arg >> eps_r, -> 0+ for arg << -eps_r. eps_r is small enough that the
    // normal-flight bias (~eps_r^2/arg, arg ~ 1e2..1e3) is negligible.
    const T arg_pos = 0.5 * (arg + sqrt(arg * arg + (4.0 * eps_r * eps_r)));
    return V * sin(gam) + sqrt(arg_pos);
  }
};

// Drag coefficient at the max-lift condition (q = 0), used by a_brk's drag term.
// The inviscid table has no true stall, so we linearly extrapolate to alpha_max
// where the rotated C_L = CLmax (one Newton-style step about alpha = 0) and
// evaluate the rotated C_D there. Cheap and approximate; the drag term itself is
// a small (~sin gamma) correction.
inline double cdAtMaxLift(const AeroLocal& a, double CLmax, double V) {
  const double mach = V / a.a_sound;
  const double CFx0 = a.off_CFx + a.dMach_CFx * mach;  // alpha = 0, q = 0
  const double CFz0 = a.off_CFz + a.dMach_CFz * mach;
  const double CL0 = CFz0;                     // -CFx0 sin0 + CFz0 cos0
  const double dCL_da0 = a.dAlpha_CFz - CFx0;  // d/dalpha[-CFx sin + CFz cos]|_0
  const double amax =
      (std::abs(dCL_da0) > 1e-9) ? (CLmax - CL0) / dCL_da0 : 0.0;
  const double CFx = a.off_CFx + a.dAlpha_CFx * amax + a.dMach_CFx * mach;
  const double CFz = a.off_CFz + a.dAlpha_CFz * amax + a.dMach_CFz * mach;
  return CFx * std::cos(amax) + CFz * std::sin(amax) + a.parasite_CD0;
}

// Build the descent barrier from the frozen aero + config, filling in the
// a_brk(V,gamma) constants (including the max-lift drag coefficient at V).
inline DescentBarrier makeDescentBarrier(const AeroLocal& a, double v_safe,
                                         double CLmax, double V) {
  DescentBarrier b;
  b.v_safe2 = v_safe * v_safe;
  b.rho = a.rho;
  b.Sref = a.Sref;
  b.mass = a.mass;
  b.g = a.g;
  b.CLmax = CLmax;
  b.CDmaxlift = cdAtMaxLift(a, CLmax, V);
  return b;
}

// Altitude-scheduled velocity limit: the surface value V_ground ramping smoothly
// up to the aloft value V_air over scale height hs:
//   Vlim(h) = V_ground + (V_air - V_ground)(1 - exp(-h/hs)).
// Vlim(0) = V_ground, -> V_air as h -> inf, and it is C-infinity so the exact-Lie
// Taylor jet stays well-posed. Set V_ground = V_air to recover a constant limit.
template <class T>
T schedVlim(const T& h, double V_air, double V_ground, double hs) {
  using std::exp;
  // Multiply by the reciprocal scale (Taylor has no operator/ by a scalar).
  return V_ground + (V_air - V_ground) * (1.0 - exp(h * (-1.0 / hs)));
}

// Lower airspeed barrier  b = V - Vmin(h), with an ALTITUDE-SCHEDULED stall floor:
// aloft the full stall-margin floor Vmin_air applies; near the surface it relaxes
// toward Vmin_ground so the aircraft may bleed airspeed to (or below) stall for a
// slow touchdown. h enters the augmented system at relative degree 4 -- one slower
// than V's degree 3 -- so the barrier stays degree 3 and the control row is
// UNCHANGED; only the drift (rhs) is reshaped by the schedule. Set Vmin_ground =
// Vmin_air (and any hs) for the original constant-floor barrier.
struct AirspeedBarrier {  // b = V - Vmin(h)
  double Vmin_air{0};     // stall-margin floor aloft [m/s]
  double Vmin_ground{0};  // relaxed floor at the surface [m/s] (<= Vmin_air)
  double h_sched{1.0};    // schedule scale height [m]
  AirspeedBarrier() = default;
  // Single-arg form keeps the original CONSTANT floor (Vmin everywhere).
  explicit AirspeedBarrier(double Vmin)
      : Vmin_air(Vmin), Vmin_ground(Vmin), h_sched(1.0) {}
  AirspeedBarrier(double v_air, double v_ground, double hs)
      : Vmin_air(v_air), Vmin_ground(v_ground), h_sched(hs) {}
  template <class T>
  T operator()(const std::array<T, NXA>& X) const {
    return X[LV] - schedVlim<T>(X[LH], Vmin_air, Vmin_ground, h_sched);
  }
};

// Upper airspeed barrier b = Vmax(h) - V >= 0 (over-speed / high-energy water-
// impact / structural-limit protection), with the same ALTITUDE-SCHEDULED ceiling:
// the over-speed limit Vmax_air aloft tightens toward Vmax_ground near the surface,
// commanding the aircraft to decelerate so it touches down at a lower horizontal
// speed. Same relative degree (3) and control-affine structure as the lower barrier
// -- signs flip -- so it reuses the machinery (barrierLie<3> -> hocbfRow). Set
// Vmax_ground = Vmax_air for the original constant ceiling.
struct AirspeedUpperBarrier {  // b = Vmax(h) - V
  double Vmax_air{0};     // over-speed ceiling aloft [m/s]
  double Vmax_ground{0};  // tightened ceiling at the surface [m/s] (<= Vmax_air)
  double h_sched{1.0};    // schedule scale height [m]
  AirspeedUpperBarrier() = default;
  // Single-arg form keeps the original CONSTANT ceiling (Vmax everywhere).
  explicit AirspeedUpperBarrier(double Vmax)
      : Vmax_air(Vmax), Vmax_ground(Vmax), h_sched(1.0) {}
  AirspeedUpperBarrier(double v_air, double v_ground, double hs)
      : Vmax_air(v_air), Vmax_ground(v_ground), h_sched(hs) {}
  template <class T>
  T operator()(const std::array<T, NXA>& X) const {
    return schedVlim<T>(X[LH], Vmax_air, Vmax_ground, h_sched) - X[LV];
  }
};

// --- Lie-derivative bundle for a relative-degree-R barrier --------------------
template <int R>
struct BarrierLie {
  std::array<double, R + 1> Lf{};            // {b, L_f b, ..., L_f^R b}
  Eigen::Matrix<double, 1, NUA> LgLf;        // L_g L_f^{R-1} b  (control row)
};

// Compute the exact Lie derivatives of barrier `b` at X0 over the augmented
// drift f and control matrix g (relative degree R for every control column).
template <int R, class Barrier>
BarrierLie<R> barrierLie(const AeroLocal& a, const Barrier& b,
                         const LonStateVec& X0) {
  const LonDrift f(a);
  const std::array<double, NXA> x0 = toArray(X0);
  BarrierLie<R> out;
  out.Lf = lieDrift<R, NXA>(f, b, x0);

  const LonGMat G = gMatrix(a, X0);
  for (int c = 0; c < NUA; ++c) {
    std::array<double, NXA> dir;
    for (int i = 0; i < NXA; ++i) dir[i] = G(i, c);
    out.LgLf(0, c) = lieAlong<R, NXA>(f, b, x0, dir);
  }
  return out;
}

// --- Generic HOCBF constraint row (linear class-K) ---------------------------
struct HocbfRow {
  Eigen::Matrix<double, 1, NUA> a{Eigen::Matrix<double, 1, NUA>::Zero()};  // a.U <= rhs
  double rhs{0};
};

// Build the degree-r constraint from the drift Lie stack `Lf` (size r+1), the
// control row `LgLf` (= L_g L_f^{r-1} b), and linear class-K gains `c` (size r).
//
//   psi_r = L_f^r b + sum_{j<r} e_{r,j} L_f^j b  +  (L_g L_f^{r-1} b) . U  >= 0
//   with e_{r,j} the elementary-symmetric polynomials in c (built by recursion).
// QP form:  -(L_g L_f^{r-1} b) . U  <=  sum_j coeff_j L_f^j b.
inline HocbfRow hocbfRow(const std::vector<double>& Lf,
                         const Eigen::Matrix<double, 1, NUA>& LgLf,
                         const std::vector<double>& c) {
  const int r = static_cast<int>(c.size());
  // psi coefficient recursion over {L_f^0 b .. L_f^r b}: coeff = shift(coeff)+c_k*coeff.
  std::vector<double> coeff(1, 1.0);  // psi_0 = b
  for (int k = 0; k < r; ++k) {
    std::vector<double> next(coeff.size() + 1, 0.0);
    for (int j = 0; j < static_cast<int>(coeff.size()); ++j) {
      next[j + 1] += coeff[j];            // L_f shift
      next[j] += c[k] * coeff[j];          // c_k * psi_{k-1}
    }
    coeff.swap(next);
  }
  HocbfRow row;
  double drift = 0.0;
  for (int j = 0; j <= r; ++j) drift += coeff[j] * Lf[j];
  row.a = -LgLf;
  row.rhs = drift;
  return row;
}

// --- Actuator barriers (degree 2 wrt Tddot; closed form, doc section 3.4) -----
// Min thrust  b1 = T >= 0:        Tddot >= -(c11+c12)Tdot - c11 c12 T
//   -> row:  -Tddot <= (c11+c12)Tdot + c11 c12 T
// Max thrust  b2 = Tmax - T >= 0: Tddot <= -(c21+c22)Tdot + c21 c22 (Tmax - T)
inline HocbfRow thrustMinRow(const LonStateVec& X, double c11, double c12) {
  HocbfRow row;
  row.a(0, LTDDOT) = -1.0;
  row.rhs = (c11 + c12) * X[LTDOT] + c11 * c12 * X[LT];
  return row;
}
inline HocbfRow thrustMaxRow(const LonStateVec& X, double Tmax, double c21,
                             double c22) {
  HocbfRow row;
  row.a(0, LTDDOT) = 1.0;
  row.rhs = -(c21 + c22) * X[LTDOT] + c21 * c22 * (Tmax - X[LT]);
  return row;
}

}  // namespace autoland
