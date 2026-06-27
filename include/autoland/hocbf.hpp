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
// and drag depend only on (V, gamma) -- already in b -- so degree 3 is preserved.
// Build via makeDescentBarrier() so the frozen-aero constants get filled in.
struct DescentBarrier {
  double v_safe2{0};
  double rho{0}, Sref{0}, mass{0}, g{0};  // environment / geometry
  double CLmax{0}, CDmaxlift{0};          // stall ceiling + its drag coefficient
  template <class T>
  T operator()(const std::array<T, NXA>& X) const {
    using std::cos;
    using std::sin;
    using std::sqrt;
    const T V = X[LV], gam = X[LGAM];
    const T a_brk = (0.5 * rho * Sref / mass) * (V * V) *
                        (CLmax * cos(gam) - CDmaxlift * sin(gam)) - g;
    return V * sin(gam) + sqrt(v_safe2 + (2.0 * a_brk) * X[LH]);
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

struct AirspeedBarrier {  // b_V = V - V_min
  double Vmin{0};
  template <class T>
  T operator()(const std::array<T, NXA>& X) const {
    return X[LV] - Vmin;
  }
};

// Upper airspeed barrier b = V_max - V >= 0 (over-speed / high-energy water-
// impact / structural-limit protection). Same relative degree (3) and control-
// affine structure as the lower airspeed barrier -- signs flip -- so it reuses
// the existing machinery (barrierLie<3> -> hocbfRow). Wired into LonCBFFilter via
// the airspeed_upper flag + Vmax_air in LonCBFConfig (mirrors the lower barrier).
struct AirspeedUpperBarrier {  // b = V_max - V
  double Vmax{0};
  template <class T>
  T operator()(const std::array<T, NXA>& X) const {
    return Vmax - X[LV];
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
