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

// --- Barriers (templated on element type so the Lie engine can autodiff them) --
// Recovery-set barriers: stall (AoA upper), nose-up (attitude lower), and a total-
// specific-energy ceiling. All are linear/polynomial in the state => C-infinity
// (no eps floors needed); relative degrees follow from the g-matrix (elevator only
// enters via G(LQ,LDE), Tddot only via G(LTDOT,LTDDOT)). The impact-load barrier
// lives in impact_barrier.hpp; the thrust actuator barriers are below.

// Stall (AoA upper bound). b = alpha_lim - (theta - gamma), alpha_lim = alpha_stall
// - margin. Active the WHOLE flight; as alpha -> alpha_lim the HOCBF drives the
// elevator nose-down (saturating at the boundary), reproducing the pilot's low-
// altitude stall recovery. Relative degree 2 (elevator: theta->q->de); Tddot is
// absent from the degree-2 control row. Class-K size 2.
struct StallBarrier {  // b = alpha_lim - (theta - gamma)
  double alpha_lim{0};
  template <class T>
  T operator()(const std::array<T, NXA>& X) const {
    return alpha_lim - (X[LTH] - X[LGAM]);
  }
};
inline StallBarrier makeStallBarrier(double alpha_stall, double margin) {
  return StallBarrier{alpha_stall - margin};
}

// Nose-up (attitude lower bound). b = theta - theta_min. THETA-based (not alpha):
// it directly keeps the impact model's gate tau = theta - tau_keel > 0 valid (an
// AoA floor would not, since a steep gamma lets theta = alpha + gamma go negative),
// and gives a clean degree 2 (theta->q->de) with no gamma coupling. Set theta_min
// >= tau_keel. Gated to the final metres in the wiring. Relative degree 2, class-K
// size 2.
struct NoseUpBarrier {  // b = theta - theta_min
  double theta_min{0};
  template <class T>
  T operator()(const std::array<T, NXA>& X) const {
    return X[LTH] - theta_min;
  }
};
inline NoseUpBarrier makeNoseUpBarrier(double theta_min) {
  return NoseUpBarrier{theta_min};
}

// Total-specific-energy ceiling. E = 1/2 V^2 + g h, capped by a height-scheduled
// ceiling E_cap(h) = 1/2 V_td_max^2 + g_eff h (g_eff sized so E_cap(h0) >= E0, so
// it starts satisfied, and E_cap(0) = 1/2 V_td_max^2 bounds the touchdown speed):
//   b = E_cap(h) - E = 1/2 (V_td_max^2 - V^2) + (g_eff - g) h.
// POLYNOMIAL in (V,h) => C-infinity, NO sqrt: the equivalent airspeed cap
// V <= sqrt(V_td_max^2 + 2 (g_eff-g) h) is b>=0 solved for V, never formed. Only
// V,h enter b (not theta/T), so the degree-3 alignment is preserved. Relative
// degree 3 (both controls enter at the 3rd derivative); class-K size 3. Set up as
// a LOOSE never-exceed ceiling (see lon_cbf_filter.hpp), independent of glide path.
struct EnergyBarrier {
  double half_Vtd2{0};  // 1/2 V_td_max^2
  double dgeff{0};      // g_eff - g
  template <class T>
  T operator()(const std::array<T, NXA>& X) const {
    const T V = X[LV];
    return half_Vtd2 - (0.5 * V * V) + dgeff * X[LH];
  }
};
inline EnergyBarrier makeEnergyBarrier(const AeroLocal& a, double V_td_max,
                                       double g_eff) {
  EnergyBarrier b;
  b.half_Vtd2 = 0.5 * V_td_max * V_td_max;
  b.dgeff = g_eff - a.g;  // freeze true g from the aero/env
  return b;
}

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
