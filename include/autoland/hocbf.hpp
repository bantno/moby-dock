#pragma once
#include <array>
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
struct DescentBarrier {  // b = V sin(gamma) + sqrt(v_safe^2 + 2 a_brk h)
  double v_safe2{0};
  double a_brk{0};
  template <class T>
  T operator()(const std::array<T, NXA>& X) const {
    using std::sin;
    using std::sqrt;
    return X[LV] * sin(X[LGAM]) + sqrt((2.0 * a_brk) * X[LH] + v_safe2);
  }
};

struct AirspeedBarrier {  // b_V = V - V_min
  double Vmin{0};
  template <class T>
  T operator()(const std::array<T, NXA>& X) const {
    return X[LV] - Vmin;
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
