#pragma once
#include <array>
#include <cmath>
#include <autodiff/forward/dual.hpp>

// =============================================================================
// Exact Lie derivatives along a (varying) vector field via the flow Taylor jet.
//
// The HOCBF constraints need L_f^k b and the control row L_g L_f^{r-1} b for a
// control-affine system Xdot = f(X) + g(X) U. These are computed EXACTLY (no
// finite differences): we build the order-r Taylor jet of the flow X(t) (Xdot =
// f(X)) by Picard iteration in truncated-Taylor arithmetic, then read
//     L_f^k b = k! * [t^k] b(X(t)).
// The control row is the directional derivative of the field L_f^{r-1} b along a
// g-column, obtained by running the same jet with the coefficient scalar
// promoted to autodiff::dual and seeding the base point along that column.
//
// Why a hand-rolled Taylor type rather than autodiff::real: autodiff's `real`
// requires an arithmetic coefficient type, so it cannot be nested over `dual`
// for the control-row directional derivative. Taylor<N,S> is templated on the
// coefficient scalar S (double for the drift, autodiff::dual for the control
// row), so the same jet code serves both. autodiff::real's `along()` would in
// any case give fixed-direction derivatives (v.grad)^k b, NOT L_f^k b, which
// differ from order 2 because f varies along the flow.
// =============================================================================
namespace autoland {

// Truncated Taylor series of order N with coefficient scalar S (the series
// variable is the flow time t). c[k] is the t^k coefficient.
template <int N, class S>
struct Taylor {
  std::array<S, N + 1> c{};
  Taylor() { for (auto& v : c) v = S(0.0); }
  Taylor(const S& v) { for (auto& x : c) x = S(0.0); c[0] = v; }  // NOLINT: implicit
  S& operator[](int i) { return c[i]; }
  const S& operator[](int i) const { return c[i]; }
};

template <int N, class S>
Taylor<N, S> operator+(const Taylor<N, S>& a, const Taylor<N, S>& b) {
  Taylor<N, S> r; for (int k = 0; k <= N; ++k) r[k] = a[k] + b[k]; return r;
}
template <int N, class S>
Taylor<N, S> operator-(const Taylor<N, S>& a, const Taylor<N, S>& b) {
  Taylor<N, S> r; for (int k = 0; k <= N; ++k) r[k] = a[k] - b[k]; return r;
}
template <int N, class S>
Taylor<N, S> operator-(const Taylor<N, S>& a) {
  Taylor<N, S> r; for (int k = 0; k <= N; ++k) r[k] = -a[k]; return r;
}
template <int N, class S>
Taylor<N, S> operator*(const Taylor<N, S>& a, const Taylor<N, S>& b) {
  Taylor<N, S> r;
  for (int k = 0; k <= N; ++k) { S s = S(0.0); for (int j = 0; j <= k; ++j) s += a[j] * b[k - j]; r[k] = s; }
  return r;
}
template <int N, class S>
Taylor<N, S> operator*(double a, const Taylor<N, S>& b) {
  Taylor<N, S> r; for (int k = 0; k <= N; ++k) r[k] = a * b[k]; return r;
}
template <int N, class S>
Taylor<N, S> operator*(const Taylor<N, S>& a, double b) { return b * a; }
template <int N, class S>
Taylor<N, S> operator+(const Taylor<N, S>& a, double b) {
  Taylor<N, S> r = a; r[0] = a[0] + b; return r;
}
template <int N, class S>
Taylor<N, S> operator+(double b, const Taylor<N, S>& a) { return a + b; }
template <int N, class S>
Taylor<N, S> operator-(const Taylor<N, S>& a, double b) {
  Taylor<N, S> r = a; r[0] = a[0] - b; return r;
}
template <int N, class S>
Taylor<N, S> operator-(double b, const Taylor<N, S>& a) {
  Taylor<N, S> r = -a; r[0] = b + r[0]; return r;
}
template <int N, class S>
Taylor<N, S> operator/(const Taylor<N, S>& a, const Taylor<N, S>& b) {
  Taylor<N, S> r;
  for (int k = 0; k <= N; ++k) { S s = a[k]; for (int j = 1; j <= k; ++j) s -= b[j] * r[k - j]; r[k] = s / b[0]; }
  return r;
}
template <int N, class S>
Taylor<N, S> sin(const Taylor<N, S>& a) {
  using std::sin; using std::cos;
  Taylor<N, S> s, co; s[0] = sin(a[0]); co[0] = cos(a[0]);
  for (int k = 1; k <= N; ++k) {
    S ds = S(0.0), dc = S(0.0);
    for (int j = 1; j <= k; ++j) { ds += double(j) * a[j] * co[k - j]; dc += double(j) * a[j] * s[k - j]; }
    s[k] = ds / double(k); co[k] = (-1.0) * dc / double(k);
  }
  return s;
}
template <int N, class S>
Taylor<N, S> cos(const Taylor<N, S>& a) {
  using std::sin; using std::cos;
  Taylor<N, S> s, co; s[0] = sin(a[0]); co[0] = cos(a[0]);
  for (int k = 1; k <= N; ++k) {
    S ds = S(0.0), dc = S(0.0);
    for (int j = 1; j <= k; ++j) { ds += double(j) * a[j] * co[k - j]; dc += double(j) * a[j] * s[k - j]; }
    s[k] = ds / double(k); co[k] = (-1.0) * dc / double(k);
  }
  return co;
}
template <int N, class S>
Taylor<N, S> sqrt(const Taylor<N, S>& a) {
  using std::sqrt;
  Taylor<N, S> r; r[0] = sqrt(a[0]);
  for (int k = 1; k <= N; ++k) {
    S s = a[k]; for (int j = 1; j < k; ++j) s -= r[j] * r[k - j];
    r[k] = s / (2.0 * r[0]);
  }
  return r;
}

// Build the order-R flow jet of X(t) (Xdot = f(X), X(0)=x0) and return
// {b, L_f b, ..., L_f^R b} as coefficient-typed scalars. `f` and `b` are
// callables templated on element type:
//   std::array<T,NX> f(const std::array<T,NX>&);   T = Taylor<R,S>
//   T               b(const std::array<T,NX>&);
template <int R, int NX, class S, class FDrift, class Barrier>
std::array<S, R + 1> lieJet(const FDrift& f, const Barrier& b,
                            const std::array<S, NX>& x0) {
  using TS = Taylor<R, S>;
  std::array<TS, NX> X;
  for (int i = 0; i < NX; ++i) { X[i] = TS(); X[i][0] = x0[i]; }
  // Picard: x_i coeff (k+1) = f_i coeff (k) / (k+1). R passes fix orders 1..R.
  for (int pass = 0; pass < R; ++pass) {
    std::array<TS, NX> F = f(X);
    for (int i = 0; i < NX; ++i)
      for (int k = 0; k < R; ++k) X[i][k + 1] = F[i][k] / double(k + 1);
  }
  TS B = b(X);
  std::array<S, R + 1> L; double fact = 1.0;
  for (int k = 0; k <= R; ++k) { L[k] = fact * B[k]; fact *= double(k + 1); }
  return L;
}

// Drift Lie derivatives {b, L_f b, ..., L_f^R b} as doubles.
template <int R, int NX, class FDrift, class Barrier>
std::array<double, R + 1> lieDrift(const FDrift& f, const Barrier& b,
                                   const std::array<double, NX>& x0) {
  return lieJet<R, NX, double>(f, b, x0);
}

// Control-row entry: L_g L_f^{R-1} b = d/ds[ L_f^{R-1} b (x0 + s*dir) ]|_0,
// with `dir` a g-column (evaluated at x0). Exact via autodiff::dual seeding.
template <int R, int NX, class FDrift, class Barrier>
double lieAlong(const FDrift& f, const Barrier& b,
                const std::array<double, NX>& x0,
                const std::array<double, NX>& dir) {
  using autodiff::dual;
  std::array<dual, NX> xs;
  for (int i = 0; i < NX; ++i) { xs[i] = x0[i]; xs[i].grad = dir[i]; }
  std::array<dual, R + 1> L = lieJet<R, NX, dual>(f, b, xs);
  return autodiff::detail::derivative(L[R - 1]);
}

}  // namespace autoland
