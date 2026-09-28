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
Taylor<N, S> operator/(const Taylor<N, S>& a, double b) {
  Taylor<N, S> r; for (int k = 0; k <= N; ++k) r[k] = a[k] / b; return r;
}
template <int N, class S>
Taylor<N, S> operator/(double a, const Taylor<N, S>& b) { return Taylor<N, S>(S(a)) / b; }
// Derivative series a'(t) (top coefficient truncated to 0), for the
// integral-form rules below: y = int q dt  <=>  y[k] = q[k-1]/k.
template <int N, class S>
Taylor<N, S> taylorDeriv(const Taylor<N, S>& a) {
  Taylor<N, S> d; for (int k = 0; k < N; ++k) d[k] = double(k + 1) * a[k + 1]; return d;
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
// e(t) = exp(a(t)), from e' = a' e:  k e[k] = sum_{j=1}^k j a[j] e[k-j].
template <int N, class S>
Taylor<N, S> exp(const Taylor<N, S>& a) {
  using std::exp;
  Taylor<N, S> r; r[0] = exp(a[0]);
  for (int k = 1; k <= N; ++k) {
    S s = S(0.0);
    for (int j = 1; j <= k; ++j) s += double(j) * a[j] * r[k - j];
    r[k] = s / double(k);
  }
  return r;
}
// t(t) = tanh(a(t)), from t' = (1 - t^2) a'. w = 1 - t*t is built progressively:
// w[i] needs t[0..i], available before t[k] for every i <= k-1 in the k t[k]
// = sum_{j=1}^k j a[j] w[k-j] recursion.
template <int N, class S>
Taylor<N, S> tanh(const Taylor<N, S>& a) {
  using std::tanh;
  Taylor<N, S> t, w;
  t[0] = tanh(a[0]);
  w[0] = 1.0 - t[0] * t[0];
  for (int k = 1; k <= N; ++k) {
    S s = S(0.0);
    for (int j = 1; j <= k; ++j) s += double(j) * a[j] * w[k - j];
    t[k] = s / double(k);
    S ww = S(0.0);
    for (int l = 0; l <= k; ++l) ww += t[l] * t[k - l];
    w[k] = -ww;
  }
  return t;
}
// tan = sin / cos (both exact above).
template <int N, class S>
Taylor<N, S> tan(const Taylor<N, S>& a) { return sin(a) / cos(a); }
// atan(a): y' = a' / (1 + a^2); asin(a): y' = a' / sqrt(1 - a^2). Integral
// form: y[k] = q[k-1]/k, y[0] from the scalar function. As with atan2 the
// top coefficient of q is never read.
template <int N, class S>
Taylor<N, S> atan(const Taylor<N, S>& a) {
  using std::atan;
  const Taylor<N, S> q = taylorDeriv(a) / (1.0 + a * a);
  Taylor<N, S> y; y[0] = atan(a[0]);
  for (int k = 1; k <= N; ++k) y[k] = q[k - 1] / double(k);
  return y;
}
template <int N, class S>
Taylor<N, S> asin(const Taylor<N, S>& a) {
  using std::asin;
  const Taylor<N, S> q = taylorDeriv(a) / sqrt(1.0 - a * a);
  Taylor<N, S> y; y[0] = asin(a[0]);
  for (int k = 1; k <= N; ++k) y[k] = q[k - 1] / double(k);
  return y;
}
// theta(t) = atan2(y(t), x(t)), from theta' = (x y' - y x')/(x^2 + y^2):
//   theta[k] = q[k-1]/k with q = (x y' - y x')/(x^2 + y^2).
// The [N] coefficient of the numerator would need x[N+1]/y[N+1]; it only feeds
// q[N], which is never read, so it is left at its (wrong) truncated value.
template <int N, class S>
Taylor<N, S> atan2(const Taylor<N, S>& y, const Taylor<N, S>& x) {
  using std::atan2;
  Taylor<N, S> xd, yd;  // derivative series (top coefficient truncated to 0)
  for (int k = 0; k < N; ++k) { xd[k] = double(k + 1) * x[k + 1]; yd[k] = double(k + 1) * y[k + 1]; }
  const Taylor<N, S> num = x * yd - y * xd;
  const Taylor<N, S> den = x * x + y * y;
  const Taylor<N, S> q = num / den;
  Taylor<N, S> th;
  th[0] = atan2(y[0], x[0]);
  for (int k = 1; k <= N; ++k) th[k] = q[k - 1] / double(k);
  return th;
}

// Value part of a (possibly nested) jet scalar -- for BRANCH DECISIONS inside
// templated barrier/drift code (e.g. a series-vs-exact switch near a removable
// singularity). Branching on the value keeps the jet exact away from the seam.
inline double scalarValue(double x) { return x; }
inline double scalarValue(const autodiff::dual& x) {
  return autodiff::detail::val(x);
}
template <int N, class S>
double scalarValue(const Taylor<N, S>& x) { return scalarValue(x.c[0]); }

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
