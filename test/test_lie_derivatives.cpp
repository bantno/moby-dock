#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <vector>

#include "autoland/hocbf.hpp"       // hocbfRow + NUA
#include "autoland/lie_taylor.hpp"  // lieDrift / lieAlong

// =============================================================================
// Closed-form verification of the Lie-derivative engine (lie_taylor.hpp) and the
// HOCBF row assembly (hocbf.hpp::hocbfRow), on textbook systems whose L_f^k b,
// L_g L_f^{r-1} b, and HOCBF constraint rows are computable by hand. This is an
// ANALYTIC cross-check (complements the finite-difference flow oracle in
// test_lon_cbf.cpp, which checks lieDrift on the real longitudinal field):
//   - exact L_f^k b incl. a state-varying field (pendulum) and a high-order jet,
//   - exact, possibly state-dependent control row L_g L_f^{r-1} b,
//   - linear class-K coefficients landing on the right Lie terms:
//       RD2:  rhs = L_f^2 b + (c1+c2) L_f b + c1 c2 b,         a = -L_g L_f b
//       RD3:  rhs = L_f^3 b + (Sum c) L_f^2 b + (Sum cc) L_f b + (prod c) b
// =============================================================================
using namespace autoland;
using Catch::Approx;

namespace {

// ---- Textbook drift / barrier functors (templated on element type) ----------
struct DoubleInt {  // x=[p,v], p'=v, v'=u  -> f=[v,0], g=[0,1]
  template <class T>
  std::array<T, 2> operator()(const std::array<T, 2>& x) const {
    return {x[1], T(0.0)};
  }
};
struct TripleInt {  // x=[p,v,a], f=[v,a,0], g=[0,0,1]
  template <class T>
  std::array<T, 3> operator()(const std::array<T, 3>& x) const {
    return {x[1], x[2], T(0.0)};
  }
};
struct Pendulum {  // x=[th,w], th'=w, w'=-sin th (+u) -> f=[w,-sin th], g=[0,1]
  template <class T>
  std::array<T, 2> operator()(const std::array<T, 2>& x) const {
    using std::sin;
    return {x[1], -sin(x[0])};
  }
};
struct ScalarSq {  // x'=x^2, b=x  -> flow x0/(1-x0 t), L_f^k b = k! x^{k+1}
  template <class T>
  std::array<T, 1> operator()(const std::array<T, 1>& x) const {
    return {x[0] * x[0]};
  }
};
struct B_p  { template <class T> T operator()(const std::array<T, 2>& x) const { return x[0]; } };
struct B_p2 { template <class T> T operator()(const std::array<T, 2>& x) const { return x[0] * x[0]; } };
struct B_p3 { template <class T> T operator()(const std::array<T, 3>& x) const { return x[0]; } };
struct B_x  { template <class T> T operator()(const std::array<T, 1>& x) const { return x[0]; } };

// Assemble a hocbfRow from a SCALAR control coefficient (single input, column 0).
HocbfRow row1(const std::vector<double>& Lf, double LgLf_scalar,
              const std::vector<double>& c) {
  Eigen::Matrix<double, 1, NUA> LgLf;
  LgLf.setZero();
  LgLf(0, 0) = LgLf_scalar;
  return hocbfRow(Lf, LgLf, c);
}

constexpr double kTol = 1e-9;

}  // namespace

TEST_CASE("Lie engine: double integrator, b=p (RD2)") {
  DoubleInt f;
  B_p b;
  const std::array<double, 2> x0{1.3, -0.7};
  const auto L = lieDrift<2, 2>(f, b, x0);
  const double LgLf = lieAlong<2, 2>(f, b, x0, {0.0, 1.0});

  CHECK(L[0] == Approx(1.3).margin(kTol));   // b   = p
  CHECK(L[1] == Approx(-0.7).margin(kTol));  // L_f b   = v
  CHECK(L[2] == Approx(0.0).margin(kTol));   // L_f^2 b = 0
  CHECK(LgLf == Approx(1.0).margin(kTol));   // L_g L_f b = 1

  const HocbfRow r = row1({L[0], L[1], L[2]}, LgLf, {1.5, 2.5});
  CHECK(r.a(0, 0) == Approx(-1.0).margin(kTol));  // a = -L_g L_f b
  // rhs = L_f^2 b + (c1+c2) L_f b + c1 c2 b
  CHECK(r.rhs == Approx((1.5 + 2.5) * (-0.7) + (1.5 * 2.5) * 1.3).margin(kTol));
}

TEST_CASE("Lie engine: double integrator, b=p^2 (RD2, state-dependent control row)") {
  DoubleInt f;
  B_p2 b;
  const std::array<double, 2> x0{1.3, -0.7};
  const auto L = lieDrift<2, 2>(f, b, x0);
  const double LgLf = lieAlong<2, 2>(f, b, x0, {0.0, 1.0});

  CHECK(L[0] == Approx(1.3 * 1.3).margin(kTol));      // p^2
  CHECK(L[1] == Approx(2 * 1.3 * -0.7).margin(kTol));  // 2pv
  CHECK(L[2] == Approx(2 * 0.7 * 0.7).margin(kTol));   // 2v^2
  CHECK(LgLf == Approx(2 * 1.3).margin(kTol));         // L_g L_f b = 2p (state-dependent)

  const HocbfRow r = row1({L[0], L[1], L[2]}, LgLf, {1.5, 2.5});
  CHECK(r.a(0, 0) == Approx(-2 * 1.3).margin(kTol));
  CHECK(r.rhs == Approx(2 * 0.49 + 4.0 * (2 * 1.3 * -0.7) + 3.75 * (1.3 * 1.3)).margin(kTol));
}

TEST_CASE("Lie engine: triple integrator, b=p (RD3 coefficient placement)") {
  TripleInt f;
  B_p3 b;
  const std::array<double, 3> x0{0.4, -0.2, 0.9};
  const auto L = lieDrift<3, 3>(f, b, x0);
  const double LgLf2 = lieAlong<3, 3>(f, b, x0, {0.0, 0.0, 1.0});

  CHECK(L[1] == Approx(-0.2).margin(kTol));   // v
  CHECK(L[2] == Approx(0.9).margin(kTol));    // a
  CHECK(L[3] == Approx(0.0).margin(kTol));    // 0
  CHECK(LgLf2 == Approx(1.0).margin(kTol));   // L_g L_f^2 b = 1

  // c=[1,2,3]: Sum c=6, Sum cc=11, prod=6 -> rhs = 6a + 11v + 6p
  const HocbfRow r = row1({L[0], L[1], L[2], L[3]}, LgLf2, {1.0, 2.0, 3.0});
  CHECK(r.a(0, 0) == Approx(-1.0).margin(kTol));
  CHECK(r.rhs == Approx(6 * 0.9 + 11 * -0.2 + 6 * 0.4).margin(kTol));
}

TEST_CASE("Lie engine: pendulum, b=theta (RD3, state-varying field, sin/cos)") {
  Pendulum f;
  B_p b;  // b = x[0] = theta
  const double th = 0.6, w = -0.3;
  const std::array<double, 2> x0{th, w};
  const auto L = lieDrift<3, 2>(f, b, x0);
  const double LgLf = lieAlong<2, 2>(f, b, x0, {0.0, 1.0});

  CHECK(L[1] == Approx(w).margin(kTol));                  // L_f b   = w
  CHECK(L[2] == Approx(-std::sin(th)).margin(kTol));      // L_f^2 b = -sin th
  CHECK(L[3] == Approx(-w * std::cos(th)).margin(kTol));  // L_f^3 b = -w cos th
  CHECK(LgLf == Approx(1.0).margin(kTol));                // L_g L_f b = 1

  const HocbfRow r = row1({L[0], L[1], L[2], L[3]}, LgLf, {1.0, 2.0, 3.0});
  CHECK(r.rhs ==
        Approx((-w * std::cos(th)) + 6 * (-std::sin(th)) + 11 * w + 6 * th).margin(kTol));
}

TEST_CASE("Lie engine: scalar x'=x^2, b=x (high-order jet, L_f^k b = k! x^{k+1})") {
  ScalarSq f;
  B_x b;
  const double x = 0.5;
  const std::array<double, 1> x0{x};
  const auto L = lieDrift<4, 1>(f, b, x0);

  double fact = 1.0, xp = x;
  for (int k = 0; k <= 4; ++k) {
    CHECK(L[k] == Approx(fact * xp).margin(kTol));
    fact *= (k + 1);
    xp *= x;
  }
}
