#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <array>
#include <cmath>
#include "autoland/lie_taylor.hpp"
#include "autoland/water_rollout.hpp"
#include "autoland/world_geometry.hpp"

using namespace autoland;
using Catch::Approx;

namespace {
constexpr double kDeg = M_PI / 180.0;

// Hand chain-rule gradient of one bump (independent of TerrainField::height).
Eigen::Vector2d handBumpGrad(const TerrainBump& b, double N, double E) {
  const double c = std::cos(b.rot), s = std::sin(b.rot);
  const double dN = N - b.N0, dE = E - b.E0;
  const double xa = (dN * c + dE * s) / b.sigma_a;
  const double xb = (dE * c - dN * s) / b.sigma_b;
  const double q = xa * xa + xb * xb;
  const double qp = std::pow(q, b.order);
  const double h = b.height * std::exp(-qp);
  // dh/dq = -h * order * q^(order-1)
  const double dh_dq = -h * b.order * std::pow(q, b.order - 1);
  // dq/dN = 2 xa (c/sa) + 2 xb (-s/sb);  dq/dE = 2 xa (s/sa) + 2 xb (c/sb)
  const double dq_dN = 2.0 * xa * c / b.sigma_a - 2.0 * xb * s / b.sigma_b;
  const double dq_dE = 2.0 * xa * s / b.sigma_a + 2.0 * xb * c / b.sigma_b;
  return {dh_dq * dq_dN, dh_dq * dq_dE};
}
}  // namespace

TEST_CASE("Corridor frame: hand values, round trip, containment", "[world]") {
  Corridor c;
  c.N_c = 100.0; c.E_c = -50.0; c.length = 600.0; c.width = 80.0;
  c.heading = 30.0 * kDeg; c.s_aim = -200.0;

  // A point 100 m ahead along the heading and 10 m to the right of it.
  const double N = c.N_c + 100.0 * std::cos(c.heading) + 10.0 * std::cos(c.heading + M_PI / 2);
  const double E = c.E_c + 100.0 * std::sin(c.heading) + 10.0 * std::sin(c.heading + M_PI / 2);
  double s, e;
  c.frame(N, E, s, e);
  CHECK(s == Approx(100.0).margin(1e-12));
  CHECK(e == Approx(10.0).margin(1e-12));
  double N2, E2;
  c.inverse(s, e, N2, E2);
  CHECK(N2 == Approx(N).margin(1e-12));
  CHECK(E2 == Approx(E).margin(1e-12));

  // Corners: just inside vs just outside.
  const double eps = 1e-6;
  for (double ss : {c.sNear(), c.sFar()})
    for (double ee : {-0.5 * c.width, 0.5 * c.width}) {
      const double sin_ = ss - std::copysign(eps, ss), ein = ee - std::copysign(eps, ee);
      const double sout = ss + std::copysign(eps, ss), eout = ee + std::copysign(eps, ee);
      double Ni, Ei, No, Eo;
      c.inverse(sin_, ein, Ni, Ei);
      c.inverse(sout, eout, No, Eo);
      CHECK(c.contains(Ni, Ei));
      CHECK_FALSE(c.contains(No, Eo));
    }
  const Eigen::Vector2d aim = c.aimPoint();
  double sa, ea;
  c.frame(aim[0], aim[1], sa, ea);
  CHECK(sa == Approx(c.s_aim).margin(1e-12));
  CHECK(ea == Approx(0.0).margin(1e-12));
  CHECK(c.contains(aim[0], aim[1]));
}

TEST_CASE("Terrain field: heights, exact gradient vs hand chain rule", "[world]") {
  TerrainField tf;
  tf.base = 0.0;
  TerrainBump mtn{"mtn", 3000.0, 900.0, 250.0, 400.0, 250.0, 30.0 * kDeg, 1};
  TerrainBump house{"house", -800.0, 250.0, 8.0, 6.0, 6.0, 0.0, 2};
  tf.bumps = {mtn, house};

  CHECK(tf.height(3000.0, 900.0) == Approx(250.0).margin(1e-9));
  CHECK(tf.height(-800.0, 250.0) == Approx(8.0).margin(1e-9));
  CHECK(tf.height(0.0, 0.0) < 1e-6);                        // lake center
  CHECK(tf.height(3000.0 + 400.0 * std::cos(mtn.rot), 900.0 + 400.0 * std::sin(mtn.rot))
        == Approx(250.0 * std::exp(-1.0)).margin(1e-9));    // one sigma along a
  // Flat-top: at 0.5 sigma the quartic is much fuller than a Gaussian.
  CHECK(tf.height(-800.0 + 3.0, 250.0) == Approx(8.0 * std::exp(-0.0625)).margin(1e-9));

  for (auto pt : {std::pair{3100.0, 800.0}, std::pair{2700.0, 1100.0},
                  std::pair{-797.0, 253.0}, std::pair{-805.0, 246.0}}) {
    const Eigen::Vector2d g = tf.gradient(pt.first, pt.second);
    const Eigen::Vector2d hand =
        handBumpGrad(mtn, pt.first, pt.second) + handBumpGrad(house, pt.first, pt.second);
    CHECK(g[0] == Approx(hand[0]).margin(1e-12));
    CHECK(g[1] == Approx(hand[1]).margin(1e-12));
  }
}

TEST_CASE("Terrain field evaluates exactly under the Taylor jet", "[world]") {
  TerrainField tf;
  tf.bumps = {TerrainBump{"g", 200.0, -100.0, 40.0, 80.0, 50.0, 20.0 * kDeg, 1},
              TerrainBump{"q", -50.0, 60.0, 10.0, 30.0, 30.0, 0.0, 2}};
  // Constant-velocity flow on (N, E): L_f h = grad h . v, L_f^2 h = v' H v.
  const std::array<double, 2> v{31.0, -12.0};
  auto f = [&](const auto& X) {
    using T = std::decay_t<decltype(X[0])>;
    return std::array<T, 2>{T(v[0]), T(v[1])};
  };
  auto b = [&](const auto& X) { return tf.height(X[0], X[1]); };
  const std::array<double, 2> x0{150.0, -60.0};
  const auto L = lieDrift<2, 2>(f, b, x0);
  const Eigen::Vector2d g = tf.gradient(x0[0], x0[1]);
  CHECK(L[0] == Approx(tf.height(x0[0], x0[1])).margin(1e-12));
  CHECK(L[1] == Approx(g[0] * v[0] + g[1] * v[1]).margin(1e-12));
  // Second order: L_f^2 h = v' H v. Oracle = central difference of the EXACT
  // gradient along v (test-side oracle only, cf. test_lie_derivatives).
  const double hh = 1e-3;
  const Eigen::Vector2d gp = tf.gradient(x0[0] + hh * v[0], x0[1] + hh * v[1]);
  const Eigen::Vector2d gm = tf.gradient(x0[0] - hh * v[0], x0[1] - hh * v[1]);
  const double vHv = ((gp - gm) / (2.0 * hh)).dot(Eigen::Vector2d(v[0], v[1]));
  CHECK(L[2] == Approx(vHv).epsilon(1e-6));
}

TEST_CASE("Terrain gradient bound and row margin", "[world]") {
  TerrainField tf;
  tf.bumps = {TerrainBump{"g", 0.0, 0.0, 100.0, 200.0, 300.0, 0.3, 1},
              TerrainBump{"q", 900.0, 0.0, 15.0, 25.0, 25.0, 0.0, 2}};
  const double gmax = terrainMaxGradient(tf);
  // Sampled gradient never exceeds the analytic bound.
  double worst = 0.0;
  for (double N = -1200; N <= 1200; N += 7.0)
    for (double E = -800; E <= 800; E += 7.0) worst = std::max(worst, tf.gradient(N, E).norm());
  CHECK(worst <= gmax);
  CHECK(worst > 0.3 * gmax);  // and the bound is not absurdly loose
  const double m = terrainRowMargin(5.0, 20.0, gmax, 14.63, 30.0 * kDeg);
  CHECK(m == Approx(5.0 + gmax * 20.0 + 0.5 * 14.63 * std::sin(30.0 * kDeg)).margin(1e-12));
}

TEST_CASE("Terrain ring builder", "[world]") {
  TerrainField tf;
  tf.addRing(0.0, 0.0, 600.0, 48, 15.0, 25.0, 2, "trees");
  REQUIRE(tf.bumps.size() == 48);
  for (const auto& b : tf.bumps) {
    CHECK(std::hypot(b.N0, b.E0) == Approx(600.0).margin(1e-9));
    CHECK(b.order == 2);
  }
  CHECK(tf.height(600.0, 0.0) == Approx(15.0).margin(1e-3));
  CHECK(tf.height(0.0, 0.0) < 1e-6);
  // Gap on the approach heading (both ends of the axis).
  TerrainField gap;
  gap.addRing(0.0, 0.0, 600.0, 48, 15.0, 25.0, 2, "trees", 0.0, 10.0 * kDeg);
  CHECK(gap.bumps.size() < 48);
  CHECK(gap.height(600.0, 0.0) < 1e-3);
  CHECK(gap.height(-600.0, 0.0) < 1e-3);
}

TEST_CASE("Water rollout closed forms and corridor containment", "[world][rollout]") {
  Corridor cor;
  cor.length = 700.0; cor.width = 80.0; cor.heading = 0.0; cor.s_aim = -250.0;

  SECTION("constant deceleration") {
    RolloutConfig c;
    c.a0 = 0.8; c.kq = 0.0; c.k_psi = 0.0; c.V_stop = 0.0; c.t_max = 200.0;
    RolloutState s;
    s.N = cor.sNear() + 50.0; s.E = 0.0; s.chi = 0.0; s.Vg = 30.0;
    const RolloutResult r = runRollout(s, 0.0, 0.01, c, cor);
    CHECK(r.stopped);
    CHECK(r.distance == Approx(30.0 * 30.0 / (2.0 * 0.8)).epsilon(1e-6));
    CHECK(r.t_stop == Approx(30.0 / 0.8).margin(0.011));
  }
  SECTION("quadratic drag") {
    RolloutConfig c;
    c.a0 = 0.0; c.kq = 0.002; c.k_psi = 0.0; c.V_stop = 0.0; c.t_max = 20.0;
    RolloutState s;
    s.N = 0.0; s.E = 0.0; s.chi = 0.0; s.Vg = 30.0;
    const RolloutResult r = runRollout(s, 0.0, 0.01, c, cor);
    CHECK_FALSE(r.stopped);  // asymptotic
    const double T = r.t_stop;
    CHECK(r.distance == Approx(std::log(1.0 + 0.002 * 30.0 * T) / 0.002).epsilon(1e-6));
  }
  SECTION("heading steers to the corridor heading") {
    Corridor rot = cor;
    rot.heading = 40.0 * kDeg;
    RolloutConfig c;
    c.a0 = 0.3; c.kq = 0.0; c.k_psi = 0.5; c.V_stop = 0.5; c.t_max = 200.0;
    RolloutState s;
    const Eigen::Vector2d p0 = rot.axisPoint(rot.s_aim);
    s.N = p0[0]; s.E = p0[1]; s.chi = rot.heading + 8.0 * kDeg; s.Vg = 25.0;
    const RolloutResult r = runRollout(s, 0.0, 0.01, c, rot);
    CHECK(std::abs(wrapAngle(r.chi_stop - rot.heading)) < 0.1 * kDeg);
  }
  SECTION("touchdown at the aim point stops inside; hot and long stops outside") {
    RolloutConfig c;  // defaults a0 0.8, kq 0.002
    RolloutState s;
    s.N = cor.s_aim; s.E = 3.0; s.chi = 2.0 * kDeg; s.Vg = 33.0;
    const RolloutResult ok = runRollout(s, 10.0, 0.01, c, cor);
    CHECK(ok.stopped);
    CHECK(ok.inside);
    CHECK(ok.stayed_inside);
    CHECK(ok.s_stop < cor.sFar());
    // Combined a0, kq closed form: d_stop = ln((a0 + kq V0^2)/a0) / (2 kq).
    const double d_ok = std::log((c.a0 + c.kq * 33.0 * 33.0) / c.a0) / (2.0 * c.kq);
    CHECK(ok.distance == Approx(d_ok).epsilon(2e-3));  // V_stop = 0.5 truncates
    CHECK(ok.s_stop == Approx(cor.s_aim + d_ok).margin(3.0));
    RolloutState hot;
    hot.N = cor.sFar() - 0.33 * cor.length; hot.E = 0.0; hot.chi = 0.0; hot.Vg = 45.0;
    const RolloutResult bad = runRollout(hot, 10.0, 0.01, c, cor);
    const double d_hot = std::log((c.a0 + c.kq * 45.0 * 45.0) / c.a0) / (2.0 * c.kq);
    CHECK(hot.N + d_hot > cor.sFar());  // the closed form predicts the overrun
    CHECK_FALSE(bad.inside);
    CHECK_FALSE(bad.stayed_inside);
  }
}
