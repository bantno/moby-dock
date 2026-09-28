#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <random>
#include "autoland/lie_taylor.hpp"
#include "autoland/sixdof_cbf.hpp"

using namespace autoland;
using Catch::Approx;

namespace {
constexpr double kDeg = M_PI / 180.0;

std::array<double, NX> arr(const StateVec& x) {
  std::array<double, NX> a;
  for (int i = 0; i < NX; ++i) a[i] = x[i];
  return a;
}
StateVec randomState(std::mt19937& rng) {
  std::uniform_real_distribution<double> un(-1.0, 1.0);
  StateVec x = StateVec::Zero();
  // NOTE: the local must NOT be named V -- that shadows the State::V index.
  const double Vt = 42.0 + 8.0 * un(rng), al = (6.0 + 5.0 * un(rng)) * kDeg,
               be = 6.0 * kDeg * un(rng);
  x[U] = Vt * std::cos(al) * std::cos(be);
  x[V] = Vt * std::sin(be);
  x[W] = Vt * std::sin(al) * std::cos(be);
  x[P] = 0.4 * un(rng); x[Q] = 0.2 * un(rng); x[R] = 0.2 * un(rng);
  x[PHI] = 25.0 * kDeg * un(rng); x[THETA] = 10.0 * kDeg * un(rng);
  x[PSI] = M_PI * un(rng); x[H] = 100.0 + 50.0 * un(rng);
  x[Y] = 100.0 * un(rng); x[XN] = 500.0 * un(rng);
  return x;
}
SixDofCbfConfig baseCfg() {
  SixDofCbfConfig c;
  c.limits.de_min = -25 * kDeg; c.limits.de_max = 25 * kDeg;
  c.limits.da_min = -25 * kDeg; c.limits.da_max = 25 * kDeg;
  c.limits.dr_min = -25 * kDeg; c.limits.dr_max = 25 * kDeg;
  c.limits.rate = 150.0 * kDeg;
  c.dt = 0.01;
  return c;
}
// RK4 on the full plant with ZOH surfaces, throttle fixed.
StateVec rk4(const BeaverDynamics& dyn, const StateVec& x, const CtrlVec& u, double dt) {
  const StateVec k1 = dyn.xdot(x, u);
  const StateVec k2 = dyn.xdot(x + 0.5 * dt * k1, u);
  const StateVec k3 = dyn.xdot(x + 0.5 * dt * k2, u);
  const StateVec k4 = dyn.xdot(x + dt * k3, u);
  return x + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
}
}  // namespace

TEST_CASE("6-DOF CBF model: f + g u reproduces the plant with frozen forces (with wind)",
          "[sixdof_cbf]") {
  BeaverDynamics dyn;
  std::mt19937 rng(21);
  for (int k = 0; k < 20; ++k) {
    const StateVec x = randomState(rng);
    std::uniform_real_distribution<double> un(-1.0, 1.0);
    const SurfVec up(0.1 * un(rng), 0.1 * un(rng), 0.1 * un(rng));
    const double dT = 0.5 + 0.4 * un(rng);
    const Eigen::Vector3d We(6.0 * un(rng), 6.0 * un(rng), 2.0 * un(rng));
    const SixDofCbfModel f(dyn, up, dT, We);
    const auto x0 = arr(x);
    const auto fx = f(x0);
    const auto g = f.gColumns(x0);
    const auto F = f.forceColumns(x0);
    CtrlVec u;
    u[DE] = up[0]; u[DA] = up[1]; u[DR] = up[2]; u[DT] = dT;
    const StateVec plant = dyn.xdot(x, u, We);
    // At U = u_prev the split is exact in every row.
    for (int i = 0; i < NX; ++i) {
      double s = fx[i];
      for (int kk = 0; kk < NUS; ++kk) s += g[kk][i] * up[kk];
      CHECK(s == Approx(plant[i]).margin(1e-9));
    }
    // At U != u_prev: moment rows exact, force rows differ by F (U - u_prev).
    const SurfVec U(up[0] + 0.05, up[1] - 0.04, up[2] + 0.03);
    CtrlVec u2 = u; u2[DE] = U[0]; u2[DA] = U[1]; u2[DR] = U[2];
    const StateVec plant2 = dyn.xdot(x, u2, We);
    for (int i = 0; i < NX; ++i) {
      double s = fx[i];
      for (int kk = 0; kk < NUS; ++kk) s += g[kk][i] * U[kk] + F[kk][i] * (U[kk] - up[kk]);
      CHECK(s == Approx(plant2[i]).margin(1e-9));
    }
  }
}

TEST_CASE("6-DOF CBF model: Taylor jets equal the plant derivative and its linearization",
          "[sixdof_cbf]") {
  BeaverDynamics dyn;
  std::mt19937 rng(22);
  for (int k = 0; k < 10; ++k) {
    const StateVec x = randomState(rng);
    const double dT = 0.6;
    // u_prev = 0: f(x) = plant(x, [0,0,0,dT]) exactly, so L_f^2 x_i = (A f)_i.
    const SixDofCbfModel f(dyn, SurfVec::Zero(), dT);
    const auto x0 = arr(x);
    CtrlVec u = CtrlVec::Zero(); u[DT] = dT;
    const StateVec fx = dyn.xdot(x, u);
    Mat A, B;
    dyn.linearize(x, u, A, B);
    const StateVec Af = A * fx;
    for (int i = 0; i < NX; ++i) {
      auto b = [i](const auto& X) { return X[i]; };
      const auto L = lieDrift<2, NX>(f, b, x0);
      CHECK(L[0] == Approx(x[i]).margin(1e-12));
      CHECK(L[1] == Approx(fx[i]).epsilon(1e-10).margin(1e-10));
      CHECK(L[2] == Approx(Af[i]).epsilon(1e-8).margin(1e-8));
    }
  }
}

TEST_CASE("6-DOF CBF rows: bank row matches an independent kinematic derivation",
          "[sixdof_cbf]") {
  BeaverDynamics dyn;
  SixDofCbfConfig cfg = baseCfg();
  cfg.row_alpha = cfg.row_V = cfg.row_beta = false;
  cfg.margin_force_terms = false;
  SixDofCbfFilter filt(dyn, cfg);
  std::mt19937 rng(23);
  for (int k = 0; k < 10; ++k) {
    const StateVec x = randomState(rng);
    const SurfVec up(0.02, -0.03, 0.01);
    const double dT = 0.5;
    const auto rs = filt.rows(x, dT, up);
    REQUIRE(rs.size() == 2);
    const SixDofCbfRow& r = rs[0];  // phi_max - phi
    // Hand: h = phi_max - phi; L_f h = -phidot; L_f^2 h = -d/dt(phidot) along
    // f, with phidot = p + (q sphi + r cphi) tan(theta). Only f's components
    // (the plant rows with frozen forces) and the kinematic chain rule are used.
    const SixDofCbfModel f(dyn, up, dT);
    const auto x0 = arr(x);
    const auto fx = f(x0);
    const double sp = std::sin(x[PHI]), cp = std::cos(x[PHI]);
    const double tt = std::tan(x[THETA]), ct = std::cos(x[THETA]);
    const double phidot = x[P] + (x[Q] * sp + x[R] * cp) * tt;
    CHECK(fx[PHI] == Approx(phidot).margin(1e-12));
    const double phiddot = fx[P] + (fx[Q] * sp + fx[R] * cp) * tt +
                           (x[Q] * cp - x[R] * sp) * fx[PHI] * tt +
                           (x[Q] * sp + x[R] * cp) * fx[THETA] / (ct * ct);
    CHECK(r.h == Approx(cfg.phi_max - x[PHI]).margin(1e-12));
    CHECK(r.Lfh == Approx(-phidot).margin(1e-10));
    CHECK(r.Lf2h == Approx(-phiddot).epsilon(1e-8).margin(1e-8));
    // Control row: a = -L_g L_f h = +(g_P + (g_Q sphi + g_R cphi) tan theta).
    const auto g = f.gColumns(x0);
    for (int kk = 0; kk < NUS; ++kk) {
      const double LgLf = -(g[kk][P] + (g[kk][Q] * sp + g[kk][R] * cp) * tt);
      CHECK(r.a[kk] == Approx(-LgLf).epsilon(1e-9).margin(1e-9));
    }
    // Sign semantics: a positive aileron (Cl_da < 0 -> rolls LEFT on the
    // Beaver) must relax the phi_max+ row (a < 0) and tighten phi_max- (a > 0).
    CHECK(rs[0].a[1] < 0.0);
    CHECK(rs[1].a[1] > 0.0);
  }
}

TEST_CASE("6-DOF CBF rows: AoA row first order by hand, second order vs flow oracle",
          "[sixdof_cbf]") {
  BeaverDynamics dyn;
  SixDofCbfConfig cfg = baseCfg();
  cfg.row_V = cfg.row_phi = cfg.row_beta = false;
  cfg.margin_force_terms = false;
  SixDofCbfFilter filt(dyn, cfg);
  std::mt19937 rng(24);
  for (int k = 0; k < 6; ++k) {
    const StateVec x = randomState(rng);
    const SurfVec up(0.0, 0.0, 0.0);
    const double dT = 0.5;
    const auto rs = filt.rows(x, dT, up);
    REQUIRE(rs.size() == 1);
    const SixDofCbfModel f(dyn, up, dT);
    const auto fx = f(arr(x));
    const double u0 = x[U], w0 = x[W];
    const double alphadot = (u0 * fx[W] - w0 * fx[U]) / (u0 * u0 + w0 * w0);
    CHECK(rs[0].Lfh == Approx(-alphadot).epsilon(1e-9).margin(1e-9));
    // L_f^2 h = d/dt (L_f h) along the flow: test-side central difference of
    // the exact first-order term along a short RK4 flow of f.
    auto Lfh_at = [&](const StateVec& xs) {
      const auto fs = f(arr(xs));
      return -(xs[U] * fs[W] - xs[W] * fs[U]) / (xs[U] * xs[U] + xs[W] * xs[W]);
    };
    const double h = 1e-4;
    CtrlVec u = CtrlVec::Zero(); u[DT] = dT;  // u_prev = 0: f == plant
    const StateVec xp = rk4(dyn, x, u, h), xm = rk4(dyn, x, u, -h);
    const double Lf2 = (Lfh_at(xp) - Lfh_at(xm)) / (2.0 * h);
    CHECK(rs[0].Lf2h == Approx(Lf2).epsilon(1e-5).margin(1e-6));
  }
}

TEST_CASE("6-DOF CBF: robustness margin equals the force-column sensitivity",
          "[sixdof_cbf]") {
  BeaverDynamics dyn;
  SixDofCbfConfig cfg = baseCfg();
  cfg.row_V = cfg.row_phi = cfg.row_beta = false;  // AoA only: elevator lift
  SixDofCbfFilter filt(dyn, cfg);
  std::mt19937 rng(25);
  const StateVec x = randomState(rng);
  const SurfVec up(0.0, 0.0, 0.0);
  const auto rs = filt.rows(x, 0.5, up);
  REQUIRE(rs.size() == 1);
  // Hand: D1_de = L_F h = grad(alpha_max - alpha) . F_de = -(u0/(u0^2+w0^2)) Fz_de/m
  // (h depends on u, w only; the elevator force column moves wdot).
  const SixDofCbfModel f(dyn, up, 0.5);
  const auto F = f.forceColumns(arr(x));
  const double u0 = x[U], w0 = x[W];
  const double D1_hand = -(u0 / (u0 * u0 + w0 * w0)) * F[0][W];
  CHECK(rs[0].D1[0] == Approx(D1_hand).epsilon(1e-9));
  CHECK(rs[0].D1[1] == Approx(0.0).margin(1e-12));   // aileron: side force only
  // D2 = L_F L_f h + L_f L_F h: the first term via lieAlong (verified
  // elsewhere), the second as a test-side flow difference of L_F h.
  const auto x0 = arr(x);
  auto b = [&](const auto& X) { using T = std::decay_t<decltype(X[0])>; using std::atan; return T(cfg.alpha_max) - atan(X[W] / X[U]); };
  const double LFLf = lieAlong<2, NX>(f, b, x0, F[0]);
  auto LFh_at = [&](const StateVec& xs) {
    const auto Fs = f.forceColumns(arr(xs));
    return -(xs[U] / (xs[U] * xs[U] + xs[W] * xs[W])) * Fs[0][W];
  };
  CtrlVec u = CtrlVec::Zero(); u[DT] = 0.5;
  const double hh = 1e-4;
  const double LfLF = (LFh_at(rk4(dyn, x, u, hh)) - LFh_at(rk4(dyn, x, u, -hh))) / (2.0 * hh);
  CHECK(rs[0].D2[0] == Approx(LFLf + LfLF).epsilon(1e-5).margin(1e-8));
  // Margin = (|D2| + (c1 + c2)|D1|) * rate dt, summed over the surfaces.
  const double du = cfg.limits.rate * cfg.dt;
  double hand = 0.0;
  for (int k = 0; k < NUS; ++k)
    hand += (std::abs(rs[0].D2[k]) + (cfg.c_alpha[0] + cfg.c_alpha[1]) * std::abs(rs[0].D1[k])) * du;
  CHECK(rs[0].margin == Approx(hand).epsilon(1e-9));
  CHECK(rs[0].margin > 0.0);
}

TEST_CASE("6-DOF CBF: hard envelope rows are forward invariant under a hostile nominal",
          "[sixdof_cbf][invariance]") {
  BeaverDynamics dyn;
  const TrimResult tr = beaverTrim(dyn, 40.0, -3.5 * kDeg);
  REQUIRE(tr.converged);
  struct Case { const char* name; int axis; double sign; };
  const Case cases[] = {{"bank", DA, -1.0},   // full aileron: Cl_da < 0 -> +da rolls left; -da rolls right
                        {"alpha", DE, -1.0},  // Cm_de < 0: -de pitches up -> alpha rises
                        {"beta", DR, -1.0}};
  for (const Case& c : cases) {
    for (double dt : {0.01, 0.0025}) {
      SixDofCbfConfig cfg = baseCfg();
      cfg.dt = dt;
      cfg.hard_phi = cfg.hard_alpha = cfg.hard_beta = true;
      cfg.hard_V = false;
      SixDofCbfFilter filt(dyn, cfg);
      StateVec x = tr.x;
      CtrlVec u_applied = tr.u;
      double min_h = 1e30, min_phi_h = 1e30, min_alpha_h = 1e30, min_beta_h = 1e30;
      int best_effort = 0;
      const int n = static_cast<int>(15.0 / dt);
      for (int k = 0; k < n; ++k) {
        SurfVec U_nom(tr.u[DE], tr.u[DA], tr.u[DR]);
        U_nom[c.axis == DE ? 0 : c.axis == DA ? 1 : 2] = c.sign * 25.0 * kDeg;
        const SurfVec up(u_applied[DE], u_applied[DA], u_applied[DR]);
        SixDofCbfDiag d;
        const SurfVec U = filt.filter(U_nom, x, tr.u[DT], up, &d);
        if (d.best_effort) ++best_effort;
        min_phi_h = std::min(min_phi_h, d.h_phi);
        min_alpha_h = std::min(min_alpha_h, d.h_alpha);
        min_beta_h = std::min(min_beta_h, d.h_beta);
        min_h = std::min(min_h, d.min_h);
        CtrlVec u = tr.u;
        u[DE] = U[0]; u[DA] = U[1]; u[DR] = U[2];
        x = rk4(dyn, x, u, dt);
        u_applied = u;
      }
      INFO(c.name << " dt=" << dt);
      // Every hard barrier stays non-negative up to the sampled-data slack of
      // one step at the envelope's own rate (bank: 0.1 deg; AoA/beta 0.05).
      CHECK(min_phi_h > -0.1 * kDeg);
      CHECK(min_alpha_h > -0.05 * kDeg);
      CHECK(min_beta_h > -0.05 * kDeg);
      // The hostile axis actually reached its boundary (the row did work).
      if (c.axis == DA) CHECK(min_phi_h < 0.5 * kDeg);
      if (c.axis == DE) CHECK(min_alpha_h < 0.5 * kDeg);
      if (c.axis == DR) CHECK(min_beta_h < 1.0 * kDeg);
      CHECK(best_effort == 0);
    }
  }
}

TEST_CASE("6-DOF CBF: feasibility sweep -- hard families with and without the rate box",
          "[sixdof_cbf]") {
  // The compatibility question a reviewer asks: is {hard rows + input box}
  // feasible everywhere inside the safe set? Measured over 200 random states
  // that are in-set at BOTH levels of every degree-2 row:
  //  * each family alone, deflection box only          -> 100 %
  //  * each family alone, with the one-step RATE box   -> reported (the AoA
  //    row is the weak one: a rate-limited surface cannot always meet a hard
  //    row within one 10 ms sample from an arbitrary in-set state -- it can
  //    over a few samples; the flight cases never hit best-effort)
  //  * all four hard at once, deflection box            -> 100 % (the rows
  //    ARE compatible; the joint shortfall seen with the rate box is the
  //    rate limit, not the rows)
  // Hence the policy: one hard row per scenario, the rest soft, best-effort
  // counted and reported, and the rate limit stated as the sampled-data
  // caveat.
  BeaverDynamics dyn;
  std::mt19937 rng(26);
  std::uniform_real_distribution<double> un(-1.0, 1.0);
  std::vector<StateVec> xs;
  std::vector<SurfVec> ups;
  {
    SixDofCbfConfig probe = baseCfg();
    SixDofCbfFilter pf(dyn, probe);
    while (xs.size() < 200) {
      const StateVec x = randomState(rng);
      const SurfVec up(0.1 * un(rng), 0.1 * un(rng), 0.1 * un(rng));
      bool in = true;
      for (const auto& r : pf.rows(x, 0.5, up)) {
        const double c1 = r.name.rfind("alpha", 0) == 0 ? probe.c_alpha[0]
                        : r.name.rfind("V_", 0) == 0 ? probe.c_V[0]
                        : r.name.rfind("phi", 0) == 0 ? probe.c_phi[0] : probe.c_beta[0];
        in = in && r.h > 0.0 && (r.Lfh + c1 * r.h) > 0.0;
      }
      if (in) { xs.push_back(x); ups.push_back(up); }
    }
  }
  auto sweep = [&](const SixDofCbfConfig& cfg) {
    SixDofCbfFilter filt(dyn, cfg);
    int ok = 0;
    for (size_t i = 0; i < xs.size(); ++i)
      if (filt.hardSetFeasible(xs[i], 0.5, ups[i])) ++ok;
    return ok;
  };
  const int n = static_cast<int>(xs.size());
  const char* names[4] = {"alpha", "V", "phi", "beta"};
  for (int fam = 0; fam < 4; ++fam) {
    SixDofCbfConfig cfg = baseCfg();
    cfg.hard_alpha = fam == 0; cfg.hard_V = fam == 1; cfg.hard_phi = fam == 2; cfg.hard_beta = fam == 3;
    cfg.rate_box = false;
    const int ok = sweep(cfg);
    INFO(names[fam] << " alone, deflection box: feasible " << ok << " / " << n);
    CHECK(ok == n);
    cfg.rate_box = true;
    const int ok_rate = sweep(cfg);
    WARN(names[fam] << " alone, with the one-step rate box: feasible " << ok_rate << " / " << n);
    CHECK(ok_rate > 0.7 * n);
  }
  SixDofCbfConfig all = baseCfg();
  all.hard_alpha = all.hard_V = all.hard_phi = all.hard_beta = true;
  all.rate_box = false;
  const int ok_all = sweep(all);
  INFO("all four families hard (deflection box): feasible " << ok_all << " / " << n);
  CHECK(ok_all == n);
}

// ---------------------------------------------------------------------------
// Phase 2: terrain keep-out row (relative degree 3).
// ---------------------------------------------------------------------------
namespace {
TerrainField testTerrain() {
  TerrainField tf;
  tf.bumps = {TerrainBump{"hill", 800.0, 0.0, 120.0, 250.0, 400.0, 0.2, 1},
              TerrainBump{"knob", 300.0, 150.0, 30.0, 40.0, 40.0, 0.0, 2}};
  return tf;
}
}  // namespace

TEST_CASE("6-DOF CBF terrain row: first order by hand, higher orders vs flow, degree-3 margin",
          "[sixdof_cbf][terrain]") {
  BeaverDynamics dyn;
  const TerrainField tf = testTerrain();
  SixDofCbfConfig cfg = baseCfg();
  cfg.row_alpha = cfg.row_V = cfg.row_phi = cfg.row_beta = false;
  cfg.row_terrain = true; cfg.terrain = &tf; cfg.terrain_margin = 4.0; cfg.terrain_gate = 1e9;
  cfg.terrain_degree3 = true;
  cfg.c_terrain = {0.8, 1.1, 1.4};
  SixDofCbfFilter filt(dyn, cfg);
  std::mt19937 rng(31);
  for (int k = 0; k < 6; ++k) {
    StateVec x = randomState(rng);
    x[XN] = 400.0 + 300.0 * (k / 5.0); x[Y] = 60.0 - 20.0 * k; x[H] = 90.0 + 10.0 * k;
    const SurfVec up(0.0, 0.0, 0.0);
    const double dT = 0.5;
    const auto rs = filt.rows(x, dT, up);
    REQUIRE(rs.size() == 1);
    const SixDofCbfRow& r = rs[0];
    CHECK(r.degree == 3);
    const SixDofCbfModel f(dyn, up, dT);
    const auto fx = f(arr(x));
    const Eigen::Vector2d gT = tf.gradient(x[XN], x[Y]);
    // h = H - h_T(N, E) - m;  L_f h = Hdot - dh_T/dN Ndot - dh_T/dE Edot.
    CHECK(r.h == Approx(x[H] - tf.height(x[XN], x[Y]) - 4.0).margin(1e-12));
    CHECK(r.Lfh == Approx(fx[H] - gT[0] * fx[XN] - gT[1] * fx[Y]).epsilon(1e-9).margin(1e-9));
    // L_f^2 h and L_f^3 h: flow differences of the exact lower orders
    // (u_prev = 0 -> f is the plant at surfaces 0).
    CtrlVec u = CtrlVec::Zero(); u[DT] = dT;
    auto Lfh_at = [&](const StateVec& xs) {
      const auto fs = f(arr(xs));
      const Eigen::Vector2d g = tf.gradient(xs[XN], xs[Y]);
      return fs[H] - g[0] * fs[XN] - g[1] * fs[Y];
    };
    const double hh = 1e-3;
    const StateVec xp = rk4(dyn, x, u, hh), xm = rk4(dyn, x, u, -hh);
    const double Lf2 = (Lfh_at(xp) - Lfh_at(xm)) / (2.0 * hh);
    CHECK(r.Lf2h == Approx(Lf2).epsilon(1e-5).margin(1e-6));
    auto Lf2_at = [&](const StateVec& xs) {
      const StateVec a = rk4(dyn, xs, u, hh), b = rk4(dyn, xs, u, -hh);
      return (Lfh_at(a) - Lfh_at(b)) / (2.0 * hh);
    };
    const double Lf3 = (Lf2_at(xp) - Lf2_at(xm)) / (2.0 * hh);
    CHECK(r.Lf3h == Approx(Lf3).epsilon(1e-3).margin(1e-4));
    // Margin structure for degree 3: sum_k (e2|D1| + e1|D2| + |D3|) du with
    // e1 = c1+c2+c3, e2 = c1c2+c1c3+c2c3.
    const double c1 = 0.8, c2 = 1.1, c3 = 1.4;
    const double e1 = c1 + c2 + c3, e2 = c1 * c2 + c1 * c3 + c2 * c3;
    const double du = cfg.limits.rate * cfg.dt;
    double hand = 0.0;
    for (int kk = 0; kk < NUS; ++kk)
      hand += (e2 * std::abs(r.D1[kk]) + e1 * std::abs(r.D2[kk]) + std::abs(r.D3[kk])) * du;
    CHECK(r.margin == Approx(hand).epsilon(1e-9));
    // Position rows do not leak at first order (h depends on H, N, E only).
    for (int kk = 0; kk < NUS; ++kk) CHECK(r.D1[kk] == Approx(0.0).margin(1e-12));
    // Relative degree 3: every surface has a nonzero control row.
    CHECK(std::abs(r.a[0]) > 1e-6);
    CHECK(std::abs(r.a[1]) > 1e-9);
  }
}

TEST_CASE("6-DOF CBF terrain braking row: closed-form gradient and first order by hand",
          "[sixdof_cbf][terrain]") {
  BeaverDynamics dyn;
  const TerrainField tf = testTerrain();
  // gradientT (closed form) == gradient (autodiff).
  for (auto pt : {std::pair{700.0, 50.0}, std::pair{320.0, 130.0}, std::pair{-100.0, 400.0}}) {
    const auto g = tf.gradientT<double>(pt.first, pt.second);
    const Eigen::Vector2d ga = tf.gradient(pt.first, pt.second);
    CHECK(g[0] == Approx(ga[0]).margin(1e-12));
    CHECK(g[1] == Approx(ga[1]).margin(1e-12));
  }
  SixDofCbfConfig cfg = baseCfg();
  cfg.row_alpha = cfg.row_V = cfg.row_phi = cfg.row_beta = false;
  cfg.row_terrain = true; cfg.terrain = &tf; cfg.terrain_margin = 4.0; cfg.terrain_gate = 1e9;
  cfg.terrain_a_brk = 3.0; cfg.terrain_v_safe = 0.5; cfg.terrain_z_eps = 1.0; cfg.terrain_k_neg = 1.0;
  SixDofCbfFilter filt(dyn, cfg);
  std::mt19937 rng(32);
  for (int k = 0; k < 5; ++k) {
    StateVec x = randomState(rng);
    x[XN] = 500.0 + 60.0 * k; x[Y] = 40.0 - 15.0 * k; x[H] = 80.0 + 12.0 * k;
    const SurfVec up(0.0, 0.0, 0.0);
    const auto rs = filt.rows(x, 0.5, up);
    REQUIRE(rs.size() == 1);
    CHECK(rs[0].degree == 2);
    const SixDofCbfModel f(dyn, up, 0.5);
    const auto fx = f(arr(x));
    const Eigen::Vector2d gT = tf.gradient(x[XN], x[Y]);
    const double z = x[H] - tf.height(x[XN], x[Y]) - 4.0;
    const double zdot = fx[H] - gT[0] * fx[XN] - gT[1] * fx[Y];
    const double zp = 0.5 * (z + std::sqrt(z * z + 1.0));
    CHECK(rs[0].h == Approx(zdot + std::sqrt(0.25 + 6.0 * zp) + 1.0 * (z - zp)).epsilon(1e-9));
    // Relative degree 2: the elevator control row is nonzero (zdot -> forces
    // -> pitch rate), and the row leaks at first order through the
    // elevator lift (D1 != 0) -- unlike the pure clearance row.
    CHECK(std::abs(rs[0].a[0]) > 1e-6);
    CHECK(std::abs(rs[0].D1[0]) > 1e-9);
  }
}

// What this test CLAIMS: against a full nose-down hostile nominal at climb
// power, flying at a 120 m hill, the CG clearance never drops below the
// margin minus 0.5 m (sampled-data slack) at 10 ms and at 2.5 ms. What it
// does NOT claim: that the barrier quantity b stays non-negative -- the
// climb-capability lookahead is conservative and the hostile dive drives b
// well below zero while the row runs best-effort; the physical clearance is
// nevertheless held, which is the point of the capability form. Nor does it
// claim the aircraft clears the hill within the window (it climbs the slope
// slowly against the nose-down nominal).
TEST_CASE("6-DOF CBF terrain braking row holds the clearance under a nose-down hostile nominal",
          "[sixdof_cbf][terrain][invariance]") {
  BeaverDynamics dyn;
  const TrimResult tr = beaverTrim(dyn, 40.0, -3.5 * kDeg);
  REQUIRE(tr.converged);
  const TerrainField tf = testTerrain();
  for (double dt : {0.01, 0.0025}) {
    SixDofCbfConfig cfg = baseCfg();
    cfg.dt = dt;
    cfg.row_terrain = true; cfg.terrain = &tf; cfg.terrain_margin = 3.0; cfg.terrain_gate = 1e9;
    cfg.hard_terrain = true;
    cfg.terrain_a_brk = 3.0; cfg.terrain_v_safe = 0.1; cfg.terrain_k_neg = 1.0; cfg.c_terrain = {1.0, 1.0, 1.0};
    cfg.terrain_lookahead = 1200.0;  // climb-capability lookahead (see the header)
    cfg.terrain_climb_grad = 0.05;
    cfg.hard_alpha = false; cfg.hard_phi = false; cfg.hard_beta = false; cfg.hard_V = false;
    SixDofCbfFilter filt(dyn, cfg);
    StateVec x = tr.x;
    x[XN] = -400.0; x[Y] = 0.0; x[H] = 90.0;   // heading north at the hill, 90 m up (inside the set)
    CtrlVec u_applied = tr.u;
    double min_h = 1e30, min_clr = 1e30;
    int be = 0, active = 0;
    const int n = static_cast<int>(30.0 / dt);
    for (int k = 0; k < n; ++k) {
      SurfVec U_nom(+15.0 * kDeg, tr.u[DA], tr.u[DR]);  // Cm_de < 0: +de = nose DOWN
      const SurfVec up(u_applied[DE], u_applied[DA], u_applied[DR]);
      SixDofCbfDiag d;
      const double dT = 0.9;   // the throttle is the nominal's: give it climb power
      const SurfVec U = filt.filter(U_nom, x, dT, up, &d);
      if (d.best_effort) ++be;
      if (d.active) ++active;
      min_h = std::min(min_h, d.h_terrain);
      min_clr = std::min(min_clr, x[H] - tf.height(x[XN], x[Y]));
      if (x[H] < -50.0) break;  // through the floor: give up
      CtrlVec u = tr.u;
      u[DE] = U[0]; u[DA] = U[1]; u[DR] = U[2]; u[DT] = dT;
      x = rk4(dyn, x, u, dt);
      u_applied = u;
      if (x[XN] > 1400.0) break;  // past the hill
    }
    WARN("dt=" << dt << " min b=" << min_h << " min CG clearance=" << min_clr
         << " m (margin 3)  best-effort steps=" << be << "  active=" << active);
    CHECK(min_clr > 2.5);         // the physical claim: never below margin - 0.5 m
    CHECK(active > 100);          // the row did the work against full nose-down
    CHECK(x[H] > 0.0);            // still airborne at the end of the window
  }
}
