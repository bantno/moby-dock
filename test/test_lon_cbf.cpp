#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "autoland/aero_table.hpp"
#include "autoland/config.hpp"
#include "autoland/hocbf.hpp"
#include "autoland/impact_barrier.hpp"
#include "autoland/lon_augmented.hpp"
#include "autoland/lon_cbf_filter.hpp"
#include "autoland/mixing.hpp"

using namespace autoland;
using Catch::Approx;

namespace {
const std::string kData = AUTOLAND_DATA_DIR;
constexpr double kDeg = M_PI / 180.0;

struct Setup {
  AeroTable table;
  AircraftConfig cfg;
  Mixing mx;
  Setup()
      : table(AeroTable::fromFile(kData + "/example.stab")),
        cfg(loadAircraftConfig(kData + "/aircraft.yaml")),
        mx(Mixing::build(cfg, table)) {}
};

// ---- Finite-difference flow oracle for the drift Lie stack ------------------
// Independent cross-check of lieDrift<3>. With drift-only dynamics Xdot = f(X)
// (control U = 0), L_f^k b = d^k/dt^k b(X(t))|_{t=0}. We integrate the SAME
// frozen-aero field the engine differentiates (AeroLocal `a` held fixed at X0,
// never re-looked-up along the flow), sample phi(t) = b(X(t)), and
// central-difference in time. This touches none of lie_taylor.hpp's
// Taylor/Picard machinery, so it catches a bug in that custom code.

// One RK4 step of the pure drift (control = 0) under the frozen aero `a`.
LonStateVec driftStep(const AeroLocal& a, const LonStateVec& X, double dt) {
  const LonCtrlVec Z = LonCtrlVec::Zero();
  const LonStateVec k1 = lonXdot(a, X, Z);
  const LonStateVec k2 = lonXdot(a, X + 0.5 * dt * k1, Z);
  const LonStateVec k3 = lonXdot(a, X + 0.5 * dt * k2, Z);
  const LonStateVec k4 = lonXdot(a, X + dt * k3, Z);
  return X + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

// Integrate the drift from X0 to time t (t may be negative) with n substeps,
// fine enough that RK4 error (O(dt^4)) is negligible vs the stencil truncation.
LonStateVec flowTo(const AeroLocal& a, const LonStateVec& X0, double t, int n) {
  LonStateVec X = X0;
  const double dt = t / n;
  for (int i = 0; i < n; ++i) X = driftStep(a, X, dt);
  return X;
}

// {b, L_f b, L_f^2 b, L_f^3 b} from central time-differences of phi(t)=b(X(t)),
// Richardson-extrapolated over (h, h/2) to cancel the O(h^2) stencil error.
template <class Barrier>
std::array<double, 4> lieDriftOracle(const AeroLocal& a, const Barrier& b,
                                     const LonStateVec& X0, double h, int nsub) {
  auto phi = [&](double t) { return b(toArray(flowTo(a, X0, t, nsub))); };
  auto stencils = [&](double hh) {
    const double f0 = phi(0.0);
    const double fp1 = phi(hh), fm1 = phi(-hh);
    const double fp2 = phi(2.0 * hh), fm2 = phi(-2.0 * hh);
    std::array<double, 4> D;
    D[0] = f0;
    D[1] = (fp1 - fm1) / (2.0 * hh);
    D[2] = (fp1 - 2.0 * f0 + fm1) / (hh * hh);
    D[3] = (fp2 - 2.0 * fp1 + 2.0 * fm1 - fm2) / (2.0 * hh * hh * hh);
    return D;
  };
  const std::array<double, 4> Dh = stencils(h);
  const std::array<double, 4> Dh2 = stencils(0.5 * h);
  std::array<double, 4> L;
  L[0] = Dh[0];
  for (int k = 1; k < 4; ++k) L[k] = (4.0 * Dh2[k] - Dh[k]) / 3.0;  // O(h^4)
  return L;
}
}  // namespace

// The recovery-set barriers reduce to their literal definitions, and their exact
// Lie stacks have the relative degree the wiring assumes: stall / nose-up are
// degree 2 (elevator only, thrust column ~0); energy is degree 3 (BOTH controls,
// since thrust reaches airspeed only through the T,Tdot,Tddot chain).
TEST_CASE("recovery barriers: closed form + relative degree", "[lon_cbf]") {
  Setup s;
  const double V = 18.0, gamma = -5.0 * kDeg, theta = -2.0 * kDeg;
  const double q = 0.3, alpha = theta - gamma;
  AeroLocal a = makeAeroLocal(s.table, s.mx, s.cfg, V, alpha);
  LonStateVec X;
  X << 30.0, V, gamma, theta, q, 5.0, 1.0;  // nonzero q, T, Tdot
  const LonStateVec xd = lonXdot(a, X, LonCtrlVec::Zero());

  SECTION("stall (deg 2, elevator only)") {
    const double alpha_stall = 11.0 * kDeg, margin = 2.0 * kDeg;
    const StallBarrier b = makeStallBarrier(alpha_stall, margin);
    CHECK(b(toArray(X)) == Approx((alpha_stall - margin) - alpha).epsilon(1e-12));
    auto lie = barrierLie<2>(a, b, X);
    // d/dt[alpha_lim - (theta - gamma)] = -thetadot + gammadot = -q + gammadot.
    CHECK(lie.Lf[1] == Approx(-q + xd[LGAM]).epsilon(1e-9));
    CHECK(std::abs(lie.LgLf(0, LTDDOT)) < 1e-9);  // thrust absent at degree 2
    CHECK(std::abs(lie.LgLf(0, LDE)) > 1e-6);     // elevator has authority
  }
  SECTION("nose-up (deg 2, elevator only)") {
    const double theta_min = 3.0 * kDeg;
    const NoseUpBarrier b = makeNoseUpBarrier(theta_min);
    CHECK(b(toArray(X)) == Approx(theta - theta_min).epsilon(1e-12));
    auto lie = barrierLie<2>(a, b, X);
    CHECK(lie.Lf[1] == Approx(q).epsilon(1e-9));  // d/dt(theta) = thetadot = q
    CHECK(std::abs(lie.LgLf(0, LTDDOT)) < 1e-9);
    CHECK(std::abs(lie.LgLf(0, LDE)) > 1e-6);
  }
  SECTION("energy (deg 3, both controls)") {
    const double V_td_max = 14.0, g_eff = 16.0;
    const EnergyBarrier b = makeEnergyBarrier(a, V_td_max, g_eff);
    CHECK(b(toArray(X)) == Approx(0.5 * (V_td_max * V_td_max - V * V) +
                                  (g_eff - a.g) * X[LH])
                               .epsilon(1e-9));
    auto lie = barrierLie<3>(a, b, X);
    // d/dt[1/2(Vtd^2 - V^2) + (g_eff-g)h] = -V*Vdot + (g_eff-g)*hdot.
    CHECK(lie.Lf[1] ==
          Approx(-V * xd[LV] + (g_eff - a.g) * xd[LH]).epsilon(1e-7));
    CHECK(std::abs(lie.LgLf(0, LDE)) > 1e-6);     // elevator enters at degree 3
    CHECK(std::abs(lie.LgLf(0, LTDDOT)) > 1e-9);  // thrust enters at degree 3
  }
}

// The generic HOCBF builder must reproduce the hand-expanded linear-class-K
// constraint, and the actuator rows the doc section 3.4 inequalities.
TEST_CASE("hocbfRow linear class-K expansion", "[lon_cbf]") {
  std::vector<double> Lf = {1.0, 2.0, 3.0, 4.0};  // b, L_f b, L_f^2 b, L_f^3 b
  Eigen::Matrix<double, 1, NUA> LgLf;
  LgLf << 0.5, -0.7;
  const double c1 = 1.5, c2 = 2.0, c3 = 0.8;
  HocbfRow row = hocbfRow(Lf, LgLf, {c1, c2, c3});

  const double drift = Lf[3] + (c1 + c2 + c3) * Lf[2] +
                       (c1 * c2 + c1 * c3 + c2 * c3) * Lf[1] + (c1 * c2 * c3) * Lf[0];
  CHECK(row.rhs == Approx(drift));
  CHECK(row.a(0, LDE) == Approx(-0.5));
  CHECK(row.a(0, LTDDOT) == Approx(0.7));
}

// A disabled filter is an exact pass-through.
TEST_CASE("disabled lon filter passes the nominal through", "[lon_cbf]") {
  Setup s;
  LonCBFConfig cfg;
  cfg.enabled = false;
  LonCBFFilter filter(cfg);
  LonStateVec X;
  X << 30.0, 18.0, -3.0 * kDeg, -1.8 * kDeg, 0.0, 2.0, 0.0;
  LonCtrlVec Un;
  Un << 0.1, -5.0;
  LonCtrlVec U = filter.filter(Un, X, s.table, s.mx, s.cfg);
  CHECK(U[LDE] == Approx(Un[LDE]));
  CHECK(U[LTDDOT] == Approx(Un[LTDDOT]));
}

// ---- Hard-barrier enforcement: after filtering, the assembled HOCBF row holds
// even when the nominal drives hard onto the violating side. Each test isolates
// one barrier (all others disabled) so the correction is unambiguous.

TEST_CASE("lon filter enforces the stall barrier when hard", "[lon_cbf]") {
  Setup s;
  // State close to the stall limit (alpha ~ 8.5 deg vs alpha_lim = 9 deg).
  const double V = 16.0, gamma = -3.0 * kDeg, theta = 5.5 * kDeg;
  AeroLocal a = makeAeroLocal(s.table, s.mx, s.cfg, V, theta - gamma);
  LonStateVec X;
  X << 20.0, V, gamma, theta, 0.0, 2.0, 0.0;

  LonCBFConfig cfg;
  cfg.stall = true; cfg.stall_hard = true;
  cfg.noseup = false; cfg.energy = false; cfg.impact = false;
  cfg.thrust_limits = false;
  LonCBFFilter filter(cfg);

  const StallBarrier b = makeStallBarrier(cfg.alpha_stall, cfg.stall_margin);
  auto lie = barrierLie<2>(a, b, X);
  std::vector<double> Lf(lie.Lf.begin(), lie.Lf.end());
  HocbfRow row = hocbfRow(Lf, lie.LgLf, {cfg.c_stall[0], cfg.c_stall[1]});
  const double a_de = row.a(0, LDE);
  REQUIRE(std::abs(a_de) > 1e-6);
  CHECK(std::abs(row.a(0, LTDDOT)) < 1e-9);  // thrust absent from the degree-2 row

  LonCtrlVec Un;
  Un << (a_de > 0 ? cfg.de_max : cfg.de_min), 0.0;  // most-violating elevator (into stall)
  REQUIRE(a_de * Un[LDE] > row.rhs + 1e-6);

  LonCtrlVec U = filter.filter(Un, X, s.table, s.mx, s.cfg);
  CHECK_FALSE(filter.lastRecovery());
  CHECK((U - Un).norm() > 1e-6);
  const double lhs = row.a(0, LDE) * U[LDE] + row.a(0, LTDDOT) * U[LTDDOT];
  CHECK(lhs <= row.rhs + 1e-6);
}

TEST_CASE("lon filter enforces the energy barrier when hard", "[lon_cbf]") {
  Setup s;
  // Hot + low: V just under the height-scheduled cap so the ceiling is active.
  const double V = 16.0, gamma = -3.0 * kDeg, theta = 0.0;
  AeroLocal a = makeAeroLocal(s.table, s.mx, s.cfg, V, theta - gamma);
  LonStateVec X;
  X << 5.0, V, gamma, theta, 0.0, 2.0, 0.0;

  LonCBFConfig cfg;
  cfg.energy = true; cfg.energy_hard = true;
  cfg.stall = false; cfg.noseup = false; cfg.impact = false;
  cfg.thrust_limits = false;
  LonCBFFilter filter(cfg);

  const EnergyBarrier b = makeEnergyBarrier(a, cfg.V_td_max, cfg.g_eff);
  auto lie = barrierLie<3>(a, b, X);
  std::vector<double> Lf(lie.Lf.begin(), lie.Lf.end());
  HocbfRow row = hocbfRow(Lf, lie.LgLf, {cfg.c_energy[0], cfg.c_energy[1], cfg.c_energy[2]});
  const double a_de = row.a(0, LDE), a_td = row.a(0, LTDDOT);
  REQUIRE((std::abs(a_de) > 1e-9 || std::abs(a_td) > 1e-9));

  // Most-violating in-box control: push both controls to the energy-raising edge.
  LonCtrlVec Un;
  Un << (a_de > 0 ? cfg.de_max : cfg.de_min),
        (a_td > 0 ? cfg.Tddot_max : cfg.Tddot_min);
  REQUIRE(a_de * Un[LDE] + a_td * Un[LTDDOT] > row.rhs + 1e-6);

  LonCtrlVec U = filter.filter(Un, X, s.table, s.mx, s.cfg);
  CHECK_FALSE(filter.lastRecovery());
  CHECK((U - Un).norm() > 1e-6);
  const double lhs = row.a(0, LDE) * U[LDE] + row.a(0, LTDDOT) * U[LTDDOT];
  CHECK(lhs <= row.rhs + 1e-6);
}

TEST_CASE("lon filter enforces the nose-up barrier when hard", "[lon_cbf]") {
  Setup s;
  // Below h_noseup (default 3 m) so the row is assembled; theta just above floor.
  const double V = 15.0, gamma = -3.0 * kDeg, theta = 3.5 * kDeg;
  AeroLocal a = makeAeroLocal(s.table, s.mx, s.cfg, V, theta - gamma);
  LonStateVec X;
  X << 2.0, V, gamma, theta, 0.0, 2.0, 0.0;

  LonCBFConfig cfg;
  cfg.noseup = true; cfg.noseup_hard = true; cfg.h_noseup = 3.0;
  cfg.stall = false; cfg.energy = false; cfg.impact = false;
  cfg.thrust_limits = false;
  LonCBFFilter filter(cfg);

  const NoseUpBarrier b = makeNoseUpBarrier(cfg.theta_min);
  auto lie = barrierLie<2>(a, b, X);
  std::vector<double> Lf(lie.Lf.begin(), lie.Lf.end());
  HocbfRow row = hocbfRow(Lf, lie.LgLf, {cfg.c_noseup[0], cfg.c_noseup[1]});
  const double a_de = row.a(0, LDE);
  REQUIRE(std::abs(a_de) > 1e-6);
  CHECK(std::abs(row.a(0, LTDDOT)) < 1e-9);

  LonCtrlVec Un;
  Un << (a_de > 0 ? cfg.de_max : cfg.de_min), 0.0;  // most-violating (nose-down)
  REQUIRE(a_de * Un[LDE] > row.rhs + 1e-6);

  LonCtrlVec U = filter.filter(Un, X, s.table, s.mx, s.cfg);
  CHECK_FALSE(filter.lastRecovery());
  CHECK((U - Un).norm() > 1e-6);
  const double lhs = row.a(0, LDE) * U[LDE] + row.a(0, LTDDOT) * U[LTDDOT];
  CHECK(lhs <= row.rhs + 1e-6);
}

// The thrust actuator / HOCBF-validity guards (always hard) hold when the nominal
// tries to drive the thrust state out of [0, Tmax].
TEST_CASE("lon filter keeps thrust within the actuator limits", "[lon_cbf]") {
  Setup s;
  LonCBFConfig cfg;
  cfg.stall = false; cfg.noseup = false; cfg.energy = false; cfg.impact = false;
  cfg.thrust_limits = true;  // min/max thrust rows (always hard)
  LonCBFFilter filter(cfg);

  LonStateVec X;
  X << 20.0, 16.0, -3.0 * kDeg, 0.0, 0.0, 0.2, 0.0;  // low thrust near the floor
  LonCtrlVec Un;
  Un << 0.0, -500.0;  // slam Tddot negative -> would drive T below 0
  LonCtrlVec U = filter.filter(Un, X, s.table, s.mx, s.cfg);

  CHECK_FALSE(filter.lastRecovery());
  // Min-thrust row (hard) must hold: -Tddot <= (c11+c12)Tdot + c11 c12 T.
  HocbfRow tmin = thrustMinRow(X, cfg.c_thrust_min[0], cfg.c_thrust_min[1]);
  CHECK(tmin.a(0, LTDDOT) * U[LTDDOT] <= tmin.rhs + 1e-6);
}

// Independently cross-check the recovery barriers' DRIFT Lie stacks against the
// finite-difference flow oracle (catches a bug in the bespoke Taylor-jet engine).
// State carries nonzero q, T, Tdot so the higher-order terms exercise the
// thrust/pitch couplings.
TEST_CASE("recovery-barrier drift Lie stacks match a finite-difference flow oracle",
          "[lon_cbf]") {
  Setup s;
  const double V = 18.0, gamma = -5.0 * kDeg, theta = -2.0 * kDeg;
  const double q = 0.3, T = 5.0, Tdot = 1.0, alpha = theta - gamma;
  AeroLocal a = makeAeroLocal(s.table, s.mx, s.cfg, V, alpha);
  LonStateVec X;
  X << 30.0, V, gamma, theta, q, T, Tdot;
  const double h = 1e-2;
  const int nsub = 200;

  // Loose 1e-3 relative guard: the finite-difference ORACLE (even Richardson-
  // extrapolated) is the accuracy floor, not the engine. Still catches any
  // structural bug (sign flip, factor of 2, missing term are all several percent).
  SECTION("stall (deg 2)") {
    const StallBarrier b = makeStallBarrier(11.0 * kDeg, 2.0 * kDeg);
    const auto lie = barrierLie<2>(a, b, X);
    const auto orc = lieDriftOracle(a, b, X, h, nsub);
    CHECK(lie.Lf[0] == Approx(orc[0]).epsilon(1e-9).margin(1e-9));
    CHECK(lie.Lf[1] == Approx(orc[1]).epsilon(1e-3).margin(1e-6));
    CHECK(lie.Lf[2] == Approx(orc[2]).epsilon(1e-3).margin(1e-6));
  }
  SECTION("nose-up (deg 2)") {
    const NoseUpBarrier b = makeNoseUpBarrier(3.0 * kDeg);
    const auto lie = barrierLie<2>(a, b, X);
    const auto orc = lieDriftOracle(a, b, X, h, nsub);
    CHECK(lie.Lf[0] == Approx(orc[0]).epsilon(1e-9).margin(1e-9));
    CHECK(lie.Lf[1] == Approx(orc[1]).epsilon(1e-3).margin(1e-6));
    CHECK(lie.Lf[2] == Approx(orc[2]).epsilon(1e-3).margin(1e-6));
  }
  SECTION("energy (deg 3)") {
    const EnergyBarrier b = makeEnergyBarrier(a, 14.0, 16.0);
    const auto lie = barrierLie<3>(a, b, X);
    const auto orc = lieDriftOracle(a, b, X, h, nsub);
    CHECK(lie.Lf[0] == Approx(orc[0]).epsilon(1e-9).margin(1e-9));
    CHECK(lie.Lf[1] == Approx(orc[1]).epsilon(1e-3).margin(1e-6));
    CHECK(lie.Lf[2] == Approx(orc[2]).epsilon(1e-3).margin(1e-6));
    CHECK(lie.Lf[3] == Approx(orc[3]).epsilon(1e-3).margin(1e-6));
  }
}

// =============================================================================
// Hydrodynamic impact-load barrier (NACA TN 1516).
// =============================================================================

// clfLookup reproduces the generated table endpoints, is monotone, and clamps
// outside the practical range. (The paper anchor Clf(0)=0.6123 is enforced by
// scripts/precompute_impact_clf.py, which generates the table; kappa=0 is below
// the table's [0.2,10] range and clamps to the first entry at run time.)
TEST_CASE("clfLookup is monotone and clamps to the practical range", "[lon_cbf]") {
  CHECK(clfLookup(kImpactClfKappa[0]).value == Approx(kImpactClfVal[0]).epsilon(1e-9));
  CHECK(clfLookup(kImpactClfKappa[kImpactClfN - 1]).value ==
        Approx(kImpactClfVal[kImpactClfN - 1]).epsilon(1e-9));
  CHECK(clfLookup(0.0).value == Approx(kImpactClfVal[0]).epsilon(1e-12));        // clamp lo
  CHECK(clfLookup(50.0).value == Approx(kImpactClfVal[kImpactClfN - 1]).epsilon(1e-12));  // clamp hi
  double prev = -1e9;
  for (int i = 0; i < kImpactClfN; ++i) { CHECK(kImpactClfVal[i] > prev); prev = kImpactClfVal[i]; }
  const double km = 0.5 * (kImpactClfKappa[3] + kImpactClfKappa[4]);  // interpolated point
  const ClfLocal m = clfLookup(km);
  CHECK(m.value > kImpactClfVal[3]);
  CHECK(m.value < kImpactClfVal[4]);
}

// The barrier value matches its closed form (n_limit - K0 ydot0^2 Clf) + Phi(z)
// exactly, using the same frozen constants (catches a sign/term error). Also
// pins f(beta) = 3.0 at beta = 22.5 deg (eq 45).
TEST_CASE("impact barrier value matches its closed form", "[lon_cbf]") {
  Setup s;
  const double V = 16.0, gamma = -4.0 * kDeg, theta = 3.0 * kDeg, h = 4.0;
  const double beta = 22.5 * kDeg, n_limit = 3.0, rho_w = 1000.0;
  const double Nb = 10.0, zs = 2.0, tau_keel = 0.0, eps_g0 = 0.02;
  AeroLocal a = makeAeroLocal(s.table, s.mx, s.cfg, V, theta - gamma);
  LonStateVec X;
  X << h, V, gamma, theta, 0.0, 2.0, 0.0;

  const ImpactLoadBarrier b =
      makeImpactLoadBarrier(a, n_limit, beta, rho_w, Nb, zs, tau_keel, eps_g0, X);

  const double gamma0 = -gamma;
  const double sg0 = std::sqrt(std::sin(gamma0) * std::sin(gamma0) + eps_g0 * eps_g0);
  const double tau = theta - tau_keel;
  const double kappa = std::sin(tau) * std::cos(tau + gamma0) / sg0;
  const double ydot0 = -V * std::sin(gamma);
  const double Clf = b.Clf0 + b.dClf_dk * (kappa - b.kappa0);
  const double n_peak = b.K0 * ydot0 * ydot0 * Clf;
  const double Phi = Nb * (1.0 - std::exp(-h / zs));
  CHECK(b(toArray(X)) == Approx((n_limit - n_peak) + Phi).epsilon(1e-12));
  CHECK((M_PI / (2.0 * beta) - 1.0) == Approx(3.0).epsilon(1e-9));  // f(beta), eq 45
}

// Relative degree: a degree-2 HOCBF. The elevator appears in L_g L_f b; thrust
// (Tddot, relative degree 3) does NOT (its column is ~0). This is exactly what
// lets the degree-2 row be flare-enforced while reusing barrierLie<2>.
TEST_CASE("impact barrier is relative degree 2 via the elevator", "[lon_cbf]") {
  Setup s;
  const double V = 16.0, gamma = -4.0 * kDeg, theta = 3.0 * kDeg;
  AeroLocal a = makeAeroLocal(s.table, s.mx, s.cfg, V, theta - gamma);
  LonStateVec X;
  X << 4.0, V, gamma, theta, 0.2, 2.0, 0.5;  // nonzero q, T, Tdot
  const ImpactLoadBarrier b =
      makeImpactLoadBarrier(a, 3.0, 22.5 * kDeg, 1000.0, 10.0, 2.0, 0.0, 0.02, X);
  auto lie = barrierLie<2>(a, b, X);
  CHECK(std::abs(lie.LgLf(0, LTDDOT)) < 1e-9);  // thrust drops out at degree 2
  CHECK(std::abs(lie.LgLf(0, LDE)) > 1e-6);     // elevator has flare authority
}

// Phi(z) makes the barrier touchdown-only: at altitude b ~ (n_limit-n_peak)+Nb
// (inactive); at z=0 Phi=0 so b = n_limit - n_peak (the true contact constraint).
// The frozen constants depend on (V,theta,gamma) only, so b(h) - b(0) = Phi(h).
TEST_CASE("impact barrier height relaxation engages only near the surface", "[lon_cbf]") {
  Setup s;
  const double V = 16.0, gamma = -4.0 * kDeg, theta = 3.0 * kDeg;
  const double Nb = 10.0, zs = 2.0;
  AeroLocal a = makeAeroLocal(s.table, s.mx, s.cfg, V, theta - gamma);
  auto barrier_at = [&](double h) {
    LonStateVec X;
    X << h, V, gamma, theta, 0.0, 2.0, 0.0;
    const ImpactLoadBarrier b =
        makeImpactLoadBarrier(a, 3.0, 22.5 * kDeg, 1000.0, Nb, zs, 0.0, 0.02, X);
    return b(toArray(X));
  };
  const double b0 = barrier_at(0.0);     // Phi = 0
  const double bhi = barrier_at(20.0);   // Phi ~ Nb
  CHECK(bhi - b0 == Approx(Nb * (1.0 - std::exp(-20.0 / zs))).epsilon(1e-9));
  CHECK(bhi > b0 + 0.9 * Nb);            // strongly relaxed aloft
}

// Independently cross-check the impact barrier's degree-2 drift stack
// {b, L_f b, L_f^2 b} against the finite-difference flow oracle. Also exercises
// the exp() Taylor overload (through Phi) along a descending flow.
TEST_CASE("impact barrier drift Lie stack matches the flow oracle", "[lon_cbf]") {
  Setup s;
  const double V = 16.0, gamma = -4.0 * kDeg, theta = 3.0 * kDeg;
  AeroLocal a = makeAeroLocal(s.table, s.mx, s.cfg, V, theta - gamma);
  LonStateVec X;
  X << 4.0, V, gamma, theta, 0.2, 2.0, 0.5;  // descending; nonzero q, T, Tdot
  const ImpactLoadBarrier b =
      makeImpactLoadBarrier(a, 3.0, 22.5 * kDeg, 1000.0, 10.0, 2.0, 0.0, 0.02, X);
  const auto lie = barrierLie<2>(a, b, X);
  const auto orc = lieDriftOracle(a, b, X, 1e-2, 200);  // returns 4; compare [0..2]
  CHECK(lie.Lf[0] == Approx(orc[0]).epsilon(1e-9).margin(1e-9));
  CHECK(lie.Lf[1] == Approx(orc[1]).epsilon(1e-3).margin(1e-6));
  CHECK(lie.Lf[2] == Approx(orc[2]).epsilon(1e-3).margin(1e-6));
}

// Safety guarantee: with the impact barrier hard and a state just inside the
// safe set, an unsafe nominal (nose-down, which drives the predicted load up) is
// corrected so the assembled degree-2 impact row a . U <= rhs holds.
TEST_CASE("lon filter enforces the impact barrier when hard", "[lon_cbf]") {
  Setup s;
  const double V = 16.0, gamma = -5.0 * kDeg, theta = 4.0 * kDeg;
  AeroLocal a = makeAeroLocal(s.table, s.mx, s.cfg, V, theta - gamma);
  LonStateVec X;
  X << 3.0, V, gamma, theta, 0.0, 1.5, 0.0;  // low altitude, steady descent

  // Measure n_peak here (Nb=0, n_limit=0 => b = -n_peak), then set n_limit just
  // above it so the state is barely inside the safe set and the row is active.
  const ImpactLoadBarrier probe =
      makeImpactLoadBarrier(a, 0.0, 22.5 * kDeg, 1000.0, 0.0, 2.0, 0.0, 0.02, X);
  const double n_peak = -probe(toArray(X));

  LonCBFConfig cfg;
  cfg.stall = false; cfg.noseup = false; cfg.energy = false;
  cfg.thrust_limits = false; cfg.impact = true; cfg.impact_hard = true;
  cfg.Nb = 0.0;                 // no height relaxation -> barrier fully active
  cfg.n_limit = n_peak + 1.0;   // ~1 g margin inside the safe set
  cfg.z_gate = 50.0;
  LonCBFFilter filter(cfg);

  // Assemble the same (hard) impact row the filter will, then drive the nominal
  // elevator hard onto the violating side of a . U <= rhs.
  const ImpactLoadBarrier b = makeImpactLoadBarrier(
      a, cfg.n_limit, cfg.beta, cfg.rho_water, cfg.Nb, cfg.zs, cfg.tau_keel, cfg.eps_g0, X);
  auto lie = barrierLie<2>(a, b, X);
  std::vector<double> Lf(lie.Lf.begin(), lie.Lf.end());
  HocbfRow row = hocbfRow(Lf, lie.LgLf, {cfg.c_impact[0], cfg.c_impact[1]});
  const double a_de = row.a(0, LDE);
  REQUIRE(std::abs(a_de) > 1e-6);
  CHECK(std::abs(row.a(0, LTDDOT)) < 1e-9);  // thrust absent from the degree-2 row

  LonCtrlVec Un;
  Un << (a_de > 0 ? cfg.de_max : cfg.de_min), 0.0;  // most-violating elevator in-box
  REQUIRE(a_de * Un[LDE] > row.rhs + 1e-6);          // it really violates the row

  LonCtrlVec U = filter.filter(Un, X, s.table, s.mx, s.cfg);
  CHECK_FALSE(filter.lastRecovery());
  CHECK((U - Un).norm() > 1e-6);  // the unsafe command was corrected
  const double lhs = row.a(0, LDE) * U[LDE] + row.a(0, LTDDOT) * U[LTDDOT];
  CHECK(lhs <= row.rhs + 1e-6);   // assembled impact row now holds
}
