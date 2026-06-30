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

// The augmented EOM's first Lie derivatives must reproduce the raw drift: e.g.
// L_f b_V = Vdot (the airspeed equation), since b_V = V - Vmin.
TEST_CASE("L_f b_V equals the airspeed equation", "[lon_cbf]") {
  Setup s;
  const double V = 18.0, gamma = -5.0 * kDeg, theta = -2.0 * kDeg;
  const double alpha = theta - gamma;
  AeroLocal a = makeAeroLocal(s.table, s.mx, s.cfg, V, alpha);

  LonStateVec X;
  X << 30.0, V, gamma, theta, 0.0, 5.0, 0.0;
  AirspeedBarrier bv{0.85 * V};
  auto lie = barrierLie<3>(a, bv, X);

  const LonStateVec xd = lonXdot(a, X, LonCtrlVec::Zero());
  CHECK(lie.Lf[1] == Approx(xd[LV]).epsilon(1e-9));  // L_f b_V == Vdot
  CHECK(lie.Lf[0] == Approx(V - 0.85 * V).epsilon(1e-12));
}

// The descent barrier now uses the speed/path-angle-dependent braking accel
//   a_brk(V,gamma) = (rho V^2 S / 2m)[CLmax cos g - CDmaxlift sin g] - g,
// built by makeDescentBarrier(). Check the barrier value matches that formula
// exactly (catches a constant-a_brk regression or a sign/term error) and that
// a_brk is positive (and not the old constant 3.0) at a normal approach state.
TEST_CASE("descent barrier uses speed-dependent a_brk(V,gamma)", "[lon_cbf]") {
  Setup s;
  const double V = 18.0, gamma = -5.0 * kDeg, theta = -2.0 * kDeg;
  const double alpha = theta - gamma;
  const double CLmax = 1.2, v_safe = 0.6;
  AeroLocal a = makeAeroLocal(s.table, s.mx, s.cfg, V, alpha);

  LonStateVec X;
  X << 30.0, V, gamma, theta, 0.0, 5.0, 0.0;
  const DescentBarrier b = makeDescentBarrier(a, v_safe, CLmax, V);

  const double CDml = cdAtMaxLift(a, CLmax, V);
  const double abrk = (0.5 * a.rho * a.Sref / a.mass) * V * V *
                          (CLmax * std::cos(gamma) - CDml * std::sin(gamma)) - a.g;
  // The barrier floors the radicand with the smooth positive map eps_r (so it can
  // never go NaN); mirror it here. At this safe state (h=30) the floor is a ~1e-9
  // bias, but reproducing it keeps the exact-match tolerance meaningful.
  const double arg = v_safe * v_safe + 2.0 * abrk * X[LH];
  const double arg_pos =
      0.5 * (arg + std::sqrt(arg * arg + 4.0 * b.eps_r * b.eps_r));
  const double b_expected = V * std::sin(gamma) + std::sqrt(arg_pos);

  CHECK(b(toArray(X)) == Approx(b_expected).epsilon(1e-12));
  CHECK(abrk > 0.0);                            // real braking authority
  CHECK(abrk != Approx(3.0).epsilon(1e-3));     // tracks dynamic pressure, not const
}

// Airspeed control authority (the coefficients of U in the 3rd derivative) must
// match the closed form in documentation/water_landing_cbf_math.md section 3.2.
// (The descent barrier now uses the speed-dependent a_brk(V,gamma), which has no
// simple closed form; it is covered by the a_brk-formula + flow-oracle tests.)
//
// NOTE: the doc's ELEVATOR authority models the aero as L(alpha,V), D(alpha,V)
// only -- it omits the aerodynamic pitch-RATE dependence (dCL/dqhat, dCD/dqhat).
// Our autodiff differentiates the full tabulated aero, so it matches the doc
// exactly only in the no-rate-aero limit (dQ = 0). The THRUST authority has no
// aero term and matches with full aero. The separate test below records that
// full aero diverges from the doc (the engine is the more complete quantity).
TEST_CASE("airspeed control authority matches the doc closed form", "[lon_cbf]") {
  Setup s;
  const double V = 18.0, gamma = -5.0 * kDeg, theta = -2.0 * kDeg;
  const double alpha = theta - gamma;
  const double T = 5.0, q = 0.3;  // nonzero q to exercise the explicit-q terms
  AeroLocal a = makeAeroLocal(s.table, s.mx, s.cfg, V, alpha);
  // No-rate-aero idealization matching the doc's L(alpha,V), D(alpha,V) model.
  a.dQ_CFx = 0.0; a.dQ_CFz = 0.0; a.dQ_CMy = 0.0;

  const double m = s.cfg.inertia.mass, Iyy = s.cfg.inertia.Iyy;
  const double Sref = s.table.Sref(), cref = s.table.Cref();
  const double qbar = 0.5 * s.cfg.env.rho * V * V;
  const double mach = V / s.cfg.env.a_sound, qhat = q * cref / (2.0 * V);

  // File-frame body coeffs and rotated drag slope at this state.
  const double CFx = a.off_CFx + a.dAlpha_CFx * alpha + a.dMach_CFx * mach + a.dQ_CFx * qhat;
  const double CFz = a.off_CFz + a.dAlpha_CFz * alpha + a.dMach_CFz * mach + a.dQ_CFz * qhat;
  const double dCD_da = a.dAlpha_CFx * std::cos(alpha) - CFx * std::sin(alpha) +
                        a.dAlpha_CFz * std::sin(alpha) + CFz * std::cos(alpha);
  const double dD_da = qbar * Sref * dCD_da;
  const double Cmde_factor = qbar * Sref * cref / Iyy;  // (rho V^2 S cbar / 2 Iyy)
  const double Cmde = a.dDe_CMy;

  LonStateVec X;
  X << 30.0, V, gamma, theta, q, T, 0.0;

  SECTION("airspeed barrier") {
    AirspeedBarrier b{0.85 * V};
    auto lie = barrierLie<3>(a, b, X);
    // Thrust authority = cos(alpha)/m
    CHECK(lie.LgLf(0, LTDDOT) == Approx(std::cos(alpha) / m).epsilon(1e-7));
    // Elevator authority = -(1/m)(T sin a + dD/da) * factor * Cmde
    const double elev = -(1.0 / m) * (T * std::sin(alpha) + dD_da) * Cmde_factor * Cmde;
    CHECK(lie.LgLf(0, LDE) == Approx(elev).epsilon(1e-7));
  }
}

// Records the finding: with the FULL tabulated aero (pitch-rate derivatives
// retained), the elevator authority differs from the doc's reduced closed form,
// while the thrust authority still matches exactly.
TEST_CASE("full aero elevator authority diverges from the doc", "[lon_cbf]") {
  Setup s;
  const double V = 18.0, gamma = -5.0 * kDeg, theta = -2.0 * kDeg;
  const double alpha = theta - gamma;
  AeroLocal a = makeAeroLocal(s.table, s.mx, s.cfg, V, alpha);
  REQUIRE(a.dQ_CMy != 0.0);  // the model really does have pitch-rate aero

  LonStateVec X;
  X << 30.0, V, gamma, theta, 0.3, 5.0, 0.0;
  AirspeedBarrier b{0.85 * V};
  auto lie = barrierLie<3>(a, b, X);

  // Thrust authority is aero-free -> still exactly cos(alpha)/m.
  CHECK(lie.LgLf(0, LTDDOT) == Approx(std::cos(alpha) / s.cfg.inertia.mass).epsilon(1e-7));

  // Elevator authority vs the doc's reduced form: a real, non-tiny gap.
  AeroLocal a0 = a; a0.dQ_CFx = 0; a0.dQ_CFz = 0; a0.dQ_CMy = 0;
  auto lie0 = barrierLie<3>(a0, b, X);
  CHECK(lie.LgLf(0, LDE) != Approx(lie0.LgLf(0, LDE)).epsilon(1e-3));
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

// The safety guarantee: after filtering, the hard barrier rows hold even when
// the nominal command is unsafe (here a strong nose-down into the surface).
TEST_CASE("lon filter enforces the hard barrier rows", "[lon_cbf]") {
  Setup s;
  LonCBFConfig cfg;             // descent + thrust hard; airspeed soft
  cfg.Vmin = 13.0;
  cfg.Vmin_ground = 13.0;       // constant floor (collapse the altitude schedule)
  LonCBFFilter filter(cfg);

  LonStateVec X;
  X << 2.0, 16.0, -6.0 * kDeg, -2.0 * kDeg, 0.0, 1.5, 0.0;  // low + steep descent
  LonCtrlVec Un;
  Un << -0.3, -100.0;           // unsafe: nose-down + chop thrust
  LonCtrlVec U = filter.filter(Un, X, s.table, s.mx, s.cfg);

  CHECK_FALSE(filter.lastRecovery());
  CHECK((U - Un).norm() > 1e-6);  // the unsafe command was corrected

  const AeroLocal a = makeAeroLocal(s.table, s.mx, s.cfg, X[LV], X[LTH] - X[LGAM]);
  // Descent HOCBF row (hard) must hold: a . U <= rhs.
  DescentBarrier bd = makeDescentBarrier(a, cfg.v_safe, cfg.CLmax, X[LV]);
  auto lie = barrierLie<3>(a, bd, X);
  std::vector<double> Lf(lie.Lf.begin(), lie.Lf.end());
  HocbfRow row = hocbfRow(Lf, lie.LgLf,
                          {cfg.c_descent[0], cfg.c_descent[1], cfg.c_descent[2]});
  const double lhs = row.a(0, LDE) * U[LDE] + row.a(0, LTDDOT) * U[LTDDOT];
  CHECK(lhs <= row.rhs + 1e-6);

  // Min-thrust row (hard) must hold: -Tddot <= (c11+c12)Tdot + c11 c12 T.
  HocbfRow tmin = thrustMinRow(X, cfg.c_thrust_min[0], cfg.c_thrust_min[1]);
  CHECK(tmin.a(0, LTDDOT) * U[LTDDOT] <= tmin.rhs + 1e-6);
}

// Independently cross-check the DRIFT Lie stack {b, L_f b, L_f^2 b, L_f^3 b}
// against the finite-difference flow oracle. The existing tests verify the
// control row L_g L_f^2 b vs the doc, but never the drift terms (the doc never
// writes L_f^3 b), so a bug in the bespoke Taylor-jet engine (lie_taylor.hpp)
// could slip through. The state carries nonzero q, T, Tdot so the higher-order
// terms genuinely exercise the thrust/pitch couplings.
TEST_CASE("drift Lie stack matches a finite-difference flow oracle", "[lon_cbf]") {
  Setup s;
  const double V = 18.0, gamma = -5.0 * kDeg, theta = -2.0 * kDeg;
  const double alpha = theta - gamma;
  const double q = 0.3, T = 5.0, Tdot = 1.0;
  AeroLocal a = makeAeroLocal(s.table, s.mx, s.cfg, V, alpha);

  LonStateVec X;
  X << 30.0, V, gamma, theta, q, T, Tdot;  // h=30 keeps the descent sqrt safe

  const double h = 1e-2;
  const int nsub = 200;

  // L_f^0 b is exact (just b(X0)). For k>=1 the accuracy floor is the
  // finite-difference ORACLE, not the engine: even Richardson-extrapolated time
  // stencils bottom out around ~1e-4 relative on L_f^3 b. 1e-3 is therefore a
  // deliberately loose guard -- it still catches any structural engine bug
  // (a sign flip, a factor of 2, or a missing term are all >= several percent),
  // which is the point of this cross-check.
  auto checkStack = [](const std::array<double, 4>& engine,
                       const std::array<double, 4>& oracle) {
    CHECK(engine[0] == Approx(oracle[0]).epsilon(1e-9).margin(1e-9));
    CHECK(engine[1] == Approx(oracle[1]).epsilon(1e-3).margin(1e-6));
    CHECK(engine[2] == Approx(oracle[2]).epsilon(1e-3).margin(1e-6));
    CHECK(engine[3] == Approx(oracle[3]).epsilon(1e-3).margin(1e-6));
  };

  SECTION("descent barrier") {
    const DescentBarrier b = makeDescentBarrier(a, 0.6, 1.2, V);
    const auto lie = barrierLie<3>(a, b, X);
    checkStack(lie.Lf, lieDriftOracle(a, b, X, h, nsub));
  }

  SECTION("airspeed barrier") {
    AirspeedBarrier b{0.85 * V};
    const auto lie = barrierLie<3>(a, b, X);
    checkStack(lie.Lf, lieDriftOracle(a, b, X, h, nsub));
  }

  SECTION("upper airspeed barrier") {
    AirspeedUpperBarrier b{1.5 * V};
    const auto lie = barrierLie<3>(a, b, X);
    checkStack(lie.Lf, lieDriftOracle(a, b, X, h, nsub));
  }

  SECTION("energy barrier") {
    // The key correctness gate for the energy-reachability barrier: the series/
    // series Dmax/|sin gamma| divide and the gamma/h coupling must reproduce the
    // finite-difference drift stack (catches any relative-degree-3 surprise).
    const EnergyBarrier b = makeEnergyBarrier(a, 12.0, 1.2, V, 0.02);
    const auto lie = barrierLie<3>(a, b, X);
    checkStack(lie.Lf, lieDriftOracle(a, b, X, h, nsub));
  }
}

// The upper airspeed barrier b = Vmax - V is the sign-flip of the lower
// b_V = V - Vmin: same state dependence (V only), opposite sign. So for k >= 1
// the drift Lie terms and the control row negate the lower barrier's, while the
// k=0 value is Vmax - V (not V - Vmin). This is what lets it reuse the same
// degree-3 machinery (the HOCBF builder handles the sign through L_g L_f^2 b).
TEST_CASE("upper airspeed barrier is the sign-flip of the lower", "[lon_cbf]") {
  Setup s;
  const double V = 18.0, gamma = -5.0 * kDeg, theta = -2.0 * kDeg;
  const double alpha = theta - gamma;
  const double Vmin = 0.85 * V, Vmax = 1.5 * V;
  AeroLocal a = makeAeroLocal(s.table, s.mx, s.cfg, V, alpha);

  LonStateVec X;
  X << 30.0, V, gamma, theta, 0.3, 5.0, 1.0;  // nonzero q, T, Tdot
  AirspeedBarrier blo{Vmin};
  AirspeedUpperBarrier bup{Vmax};
  auto lielo = barrierLie<3>(a, blo, X);
  auto lieup = barrierLie<3>(a, bup, X);

  // k=0 values are the literal barrier definitions.
  CHECK(bup(toArray(X)) == Approx(Vmax - V).epsilon(1e-12));
  const LonStateVec xd = lonXdot(a, X, LonCtrlVec::Zero());
  CHECK(lieup.Lf[1] == Approx(-xd[LV]).epsilon(1e-9));  // L_f b_up == -Vdot

  // k>=1 drift terms and the control row negate the lower barrier's.
  for (int k = 1; k <= 3; ++k)
    CHECK(lieup.Lf[k] == Approx(-lielo.Lf[k]).epsilon(1e-9).margin(1e-12));
  CHECK(lieup.LgLf(0, LDE) == Approx(-lielo.LgLf(0, LDE)).epsilon(1e-9).margin(1e-12));
  CHECK(lieup.LgLf(0, LTDDOT) == Approx(-lielo.LgLf(0, LTDDOT)).epsilon(1e-9).margin(1e-12));
}

// The safety guarantee for the upper barrier: with it hard and a low Vmax, an
// over-speed nominal (nose-down + max thrust to accelerate) must be corrected so
// the assembled upper airspeed HOCBF row a . U <= rhs holds.
TEST_CASE("lon filter enforces the upper airspeed barrier when hard", "[lon_cbf]") {
  Setup s;
  LonCBFConfig cfg;
  cfg.descent = false;          // isolate the upper airspeed barrier
  cfg.airspeed = false;
  cfg.thrust_limits = false;
  cfg.airspeed_upper = true;
  cfg.airspeed_upper_hard = true;
  cfg.Vmax_air = 18.5;          // just above the state's V -> barrier is active
  cfg.Vmax_ground = 18.5;       // constant ceiling (collapse the altitude schedule)
  LonCBFFilter filter(cfg);

  LonStateVec X;
  X << 30.0, 18.0, -3.0 * kDeg, -1.0 * kDeg, 0.0, 5.0, 0.0;  // near Vmax
  LonCtrlVec Un;
  Un << -0.3, 500.0;            // unsafe: nose-down + slam thrust (accelerate)
  LonCtrlVec U = filter.filter(Un, X, s.table, s.mx, s.cfg);

  CHECK_FALSE(filter.lastRecovery());
  CHECK((U - Un).norm() > 1e-6);  // the over-speed command was corrected

  const AeroLocal a = makeAeroLocal(s.table, s.mx, s.cfg, X[LV], X[LTH] - X[LGAM]);
  AirspeedUpperBarrier b{cfg.Vmax_air};
  auto lie = barrierLie<3>(a, b, X);
  std::vector<double> Lf(lie.Lf.begin(), lie.Lf.end());
  HocbfRow row = hocbfRow(Lf, lie.LgLf,
                          {cfg.c_airspeed_upper[0], cfg.c_airspeed_upper[1],
                           cfg.c_airspeed_upper[2]});
  const double lhs = row.a(0, LDE) * U[LDE] + row.a(0, LTDDOT) * U[LTDDOT];
  CHECK(lhs <= row.rhs + 1e-6);
}

// =============================================================================
// Energy-reachability upper (anti-bounce) barrier.
// =============================================================================

// b = [E_land + (Dmax/|sin gamma|) h] - (1/2 m V^2 + m g h) matches its closed
// form aloft, and at h=0 reduces to the touchdown-speed bound b >= 0 <=> V <= V_land.
TEST_CASE("energy barrier value and touchdown-speed bound", "[lon_cbf]") {
  Setup s;
  const double V = 18.0, gamma = -5.0 * kDeg, theta = -2.0 * kDeg;
  const double alpha = theta - gamma;
  const double V_land = 12.0, eps_sg = 0.02, CLmax = 1.2;
  AeroLocal a = makeAeroLocal(s.table, s.mx, s.cfg, V, alpha);
  const EnergyBarrier b = makeEnergyBarrier(a, V_land, CLmax, V, eps_sg);

  // (1) Value matches the closed form at an aloft state (h = 30). SPECIFIC energy:
  // e = 1/2 V^2 + g h, dmax = (rho Sref CDml / 2m) V^2, e_land = 1/2 V_land^2.
  LonStateVec X;
  X << 30.0, V, gamma, theta, 0.3, 5.0, 1.0;
  const double CDml = cdAtMaxLift(a, CLmax, V);
  const double dmax = 0.5 * a.rho * a.Sref * CDml / a.mass * V * V;
  const double sg = std::sqrt(std::sin(gamma) * std::sin(gamma) + eps_sg * eps_sg);
  const double e = 0.5 * V * V + a.g * 30.0;
  const double e_land = 0.5 * V_land * V_land;
  const double expect = (e_land + (dmax / sg) * 30.0) - e;
  CHECK(b(toArray(X)) == Approx(expect).epsilon(1e-12));

  // (2) At h = 0 the budget term vanishes -> b = 1/2 (V_land^2 - V^2):
  // zero at V = V_land, negative above, positive below.
  LonStateVec X0 = X; X0[LH] = 0.0;
  LonStateVec Xat = X0;   Xat[LV] = V_land;
  CHECK(b(toArray(Xat)) == Approx(0.0).margin(1e-9));
  LonStateVec Xhot = X0;  Xhot[LV] = V_land + 2.0;
  CHECK(b(toArray(Xhot)) < 0.0);
  LonStateVec Xslow = X0; Xslow[LV] = V_land - 2.0;
  CHECK(b(toArray(Xslow)) > 0.0);
}

// With the energy barrier hard and a low V_land, an energy-adding nominal
// (nose-down + max thrust) must be corrected so the assembled HOCBF row holds.
TEST_CASE("lon filter enforces the energy barrier when hard", "[lon_cbf]") {
  Setup s;
  LonCBFConfig cfg;
  cfg.descent = false;          // isolate the energy barrier
  cfg.airspeed = false;
  cfg.thrust_limits = false;
  cfg.airspeed_upper = false;
  cfg.impact = false;
  cfg.energy = true;
  cfg.energy_hard = true;
  cfg.V_land = 12.0;            // low target; the state is well above -> active
  LonCBFFilter filter(cfg);

  LonStateVec X;
  X << 3.0, 16.0, -3.0 * kDeg, -1.0 * kDeg, 0.0, 5.0, 0.0;  // low + fast -> hot
  LonCtrlVec Un;
  Un << -0.3, 500.0;            // unsafe: nose-down + slam thrust (add energy)
  LonCtrlVec U = filter.filter(Un, X, s.table, s.mx, s.cfg);

  CHECK_FALSE(filter.lastRecovery());
  CHECK((U - Un).norm() > 1e-6);  // the energy-adding command was corrected

  const AeroLocal a = makeAeroLocal(s.table, s.mx, s.cfg, X[LV], X[LTH] - X[LGAM]);
  const EnergyBarrier b = makeEnergyBarrier(a, cfg.V_land, cfg.CLmax, X[LV], cfg.eps_sg);
  auto lie = barrierLie<3>(a, b, X);
  std::vector<double> Lf(lie.Lf.begin(), lie.Lf.end());
  HocbfRow row = hocbfRow(Lf, lie.LgLf,
                          {cfg.c_energy[0], cfg.c_energy[1], cfg.c_energy[2]});
  const double lhs = row.a(0, LDE) * U[LDE] + row.a(0, LTDDOT) * U[LTDDOT];
  CHECK(lhs <= row.rhs + 1e-6);
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
// the new exp() Taylor overload (through Phi) along a descending flow.
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
  cfg.descent = false; cfg.airspeed = false; cfg.airspeed_upper = false;
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
