#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "autoland/aero_table.hpp"
#include "autoland/config.hpp"
#include "autoland/hocbf.hpp"
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

// Control authorities (the coefficients of U in the 3rd derivative) must match
// the closed forms in documentation/water_landing_cbf_math.md sections 3.1-3.2.
//
// NOTE: the doc's ELEVATOR authority models the aero as L(alpha,V), D(alpha,V)
// only -- it omits the aerodynamic pitch-RATE dependence (dCL/dqhat, dCD/dqhat).
// Our autodiff differentiates the full tabulated aero, so it matches the doc
// exactly only in the no-rate-aero limit (dQ = 0). The THRUST authority has no
// aero term and matches with full aero. The separate test below records that
// full aero diverges from the doc (the engine is the more complete quantity).
TEST_CASE("autodiff control authorities match the doc closed forms", "[lon_cbf]") {
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
  const double dCL_da = -a.dAlpha_CFx * std::sin(alpha) - CFx * std::cos(alpha) +
                        a.dAlpha_CFz * std::cos(alpha) - CFz * std::sin(alpha);
  const double dCD_da = a.dAlpha_CFx * std::cos(alpha) - CFx * std::sin(alpha) +
                        a.dAlpha_CFz * std::sin(alpha) + CFz * std::cos(alpha);
  const double dL_da = qbar * Sref * dCL_da;
  const double dD_da = qbar * Sref * dCD_da;
  const double Cmde_factor = qbar * Sref * cref / Iyy;  // (rho V^2 S cbar / 2 Iyy)
  const double Cmde = a.dDe_CMy;

  LonStateVec X;
  X << 30.0, V, gamma, theta, q, T, 0.0;

  SECTION("descent barrier") {
    DescentBarrier b{0.6 * 0.6, 3.0};
    auto lie = barrierLie<3>(a, b, X);
    // Thrust authority = sin(theta)/m
    CHECK(lie.LgLf(0, LTDDOT) == Approx(std::sin(theta) / m).epsilon(1e-7));
    // Elevator authority = (1/m)(T cos th + dL/da cos g - dD/da sin g) * factor * Cmde
    const double elev = (1.0 / m) *
                        (T * std::cos(theta) + dL_da * std::cos(gamma) - dD_da * std::sin(gamma)) *
                        Cmde_factor * Cmde;
    CHECK(lie.LgLf(0, LDE) == Approx(elev).epsilon(1e-7));
  }

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
  DescentBarrier bd{cfg.v_safe * cfg.v_safe, cfg.a_brk};
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
