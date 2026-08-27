#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <Eigen/Eigenvalues>
#include <cmath>
#include <fstream>
#include <string>

#include "autoland/beaver_dynamics.hpp"
#include "autoland/sixdof_sim.hpp"

using namespace autoland;
using Catch::Approx;

namespace {
constexpr double kDeg = M_PI / 180.0;

// The FDC 1.2 manual's ACTRIM check case (figs. 10.18/10.19): the complete
// trimmed state, inputs, and state derivative printed for V=35 m/s, n=1800
// RPM, pz=20 "Hg. The printed trim state carries H = 0 and its accelerations
// are only reproduced with the SEA-LEVEL atmosphere (the 2000 ft prompt seeds
// the simulation IC, not the trim atmosphere) -- verified empirically: at
// rho(2000 ft) the residuals are ~0.1 m/s^2, at rho(0) they match to ~3e-6.
struct FdcCase {
  double V = 35.0, alpha = 2.1131e-1, beta = -2.0667e-2, theta = 1.9190e-1;
  double de = -9.3083e-2, da = 9.6242e-3, dr = -4.9506e-2;
  double pz = 20.0, n = 1800.0;
  double Vdot = -1.8871e-4, alphadot = -1.2348e-5, betadot = 4.6356e-4;
  double pdot = -2.5027e-5, qdot = -2.0660e-5, rdot = -5.2604e-5;
  double yedot = -7.2330e-1, Hdot = -6.7909e-1;
};

StateVec caseState(const FdcCase& c) {
  StateVec x = StateVec::Zero();
  x[U] = c.V * std::cos(c.alpha) * std::cos(c.beta);
  x[V] = c.V * std::sin(c.beta);
  x[W] = c.V * std::sin(c.alpha) * std::cos(c.beta);
  x[THETA] = c.theta;
  return x;
}

CtrlVec caseCtrl(const FdcCase& c, const BeaverDynamics& dyn) {
  CtrlVec u = CtrlVec::Zero();
  u[DE] = c.de;
  u[DA] = c.da;
  u[DR] = c.dr;
  u[DT] = (c.pz - dyn.config().pz_idle) /
          (dyn.config().pz_max - dyn.config().pz_idle);
  return u;
}
}  // namespace

// External oracle: our 6-DOF xdot evaluated at the FDC-printed trim point must
// reproduce the FDC-printed xdot. This validates the whole nonlinear chain --
// aero polynomials, engine model, atmosphere, gravity, EOM assembly, and
// kinematics -- against an independent published implementation.
TEST_CASE("Beaver 6-DOF xdot reproduces the FDC ACTRIM check case",
          "[beaver6]") {
  FdcCase c;
  BeaverPlantConfig pc;
  pc.n_rpm = c.n;
  pc.h_ref = 0.0;
  BeaverDynamics dyn(pc);
  REQUIRE(dyn.rho() == Approx(1.225).epsilon(1e-4));

  const StateVec x = caseState(c);
  const StateVec xd = dyn.xdot(x, caseCtrl(c, dyn));

  // Convert body-axis rates to the FDC (Vdot, alphadot, betadot) coordinates.
  const double u = x[U], v = x[V], w = x[W];
  const double Vt = std::sqrt(u * u + v * v + w * w);
  const double Vdot = (u * xd[U] + v * xd[V] + w * xd[W]) / Vt;
  const double alphadot = (u * xd[W] - w * xd[U]) / (u * u + w * w);
  const double betadot =
      (xd[V] * Vt - v * Vdot) / (Vt * std::sqrt(u * u + w * w));

  // The FDC printout has 5 significant digits; 1e-5 absolute on the
  // accelerations is print precision (forces are O(g), so this is ~1e-6 rel).
  CHECK(Vdot == Approx(c.Vdot).margin(1e-5));
  CHECK(alphadot == Approx(c.alphadot).margin(1e-5));
  CHECK(betadot == Approx(c.betadot).margin(1e-5));
  CHECK(xd[P] == Approx(c.pdot).margin(1e-5));
  CHECK(xd[Q] == Approx(c.qdot).margin(1e-5));
  CHECK(xd[R] == Approx(c.rdot).margin(1e-5));
  // Kinematic rows to the printout's own precision.
  CHECK(xd[Y] == Approx(c.yedot).margin(1e-4));
  CHECK(xd[H] == Approx(c.Hdot).margin(1e-3));
}

// Our 6-axis Newton trim at the FDC descent angle must recover the printed
// longitudinal trim to ~print precision. (The lateral unknowns beta/da/dr are
// compared loosely: FDC's fmins stopped with a beta-dot residual of 4.6e-4,
// so its printed lateral trim is only converged to ~0.2 deg.)
TEST_CASE("Beaver trim recovers the FDC check-case trim", "[beaver6]") {
  FdcCase c;
  BeaverPlantConfig pc;
  pc.n_rpm = c.n;
  pc.h_ref = 0.0;
  BeaverDynamics dyn(pc);

  const double gamma = std::asin(c.Hdot / c.V);
  const TrimResult t = beaverTrim(dyn, c.V, gamma);
  REQUIRE(t.converged);
  CHECK(t.residual < 1e-9);

  CHECK(t.alpha == Approx(c.alpha).margin(2e-4));               // 0.01 deg
  CHECK(t.theta == Approx(c.theta).margin(2e-4));
  CHECK(t.u[DE] == Approx(c.de).margin(2e-4));
  CHECK(dyn.pzFromThrottle(t.u[DT]) == Approx(c.pz).margin(0.01));
  CHECK(std::asin(t.x[V] / c.V) == Approx(c.beta).margin(4e-3));  // ~0.2 deg
  CHECK(t.u[DA] == Approx(c.da).margin(4e-3));
  CHECK(t.u[DR] == Approx(c.dr).margin(4e-3));
}

// The exact autodiff linearization at cruise must show the classic Beaver
// modes with sane frequencies/damping: stable short period, lightly-damped
// stable phugoid, fast stable roll, damped Dutch roll, near-neutral spiral.
TEST_CASE("Beaver linearized modes are the classic set at cruise",
          "[beaver6]") {
  BeaverPlantConfig pc;
  pc.h_ref = 1828.8;  // the LR-556 aero-basis altitude (Table C.2)
  BeaverDynamics dyn(pc);
  const TrimResult t = beaverTrim(dyn, 45.0, 0.0);
  REQUIRE(t.converged);

  Mat A, B;
  dyn.linearize(t.x, t.u, A, B);

  auto sub = [&](std::array<int, 4> idx) {
    Eigen::Matrix4d S;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) S(i, j) = A(idx[i], idx[j]);
    return Eigen::EigenSolver<Eigen::Matrix4d>(S).eigenvalues();
  };

  const Eigen::Vector4cd lon = sub({U, W, Q, THETA});
  double wn_sp = 0, z_sp = 0, wn_ph = 0, z_ph = 0;
  for (int i = 0; i < 4; ++i) {
    if (lon[i].imag() <= 0) continue;  // one of each conjugate pair
    const double wn = std::abs(lon[i]);
    const double z = -lon[i].real() / wn;
    if (wn > 1.0) { wn_sp = wn; z_sp = z; }
    else { wn_ph = wn; z_ph = z; }
  }
  CHECK(wn_sp == Approx(3.2).margin(0.6));   // short period ~3.2 rad/s
  CHECK(z_sp == Approx(0.66).margin(0.15));
  CHECK(wn_ph == Approx(0.26).margin(0.08)); // phugoid ~0.26 rad/s
  CHECK(z_ph > 0.0);
  CHECK(z_ph < 0.25);

  const Eigen::Vector4cd lat = sub({V, P, R, PHI});
  double roll = 0, spiral = -1, wn_dr = 0, z_dr = 0;
  for (int i = 0; i < 4; ++i) {
    if (lat[i].imag() > 0) {
      wn_dr = std::abs(lat[i]);
      z_dr = -lat[i].real() / wn_dr;
    } else if (lat[i].imag() == 0) {
      if (std::abs(lat[i].real()) > 1.0) roll = lat[i].real();
      else spiral = lat[i].real();
    }
  }
  CHECK(roll < -3.0);                        // roll subsidence, tau < 0.33 s
  CHECK(wn_dr == Approx(1.1).margin(0.3));   // Dutch roll ~1.1 rad/s
  CHECK(z_dr == Approx(0.44).margin(0.15));
  CHECK(std::abs(spiral) < 0.15);            // near-neutral (here: stable)
}

// A steady coordinated 30-deg-bank level turn at 40 m/s IS an equilibrium of
// the plant. The trim state below was solved to ~1e-12 on the INDEPENDENT
// Python implementation (scripts/validate_beaver_sixdof.py turn_trim), so
// this is a cross-implementation check at NONZERO body rates -- the one
// regime the wings-level oracles cannot see: it exercises the quadratic
// gyroscopic terms ((Izz-Iyy)qr, Ixz(p^2-r^2), ...), the omega x v Coriolis
// terms, and the banked-attitude kinematics. The physics anchor: n_z matches
// 1/cos(phi) and psidot = g tan(phi)/V by construction of the solve.
TEST_CASE("Steady 30-deg coordinated turn is an equilibrium (gyroscopic terms)",
          "[beaver6]") {
  const double V0 = 40.0, phi = 30.0 * kDeg;
  const double alpha = 0.180149730674244, beta = -0.022711703564002;
  const double theta = 0.145151507839705;
  const double de = -0.099018121327235, da = 0.058757054469476,
               dr = -0.073000017611803, pz = 24.563634304513656;
  const double p = -0.020473661276978, q = 0.070029147336906,
               r = 0.121294041198248;

  BeaverDynamics dyn;  // sea level, n=1800, flaps up
  StateVec x = StateVec::Zero();
  x[U] = V0 * std::cos(alpha) * std::cos(beta);
  x[V] = V0 * std::sin(beta);
  x[W] = V0 * std::sin(alpha) * std::cos(beta);
  x[P] = p;
  x[Q] = q;
  x[R] = r;
  x[PHI] = phi;
  x[THETA] = theta;
  CtrlVec u = CtrlVec::Zero();
  u[DE] = de;
  u[DA] = da;
  u[DR] = dr;
  u[DT] = (pz - dyn.config().pz_idle) /
          (dyn.config().pz_max - dyn.config().pz_idle);

  const StateVec xd = dyn.xdot(x, u);
  for (int i : {U, V, W, P, Q, R})
    CHECK(std::abs(xd[i]) < 1e-8);  // dynamic rows vanish
  // Euler kinematics of the steady turn: phidot = thetadot = 0, psidot = Om.
  const double Om = 9.80665 * std::tan(phi) / V0;
  CHECK(std::abs(xd[PHI]) < 1e-9);
  CHECK(std::abs(xd[THETA]) < 1e-9);
  CHECK(xd[PSI] == Approx(Om).epsilon(1e-9));
  CHECK(std::abs(xd[H]) < 1e-8);  // level turn
}

// Wind enters only through the aerodynamics: zero wind is bit-identical, the
// kinematic rows are wind-invariant, and the physical gust signs are right.
TEST_CASE("Beaver wind coupling: zero-wind identity and aero signs",
          "[beaver6]") {
  BeaverDynamics dyn;
  const TrimResult t = beaverTrim(dyn, 40.0, 0.0);
  const StateVec x = t.x;
  const CtrlVec u = t.u;

  const StateVec still = dyn.xdot(x, u);
  const StateVec still3 = dyn.xdot(x, u, Eigen::Vector3d::Zero());
  for (int j = 0; j < NX; ++j) CHECK(still[j] == still3[j]);

  const StateVec updraft = dyn.xdot(x, u, Eigen::Vector3d(0, 0, 2.0));
  CHECK(updraft[W] < still[W]);  // more lift (z-down)
  const StateVec tail = dyn.xdot(x, u, Eigen::Vector3d(5.0, 0, 0));
  CHECK(tail[W] > still[W]);     // less airspeed -> less lift
  for (int j : {PHI, THETA, PSI, H, Y}) CHECK(tail[j] == still[j]);
}

// Closed-loop calm straight-in on the Beaver plant (the default scenario):
// the loop must fly the glideslope to touchdown at the trim sink with the
// lateral axes quiet, all inside the LR-556 validity band.
TEST_CASE("Beaver calm 6-DOF straight-in lands on the glideslope",
          "[beaver6][sim]") {
  SixDofSim sim("", "", std::string(AUTOLAND_DATA_DIR) +
                            "/beaver_landing_calm.yaml");
  REQUIRE(sim.trimResult().converged);
  const SixDofTouchdown td = sim.run("test_beaver_calm.csv");

  REQUIRE(td.reached);
  CHECK(td.sink == Approx(40.0 * std::sin(3.5 * kDeg)).margin(0.15));
  CHECK(td.V == Approx(40.0).margin(0.3));
  CHECK(std::abs(td.gamma + 3.5 * kDeg) < 0.3 * kDeg);
  CHECK(std::abs(td.phi) < 1.0 * kDeg);
  CHECK(std::abs(td.y) < 1.0);
  CHECK(sim.stats().max_alpha < 12.0 * kDeg);  // stall margin held
  CHECK(sim.stats().min_V > 35.0);             // inside the validity band
}
