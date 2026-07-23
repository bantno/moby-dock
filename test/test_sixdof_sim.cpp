#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <fstream>
#include <string>

#include "autoland/config.hpp"
#include "autoland/dynamics.hpp"
#include "autoland/mixing.hpp"
#include "autoland/sixdof_sim.hpp"
#include "autoland/wind_gust.hpp"

using namespace autoland;
using Catch::Approx;

namespace {
constexpr double kDeg = M_PI / 180.0;

struct Setup {
  AeroTable table;
  AircraftConfig cfg;
  Mixing mx;
  Dynamics dyn;
  explicit Setup(const std::string& stab = "/example.stab")
      : table(AeroTable::fromFile(std::string(AUTOLAND_DATA_DIR) + stab)),
        cfg(loadAircraftConfig(std::string(AUTOLAND_DATA_DIR) +
                               "/aircraft.yaml")),
        mx(Mixing::build(cfg, table)),
        dyn(table, mx, cfg) {}
};

StateVec sampleState(int i) {
  StateVec x = StateVec::Zero();
  x[U] = 16.0 + i;
  x[V] = 0.4 * i - 0.6;
  x[W] = 1.0 + 0.3 * i;
  x[P] = 0.1 * i - 0.15;
  x[Q] = -0.05 * i;
  x[R] = 0.07 * i - 0.1;
  x[PHI] = 0.06 * i - 0.1;
  x[THETA] = 0.04 + 0.02 * i;
  x[PSI] = 0.3 * i - 0.5;
  x[H] = 30.0 + 5.0 * i;
  x[Y] = 2.0 * i - 3.0;
  return x;
}

// Minimal scenario writer for the closed-loop tests: calm base + optional
// wind/waves blocks appended verbatim.
std::string writeScenario(const std::string& name, const std::string& extra) {
  const std::string path = name;
  std::ofstream f(path);
  f << "V_app: 18.0\n"
       "gamma_app_deg: -3.0\n"
       "dt: 0.01\n"
       "t_max: 120.0\n"
       "initial:\n"
       "  h0: 40.0\n"
    << extra;
  return path;
}
}  // namespace

// The wind-aware overload with zero wind must reproduce xdot(x, u) EXACTLY
// (bit-identical), so every existing consumer of the 2-arg EOM is untouched.
TEST_CASE("Zero-wind xdot overload is bit-identical", "[sixdof]") {
  Setup s;
  CtrlVec u;
  u << 0.05, -0.03, 0.02, 0.4;
  for (int i = 0; i < 4; ++i) {
    const StateVec x = sampleState(i);
    const StateVec a = s.dyn.xdot(x, u);
    const StateVec b = s.dyn.xdot(x, u, Eigen::Vector3d::Zero());
    for (int j = 0; j < NX; ++j) CHECK(a[j] == b[j]);
  }
}

// Wind enters ONLY through the aerodynamics/thrust: the attitude and position
// kinematics rows depend on the inertial state alone and must be unchanged,
// while the physical responses carry the right signs (updraft -> more lift;
// lateral wind from the left -> pushed right + weathervane nose-left).
TEST_CASE("Wind coupling: kinematics untouched, aero signs correct",
          "[sixdof]") {
  Setup s;
  StateVec x = StateVec::Zero();
  x[U] = 18.0;  // level flight straight north, wings level
  x[W] = 1.0;
  x[H] = 30.0;
  CtrlVec u = CtrlVec::Zero();
  u[DT] = 0.3;

  const StateVec still = s.dyn.xdot(x, u);

  SECTION("kinematic rows are wind-invariant") {
    const StateVec wind = s.dyn.xdot(x, u, Eigen::Vector3d(3.0, 2.0, 1.0));
    for (int j : {PHI, THETA, PSI, H, Y}) CHECK(wind[j] == still[j]);
  }
  SECTION("updraft raises alpha -> more lift (z-down: wdot decreases)") {
    const StateVec wind = s.dyn.xdot(x, u, Eigen::Vector3d(0.0, 0.0, 2.0));
    CHECK(wind[W] < still[W]);
  }
  SECTION("tailwind cuts airspeed -> less lift (wdot increases)") {
    const StateVec wind = s.dyn.xdot(x, u, Eigen::Vector3d(5.0, 0.0, 0.0));
    CHECK(wind[W] > still[W]);
  }
  SECTION("crosswind toward +y: pushed downwind, weathervanes into wind") {
    // Air-relative v_air = -W_v < 0 => beta < 0 => side force +y (Cy_beta<0)
    // and yaw moment nose-left (Cn_beta>0).
    const StateVec wind = s.dyn.xdot(x, u, Eigen::Vector3d(0.0, 3.0, 0.0));
    CHECK(wind[V] > still[V]);
    CHECK(wind[R] < still[R]);
  }
}

// The new lateral gust axis follows the same 1-cosine shape as u/w and rides
// through gustWind; len <= 0 degenerates to a step (steady crosswind).
TEST_CASE("Lateral gust axis shape and step limit", "[sixdof][wind]") {
  DiscreteGustConfig g;
  g.enabled = true;
  g.v.amp = 3.0;
  g.v.len = 100.0;

  CHECK(gustWind(g, 0.0).v == 0.0);
  CHECK(gustWind(g, 50.0).v == Approx(1.5));   // half amplitude at midpoint
  CHECK(gustWind(g, 100.0).v == Approx(3.0));  // full amplitude at dm
  CHECK(gustWind(g, 500.0).v == Approx(3.0));  // holds beyond
  CHECK(gustWindRate(g, 50.0, 18.0).v ==
        Approx(0.5 * 3.0 * (M_PI / 100.0) * 18.0));  // peak slope * xdot

  g.v.len = 0.0;  // step gust = steady crosswind from onset
  CHECK(gustWind(g, 1e-3).v == Approx(3.0));
  CHECK(gustWind(g, 0.0).v == 0.0);
  CHECK(gustWind(g, 300.0).u == 0.0);  // other axes untouched
}

// Calm straight-in from trim: the closed loop must fly the glideslope to
// touchdown at the trim sink rate with the lateral axes quiet.
TEST_CASE("Calm 6-DOF straight-in reaches touchdown near trim sink",
          "[sixdof][sim]") {
  const std::string scenario = writeScenario("test_sixdof_calm.yaml", "");
  SixDofSim sim(std::string(AUTOLAND_DATA_DIR) + "/AHAB_combined_betasym.stab",
                std::string(AUTOLAND_DATA_DIR) + "/aircraft.yaml", scenario);
  REQUIRE(sim.trimResult().converged);
  const SixDofTouchdown td = sim.run("test_sixdof_calm.csv");

  REQUIRE(td.reached);
  // Nominal sink on the 18 m/s / -3 deg glideslope is ~0.94 m/s.
  CHECK(td.sink == Approx(18.0 * std::sin(3.0 * kDeg)).margin(0.3));
  CHECK(std::abs(td.gamma + 3.0 * kDeg) < 1.0 * kDeg);
  CHECK(sim.stats().max_abs_phi < 2.0 * kDeg);
  CHECK(sim.stats().max_abs_y < 0.5);
  CHECK(sim.stats().min_V > 15.0);
}

// Hot-and-high entry (fast + pitched up + off-heading): the speed axis must
// converge back to V_ref before touchdown. Regression for the speed runaway
// found during development: with the throttle railed at idle and (then) the
// CFx sign bug zeroing the deck's drag, a +2 m/s entry accelerated all the
// way to a ~29 m/s touchdown while gamma tracked its reference. The guards
// are the Kv_gamma speed->path reference shift (idle-rail robustness) and
// the corrected drag polar; this test keeps both honest.
TEST_CASE("Hot entry bleeds back to V_ref before touchdown",
          "[sixdof][sim]") {
  const std::string scenario = writeScenario("test_sixdof_hot.yaml",
                                             "  dV: 2.0\n"
                                             "  dtheta_deg: 2.0\n"
                                             "  dpsi_deg: 5.0\n");
  SixDofSim sim(std::string(AUTOLAND_DATA_DIR) + "/AHAB_combined_betasym.stab",
                std::string(AUTOLAND_DATA_DIR) + "/aircraft.yaml", scenario);
  const SixDofTouchdown td = sim.run("test_sixdof_hot.csv");

  REQUIRE(td.reached);
  CHECK(td.V == Approx(18.0).margin(0.5));
  CHECK(std::abs(td.gamma + 3.0 * kDeg) < 0.5 * kDeg);
  CHECK(td.sink < 1.3);
}

// Steady crosswind (step gust): the loop must still land, holding the
// centerline while the nose crabs into the wind.
TEST_CASE("Crosswind straight-in: bounded cross-track, crab into wind",
          "[sixdof][sim]") {
  const std::string scenario = writeScenario("test_sixdof_xwind.yaml",
                                             "wind:\n"
                                             "  enabled: true\n"
                                             "  t_start: 5.0\n"
                                             "  v_amp: 3.0\n"
                                             "  v_len: 0.0\n");
  SixDofSim sim(std::string(AUTOLAND_DATA_DIR) + "/AHAB_combined_betasym.stab",
                std::string(AUTOLAND_DATA_DIR) + "/aircraft.yaml", scenario);
  const SixDofTouchdown td = sim.run("test_sixdof_xwind.csv");

  REQUIRE(td.reached);
  CHECK(sim.stats().max_abs_y < 15.0);   // blown off but recovered/bounded
  CHECK(std::abs(td.y) < 5.0);           // near centerline at contact
  // Wind toward +y => crab nose-left (psi < 0) to hold the ground track.
  CHECK(td.psi < -0.5 * kDeg);
  CHECK(td.sink < 2.0);
}

// Waves on: touchdown happens at the instantaneous surface h = eta, not at
// h = 0, and the contact record carries the surface state.
TEST_CASE("Wave-field touchdown lands on eta, not on h=0", "[sixdof][sim]") {
  const std::string scenario = writeScenario("test_sixdof_waves.yaml",
                                             "waves:\n"
                                             "  enabled: true\n"
                                             "  regular: true\n"
                                             "  Hs: 0.4\n"
                                             "  Tp: 3.0\n"
                                             "  phase_deg: 45.0\n");
  SixDofSim sim(std::string(AUTOLAND_DATA_DIR) + "/AHAB_combined_betasym.stab",
                std::string(AUTOLAND_DATA_DIR) + "/aircraft.yaml", scenario);
  const SixDofTouchdown td = sim.run("test_sixdof_waves.csv");

  REQUIRE(td.reached);
  CHECK(td.h <= td.eta + 1e-12);          // contact test is h <= eta
  CHECK(td.h == Approx(td.eta).margin(0.05));  // within one step's descent
  CHECK(std::abs(td.eta) <= 0.2 + 1e-9);  // |eta| bounded by the amplitude
  CHECK(td.n_peak_flat >= 0.0);
}
