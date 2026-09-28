#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <fstream>
#include <random>
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
  x[XN] = 50.0 * i;
  return x;
}

// Minimal scenario writer for the closed-loop tests: calm base + optional
// wind/waves blocks appended verbatim.
std::string writeScenario(const std::string& name, const std::string& extra) {
  const std::string path = name;
  std::ofstream f(path);
  f << "plant: vspaero\n"
       "V_app: 18.0\n"
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
    for (int j : {PHI, THETA, PSI, H, Y, XN}) CHECK(wind[j] == still[j]);
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

// ---------------------------------------------------------------------------
// Corridor landing (world: block, pattern guidance, terrain, rollout).
// ---------------------------------------------------------------------------

// Straight-in on the corridor axis: the planner takes the DIRECT line (the
// aircraft is on the axis ahead of the FAF; a pure Dubins CSC plan would loop
// here, see test_pattern_guidance) and the run must reproduce the legacy calm
// touchdown (same trim sink / speed) through the guidance + rollout path,
// landing at the aim point.
TEST_CASE("Corridor straight-in lands at the aim point and rolls out inside",
          "[sixdof][corridor]") {
  SixDofSim sim(std::string(AUTOLAND_DATA_DIR) + "/AHAB_combined_betasym.stab",
                std::string(AUTOLAND_DATA_DIR) + "/aircraft.yaml",
                std::string(AUTOLAND_DATA_DIR) + "/beaver_corridor_straight.yaml");
  const SixDofTouchdown td = sim.run("test_sixdof_corridor_straight.csv");
  const SixDofRunStats& st = sim.stats();
  REQUIRE(st.trim_converged);
  REQUIRE(td.reached);
  CHECK_FALSE(td.terrain_strike);
  CHECK(td.in_corridor);
  const Corridor& cor = sim.scenario().world.corridor;
  CHECK(td.s_corr == Approx(cor.s_aim).margin(3.0));
  CHECK(std::abs(td.e_corr) < 3.0);
  CHECK(std::abs(td.dchi_corr) < 1.0 * kDeg);
  CHECK(td.sink == Approx(40.0 * std::sin(3.5 * kDeg)).margin(0.15));
  CHECK(td.V == Approx(40.0).margin(0.5));
  // Rollout: stopped, inside, never left the corridor.
  REQUIRE(td.rollout.ran);
  CHECK(td.rollout.stopped);
  CHECK(td.rollout.inside);
  CHECK(td.rollout.stayed_inside);
  CHECK(td.rollout.s_stop < cor.sFar());
  // The whole path is the direct line to the aim point.
  CHECK(st.plan_word == "DIRECT");
  CHECK(st.plan_length == Approx(3000.0).margin(30.0));
  CHECK(st.phi_sat_steps == 0);
  CHECK(st.max_abs_e_final < 3.0);
  // Shore trees under the final: cleared.
  CHECK(st.min_terrain_clearance > 5.0);
  CHECK(st.min_terrain_clearance < 60.0);  // and actually flown over them
}

// Landing-accuracy REQUIREMENTS shared by the pattern cases (not tuned to
// the observed values): touchdown inside the corridor with the ground
// course within 3 deg of the axis (a float touchdown crabbed more than that
// is a side-load case), lateral offset within a quarter of the half-width,
// wings within 3 deg of level at contact (no level-off exists yet, so the
// tracker must simply be converged), the bank command never beyond the
// guidance limit, airspeed inside the LR-556 validity band, and the stop
// point inside the corridor with a stopping-run margin.
static void checkPatternLanding(const SixDofSim& sim, const SixDofTouchdown& td) {
  const SixDofRunStats& st = sim.stats();
  const Corridor& cor = sim.scenario().world.corridor;
  const double phi_max = sim.scenario().guidance.phi_max;
  REQUIRE(st.trim_converged);
  REQUIRE(td.reached);
  CHECK_FALSE(td.terrain_strike);
  CHECK(td.in_corridor);
  CHECK(std::abs(td.e_corr) < 0.25 * 0.5 * cor.width);
  CHECK(std::abs(td.dchi_corr) < 3.0 * kDeg);
  CHECK(std::abs(td.phi) < 3.0 * kDeg);
  CHECK(td.sink < 3.0);
  CHECK(td.V == Approx(sim.scenario().V_app).margin(1.0));
  REQUIRE(td.rollout.ran);
  CHECK(td.rollout.stopped);
  CHECK(td.rollout.inside);
  CHECK(td.rollout.stayed_inside);
  CHECK(td.rollout.s_stop < cor.sFar() - 50.0);
  CHECK(st.max_abs_phi <= phi_max + 0.5 * kDeg);  // roll-loop overshoot only
  CHECK(st.phi_sat_steps == 0);
  CHECK(st.min_V > 35.0);       // LR-556 validity band
  CHECK(st.max_alpha < 12.0 * kDeg);
  // Shore trees under the final: actually overflown (sentinel would be 1e30)
  // and cleared by the CG by more than 5 m.
  CHECK(st.min_terrain_clearance < 60.0);
  CHECK(st.min_terrain_clearance > 5.0);
}

// Overhead entry, opposite heading: the pattern (Dubins CSC + final) must
// bring the aircraft round, down the glideslope, into the corridor with the
// course aligned, over the shore trees, and the rollout must stop inside.
TEST_CASE("Corridor overhead entry flies the pattern into the corridor",
          "[sixdof][corridor][slow]") {
  SixDofSim sim(std::string(AUTOLAND_DATA_DIR) + "/AHAB_combined_betasym.stab",
                std::string(AUTOLAND_DATA_DIR) + "/aircraft.yaml",
                std::string(AUTOLAND_DATA_DIR) + "/beaver_corridor_overhead.yaml");
  const SixDofTouchdown td = sim.run("test_sixdof_corridor_overhead.csv");
  checkPatternLanding(sim, td);
  const SixDofRunStats& st = sim.stats();
  // A real pattern was flown: a turn of > 90 deg at a substantial bank.
  CHECK(st.max_abs_phi > 15.0 * kDeg);
  CHECK(st.plan_length > 3500.0);
  CHECK(std::abs(td.crab) < 2.0 * kDeg);  // no wind: heading ~ course
  CHECK(st.max_xte < 20.0);               // tracking, not just arrival
}

// Same entry with a 5 m/s steady crosswind (plant-side step at t = 0, so the
// first seconds are a gust response, not an equilibrium): the ground course
// is tracked, the touchdown is crabbed into wind, everything else holds.
TEST_CASE("Corridor overhead entry in a crosswind tracks the ground course",
          "[sixdof][corridor][slow]") {
  SixDofSim sim(std::string(AUTOLAND_DATA_DIR) + "/AHAB_combined_betasym.stab",
                std::string(AUTOLAND_DATA_DIR) + "/aircraft.yaml",
                std::string(AUTOLAND_DATA_DIR) + "/beaver_corridor_crosswind.yaml");
  const SixDofTouchdown td = sim.run("test_sixdof_corridor_crosswind.csv");
  checkPatternLanding(sim, td);
  // Crab into a 5 m/s east wind at ~40 m/s: asin(5/40) = 7.2 deg (the
  // heading error is crab-dominated; the course error is what the corridor
  // requirement bounds).
  CHECK(td.crab < -4.0 * kDeg);
  CHECK(td.crab > -9.0 * kDeg);
  CHECK(std::abs(td.dpsi_corr - td.crab) < 3.0 * kDeg);
  CHECK(sim.stats().max_xte < 25.0);
}

// The guidance outer-loop refactor: stepOuter() must be bit-identical to
// stepOuterRef() at the config references, for both longitudinal modes,
// over a long random-input sequence (integrator states included).
TEST_CASE("stepOuter is bit-identical to stepOuterRef at the config references",
          "[sixdof][corridor]") {
  for (LonMode mode : {LonMode::Cascade, LonMode::Tecs}) {
    SixDofNominalConfig cfg;
    cfg.V_ref = 40.0; cfg.gamma_ref = -3.5 * kDeg; cfg.lon_mode = mode;
    cfg.theta_trim = 6.0 * kDeg; cfg.dT_trim = 0.4;
    cfg.tecs.tas_min = 30.0; cfg.tecs.tas_max = 55.0;
    SixDofNominal a(cfg), b(cfg);
    std::mt19937 rng(5);
    std::uniform_real_distribution<double> un(-1.0, 1.0);
    for (int k = 0; k < 500; ++k) {
      StateVec x = sampleState(k % 4);
      x[U] += 20.0 + 5.0 * un(rng); x[THETA] += 0.1 * un(rng); x[H] = 100.0 + 50.0 * un(rng);
      const double V = 38.0 + 4.0 * un(rng), Vd = 0.5 * un(rng);
      const SixDofNominal::LonOuterCmd oa = a.stepOuter(x, V, Vd, 0.01);
      const SixDofNominal::LonOuterCmd ob = b.stepOuterRef(x, V, Vd, cfg.gamma_ref, cfg.V_ref, 0.01);
      CHECK(oa.theta_cmd == ob.theta_cmd);
      CHECK(oa.dT == ob.dT);
    }
  }
}

// ---------------------------------------------------------------------------
// NOMINAL-FAILURE BASELINES (data/beaver_corridor_fail_*.yaml): cases built so
// the pattern nominal alone does NOT land safely. Each pins the failure
// signature the corresponding barrier row (roadmap Phases 1-5) must flip;
// when a phase lands, its case moves out of this test into a filter:on/off
// A/B pair. If a nominal change makes one of these pass, revisit the case --
// the point is to keep a ground truth where a barrier is actually needed.
// ---------------------------------------------------------------------------
namespace {
struct FailRun {
  SixDofTouchdown td;
  SixDofRunStats st;
  double V_app{0};
  Corridor cor;
};
FailRun runFail(const std::string& name) {
  SixDofSim sim(std::string(AUTOLAND_DATA_DIR) + "/AHAB_combined_betasym.stab",
                std::string(AUTOLAND_DATA_DIR) + "/aircraft.yaml",
                std::string(AUTOLAND_DATA_DIR) + "/beaver_corridor_fail_" + name + ".yaml");
  FailRun r;
  r.td = sim.run("test_sixdof_fail_" + name + ".csv");
  r.st = sim.stats();
  r.V_app = sim.scenario().V_app;
  r.cor = sim.scenario().world.corridor;
  REQUIRE(r.st.trim_converged);
  REQUIRE(r.td.reached);
  return r;
}
}  // namespace

TEST_CASE("Nominal-failure baselines: where a barrier row is actually needed",
          "[sixdof][corridor][baseline][slow]") {
  SECTION("ridge: tree line on the final -> terrain strike (Phase 2 terrain row)") {
    const FailRun r = runFail("ridge");
    CHECK(r.td.terrain_strike);
    CHECK(r.td.hT > 20.0);                 // hit the tree line, not the water
    CHECK(r.st.min_terrain_clearance < 0.0);
  }
  SECTION("mountain: under the pattern loop -> terrain strike (Phase 2 terrain row + planner)") {
    const FailRun r = runFail("mountain");
    CHECK(r.td.terrain_strike);
    CHECK(r.td.hT > 100.0);
    CHECK(r.td.h > 100.0);                 // mid-pattern, not on final
  }
  SECTION("short lake, hot: on-speed touchdown overruns the far edge (Phase 3 energy ceiling)") {
    const FailRun r = runFail("short_lake");
    CHECK(r.td.in_corridor);               // the touchdown itself is fine...
    CHECK(r.td.V > r.V_app - 1.0);         // ...at the hot approach speed
    REQUIRE(r.td.rollout.ran);
    CHECK_FALSE(r.td.rollout.inside);      // ...and the water run does not fit
    CHECK(r.td.rollout.s_stop > r.cor.sFar());
  }
  SECTION("engine out: pattern infeasible on the idle glide -> lands short (Phase 3 energy floor)") {
    const FailRun r = runFail("engine_out");
    CHECK_FALSE(r.td.in_corridor);
    CHECK(r.td.s_corr < r.cor.sNear() - 500.0);   // far short
    CHECK(r.st.max_alpha > 16.0 * kDeg);          // and through the stall angle
    CHECK(r.st.min_V < 30.0);                     // below the LR-556 band
  }
  SECTION("gust on final: ground track held, floats touch with large crab (Phase 4 decrab row)") {
    const FailRun r = runFail("gust_final");
    CHECK(std::abs(r.td.e_corr) < 2.0);           // cross-track is NOT the failure
    CHECK(std::abs(r.td.dchi_corr) < 1.0 * kDeg);
    CHECK(std::abs(r.td.crab) > 10.0 * kDeg);     // the crab is
  }
  SECTION("slow approach: decelerates into the stall region (Phase 1 AoA row)") {
    const FailRun r = runFail("slow");
    CHECK(r.st.stall_entered);
    CHECK(r.st.max_alpha > 16.0 * kDeg);
    CHECK(r.st.min_V < 33.0);
    CHECK_FALSE(r.td.success);
  }
  SECTION("gust in the turn: AoA transient into the stall region (Phase 1 AoA row -- DONE, see cbf_gust_turn)") {
    const FailRun r = runFail("gust_turn");
    CHECK(r.st.stall_entered);
    CHECK(r.st.max_alpha > 16.0 * kDeg);
    CHECK_FALSE(r.td.success);
  }
  SECTION("gap between two mountains on the final -> terrain strike (Phase 2; the terrain row does NOT yet flip it)") {
    const FailRun r = runFail("gap_final");
    CHECK(r.td.terrain_strike);
    CHECK(r.td.hT > 50.0);
    CHECK(std::abs(r.td.e_corr) < 5.0);   // flew the axis straight into the flank
  }
  SECTION("tailwind on a short lake: airspeed on target, ground run overruns (Phase 3 energy ceiling)") {
    const FailRun r = runFail("tailwind");
    CHECK(r.td.in_corridor);
    CHECK(r.td.V == Approx(r.V_app).margin(1.5));
    REQUIRE(r.td.rollout.ran);
    CHECK_FALSE(r.td.rollout.inside);
    CHECK(r.td.rollout.s_stop > r.cor.sFar());
    CHECK_FALSE(r.td.success);
  }
  SECTION("tight pattern: bank command saturated, flown bank over the limit (Phase 1 bank row)") {
    const FailRun r = runFail("tight_pattern");
    CHECK(r.st.phi_sat_steps > 0.15 * r.st.steps);
    CHECK(r.st.max_abs_phi > 25.0 * kDeg);
    CHECK(r.st.max_xte > 30.0);
  }
}

// ---------------------------------------------------------------------------
// Phase 1: the 6-DOF surfaces-only CBF filter on the corridor cases.
// ---------------------------------------------------------------------------
namespace {
struct CbfRun {
  SixDofTouchdown td;
  SixDofRunStats st;
  SixDofScenario sc;
};
CbfRun runCbf(const std::string& name) {
  SixDofSim sim(std::string(AUTOLAND_DATA_DIR) + "/AHAB_combined_betasym.stab",
                std::string(AUTOLAND_DATA_DIR) + "/aircraft.yaml",
                std::string(AUTOLAND_DATA_DIR) + "/beaver_corridor_cbf_" + name + ".yaml");
  CbfRun r;
  r.td = sim.run("test_sixdof_cbf_" + name + ".csv");
  r.st = sim.stats();
  r.sc = sim.scenario();
  REQUIRE(r.st.trim_converged);
  REQUIRE(r.td.reached);
  return r;
}
}  // namespace

// Rare intervention by construction: with the envelope limits a margin above
// the nominal's clamps, the filter must stay quiet on every well-flown pattern
// (at most a handful of steps, never best-effort) and the landing verdict is
// unchanged.
TEST_CASE("6-DOF CBF stays quiet on the landing cases", "[sixdof][corridor][cbf][slow]") {
  for (const char* name : {"overhead", "crosswind", "straight"}) {
    const CbfRun r = runCbf(name);
    INFO(name);
    CHECK(r.td.success);
    CHECK(r.st.cbf6_active_steps <= 5);
    CHECK(r.st.cbf6_best_effort_steps == 0);
    CHECK_FALSE(r.st.stall_entered);
    CHECK(r.st.cbf6_time_mean_us < 5000.0);  // compute budget: well under a 10 ms step
  }
}

// Tight pattern A/B: the tracker may ask for 35 deg; with the hard bank row
// the FLOWN bank never exceeds 25 deg (sampled-data slack of 0.5 deg), the
// row never goes best-effort, and the landing still succeeds. Without the
// filter the same nominal banks past 25.
TEST_CASE("6-DOF CBF hard bank row holds the flown bank (tight pattern A/B)",
          "[sixdof][corridor][cbf][slow]") {
  const CbfRun on = runCbf("tight_pattern");
  const CbfRun off = runCbf("tight_pattern_off");
  const double lim = on.sc.cbf.phi_max;
  CHECK(off.st.max_abs_phi > lim + 1.0 * kDeg);       // the nominal alone exceeds it (26.4 vs 25)
  CHECK(on.st.max_abs_phi <= lim + 0.5 * kDeg);       // the row holds it
  CHECK(on.st.cbf6_min_h_phi > -0.5 * kDeg);
  CHECK(on.st.cbf6_best_effort_steps == 0);
  CHECK(on.st.cbf6_active_steps > 1000);              // it actually worked
  CHECK(on.td.success);
  CHECK_FALSE(on.st.stall_entered);
}

// Slow approach A/B: the nominal decelerates into the stall; with the hard
// AoA ceiling the flown alpha never exceeds the limit and the stall region is
// never entered. The landing is NOT claimed: an approach at 31 m/s cannot be
// flown inside the band, so the aircraft runs out of altitude elsewhere --
// keeping it out of the stall is the whole of Phase 1's claim here.
TEST_CASE("6-DOF CBF hard AoA row keeps the slow approach out of the stall (A/B)",
          "[sixdof][corridor][cbf][slow]") {
  const CbfRun on = runCbf("slow");
  const CbfRun off = runCbf("slow_off");
  CHECK(off.st.stall_entered);
  CHECK(off.st.max_alpha > 16.0 * kDeg);
  CHECK_FALSE(on.st.stall_entered);
  CHECK(on.st.max_alpha <= on.sc.cbf.alpha_max + 0.3 * kDeg);
  CHECK(on.st.cbf6_min_h_alpha > -0.3 * kDeg);
  CHECK(on.st.cbf6_best_effort_steps < 10);
  CHECK(on.st.min_V > 33.0);
}

// Gust in the turn A/B: the nominal alone lets the sideslip-induced pitch-up
// carry alpha through the stall line; with the hard AoA row the transient is
// clipped and the landing proceeds unchanged. The first case where a barrier
// row turns a FAILED run into a landing.
TEST_CASE("6-DOF CBF hard AoA row clips the gust-in-turn transient and the landing proceeds (A/B)",
          "[sixdof][corridor][cbf][slow]") {
  const CbfRun on = runCbf("gust_turn");
  const CbfRun off = runCbf("gust_turn_off");
  CHECK(off.st.stall_entered);
  CHECK_FALSE(off.td.success);
  CHECK_FALSE(on.st.stall_entered);
  CHECK(on.st.max_alpha <= on.sc.cbf.alpha_max + 0.3 * kDeg);
  CHECK(on.st.cbf6_best_effort_steps == 0);
  CHECK(on.st.cbf6_active_steps > 100);
  CHECK(on.td.success);
  CHECK(std::abs(on.td.e_corr) < 0.25 * 0.5 * on.sc.world.corridor.width);
}

// NEGATIVE RESULT, pinned: the terrain row (braking form, climb-capability
// lookahead fan, hard, with the envelope rows hard) does not steer the
// straight-in through the gap between two mountains -- the run ends in a
// terrain strike and a stall. See corridor_landing_roadmap.md "Phase 2 as
// built" for the analysis: within a fan wide enough to keep the barrier
// continuous the safe set along the final is EMPTY (some direction in the
// fan always sees a mountain), and a single ray narrow enough to fit the gap
// makes the barrier discontinuous under yaw. The gap is a corridor problem
// (Phase 3 cross-track funnel with the gap as an intermediate gate), not a
// terrain-clearance problem. Flip this test when a formulation lands.
TEST_CASE("6-DOF CBF terrain row does not (yet) steer through the mountain gap",
          "[sixdof][corridor][cbf][baseline][slow]") {
  const CbfRun r = runCbf("gap_final");
  CHECK_FALSE(r.td.success);
  CHECK(r.st.cbf6_best_effort_steps > 100);   // the hard set is infeasible along the final
}
