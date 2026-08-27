#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <fstream>
#include <string>

#include "autoland/px4_tecs.hpp"
#include "autoland/sixdof_sim.hpp"

using namespace autoland;
using Catch::Approx;

namespace {
constexpr double kDeg = M_PI / 180.0;

// PX4's own fixture (src/lib/tecs/TECSTest.cpp makeParam/makeInput/
// makeSetpoint, main @ a906b728), values verbatim, so the unit checks run the
// port on the same numbers upstream tests against.
px4::TecsParam makeParam() {
  px4::TecsParam p;
  p.max_sink_rate = 5.0;
  p.min_sink_rate = 2.0;
  p.max_climb_rate = 5.0;
  p.vert_accel_limit = 10.0;
  p.equivalent_airspeed_trim = 15.0;
  p.tas_min = 10.0;
  p.tas_max = 30.0;
  p.pitch_max = 15.0 * kDeg;
  p.pitch_min = -15.0 * kDeg;
  p.throttle_trim = 0.5;
  p.throttle_max = 1.0;
  p.throttle_min = 0.0;
  p.altitude_error_gain = 0.2;
  p.altitude_setpoint_gain_ff = 0.0;
  p.tas_error_percentage = 0.15;
  p.airspeed_error_gain = 0.1;
  p.ste_rate_time_const = 0.1;
  p.seb_rate_ff = 1.0;
  p.pitch_speed_weight = 1.0;
  p.integrator_gain_pitch = 0.4;
  p.pitch_damping_gain = 0.1;
  p.integrator_gain_throttle = 0.3;
  p.throttle_damping_gain = 0.1;
  p.throttle_slewrate = 0.0;
  p.load_factor_correction = 0.0;
  p.load_factor = 1.0;
  return p;
}
px4::TecsFlag makeFlag() { return {true, true}; }
px4::TecsInput makeInput() { return {100.0, 0.0, 15.0, 0.0}; }
px4::TecsSetpoint makeSetpoint() {
  px4::TecsSetpoint s;
  s.altitude_reference.alt = 100.0;
  s.altitude_reference.alt_rate = 0.0;
  s.altitude_rate_setpoint_direct = NAN;
  s.tas_setpoint = 15.0;
  return s;
}
constexpr double kDt = 0.02;  // upstream test step

// Beaver scenario writer for the closed-loop TECS cases (calm base +
// optional blocks appended verbatim).
std::string writeBeaverTecs(const std::string& name, const std::string& initial_extra,
                            const std::string& extra = "") {
  std::ofstream f(name);
  f << "plant: beaver\n"
       "V_app: 40.0\n"
       "gamma_app_deg: -3.5\n"
       "beaver:\n"
       "  n_rpm: 1800.0\n"
       "dt: 0.01\n"
       "t_max: 120.0\n"
       "initial:\n"
       "  h0: 60.0\n"
    << initial_extra
    << "nominal:\n"
       "  type: tecs\n"
    << extra;
  return name;
}
}  // namespace

// At the setpoint (on altitude, on speed, zero rates) the law must return the
// trim throttle and zero pitch-above-trim, and hold them: no integrator drift.
TEST_CASE("TECS at the setpoint holds trim throttle and zero pitch", "[tecs]") {
  px4::TecsControl c;
  const px4::TecsParam p = makeParam();
  const px4::TecsFlag fl = makeFlag();
  const px4::TecsInput in = makeInput();
  const px4::TecsSetpoint sp = makeSetpoint();

  c.initialize(sp, in, p, fl);
  CHECK(c.getThrottleSetpoint() == Approx(0.5).margin(1e-12));
  CHECK(c.getPitchSetpoint() == Approx(0.0).margin(1e-12));
  CHECK(c.getRatioUndersped() == 0.0);

  for (int i = 0; i < 500; ++i) c.update(kDt, sp, in, p, fl);
  CHECK(c.getThrottleSetpoint() == Approx(0.5).margin(1e-12));
  CHECK(c.getPitchSetpoint() == Approx(0.0).margin(1e-12));
  CHECK(c.getDebugOutput().pitch_integrator == Approx(0.0).margin(1e-12));
  CHECK(c.getDebugOutput().throttle_integrator == Approx(0.0).margin(1e-12));
}

// Energy bookkeeping signs: total-energy errors go to the throttle, energy-
// BALANCE errors go to the pitch (trade speed for height and vice versa).
TEST_CASE("TECS energy signs: throttle on total energy, pitch on balance", "[tecs]") {
  const px4::TecsParam p = makeParam();
  const px4::TecsFlag fl = makeFlag();
  const px4::TecsSetpoint sp = makeSetpoint();
  auto run = [&](double alt, double tas) {
    px4::TecsControl c;
    px4::TecsInput in = makeInput();
    in.altitude = alt;
    in.tas = tas;
    c.initialize(sp, in, p, fl);
    c.update(kDt, sp, in, p, fl);
    return c;
  };

  SECTION("low and slow -> more throttle, pitch from the balance only") {
    const auto c = run(90.0, 13.0);
    CHECK(c.getThrottleSetpoint() > 0.5);
    // SPE demand +g*2 (alt err 10 m * 0.2), SKE demand 15*(-0.2)... net
    // balance = spe_sp - ske_sp = 19.6 + 3 > 0 -> nose up (height priority
    // wins at equal weighting because the altitude term is larger here).
    CHECK(c.getDebugOutput().energy_balance_rate_sp > 0.0);
  }
  SECTION("high and fast -> less throttle") {
    CHECK(run(110.0, 17.0).getThrottleSetpoint() < 0.5);
  }
  SECTION("fast on altitude -> nose up (speed into height), throttle back") {
    const auto c = run(100.0, 17.0);
    CHECK(c.getPitchSetpoint() > 0.0);
    CHECK(c.getThrottleSetpoint() < 0.5);
  }
  SECTION("slow on altitude -> nose down (height into speed), throttle up") {
    const auto c = run(100.0, 13.0);
    CHECK(c.getPitchSetpoint() < 0.0);
    CHECK(c.getThrottleSetpoint() > 0.5);
  }
  SECTION("low on speed -> nose up and throttle up") {
    const auto c = run(90.0, 15.0);
    CHECK(c.getPitchSetpoint() > 0.0);
    CHECK(c.getThrottleSetpoint() > 0.5);
  }
}

// Altitude-loop numbers against the upstream formulas: hdot_sp = err/TC (+ff),
// clamped to the climb/sink limits; SPE/SKE rates as g*hdot and V*Vdot.
TEST_CASE("TECS altitude loop and specific-energy rates match the PX4 formulas",
          "[tecs]") {
  px4::TecsControl c;
  const px4::TecsParam p = makeParam();
  const px4::TecsFlag fl = makeFlag();
  px4::TecsSetpoint sp = makeSetpoint();
  px4::TecsInput in = makeInput();

  in.altitude = 95.0;   // 5 m low -> 5 * 0.2 = 1 m/s climb demand
  in.altitude_rate = 0.3;
  in.tas = 14.0;        // 1 m/s slow -> 0.1 m/s^2 demand (within 0.5*STE limits / V)
  in.tas_rate = 0.05;
  c.initialize(sp, in, p, fl);
  c.update(kDt, sp, in, p, fl);
  const px4::TecsDebugOutput& d = c.getDebugOutput();
  CHECK(d.altitude_rate_control == Approx(1.0));
  CHECK(d.true_airspeed_derivative_control == Approx(0.1));
  // STE_sp = g*1 + 15*0.1 ; SEB_sp = g*1 - 15*0.1 (weight 1 / 1)
  CHECK(d.total_energy_rate_sp == Approx(px4::kOneG * 1.0 + 15.0 * 0.1));
  CHECK(d.energy_balance_rate_sp == Approx(px4::kOneG * 1.0 - 15.0 * 0.1));
  CHECK(d.energy_balance_rate_estimate == Approx(px4::kOneG * 0.3 - 14.0 * 0.05));

  // Clamp: 50 m low asks 10 m/s but max_climb_rate is 5.
  in.altitude = 50.0;
  c.update(kDt, sp, in, p, fl);
  CHECK(c.getDebugOutput().altitude_rate_control == Approx(5.0));
}

// The pitch setpoint is rate-limited by vert_accel_limit / V per step and the
// integrator is frozen in the saturating direction at the pitch limit.
TEST_CASE("TECS pitch: vertical-acceleration rate limit and anti-windup", "[tecs]") {
  px4::TecsControl c;
  const px4::TecsParam p = makeParam();
  const px4::TecsFlag fl = makeFlag();
  const px4::TecsSetpoint sp = makeSetpoint();
  px4::TecsInput in = makeInput();

  c.initialize(sp, in, p, fl);  // at the setpoint: pitch 0
  in.altitude = 50.0;           // large step -> pitch demand saturates
  const double max_step = kDt * p.vert_accel_limit / in.tas;
  double prev = c.getPitchSetpoint();
  for (int i = 0; i < 5; ++i) {
    c.update(kDt, sp, in, p, fl);
    CHECK(std::abs(c.getPitchSetpoint() - prev) <= max_step + 1e-12);
    prev = c.getPitchSetpoint();
  }
  for (int i = 0; i < 200; ++i) c.update(kDt, sp, in, p, fl);
  CHECK(c.getPitchSetpoint() == Approx(p.pitch_max));
  // Once saturated at pitch_max the integrator input is clipped to <= 0, so
  // the state cannot keep growing; record and confirm.
  const double integ_at_sat = c.getDebugOutput().pitch_integrator;
  for (int i = 0; i < 200; ++i) c.update(kDt, sp, in, p, fl);
  CHECK(c.getDebugOutput().pitch_integrator <= integ_at_sat + 1e-12);
  CHECK(c.getPitchSetpoint() == Approx(p.pitch_max));
}

// Underspeed: below tas_min - 2 * 0.15 * V_trim the ratio is 1 -> full
// throttle, speed-only weighting (nose DOWN to recover speed), and the
// throttle integrator is frozen.
TEST_CASE("TECS underspeed ramps in full throttle and speed priority", "[tecs]") {
  px4::TecsControl c;
  const px4::TecsParam p = makeParam();
  const px4::TecsFlag fl = makeFlag();
  const px4::TecsSetpoint sp = makeSetpoint();
  px4::TecsInput in = makeInput();
  in.tas = 5.0;  // fully undersped: 10 - 0.15*15 - 0.15*15 = 5.5
  c.initialize(sp, in, p, fl);
  for (int i = 0; i < 20; ++i) c.update(kDt, sp, in, p, fl);
  CHECK(c.getRatioUndersped() == Approx(1.0));
  CHECK(c.getThrottleSetpoint() == Approx(p.throttle_max));
  CHECK(c.getPitchSetpoint() < 0.0);
  CHECK(c.getDebugOutput().throttle_integrator == Approx(0.0).margin(1e-12));

  // Half-way through the soft band (tas = 6.625 -> ratio 0.5).
  in.tas = 0.5 * (5.5 + 7.75);
  c.update(kDt, sp, in, p, fl);
  CHECK(c.getRatioUndersped() == Approx(0.5));

  // Detection disabled: ratio 0 regardless.
  px4::TecsFlag no_det = fl;
  no_det.detect_underspeed_enabled = false;
  in.tas = 5.0;
  c.update(kDt, sp, in, p, no_det);
  CHECK(c.getRatioUndersped() == 0.0);
}

// Throttle slew: |dT| per step <= dt * (thr_max - thr_min) * slewrate.
TEST_CASE("TECS throttle slew limit bounds the per-step change", "[tecs]") {
  px4::TecsControl c;
  px4::TecsParam p = makeParam();
  p.throttle_slewrate = 0.5;
  const px4::TecsFlag fl = makeFlag();
  const px4::TecsSetpoint sp = makeSetpoint();
  px4::TecsInput in = makeInput();
  c.initialize(sp, in, p, fl);
  in.altitude = 50.0;  // demands much more throttle
  const double inc = kDt * (p.throttle_max - p.throttle_min) * p.throttle_slewrate;
  double prev = c.getThrottleSetpoint();
  for (int i = 0; i < 20; ++i) {
    c.update(kDt, sp, in, p, fl);
    CHECK(std::abs(c.getThrottleSetpoint() - prev) <= inc + 1e-12);
    prev = c.getThrottleSetpoint();
  }
  CHECK(c.getThrottleSetpoint() > 0.5);
}

// A finite direct height-rate setpoint bypasses the altitude loop entirely
// (the sim's glideslope path): a huge altitude error must not act.
TEST_CASE("TECS direct height-rate setpoint bypasses the altitude loop", "[tecs]") {
  px4::TecsControl c;
  const px4::TecsParam p = makeParam();
  const px4::TecsFlag fl = makeFlag();
  px4::TecsSetpoint sp = makeSetpoint();
  px4::TecsInput in = makeInput();
  c.initialize(sp, in, p, fl);
  in.altitude = 50.0;
  sp.altitude_rate_setpoint_direct = 0.0;
  for (int i = 0; i < 50; ++i) c.update(kDt, sp, in, p, fl);
  CHECK(c.getDebugOutput().altitude_rate_control == 0.0);
  CHECK(c.getThrottleSetpoint() == Approx(0.5).margin(1e-12));
  CHECK(c.getPitchSetpoint() == Approx(0.0).margin(1e-12));

  // Descent demand of -2 m/s at constant speed: the predicted throttle sits
  // between trim and idle (STE_sp = -2g = STE_rate_min -> exactly thr_min +
  // damping/integrator on the estimate error).
  sp.altitude_rate_setpoint_direct = -2.0;
  in.altitude_rate = -2.0;  // already sinking at the demand: zero error
  // The pitch demand walks down under the vert_accel_limit rate limit
  // (10/15 * 0.02 = 0.0133 rad per step) -- 10 steps to reach the FF value.
  for (int i = 0; i < 20; ++i) c.update(kDt, sp, in, p, fl);
  CHECK(c.getThrottleSetpoint() == Approx(p.throttle_min).margin(1e-9));
  // Pitch FF = -2g / (V g) = -2/15 rad above trim (zero error: no integrator).
  CHECK(c.getPitchSetpoint() == Approx(-2.0 / 15.0).margin(1e-9));
}

// Bank-angle induced-drag compensation raises the total-energy demand by
// load_factor_correction * (n - 1).
TEST_CASE("TECS load-factor correction adds to the total-energy demand", "[tecs]") {
  px4::TecsParam p = makeParam();
  p.load_factor_correction = 15.0;
  p.load_factor = 1.0 / std::cos(30.0 * kDeg);
  const px4::TecsFlag fl = makeFlag();
  const px4::TecsSetpoint sp = makeSetpoint();
  const px4::TecsInput in = makeInput();
  px4::TecsControl c;
  c.initialize(sp, in, p, fl);
  c.update(kDt, sp, in, p, fl);
  CHECK(c.getDebugOutput().total_energy_rate_sp ==
        Approx(15.0 * (p.load_factor - 1.0)));
  CHECK(c.getThrottleSetpoint() > 0.5);
}

// Closed loop on the Beaver plant: the TECS nominal must fly the calm
// straight-in to touchdown on the glideslope at V_app, to the same tolerances
// the cascaded PID meets (test_beaver_dynamics.cpp), with the plant-derived
// anchors physically ordered.
TEST_CASE("Beaver TECS calm straight-in lands on the glideslope", "[tecs][beaver6][sim]") {
  SixDofSim sim("", "", std::string(AUTOLAND_DATA_DIR) + "/beaver_landing_tecs.yaml");
  REQUIRE(sim.trimResult().converged);
  const SixDofNominalConfig& nom = sim.scenario().nominal;
  REQUIRE(nom.lon_mode == LonMode::Tecs);
  const px4::TecsParam& tp = nom.tecs;
  // Anchors: idle sink < 0 < level < full-throttle climb, in throttle order;
  // the level-flight pitch sits above the descent trim pitch at the same V.
  CHECK(tp.max_climb_rate > 0.5);
  CHECK(tp.min_sink_rate > 0.5);
  CHECK(tp.max_sink_rate >= tp.min_sink_rate);
  CHECK(tp.throttle_trim > nom.limits.dT_min);
  CHECK(tp.throttle_trim < nom.limits.dT_max);
  CHECK(tp.throttle_trim > nom.dT_trim);  // level needs more power than the descent
  CHECK(nom.tecs_pitch_offset > nom.theta_trim);

  const SixDofTouchdown td = sim.run("test_beaver_tecs_calm.csv");
  REQUIRE(td.reached);
  CHECK(td.sink == Approx(40.0 * std::sin(3.5 * kDeg)).margin(0.15));
  CHECK(td.V == Approx(40.0).margin(0.3));
  CHECK(std::abs(td.gamma + 3.5 * kDeg) < 0.3 * kDeg);
  CHECK(std::abs(td.phi) < 1.0 * kDeg);
  CHECK(std::abs(td.y) < 1.0);
  CHECK(sim.stats().max_alpha < 12.0 * kDeg);
  CHECK(sim.stats().min_V > 35.0);
}

// Hot-and-high entry (+3 m/s, +2 deg pitch): TECS must trade the excess
// energy away and settle on V_app and the glideslope before touchdown.
TEST_CASE("Beaver TECS hot entry settles on V_app before touchdown", "[tecs][beaver6][sim]") {
  const std::string scenario = writeBeaverTecs("test_beaver_tecs_hot.yaml",
                                               "  dV: 3.0\n"
                                               "  dtheta_deg: 2.0\n");
  SixDofSim sim("", "", scenario);
  REQUIRE(sim.trimResult().converged);
  const SixDofTouchdown td = sim.run("test_beaver_tecs_hot.csv");
  REQUIRE(td.reached);
  CHECK(td.V == Approx(40.0).margin(0.5));
  CHECK(std::abs(td.gamma + 3.5 * kDeg) < 0.5 * kDeg);
  CHECK(td.sink < 2.8);
  CHECK(sim.stats().max_alpha < 12.0 * kDeg);
}
