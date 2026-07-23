#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <string>

#include "autoland/lon_sim.hpp"

using namespace autoland;

namespace {
const std::string kData = AUTOLAND_DATA_DIR;
constexpr double kDeg = M_PI / 180.0;
}  // namespace

// End-to-end closed-loop landing: the bundled nominal scenario -- a powered
// -3 deg approach from 40 m -- flown to touchdown on the real AHAB aero deck
// under the full recovery-set CBF-QP (impact HARD, stall/nose-up/energy soft,
// thrust guards). Deterministic: the altitude-sensor noise uses a fixed seed.
//
// This is the flagship behavior every unit above serves, pinned as a test:
// the aircraft must actually reach the water (no float-forever / divergence),
// inside the touchdown envelope, with the filter healthy the whole way down
// (no best-effort recoveries, no dropped HARD rows). The bounds are envelope
// limits and deliberately looser than the current baseline (t=68.5 s,
// sink=0.21 m/s, V=13.37 m/s, theta=3.0 deg) so retuning has headroom, but a
// real regression -- a hot, fast, nose-down, or non-arriving touchdown --
// fails loudly.
TEST_CASE("nominal -3 deg approach lands within the touchdown envelope",
          "[lon_sim]") {
  LonSim sim(kData + "/AHAB_combined.stab", kData + "/aircraft.yaml",
             kData + "/lon_scenario.yaml");
  REQUIRE(sim.scenario().nominal.gamma_ref == Catch::Approx(-3.0 * kDeg));

  const LonTouchdown td = sim.run("test_nominal_landing_trace.csv");

  REQUIRE(td.reached);              // descends all the way to touchdown
  CHECK(td.t > 20.0);               // ... on a real approach timescale
  CHECK(td.t < 150.0);
  CHECK(td.V <= 14.0 + 1e-6);       // V_td_max energy ceiling holds at contact
  CHECK(td.sink > 0.0);             // still descending at contact
  CHECK(td.sink < 0.5);             // soft touchdown (envelope, not the tune)
  CHECK(td.theta > 0.0);            // nose-up attitude at contact
  CHECK(td.gamma < 0.0);            // on a descending path ...
  CHECK(td.gamma > -3.0 * kDeg);    // ... flared shallower than the approach
  CHECK(td.recoveries == 0);        // hard set stayed feasible every step
  CHECK(td.steps_hard_dropped == 0);  // no HARD row lost to bad assembly
}
