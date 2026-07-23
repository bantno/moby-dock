#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <filesystem>
#include <string>

#include "autoland/lon_sim.hpp"

// =============================================================================
// Emergent stall-recovery integration tests (end-to-end LonSim runs).
//
// The nominal controller commands only a constant gamma_ref + constant T_set
// -- it has no flare and no stall-recovery logic. These tests pin the claim
// that on approaches which stall near the water, the CBF-QP filter recovers
// the aircraft and lands it safely anyway. The scenarios are LOCKED yamls
// (data/lon_stall_{pull,entry}_cbf.yaml) shared with the Python suite
// (scripts/stall_recovery_suite.py); edit numbers only together.
//
// Threshold philosophy: the hard-row checks (impact psi minima, no dropped
// HARD rows, thrust bounds) are tuning-independent safety invariants; the
// sink/V/theta bounds are per-scenario envelopes frozen with >= 1.5x margin
// over the calibrated values so barrier retuning doesn't churn them.
// =============================================================================

using namespace autoland;

namespace {
const std::string kData = AUTOLAND_DATA_DIR;
constexpr double kDeg = M_PI / 180.0;

struct RunResult {
  LonTouchdown td;
  LonRunStats stats;
};

// Run a locked scenario end-to-end; traces go under the ctest working dir.
RunResult runScenario(const std::string& scenario_yaml, const std::string& name) {
  std::filesystem::create_directories("stall_recovery_test_out");
  LonSim sim(kData + "/AHAB_combined.stab", kData + "/aircraft.yaml",
             kData + "/" + scenario_yaml);
  LonTouchdown td = sim.run("stall_recovery_test_out/" + name + ".csv");
  return {td, sim.stats()};
}

// The tuning-independent safety invariants every CBF-on run must satisfy.
void checkHardInvariants(const RunResult& r) {
  // A dropped HARD row is a silent loss of the safety guarantee.
  CHECK(r.stats.steps_hard_dropped == 0);
  // Impact barrier (the only HARD safety-envelope row) held over its active
  // window. psi1 is the state-only forward-invariance condition for the
  // degree-2 row; the residual confirms the QP row held with the APPLIED
  // control (psi2 drift-only goes negative exactly when authority is spent,
  // so it is deliberately NOT asserted). 1e-3 covers the OSQP tolerance.
  CHECK(r.stats.min_psi1_impact >= -1e-6);
  CHECK(r.stats.min_b_impact_active >= -1e-6);
  CHECK(r.stats.min_res_impact_active >= -1e-3);
  // Thrust actuator/HOCBF-validity guards (HARD): 0 <= T <= T_static.
  CHECK(r.stats.min_T >= -1e-6);
  CHECK(r.stats.max_T <= 50.0 + 1e-6);
  // The trajectory stayed inside the stall-table validity range; past
  // kAlphaModelLimit any "recovery" is not physical.
  CHECK_FALSE(r.stats.out_of_model);
}
}  // namespace

TEST_CASE("CBF-on pull-into-stall recovers and lands safely",
          "[stall_recovery][integration]") {
  // Nominal commands an unsustainable +22 deg climb at 18 m and would depart
  // (that raw departure is the CBF-off demo); the filter must cap alpha and
  // land within the envelope with no recovery logic in the nominal.
  const RunResult r = runScenario("lon_stall_pull_cbf.yaml", "pull_cbf_on");
  REQUIRE(r.stats.ic_trim_converged);
  REQUIRE(r.td.reached);
  CHECK(r.td.sink < 2.5);
  CHECK(r.td.V < 14.5);                       // energy cap 14 + margin
  CHECK(r.td.theta > -5.0 * kDeg);
  // Soft stall row may dip past alpha_lim = 9 deg transiently, but the wing
  // must not depart: the CBF-off twin exceeds 30 deg here.
  CHECK(r.stats.max_alpha < 14.0 * kDeg);
  checkHardInvariants(r);
}

TEST_CASE("CBF-on stalled entry at 12 m arrests the departure and lands",
          "[stall_recovery][integration]") {
  // Starts already past stall (alpha0 ~ 22 deg) at 12 m; the entry excursion
  // is by design, so the alpha bound only excludes a full departure, and the
  // sink bound admits the deliberate hard-nose-down recovery landing the
  // design doc calls the honest tradeoff (~2.7 m/s cited).
  const RunResult r = runScenario("lon_stall_entry_cbf.yaml", "entry_cbf_on");
  REQUIRE(r.stats.ic_trim_converged);
  REQUIRE(r.td.reached);
  CHECK(r.td.sink < 3.5);
  CHECK(r.td.V < 14.5);
  CHECK(r.td.theta > -8.0 * kDeg);
  CHECK(r.stats.max_alpha < 45.0 * kDeg);
  checkHardInvariants(r);
}

TEST_CASE("Nominal alone departs -- the recovery is emergent, not commanded",
          "[stall_recovery][integration]") {
  // Premise guard: the existing CBF-off demo scenario genuinely stalls and
  // departs within its 8 s window. If a future nominal gains recovery logic,
  // this fails and the "emergent" claim of the two tests above needs revisit.
  const RunResult r = runScenario("lon_stall_recovery.yaml", "demo_cbf_off");
  CHECK_FALSE(r.td.reached);
  CHECK(r.stats.max_alpha > 12.0 * kDeg);     // past the 11 deg wing stall
}
