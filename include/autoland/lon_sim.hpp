#pragma once
#include <memory>
#include <string>
#include "autoland/aero_table.hpp"
#include "autoland/config.hpp"
#include "autoland/lon_cbf_filter.hpp"
#include "autoland/lon_nominal.hpp"
#include "autoland/mixing.hpp"
#include "autoland/water_waves.hpp"
#include "autoland/wind_gust.hpp"

// =============================================================================
// Self-consistent augmented-longitudinal landing simulation: the plant IS the
// CBF model (Xdot = f(X) + g(X) U*), so there is zero model mismatch -- a clean
// test of the CBF/QP math. RK4 integration with a zero-order-hold control.
// =============================================================================
namespace autoland {

// Scenario bundle. Defaults are derived from a trim solve at (V_app, gamma_app);
// the YAML may override any field. Held together here to avoid coupling the
// generic config.hpp to the longitudinal control headers.
struct LonScenario {
  LonNominalConfig nominal;
  LonCBFConfig cbf;
  DiscreteGustConfig wind;  // MIL-F-8785C discrete gust, plant-side only
  WaveConfig waves;         // Airy/JONSWAP surface-wave field, plant-side only
  LonStateVec X0{LonStateVec::Zero()};
  double dt{0.01};
  double t_max{60.0};
  bool cbf_enabled{true};
};

struct LonTouchdown {
  bool reached{false};
  double t{0}, sink{0}, V{0}, theta{0}, gamma{0};
  // Wave-field contact record (all zero on flat water): surface elevation /
  // slope under the aircraft at contact, keel altitude h and ground position x
  // there, the SURFACE-relative closure rate -d/dt(h - eta) (the face may be
  // rising into the path), and the exact TN 1516 peak-load truth evaluated
  // flat-water- vs wave-referenced (tau / gamma0 tilted by the local surface
  // angle, closure onto the moving surface). The filter saw only the
  // flat-water quantities, so n_peak_wave / n_peak_flat is the factor the
  // smooth-water assumption missed on this landing.
  double h{0}, x{0}, eta{0}, slope{0}, sink_rel{0};
  double n_peak_flat{0}, n_peak_wave{0};
  // Filter-health tallies over the whole run (also printed to the console):
  // best-effort feasibility recoveries and steps that dropped a HARD row. Both
  // are 0 on a healthy run; exposed so an end-to-end test can assert it.
  int recoveries{0};
  int steps_hard_dropped{0};
};

// Plant validity ceiling for the stall table: the NACA 4414 Viterna curve is
// held flat past ~90 deg and a fully-departed tumble there is not physical, so
// any run whose |alpha| exceeds this is flagged out-of-model (85 deg leaves a
// conservative margin before the clamp).
constexpr double kAlphaModelLimit = 85.0 * M_PI / 180.0;

// Full-run record for end-to-end tests: everything LonSim::run tallies for the
// console report (HOCBF psi minima over each row's ACTIVE window, filter-health
// counters) plus flight-envelope extrema and validity flags, exposed so tests
// can assert on them without parsing the CSV or stdout. Populated by run();
// valid via stats() afterwards. Minima keep the 1e30 "never active" sentinel.
struct LonRunStats {
  // CBF/QP health (recoveries/steps_hard_dropped mirror LonTouchdown).
  int recoveries{0};
  int steps_rows_dropped{0};
  int steps_hard_dropped{0};
  // HOCBF nested-function minima over each row's active window: stall/energy
  // whole flight, nose-up over h < h_noseup, energy-floor only when opt-in,
  // impact over (h < z_gate && descending && theta > tau_keel).
  double min_psi1_stall{1e30}, min_psi2_stall{1e30};
  double min_psi1_noseup{1e30}, min_psi2_noseup{1e30};
  double min_psi1_energy{1e30}, min_psi2_energy{1e30}, min_psi3_energy{1e30};
  double min_psi1_efloor{1e30}, min_psi2_efloor{1e30}, min_psi3_efloor{1e30};
  // NOTE: psi2 (and the row residual) involve the CONTROL; the drift-only psi2
  // minima here go legitimately negative exactly when the QP has to spend
  // authority. The state-only forward-invariance condition for the degree-2
  // impact row is min_psi1_impact >= 0; min_res_impact_active tracks the
  // enforced row residual (rhs - a.u, >= 0 up to QP tolerance) as the check
  // that the applied control actually satisfied the HARD constraint.
  double min_psi1_impact{1e30}, min_psi2_impact{1e30};
  double min_b_impact_active{1e30};    // b_impact itself over the same window
  double min_res_impact_active{1e30};  // enforced-row residual, same window
  int nan_stall{0}, nan_noseup{0}, nan_energy{0}, nan_efloor{0}, nan_impact{0};
  // Flight-envelope extrema [rad / SI] for safe-recovery verdicts.
  double max_alpha{-1e30};
  double min_V{1e30};
  double min_T{1e30}, max_T{-1e30};
  bool out_of_model{false};      // |alpha| ever exceeded kAlphaModelLimit
  bool ic_trim_converged{true};  // level-flight trim at V0 (constructor)
  int steps{0};
  double t_end{0.0};
};

class LonSim {
 public:
  LonSim(const std::string& stab_path, const std::string& aircraft_yaml,
         const std::string& scenario_yaml);

  // Run to touchdown (h <= eta(x,t); flat-water 0) or t_max. Writes a CSV
  // trace to csv_path.
  LonTouchdown run(const std::string& csv_path);

  // Full-run tallies (psi minima, envelope extrema, health counters); the
  // ic_trim_converged flag is set by the constructor, the rest by run().
  const LonRunStats& stats() const { return stats_; }

  const LonScenario& scenario() const { return sc_; }

 private:
  AeroTable table_;
  AircraftConfig ac_;
  std::unique_ptr<Mixing> mixing_;
  LonScenario sc_;
  LonRunStats stats_;
};

// Full augmented RHS with aero rebuilt from the state (plant + control).
LonStateVec lonXdotFull(const AeroTable& table, const Mixing& mixing,
                        const AircraftConfig& cfg, const LonStateVec& X,
                        const LonCtrlVec& U);

// Wind-perturbed plant RHS. The state keeps its air-relative meaning (V =
// airspeed, gamma = air-path angle, so alpha = theta - gamma is still the true
// aerodynamic AoA); the earth-frame gust W = (u tailwind+, w updraft+) enters
// exactly, via
//   hdot     += W.w                                     (kinematic transport)
//   Vdot     += -( Wdot.u cos(gamma) + Wdot.w sin(gamma) )
//   gammadot += ( Wdot.u sin(gamma) - Wdot.w cos(gamma) ) / V
// i.e. the -m Wdot inertial forcing of m d(v_air)/dt = F - m Wdot projected on
// the air-path axes. A steady wind (Wdot = 0) therefore only transports the
// aircraft; only the gust RAMP forces V/gamma -- the classic airspeed loss in
// a tailwind ramp and the AoA rise in an updraft ramp.
LonStateVec lonXdotFullWind(const AeroTable& table, const Mixing& mixing,
                            const AircraftConfig& cfg, const LonStateVec& X,
                            const LonCtrlVec& U, const GustWind& W,
                            const GustWind& Wdot);

}  // namespace autoland
