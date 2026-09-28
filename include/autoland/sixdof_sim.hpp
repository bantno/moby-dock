#pragma once
#include <functional>
#include <memory>
#include <string>
#include "autoland/aero_table.hpp"
#include "autoland/beaver_dynamics.hpp"
#include "autoland/config.hpp"
#include "autoland/dynamics.hpp"
#include "autoland/energy_reach_cbf.hpp"
#include "autoland/mixing.hpp"
#include "autoland/pattern_guidance.hpp"
#include "autoland/sixdof_cbf.hpp"
#include "autoland/sixdof_nominal.hpp"
#include "autoland/trim.hpp"
#include "autoland/water_rollout.hpp"
#include "autoland/water_waves.hpp"
#include "autoland/world_geometry.hpp"
#include "autoland/wind_gust.hpp"

// =============================================================================
// 6-DOF straight-in water-landing simulation: a full nonlinear body-axis EOM
// plant closed with the cascaded PID nominal (sixdof_nominal.hpp), against the
// plant-side wind gust and surface-wave models. RK4 with zero-order-hold
// controls; runs to touchdown (h <= eta(x, t); flat water when waves are
// disabled) or t_max. No CBF filter yet -- the nominal command is applied
// directly, but the loop is shaped so the filter drops in between step() and
// the integrator later.
//
// PLANT SELECTION (scenario `plant:` key):
//   * "beaver" (DEFAULT): the flight-validated DHC-2 Beaver polynomial model
//     (beaver_dynamics.hpp, LR-556/FDC, validated against the FDC check case
//     -- see apps/beaver_validation.cpp). Optional `beaver:` block sets
//     n_rpm / pz_idle / pz_max / flap_deg / h_ref.
//   * "vspaero": the original AHAB VSPAERO-table plant (Dynamics::xdot); needs
//     the stab/aircraft.yaml paths passed to the constructor.
// =============================================================================
namespace autoland {

// Hull/contact parameters for the TN 1516 slam-load truth reported at
// touchdown. Purely diagnostic in this sim (no impact barrier).
struct HullParams {
  double beta{22.5 * M_PI / 180.0};  // dead-rise angle [rad]
  double rho_water{1000.0};          // [kg/m^3]
  double tau_keel{0.0};              // keel incidence: tau = theta - tau_keel
  double eps_g0{0.02};               // smooth floor on sin(gamma0)
  int n_surfaces{1};                 // planing surfaces: 1 hull / 2 floats -> W/n
};

// Energy-reach CBF layer configuration (scenario `energy_cbf:` block). When
// enabled the sim runs the nominal OUTER loop, maps theta_cmd into the
// point-mass (alpha_dot, beta) space, filters through EnergyReachFilter, and
// closes the INNER loops on the filtered commands (wings-level + rudder
// sideslip lateral mode). `filter` = false keeps the whole layer (idle
// throttle, beta lateral, h_E logging against the same target) but passes the
// nominal through unfiltered -- the A/B baseline.
struct ReachSimCfg {
  bool enabled{false};
  bool filter{true};        // false -> monitor-only baseline
  bool throttle_idle{true}; // force dT = dT_min (glide, matches the model)
  double tau_alpha{1.0};    // alpha_dot_nom = (alpha_des - alpha_cmd)/tau [s]
  double target_N{NAN};     // NaN -> glideslope aim point h0/tan(-gamma_app)
  double target_E{NAN};     // NaN -> 0 (on the centerline)
  EnergyReachCfg cbf;       // filter config; prm's physical/polar/target
                            // fields are filled by the sim at construction
};

// Scenario bundle. Nominal feedforwards default to the trim solve at
// (V_app, gamma_app); the YAML may override any field. Wind and waves each
// carry their own `enabled` flag -- one YAML line toggles them.
//
// `world:` (corridor + terrain + rollout, world_geometry.hpp) switches the
// sim from the legacy straight-in (centerline along north through the
// origin, run ends at touchdown) to the corridor landing: the pattern
// guidance (pattern_guidance.hpp) plans a Dubins path to the corridor's
// final-approach fix and tracks it, touchdown is tested against
// max(waves, terrain), and the water rollout runs to a stop. The initial
// condition is then ABSOLUTE (N0, E0, psi0) instead of centerline offsets.
struct SixDofScenario {
  SixDofNominalConfig nominal;
  DiscreteGustConfig wind;  // MIL-F-8785C discrete gust, plant-side only
  WaveConfig waves;         // Airy/JONSWAP surface waves, plant-side only
  HullParams hull;
  ReachSimCfg reach;        // reach-the-target energy CBF layer
  WorldConfig world;        // corridor / terrain / rollout (corridor landing)
  PatternGuidanceConfig guidance;  // pattern nominal (world.enabled only)
  SixDofCbfConfig cbf;      // 6-DOF surfaces-only CBF-QP filter (`sixdof_cbf:`)
  // Envelope verdict: any step with alpha above this is a FAILED run (the
  // plant has no post-stall aero; numbers past it are not physical).
  double alpha_stall{16.0 * M_PI / 180.0};
  double V_app{18.0};
  double gamma_app{-3.0 * M_PI / 180.0};
  // Initial condition = trim state with these offsets applied.
  double h0{40.0};    // start altitude [m]
  double dV{0.0};     // airspeed offset [m/s]
  double dy{0.0};     // cross-track offset [m]
  double dpsi{0.0};   // heading offset [rad]
  double dtheta{0.0}; // pitch offset [rad]
  // world.enabled: absolute start position / heading (dy, dpsi ignored).
  double N0{0.0}, E0{0.0}, psi0{0.0};
  double dt{0.01};
  double t_max{200.0};
};

struct SixDofTouchdown {
  bool reached{false};
  double t{0};
  double sink{0};   // inertial descent rate, positive down [m/s]
  double V{0};      // airspeed [m/s]
  double gamma{0};  // inertial flight-path angle [rad]
  double alpha{0}, beta{0};          // air-relative [rad]
  double theta{0}, phi{0}, psi{0};   // attitude [rad]
  double y{0};                       // cross-track [m]
  // Wave-field contact record (zero on flat water) -- mirrors LonTouchdown:
  // surface elevation/slope under the keel, ground position x, the
  // surface-relative closure rate, and the TN 1516 peak-load truth evaluated
  // flat- vs wave-referenced.
  double x{0}, h{0}, eta{0}, slope{0}, sink_rel{0};
  double n_peak_flat{0}, n_peak_wave{0};
  // Energy-reach layer (zero unless the scenario enables it): barrier value
  // and NE distance to the target at touchdown.
  double hE{0}, dist_target{0};
  // Corridor landing (world.enabled): terrain height under the contact
  // point, whether contact was terrain (crash) rather than water, corridor
  // frame position / heading error / crab at touchdown, and the rollout.
  double hT{0};
  bool terrain_strike{false};
  bool in_corridor{false};
  double s_corr{0}, e_corr{0};   // corridor frame [m]
  double dpsi_corr{0};           // psi - psi_c [rad] (heading; crab-dominated in wind)
  double dchi_corr{0};           // chi - psi_c [rad] (ground course: the alignment metric)
  double crab{0};                // psi - ground course [rad]
  RolloutResult rollout;
  // Verdict (world.enabled): water contact inside the corridor, rollout
  // stopped inside, never a terrain strike, never in the stall region.
  bool success{false};
  std::string fail_reasons;      // "" when success
};

// Full-run tallies for end-to-end tests (no CSV/stdout parsing).
struct SixDofRunStats {
  bool trim_converged{true};  // approach trim at (V_app, gamma_app)
  int steps{0};
  double t_end{0.0};
  double max_abs_phi{0.0};   // [rad]
  double max_abs_beta{0.0};  // air-relative [rad]
  double max_alpha{-1e30};   // air-relative [rad]
  double min_V{1e30};        // airspeed [m/s]
  double max_abs_y{0.0};     // [m]
  // Energy-reach layer (meaningful only when the scenario enables it).
  double min_hE{1e30};       // barrier minimum over the run [m^2/s^2]
  int cbf_active_steps{0};   // steps where the filter moved the command
  int cbf_best_effort_steps{0};  // hard row infeasible -> max-hdot fallback
  // Corridor landing (world.enabled).
  double min_terrain_clearance{1e30};  // min over the run of h - h_T(N,E) [m] (CG)
  double min_wingtip_clearance{1e30};  // ... minus (b/2)|sin phi| for the low wingtip [m]
  double max_abs_e_final{0.0};         // max |e| while on final [m]
  double max_xte{0.0};                 // max cross-track distance to the planned path [m]
  int phi_sat_steps{0};                // steps with |phi_cmd| at the guidance bank limit
  double plan_length{0.0};             // Dubins + final leg to the aim point [m]
  double d_final{0.0};                 // FAF-to-aim distance [m]
  std::string plan_word;               // Dubins word
  // Envelope verdict.
  bool stall_entered{false};
  int stall_steps{0};
  // 6-DOF CBF filter (sixdof_cbf.enabled).
  int cbf6_active_steps{0};
  int cbf6_best_effort_steps{0};
  double cbf6_min_h_alpha{1e30}, cbf6_min_h_V{1e30}, cbf6_min_h_phi{1e30}, cbf6_min_h_beta{1e30};
  double cbf6_min_h_terrain{1e30};     // only while the terrain row is gated in
  double cbf6_margin_max{0.0};
  double cbf6_time_mean_us{0.0}, cbf6_time_max_us{0.0};
};

class SixDofSim {
 public:
  SixDofSim(const std::string& stab_path, const std::string& aircraft_yaml,
            const std::string& scenario_yaml);

  // Run to touchdown (h <= eta) or t_max. Writes a CSV trace to csv_path.
  SixDofTouchdown run(const std::string& csv_path);

  const SixDofRunStats& stats() const { return stats_; }
  const SixDofScenario& scenario() const { return sc_; }
  const TrimResult& trimResult() const { return trim_; }
  // Idle-glide polar fitted at construction (energy_cbf scenarios only).
  const ReachPolarFit& reachPolar() const { return polar_; }
  // Initial state (after the scenario offsets / absolute placement).
  const StateVec& initialState() const { return x0_; }

  // Plant state derivative (wind-aware); bound to whichever plant is active.
  using PlantFn =
      std::function<StateVec(const StateVec&, const CtrlVec&,
                             const Eigen::Vector3d&)>;

 private:
  std::unique_ptr<AeroTable> table_;   // vspaero plant only
  AircraftConfig ac_;                  // mass/env/limits for diagnostics
  std::unique_ptr<Mixing> mixing_;     // vspaero plant only
  std::unique_ptr<Dynamics> dyn_;      // vspaero plant only
  std::unique_ptr<BeaverDynamics> bdyn_;  // beaver plant only
  PlantFn plant_;
  TrimResult trim_;
  SixDofScenario sc_;
  StateVec x0_{StateVec::Zero()};
  SixDofRunStats stats_;
  ReachPolarFit polar_{};
};

}  // namespace autoland
