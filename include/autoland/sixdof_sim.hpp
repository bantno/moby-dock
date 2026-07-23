#pragma once
#include <memory>
#include <string>
#include "autoland/aero_table.hpp"
#include "autoland/config.hpp"
#include "autoland/dynamics.hpp"
#include "autoland/mixing.hpp"
#include "autoland/sixdof_nominal.hpp"
#include "autoland/trim.hpp"
#include "autoland/water_waves.hpp"
#include "autoland/wind_gust.hpp"

// =============================================================================
// 6-DOF straight-in water-landing simulation: the full nonlinear body-axis EOM
// (Dynamics::xdot, the repo's single source of truth) closed with the cascaded
// PID nominal (sixdof_nominal.hpp), against the plant-side wind gust and
// surface-wave models. RK4 with zero-order-hold controls; runs to touchdown
// (h <= eta(x, t); flat water when waves are disabled) or t_max. No CBF filter
// yet -- the nominal command is applied directly, but the loop is shaped so
// the filter drops in between step() and the integrator later.
// =============================================================================
namespace autoland {

// Hull/contact parameters for the TN 1516 slam-load truth reported at
// touchdown. Purely diagnostic in this sim (no impact barrier).
struct HullParams {
  double beta{22.5 * M_PI / 180.0};  // dead-rise angle [rad]
  double rho_water{1000.0};          // [kg/m^3]
  double tau_keel{0.0};              // keel incidence: tau = theta - tau_keel
  double eps_g0{0.02};               // smooth floor on sin(gamma0)
};

// Scenario bundle. Nominal feedforwards default to the trim solve at
// (V_app, gamma_app); the YAML may override any field. Wind and waves each
// carry their own `enabled` flag -- one YAML line toggles them.
struct SixDofScenario {
  SixDofNominalConfig nominal;
  DiscreteGustConfig wind;  // MIL-F-8785C discrete gust, plant-side only
  WaveConfig waves;         // Airy/JONSWAP surface waves, plant-side only
  HullParams hull;
  double V_app{18.0};
  double gamma_app{-3.0 * M_PI / 180.0};
  // Initial condition = trim state with these offsets applied.
  double h0{40.0};    // start altitude [m]
  double dV{0.0};     // airspeed offset [m/s]
  double dy{0.0};     // cross-track offset [m]
  double dpsi{0.0};   // heading offset [rad]
  double dtheta{0.0}; // pitch offset [rad]
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

 private:
  AeroTable table_;
  AircraftConfig ac_;
  std::unique_ptr<Mixing> mixing_;
  std::unique_ptr<Dynamics> dyn_;
  TrimResult trim_;
  SixDofScenario sc_;
  StateVec x0_{StateVec::Zero()};
  SixDofRunStats stats_;
};

}  // namespace autoland
