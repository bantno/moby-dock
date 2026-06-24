#pragma once
#include <memory>
#include <string>
#include <vector>
#include "autoland/aero_table.hpp"
#include "autoland/cbf.hpp"
#include "autoland/config.hpp"
#include "autoland/controller.hpp"
#include "autoland/dynamics.hpp"
#include "autoland/linear_model.hpp"
#include "autoland/mixing.hpp"
#include "autoland/trim.hpp"

// =============================================================================
// Closed-loop autoland simulation driver.
//
// Pipeline: load config + .stab -> trim -> linearize about trim -> initialize
// on a stabilized approach a configurable distance out -> run the closed loop
// to touchdown (h crossing zero) -> log states/controls/references to CSV.
//
// The closed loop currently runs against the LINEAR plant (LinearModel::xdot),
// but the nonlinear EOM is the foundation: swap plant_ to Dynamics::xdot and the
// same loop runs the full nonlinear plant with no other change.
// =============================================================================
namespace autoland {

struct TouchdownSummary {
  bool reached{false};
  double t{0};
  double sink_rate{0};  // positive down [m/s]
  double V{0};          // [m/s]
  double theta{0};      // [rad]
  double phi{0};        // [rad]
  double beta{0};       // [rad]
  double y{0};          // cross-track [m]
};

class Sim {
 public:
  Sim(const std::string& stab_path, const std::string& aircraft_yaml,
      const std::string& scenario_yaml);

  // Run the closed loop, writing the CSV log to csv_path. Returns the summary.
  TouchdownSummary run(const std::string& csv_path);

  const TrimResult& trimResult() const { return trim_; }
  const LinearModel& linearModel() const { return lin_; }

 private:
  AeroTable table_;
  AircraftConfig ac_;
  ScenarioConfig sc_;
  std::unique_ptr<Mixing> mixing_;
  std::unique_ptr<Dynamics> dyn_;
  TrimResult trim_;
  LinearModel lin_;
  CBFFilter cbf_;

  // Candidate CBF barriers (currently passed to the pass-through filter).
  std::vector<Barrier> barriers_;
  void buildCandidateBarriers();
};

}  // namespace autoland
