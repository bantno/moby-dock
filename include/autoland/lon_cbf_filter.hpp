#pragma once
#include <array>
#include <memory>
#include "autoland/aero_table.hpp"
#include "autoland/config.hpp"
#include "autoland/hocbf.hpp"
#include "autoland/lon_augmented.hpp"
#include "autoland/mixing.hpp"
#include "autoland/qp_solver.hpp"

// =============================================================================
// Augmented-longitudinal CBF-QP safety filter (documentation section 4).
//
//   U* = argmin 1/2 ||U - U_nom||_W^2   s.t.  the HOCBF + actuator rows
//
// Barriers: descent-rate (deg 3), airspeed (deg 3), min/max thrust (deg 2).
// The contact-force barrier (section 3.3) is intentionally deferred. The QP is
// solved with OSQP; hard rows are honored first, and if that set is infeasible
// the filter recovers by softening every row for that step.
// =============================================================================
namespace autoland {

struct LonCBFConfig {
  bool enabled{true};
  bool descent{true};
  bool airspeed{true};
  bool thrust_limits{true};
  bool descent_hard{true};
  bool airspeed_hard{false};

  double v_safe{0.6};     // hull-safe touchdown sink rate [m/s]
  double a_brk{3.0};      // braking deceleration for the descent envelope [m/s^2]
  double Vmin{15.0};      // stall-margin airspeed [m/s]
  double Tmax{12.0};      // max thrust [N]

  std::array<double, 3> c_descent{2.0, 2.0, 2.0};   // class-K gains (deg 3)
  std::array<double, 3> c_airspeed{2.0, 2.0, 2.0};
  std::array<double, 2> c_thrust_min{4.0, 4.0};     // {c11, c12} (deg 2)
  std::array<double, 2> c_thrust_max{4.0, 4.0};     // {c21, c22}

  double w_de{1.0};        // QP weight on elevator deviation
  double w_Tddot{1.0};     // QP weight on Tddot deviation
  double slack_penalty{1.0e4};

  double de_min{-0.5}, de_max{0.5};            // elevator bounds [rad]
  double Tddot_min{-500.0}, Tddot_max{500.0};  // Tddot bounds [N/s^2]
};

class LonCBFFilter {
 public:
  LonCBFFilter() : solver_(std::make_unique<OsqpSolver>()) {}
  explicit LonCBFFilter(LonCBFConfig cfg,
                        std::unique_ptr<QPSolver> solver = nullptr)
      : cfg_(cfg),
        solver_(solver ? std::move(solver)
                       : std::unique_ptr<QPSolver>(std::make_unique<OsqpSolver>())) {}

  // Filter U_nom at state X. Aero is sourced from (table, mixing, cfg).
  LonCtrlVec filter(const LonCtrlVec& U_nom, const LonStateVec& X,
                    const AeroTable& table, const Mixing& mixing,
                    const AircraftConfig& cfg) const;

  bool enabled() const { return cfg_.enabled; }
  bool lastRecovery() const { return recovered_; }
  const LonCBFConfig& config() const { return cfg_; }

 private:
  LonCBFConfig cfg_;
  std::unique_ptr<QPSolver> solver_;
  mutable bool recovered_{false};
};

}  // namespace autoland
