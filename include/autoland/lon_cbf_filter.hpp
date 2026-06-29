#pragma once
#include <array>
#include <memory>
#include "autoland/aero_table.hpp"
#include "autoland/config.hpp"
#include "autoland/hocbf.hpp"
#include "autoland/impact_barrier.hpp"
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
  bool airspeed_upper{true};
  bool thrust_limits{true};
  bool impact{true};
  bool descent_hard{true};
  bool airspeed_hard{false};
  bool airspeed_upper_hard{false};
  bool impact_hard{false};

  double v_safe{0.6};     // hull-safe touchdown sink rate [m/s]
  double a_brk{3.0};      // [DEPRECATED] constant braking accel. The descent
                          // barrier now computes a_brk(V,gamma) from CLmax; this
                          // is kept only for logging / plot-script compatibility.
  double CLmax{1.2};      // max lift coeff (stall ceiling) for a_brk(V,gamma).
                          // PLACEHOLDER -- calibrate to the real airframe.
  double Vmin{15.0};      // stall-margin airspeed [m/s]
  double Vmax_air{30.0};  // never-exceed airspeed [m/s] for the upper airspeed
                          // barrier (over-speed / high-energy impact). NOTE: this
                          // is airspeed -- distinct from Tmax (thrust). PLACEHOLDER.
  double Tmax{12.0};      // max thrust [N]

  std::array<double, 3> c_descent{2.0, 2.0, 2.0};   // class-K gains (deg 3)
  std::array<double, 3> c_airspeed{2.0, 2.0, 2.0};
  std::array<double, 3> c_airspeed_upper{2.0, 2.0, 2.0};
  std::array<double, 2> c_thrust_min{4.0, 4.0};     // {c11, c12} (deg 2)
  std::array<double, 2> c_thrust_max{4.0, 4.0};     // {c21, c22}

  // --- Impact-load barrier (NACA TN 1516; degree-2 HOCBF, elevator-enforced) --
  // Bounds the peak CG load factor at water touchdown. n_limit/beta/Nb/zs are
  // PLACEHOLDERS to calibrate -- see TODO.md. NOTE: n_limit (load factor) and
  // beta/rho_water are HYDRODYNAMIC; do not confuse with the aero CLmax above.
  double n_limit{3.0};        // structural CG load-factor limit [g] (normal to water)
  double beta{22.5 * M_PI / 180.0};  // hull dead-rise [rad]
  double rho_water{1000.0};   // water density [kg/m^3] (1000 fresh / 1025 sea)
  double Nb{10.0};            // Phi(z) budget [g]: counterfactual excess load tolerated
  double zs{2.0};             // Phi(z) altitude scale [m] (flare-authority height)
  double tau_keel{0.0};       // keel incidence: tau = theta - tau_keel [rad]
  double z_gate{10.0};        // assemble the impact row only below this height [m]
  double eps_g0{0.02};        // smooth floor on sin(gamma0) (planing-singularity guard)
  double impact_slack_lo{1.0e2};  // Option C height-scheduled slack penalty (high z)
  double impact_slack_hi{1.0e4};  //   "   (z -> 0): cheap to relax high, firm near water
  std::array<double, 2> c_impact{2.0, 2.0};  // class-K gains (deg 2)

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
