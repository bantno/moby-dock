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
  double h_flare{0.0};    // flare height [m]: the descent barrier levels the sink
                          // toward v_safe at h = h_flare instead of at the surface
                          // (a knob on where the flare begins). 0 => flare at h=0.
                          // IMPLEMENTED but needs tuning -- see TODO / CHANGELOG.
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
  std::array<double, 2> c_impact{2.0, 2.0};  // class-K gains (deg 2)

  // --- QP objective ----------------------------------------------------------
  // J = 1/2 ||u - u_nom||^2  +  sum_i 1/2 w_slack_i * delta_i^2
  // The control term is UNWEIGHTED (identity) -- it only regularizes u toward the
  // nominal. Each SOFT barrier i carries a free slack delta_i >= 0 (constraint
  // a_i.u - delta_i <= rhs_i) penalized QUADRATICALLY by w_slack_i: larger w =>
  // firmer constraint. (Quadratic, not linear w_i*delta_i: a pure-linear penalty
  // leaves P singular on the slacks and OSQP then fails to converge on many
  // feasible steps, firing best-effort spuriously.) HARD barriers (the *_hard
  // flags, and the thrust-limit rows) get NO slack; if the hard set is infeasible
  // the filter falls back to a best-effort minimum-violation solve.
  double w_slack_descent{1.0e4};
  double w_slack_airspeed{1.0e4};
  double w_slack_airspeed_upper{1.0e4};
  double w_slack_impact{1.0e4};

  double de_min{-0.5}, de_max{0.5};            // elevator bounds [rad]
  double Tddot_min{-500.0}, Tddot_max{500.0};  // Tddot bounds [N/s^2]

  // Altitude-measurement noise standard deviation [m]. The sim feeds the filter
  // a noisy AGL measurement h_meas = h + N(0, h_meas_stddev^2); the QP math is
  // unchanged and the plant still integrates the TRUE h (the noise is a sensor
  // model applied in LonSim::run, not in the barrier math). 0 => perfect
  // altitude (default). Tunable; ~1 m models a noisy AGL/radar-altimeter sensor.
  double h_meas_stddev{0.0};

  // RNG seed for the altitude-measurement noise. Exposed so a Monte-Carlo sweep
  // can evaluate a gain choice across many independent noise realizations rather
  // than a single fixed draw. Fixed by default for run-to-run reproducibility.
  unsigned int h_meas_seed{0xA17B0A11u};

  // First-order low-pass time constant [s] on the altitude measurement, applied
  // before the CBF sees it (h_filt += alpha*(h_meas - h_filt), alpha=dt/(tau+dt)).
  // Attenuates the sensor noise that drives near-surface barrier/elevator chatter,
  // at the cost of ~tau lag. 0 => no filtering / pass-through (default). Cutoff
  // f_c = 1/(2*pi*tau); e.g. tau=0.2 s ~ 0.8 Hz.
  double h_lpf_tau{0.0};
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
  // Per-call row-assembly health, set by the most recent filter() call.
  //   lastDroppedRows()      : barrier rows that came out non-finite and were
  //                            dropped from the QP (coeffs zeroed, rhs -> +inf).
  //   lastHardDropped()      : at least one dropped row was a HARD constraint --
  //                            an unannunciated loss of a safety guarantee.
  //   lastDescentInfeasible(): a_brk(V,gamma) <= 0 -- the soft-landing envelope
  //                            had no real braking solution (root-cause signal).
  int lastDroppedRows() const { return dropped_rows_; }
  bool lastHardDropped() const { return hard_dropped_; }
  bool lastDescentInfeasible() const { return descent_infeasible_; }
  const LonCBFConfig& config() const { return cfg_; }

 private:
  LonCBFConfig cfg_;
  std::unique_ptr<QPSolver> solver_;
  mutable bool recovered_{false};
  mutable int dropped_rows_{0};
  mutable bool hard_dropped_{false};
  mutable bool descent_infeasible_{false};
};

}  // namespace autoland
