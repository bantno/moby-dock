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
// Barriers: impact-load (deg 2, the only HARD safety row), stall/AoA (deg 2),
// nose-up attitude floor (deg 2, gated near the surface), total-energy ceiling
// (deg 3), and min/max thrust actuator guards (deg 2, HARD). There is NO airspeed
// floor -- the AoA barrier guards stall directly; the thrust CBFs are the
// actuator-effectiveness / HOCBF-validity guards. The QP is solved with OSQP; the
// hard rows are honored first, and if that set is infeasible the filter falls back
// to a best-effort minimum-violation solve.
// =============================================================================
namespace autoland {

struct LonCBFConfig {
  bool enabled{true};
  bool thrust_limits{true};
  bool impact{true};
  bool stall{true};
  bool noseup{true};
  bool energy{true};
  bool impact_hard{true};    // impact is the only hard SAFETY-envelope barrier
  bool stall_hard{false};    // stall / nose-up / energy are all SOFT (slacked)
  bool noseup_hard{false};
  bool energy_hard{false};
  // Energy FLOOR (powered stall recovery; see hocbf.hpp). OPT-IN: default off
  // preserves the recovery set's no-airspeed-floor design.
  bool energy_floor{false};
  bool energy_floor_hard{false};

  // Stall / nose-up / energy-corridor parameters. Angles stored in RADIANS (the
  // YAML gives degrees). The AoA barrier replaces the old airspeed floor as the
  // stall guard; there is deliberately NO airspeed-floor CBF (see class doc).
  double alpha_stall{11.0 * M_PI / 180.0};  // NACA-4414 wing-stall AoA [rad]
  double stall_margin{2.0 * M_PI / 180.0};  // AoA buffer below stall [rad]
  double theta_min{3.0 * M_PI / 180.0};     // nose-up attitude floor [rad]
                                            // (keep >= tau_keel for impact validity)
  double h_noseup{3.0};   // nose-up barrier active only below this height [m]
  double V_td_max{14.0};  // never-exceed touchdown airspeed [m/s]; E_cap(0) =
                          // 1/2 V_td_max^2. From the hull/structural limit
                          // (~1.4 V_stall), NOT the approach speed.
  double g_eff{16.0};     // energy-ceiling loosening rate [m/s^2]. LOOSE cap:
                          //   V <= sqrt(V_td_max^2 + 2(g_eff-g)h).
                          // Size so E_cap(h0) >= E0 (starts satisfied):
                          //   g_eff >= g + (V0^2 - V_td_max^2)/(2 h0).
  double Tmax{12.0};      // max thrust [N]
  double V_floor{11.0};   // energy-floor airspeed at h = 0 [m/s] (~1.1 V_stall)
  double g_floor{9.80665};  // floor slope [m/s^2]: = g -> pure airspeed floor
                            // V >= V_floor everywhere; < g -> relaxes aloft as
                            // V >= sqrt(V_floor^2 + 2(g_floor - g)h).

  std::array<double, 2> c_stall{4.0, 4.0};          // class-K gains (deg 2)
  std::array<double, 2> c_noseup{2.0, 2.0};         // (deg 2, gated near surface)
  std::array<double, 3> c_energy{2.0, 2.0, 2.0};    // (deg 3)
  std::array<double, 3> c_energy_floor{2.0, 2.0, 2.0};  // (deg 3, opt-in floor)
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
  // Per-control tracking weights, J_u = 1/2 [w_de (de - de_nom)^2 +
  // w_Tddot (Tddot - Tddot_nom)^2]. DEFAULT IDENTITY (backward compatible).
  // The identity cost prices 1 rad of elevator equal to 1 N/s^2 of Tddot, so on
  // any mixed-authority row the QP corrects with the elevator alone (its row
  // coefficients run ~1e3-1e4x thrust's) and thrust only moves once the
  // elevator saturates -- late and slammed. Range-scaled weights (w = 1/u_max^2:
  // w_de ~ 4, w_Tddot ~ 4e-6) price a full-scale deflection of each actuator
  // equally, letting thrust participate early and smoothly (the powered half of
  // the stall recovery, with the opt-in energy floor).
  double w_de{1.0}, w_Tddot{1.0};

  double w_slack_stall{1.0e5};   // firmest soft row (protect the stall recovery)
  double w_slack_energy{1.0e4};  // medium
  double w_slack_energy_floor{1.0e4};  // same tier as the ceiling (opt-in row)
  double w_slack_noseup{1.0e3};  // softest: yields to stall on the shared elevator
  double w_slack_impact{1.0e4};  // unused while impact_hard=true; kept for testing

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
  int lastDroppedRows() const { return dropped_rows_; }
  bool lastHardDropped() const { return hard_dropped_; }
  const LonCBFConfig& config() const { return cfg_; }

 private:
  LonCBFConfig cfg_;
  std::unique_ptr<QPSolver> solver_;
  mutable bool recovered_{false};
  mutable int dropped_rows_{0};
  mutable bool hard_dropped_{false};
};

}  // namespace autoland
