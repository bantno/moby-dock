#include "autoland/lon_cbf_filter.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace autoland {
namespace {

constexpr double kInf = 1.0e30;  // OSQP +/-inf sentinel
constexpr double kCap = 1.0e20;  // finite cap kept strictly below kInf

// Best-effort fallback weight. When the hard-constrained QP is infeasible we do
// NOT relax the hard rows to regain feasibility -- that would simply delete the
// constraint's effect. Instead we solve a minimum-violation QP whose hard-row
// slack (= the constraint violation) is penalized far above the soft rows and
// the nominal-tracking term, so the optimizer spends all available actuator
// authority driving the hard-barrier violation toward zero. The result is the
// best control we can apply to push the state back toward the safe set, even
// though no control fully satisfies the (infeasible) hard constraints.
constexpr double kBestEffortHardPenalty = 1.0e6;

double clampd(double v, double lo, double hi) { return std::max(lo, std::min(hi, v)); }

// One assembled barrier row in QP <= form: a . U <= rhs (with optional slack).
struct QRow {
  Eigen::Matrix<double, 1, NUA> a;
  double rhs{0};
  bool hard{true};
  bool finite{true};
  double slack_weight{0};  // per-row quadratic slack penalty w_i (soft rows only)
};

}  // namespace

LonCtrlVec LonCBFFilter::filter(const LonCtrlVec& U_nom, const LonStateVec& X,
                                const AeroTable& table, const Mixing& mixing,
                                const AircraftConfig& cfg) const {
  recovered_ = false;
  dropped_rows_ = 0;
  hard_dropped_ = false;
  if (!cfg_.enabled) return U_nom;

  const double V = X[LV];
  const double alpha = X[LTH] - X[LGAM];
  const AeroLocal aero = makeAeroLocal(table, mixing, cfg, V, alpha);

  // ---- Assemble the barrier rows -------------------------------------------
  std::vector<QRow> rows;
  auto pushHocbf = [&](const HocbfRow& h, bool hard, double slack_w) {
    QRow r;
    r.a = h.a;
    r.rhs = h.rhs;
    r.hard = hard;
    r.finite = std::isfinite(h.rhs) && std::isfinite(h.a(0, LDE)) &&
               std::isfinite(h.a(0, LTDDOT));
    r.slack_weight = slack_w;  // ignored for hard rows (no slack column)
    rows.push_back(r);
  };

  if (cfg_.stall) {
    const StallBarrier b = makeStallBarrier(cfg_.alpha_stall, cfg_.stall_margin);
    auto lie = barrierLie<2>(aero, b, X);
    std::vector<double> Lf(lie.Lf.begin(), lie.Lf.end());
    const std::vector<double> c(cfg_.c_stall.begin(), cfg_.c_stall.end());
    pushHocbf(hocbfRow(Lf, lie.LgLf, c), cfg_.stall_hard, cfg_.w_slack_stall);
  }
  if (cfg_.energy) {
    const EnergyBarrier b = makeEnergyBarrier(aero, cfg_.V_td_max, cfg_.g_eff);
    auto lie = barrierLie<3>(aero, b, X);
    std::vector<double> Lf(lie.Lf.begin(), lie.Lf.end());
    const std::vector<double> c(cfg_.c_energy.begin(), cfg_.c_energy.end());
    pushHocbf(hocbfRow(Lf, lie.LgLf, c), cfg_.energy_hard, cfg_.w_slack_energy);
  }
  // Nose-up attitude floor, gated to the final h_noseup metres (mirrors the impact
  // z-gate below). Keeps theta above tau_keel so the impact model stays valid and
  // the touchdown is nose-up; SOFT with the lowest weight, so it yields to the
  // stall guard on the shared elevator (the recovery behavior).
  if (cfg_.noseup && X[LH] < cfg_.h_noseup) {
    const NoseUpBarrier b = makeNoseUpBarrier(cfg_.theta_min);
    auto lie = barrierLie<2>(aero, b, X);
    std::vector<double> Lf(lie.Lf.begin(), lie.Lf.end());
    const std::vector<double> c(cfg_.c_noseup.begin(), cfg_.c_noseup.end());
    pushHocbf(hocbfRow(Lf, lie.LgLf, c), cfg_.noseup_hard, cfg_.w_slack_noseup);
  }
  // Hydrodynamic impact-load barrier (degree-2 HOCBF; elevator-enforced). Gated
  // on (a) height for efficiency -- inactive above z_gate where Phi(z) ~ Nb -- and
  // (b) model validity: NACA TN 1516 assumes a descending, positive-trim contact
  // (gamma0 > 0, tau > 0). Outside that the prediction (kappa < 0, n_peak < 0) is
  // meaningless, so the row is skipped rather than fed to the QP as a pathological
  // (huge-coefficient) constraint. In a normal landing tau/gamma0 stay positive
  // through the flare, so this coincides with the active window.
  const double tau_gate = X[LTH] - cfg_.tau_keel;
  if (cfg_.impact && X[LH] < cfg_.z_gate && -X[LGAM] > 0.0 && tau_gate > 0.0) {
    const ImpactLoadBarrier b =
        makeImpactLoadBarrier(aero, cfg_.n_limit, cfg_.beta, cfg_.rho_water,
                              cfg_.Nb, cfg_.zs, cfg_.tau_keel, cfg_.eps_g0, X);
    auto lie = barrierLie<2>(aero, b, X);  // relative degree 2 (elevator)
    std::vector<double> Lf(lie.Lf.begin(), lie.Lf.end());
    const std::vector<double> c(cfg_.c_impact.begin(), cfg_.c_impact.end());
    pushHocbf(hocbfRow(Lf, lie.LgLf, c), cfg_.impact_hard, cfg_.w_slack_impact);
  }
  if (cfg_.thrust_limits) {
    pushHocbf(thrustMinRow(X, cfg_.c_thrust_min[0], cfg_.c_thrust_min[1]), true, 0.0);
    pushHocbf(thrustMaxRow(X, cfg_.Tmax, cfg_.c_thrust_max[0], cfg_.c_thrust_max[1]), true,
              0.0);
  }

  // Row-assembly health: any barrier whose Lie row came out non-finite is dropped
  // from the QP below (coefficients zeroed, rhs -> +inf). A dropped HARD row is a
  // lost guarantee, so record it for the caller rather than failing silently.
  for (const QRow& r : rows)
    if (!r.finite) { ++dropped_rows_; if (r.hard) hard_dropped_ = true; }

  const int nb = static_cast<int>(rows.size());
  if (nb == 0) return U_nom;

  // ---- Assemble + solve the QP --------------------------------------------
  //   z = [de, Tddot, slack(one per slacked row)]
  // best_effort=false (primary): hard rows are HARD (no slack), soft rows carry
  //   penalized slack, and the objective tracks U_nom. Used whenever feasible, so
  //   the hard barriers are satisfied exactly with minimal deviation from nominal.
  // best_effort=true (infeasibility fallback): every row gets a slack so the QP
  //   is always solvable, but each HARD row's slack (its constraint violation) is
  //   penalized at kBestEffortHardPenalty -- far above the soft rows and tracking
  //   -- so the solution is the actuator-box-limited control that MINIMIZES the
  //   hard-barrier violation. We never soften a hard constraint to dodge it; we
  //   apply the most-safety-driving control the actuators allow.
  auto solveWith = [&](bool best_effort) -> QPResult {
    std::vector<int> slack_col(nb, -1);
    int nslack = 0;
    for (int i = 0; i < nb; ++i)
      if (best_effort || !rows[i].hard) slack_col[i] = NUA + nslack++;

    const int n = NUA + nslack;
    const int m = nb + NUA + nslack;  // barrier rows + control box + slack >= 0

    Mat P = Mat::Zero(n, n);
    Vec q = Vec::Zero(n);
    // Control: identity-weighted tracking, 1/2 ||u - u_nom||^2 (regularizer only,
    // no per-control weights). In the best-effort solve this is dominated by the
    // hard-violation penalty below, so it only breaks ties among equally-safe
    // controls (it never trades away safety).
    P(LDE, LDE) = 1.0;
    P(LTDDOT, LTDDOT) = 1.0;
    q[LDE] = -U_nom[LDE];
    q[LTDDOT] = -U_nom[LTDDOT];
    // Soft constraints: per-constraint QUADRATIC slack penalty 1/2 w_i delta_i^2
    // (P-term; q stays 0 on slacks). Larger w_i => firmer constraint. In the
    // best-effort solve, hard rows get a slack at the dominant kBestEffortHardPenalty
    // so the solution MINIMIZES their violation.
    for (int i = 0; i < nb; ++i)
      if (slack_col[i] >= 0)
        P(slack_col[i], slack_col[i]) =
            (best_effort && rows[i].hard) ? kBestEffortHardPenalty : rows[i].slack_weight;

    Mat A = Mat::Zero(m, n);
    Vec lo = Vec::Zero(m);
    Vec hi = Vec::Zero(m);
    int row = 0;

    // (1) barrier rows:  a . U (- slack) <= rhs   [lo = -inf, hi = rhs]
    for (int i = 0; i < nb; ++i) {
      const QRow& r = rows[i];
      A(row, LDE) = r.finite ? r.a(0, LDE) : 0.0;
      A(row, LTDDOT) = r.finite ? r.a(0, LTDDOT) : 0.0;
      if (slack_col[i] >= 0) A(row, slack_col[i]) = -1.0;
      lo[row] = -kInf;
      hi[row] = r.finite ? clampd(r.rhs, -kCap, kCap) : kInf;
      ++row;
    }
    // (2) control box.
    A(row, LDE) = 1.0; lo[row] = cfg_.de_min; hi[row] = cfg_.de_max; ++row;
    A(row, LTDDOT) = 1.0; lo[row] = cfg_.Tddot_min; hi[row] = cfg_.Tddot_max; ++row;
    // (3) non-negative slacks.
    for (int i = 0; i < nb; ++i)
      if (slack_col[i] >= 0) { A(row, slack_col[i]) = 1.0; lo[row] = 0.0; hi[row] = kInf; ++row; }

    return solver_->solve(P, q, A, lo, hi);
  };

  bool any_hard = false;
  for (const QRow& r : rows) any_hard = any_hard || r.hard;

  // Primary: enforce the hard rows exactly. If that is infeasible, switch to the
  // best-effort (minimum hard-violation) control -- NOT a softened relaxation.
  // Relaxing the hard rows would delete their influence; best-effort still spends
  // the actuators driving the state toward the safe set as hard as possible.
  QPResult res = solveWith(/*best_effort=*/false);
  if (!res.success && any_hard) {
    recovered_ = true;
    res = solveWith(/*best_effort=*/true);
  }

  // Last resort: if even the best-effort QP fails to solve (solver error -- it is
  // structurally feasible, so this is numerical), apply the box-clamped nominal.
  LonCtrlVec u_out = U_nom;
  if (res.success) {
    u_out[LDE] = clampd(res.z[LDE], cfg_.de_min, cfg_.de_max);
    u_out[LTDDOT] = clampd(res.z[LTDDOT], cfg_.Tddot_min, cfg_.Tddot_max);
  } else {
    u_out[LDE] = clampd(U_nom[LDE], cfg_.de_min, cfg_.de_max);
    u_out[LTDDOT] = clampd(U_nom[LTDDOT], cfg_.Tddot_min, cfg_.Tddot_max);
  }
  return u_out;
}

}  // namespace autoland
