#include "autoland/cbf.hpp"

#include <algorithm>
#include <cmath>

namespace autoland {
namespace {

constexpr double kInf = 1.0e30;   // OSQP +/-inf sentinel
constexpr double kCap = 1.0e20;   // finite cap kept strictly below kInf, so a
                                  // pathological barrier row never produces
                                  // lo > hi (which OSQP rejects).

double clampd(double v, double lo, double hi) {
  return std::max(lo, std::min(hi, v));
}

}  // namespace

CBFFilter::CBFFilter()
    : cfg_(CBFConfig{}), lim_(SurfaceLimits{}),
      solver_(std::make_unique<OsqpSolver>()) {}

CBFFilter::CBFFilter(CBFConfig cfg, SurfaceLimits limits,
                     std::unique_ptr<QPSolver> solver)
    : cfg_(std::move(cfg)),
      lim_(limits),
      solver_(solver ? std::move(solver)
                     : std::unique_ptr<QPSolver>(std::make_unique<OsqpSolver>())) {}

std::pair<double, double> CBFFilter::ctrlBounds(int ctrl) const {
  switch (ctrl) {
    case DE: return {lim_.de_min, lim_.de_max};
    case DA: return {lim_.da_min, lim_.da_max};
    case DR: return {lim_.dr_min, lim_.dr_max};
    case DT: return {lim_.dT_min, lim_.dT_max};
    default: return {-kInf, kInf};
  }
}

CtrlVec CBFFilter::filter(const CtrlVec& u_nom, const StateVec& x,
                          const ControlAffineModel& model,
                          const std::vector<Barrier>& barriers) const {
  recovered_ = false;
  if (!enabled_ || barriers.empty()) return u_nom;

  const int na = static_cast<int>(cfg_.ctrl_idx.size());  // active controls
  const int nb = static_cast<int>(barriers.size());

  // ---- Control-affine pieces at the current state ---------------------------
  const StateVec f = model.drift(x);  // f(x)        (NX)
  const Mat g = model.ctrlMatrix(x);  // g(x) = B    (NX x NU)

  std::vector<bool> active(NU, false);
  for (int idx : cfg_.ctrl_idx) active[idx] = true;

  // ---- Per-barrier CBF row, precomputed once --------------------------------
  //   Lgh_active . z (+ slack) >= rhs,   rhs = -(Lf h + alpha(h) + Lgh_fixed.u_fixed)
  struct Row {
    Eigen::RowVectorXd Lgh;  // 1 x NU
    double rhs{0.0};
    bool finite{true};
  };
  std::vector<Row> rows(nb);
  for (int i = 0; i < nb; ++i) {
    const Barrier& b = barriers[i];
    const Eigen::RowVectorXd grad = b.grad(x);  // 1 x NX
    Row& r = rows[i];
    r.Lgh = grad * g;                           // 1 x NU
    const double Lfh = (grad * f)(0, 0);        // grad . f
    const double alpha = b.alpha(b.h(x));
    double fixed_term = 0.0;
    for (int j = 0; j < NU; ++j)
      if (!active[j]) fixed_term += r.Lgh(j) * u_nom[j];
    r.rhs = -(Lfh + alpha + fixed_term);
    r.finite = std::isfinite(r.rhs);
    for (int k = 0; k < na; ++k)
      r.finite = r.finite && std::isfinite(r.Lgh(cfg_.ctrl_idx[k]));
  }

  // ---- Assemble + solve the QP, optionally forcing every barrier soft -------
  //   z = [u_active (na) ; slack (one per SOFT barrier)]
  auto solveWith = [&](bool force_soft) -> QPResult {
    std::vector<int> slack_col(nb, -1);  // QP column of each barrier's slack
    int nslack = 0;
    for (int i = 0; i < nb; ++i)
      if (force_soft || !barriers[i].hard) slack_col[i] = na + nslack++;

    const int n = na + nslack;
    const int m = nb + na + nslack;  // barrier rows + control box + slack >= 0

    Mat P = Mat::Zero(n, n);
    Vec q = Vec::Zero(n);
    for (int k = 0; k < na; ++k) {
      const double w = (k < static_cast<int>(cfg_.ctrl_weights.size()))
                           ? cfg_.ctrl_weights[k]
                           : 1.0;
      P(k, k) = w;
      q[k] = -w * u_nom[cfg_.ctrl_idx[k]];
    }
    for (int i = 0; i < nb; ++i)
      if (slack_col[i] >= 0) P(slack_col[i], slack_col[i]) = cfg_.slack_penalty;

    Mat A = Mat::Zero(m, n);
    Vec lo = Vec::Zero(m);
    Vec hi = Vec::Zero(m);
    int row = 0;

    // (1) One CBF inequality per barrier. Hard barriers get no slack column, so
    //     their row must hold exactly; soft barriers can absorb violation.
    for (int i = 0; i < nb; ++i) {
      const Row& r = rows[i];
      for (int k = 0; k < na; ++k)
        A(row, k) = r.finite ? r.Lgh(cfg_.ctrl_idx[k]) : 0.0;
      if (slack_col[i] >= 0) A(row, slack_col[i]) = 1.0;
      // A non-finite row (bad linearization) is dropped to inactive so it can
      // never make even a HARD constraint spuriously infeasible.
      const double lo_val = r.finite ? r.rhs : -kInf;
      lo[row] = clampd(lo_val, -kInf, kCap);
      hi[row] = kInf;
      ++row;
    }
    // (2) Box constraints on the active controls.
    for (int k = 0; k < na; ++k) {
      const auto [clo, chi] = ctrlBounds(cfg_.ctrl_idx[k]);
      A(row, k) = 1.0;
      lo[row] = clo;
      hi[row] = chi;
      ++row;
    }
    // (3) Non-negative slacks.
    for (int i = 0; i < nb; ++i)
      if (slack_col[i] >= 0) {
        A(row, slack_col[i]) = 1.0;
        lo[row] = 0.0;
        hi[row] = kInf;
        ++row;
      }

    return solver_->solve(P, q, A, lo, hi);
  };

  bool any_hard = false;
  for (const Barrier& b : barriers) any_hard = any_hard || b.hard;

  // Honor hardness first; if the hard set is infeasible, recover by softening
  // every barrier so the filter still returns a sane, in-box control.
  QPResult res = solveWith(/*force_soft=*/false);
  if (!res.success && any_hard) {
    recovered_ = true;
    res = solveWith(/*force_soft=*/true);
  }

  CtrlVec u_out = u_nom;
  for (int k = 0; k < na; ++k) {
    const int idx = cfg_.ctrl_idx[k];
    const auto [clo, chi] = ctrlBounds(idx);
    // On total failure, clamp the nominal to the box; otherwise use the solution.
    u_out[idx] = res.success ? clampd(res.z[k], clo, chi)
                             : clampd(u_nom[idx], clo, chi);
  }
  return u_out;
}

}  // namespace autoland
