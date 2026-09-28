#include "autoland/energy_reach_cbf.hpp"

#include <algorithm>
#include <cmath>

namespace autoland {
namespace {
constexpr double kInf = 1.0e30;  // QPSolver one-sided sentinel
}

Eigen::Vector2d EnergyReachFilter::filter(double ua_nom, double be_nom,
                                          const ReachStateArr& X,
                                          EnergyReachDiag* diag) const {
  const ReachParams& p = cfg_.prm;
  const ReachLie lie = reachLie(p, X);

  // AoA band folded into the u_alpha variable bounds (each is itself a valid
  // RD-1 CBF row: alphadot <= c_alpha*(alpha_max - alpha) etc.). If alpha has
  // escaped the band beyond the rate the rows can ask for, command the full
  // rate back toward it.
  const double alpha = X[ERAL];
  double ua_lo = std::max(-cfg_.ua_max, -cfg_.c_alpha * (alpha - cfg_.alpha_cmd_min));
  double ua_hi = std::min(cfg_.ua_max, cfg_.c_alpha * (cfg_.alpha_cmd_max - alpha));
  if (ua_lo > ua_hi) {
    const double mid = 0.5 * (cfg_.alpha_cmd_min + cfg_.alpha_cmd_max);
    ua_lo = ua_hi = (alpha > mid) ? -cfg_.ua_max : cfg_.ua_max;
  }

  // Diagnostics (recomputed with the same guarded formulas as the barrier).
  EnergyReachDiag d;
  d.h = lie.h;
  d.Lfh = lie.Lfh;
  d.A_ua = lie.A_ua;
  d.A_be = lie.A_be;
  {
    const ReachEnergyBarrier b(p);
    d.d = b.arcDistance(X[ERN], X[ERE], X[ERCHI]);
    const double ats = b.alphaFloored(alpha);
    d.LD = p.CLa * ats / (p.CD0 + p.CDa2 * ats * ats);
    const double dN = p.Nt - X[ERN], dE = p.Et - X[ERE];
    d.phi = (dN * dN + dE * dE < p.d_near * p.d_near)
                ? 0.0
                : std::atan2(dN * std::sin(X[ERCHI]) - dE * std::cos(X[ERCHI]),
                             dN * std::cos(X[ERCHI]) + dE * std::sin(X[ERCHI]));
  }

  const double ua_pass = std::clamp(ua_nom, ua_lo, ua_hi);
  const double be_pass = std::clamp(be_nom, -cfg_.beta_max, cfg_.beta_max);
  const double rhs = lie.Lfh + cfg_.ck * lie.h;  // hdot >= -ck h  <=>  A.u >= -rhs

  const bool finite = std::isfinite(lie.h) && std::isfinite(rhs) &&
                      std::isfinite(lie.A_ua) && std::isfinite(lie.A_be);
  if (!cfg_.enabled || !finite) {
    d.finite = finite;
    d.ua = ua_pass;
    d.be = be_pass;
    if (diag) *diag = d;
    return {ua_pass, be_pass};
  }

  // Feasibility against the input bounds, checked in closed form: the max of
  // A.u over the box is attained railing each control toward its A-sign.
  const double ua_best = (lie.A_ua >= 0.0) ? ua_hi : ua_lo;
  const double be_best = (lie.A_be >= 0.0) ? cfg_.beta_max : -cfg_.beta_max;
  const double Au_max = lie.A_ua * ua_best + lie.A_be * be_best;
  if (Au_max < -rhs) {
    // The hard row cannot be met: minimum-violation (max-hdot) action.
    d.best_effort = true;
    d.active = true;
    d.ua = ua_best;
    d.be = be_best;
    d.viol = -rhs - Au_max;
    if (diag) *diag = d;
    return {d.ua, d.be};
  }

  // Range-scaled tracking weights: full scale of each channel priced equally.
  const double w_ua = 1.0 / (cfg_.ua_max * cfg_.ua_max);
  const double w_be = 1.0 / (cfg_.beta_max * cfg_.beta_max);

  Mat P = Mat::Zero(2, 2);
  P(0, 0) = w_ua;
  P(1, 1) = w_be;
  Vec q(2);
  q << -w_ua * ua_nom, -w_be * be_nom;

  // Energy-reach row, NORMALIZED by the largest control coefficient: the raw
  // A_ua runs O(1e4) (g*d*dPsi) against the O(1) variable-bound rows, and
  // that imbalance makes OSQP return spurious infeasibility certificates on
  // feasible steps.
  const double s = std::max({1.0, std::abs(lie.A_ua), std::abs(lie.A_be)});
  Mat A = Mat::Zero(3, 2);
  Vec l(3), u(3);
  A(0, 0) = -lie.A_ua / s;
  A(0, 1) = -lie.A_be / s;
  l(0) = -kInf;
  u(0) = rhs / s;
  A(1, 0) = 1.0; l(1) = ua_lo; u(1) = ua_hi;
  A(2, 1) = 1.0; l(2) = -cfg_.beta_max; u(2) = cfg_.beta_max;

  const QPResult r = solver_->solve(P, q, A, l, u);
  if (!r.success) {
    // Should not happen (feasibility was certified above); fail toward the
    // max-hdot action rather than the nominal, and record it.
    d.best_effort = true;
    d.active = true;
    d.ua = ua_best;
    d.be = be_best;
    if (diag) *diag = d;
    return {d.ua, d.be};
  }

  d.ua = std::clamp(r.z(0), ua_lo, ua_hi);  // shave OSQP's solution tolerance
  d.be = std::clamp(r.z(1), -cfg_.beta_max, cfg_.beta_max);
  d.viol = std::max(0.0, -rhs - (lie.A_ua * d.ua + lie.A_be * d.be));
  d.active = std::abs(d.ua - ua_pass) > 1e-9 || std::abs(d.be - be_pass) > 1e-9;
  if (diag) *diag = d;
  return {d.ua, d.be};
}

ReachPolarFit fitBeaverReachPolar(const BeaverDynamics& dyn, double V_ref,
                                  double dT, double alpha_min, double alpha_max,
                                  int n) {
  const double dpt = dyn.dpt(dT, V_ref);
  const double flap = dyn.config().flap;
  const BeaverAeroCoef& k = dyn.config().aero;

  Eigen::VectorXd al(n), CL(n), CD(n);
  for (int i = 0; i < n; ++i) {
    const double a = alpha_min + (alpha_max - alpha_min) * i / (n - 1);
    const std::array<double, 6> C =
        beaverAeroCoeffs(a, 0.0, 0.0, 0.0, 0.0, V_ref, 0.0, 0.0, 0.0, flap, dpt, k);
    al(i) = a;
    CL(i) = C[0] * std::sin(a) - C[2] * std::cos(a);
    CD(i) = -(C[0] * std::cos(a) + C[2] * std::sin(a));
  }

  ReachPolarFit f;
  // CL ~ c0 + CLa*alpha  ->  alpha0L = -c0/CLa.
  {
    Eigen::MatrixXd B(n, 2);
    B.col(0).setOnes();
    B.col(1) = al;
    const Eigen::Vector2d c = (B.transpose() * B).ldlt().solve(B.transpose() * CL);
    f.CLa = c(1);
    f.alpha0L = -c(0) / c(1);
    f.rms_CL = std::sqrt((CL - B * c).squaredNorm() / n);
  }
  // CD ~ CD0 + CDa2*(alpha - alpha0L)^2.
  {
    Eigen::MatrixXd B(n, 2);
    B.col(0).setOnes();
    B.col(1) = (al.array() - f.alpha0L).square();
    const Eigen::Vector2d c = (B.transpose() * B).ldlt().solve(B.transpose() * CD);
    f.CD0 = c(0);
    f.CDa2 = c(1);
    f.rms_CD = std::sqrt((CD - B * c).squaredNorm() / n);
  }
  f.alpha_star = std::sqrt(std::max(0.0, f.CD0 / f.CDa2));
  f.LD_max = f.CLa / (2.0 * std::sqrt(std::max(1e-12, f.CD0 * f.CDa2)));
  return f;
}

}  // namespace autoland
