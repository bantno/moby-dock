#include "autoland/lon_cbf_filter.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace autoland {
namespace {

constexpr double kInf = 1.0e30;  // OSQP +/-inf sentinel
constexpr double kCap = 1.0e20;  // finite cap kept strictly below kInf

double clampd(double v, double lo, double hi) { return std::max(lo, std::min(hi, v)); }

// One assembled barrier row in QP <= form: a . U <= rhs (with optional slack).
struct QRow {
  Eigen::Matrix<double, 1, NUA> a;
  double rhs{0};
  bool hard{true};
  bool finite{true};
};

}  // namespace

LonCtrlVec LonCBFFilter::filter(const LonCtrlVec& U_nom, const LonStateVec& X,
                                const AeroTable& table, const Mixing& mixing,
                                const AircraftConfig& cfg) const {
  recovered_ = false;
  if (!cfg_.enabled) return U_nom;

  const double V = X[LV];
  const double alpha = X[LTH] - X[LGAM];
  const AeroLocal aero = makeAeroLocal(table, mixing, cfg, V, alpha);

  // ---- Assemble the barrier rows -------------------------------------------
  std::vector<QRow> rows;
  auto pushHocbf = [&](const HocbfRow& h, bool hard) {
    QRow r;
    r.a = h.a;
    r.rhs = h.rhs;
    r.hard = hard;
    r.finite = std::isfinite(h.rhs) && std::isfinite(h.a(0, LDE)) &&
               std::isfinite(h.a(0, LTDDOT));
    rows.push_back(r);
  };

  if (cfg_.descent) {
    const DescentBarrier b = makeDescentBarrier(aero, cfg_.v_safe, cfg_.CLmax, V);
    // a_brk(V,gamma) must stay positive for the soft-landing envelope to be
    // well-posed; it is, as long as V stays above ~stall (the airspeed barrier).
    // Warn once if it ever goes non-positive (the descent row then turns
    // non-finite and is dropped below -- a silent loss of the guarantee).
    const double abrk = (0.5 * b.rho * b.Sref / b.mass) * V * V *
                            (b.CLmax * std::cos(X[LGAM]) -
                             b.CDmaxlift * std::sin(X[LGAM])) -
                        b.g;
    if (abrk <= 0.0) {
      static bool warned = false;
      if (!warned) {
        std::cerr << "[lon_cbf] warning: a_brk(V,gamma) <= 0 at V=" << V
                  << " m/s -- descent barrier may be ill-posed\n";
        warned = true;
      }
    }
    auto lie = barrierLie<3>(aero, b, X);
    std::vector<double> Lf(lie.Lf.begin(), lie.Lf.end());
    const std::vector<double> c(cfg_.c_descent.begin(), cfg_.c_descent.end());
    pushHocbf(hocbfRow(Lf, lie.LgLf, c), cfg_.descent_hard);
  }
  if (cfg_.airspeed) {
    AirspeedBarrier b{cfg_.Vmin};
    auto lie = barrierLie<3>(aero, b, X);
    std::vector<double> Lf(lie.Lf.begin(), lie.Lf.end());
    const std::vector<double> c(cfg_.c_airspeed.begin(), cfg_.c_airspeed.end());
    pushHocbf(hocbfRow(Lf, lie.LgLf, c), cfg_.airspeed_hard);
  }
  if (cfg_.thrust_limits) {
    pushHocbf(thrustMinRow(X, cfg_.c_thrust_min[0], cfg_.c_thrust_min[1]), true);
    pushHocbf(thrustMaxRow(X, cfg_.Tmax, cfg_.c_thrust_max[0], cfg_.c_thrust_max[1]), true);
  }

  const int nb = static_cast<int>(rows.size());
  if (nb == 0) return U_nom;

  // ---- Assemble + solve the QP, optionally forcing every row soft ----------
  //   z = [de, Tddot, slack(one per SOFT row)]
  auto solveWith = [&](bool force_soft) -> QPResult {
    std::vector<int> slack_col(nb, -1);
    int nslack = 0;
    for (int i = 0; i < nb; ++i)
      if (force_soft || !rows[i].hard) slack_col[i] = NUA + nslack++;

    const int n = NUA + nslack;
    const int m = nb + NUA + nslack;  // barrier rows + control box + slack >= 0

    Mat P = Mat::Zero(n, n);
    Vec q = Vec::Zero(n);
    P(LDE, LDE) = cfg_.w_de;
    P(LTDDOT, LTDDOT) = cfg_.w_Tddot;
    q[LDE] = -cfg_.w_de * U_nom[LDE];
    q[LTDDOT] = -cfg_.w_Tddot * U_nom[LTDDOT];
    for (int i = 0; i < nb; ++i)
      if (slack_col[i] >= 0) P(slack_col[i], slack_col[i]) = cfg_.slack_penalty;

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

  QPResult res = solveWith(/*force_soft=*/false);
  if (!res.success && any_hard) {
    recovered_ = true;
    res = solveWith(/*force_soft=*/true);
  }

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
