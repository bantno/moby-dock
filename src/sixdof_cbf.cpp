#include "autoland/sixdof_cbf.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <type_traits>
#include "autoland/hocbf.hpp"

namespace autoland {
namespace {
constexpr double kInf = 1e30;
constexpr double kCap = 1e12;

std::array<double, NX> toArr(const StateVec& x) {
  std::array<double, NX> a;
  for (int i = 0; i < NX; ++i) a[i] = x[i];
  return a;
}

// The four barrier families as scalar-generic functors (air-relative
// quantities through the model's wind rotation).
struct AlphaCeil {
  const SixDofCbfModel* m; double amax;
  template <class T> T operator()(const std::array<T, NX>& X) const {
    using std::atan;
    const auto ar = m->airRel<T>(X);
    return T(amax) - atan(ar[2] / ar[0]);
  }
};
struct AirspeedRow {
  const SixDofCbfModel* m; double lim; double sign;  // +1: V - lim; -1: lim - V
  template <class T> T operator()(const std::array<T, NX>& X) const {
    using std::sqrt;
    const auto ar = m->airRel<T>(X);
    const T Vt = sqrt(ar[0] * ar[0] + ar[1] * ar[1] + ar[2] * ar[2] + T(1e-12));
    return sign > 0 ? Vt - lim : T(lim) - Vt;
  }
};
struct BankRow {
  double lim; double sign;  // +1: lim - phi; -1: lim + phi
  template <class T> T operator()(const std::array<T, NX>& X) const {
    return sign > 0 ? T(lim) - X[PHI] : T(lim) + X[PHI];
  }
};
struct BetaRow {
  const SixDofCbfModel* m; double lim; double sign;
  template <class T> T operator()(const std::array<T, NX>& X) const {
    using std::asin;
    using std::sqrt;
    const auto ar = m->airRel<T>(X);
    const T Vt = sqrt(ar[0] * ar[0] + ar[1] * ar[1] + ar[2] * ar[2] + T(1e-12));
    const T beta = asin(ar[1] / Vt);
    return sign > 0 ? T(lim) - beta : T(lim) + beta;
  }
};
struct TerrainRow {  // pure clearance, degree 3
  const TerrainField* tf; double margin;
  template <class T> T operator()(const std::array<T, NX>& X) const {
    return X[H] - tf->height<T>(X[XN], X[Y]) - margin;
  }
};
struct TerrainBrakingRow {  // closure-rate braking form, degree 2
  const SixDofCbfModel* m; const TerrainField* tf;
  double margin, a_brk, v_safe, z_eps, k_neg, lookahead, smax_eps, climb_grad;
  double chi0;         // ground course FROZEN at the current step [rad]
  double fan_half;     // half-angle of the direction fan [rad]
  int fan_n;           // number of fan directions (odd; 1 = single direction)
  // Lookahead terrain height and its gradient: smooth max over the field at
  // the aircraft and at lookahead/2, lookahead ahead along the frozen course.
  // With the direction frozen the height depends on (N, E) only, so the row
  // keeps relative degree 2; the smooth-max partials chain the gradients.
  template <class T>
  void terrainL(const T& N, const T& E, T& h, T& gN, T& gE) const {
    using std::sqrt;
    h = tf->height<T>(N, E);
    auto g = tf->gradientT<T>(N, E);
    gN = g[0]; gE = g[1];
    if (lookahead <= 0.0) return;
    // Height REQUIRED NOW to clear the terrain at s ahead on a sustainable
    // climb gradient: h_T(s) - s * climb_grad, over a FAN of directions
    // about the frozen course (so a small course change moves the smooth
    // maximum continuously instead of sweeping a single ray across a peak).
    // The max over s and direction is the climb-capability lookahead (a
    // slope cannot be "arrested" by a pull-up; it has to be climbed from far
    // enough back, or turned away from -- the fan's lateral spread is what
    // gives the closure rate a smooth lateral gradient to steer with).
    const int nf = std::max(1, fan_n);
    for (int i = 0; i < nf; ++i) {
      const double chi = nf == 1 ? chi0 : chi0 + fan_half * (2.0 * i / (nf - 1) - 1.0);
      const double dN = std::cos(chi), dE = std::sin(chi);
      for (double frac : {0.25, 0.5, 0.75, 1.0}) {
        const double s_ahead = lookahead * frac;
        const T Na = N + (s_ahead * dN), Ea = E + (s_ahead * dE);
        const T ha = tf->height<T>(Na, Ea) - (s_ahead * climb_grad);
        const auto ga = tf->gradientT<T>(Na, Ea);
        const T d = h - ha;
        const T r = sqrt(d * d + T(smax_eps * smax_eps));
        const T hn = 0.5 * (h + ha + r);
        const T wa = 0.5 * (1.0 - d / r);   // d smax / d ha
        const T wh = 1.0 - wa;              // d smax / d h
        gN = wh * gN + wa * ga[0];
        gE = wh * gE + wa * ga[1];
        h = hn;
      }
    }
  }
  template <class T> T operator()(const std::array<T, NX>& X) const {
    using std::cos;
    using std::sin;
    using std::sqrt;
    const T ct = cos(X[THETA]), st = sin(X[THETA]);
    const T cp = cos(X[PHI]), sp = sin(X[PHI]);
    const T cy = cos(X[PSI]), sy = sin(X[PSI]);
    const T Ndot = X[U] * ct * cy + X[V] * (sp * st * cy - cp * sy) + X[W] * (cp * st * cy + sp * sy);
    const T Edot = X[U] * ct * sy + X[V] * (sp * st * sy + cp * cy) + X[W] * (cp * st * sy - sp * cy);
    const T Hdot = X[U] * st - X[V] * sp * ct - X[W] * cp * ct;
    T hT, gN, gE;
    terrainL<T>(X[XN], X[Y], hT, gN, gE);
    const T z = X[H] - hT - margin;
    const T zdot = Hdot - gN * Ndot - gE * Edot;
    const T zp = 0.5 * (z + sqrt(z * z + T(z_eps * z_eps)));   // smooth max(z, 0)
    const T zn = z - zp;                                        // smooth min(z, 0)
    return zdot + sqrt(T(v_safe * v_safe) + 2.0 * a_brk * zp) + k_neg * zn;
  }
};
}  // namespace

std::array<std::array<double, NX>, NUS> SixDofCbfModel::gColumns(
    const std::array<double, NX>& x0) const {
  const auto gm = momentColumns<double>(x0);
  std::array<std::array<double, NX>, NUS> g{};
  for (int k = 0; k < NUS; ++k) {
    g[k].fill(0.0);
    g[k][P] = gm[k][0]; g[k][Q] = gm[k][1]; g[k][R] = gm[k][2];
  }
  return g;
}

// Drift perturbed by epsilon * F_k (epsilon an autodiff dual seeded with
// grad = 1): the first-order-in-epsilon coefficients of its Lie stack are
// L_F h (order 1) and L_F L_f h + L_f L_F h (order 2).
struct PerturbedDrift {
  const SixDofCbfModel& f;
  int k;
  autodiff::dual eps;
  template <class T>
  std::array<T, NX> operator()(const std::array<T, NX>& X) const {
    std::array<T, NX> xd = f(X);
    const auto F = f.forceColumnsT<T>(X);
    for (int i = 0; i < NX; ++i) xd[i] = xd[i] + T(eps) * F[k][i];
    return xd;
  }
};

SixDofCbfFilter::SixDofCbfFilter(const BeaverDynamics& dyn,
                                 const SixDofCbfConfig& cfg)
    : dyn_(dyn), cfg_(cfg), solver_(std::make_unique<OsqpSolver>()) {}
SixDofCbfFilter::~SixDofCbfFilter() = default;

std::vector<SixDofCbfRow> SixDofCbfFilter::rows(const StateVec& x, double dT,
                                                const SurfVec& u_prev,
                                                const Eigen::Vector3d& W_earth) const {
  const SixDofCbfModel f(dyn_, u_prev, dT, W_earth);
  const std::array<double, NX> x0 = toArr(x);
  const auto g = f.gColumns(x0);
  const SurfaceLimits& L = cfg_.limits;
  const double du_step = L.rate * cfg_.dt;
  const double range[NUS] = {L.de_max - L.de_min, L.da_max - L.da_min,
                             L.dr_max - L.dr_min};

  std::vector<SixDofCbfRow> out;
  // Degree-R row: Lie stack, control row, and the force-term margin
  // sum_k sum_{j=1..R} coeff_j |D_j,k| dU_k (see the header).
  auto pushR = [&](auto Rtag, const std::string& name, const auto& b,
                   const std::vector<double>& c, bool hard) {
    constexpr int R = decltype(Rtag)::value;
    SixDofCbfRow r;
    r.name = name;
    r.hard = hard;
    r.degree = R;
    const std::array<double, R + 1> Lf = lieDrift<R, NX>(f, b, x0);
    r.h = Lf[0]; r.Lfh = Lf[1]; r.Lf2h = Lf[2];
    if constexpr (R >= 3) r.Lf3h = Lf[3];
    const std::vector<double> coeff = hocbfPsiCoeffs(c);
    double drift = 0.0;
    for (int j = 0; j <= R; ++j) drift += coeff[j] * Lf[j];
    for (int k = 0; k < NUS; ++k) r.a[k] = -lieAlong<R, NX>(f, b, x0, g[k]);
    double margin = cfg_.margin_extra;
    if (cfg_.margin_force_terms) {
      using autodiff::dual;
      std::array<dual, NX> xs;
      for (int i = 0; i < NX; ++i) xs[i] = x0[i];
      for (int k = 0; k < NUS; ++k) {
        PerturbedDrift fe{f, k, dual(0.0)};
        fe.eps.grad = 1.0;
        const std::array<dual, R + 1> Le = lieJet<R, NX, dual>(fe, b, xs);
        const double du = std::min(du_step, range[k]);
        r.D1[k] = autodiff::detail::derivative(Le[1]);
        r.D2[k] = autodiff::detail::derivative(Le[2]);
        double m = coeff[1] * std::abs(r.D1[k]) + coeff[2] * std::abs(r.D2[k]);
        if constexpr (R >= 3) {
          r.D3[k] = autodiff::detail::derivative(Le[3]);
          m += coeff[3] * std::abs(r.D3[k]);
        }
        margin += m * du;
      }
    }
    r.margin = margin;
    r.rhs = drift - margin;
    r.finite = std::isfinite(r.rhs) && r.a.allFinite() && std::isfinite(r.h);
    out.push_back(r);
  };
  auto push = [&](const std::string& name, const auto& b,
                  const std::array<double, 2>& c, bool hard) {
    pushR(std::integral_constant<int, 2>{}, name, b, {c[0], c[1]}, hard);
  };

  if (cfg_.row_alpha) push("alpha_max", AlphaCeil{&f, cfg_.alpha_max}, cfg_.c_alpha, cfg_.hard_alpha);
  if (cfg_.row_V) {
    push("V_min", AirspeedRow{&f, cfg_.V_min, +1.0}, cfg_.c_V, cfg_.hard_V);
    push("V_max", AirspeedRow{&f, cfg_.V_max, -1.0}, cfg_.c_V, cfg_.hard_V);
  }
  if (cfg_.row_phi) {
    push("phi_max+", BankRow{cfg_.phi_max, +1.0}, cfg_.c_phi, cfg_.hard_phi);
    push("phi_max-", BankRow{cfg_.phi_max, -1.0}, cfg_.c_phi, cfg_.hard_phi);
  }
  if (cfg_.row_beta) {
    push("beta_max+", BetaRow{&f, cfg_.beta_max, +1.0}, cfg_.c_beta, cfg_.hard_beta);
    push("beta_max-", BetaRow{&f, cfg_.beta_max, -1.0}, cfg_.c_beta, cfg_.hard_beta);
  }
  if (cfg_.row_terrain && cfg_.terrain) {
    const double clr = x[H] - cfg_.terrain->height(x[XN], x[Y]);
    if (clr < cfg_.terrain_gate) {
      if (cfg_.terrain_degree3)
        pushR(std::integral_constant<int, 3>{}, "terrain",
              TerrainRow{cfg_.terrain, cfg_.terrain_margin},
              {cfg_.c_terrain[0], cfg_.c_terrain[1], cfg_.c_terrain[2]}, cfg_.hard_terrain);
      else {
        // Ground course at x0 (frozen for the step; see TerrainBrakingRow).
        const double ct = std::cos(x[THETA]), st = std::sin(x[THETA]);
        const double cp = std::cos(x[PHI]), sp = std::sin(x[PHI]);
        const double cy = std::cos(x[PSI]), sy = std::sin(x[PSI]);
        const double Ndot = x[U] * ct * cy + x[V] * (sp * st * cy - cp * sy) + x[W] * (cp * st * cy + sp * sy);
        const double Edot = x[U] * ct * sy + x[V] * (sp * st * sy + cp * cy) + x[W] * (cp * st * sy - sp * cy);
        (void)Ndot; (void)Edot;
        pushR(std::integral_constant<int, 2>{}, "terrain",
              TerrainBrakingRow{&f, cfg_.terrain, cfg_.terrain_margin, cfg_.terrain_a_brk,
                                cfg_.terrain_v_safe, cfg_.terrain_z_eps, cfg_.terrain_k_neg,
                                cfg_.terrain_lookahead, cfg_.terrain_smax_eps,
                                cfg_.terrain_climb_grad, std::atan2(Edot, Ndot),
                                cfg_.terrain_fan_half, cfg_.terrain_fan_n},
              {cfg_.c_terrain[0], cfg_.c_terrain[1]}, cfg_.hard_terrain);
      }
    }
  }
  return out;
}

namespace {
struct Box {
  double lo[NUS], hi[NUS];
};
Box inputBox(const SurfaceLimits& L, double dt, const SurfVec& u_prev, bool rate_box) {
  const double du = rate_box ? L.rate * dt : 1e9;
  Box b;
  const double mn[NUS] = {L.de_min, L.da_min, L.dr_min};
  const double mx[NUS] = {L.de_max, L.da_max, L.dr_max};
  for (int k = 0; k < NUS; ++k) {
    b.lo[k] = std::max(mn[k], u_prev[k] - du);
    b.hi[k] = std::min(mx[k], u_prev[k] + du);
    if (b.hi[k] < b.lo[k]) b.lo[k] = b.hi[k] = std::clamp(u_prev[k], mn[k], mx[k]);
  }
  return b;
}
double wRow(const SixDofCbfConfig& c, const std::string& n) {
  if (n.rfind("alpha", 0) == 0) return c.w_alpha;
  if (n.rfind("V_", 0) == 0) return c.w_V;
  if (n.rfind("phi", 0) == 0) return c.w_phi;
  if (n.rfind("terrain", 0) == 0) return c.w_terrain;
  return c.w_beta;
}
}  // namespace

bool SixDofCbfFilter::hardSetFeasible(const StateVec& x, double dT,
                                      const SurfVec& u_prev,
                                      const Eigen::Vector3d& W_earth) const {
  const std::vector<SixDofCbfRow> rs = rows(x, dT, u_prev, W_earth);
  const Box box = inputBox(cfg_.limits, cfg_.dt, u_prev, cfg_.rate_box);
  int nh = 0;
  for (const auto& r : rs) if (r.hard && r.finite) ++nh;
  if (nh == 0) return true;
  const int m = nh + NUS;
  Mat P = Mat::Identity(NUS, NUS) * 1e-6;
  Vec q = Vec::Zero(NUS);
  Mat A = Mat::Zero(m, NUS);
  Vec lo = Vec::Zero(m), hi = Vec::Zero(m);
  int row = 0;
  for (const auto& r : rs) {
    if (!(r.hard && r.finite)) continue;
    const double s = std::max(1.0, r.a.cwiseAbs().maxCoeff());
    for (int k = 0; k < NUS; ++k) A(row, k) = r.a[k] / s;
    lo[row] = -kInf; hi[row] = std::clamp(r.rhs / s, -kCap, kCap);
    ++row;
  }
  for (int k = 0; k < NUS; ++k) { A(row, k) = 1.0; lo[row] = box.lo[k]; hi[row] = box.hi[k]; ++row; }
  return solver_->solve(P, q, A, lo, hi).success;
}

SurfVec SixDofCbfFilter::filter(const SurfVec& U_nom, const StateVec& x,
                                double dT, const SurfVec& u_prev,
                                SixDofCbfDiag* diag,
                                const Eigen::Vector3d& W_earth) const {
  const auto t0 = std::chrono::steady_clock::now();
  const std::vector<SixDofCbfRow> rs = rows(x, dT, u_prev, W_earth);
  const Box box = inputBox(cfg_.limits, cfg_.dt, u_prev, cfg_.rate_box);
  const SurfaceLimits& L = cfg_.limits;
  const double range[NUS] = {L.de_max - L.de_min, L.da_max - L.da_min,
                             L.dr_max - L.dr_min};
  const double w[NUS] = {cfg_.w_de, cfg_.w_da, cfg_.w_dr};

  SixDofCbfDiag d;
  d.rows = rs;
  auto fam = [&](const std::string& prefix) {
    double m = 1e30;
    for (const auto& r : rs) if (r.name.rfind(prefix, 0) == 0) m = std::min(m, r.h);
    return m < 1e30 ? m : 0.0;
  };
  d.h_alpha = fam("alpha"); d.h_V = fam("V_"); d.h_phi = fam("phi"); d.h_beta = fam("beta");
  d.h_terrain = fam("terrain");
  for (const auto& r : rs) { d.min_h = std::min(d.min_h, r.h); d.margin_max = std::max(d.margin_max, r.margin); }

  const int nb = static_cast<int>(rs.size());
  bool any_hard = false;
  for (const auto& r : rs) any_hard = any_hard || (r.hard && r.finite);

  auto solveWith = [&](bool best_effort) -> QPResult {
    std::vector<int> slack_col(nb, -1);
    int nslack = 0;
    for (int i = 0; i < nb; ++i)
      if (rs[i].finite && (best_effort || !rs[i].hard)) slack_col[i] = NUS + nslack++;
    const int n = NUS + nslack;
    const int m = nb + NUS + nslack;
    Mat P = Mat::Zero(n, n);
    Vec q = Vec::Zero(n);
    for (int k = 0; k < NUS; ++k) {
      const double wk = w[k] / (range[k] * range[k]);
      P(k, k) = wk;
      q[k] = -wk * U_nom[k];
    }
    for (int i = 0; i < nb; ++i)
      if (slack_col[i] >= 0)
        P(slack_col[i], slack_col[i]) =
            (best_effort && rs[i].hard) ? cfg_.best_effort_penalty : wRow(cfg_, rs[i].name);
    Mat A = Mat::Zero(m, n);
    Vec lo = Vec::Zero(m), hi = Vec::Zero(m);
    int row = 0;
    for (int i = 0; i < nb; ++i) {
      const SixDofCbfRow& r = rs[i];
      const double s = r.finite ? std::max(1.0, r.a.cwiseAbs().maxCoeff()) : 1.0;
      for (int k = 0; k < NUS; ++k) A(row, k) = r.finite ? r.a[k] / s : 0.0;
      if (slack_col[i] >= 0) A(row, slack_col[i]) = -1.0;
      lo[row] = -kInf;
      hi[row] = r.finite ? std::clamp(r.rhs / s, -kCap, kCap) : kInf;
      ++row;
    }
    for (int k = 0; k < NUS; ++k) { A(row, k) = 1.0; lo[row] = box.lo[k]; hi[row] = box.hi[k]; ++row; }
    for (int i = 0; i < nb; ++i)
      if (slack_col[i] >= 0) { A(row, slack_col[i]) = 1.0; lo[row] = 0.0; hi[row] = kInf; ++row; }
    return solver_->solve(P, q, A, lo, hi);
  };

  QPResult res = solveWith(false);
  if (!res.success && any_hard) {
    d.best_effort = true;
    res = solveWith(true);
  }
  d.qp_status = res.status_val;
  SurfVec out;
  if (res.success) {
    for (int k = 0; k < NUS; ++k) out[k] = std::clamp(res.z[k], box.lo[k], box.hi[k]);
  } else {
    d.solved = false;
    for (int k = 0; k < NUS; ++k) out[k] = std::clamp(U_nom[k], box.lo[k], box.hi[k]);
  }
  d.active = (out - U_nom).cwiseAbs().maxCoeff() > 1e-6;
  d.t_us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count();
  if (diag) *diag = d;
  return out;
}

}  // namespace autoland
