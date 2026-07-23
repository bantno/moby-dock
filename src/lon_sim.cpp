#include "autoland/lon_sim.hpp"

#include <yaml-cpp/yaml.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

#include "autoland/dynamics.hpp"
#include "autoland/hocbf.hpp"
#include "autoland/trim.hpp"

namespace autoland {
namespace {
constexpr double kDeg = M_PI / 180.0;

double getOr(const YAML::Node& n, const std::string& key, double def) {
  return (n && n[key]) ? n[key].as<double>() : def;
}
bool getOrB(const YAML::Node& n, const std::string& key, bool def) {
  return (n && n[key]) ? n[key].as<bool>() : def;
}
std::array<double, 3> getArr3(const YAML::Node& n, const std::string& key,
                              std::array<double, 3> def) {
  if (n && n[key] && n[key].IsSequence() && n[key].size() == 3)
    return {n[key][0].as<double>(), n[key][1].as<double>(), n[key][2].as<double>()};
  return def;
}
std::array<double, 2> getArr2(const YAML::Node& n, const std::string& key,
                              std::array<double, 2> def) {
  if (n && n[key] && n[key].IsSequence() && n[key].size() == 2)
    return {n[key][0].as<double>(), n[key][1].as<double>()};
  return def;
}

double thrustOf(const AircraftConfig& cfg, double throttle, double V) {
  return std::max(0.0, throttle * (cfg.thrust.T_static - cfg.thrust.k_v * V));
}
}  // namespace

LonStateVec lonXdotFull(const AeroTable& table, const Mixing& mixing,
                        const AircraftConfig& cfg, const LonStateVec& X,
                        const LonCtrlVec& U) {
  const AeroLocal a = makeAeroLocal(table, mixing, cfg, X[LV], X[LTH] - X[LGAM]);
  return lonXdot(a, X, U);
}

LonStateVec lonXdotFullWind(const AeroTable& table, const Mixing& mixing,
                            const AircraftConfig& cfg, const LonStateVec& X,
                            const LonCtrlVec& U, const GustWind& W,
                            const GustWind& Wdot) {
  LonStateVec xd = lonXdotFull(table, mixing, cfg, X, U);
  const double sg = std::sin(X[LGAM]), cg = std::cos(X[LGAM]);
  xd[LH] += W.w;
  xd[LV] += -(Wdot.u * cg + Wdot.w * sg);
  xd[LGAM] += (Wdot.u * sg - Wdot.w * cg) / X[LV];
  return xd;
}

namespace {
LonStateVec rk4(const AeroTable& t, const Mixing& mx, const AircraftConfig& cfg,
                const LonStateVec& X, const LonCtrlVec& U, double dt) {
  const LonStateVec k1 = lonXdotFull(t, mx, cfg, X, U);
  const LonStateVec k2 = lonXdotFull(t, mx, cfg, X + 0.5 * dt * k1, U);
  const LonStateVec k3 = lonXdotFull(t, mx, cfg, X + 0.5 * dt * k2, U);
  const LonStateVec k4 = lonXdotFull(t, mx, cfg, X + dt * k3, U);
  return X + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

// RK4 step of the wind-augmented system [X; x_gust]. The gust penetration
// distance x_gust is a genuine state (xdot = V after onset, like the Simulink
// block's airspeed integrator), so it is integrated jointly with X rather than
// held over the step. t is the step start time (gates the onset within
// substeps).
struct WindStep {
  LonStateVec X;
  double xg;
};

WindStep rk4Wind(const AeroTable& tab, const Mixing& mx, const AircraftConfig& cfg,
                 const DiscreteGustConfig& g, double t, const LonStateVec& X,
                 double xg, const LonCtrlVec& U, double dt) {
  auto deriv = [&](double tl, const LonStateVec& Xl, double xgl, LonStateVec& xd,
                   double& xgd) {
    xgd = gustXdot(g, tl, Xl[LV]);
    xd = lonXdotFullWind(tab, mx, cfg, Xl, U, gustWind(g, xgl),
                         gustWindRate(g, xgl, xgd));
  };
  LonStateVec k1, k2, k3, k4;
  double x1, x2, x3, x4;
  deriv(t, X, xg, k1, x1);
  deriv(t + 0.5 * dt, X + 0.5 * dt * k1, xg + 0.5 * dt * x1, k2, x2);
  deriv(t + 0.5 * dt, X + 0.5 * dt * k2, xg + 0.5 * dt * x2, k3, x3);
  deriv(t + dt, X + dt * k3, xg + dt * x3, k4, x4);
  return {X + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4),
          xg + (dt / 6.0) * (x1 + 2.0 * x2 + 2.0 * x3 + x4)};
}

// Steady-flight longitudinal trim of the lon model itself: solve for
// [alpha, delta_e, T] (with q = Tdot = Tddot = 0) so that [Vdot, gammadot, qdot]
// vanish at (V, gamma). Mirrors the body-axis trim() Newton scheme but on the
// lon EOM, so the result is a true equilibrium of THIS model (no startup
// transient). Thrust T is a state here (Newtons), so unlike the body-axis trim
// there is no throttle clamp to make the Jacobian singular.
struct LonTrimResult { double alpha, theta, de, T, residual; bool converged; };

LonTrimResult lonTrim(const AeroTable& table, const Mixing& mx,
                      const AircraftConfig& cfg, double V, double gamma) {
  auto residual = [&](const Eigen::Vector3d& v) -> Eigen::Vector3d {
    LonStateVec X;
    X << 0.0, V, gamma, gamma + v[0], 0.0, v[2], 0.0;  // alpha=v0, T=v2, q=Tdot=0
    LonCtrlVec U;
    U << v[1], 0.0;                                     // de=v1, Tddot=0
    const LonStateVec xd = lonXdotFull(table, mx, cfg, X, U);
    return Eigen::Vector3d(xd[LV], xd[LGAM], xd[LQ]);   // Vdot, gammadot, qdot
  };

  Eigen::Vector3d v(0.05, 0.0, 0.3 * cfg.thrust.T_static);  // alpha, de, T guess
  const double hh = 1e-6;
  Eigen::Vector3d R = residual(v);
  int it = 0;
  for (; it < 100; ++it) {
    if (R.norm() < 1e-10) break;
    Eigen::Matrix3d J;
    for (int j = 0; j < 3; ++j) {
      Eigen::Vector3d vp = v, vm = v;
      vp[j] += hh;
      vm[j] -= hh;
      J.col(j) = (residual(vp) - residual(vm)) / (2.0 * hh);
    }
    Eigen::Vector3d dv = J.fullPivLu().solve(-R);
    double step = 1.0;
    Eigen::Vector3d vn = v + step * dv;
    Eigen::Vector3d Rn = residual(vn);
    for (int bt = 0; bt < 20 && Rn.norm() > R.norm(); ++bt) {
      step *= 0.5;
      vn = v + step * dv;
      Rn = residual(vn);
    }
    v = vn;
    R = Rn;
  }
  return {v[0], gamma + v[0], v[1], v[2], R.norm(), R.norm() < 1e-6};
}
}  // namespace

LonSim::LonSim(const std::string& stab_path, const std::string& aircraft_yaml,
               const std::string& scenario_yaml)
    : table_(AeroTable::fromFile(stab_path)),
      ac_(loadAircraftConfig(aircraft_yaml)) {
  mixing_ = std::make_unique<Mixing>(Mixing::build(ac_, table_));

  YAML::Node root = YAML::LoadFile(scenario_yaml);
  const double V_app = getOr(root, "V_app", 18.0);
  const double gamma_app = getOr(root, "gamma_app_deg", -3.0) * kDeg;

  // Trim the full 6-DOF model for nominal feedforward + initial condition.
  Dynamics dyn(table_, *mixing_, ac_);
  TrimResult tr = trim(dyn, V_app, gamma_app);
  const double T_trim = thrustOf(ac_, tr.u[DT], V_app);

  // --- Nominal defaults (then YAML overrides) ---
  LonNominalConfig& nom = sc_.nominal;
  nom.theta_trim = tr.theta;
  nom.gamma_ref = gamma_app;
  nom.T_set = T_trim;
  YAML::Node yn = root["nominal"];
  nom.T_set = getOr(yn, "T_set", nom.T_set);
  nom.Kp_T = getOr(yn, "Kp_T", 4.0);
  nom.Kd_T = getOr(yn, "Kd_T", 4.0);
  nom.Kp_gamma = getOr(yn, "Kp_gamma", 2.0);
  nom.Ki_gamma = getOr(yn, "Ki_gamma", 0.5);
  nom.theta_cmd_max = getOr(yn, "theta_cmd_max_deg", 17.0) * kDeg;
  nom.Kp_theta = getOr(yn, "Kp_theta", 6.0);
  nom.Kq = getOr(yn, "Kq", 1.5);

  // --- CBF defaults (then YAML overrides) ---
  LonCBFConfig& cb = sc_.cbf;
  YAML::Node yc = root["cbf"];
  cb.enabled = getOrB(yc, "enabled", true);
  cb.thrust_limits = getOrB(yc, "thrust_limits", true);
  cb.impact = getOrB(yc, "impact", true);
  cb.stall = getOrB(yc, "stall", true);
  cb.noseup = getOrB(yc, "noseup", true);
  cb.energy = getOrB(yc, "energy", true);
  cb.impact_hard = getOrB(yc, "impact_hard", true);
  cb.stall_hard = getOrB(yc, "stall_hard", false);
  cb.noseup_hard = getOrB(yc, "noseup_hard", false);
  cb.energy_hard = getOrB(yc, "energy_hard", false);
  cb.alpha_stall = getOr(yc, "alpha_stall_deg", 11.0) * kDeg;
  cb.stall_margin = getOr(yc, "stall_margin_deg", 2.0) * kDeg;
  cb.theta_min = getOr(yc, "theta_min_deg", 3.0) * kDeg;
  cb.h_noseup = getOr(yc, "h_noseup", 3.0);
  cb.V_td_max = getOr(yc, "V_td_max", 14.0);
  cb.g_eff = getOr(yc, "g_eff", 16.0);
  cb.Tmax = getOr(yc, "Tmax", ac_.thrust.T_static);
  // Opt-in energy floor (powered stall recovery; default OFF -> set unchanged).
  cb.energy_floor = getOrB(yc, "energy_floor", false);
  cb.energy_floor_hard = getOrB(yc, "energy_floor_hard", false);
  cb.V_floor = getOr(yc, "V_floor", cb.V_floor);
  cb.g_floor = getOr(yc, "g_floor", cb.g_floor);
  // Impact-load barrier (NACA TN 1516). beta given in degrees in YAML.
  cb.n_limit = getOr(yc, "n_limit", 3.0);
  cb.beta = getOr(yc, "beta_deg", 22.5) * kDeg;
  cb.rho_water = getOr(yc, "rho_water", 1000.0);
  cb.Nb = getOr(yc, "Nb", 10.0);
  cb.zs = getOr(yc, "zs", 2.0);
  cb.tau_keel = getOr(yc, "tau_keel_deg", 0.0) * kDeg;
  cb.z_gate = getOr(yc, "z_gate", 10.0);
  cb.eps_g0 = getOr(yc, "eps_g0", 0.02);
  cb.c_stall = getArr2(yc, "c_stall", {4.0, 4.0});
  cb.c_noseup = getArr2(yc, "c_noseup", {2.0, 2.0});
  cb.c_energy = getArr3(yc, "c_energy", {2.0, 2.0, 2.0});
  cb.c_energy_floor = getArr3(yc, "c_energy_floor", cb.c_energy_floor);
  cb.c_thrust_min = getArr2(yc, "c_thrust_min", {4.0, 4.0});
  cb.c_thrust_max = getArr2(yc, "c_thrust_max", {4.0, 4.0});
  cb.c_impact = getArr2(yc, "c_impact", {2.0, 2.0});
  // Control tracking weights (identity default; range-scale to let thrust act).
  cb.w_de = getOr(yc, "w_de", cb.w_de);
  cb.w_Tddot = getOr(yc, "w_Tddot", cb.w_Tddot);
  // Per-constraint quadratic slack weights (soft rows only).
  cb.w_slack_stall = getOr(yc, "w_slack_stall", 1.0e5);
  cb.w_slack_energy = getOr(yc, "w_slack_energy", 1.0e4);
  cb.w_slack_energy_floor = getOr(yc, "w_slack_energy_floor", cb.w_slack_energy_floor);
  cb.w_slack_noseup = getOr(yc, "w_slack_noseup", 1.0e3);
  cb.w_slack_impact = getOr(yc, "w_slack_impact", 1.0e4);
  cb.h_meas_seed = static_cast<unsigned int>(
      getOr(yc, "h_meas_seed", static_cast<double>(cb.h_meas_seed)));
  cb.de_min = getOr(yc, "de_min_deg", -28.6) * kDeg;
  cb.de_max = getOr(yc, "de_max_deg", 28.6) * kDeg;
  cb.Tddot_min = getOr(yc, "Tddot_min", -500.0);
  cb.Tddot_max = getOr(yc, "Tddot_max", 500.0);
  cb.h_meas_stddev = getOr(yc, "h_meas_stddev", 0.0);
  cb.h_lpf_tau = getOr(yc, "h_lpf_tau", 0.0);

  sc_.cbf_enabled = getOrB(root, "cbf_enabled", true);
  cb.enabled = cb.enabled && sc_.cbf_enabled;

  // Optional scenario-level viscous-stall override (the plant default lives in
  // aircraft.yaml, OFF). Lets a stall-demo scenario enable the NACA 4414 stall
  // model without a duplicate aircraft config.
  YAML::Node yst = root["stall"];
  if (yst) {
    ac_.stall.enabled = getOrB(yst, "enabled", ac_.stall.enabled);
    ac_.stall.severity = getOr(yst, "severity", ac_.stall.severity);
  }

  // --- MIL-F-8785C discrete wind gust ("1-cosine"), plant disturbance only.
  // The nominal controller and the CBF-QP keep the still-air model, so an
  // enabled gust probes the filter's robustness to unmeasured wind.
  YAML::Node yw = root["wind"];
  DiscreteGustConfig& wd = sc_.wind;
  wd.enabled = getOrB(yw, "enabled", false);
  wd.t_start = getOr(yw, "t_start", 0.0);
  wd.u.amp = getOr(yw, "u_amp", 0.0);
  wd.u.len = getOr(yw, "u_len", 120.0);
  wd.w.amp = getOr(yw, "w_amp", 0.0);
  wd.w.len = getOr(yw, "w_len", 120.0);

  // --- Airy/JONSWAP surface-wave field (see water_waves.hpp), plant-side
  // truth only: touchdown surface, radar-altimeter reference, and slam-load
  // truth diagnostics. The nominal and the CBF-QP keep the flat-water model,
  // so an enabled seaway probes the filter's smooth-water assumption.
  YAML::Node yv = root["waves"];
  WaveConfig& wv = sc_.waves;
  wv.enabled = getOrB(yv, "enabled", false);
  wv.regular = getOrB(yv, "regular", false);
  wv.Hs = getOr(yv, "Hs", wv.Hs);
  wv.Tp = getOr(yv, "Tp", wv.Tp);
  wv.gamma = getOr(yv, "gamma", wv.gamma);
  wv.n = static_cast<int>(getOr(yv, "n", wv.n));
  wv.seed = static_cast<unsigned>(getOr(yv, "seed", wv.seed));
  wv.phase_deg = getOr(yv, "phase_deg", wv.phase_deg);
  wv.contact_len = getOr(yv, "contact_len", wv.contact_len);
  if (yv && yv["direction"])
    wv.dir = yv["direction"].as<std::string>() == "following" ? 1.0 : -1.0;

  sc_.dt = getOr(root, "dt", 0.01);
  sc_.t_max = getOr(root, "t_max", 60.0);

  // --- Initial condition: steady LEVEL-flight trim of the lon model itself.
  // Starting from a true equilibrium (gamma=0, Vdot=gammadot=qdot=0) reflects an
  // aircraft cruising level before it captures the approach, and avoids the
  // startup transient of seeding the IC from the (different) body-axis trim. The
  // nominal then commands gamma_ref = gamma_app and the aircraft pushes over
  // into the approach.
  YAML::Node yi = root["initial"];
  const double h0 = getOr(yi, "h0", 40.0);
  const double V0 = V_app + getOr(yi, "dV", 0.0);
  const LonTrimResult lt = lonTrim(table_, *mixing_, ac_, V0, 0.0);
  stats_.ic_trim_converged = lt.converged;
  if (!lt.converged)
    std::cerr << "[lon_sim] warning: level-flight trim did not converge (res="
              << lt.residual << ")\n";
  const double th0 = lt.theta + getOr(yi, "dtheta_deg", 0.0) * kDeg;
  sc_.X0 << h0, V0, 0.0, th0, 0.0, lt.T, 0.0;
}

LonTouchdown LonSim::run(const std::string& csv_path) {
  const double dt = sc_.dt;
  const int nsteps = static_cast<int>(sc_.t_max / dt);

  // Fresh stats per run; the IC-trim flag was set by the constructor and holds
  // for every run of this sim instance.
  const bool ic_trim_ok = stats_.ic_trim_converged;
  stats_ = LonRunStats{};
  stats_.ic_trim_converged = ic_trim_ok;

  LonNominal nominal(sc_.nominal);
  LonCBFFilter filter(sc_.cbf);
  // Diagnostic barriers (stall / nose-up / energy) are rebuilt per-step below,
  // where the frozen aero is available (mirrors the impact barrier).

  // Altitude-sensor model: the CBF acts on a noisy, low-pass-filtered measurement
  // h_filt(h + N(0,sigma^2)), while the plant (RK4) and the CSV diagnostics keep
  // the true state. Fixed seed so a given scenario is reproducible run-to-run
  // (deterministic plant + noise).
  std::mt19937 h_rng(sc_.cbf.h_meas_seed);
  std::normal_distribution<double> h_noise(0.0, sc_.cbf.h_meas_stddev);
  // First-order low-pass: alpha = dt/(tau+dt); tau=0 => alpha=1 => pass-through.
  // Seed the filter state at the true initial altitude (no startup transient).
  const double h_lpf_alpha =
      sc_.cbf.h_lpf_tau > 0.0 ? dt / (sc_.cbf.h_lpf_tau + dt) : 1.0;

  // Surface-wave field (one locked, seeded sea; empty/zero when disabled). The
  // altimeter chain below measures clearance over the INSTANTANEOUS surface
  // (radar-altimeter physics), so h_filt is seeded at the initial clearance.
  const WaveField wf = makeWaveField(sc_.waves);
  double h_filt = sc_.X0[LH] - wf.eta(0.0, 0.0);

  // Gust penetration distance [m] (MIL-F-8785C discrete gust; 0 until onset).
  double x_gust = 0.0;
  // Earth-frame ground track [m]: the wave-field abscissa under the aircraft.
  // Feeds only eta(x, t) and the diagnostics, never the dynamics, so it is
  // advanced by trapezoid outside the RK4 (error ~ mm against wavelengths of
  // metres). xdot = V cos(gamma) + W_u (ground speed incl. gust transport).
  double x_pos = 0.0;

  std::ofstream csv(csv_path);
  csv << "t,h,V,gamma_deg,theta_deg,q,T,Tdot,alpha_deg,sink,de,Tddot,"
         "de_nom,Tddot_nom,theta_cmd_deg,recovered,"
         "b_stall,psi1_stall,psi2_stall,res_stall,"
         "b_noseup,psi1_noseup,psi2_noseup,res_noseup,"
         "b_energy,psi1_energy,psi2_energy,psi3_energy,res_energy,"
         "res_tmin,res_tmax,"
         "b_impact,n_peak,kappa_imp,psi1_imp,psi2_imp,res_imp,h_meas,h_filt,"
         "n_rows_dropped,hard_dropped,CL,CD,dCL_stall,W_u,W_h,x_gust,"
         "b_efloor,psi1_efloor,psi2_efloor,psi3_efloor,res_efloor,"
         "x,eta,eta_slope\n";

  LonStateVec X = sc_.X0;
  LonTouchdown td;
  double t = 0.0;
  int recoveries = 0;
  // Track the minima of the HOCBF nested functions psi_i (the real
  // forward-invariance condition is psi_i >= 0 for all i, not just b = psi_0).
  double min_psi1_st = 1e30, min_psi2_st = 1e30;              // stall (whole flight)
  double min_psi1_nu = 1e30, min_psi2_nu = 1e30;              // nose-up (h<h_noseup)
  double min_psi1_en = 1e30, min_psi2_en = 1e30, min_psi3_en = 1e30;  // energy (deg 3)
  double min_psi1_ef = 1e30, min_psi2_ef = 1e30, min_psi3_ef = 1e30;  // energy floor (opt-in)
  double min_psi1i = 1e30, min_psi2i = 1e30;  // impact (over the active z<z_gate window)
  // Diagnostic-integrity + filter-health accounting: non-finite psi samples must
  // not be silently swallowed by std::min (NaN < x is false), and dropped HARD
  // barrier rows must be visible rather than masquerading as recoveries=0.
  int nan_stall = 0, nan_noseup = 0, nan_energy = 0, nan_efloor = 0, nan_imp = 0;
  int steps_rows_dropped = 0, steps_hard_dropped = 0;
  auto gmin = [](double& mn, double v) { if (std::isfinite(v)) mn = std::min(mn, v); };

  for (int k = 0; k <= nsteps; ++k) {
    const LonCtrlVec U_nom = nominal.step(X, dt);
    // Wave surface under the aircraft this step (identically 0 on flat water).
    // The slope is the mean tilt across the hull contact length -- the slam
    // geometry a keel of that length meets -- not the ripple-dominated point
    // slope (see water_waves.hpp::slopeMean).
    const double eta_now = wf.eta(x_pos, t);
    const double eta_x = wf.slopeMean(x_pos, t, sc_.waves.contact_len);
    // Sensor chain: radar-altimeter clearance over the instantaneous surface
    // (true h - eta) -> add noise -> low-pass -> CBF. The filter therefore
    // flies the surface-relative altitude but stays wave-BLIND otherwise (its
    // barrier model keeps flat water; it never sees the slope). The true X is
    // still integrated below and used for the diagnostics/touchdown test.
    double h_meas = X[LH] - eta_now;
    if (sc_.cbf.h_meas_stddev > 0.0) h_meas += h_noise(h_rng);
    h_filt += h_lpf_alpha * (h_meas - h_filt);
    LonStateVec X_meas = X;
    X_meas[LH] = h_filt;
    const LonCtrlVec U = filter.filter(U_nom, X_meas, table_, *mixing_, ac_);
    if (filter.lastRecovery()) ++recoveries;
    const int n_rows_dropped = filter.lastDroppedRows();
    const bool hard_dropped = filter.lastHardDropped();
    if (n_rows_dropped > 0) ++steps_rows_dropped;
    if (hard_dropped) ++steps_hard_dropped;

    const double alpha = X[LTH] - X[LGAM];
    // True inertial descent rate -hdot (positive down): air-relative sink plus
    // the vertical gust transport. Identical to -V sin(gamma) with wind off.
    const GustWind Wg = gustWind(sc_.wind, x_gust);
    const double sink = -(X[LV] * std::sin(X[LGAM]) + Wg.w);
    // Ground speed (earth x rate): the wave-field advance under the aircraft.
    const double xdot_gnd = X[LV] * std::cos(X[LGAM]) + Wg.u;
    const auto xa = toArray(X);

    // Flight-envelope extrema + plant-validity flag for the run stats.
    stats_.max_alpha = std::max(stats_.max_alpha, alpha);
    stats_.min_V = std::min(stats_.min_V, X[LV]);
    stats_.min_T = std::min(stats_.min_T, X[LT]);
    stats_.max_T = std::max(stats_.max_T, X[LT]);
    if (std::abs(alpha) > kAlphaModelLimit) stats_.out_of_model = true;

    // HOCBF nested functions psi_1, psi_2 (linear class-K): the true safe-set
    // membership condition. psi_1 = L_f b + c1 b; psi_2 = L_f^2 b + (c1+c2)L_f b
    // + c1 c2 b.
    const AeroLocal aero = makeAeroLocal(table_, *mixing_, ac_, X[LV], alpha);
    // Realized wind-axis aero coefficients (mirrors LonDrift) for the trace,
    // including the viscous-stall lift delta dCL_stall (0 when stall disabled).
    const double mach_d = X[LV] / aero.a_sound;
    const double qhat_d = (X[LQ] * aero.cref / 2.0) / X[LV];
    const double CFx_d = aero.off_CFx + aero.dAlpha_CFx * alpha +
                         aero.dMach_CFx * mach_d + aero.dQ_CFx * qhat_d;
    const double CFz_d = aero.off_CFz + aero.dAlpha_CFz * alpha +
                         aero.dMach_CFz * mach_d + aero.dQ_CFz * qhat_d;
    double CL_d = -CFx_d * std::sin(alpha) + CFz_d * std::cos(alpha);
    double CD_d = CFx_d * std::cos(alpha) + CFz_d * std::sin(alpha) + aero.parasite_CD0;
    const double CL_att = CL_d;  // pre-stall (VSPAERO) lift, for the dCL_stall trace
    if (aero.stall_on) {         // blend toward the Viterna post-stall curve
      const double w = aero.off_w + aero.dAlpha_w * alpha;
      CL_d = (1.0 - w) * CL_d + w * (aero.off_CLp + aero.dAlpha_CLp * alpha);
      CD_d = (1.0 - w) * CD_d + w * (aero.off_CDp + aero.dAlpha_CDp * alpha);
    }
    const double dCL_stall = CL_d - CL_att;  // lift lost to stall (<= 0)
    const StallBarrier bstall = makeStallBarrier(sc_.cbf.alpha_stall, sc_.cbf.stall_margin);
    const NoseUpBarrier bnose = makeNoseUpBarrier(sc_.cbf.theta_min);
    const EnergyBarrier ben = makeEnergyBarrier(aero, sc_.cbf.V_td_max, sc_.cbf.g_eff);
    const EnergyFloorBarrier befl =
        makeEnergyFloorBarrier(aero, sc_.cbf.V_floor, sc_.cbf.g_floor);
    const auto lst = barrierLie<2>(aero, bstall, X);
    const auto lnu = barrierLie<2>(aero, bnose, X);
    const auto len = barrierLie<3>(aero, ben, X);
    const auto lef = barrierLie<3>(aero, befl, X);
    const auto& cs = sc_.cbf.c_stall;
    const auto& cn = sc_.cbf.c_noseup;
    const auto& ce = sc_.cbf.c_energy;
    const auto& cf = sc_.cbf.c_energy_floor;
    const double psi1_st = lst.Lf[1] + cs[0] * lst.Lf[0];
    const double psi2_st = lst.Lf[2] + (cs[0] + cs[1]) * lst.Lf[1] + cs[0] * cs[1] * lst.Lf[0];
    const double psi1_nu = lnu.Lf[1] + cn[0] * lnu.Lf[0];
    const double psi2_nu = lnu.Lf[2] + (cn[0] + cn[1]) * lnu.Lf[1] + cn[0] * cn[1] * lnu.Lf[0];
    const double psi1_en = len.Lf[1] + ce[0] * len.Lf[0];
    const double psi2_en = len.Lf[2] + (ce[0] + ce[1]) * len.Lf[1] + ce[0] * ce[1] * len.Lf[0];
    const double psi3_en = len.Lf[3] + (ce[0] + ce[1] + ce[2]) * len.Lf[2]
                         + (ce[0] * ce[1] + ce[0] * ce[2] + ce[1] * ce[2]) * len.Lf[1]
                         + ce[0] * ce[1] * ce[2] * len.Lf[0];
    const double psi1_ef = lef.Lf[1] + cf[0] * lef.Lf[0];
    const double psi2_ef = lef.Lf[2] + (cf[0] + cf[1]) * lef.Lf[1] + cf[0] * cf[1] * lef.Lf[0];
    const double psi3_ef = lef.Lf[3] + (cf[0] + cf[1] + cf[2]) * lef.Lf[2]
                         + (cf[0] * cf[1] + cf[0] * cf[2] + cf[1] * cf[2]) * lef.Lf[1]
                         + cf[0] * cf[1] * cf[2] * lef.Lf[0];
    gmin(min_psi1_st, psi1_st); gmin(min_psi2_st, psi2_st);
    gmin(min_psi1_en, psi1_en); gmin(min_psi2_en, psi2_en); gmin(min_psi3_en, psi3_en);
    if (!std::isfinite(psi1_st) || !std::isfinite(psi2_st)) ++nan_stall;
    if (!std::isfinite(psi1_en) || !std::isfinite(psi2_en) || !std::isfinite(psi3_en)) ++nan_energy;
    // Energy-floor psi minima only when its row is actually assembled (opt-in).
    if (sc_.cbf.energy_floor) {
      gmin(min_psi1_ef, psi1_ef); gmin(min_psi2_ef, psi2_ef); gmin(min_psi3_ef, psi3_ef);
      if (!std::isfinite(psi1_ef) || !std::isfinite(psi2_ef) || !std::isfinite(psi3_ef))
        ++nan_efloor;
    }
    // Nose-up psi minima only over its active window (h < h_noseup), where the row
    // is actually assembled.
    if (X[LH] < sc_.cbf.h_noseup) {
      gmin(min_psi1_nu, psi1_nu); gmin(min_psi2_nu, psi2_nu);
      if (!std::isfinite(psi1_nu) || !std::isfinite(psi2_nu)) ++nan_noseup;
    }

    // Impact-load barrier (degree 2). psi minima are tracked only over the active
    // window (z < z_gate), where the row is actually assembled/enforced.
    const ImpactLoadBarrier bimp =
        makeImpactLoadBarrier(aero, sc_.cbf.n_limit, sc_.cbf.beta, sc_.cbf.rho_water,
                              sc_.cbf.Nb, sc_.cbf.zs, sc_.cbf.tau_keel, sc_.cbf.eps_g0, X);
    const auto lii = barrierLie<2>(aero, bimp, X);
    const auto& ci = sc_.cbf.c_impact;
    const double psi1i = lii.Lf[1] + ci[0] * lii.Lf[0];
    const double psi2i = lii.Lf[2] + (ci[0] + ci[1]) * lii.Lf[1] + ci[0] * ci[1] * lii.Lf[0];
    // Track minima only where the filter actually assembles the row (the active,
    // model-valid window): below z_gate, descending, positive trim.
    const bool impact_active = X[LH] < sc_.cbf.z_gate && -X[LGAM] > 0.0 &&
                               (X[LTH] - sc_.cbf.tau_keel) > 0.0;
    const double b_imp_val = bimp(xa);
    if (impact_active) {
      gmin(min_psi1i, psi1i); gmin(min_psi2i, psi2i);
      gmin(stats_.min_b_impact_active, b_imp_val);
      if (!std::isfinite(psi1i) || !std::isfinite(psi2i)) ++nan_imp;
    }
    // n_peak and kappa for the trace, recomputed from the frozen barrier.
    const double gamma0_l = -X[LGAM];
    const double tau_l = X[LTH] - sc_.cbf.tau_keel;
    const double sg0_l = std::sqrt(std::sin(gamma0_l) * std::sin(gamma0_l) +
                                   sc_.cbf.eps_g0 * sc_.cbf.eps_g0);
    const double kappa_l = std::sin(tau_l) * std::cos(tau_l + gamma0_l) / sg0_l;
    const double ydot0_l = -X[LV] * std::sin(X[LGAM]);
    const double n_peak = (bimp.K00 + bimp.dK0_dtau * (tau_l - bimp.tau0)) *
                          ydot0_l * ydot0_l *
                          (bimp.Clf0 + bimp.dClf_dk * (kappa_l - bimp.kappa0));

    // Active-set diagnostics: constraint residual rhs - a.U for each barrier row
    // (residual ~ 0 => that barrier is binding/active and shaping the control).
    std::vector<double> Lfst(lst.Lf.begin(), lst.Lf.end());
    std::vector<double> Lfnu(lnu.Lf.begin(), lnu.Lf.end());
    std::vector<double> Lfen(len.Lf.begin(), len.Lf.end());
    const HocbfRow rst = hocbfRow(Lfst, lst.LgLf, {cs[0], cs[1]});
    const HocbfRow rnu = hocbfRow(Lfnu, lnu.LgLf, {cn[0], cn[1]});
    const HocbfRow ren = hocbfRow(Lfen, len.LgLf, {ce[0], ce[1], ce[2]});
    std::vector<double> Lfef(lef.Lf.begin(), lef.Lf.end());
    const HocbfRow ref = hocbfRow(Lfef, lef.LgLf, {cf[0], cf[1], cf[2]});
    const HocbfRow rtl = thrustMinRow(X, sc_.cbf.c_thrust_min[0], sc_.cbf.c_thrust_min[1]);
    const HocbfRow rtu = thrustMaxRow(X, sc_.cbf.Tmax, sc_.cbf.c_thrust_max[0], sc_.cbf.c_thrust_max[1]);
    std::vector<double> Lfi(lii.Lf.begin(), lii.Lf.end());
    const HocbfRow ri = hocbfRow(Lfi, lii.LgLf, {ci[0], ci[1]});
    auto resid = [&](const HocbfRow& r) {
      return r.rhs - (r.a(0, LDE) * U[LDE] + r.a(0, LTDDOT) * U[LTDDOT]);
    };
    const double res_st = resid(rst), res_nu = resid(rnu), res_en = resid(ren);
    const double res_ef = resid(ref);
    const double res_tl = resid(rtl), res_tu = resid(rtu);
    const double res_i = resid(ri);
    // Enforced impact-row residual with the APPLIED control (>= 0 up to QP
    // tolerance when the hard row held), tracked over the active window.
    if (impact_active) gmin(stats_.min_res_impact_active, res_i);

    csv << t << ',' << X[LH] << ',' << X[LV] << ',' << X[LGAM] / kDeg << ','
        << X[LTH] / kDeg << ',' << X[LQ] << ',' << X[LT] << ',' << X[LTDOT] << ','
        << alpha / kDeg << ',' << sink << ',' << U[LDE] << ',' << U[LTDDOT] << ','
        << U_nom[LDE] << ',' << U_nom[LTDDOT] << ',' << nominal.thetaCmd() / kDeg << ','
        << (filter.lastRecovery() ? 1 : 0) << ','
        << bstall(xa) << ',' << psi1_st << ',' << psi2_st << ',' << res_st << ','
        << bnose(xa) << ',' << psi1_nu << ',' << psi2_nu << ',' << res_nu << ','
        << ben(xa) << ',' << psi1_en << ',' << psi2_en << ',' << psi3_en << ',' << res_en << ','
        << res_tl << ',' << res_tu << ','
        << b_imp_val << ',' << n_peak << ',' << kappa_l << ',' << psi1i << ',' << psi2i
        << ',' << res_i << ',' << h_meas << ',' << h_filt << ',' << n_rows_dropped << ','
        << (hard_dropped ? 1 : 0) << ',' << CL_d << ',' << CD_d << ',' << dCL_stall << ','
        << Wg.u << ',' << Wg.w << ',' << x_gust << ','
        << befl(xa) << ',' << psi1_ef << ',' << psi2_ef << ',' << psi3_ef << ','
        << res_ef << ',' << x_pos << ',' << eta_now << ',' << eta_x << '\n';

    if (X[LH] <= eta_now && k > 0) {
      td.reached = true;
      td.t = t; td.sink = sink; td.V = X[LV]; td.theta = X[LTH]; td.gamma = X[LGAM];
      td.h = X[LH]; td.x = x_pos; td.eta = eta_now; td.slope = eta_x;
      // Surface-relative closure -d/dt(h - eta) = sink + eta_x xdot + eta_t:
      // what the keel actually closes at when the face is rising into the path.
      td.sink_rel = sink + eta_x * xdot_gnd + wf.etaDot(x_pos, t);
      // TN 1516 slam-load truth at contact, flat-water- vs wave-referenced
      // (tau/gamma0 tilted by the local surface angle alpha_s, closure onto
      // the moving surface -- the spec's rough-water referencing). The filter
      // enforced only the flat-water model; the ratio is what it missed.
      const double alpha_s = std::atan(eta_x);
      td.n_peak_flat = impactNPeakExact(
          aero.mass, aero.g, sc_.cbf.beta, sc_.cbf.rho_water, sc_.cbf.eps_g0,
          X[LTH] - sc_.cbf.tau_keel, -X[LGAM], std::max(0.0, sink));
      td.n_peak_wave = impactNPeakExact(
          aero.mass, aero.g, sc_.cbf.beta, sc_.cbf.rho_water, sc_.cbf.eps_g0,
          X[LTH] - alpha_s - sc_.cbf.tau_keel, alpha_s - X[LGAM],
          std::max(0.0, td.sink_rel));
      break;
    }
    if (sc_.wind.enabled) {
      const WindStep ws = rk4Wind(table_, *mixing_, ac_, sc_.wind, t, X, x_gust, U, dt);
      X = ws.X;
      x_gust = ws.xg;
    } else {
      // Still-air path kept separate so wind-off runs are bit-identical to the
      // pre-gust plant (mirrors the stall_on pattern).
      X = rk4(table_, *mixing_, ac_, X, U, dt);
    }
    // Ground-track trapezoid across the step (endpoint ground speeds).
    x_pos += 0.5 * dt *
             (xdot_gnd + X[LV] * std::cos(X[LGAM]) + gustWind(sc_.wind, x_gust).u);
    t += dt;
    ++stats_.steps;
  }
  td.recoveries = recoveries;
  td.steps_hard_dropped = steps_hard_dropped;

  // Copy the loop tallies into the run record (the locals stay the print
  // sources below, so console output is unchanged).
  stats_.recoveries = recoveries;
  stats_.steps_rows_dropped = steps_rows_dropped;
  stats_.steps_hard_dropped = steps_hard_dropped;
  stats_.min_psi1_stall = min_psi1_st;   stats_.min_psi2_stall = min_psi2_st;
  stats_.min_psi1_noseup = min_psi1_nu;  stats_.min_psi2_noseup = min_psi2_nu;
  stats_.min_psi1_energy = min_psi1_en;  stats_.min_psi2_energy = min_psi2_en;
  stats_.min_psi3_energy = min_psi3_en;
  stats_.min_psi1_efloor = min_psi1_ef;  stats_.min_psi2_efloor = min_psi2_ef;
  stats_.min_psi3_efloor = min_psi3_ef;
  stats_.min_psi1_impact = min_psi1i;    stats_.min_psi2_impact = min_psi2i;
  stats_.nan_stall = nan_stall;   stats_.nan_noseup = nan_noseup;
  stats_.nan_energy = nan_energy; stats_.nan_efloor = nan_efloor;
  stats_.nan_impact = nan_imp;
  stats_.t_end = t;

  std::cout << "=== Augmented longitudinal landing ===\n";
  std::cout << "trim theta_cmd seed, T_set=" << sc_.nominal.T_set << " N, gamma_ref="
            << sc_.nominal.gamma_ref / kDeg << " deg\n";
  if (sc_.wind.enabled) {
    std::cout << "  MIL-F-8785C discrete gust (plant-only, unmeasured): t_start="
              << sc_.wind.t_start << " s  u: Vm=" << sc_.wind.u.amp << " m/s dm="
              << sc_.wind.u.len << " m (+tailwind)  w: Vm=" << sc_.wind.w.amp
              << " m/s dm=" << sc_.wind.w.len << " m (+updraft)\n";
  }
  if (sc_.waves.enabled) {
    std::cout << "  Surface waves (plant-only, wave-blind filter): ";
    if (sc_.waves.regular)
      std::cout << "regular";
    else if (sc_.waves.gamma <= 1.0)
      std::cout << "Bretschneider (STANAG 4194)";
    else
      std::cout << "JONSWAP gamma=" << sc_.waves.gamma;
    std::cout << "  Hs=" << sc_.waves.Hs << " m  Tp=" << sc_.waves.Tp << " s  "
              << (sc_.waves.dir < 0 ? "head" : "following") << " seas\n";
  }
  if (td.reached) {
    std::cout << "TOUCHDOWN  t=" << td.t << " s  sink=" << td.sink << " m/s"
              << "  V=" << td.V << " m/s  theta=" << td.theta / kDeg << " deg\n";
    std::cout << "  V_td_max cap=" << sc_.cbf.V_td_max << " m/s  -> "
              << (td.V <= sc_.cbf.V_td_max + 1e-6 ? "WITHIN" : "EXCEEDED") << "\n";
    if (sc_.waves.enabled) {
      std::cout << "  wave contact: eta=" << td.eta << " m  surface slope="
                << std::atan(td.slope) / kDeg << " deg  sink_rel=" << td.sink_rel
                << " m/s (flat-ref " << td.sink << ")\n";
      std::cout << "  TN1516 n_peak truth: flat-ref=" << td.n_peak_flat
                << " g  wave-ref=" << td.n_peak_wave << " g";
      if (td.n_peak_flat > 1e-9)
        std::cout << "  (x" << td.n_peak_wave / td.n_peak_flat << ")";
      std::cout << "\n";
    }
  } else {
    std::cout << "No touchdown within t_max=" << sc_.t_max << " s (final h="
              << X[LH] << ")\n";
  }
  std::cout << "  CBF feasibility recoveries: " << recoveries << "\n";
  // Row-assembly faults: a dropped HARD row is an unannunciated loss of a safety
  // guarantee, so surface it explicitly -- recoveries=0 alone does NOT imply health
  // (dropping the binding row makes the QP trivially feasible).
  std::cout << "  CBF row-assembly faults: " << steps_rows_dropped
            << " step(s) dropped >=1 barrier row (" << steps_hard_dropped
            << " dropped a HARD row)\n";
  auto nanTag = [](int k) {
    return k > 0 ? "  [" + std::to_string(k) + " non-finite sample(s) excluded]"
                 : std::string();
  };
  std::cout << "  HOCBF nested-function minima (impact HARD must be >= 0; soft rows may dip):\n";
  std::cout << "    stall  : min psi1=" << min_psi1_st << "  min psi2=" << min_psi2_st
            << nanTag(nan_stall) << "\n";
  std::cout << "    nose-up: min psi1=" << min_psi1_nu << "  min psi2=" << min_psi2_nu
            << "  (h < h_noseup=" << sc_.cbf.h_noseup << " m)" << nanTag(nan_noseup) << "\n";
  std::cout << "    energy : min psi1=" << min_psi1_en << "  min psi2=" << min_psi2_en
            << "  min psi3=" << min_psi3_en << nanTag(nan_energy) << "\n";
  if (sc_.cbf.energy_floor) {
    std::cout << "    e-floor: min psi1=" << min_psi1_ef << "  min psi2=" << min_psi2_ef
              << "  min psi3=" << min_psi3_ef << "  (V_floor=" << sc_.cbf.V_floor
              << " m/s, g_floor=" << sc_.cbf.g_floor << ")" << nanTag(nan_efloor) << "\n";
  }
  std::cout << "    impact : min psi1=" << min_psi1i << "  min psi2=" << min_psi2i
            << "  (z < z_gate=" << sc_.cbf.z_gate << " m)" << nanTag(nan_imp) << "\n";
  std::cout << "  trace -> " << csv_path << "\n";
  return td;
}

}  // namespace autoland
