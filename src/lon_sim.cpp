#include "autoland/lon_sim.hpp"

#include <yaml-cpp/yaml.h>

#include <cmath>
#include <fstream>
#include <iostream>

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

namespace {
LonStateVec rk4(const AeroTable& t, const Mixing& mx, const AircraftConfig& cfg,
                const LonStateVec& X, const LonCtrlVec& U, double dt) {
  const LonStateVec k1 = lonXdotFull(t, mx, cfg, X, U);
  const LonStateVec k2 = lonXdotFull(t, mx, cfg, X + 0.5 * dt * k1, U);
  const LonStateVec k3 = lonXdotFull(t, mx, cfg, X + 0.5 * dt * k2, U);
  const LonStateVec k4 = lonXdotFull(t, mx, cfg, X + dt * k3, U);
  return X + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
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
  cb.descent = getOrB(yc, "descent", true);
  cb.airspeed = getOrB(yc, "airspeed", true);
  cb.thrust_limits = getOrB(yc, "thrust_limits", true);
  cb.descent_hard = getOrB(yc, "descent_hard", true);
  cb.airspeed_hard = getOrB(yc, "airspeed_hard", false);
  cb.v_safe = getOr(yc, "v_safe", 0.6);
  cb.a_brk = getOr(yc, "a_brk", 3.0);
  cb.Vmin = getOr(yc, "Vmin", 0.85 * V_app);
  cb.Tmax = getOr(yc, "Tmax", ac_.thrust.T_static);
  cb.c_descent = getArr3(yc, "c_descent", {2.0, 2.0, 2.0});
  cb.c_airspeed = getArr3(yc, "c_airspeed", {2.0, 2.0, 2.0});
  cb.c_thrust_min = getArr2(yc, "c_thrust_min", {4.0, 4.0});
  cb.c_thrust_max = getArr2(yc, "c_thrust_max", {4.0, 4.0});
  cb.w_de = getOr(yc, "w_de", 1.0);
  cb.w_Tddot = getOr(yc, "w_Tddot", 1.0);
  cb.slack_penalty = getOr(yc, "slack_penalty", 1.0e4);
  cb.de_min = getOr(yc, "de_min_deg", -28.6) * kDeg;
  cb.de_max = getOr(yc, "de_max_deg", 28.6) * kDeg;
  cb.Tddot_min = getOr(yc, "Tddot_min", -500.0);
  cb.Tddot_max = getOr(yc, "Tddot_max", 500.0);

  sc_.cbf_enabled = getOrB(root, "cbf_enabled", true);
  cb.enabled = cb.enabled && sc_.cbf_enabled;

  sc_.dt = getOr(root, "dt", 0.01);
  sc_.t_max = getOr(root, "t_max", 60.0);

  // --- Initial condition ---
  YAML::Node yi = root["initial"];
  const double h0 = getOr(yi, "h0", 40.0);
  const double V0 = V_app + getOr(yi, "dV", 0.0);
  const double g0 = gamma_app + getOr(yi, "dgamma_deg", 0.0) * kDeg;
  const double th0 = tr.theta + getOr(yi, "dtheta_deg", 0.0) * kDeg;
  sc_.X0 << h0, V0, g0, th0, 0.0, nom.T_set, 0.0;
}

LonTouchdown LonSim::run(const std::string& csv_path) {
  const double dt = sc_.dt;
  const int nsteps = static_cast<int>(sc_.t_max / dt);

  LonNominal nominal(sc_.nominal);
  LonCBFFilter filter(sc_.cbf);
  DescentBarrier bdesc{sc_.cbf.v_safe * sc_.cbf.v_safe, sc_.cbf.a_brk};
  AirspeedBarrier bair{sc_.cbf.Vmin};

  std::ofstream csv(csv_path);
  csv << "t,h,V,gamma_deg,theta_deg,q,T,Tdot,alpha_deg,sink,de,Tddot,"
         "de_nom,Tddot_nom,b_descent,b_airspeed,theta_cmd_deg,recovered,"
         "psi1_desc,psi2_desc,psi1_air,psi2_air,"
         "res_desc,res_air,res_tmin,res_tmax\n";

  LonStateVec X = sc_.X0;
  LonTouchdown td;
  double t = 0.0;
  int recoveries = 0;
  // Track the minima of the HOCBF nested functions psi_i (the real
  // forward-invariance condition is psi_i >= 0 for all i, not just b = psi_0).
  double min_psi1d = 1e30, min_psi2d = 1e30, min_psi1a = 1e30, min_psi2a = 1e30;

  for (int k = 0; k <= nsteps; ++k) {
    const LonCtrlVec U_nom = nominal.step(X, dt);
    const LonCtrlVec U = filter.filter(U_nom, X, table_, *mixing_, ac_);
    if (filter.lastRecovery()) ++recoveries;

    const double alpha = X[LTH] - X[LGAM];
    const double sink = -X[LV] * std::sin(X[LGAM]);  // positive down
    const auto xa = toArray(X);

    // HOCBF nested functions psi_1, psi_2 (linear class-K): the true safe-set
    // membership condition. psi_1 = L_f b + c1 b; psi_2 = L_f^2 b + (c1+c2)L_f b
    // + c1 c2 b.
    const AeroLocal aero = makeAeroLocal(table_, *mixing_, ac_, X[LV], alpha);
    const auto ldd = barrierLie<3>(aero, bdesc, X);
    const auto lda = barrierLie<3>(aero, bair, X);
    const auto& cd = sc_.cbf.c_descent;
    const auto& ca = sc_.cbf.c_airspeed;
    const double psi1d = ldd.Lf[1] + cd[0] * ldd.Lf[0];
    const double psi2d = ldd.Lf[2] + (cd[0] + cd[1]) * ldd.Lf[1] + cd[0] * cd[1] * ldd.Lf[0];
    const double psi1a = lda.Lf[1] + ca[0] * lda.Lf[0];
    const double psi2a = lda.Lf[2] + (ca[0] + ca[1]) * lda.Lf[1] + ca[0] * ca[1] * lda.Lf[0];
    min_psi1d = std::min(min_psi1d, psi1d);
    min_psi2d = std::min(min_psi2d, psi2d);
    min_psi1a = std::min(min_psi1a, psi1a);
    min_psi2a = std::min(min_psi2a, psi2a);

    // Active-set diagnostics: constraint residual rhs - a.U for each barrier row
    // (residual ~ 0 => that barrier is binding/active and shaping the control).
    std::vector<double> Lfd(ldd.Lf.begin(), ldd.Lf.end());
    std::vector<double> Lfa(lda.Lf.begin(), lda.Lf.end());
    const HocbfRow rd = hocbfRow(Lfd, ldd.LgLf, {cd[0], cd[1], cd[2]});
    const HocbfRow ra = hocbfRow(Lfa, lda.LgLf, {ca[0], ca[1], ca[2]});
    const HocbfRow rtl = thrustMinRow(X, sc_.cbf.c_thrust_min[0], sc_.cbf.c_thrust_min[1]);
    const HocbfRow rtu = thrustMaxRow(X, sc_.cbf.Tmax, sc_.cbf.c_thrust_max[0], sc_.cbf.c_thrust_max[1]);
    auto resid = [&](const HocbfRow& r) {
      return r.rhs - (r.a(0, LDE) * U[LDE] + r.a(0, LTDDOT) * U[LTDDOT]);
    };
    const double res_d = resid(rd), res_a = resid(ra), res_tl = resid(rtl), res_tu = resid(rtu);

    csv << t << ',' << X[LH] << ',' << X[LV] << ',' << X[LGAM] / kDeg << ','
        << X[LTH] / kDeg << ',' << X[LQ] << ',' << X[LT] << ',' << X[LTDOT] << ','
        << alpha / kDeg << ',' << sink << ',' << U[LDE] << ',' << U[LTDDOT] << ','
        << U_nom[LDE] << ',' << U_nom[LTDDOT] << ',' << bdesc(xa) << ','
        << bair(xa) << ',' << nominal.thetaCmd() / kDeg << ','
        << (filter.lastRecovery() ? 1 : 0) << ',' << psi1d << ',' << psi2d << ','
        << psi1a << ',' << psi2a << ',' << res_d << ',' << res_a << ',' << res_tl
        << ',' << res_tu << '\n';

    if (X[LH] <= 0.0 && k > 0) {
      td.reached = true;
      td.t = t; td.sink = sink; td.V = X[LV]; td.theta = X[LTH]; td.gamma = X[LGAM];
      break;
    }
    X = rk4(table_, *mixing_, ac_, X, U, dt);
    t += dt;
  }

  std::cout << "=== Augmented longitudinal landing ===\n";
  std::cout << "trim theta_cmd seed, T_set=" << sc_.nominal.T_set << " N, gamma_ref="
            << sc_.nominal.gamma_ref / kDeg << " deg\n";
  if (td.reached) {
    std::cout << "TOUCHDOWN  t=" << td.t << " s  sink=" << td.sink << " m/s"
              << "  V=" << td.V << " m/s  theta=" << td.theta / kDeg << " deg\n";
    std::cout << "  v_safe budget=" << sc_.cbf.v_safe << " m/s  -> "
              << (td.sink <= sc_.cbf.v_safe + 1e-6 ? "WITHIN" : "EXCEEDED") << "\n";
  } else {
    std::cout << "No touchdown within t_max=" << sc_.t_max << " s (final h="
              << X[LH] << ")\n";
  }
  std::cout << "  CBF feasibility recoveries: " << recoveries << "\n";
  std::cout << "  HOCBF nested-function minima (must be >= 0 for the guarantee):\n";
  std::cout << "    descent : min psi1=" << min_psi1d << "  min psi2=" << min_psi2d << "\n";
  std::cout << "    airspeed: min psi1=" << min_psi1a << "  min psi2=" << min_psi2a << "\n";
  std::cout << "  trace -> " << csv_path << "\n";
  return td;
}

}  // namespace autoland
