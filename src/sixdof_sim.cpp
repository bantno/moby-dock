#include "autoland/sixdof_sim.hpp"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "autoland/impact_barrier.hpp"

namespace autoland {
namespace {
constexpr double kDeg = M_PI / 180.0;

double getOr(const YAML::Node& n, const std::string& key, double def) {
  return (n && n[key]) ? n[key].as<double>() : def;
}
bool getOrB(const YAML::Node& n, const std::string& key, bool def) {
  return (n && n[key]) ? n[key].as<bool>() : def;
}

// Earth-frame gust (x north tailwind+, y east+, h updraft+) rotated into body
// axes: W_b = R_nb^T * (W_u, W_v, -W_w). Same DCM rows as Dynamics::xdot.
// NOTE: the gust parameter must NOT be named 'W' -- that shadows the State::W
// enum index (cf. the 'Vt' note in dynamics.cpp).
Eigen::Vector3d windBody(const StateVec& x, const GustWind& Wg) {
  const double ct = std::cos(x[THETA]), st = std::sin(x[THETA]);
  const double cp = std::cos(x[PHI]), sp = std::sin(x[PHI]);
  const double cy = std::cos(x[PSI]), sy = std::sin(x[PSI]);
  const double Wn = Wg.u, We = Wg.v, Wd = -Wg.w;
  return {ct * cy * Wn + ct * sy * We - st * Wd,
          (sp * st * cy - cp * sy) * Wn + (sp * st * sy + cp * cy) * We +
              sp * ct * Wd,
          (cp * st * cy + sp * sy) * Wn + (cp * st * sy - sp * cy) * We +
              cp * ct * Wd};
}

// Air-relative flight condition (airspeed, alpha, beta) under the gust field.
struct AirData {
  double V, alpha, beta;
};
AirData airData(const StateVec& x, const GustWind& Wg) {
  const Eigen::Vector3d Wb = windBody(x, Wg);
  const double ua = x[U] - Wb[0], va = x[V] - Wb[1], wa = x[W] - Wb[2];
  const double V = std::max(1e-3, std::sqrt(ua * ua + va * va + wa * wa));
  return {V, std::atan2(wa, ua), std::asin(std::clamp(va / V, -1.0, 1.0))};
}

// Inertial climb rate (positive up) and north/east ground velocities from the
// state kinematics (rows of R_nb; identical to Dynamics::xdot's hdot/ydot).
struct GroundKinematics {
  double hdot, xdot_n, ydot_e;
};
GroundKinematics groundKinematics(const StateVec& x) {
  const double ct = std::cos(x[THETA]), st = std::sin(x[THETA]);
  const double cp = std::cos(x[PHI]), sp = std::sin(x[PHI]);
  const double cy = std::cos(x[PSI]), sy = std::sin(x[PSI]);
  const double u = x[U], v = x[V], w = x[W];
  return {u * st - v * sp * ct - w * cp * ct,
          u * ct * cy + v * (sp * st * cy - cp * sy) +
              w * (cp * st * cy + sp * sy),
          u * ct * sy + v * (sp * st * sy + cp * cy) +
              w * (cp * st * sy - sp * cy)};
}

// RK4 step of the joint system [x; x_gust]: the gust penetration distance is
// a genuine state (xdot = airspeed after onset), integrated with the aircraft
// like lon_sim's rk4Wind. With wind disabled every gust term is exactly zero,
// so this single path is bit-identical to a still-air RK4 of Dynamics::xdot.
struct JointState {
  StateVec x;
  double xg;
};
JointState rk4WindStep(const Dynamics& dyn, const DiscreteGustConfig& g,
                       double t, const StateVec& x, double xg,
                       const CtrlVec& u, double dt) {
  auto deriv = [&](double tl, const StateVec& xl, double xgl, StateVec& xd,
                   double& xgd) {
    const GustWind Wg = gustWind(g, xgl);
    xgd = gustXdot(g, tl, airData(xl, Wg).V);
    xd = dyn.xdot(xl, u, Eigen::Vector3d(Wg.u, Wg.v, Wg.w));
  };
  StateVec k1, k2, k3, k4;
  double g1, g2, g3, g4;
  deriv(t, x, xg, k1, g1);
  deriv(t + 0.5 * dt, x + 0.5 * dt * k1, xg + 0.5 * dt * g1, k2, g2);
  deriv(t + 0.5 * dt, x + 0.5 * dt * k2, xg + 0.5 * dt * g2, k3, g3);
  deriv(t + dt, x + dt * k3, xg + dt * g3, k4, g4);
  return {x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4),
          xg + (dt / 6.0) * (g1 + 2.0 * g2 + 2.0 * g3 + g4)};
}
}  // namespace

SixDofSim::SixDofSim(const std::string& stab_path,
                     const std::string& aircraft_yaml,
                     const std::string& scenario_yaml)
    : table_(AeroTable::fromFile(stab_path)),
      ac_(loadAircraftConfig(aircraft_yaml)) {
  mixing_ = std::make_unique<Mixing>(Mixing::build(ac_, table_));
  dyn_ = std::make_unique<Dynamics>(table_, *mixing_, ac_);

  YAML::Node root = YAML::LoadFile(scenario_yaml);
  sc_.V_app = getOr(root, "V_app", 18.0);
  sc_.gamma_app = getOr(root, "gamma_app_deg", -3.0) * kDeg;
  sc_.dt = getOr(root, "dt", 0.01);
  sc_.t_max = getOr(root, "t_max", 200.0);

  // Approach trim: nominal feedforwards + the initial condition.
  trim_ = trim(*dyn_, sc_.V_app, sc_.gamma_app);
  stats_.trim_converged = trim_.converged;
  if (!trim_.converged)
    std::cerr << "[sixdof_sim] warning: approach trim did not converge (res="
              << trim_.residual << ")\n";

  // --- Nominal defaults from trim (then YAML overrides) ---
  SixDofNominalConfig& nom = sc_.nominal;
  nom.V_ref = sc_.V_app;
  nom.gamma_ref = sc_.gamma_app;
  nom.theta_trim = trim_.theta;
  nom.de_trim = trim_.u[DE];
  nom.dT_trim = trim_.u[DT];
  nom.limits = ac_.limits;
  YAML::Node yn = root["nominal"];
  nom.V_ref = getOr(yn, "V_ref", nom.V_ref);
  if (yn && yn["gamma_ref_deg"]) nom.gamma_ref = yn["gamma_ref_deg"].as<double>() * kDeg;
  nom.Kp_V = getOr(yn, "Kp_V", nom.Kp_V);
  nom.Ki_V = getOr(yn, "Ki_V", nom.Ki_V);
  nom.Kp_gamma = getOr(yn, "Kp_gamma", nom.Kp_gamma);
  nom.Ki_gamma = getOr(yn, "Ki_gamma", nom.Ki_gamma);
  nom.Kv_gamma = getOr(yn, "Kv_gamma", nom.Kv_gamma);
  nom.dgamma_V_max = getOr(yn, "dgamma_V_max_deg", nom.dgamma_V_max / kDeg) * kDeg;
  nom.theta_cmd_max = getOr(yn, "theta_cmd_max_deg", nom.theta_cmd_max / kDeg) * kDeg;
  nom.Kp_theta = getOr(yn, "Kp_theta", nom.Kp_theta);
  nom.Ki_theta = getOr(yn, "Ki_theta", nom.Ki_theta);
  nom.Kq = getOr(yn, "Kq", nom.Kq);
  nom.Kp_y = getOr(yn, "Kp_y", nom.Kp_y);
  nom.Kd_y = getOr(yn, "Kd_y", nom.Kd_y);
  nom.phi_max = getOr(yn, "phi_max_deg", nom.phi_max / kDeg) * kDeg;
  nom.Kp_phi = getOr(yn, "Kp_phi", nom.Kp_phi);
  nom.Kp_p = getOr(yn, "Kp_p", nom.Kp_p);
  nom.Kr = getOr(yn, "Kr", nom.Kr);

  // --- MIL-F-8785C discrete gust (plant-side, unmeasured). enabled: one-line
  // toggle; all-zero amplitudes are equivalent to off.
  YAML::Node yw = root["wind"];
  DiscreteGustConfig& wd = sc_.wind;
  wd.enabled = getOrB(yw, "enabled", false);
  wd.t_start = getOr(yw, "t_start", 0.0);
  wd.u.amp = getOr(yw, "u_amp", 0.0);
  wd.u.len = getOr(yw, "u_len", 120.0);
  wd.v.amp = getOr(yw, "v_amp", 0.0);
  wd.v.len = getOr(yw, "v_len", 120.0);
  wd.w.amp = getOr(yw, "w_amp", 0.0);
  wd.w.len = getOr(yw, "w_len", 120.0);

  // --- Airy/JONSWAP surface waves (plant-side truth: touchdown surface +
  // contact diagnostics). enabled: one-line toggle; off => flat water.
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

  // --- Hull contact diagnostics (TN 1516 truth at touchdown). ---
  YAML::Node yh = root["hull"];
  sc_.hull.beta = getOr(yh, "beta_deg", sc_.hull.beta / kDeg) * kDeg;
  sc_.hull.rho_water = getOr(yh, "rho_water", sc_.hull.rho_water);
  sc_.hull.tau_keel = getOr(yh, "tau_keel_deg", sc_.hull.tau_keel / kDeg) * kDeg;
  sc_.hull.eps_g0 = getOr(yh, "eps_g0", sc_.hull.eps_g0);

  // --- Initial condition: approach trim + scenario offsets. ---
  YAML::Node yi = root["initial"];
  sc_.h0 = getOr(yi, "h0", 40.0);
  sc_.dV = getOr(yi, "dV", 0.0);
  sc_.dy = getOr(yi, "dy", 0.0);
  sc_.dpsi = getOr(yi, "dpsi_deg", 0.0) * kDeg;
  sc_.dtheta = getOr(yi, "dtheta_deg", 0.0) * kDeg;

  x0_ = trim_.x;
  const double vscale = (sc_.V_app + sc_.dV) / std::max(1e-6, sc_.V_app);
  x0_[U] *= vscale;
  x0_[W] *= vscale;  // v = 0 at trim (wings-level, zero-sideslip)
  x0_[THETA] += sc_.dtheta;
  x0_[PSI] += sc_.dpsi;
  x0_[H] = sc_.h0;
  x0_[Y] = sc_.dy;
}

SixDofTouchdown SixDofSim::run(const std::string& csv_path) {
  const double dt = sc_.dt;
  const int nsteps = static_cast<int>(sc_.t_max / dt);

  const bool trim_ok = stats_.trim_converged;
  stats_ = SixDofRunStats{};
  stats_.trim_converged = trim_ok;

  SixDofNominal nominal(sc_.nominal);
  const WaveField wf = makeWaveField(sc_.waves);

  // Gust penetration distance [m] and earth-frame north position [m]. x_pos
  // feeds only the wave field eta(x, t) and the diagnostics, never the
  // dynamics, so it is advanced by trapezoid outside the RK4 (cf. lon_sim).
  double x_gust = 0.0;
  double x_pos = 0.0;

  std::ofstream csv(csv_path);
  csv << "t,x,y,h,u,v,w,V_air,alpha_deg,beta_deg,p,q,r,"
         "phi_deg,theta_deg,psi_deg,gamma_deg,sink,"
         "de,da,dr,dT,theta_cmd_deg,phi_cmd_deg,"
         "W_u,W_v,W_h,x_gust,eta,eta_slope\n";

  StateVec x = x0_;
  SixDofTouchdown td;
  double t = 0.0;

  for (int k = 0; k <= nsteps; ++k) {
    const GustWind Wg = gustWind(sc_.wind, x_gust);
    const AirData ad = airData(x, Wg);
    const GroundKinematics gk = groundKinematics(x);
    const double gamma =
        std::asin(std::clamp(gk.hdot / std::max(1e-3, x.head<3>().norm()), -1.0, 1.0));
    const double sink = -gk.hdot;

    // Wave surface under the keel this step (identically 0 on flat water).
    const double eta_now = wf.eta(x_pos, t);
    const double eta_x = wf.slopeMean(x_pos, t, sc_.waves.contact_len);

    const CtrlVec u = nominal.step(x, ad.V, dt);

    stats_.max_abs_phi = std::max(stats_.max_abs_phi, std::abs(x[PHI]));
    stats_.max_abs_beta = std::max(stats_.max_abs_beta, std::abs(ad.beta));
    stats_.max_alpha = std::max(stats_.max_alpha, ad.alpha);
    stats_.min_V = std::min(stats_.min_V, ad.V);
    stats_.max_abs_y = std::max(stats_.max_abs_y, std::abs(x[Y]));

    csv << t << ',' << x_pos << ',' << x[Y] << ',' << x[H] << ',' << x[U] << ','
        << x[V] << ',' << x[W] << ',' << ad.V << ',' << ad.alpha / kDeg << ','
        << ad.beta / kDeg << ',' << x[P] << ',' << x[Q] << ',' << x[R] << ','
        << x[PHI] / kDeg << ',' << x[THETA] / kDeg << ',' << x[PSI] / kDeg << ','
        << gamma / kDeg << ',' << sink << ',' << u[DE] << ',' << u[DA] << ','
        << u[DR] << ',' << u[DT] << ',' << nominal.thetaCmd() / kDeg << ','
        << nominal.phiCmd() / kDeg << ',' << Wg.u << ',' << Wg.v << ',' << Wg.w
        << ',' << x_gust << ',' << eta_now << ',' << eta_x << '\n';

    if (x[H] <= eta_now && k > 0) {
      td.reached = true;
      td.t = t;
      td.sink = sink;
      td.V = ad.V;
      td.gamma = gamma;
      td.alpha = ad.alpha;
      td.beta = ad.beta;
      td.theta = x[THETA];
      td.phi = x[PHI];
      td.psi = x[PSI];
      td.y = x[Y];
      td.x = x_pos;
      td.h = x[H];
      td.eta = eta_now;
      td.slope = eta_x;
      // Surface-relative closure -d/dt(h - eta) = sink + eta_x xdot + eta_t.
      td.sink_rel = sink + eta_x * gk.xdot_n + wf.etaDot(x_pos, t);
      // TN 1516 slam-load truth, flat- vs wave-referenced (tau / gamma0
      // tilted by the local surface angle, closure onto the moving surface).
      const double m = ac_.inertia.mass, g = ac_.env.g;
      const double alpha_s = std::atan(eta_x);
      td.n_peak_flat = impactNPeakExact(
          m, g, sc_.hull.beta, sc_.hull.rho_water, sc_.hull.eps_g0,
          x[THETA] - sc_.hull.tau_keel, -gamma, std::max(0.0, sink));
      td.n_peak_wave = impactNPeakExact(
          m, g, sc_.hull.beta, sc_.hull.rho_water, sc_.hull.eps_g0,
          x[THETA] - alpha_s - sc_.hull.tau_keel, alpha_s - gamma,
          std::max(0.0, td.sink_rel));
      break;
    }

    const JointState js = rk4WindStep(*dyn_, sc_.wind, t, x, x_gust, u, dt);
    // Ground-track trapezoid across the step (endpoint north velocities).
    x_pos += 0.5 * dt * (gk.xdot_n + groundKinematics(js.x).xdot_n);
    x = js.x;
    x_gust = js.xg;
    t += dt;
    ++stats_.steps;
  }
  stats_.t_end = t;

  std::cout << "=== 6-DOF straight-in landing ===\n";
  std::cout << "trim: theta=" << trim_.theta / kDeg << " deg  de="
            << trim_.u[DE] / kDeg << " deg  dT=" << trim_.u[DT]
            << "  (V_app=" << sc_.V_app << " m/s, gamma_app="
            << sc_.gamma_app / kDeg << " deg)\n";
  if (sc_.wind.enabled) {
    std::cout << "  MIL-F-8785C discrete gust (plant-only, unmeasured): t_start="
              << sc_.wind.t_start << " s  u: " << sc_.wind.u.amp << " m/s/"
              << sc_.wind.u.len << " m (+tailwind)  v: " << sc_.wind.v.amp
              << " m/s/" << sc_.wind.v.len << " m (+east)  w: "
              << sc_.wind.w.amp << " m/s/" << sc_.wind.w.len
              << " m (+updraft)\n";
  }
  if (sc_.waves.enabled) {
    std::cout << "  Surface waves (plant-only): ";
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
    std::cout << "TOUCHDOWN  t=" << td.t << " s  sink=" << td.sink
              << " m/s  V=" << td.V << " m/s  gamma=" << td.gamma / kDeg
              << " deg\n";
    std::cout << "  attitude: theta=" << td.theta / kDeg << " deg  phi="
              << td.phi / kDeg << " deg  psi=" << td.psi / kDeg
              << " deg  beta=" << td.beta / kDeg << " deg  y=" << td.y
              << " m\n";
    if (sc_.waves.enabled) {
      std::cout << "  wave contact: eta=" << td.eta << " m  surface slope="
                << std::atan(td.slope) / kDeg << " deg  sink_rel="
                << td.sink_rel << " m/s (flat-ref " << td.sink << ")\n";
      std::cout << "  TN1516 n_peak truth: flat-ref=" << td.n_peak_flat
                << " g  wave-ref=" << td.n_peak_wave << " g";
      if (td.n_peak_flat > 1e-9)
        std::cout << "  (x" << td.n_peak_wave / td.n_peak_flat << ")";
      std::cout << "\n";
    }
  } else {
    std::cout << "No touchdown within t_max=" << sc_.t_max << " s (final h="
              << x[H] << ")\n";
  }
  std::cout << "  trace -> " << csv_path << "\n";
  return td;
}

}  // namespace autoland
