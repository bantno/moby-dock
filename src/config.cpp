#include "autoland/config.hpp"

#include <yaml-cpp/yaml.h>

#include <cmath>
#include <stdexcept>

namespace autoland {
namespace {

constexpr double kDeg2Rad = M_PI / 180.0;

// Fetch a scalar with a default if the key is missing.
template <typename T>
T get(const YAML::Node& n, const std::string& key, T def) {
  return n[key] ? n[key].as<T>() : def;
}

// Fetch a required scalar; throw if missing.
template <typename T>
T req(const YAML::Node& n, const std::string& key, const std::string& ctx) {
  if (!n[key]) throw std::runtime_error("config: missing required key '" + key +
                                        "' in " + ctx);
  return n[key].as<T>();
}

}  // namespace

AircraftConfig loadAircraftConfig(const std::string& path) {
  YAML::Node root = YAML::LoadFile(path);
  AircraftConfig c;

  const YAML::Node in = root["inertia"];
  if (in) {
    c.inertia.mass = get<double>(in, "mass", 0.0);
    c.inertia.Ixx = get<double>(in, "Ixx", 0.0);
    c.inertia.Iyy = get<double>(in, "Iyy", 0.0);
    c.inertia.Izz = get<double>(in, "Izz", 0.0);
    c.inertia.Ixz = get<double>(in, "Ixz", 0.0);
  }

  const YAML::Node cg = root["cg_offset"];
  if (cg) {
    c.cg.dx = get<double>(cg, "dx", 0.0);
    c.cg.dy = get<double>(cg, "dy", 0.0);
    c.cg.dz = get<double>(cg, "dz", 0.0);
  }

  const YAML::Node th = root["thrust"];
  if (th) {
    c.thrust.T_static = get<double>(th, "T_static", 0.0);
    c.thrust.k_v = get<double>(th, "k_v", 0.0);
    c.thrust.zcp = get<double>(th, "zcp", 0.0);
  }

  const YAML::Node aero = root["aero"];
  if (aero) c.parasite_CD0 = get<double>(aero, "parasite_CD0", 0.0);

  const YAML::Node env = root["environment"];
  if (env) {
    c.env.rho = get<double>(env, "rho", 1.225);
    c.env.a_sound = get<double>(env, "a_sound", 340.3);
    c.env.g = get<double>(env, "g", 9.80665);
  }

  const YAML::Node lim = root["limits"];
  if (lim) {
    // limits given in degrees in YAML for readability; convert to radians.
    c.limits.de_min = get<double>(lim, "de_min_deg", -30.0) * kDeg2Rad;
    c.limits.de_max = get<double>(lim, "de_max_deg", 30.0) * kDeg2Rad;
    c.limits.da_min = get<double>(lim, "da_min_deg", -30.0) * kDeg2Rad;
    c.limits.da_max = get<double>(lim, "da_max_deg", 30.0) * kDeg2Rad;
    c.limits.dr_min = get<double>(lim, "dr_min_deg", -30.0) * kDeg2Rad;
    c.limits.dr_max = get<double>(lim, "dr_max_deg", 30.0) * kDeg2Rad;
    c.limits.rate = get<double>(lim, "rate_deg_s", 120.0) * kDeg2Rad;
    c.limits.dT_min = get<double>(lim, "dT_min", 0.0);
    c.limits.dT_max = get<double>(lim, "dT_max", 1.0);
  }

  // Mixing matrix: rows = control groups, columns = [delta_e, delta_a, delta_r].
  const YAML::Node mx = root["mixing"];
  if (mx && mx["matrix"]) {
    const YAML::Node rows = mx["matrix"];
    const int nrows = static_cast<int>(rows.size());
    c.mixing = Mat::Zero(nrows, 3);
    for (int i = 0; i < nrows; ++i) {
      const YAML::Node row = rows[i];
      if (row.size() != 3)
        throw std::runtime_error("config: each mixing matrix row needs 3 "
                                 "entries [de, da, dr]");
      for (int j = 0; j < 3; ++j) c.mixing(i, j) = row[j].as<double>();
    }
    if (mx["group_order"])
      c.mixing_group_order = mx["group_order"].as<std::vector<std::string>>();
  }
  // If absent, mixing stays empty and Mixing::build() constructs the default.

  return c;
}

ScenarioConfig loadScenarioConfig(const std::string& path) {
  YAML::Node root = YAML::LoadFile(path);
  ScenarioConfig s;

  s.V_app = req<double>(root, "V_app", "scenario");
  s.gamma_app = req<double>(root, "gamma_app_deg", "scenario") * kDeg2Rad;
  s.h_decrab = get<double>(root, "h_decrab", 0.0);
  s.dt = get<double>(root, "dt", 0.01);
  s.t_max = get<double>(root, "t_max", 120.0);

  const YAML::Node g = root["gains"];
  if (g) {
    s.gains.Kp_V = get<double>(g, "Kp_V", 0.0);
    s.gains.Ki_V = get<double>(g, "Ki_V", 0.0);
    s.gains.Kp_gs = get<double>(g, "Kp_gs", 0.0);
    s.gains.Kd_gs = get<double>(g, "Kd_gs", 0.0);
    s.gains.Kp_theta = get<double>(g, "Kp_theta", 0.0);
    s.gains.Ki_theta = get<double>(g, "Ki_theta", 0.0);
    s.gains.Kq = get<double>(g, "Kq", 0.0);
    s.gains.Kp_y = get<double>(g, "Kp_y", 0.0);
    s.gains.Kd_y = get<double>(g, "Kd_y", 0.0);
    s.gains.Ki_y = get<double>(g, "Ki_y", 0.0);
    s.gains.phi_max = get<double>(g, "phi_max_deg", 30.0) * kDeg2Rad;
    s.gains.Kp_phi = get<double>(g, "Kp_phi", 0.0);
    s.gains.Kp_p = get<double>(g, "Kp_p", 0.0);
    s.gains.Kr = get<double>(g, "Kr", 0.0);
    s.gains.Kp_beta = get<double>(g, "Kp_beta", 0.0);
  }

  const YAML::Node fl = root["flare"];
  if (fl) {
    s.flare.h_flare = get<double>(fl, "h_flare", 0.0);
    s.flare.tau = get<double>(fl, "tau", 1.0);
    s.flare.w_td = get<double>(fl, "w_td", 0.5);
  }

  const YAML::Node ic = root["initial"];
  if (ic) {
    s.init.distance = get<double>(ic, "distance", 300.0);
    s.init.dh = get<double>(ic, "dh", 0.0);
    s.init.dy = get<double>(ic, "dy", 0.0);
    s.init.dpsi = get<double>(ic, "dpsi_deg", 0.0) * kDeg2Rad;
    s.init.dV = get<double>(ic, "dV", 0.0);
  }

  // --- Nominal controller selection -----------------------------------------
  const std::string nom = get<std::string>(root, "nominal", "cascaded_pid");
  if (nom == "lon_glideslope")
    s.nominal = NominalKind::LonGlideslope;
  else if (nom == "cascaded_pid")
    s.nominal = NominalKind::CascadedPID;
  else
    throw std::runtime_error("config: unknown 'nominal' value '" + nom +
                             "' (expected cascaded_pid | lon_glideslope)");

  // --- CBF safety-filter parameters -----------------------------------------
  const YAML::Node cbf = root["cbf"];
  if (cbf) {
    s.cbf.enabled = get<bool>(cbf, "enabled", true);
    s.cbf.airspeed_barrier = get<bool>(cbf, "airspeed_barrier", true);
    s.cbf.descent_barrier = get<bool>(cbf, "descent_barrier", true);
    s.cbf.airspeed_hard = get<bool>(cbf, "airspeed_hard", false);
    s.cbf.descent_hard = get<bool>(cbf, "descent_hard", true);
    s.cbf.V_min = get<double>(cbf, "V_min", 0.0);
    s.cbf.alpha_V = get<double>(cbf, "alpha_V", 1.0);
    s.cbf.v_safe = get<double>(cbf, "v_safe", 0.6);
    s.cbf.a_brk = get<double>(cbf, "a_brk", 3.0);
    s.cbf.alpha_h = get<double>(cbf, "alpha_h", 1.0);
    s.cbf.w_de = get<double>(cbf, "w_de", 1.0);
    s.cbf.w_dT = get<double>(cbf, "w_dT", 1.0);
    s.cbf.slack_penalty = get<double>(cbf, "slack_penalty", 1.0e4);
  }

  return s;
}

}  // namespace autoland
