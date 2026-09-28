#include "autoland/sixdof_sim.hpp"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
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
// so this single path is bit-identical to a still-air RK4 of the plant xdot.
struct JointState {
  StateVec x;
  double xg;
};
JointState rk4WindStep(const SixDofSim::PlantFn& plant,
                       const DiscreteGustConfig& g, double t, const StateVec& x,
                       double xg, const CtrlVec& u, double dt) {
  auto deriv = [&](double tl, const StateVec& xl, double xgl, StateVec& xd,
                   double& xgd) {
    const GustWind Wg = gustWind(g, xgl);
    xgd = gustXdot(g, tl, airData(xl, Wg).V);
    xd = plant(xl, u, Eigen::Vector3d(Wg.u, Wg.v, Wg.w));
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

// Exact airspeed rate for the nominal (stands in for the accelerometer-fed
// airspeed filter of PX4's TECS): V_a = |v_b - W_b| with W_b the gust in body
// axes, so
//   Vdot_a = v_a . (vdot_b - Wdot_b) / V_a,   Wdot_b = R_nb^T Wdot_e - omega x W_b
// (body-resolved rate of an earth-frame vector). vdot_b is the plant xdot
// under the control held over the LAST step (zero-order hold), i.e. the
// acceleration the aircraft is actually undergoing when the new command is
// computed -- no algebraic loop with that command. Wind-off this is v.vdot/|v|.
double airspeedRate(const StateVec& x, const StateVec& xd, const GustWind& Wg,
                    const GustWind& Wgdot) {
  const Eigen::Vector3d Wb = windBody(x, Wg);
  const Eigen::Vector3d omega(x[P], x[Q], x[R]);
  const Eigen::Vector3d Wbdot = windBody(x, Wgdot) - omega.cross(Wb);  // linear in W
  const Eigen::Vector3d va(x[U] - Wb[0], x[V] - Wb[1], x[W] - Wb[2]);
  const Eigen::Vector3d vadot(xd[U] - Wbdot[0], xd[V] - Wbdot[1], xd[W] - Wbdot[2]);
  const double Va = std::max(1e-3, va.norm());
  return va.dot(vadot) / Va;
}

// Steady wings-level trim at (V, gamma) on whichever plant is active.
using TrimFn = std::function<TrimResult(double V, double gamma)>;

// Steady-flight gamma at which the trim throttle equals dT_target at airspeed
// V (NaN if a trim fails or no crossing exists within 40 deg of gamma0). The
// trim throttle is monotone in gamma, so: march 1 deg from gamma0 until the
// throttle residual changes sign, then bisect. Pure trim evaluations -- no
// derivatives.
double gammaForThrottle(const TrimFn& trimAt, double V, double dT_target,
                        double gamma0) {
  auto residual = [&](double gam, bool& ok) {
    const TrimResult r = trimAt(V, gam);
    ok = r.converged;
    return r.u[DT] - dT_target;
  };
  bool ok = false;
  double g0 = gamma0, f0 = residual(g0, ok);
  if (!ok) return NAN;
  const double step = (f0 < 0.0) ? kDeg : -kDeg;  // short of throttle -> steeper climb
  double g1 = g0, f1 = f0;
  bool bracketed = false;
  for (int i = 0; i < 40 && !bracketed; ++i) {
    g1 = g0 + step;
    f1 = residual(g1, ok);
    if (!ok) return NAN;
    if (f0 * f1 <= 0.0) bracketed = true;
    else { g0 = g1; f0 = f1; }
  }
  if (!bracketed) return NAN;
  for (int i = 0; i < 30; ++i) {
    const double gm = 0.5 * (g0 + g1);
    const double fm = residual(gm, ok);
    if (!ok) return NAN;
    if (f0 * fm <= 0.0) { g1 = gm; f1 = fm; } else { g0 = gm; f0 = fm; }
  }
  return 0.5 * (g0 + g1);
}

// Fill the TECS vehicle anchors from steady-flight trims on the active plant.
// PX4 takes these as parameters (FW_THR_TRIM, FW_PSP_OFF, FW_T_CLMB_MAX,
// FW_T_SINK_MIN, FW_T_SINK_MAX) scaled by a performance model; here they are
// the plant's own equilibria, so the throttle feedforward's three anchors
// (dT_min -> -min_sink_rate, throttle_trim -> level, dT_max -> max_climb_rate)
// are exact and the pitch offset is the level-flight theta at V_ref:
//   throttle_trim, pitch_offset : level trim at V_ref
//   max_climb_rate              : steady climb at dT_max, V_ref
//   min_sink_rate               : steady sink at dT_min, V_ref
//   max_sink_rate               : steady sink at dT_min, tas_max (demand limit)
// Any anchor that cannot be solved keeps its PX4 default, with a warning.
void deriveTecsAnchors(const TrimFn& trimAt, SixDofNominalConfig& nom,
                       double gamma_app) {
  px4::TecsParam& tp = nom.tecs;
  const double V = nom.V_ref;
  const SurfaceLimits& L = nom.limits;
  auto warn = [](const std::string& what) {
    std::cerr << "[sixdof_sim] warning: TECS anchor '" << what
              << "' not derivable from the plant; keeping the default\n";
  };

  const TrimResult level = trimAt(V, 0.0);
  if (level.converged) {
    tp.throttle_trim = std::clamp(level.u[DT], L.dT_min, L.dT_max);
    nom.tecs_pitch_offset = level.theta;
  } else {
    warn("throttle_trim / pitch_offset (level trim)");
    nom.tecs_pitch_offset = nom.theta_trim;
  }

  const double g_climb = gammaForThrottle(trimAt, V, L.dT_max, gamma_app);
  if (std::isfinite(g_climb) && V * std::sin(g_climb) > 0.0)
    tp.max_climb_rate = V * std::sin(g_climb);
  else
    warn("max_climb_rate (full-throttle climb)");

  const double g_sink = gammaForThrottle(trimAt, V, L.dT_min, gamma_app);
  if (std::isfinite(g_sink) && -V * std::sin(g_sink) > 0.0)
    tp.min_sink_rate = -V * std::sin(g_sink);
  else
    warn("min_sink_rate (idle sink at V_ref)");

  const double g_sink_fast = gammaForThrottle(trimAt, tp.tas_max, L.dT_min, gamma_app);
  if (std::isfinite(g_sink_fast) && -tp.tas_max * std::sin(g_sink_fast) > 0.0)
    tp.max_sink_rate = -tp.tas_max * std::sin(g_sink_fast);
  else
    warn("max_sink_rate (idle sink at tas_max)");
}
}  // namespace

SixDofSim::SixDofSim(const std::string& stab_path,
                     const std::string& aircraft_yaml,
                     const std::string& scenario_yaml) {
  YAML::Node root = YAML::LoadFile(scenario_yaml);
  const std::string plant =
      (root && root["plant"]) ? root["plant"].as<std::string>() : "beaver";
  const bool beaver = (plant == "beaver");
  if (!beaver && plant != "vspaero")
    throw std::runtime_error("sixdof_sim: unknown plant '" + plant + "'");

  SixDofNominalConfig& nom = sc_.nominal;
  if (beaver) {
    // --- DHC-2 Beaver plant (validated LR-556/FDC polynomials). -------------
    BeaverPlantConfig bpc;
    YAML::Node yb = root["beaver"];
    bpc.n_rpm = getOr(yb, "n_rpm", bpc.n_rpm);
    bpc.pz_idle = getOr(yb, "pz_idle", bpc.pz_idle);
    bpc.pz_max = getOr(yb, "pz_max", bpc.pz_max);
    bpc.flap = getOr(yb, "flap_deg", 0.0) * kDeg;
    bpc.h_ref = getOr(yb, "h_ref", 0.0);
    bdyn_ = std::make_unique<BeaverDynamics>(bpc);
    plant_ = [this](const StateVec& x, const CtrlVec& u,
                    const Eigen::Vector3d& W) { return bdyn_->xdot(x, u, W); };

    // AircraftConfig mirror for the sim's diagnostics (TN 1516 mass/g) and
    // the nominal's surface limits.
    ac_.inertia.mass = BeaverGeom::mass;
    ac_.inertia.Ixx = BeaverGeom::Ix;
    ac_.inertia.Iyy = BeaverGeom::Iy;
    ac_.inertia.Izz = BeaverGeom::Iz;
    ac_.inertia.Ixz = BeaverGeom::Ixz;
    ac_.env.rho = bdyn_->rho();
    ac_.env.g = bdyn_->g();
    ac_.limits = bpc.limits;

    // Beaver-scale defaults (2288 kg): approach inside the LR-556 30-55 m/s
    // validity band; gains sized to the Beaver's control authorities
    // (M_de ~ 7.7 rad/s^2 per rad at 35 m/s, L_da ~ 4.7, throttle ~ 1.5
    // m/s^2 per unit dT) -- the AHAB defaults are for a 3.6 kg airframe.
    sc_.V_app = 35.0;
    sc_.gamma_app = -3.5 * kDeg;
    nom.Kp_V = 0.15; nom.Ki_V = 0.03;
    nom.Kp_gamma = 1.5; nom.Ki_gamma = 0.3;
    nom.Kv_gamma = 0.02; nom.dgamma_V_max = 3.0 * kDeg;
    nom.theta_cmd_max = 10.0 * kDeg;
    nom.Kp_theta = 0.8; nom.Ki_theta = 0.3; nom.Kq = 0.5;
    nom.Kp_y = 0.01; nom.Kd_y = 0.05; nom.phi_max = 15.0 * kDeg;
    nom.Kp_phi = 2.0; nom.Kp_p = 0.8; nom.Kr = 0.5;
    // Standard Delft/FDC deflection signs (Cm_de, Cl_da, Cn_dr all < 0).
    nom.de_sign = -1.0; nom.da_sign = -1.0; nom.dr_sign = -1.0;
    // TECS airspeed demand band = the LR-556 validity band (underspeed
    // mitigation ramps in below tas_min - 0.15 V_ref; px4_tecs.hpp).
    nom.tecs.tas_min = 30.0; nom.tecs.tas_max = 55.0;
    // TECS gains tuned for the Beaver (scripts/tune_tecs.py; protocol and
    // tables in documentation/px4_tecs_port.md). Three of PX4's flown
    // defaults change, the rest stay upstream:
    //   FW_T_PTCH_DAMP  0.1 -> 1.0   pitch damping is exactly a pitch-per-
    //                                gamma gain (d theta = damp * d hdot / V);
    //                                0.1 leaves the SEB loop integrator-
    //                                dominated and ringing (~15 s). 1.0-1.5 is
    //                                a flat optimum; the cascade's hand-tuned
    //                                Kp_gamma is 1.5.
    //   FW_T_I_GAIN_PIT 0.1 -> 0.4   pitch-rate-per-gamma integral gain;
    //                                optimum 0.4 (degrades above 0.8); the
    //                                cascade's Ki_gamma is 0.3.
    //   FW_T_THR_INTEG  0.02 -> 0.3  at 0.02 the linear-about-trim throttle
    //                                map's ~0.05 under-prediction of the
    //                                approach throttle takes ~100 s to trim
    //                                out (5% steep touchdown); flat 0.3..1.0.
    // Picked from the flat region of a 2-pass grid on four longitudinal
    // training cases, then confirmed on seven held-out cases (POH flaps-35
    // config, gusts, crosswind, long/steep/fast approaches) -- every case
    // improves vs the PX4 defaults. FW_T_THR_DAMPING (0.05) is insensitive
    // (0.2..1.0 gains 1% on validation) and stays at the PX4 default.
    nom.tecs.pitch_damping_gain = 1.0;
    nom.tecs.integrator_gain_pitch = 0.4;
    nom.tecs.integrator_gain_throttle = 0.3;
  } else {
    // --- AHAB VSPAERO-table plant (original path). --------------------------
    table_ = std::make_unique<AeroTable>(AeroTable::fromFile(stab_path));
    ac_ = loadAircraftConfig(aircraft_yaml);
    mixing_ = std::make_unique<Mixing>(Mixing::build(ac_, *table_));
    dyn_ = std::make_unique<Dynamics>(*table_, *mixing_, ac_);
    plant_ = [this](const StateVec& x, const CtrlVec& u,
                    const Eigen::Vector3d& W) { return dyn_->xdot(x, u, W); };
    sc_.V_app = 18.0;
    sc_.gamma_app = -3.0 * kDeg;
    nom.tecs.tas_min = 12.0; nom.tecs.tas_max = 25.0;  // placeholder band
  }

  sc_.V_app = getOr(root, "V_app", sc_.V_app);
  sc_.gamma_app = getOr(root, "gamma_app_deg", sc_.gamma_app / kDeg) * kDeg;
  sc_.dt = getOr(root, "dt", 0.01);
  sc_.t_max = getOr(root, "t_max", 200.0);

  // --- Reach-the-target energy CBF layer (`energy_cbf:` block). -------------
  // Parsed BEFORE the approach trim: in idle-throttle mode the approach path
  // is the plant's own idle glide at V_app (gamma_app is REPLACED by the
  // steady idle-glide flight-path angle), so the trim, the nominal gamma_ref,
  // and the CBF's unpowered drift model all describe the same flight.
  {
    YAML::Node yr = root["energy_cbf"];
    ReachSimCfg& rc = sc_.reach;
    rc.enabled = getOrB(yr, "enabled", false);
    if (rc.enabled) {
      if (!beaver)
        throw std::runtime_error(
            "sixdof_sim: energy_cbf requires the beaver plant (the glide "
            "polar is fitted from the Beaver polynomials)");
      rc.filter = getOrB(yr, "filter", true);
      if (yr && yr["throttle"]) {
        const std::string thr = yr["throttle"].as<std::string>();
        if (thr != "idle" && thr != "nominal")
          throw std::runtime_error("sixdof_sim: energy_cbf.throttle must be "
                                   "'idle' or 'nominal'");
        rc.throttle_idle = (thr == "idle");
      }
      rc.tau_alpha = getOr(yr, "tau_alpha", rc.tau_alpha);
      rc.target_N = getOr(yr, "target_N", rc.target_N);
      rc.target_E = getOr(yr, "target_E", rc.target_E);
      EnergyReachCfg& ec = rc.cbf;
      ec.enabled = rc.filter;
      ec.ck = getOr(yr, "ck", ec.ck);
      ec.c_alpha = getOr(yr, "c_alpha", ec.c_alpha);
      ec.ua_max = getOr(yr, "ua_max_deg", ec.ua_max / kDeg) * kDeg;
      ec.beta_max = getOr(yr, "beta_max_deg", ec.beta_max / kDeg) * kDeg;
      ec.alpha_cmd_min =
          getOr(yr, "alpha_cmd_min_deg", ec.alpha_cmd_min / kDeg) * kDeg;
      ec.alpha_cmd_max =
          getOr(yr, "alpha_cmd_max_deg", ec.alpha_cmd_max / kDeg) * kDeg;
      ec.prm.V_min = getOr(yr, "V_min", ec.prm.V_min);
      ec.prm.phi_sat = getOr(yr, "phi_sat", ec.prm.phi_sat);

      // Physical constants + the idle-glide polar from the Beaver plant. The
      // fit is at the IDLE throttle so C_D0 absorbs the windmilling-propeller
      // drag -- under nominal power the barrier is then conservative.
      ReachParams& pp = ec.prm;
      pp.m = BeaverGeom::mass;
      pp.S = BeaverGeom::S;
      pp.rho = bdyn_->rho();
      pp.g = bdyn_->g();
      pp.CYb = bdyn_->config().aero.Cy_b;
      polar_ = fitBeaverReachPolar(*bdyn_, sc_.V_app, ac_.limits.dT_min);
      pp.CLa = polar_.CLa;
      pp.CD0 = polar_.CD0;
      pp.CDa2 = polar_.CDa2;
      pp.alpha0L = polar_.alpha0L;

      if (rc.throttle_idle) {
        const TrimFn trimAt = [this](double V, double gam) {
          return beaverTrim(*bdyn_, V, gam);  // beaver guaranteed above
        };
        const double g_idle = gammaForThrottle(trimAt, sc_.V_app,
                                               ac_.limits.dT_min, sc_.gamma_app);
        if (std::isfinite(g_idle))
          sc_.gamma_app = g_idle;
        else
          std::cerr << "[sixdof_sim] warning: idle-glide gamma solve failed; "
                       "keeping gamma_app = " << sc_.gamma_app / kDeg
                    << " deg\n";
      }
    }
  }

  // Approach trim: nominal feedforwards + the initial condition. The Beaver
  // trims all six axes (nonzero beta/da/dr from slipstream asymmetry).
  trim_ = beaver ? beaverTrim(*bdyn_, sc_.V_app, sc_.gamma_app)
                 : trim(*dyn_, sc_.V_app, sc_.gamma_app);
  stats_.trim_converged = trim_.converged;
  if (!trim_.converged)
    std::cerr << "[sixdof_sim] warning: approach trim did not converge (res="
              << trim_.residual << ")\n";

  // --- Nominal defaults from trim (then YAML overrides) ---
  nom.V_ref = sc_.V_app;
  nom.gamma_ref = sc_.gamma_app;
  nom.theta_trim = trim_.theta;
  nom.de_trim = trim_.u[DE];
  nom.dT_trim = trim_.u[DT];
  nom.da_trim = trim_.u[DA];
  nom.dr_trim = trim_.u[DR];
  nom.limits = ac_.limits;
  YAML::Node yn = root["nominal"];
  nom.V_ref = getOr(yn, "V_ref", nom.V_ref);
  // Throttle authority overrides (engine-out / fixed-power cases): the
  // nominal's throttle loop clamps to these; the plant is unchanged.
  nom.limits.dT_min = getOr(yn, "dT_min", nom.limits.dT_min);
  nom.limits.dT_max = getOr(yn, "dT_max", nom.limits.dT_max);
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
  nom.Kp_beta = getOr(yn, "Kp_beta", nom.Kp_beta);
  nom.Ki_beta = getOr(yn, "Ki_beta", nom.Ki_beta);

  // --- Longitudinal outer loop: cascaded PID (default) or the PX4 TECS port.
  if (yn && yn["type"]) {
    const std::string type = yn["type"].as<std::string>();
    if (type == "tecs") nom.lon_mode = LonMode::Tecs;
    else if (type != "cascade")
      throw std::runtime_error("sixdof_sim: unknown nominal.type '" + type + "'");
  }
  if (nom.lon_mode == LonMode::Tecs) {
    px4::TecsParam& tp = nom.tecs;
    YAML::Node yt = root["tecs"];
    // Vehicle limits from the nominal/aircraft config (PX4: FW_AIRSPD_TRIM,
    // FW_THR_MIN/MAX, FW_P_LIM_MIN/MAX - FW_PSP_OFF, FW_AIRSPD_MIN/MAX).
    tp.equivalent_airspeed_trim = nom.V_ref;
    tp.throttle_min = nom.limits.dT_min;
    tp.throttle_max = nom.limits.dT_max;
    tp.pitch_max = nom.theta_cmd_max;
    tp.pitch_min = -nom.theta_cmd_max;
    tp.tas_min = getOr(yt, "tas_min", tp.tas_min);
    tp.tas_max = getOr(yt, "tas_max", tp.tas_max);
    tp.pitch_max = getOr(yt, "pitch_max_deg", tp.pitch_max / kDeg) * kDeg;
    tp.pitch_min = getOr(yt, "pitch_min_deg", tp.pitch_min / kDeg) * kDeg;
    // Plant-derived anchors (level trim, full-throttle climb, idle sink).
    const TrimFn trimAt = [this, beaver](double V, double gam) {
      return beaver ? beaverTrim(*bdyn_, V, gam) : trim(*dyn_, V, gam);
    };
    deriveTecsAnchors(trimAt, nom, sc_.gamma_app);
    // YAML overrides. Keys are the PX4 parameter names, lower-case, without
    // the FW_T_ / FW_ prefix; time constants as in PX4 (gain = 1 / max(tc, 0.1)).
    tp.throttle_trim = getOr(yt, "thr_trim", tp.throttle_trim);
    nom.tecs_pitch_offset =
        getOr(yt, "pitch_offset_deg", nom.tecs_pitch_offset / kDeg) * kDeg;
    tp.max_climb_rate = getOr(yt, "clmb_max", tp.max_climb_rate);
    tp.min_sink_rate = getOr(yt, "sink_min", tp.min_sink_rate);
    tp.max_sink_rate = getOr(yt, "sink_max", tp.max_sink_rate);
    tp.vert_accel_limit = getOr(yt, "vert_acc", tp.vert_accel_limit);
    if (yt && yt["alt_tc"])
      tp.altitude_error_gain = 1.0 / std::max(yt["alt_tc"].as<double>(), 0.1);
    tp.altitude_setpoint_gain_ff = getOr(yt, "hrate_ff", tp.altitude_setpoint_gain_ff);
    if (yt && yt["tas_tc"])
      tp.airspeed_error_gain = 1.0 / std::max(yt["tas_tc"].as<double>(), 0.1);
    tp.tas_error_percentage = getOr(yt, "tas_error_percentage", tp.tas_error_percentage);
    tp.ste_rate_time_const = getOr(yt, "ste_r_tc", tp.ste_rate_time_const);
    tp.seb_rate_ff = getOr(yt, "seb_r_ff", tp.seb_rate_ff);
    tp.pitch_speed_weight = getOr(yt, "spdweight", tp.pitch_speed_weight);
    tp.integrator_gain_pitch = getOr(yt, "i_gain_pit", tp.integrator_gain_pitch);
    tp.pitch_damping_gain = getOr(yt, "ptch_damp", tp.pitch_damping_gain);
    tp.integrator_gain_throttle = getOr(yt, "thr_integ", tp.integrator_gain_throttle);
    tp.throttle_damping_gain = getOr(yt, "thr_damping", tp.throttle_damping_gain);
    tp.throttle_slewrate = getOr(yt, "thr_slew_max", tp.throttle_slewrate);
    tp.load_factor_correction = getOr(yt, "rll2thr", tp.load_factor_correction);
    nom.tecs_detect_underspeed = getOrB(yt, "detect_underspeed", true);
  }

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
  sc_.hull.n_surfaces = static_cast<int>(
      getOr(yh, "n_surfaces", static_cast<double>(sc_.hull.n_surfaces)));

  // --- Initial condition: approach trim + scenario offsets. ---
  YAML::Node yi = root["initial"];
  sc_.h0 = getOr(yi, "h0", 40.0);
  sc_.dV = getOr(yi, "dV", 0.0);
  sc_.dy = getOr(yi, "dy", 0.0);
  sc_.dpsi = getOr(yi, "dpsi_deg", 0.0) * kDeg;
  sc_.dtheta = getOr(yi, "dtheta_deg", 0.0) * kDeg;
  sc_.N0 = getOr(yi, "N0", 0.0);
  sc_.E0 = getOr(yi, "E0", 0.0);
  sc_.psi0 = getOr(yi, "psi_deg", 0.0) * kDeg;

  // --- Corridor landing world (`world:` block): corridor, terrain, rollout. -
  {
    YAML::Node yw = root["world"];
    WorldConfig& w = sc_.world;
    w.enabled = getOrB(yw, "enabled", false);
    if (w.enabled) {
      if (sc_.reach.enabled)
        throw std::runtime_error(
            "sixdof_sim: world.enabled and energy_cbf.enabled are mutually "
            "exclusive (the reach layer forces wings-level skid-to-turn)");
      YAML::Node yc = yw["corridor"];
      w.corridor.N_c = getOr(yc, "center_N", 0.0);
      w.corridor.E_c = getOr(yc, "center_E", 0.0);
      w.corridor.length = getOr(yc, "length", 700.0);
      w.corridor.width = getOr(yc, "width", 80.0);
      w.corridor.heading = getOr(yc, "heading_deg", 0.0) * kDeg;
      w.corridor.s_aim = getOr(yc, "aim_s", -0.5 * w.corridor.length + 100.0);
      YAML::Node yt = yw["terrain"];
      w.terrain.base = getOr(yt, "base", 0.0);
      if (yt && yt["bumps"]) {
        for (const YAML::Node& yb : yt["bumps"]) {
          TerrainBump b;
          b.name = (yb["name"]) ? yb["name"].as<std::string>() : "bump";
          b.N0 = getOr(yb, "N", 0.0);
          b.E0 = getOr(yb, "E", 0.0);
          b.height = getOr(yb, "height", 0.0);
          b.sigma_a = getOr(yb, "sigma_a", 10.0);
          b.sigma_b = getOr(yb, "sigma_b", b.sigma_a);
          b.rot = getOr(yb, "rot_deg", 0.0) * kDeg;
          b.order = static_cast<int>(getOr(yb, "order", 1.0));
          w.terrain.bumps.push_back(b);
        }
      }
      if (yt && yt["rings"]) {
        for (const YAML::Node& yr : yt["rings"]) {
          const std::string name =
              (yr["name"]) ? yr["name"].as<std::string>() : "ring";
          w.terrain.addRing(getOr(yr, "center_N", 0.0), getOr(yr, "center_E", 0.0),
                            getOr(yr, "radius", 500.0),
                            static_cast<int>(getOr(yr, "count", 32.0)),
                            getOr(yr, "height", 15.0), getOr(yr, "sigma", 25.0),
                            static_cast<int>(getOr(yr, "order", 2.0)), name,
                            getOr(yr, "gap_heading_deg", NAN) * kDeg,
                            getOr(yr, "gap_halfangle_deg", 0.0) * kDeg);
        }
      }
      YAML::Node yro = yw["rollout"];
      w.rollout.enabled = getOrB(yro, "enabled", true);
      w.rollout.a0 = getOr(yro, "a0", w.rollout.a0);
      w.rollout.kq = getOr(yro, "kq", w.rollout.kq);
      w.rollout.k_psi = getOr(yro, "k_psi", w.rollout.k_psi);
      w.rollout.V_stop = getOr(yro, "V_stop", w.rollout.V_stop);
      w.rollout.t_max = getOr(yro, "t_max", w.rollout.t_max);

      // Pattern guidance (`guidance:` block): the approach speed / slope are
      // the scenario's; everything else has header defaults.
      YAML::Node yg = root["guidance"];
      PatternGuidanceConfig& g = sc_.guidance;
      g.V_app = sc_.V_app;
      g.gamma_app = sc_.gamma_app;
      g.g = ac_.env.g;
      sc_.nominal.g = ac_.env.g;
      g.phi_plan = getOr(yg, "phi_plan_deg", g.phi_plan / kDeg) * kDeg;
      g.k_R = getOr(yg, "k_R", g.k_R);
      g.d_final_min = getOr(yg, "d_final_min", g.d_final_min);
      g.d_final_max = getOr(yg, "d_final_max", g.d_final_max);
      g.axis_overrun = getOr(yg, "axis_overrun", g.axis_overrun);
      g.K_h = getOr(yg, "K_h", g.K_h);
      g.gamma_min = getOr(yg, "gamma_min_deg", g.gamma_min / kDeg) * kDeg;
      g.gamma_max = getOr(yg, "gamma_max_deg", g.gamma_max / kDeg) * kDeg;
      g.L1_period = getOr(yg, "L1_period", g.L1_period);
      g.L1_scale = getOr(yg, "L1_scale", g.L1_scale);
      g.direct_e_tol = getOr(yg, "direct_e_tol", g.direct_e_tol);
      g.direct_chi_tol = getOr(yg, "direct_chi_tol_deg", g.direct_chi_tol / kDeg) * kDeg;
      if (yi && (yi["dy"] || yi["dpsi_deg"]))
        std::cerr << "[sixdof_sim] warning: initial.dy / initial.dpsi_deg are "
                     "ignored when world.enabled (use N0 / E0 / psi_deg)\n";
      if (!(sc_.gamma_app < 0.0))
        throw std::runtime_error(
            "sixdof_sim: world.enabled needs a descending gamma_app_deg");
      g.L1_min = getOr(yg, "L1_min", g.L1_min);
      g.phi_max = getOr(yg, "phi_max_deg", sc_.nominal.phi_max / kDeg) * kDeg;
    }
  }

  // The energy-CBF layer commands alpha inside its trusted band; if the
  // approach TRIM alpha already sits outside it the band rows will fight the
  // glide from the first step -- surface that loudly.
  if (sc_.reach.enabled && trim_.converged) {
    const double a_tr = std::atan2(trim_.x[W], trim_.x[U]);
    if (a_tr > sc_.reach.cbf.alpha_cmd_max || a_tr < sc_.reach.cbf.alpha_cmd_min)
      std::cerr << "[sixdof_sim] warning: trim alpha " << a_tr / kDeg
                << " deg is outside the energy_cbf AoA band ["
                << sc_.reach.cbf.alpha_cmd_min / kDeg << ", "
                << sc_.reach.cbf.alpha_cmd_max / kDeg << "] deg\n";
  }

  // Energy-CBF target: explicit YAML point, else the glideslope aim point --
  // the ground intercept of the approach path from the initial condition
  // (north position starts at 0, the target altitude datum is the water at 0).
  if (sc_.reach.enabled) {
    ReachParams& pp = sc_.reach.cbf.prm;
    const double gs = std::max(1e-3, std::tan(-sc_.gamma_app));
    pp.Nt = std::isfinite(sc_.reach.target_N) ? sc_.reach.target_N
                                              : sc_.h0 / gs;
    pp.Et = std::isfinite(sc_.reach.target_E) ? sc_.reach.target_E : 0.0;
  }

  // --- 6-DOF surfaces-only CBF filter (`sixdof_cbf:` block). ---------------
  {
    YAML::Node yc = root["sixdof_cbf"];
    SixDofCbfConfig& c = sc_.cbf;
    c.enabled = getOrB(yc, "enabled", false);
    if (c.enabled && !beaver)
      throw std::runtime_error("sixdof_sim: sixdof_cbf needs plant: beaver");
    c.filter = getOrB(yc, "filter", true);
    c.row_alpha = getOrB(yc, "row_alpha", true);
    c.row_V = getOrB(yc, "row_V", true);
    c.row_phi = getOrB(yc, "row_phi", true);
    c.row_beta = getOrB(yc, "row_beta", true);
    c.alpha_max = getOr(yc, "alpha_max_deg", c.alpha_max / kDeg) * kDeg;
    c.V_min = getOr(yc, "V_min", c.V_min);
    c.V_max = getOr(yc, "V_max", c.V_max);
    c.phi_max = getOr(yc, "phi_max_deg", c.phi_max / kDeg) * kDeg;
    c.beta_max = getOr(yc, "beta_max_deg", c.beta_max / kDeg) * kDeg;
    auto gains = [&](const char* key, std::array<double, 2>& g) {
      if (yc && yc[key]) { g[0] = yc[key][0].as<double>(); g[1] = yc[key][1].as<double>(); }
    };
    gains("c_alpha", c.c_alpha); gains("c_V", c.c_V); gains("c_phi", c.c_phi); gains("c_beta", c.c_beta);
    c.hard_alpha = getOrB(yc, "hard_alpha", false);
    c.hard_V = getOrB(yc, "hard_V", false);
    c.hard_phi = getOrB(yc, "hard_phi", false);
    c.hard_beta = getOrB(yc, "hard_beta", false);
    c.w_alpha = getOr(yc, "w_alpha", c.w_alpha);
    c.w_V = getOr(yc, "w_V", c.w_V);
    c.w_phi = getOr(yc, "w_phi", c.w_phi);
    c.w_beta = getOr(yc, "w_beta", c.w_beta);
    c.w_de = getOr(yc, "w_de", c.w_de);
    c.w_da = getOr(yc, "w_da", c.w_da);
    c.w_dr = getOr(yc, "w_dr", c.w_dr);
    c.margin_force_terms = getOrB(yc, "margin_force_terms", true);
    c.margin_extra = getOr(yc, "margin_extra", 0.0);
    c.rate_box = getOrB(yc, "rate_box", true);
    c.row_terrain = getOrB(yc, "row_terrain", false);
    c.hard_terrain = getOrB(yc, "hard_terrain", true);
    c.terrain_margin = getOr(yc, "terrain_margin", c.terrain_margin);
    c.terrain_gate = getOr(yc, "terrain_gate", c.terrain_gate);
    c.w_terrain = getOr(yc, "w_terrain", c.w_terrain);
    c.terrain_degree3 = getOrB(yc, "terrain_degree3", false);
    c.terrain_a_brk = getOr(yc, "terrain_a_brk", c.terrain_a_brk);
    c.terrain_v_safe = getOr(yc, "terrain_v_safe", c.terrain_v_safe);
    c.terrain_z_eps = getOr(yc, "terrain_z_eps", c.terrain_z_eps);
    c.terrain_k_neg = getOr(yc, "terrain_k_neg", c.terrain_k_neg);
    c.terrain_lookahead = getOr(yc, "terrain_lookahead", c.terrain_lookahead);
    c.terrain_smax_eps = getOr(yc, "terrain_smax_eps", c.terrain_smax_eps);
    c.terrain_climb_grad = getOr(yc, "terrain_climb_grad", c.terrain_climb_grad);
    c.terrain_fan_half = getOr(yc, "terrain_fan_half_deg", c.terrain_fan_half / kDeg) * kDeg;
    c.terrain_fan_n = static_cast<int>(getOr(yc, "terrain_fan_n", c.terrain_fan_n));
    if (yc && yc["c_terrain"])
      for (int i = 0; i < 3; ++i) c.c_terrain[i] = yc["c_terrain"][i].as<double>();
    if (c.row_terrain && !sc_.world.enabled)
      throw std::runtime_error("sixdof_sim: sixdof_cbf.row_terrain needs world.enabled");
    c.terrain = &sc_.world.terrain;
    c.best_effort_penalty = getOr(yc, "best_effort_penalty", c.best_effort_penalty);
    c.limits = nom.limits;
    c.dt = sc_.dt;
  }
  sc_.alpha_stall = getOr(root["envelope"], "alpha_stall_deg", sc_.alpha_stall / kDeg) * kDeg;


  x0_ = trim_.x;
  const double vscale = (sc_.V_app + sc_.dV) / std::max(1e-6, sc_.V_app);
  // Scale the whole velocity triad (the Beaver trims with small nonzero v =
  // V sin(beta); scaling preserves the trim alpha/beta).
  x0_[U] *= vscale;
  x0_[V] *= vscale;
  x0_[W] *= vscale;
  x0_[THETA] += sc_.dtheta;
  x0_[H] = sc_.h0;
  if (sc_.world.enabled) {
    // Absolute placement: the trim heading is 0 (north), so psi0 is absolute.
    x0_[PSI] += sc_.psi0;
    x0_[Y] = sc_.E0;
    x0_[XN] = sc_.N0;
  } else {
    x0_[PSI] += sc_.dpsi;
    x0_[Y] = sc_.dy;
  }
}

SixDofTouchdown SixDofSim::run(const std::string& csv_path) {
  const double dt = sc_.dt;
  const int nsteps = static_cast<int>(sc_.t_max / dt);

  const bool trim_ok = stats_.trim_converged;
  stats_ = SixDofRunStats{};
  stats_.trim_converged = trim_ok;

  SixDofNominal nominal(sc_.nominal);
  const WaveField wf = makeWaveField(sc_.waves);

  // Energy-reach CBF layer state (see the ReachSimCfg doc in sixdof_sim.hpp).
  // alpha_cmd is the filter's own AoA command integrator, seeded from the
  // measured alpha on the first step; beta commands are offsets from the trim
  // sideslip (the plant's straight-flight zero, nonzero on the Beaver).
  const EnergyReachFilter rfilter(sc_.reach.cbf);
  double alpha_cmd = NAN;
  const double Vt_trim = std::max(1e-6, trim_.x.head<3>().norm());
  const double beta_trim =
      std::asin(std::clamp(trim_.x[V] / Vt_trim, -1.0, 1.0));
  double beta_cmd_log = 0.0;

  // Corridor landing: pattern guidance planned ONCE from the initial pose
  // (ground course from the initial velocity), then tracked.
  const bool world_on = sc_.world.enabled;
  const Corridor& cor = sc_.world.corridor;
  PatternGuidance guidance(sc_.guidance, cor);
  if (world_on) {
    const GroundKinematics gk0 = groundKinematics(x0_);
    const double chi0 = std::atan2(gk0.ydot_e, gk0.xdot_n);
    if (!guidance.plan(x0_[XN], x0_[Y], chi0, x0_[H]))
      throw std::runtime_error("sixdof_sim: pattern guidance found no path");
    stats_.plan_length = guidance.sAimPath();
    stats_.d_final = guidance.dFinal();
    stats_.plan_word = guidance.word();
  }
  GuidanceCmd gc;  // zeros unless world_on

  // 6-DOF surfaces-only CBF filter (sixdof_cbf.enabled): sits between the
  // nominal's inner loops and the integrator; the throttle is the nominal's.
  const bool cbf_on = sc_.cbf.enabled;
  if (cbf_on && !bdyn_)
    throw std::runtime_error("sixdof_sim: sixdof_cbf needs the Beaver plant");
  std::unique_ptr<SixDofCbfFilter> cbf;
  if (cbf_on) cbf = std::make_unique<SixDofCbfFilter>(*bdyn_, sc_.cbf);
  SixDofCbfDiag cd;   // zeros unless the filter runs
  CtrlVec u_nom_log = CtrlVec::Zero();
  double cbf_time_sum = 0.0;

  // Gust penetration distance [m]. The earth-frame north position is the
  // plant state x[XN] (RK4 with everything else).
  double x_gust = 0.0;
  // Control held over the previous step (ZOH); seeds the exact airspeed-rate
  // input at t = 0 with the trim command the plant starts under.
  CtrlVec u_applied = trim_.u;

  std::ofstream csv(csv_path);
  // 10 significant digits: at the default 6, a 40 m/s airspeed quantizes to
  // 0.1 mm/s steps, which renders as a staircase on tightly-held channels.
  csv << std::setprecision(10);
  csv << "t,x,y,h,u,v,w,V_air,alpha_deg,beta_deg,p,q,r,"
         "phi_deg,theta_deg,psi_deg,gamma_deg,sink,"
         "de,da,dr,dT,theta_cmd_deg,phi_cmd_deg,"
         "W_u,W_v,W_h,x_gust,eta,eta_slope,"
         "Vdot_air,hdot_sp,ste_rate_sp,ste_rate_est,seb_rate_sp,seb_rate_est,"
         "tecs_pitch_int,tecs_thr_int,"
         "hE,d_tgt,phi_arc_deg,LD_inst,alpha_cmd_deg,beta_cmd_deg,"
         "ua_cbf,be_cbf_deg,cbf_viol,cbf_active,"
         "phase,chi_cmd_deg,N_ref,E_ref,h_gs,L_to_go,seg,s_corr,e_corr,hT,Vg,xte,s_path,"
         "cbf6_active,cbf6_be,h6_alpha,h6_V,h6_phi,h6_beta,de_nom,da_nom,dr_nom,cbf6_margin,cbf6_us,h6_terrain\n";

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
    const double eta_now = wf.eta(x[XN], t);
    const double eta_x = wf.slopeMean(x[XN], t, sc_.waves.contact_len);

    // Exact airspeed rate under the held control (see airspeedRate).
    const GustWind Wgdot =
        gustWindRate(sc_.wind, x_gust, gustXdot(sc_.wind, t, ad.V));
    const double Vdot_air = airspeedRate(
        x, plant_(x, u_applied, Eigen::Vector3d(Wg.u, Wg.v, Wg.w)), Wg, Wgdot);

    CtrlVec u;
    EnergyReachDiag rd;  // zeros unless the reach layer runs
    if (sc_.reach.enabled) {
      // Outer loop -> point-mass (alpha_dot, beta) space -> CBF-QP -> inner
      // loops on the filtered commands (wings-level + rudder sideslip).
      const SixDofNominal::LonOuterCmd oc =
          nominal.stepOuter(x, ad.V, Vdot_air, dt);
      const double dT =
          sc_.reach.throttle_idle ? sc_.nominal.limits.dT_min : oc.dT;
      const double chi = std::atan2(gk.ydot_e, gk.xdot_n);
      const double Vg = x.head<3>().norm();
      const ReachStateArr Xr{Vg, gamma, chi, x[XN], x[Y], x[H], ad.alpha};
      if (!std::isfinite(alpha_cmd)) alpha_cmd = ad.alpha;
      const EnergyReachCfg& ec = sc_.reach.cbf;
      // theta = alpha + gamma (wings level): the nominal theta_cmd maps to an
      // AoA request, rate-shaped through tau_alpha into alpha_dot space.
      const double alpha_des = std::clamp(oc.theta_cmd - gamma,
                                          ec.alpha_cmd_min, ec.alpha_cmd_max);
      const double ua_nom = std::clamp(
          (alpha_des - alpha_cmd) / sc_.reach.tau_alpha, -ec.ua_max, ec.ua_max);
      const Eigen::Vector2d uf = rfilter.filter(ua_nom, 0.0, Xr, &rd);
      alpha_cmd = std::clamp(alpha_cmd + uf[0] * dt, ec.alpha_cmd_min,
                             ec.alpha_cmd_max);
      const double beta_cmd = beta_trim + uf[1];
      beta_cmd_log = beta_cmd;
      u = nominal.stepInnerBeta(x, alpha_cmd + gamma, dT, ad.beta, beta_cmd, dt);
      stats_.min_hE = std::min(stats_.min_hE, rd.h);
      if (rd.active) ++stats_.cbf_active_steps;
      if (rd.best_effort) ++stats_.cbf_best_effort_steps;
    } else if (world_on) {
      // Pattern guidance -> (gamma_ref, phi_cmd) -> the same outer / inner
      // loops at those references.
      gc = guidance.step(x[XN], x[Y], x[H], gk.xdot_n, gk.ydot_e);
      const SixDofNominal::LonOuterCmd oc = nominal.stepOuterRef(
          x, ad.V, Vdot_air, gc.gamma_ref, sc_.nominal.V_ref, dt);
      u = nominal.stepInnerPhi(x, ad.V, oc.theta_cmd, oc.dT, gc.phi_cmd,
                               ad.beta, beta_trim, dt);
    } else {
      u = nominal.step(x, ad.V, Vdot_air, dt);
    }
    u_nom_log = u;
    if (cbf_on) {
      const SurfVec U_nom(u[DE], u[DA], u[DR]);
      const SurfVec up(u_applied[DE], u_applied[DA], u_applied[DR]);
      // Wind: the plant's own vector -- an idealized perfect wind estimate
      // (see sixdof_cbf.hpp); the nominal still sees only pitot air data.
      const SurfVec U = cbf->filter(U_nom, x, u[DT], up, &cd,
                                    Eigen::Vector3d(Wg.u, Wg.v, Wg.w));
      if (sc_.cbf.filter) { u[DE] = U[0]; u[DA] = U[1]; u[DR] = U[2]; }
      nominal.commitApplied(u);
      if (cd.active && sc_.cbf.filter) ++stats_.cbf6_active_steps;
      if (cd.best_effort) ++stats_.cbf6_best_effort_steps;
      stats_.cbf6_min_h_alpha = std::min(stats_.cbf6_min_h_alpha, cd.h_alpha);
      stats_.cbf6_min_h_V = std::min(stats_.cbf6_min_h_V, cd.h_V);
      stats_.cbf6_min_h_phi = std::min(stats_.cbf6_min_h_phi, cd.h_phi);
      stats_.cbf6_min_h_beta = std::min(stats_.cbf6_min_h_beta, cd.h_beta);
      for (const auto& r : cd.rows)
        if (r.name == "terrain") stats_.cbf6_min_h_terrain = std::min(stats_.cbf6_min_h_terrain, r.h);
      stats_.cbf6_margin_max = std::max(stats_.cbf6_margin_max, cd.margin_max);
      stats_.cbf6_time_max_us = std::max(stats_.cbf6_time_max_us, cd.t_us);
      cbf_time_sum += cd.t_us;
    }
    // Envelope verdict: the stall region is a failed run whatever else happens.
    if (ad.alpha > sc_.alpha_stall) { stats_.stall_entered = true; ++stats_.stall_steps; }
    const px4::TecsDebugOutput& tdb = nominal.tecsDebug();

    // Terrain under the aircraft and the corridor-frame position.
    double hT = 0.0, s_corr = 0.0, e_corr = 0.0;
    if (world_on) {
      hT = sc_.world.terrain.height(x[XN], x[Y]);
      cor.frame(x[XN], x[Y], s_corr, e_corr);
      // Clearance is only meaningful OVER terrain (over the water it is the
      // altitude, trivially 0 at touchdown).
      if (hT > sc_.world.terrain.base + 0.5) {
        stats_.min_terrain_clearance =
            std::min(stats_.min_terrain_clearance, x[H] - hT);
        stats_.min_wingtip_clearance = std::min(
            stats_.min_wingtip_clearance,
            x[H] - hT - 0.5 * BeaverGeom::b * std::abs(std::sin(x[PHI])));
      }
      if (gc.on_final)
        stats_.max_abs_e_final =
            std::max(stats_.max_abs_e_final, std::abs(e_corr));
      stats_.max_xte = std::max(stats_.max_xte, gc.xte);
      if (std::abs(gc.phi_cmd) >= sc_.guidance.phi_max - 1e-9)
        ++stats_.phi_sat_steps;
    }
    const double Vg_ground = std::hypot(gk.xdot_n, gk.ydot_e);

    stats_.max_abs_phi = std::max(stats_.max_abs_phi, std::abs(x[PHI]));
    stats_.max_abs_beta = std::max(stats_.max_abs_beta, std::abs(ad.beta));
    stats_.max_alpha = std::max(stats_.max_alpha, ad.alpha);
    stats_.min_V = std::min(stats_.min_V, ad.V);
    stats_.max_abs_y = std::max(stats_.max_abs_y, std::abs(x[Y]));

    csv << t << ',' << x[XN] << ',' << x[Y] << ',' << x[H] << ',' << x[U] << ','
        << x[V] << ',' << x[W] << ',' << ad.V << ',' << ad.alpha / kDeg << ','
        << ad.beta / kDeg << ',' << x[P] << ',' << x[Q] << ',' << x[R] << ','
        << x[PHI] / kDeg << ',' << x[THETA] / kDeg << ',' << x[PSI] / kDeg << ','
        << gamma / kDeg << ',' << sink << ',' << u[DE] << ',' << u[DA] << ','
        << u[DR] << ',' << u[DT] << ',' << nominal.thetaCmd() / kDeg << ','
        << nominal.phiCmd() / kDeg << ',' << Wg.u << ',' << Wg.v << ',' << Wg.w
        << ',' << x_gust << ',' << eta_now << ',' << eta_x << ',' << Vdot_air
        << ',' << nominal.hdotSp() << ',' << tdb.total_energy_rate_sp << ','
        << tdb.total_energy_rate_estimate << ',' << tdb.energy_balance_rate_sp
        << ',' << tdb.energy_balance_rate_estimate << ','
        << tdb.pitch_integrator << ',' << tdb.throttle_integrator << ','
        << rd.h << ',' << rd.d << ',' << rd.phi / kDeg << ',' << rd.LD << ','
        << (std::isfinite(alpha_cmd) ? alpha_cmd / kDeg : 0.0) << ','
        << beta_cmd_log / kDeg << ',' << rd.ua << ',' << rd.be / kDeg << ','
        << rd.viol << ',' << (rd.active ? 1 : 0) << ','
        << 0 << ',' << gc.chi_cmd / kDeg << ',' << gc.N_ref << ',' << gc.E_ref
        << ',' << gc.h_gs << ',' << gc.L_to_go << ',' << gc.seg << ','
        << s_corr << ',' << e_corr << ',' << hT << ',' << Vg_ground << ','
        << gc.xte << ',' << gc.s_path << ','
        << (cd.active ? 1 : 0) << ',' << (cd.best_effort ? 1 : 0) << ','
        << cd.h_alpha << ',' << cd.h_V << ',' << cd.h_phi << ',' << cd.h_beta << ','
        << u_nom_log[DE] << ',' << u_nom_log[DA] << ',' << u_nom_log[DR] << ','
        << cd.margin_max << ',' << cd.t_us << ',' << cd.h_terrain << '\n';

    // Contact surface: the wave elevation, or the terrain where it stands
    // above the water (a terrain strike is a crash, not a landing).
    const double surface = world_on ? std::max(eta_now, hT) : eta_now;
    if (x[H] <= surface && k > 0) {
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
      td.x = x[XN];
      td.h = x[H];
      td.eta = eta_now;
      td.slope = eta_x;
      if (sc_.reach.enabled) {
        const ReachParams& pp = sc_.reach.cbf.prm;
        td.hE = rd.h;
        td.dist_target = std::hypot(pp.Nt - x[XN], pp.Et - x[Y]);
      }
      // Surface-relative closure -d/dt(h - eta) = sink + eta_x xdot + eta_t.
      td.sink_rel = sink + eta_x * gk.xdot_n + wf.etaDot(x[XN], t);
      // TN 1516 slam-load truth, flat- vs wave-referenced (tau / gamma0
      // tilted by the local surface angle, closure onto the moving surface).
      const double m = ac_.inertia.mass, g = ac_.env.g;
      const double alpha_s = std::atan(eta_x);
      td.n_peak_flat = impactNPeakExact(
          m, g, sc_.hull.beta, sc_.hull.rho_water, sc_.hull.eps_g0,
          x[THETA] - sc_.hull.tau_keel, -gamma, std::max(0.0, sink),
          sc_.hull.n_surfaces);
      td.n_peak_wave = impactNPeakExact(
          m, g, sc_.hull.beta, sc_.hull.rho_water, sc_.hull.eps_g0,
          x[THETA] - alpha_s - sc_.hull.tau_keel, alpha_s - gamma,
          std::max(0.0, td.sink_rel), sc_.hull.n_surfaces);
      if (world_on) {
        const double chi = std::atan2(gk.ydot_e, gk.xdot_n);
        td.hT = hT;
        td.terrain_strike = hT > eta_now + 0.05;
        td.s_corr = s_corr;
        td.e_corr = e_corr;
        td.dpsi_corr = wrapAngle(x[PSI] - cor.heading);
        td.dchi_corr = wrapAngle(chi - cor.heading);
        td.crab = wrapAngle(x[PSI] - chi);
        td.in_corridor = !td.terrain_strike && cor.contains(x[XN], x[Y]);
        if (!td.terrain_strike && sc_.world.rollout.enabled) {
          // Water rollout: continue the CSV with phase = 1 rows (position,
          // course as psi, ground speed, corridor frame; the rest zero).
          RolloutState rs;
          rs.N = x[XN]; rs.E = x[Y]; rs.chi = chi; rs.Vg = Vg_ground;
          auto row = [&](double tt, const RolloutState& st) {
            if (tt == td.t) return;  // the airborne row at t already exists
            double sc, ec;
            cor.frame(st.N, st.E, sc, ec);
            csv << tt << ',' << st.N << ',' << st.E << ",0,0,0,0,0,0,0,0,0,0,"
                << "0,0," << st.chi / kDeg << ",0,0,0,0,0,0,0,0,"
                << "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"
                << "1,0,0,0,0,0,0," << sc << ',' << ec << ",0," << st.Vg
                << ",0,0,0,0,0,0,0,0,0,0,0,0,0,0\n";
          };
          td.rollout = runRollout(rs, td.t, dt, sc_.world.rollout, cor, row);
        }
        std::string why;
        if (td.terrain_strike) why += "terrain strike; ";
        else if (!td.in_corridor) why += "touchdown outside corridor; ";
        if (sc_.world.rollout.enabled && !td.terrain_strike) {
          if (!td.rollout.stopped) why += "rollout timed out; ";
          else if (!td.rollout.inside) why += "rollout stopped outside corridor; ";
          else if (!td.rollout.stayed_inside) why += "rollout left the corridor; ";
        }
        if (stats_.stall_entered) why += "entered stall region; ";
        td.success = why.empty();
        td.fail_reasons = why;
      }
      break;
    }

    const JointState js = rk4WindStep(plant_, sc_.wind, t, x, x_gust, u, dt);
    x = js.x;
    x_gust = js.xg;
    u_applied = u;
    t += dt;
    ++stats_.steps;
  }
  stats_.t_end = t;
  if (cbf_on && stats_.steps > 0) stats_.cbf6_time_mean_us = cbf_time_sum / (stats_.steps + 1);

  std::cout << "=== 6-DOF straight-in landing ===\n";
  std::cout << "plant: "
            << (bdyn_ ? "DHC-2 Beaver (LR-556/FDC polynomials)"
                      : "AHAB VSPAERO deck")
            << "\n";
  std::cout << "trim: theta=" << trim_.theta / kDeg << " deg  de="
            << trim_.u[DE] / kDeg << " deg  dT=" << trim_.u[DT];
  if (bdyn_)
    std::cout << " (pz=" << bdyn_->pzFromThrottle(trim_.u[DT]) << " \"Hg, P="
              << bdyn_->power(trim_.u[DT]) << " kW)  da="
              << trim_.u[DA] / kDeg << " deg  dr=" << trim_.u[DR] / kDeg
              << " deg";
  std::cout << "  (V_app=" << sc_.V_app << " m/s, gamma_app="
            << sc_.gamma_app / kDeg << " deg)\n";
  if (sc_.nominal.lon_mode == LonMode::Tecs) {
    const px4::TecsParam& tp = sc_.nominal.tecs;
    std::cout << "nominal: PX4 TECS (direct height-rate hdot_sp="
              << sc_.nominal.V_ref * std::sin(sc_.nominal.gamma_ref)
              << " m/s, TAS_sp=" << sc_.nominal.V_ref << " m/s)\n"
              << "  anchors: thr_trim=" << tp.throttle_trim
              << "  pitch_offset=" << sc_.nominal.tecs_pitch_offset / kDeg
              << " deg  clmb_max=" << tp.max_climb_rate
              << " m/s  sink_min=" << tp.min_sink_rate
              << " m/s  sink_max=" << tp.max_sink_rate << " m/s (at tas_max="
              << tp.tas_max << ")  tas_min=" << tp.tas_min << "\n";
  } else {
    std::cout << "nominal: cascaded PID (gamma -> theta, V -> throttle)\n";
  }
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
  if (stats_.stall_entered)
    std::cout << "ENVELOPE: entered the stall region (alpha > " << sc_.alpha_stall / kDeg
              << " deg) for " << stats_.stall_steps << " steps -- FAILED run; the plant has "
                 "no post-stall aero, numbers past it are not physical (max alpha="
              << stats_.max_alpha / kDeg << " deg)\n";
  if (cbf_on) {
    const SixDofCbfConfig& c = sc_.cbf;
    std::cout << "6-DOF CBF (surfaces only, " << (c.filter ? "FILTERING" : "monitor only")
              << "): rows alpha<=" << c.alpha_max / kDeg << " deg" << (c.hard_alpha ? "(hard)" : "")
              << "  V in [" << c.V_min << "," << c.V_max << "]" << (c.hard_V ? "(hard)" : "")
              << "  |phi|<=" << c.phi_max / kDeg << " deg" << (c.hard_phi ? "(hard)" : "")
              << "  |beta|<=" << c.beta_max / kDeg << " deg" << (c.hard_beta ? "(hard)" : "") << "\n"
              << "  active " << stats_.cbf6_active_steps << "/" << stats_.steps + 1
              << " steps  best-effort=" << stats_.cbf6_best_effort_steps
              << "  min h: alpha=" << stats_.cbf6_min_h_alpha / kDeg << " deg  V="
              << stats_.cbf6_min_h_V << " m/s  phi=" << stats_.cbf6_min_h_phi / kDeg
              << " deg  beta=" << stats_.cbf6_min_h_beta / kDeg << " deg";
    if (c.row_terrain)
      std::cout << "  terrain=" << (stats_.cbf6_min_h_terrain < 1e29 ? stats_.cbf6_min_h_terrain : NAN)
                << " m (margin " << c.terrain_margin << ", " << (c.hard_terrain ? "hard" : "soft") << ")";
    std::cout << "\n"
              << "  robustness margin max=" << stats_.cbf6_margin_max
              << "  compute mean=" << stats_.cbf6_time_mean_us << " us  max="
              << stats_.cbf6_time_max_us << " us per step\n";
  }
  if (world_on) {
    std::cout << "corridor landing: corridor " << cor.length << " x "
              << cor.width << " m at heading " << cor.heading / kDeg
              << " deg, aim s=" << cor.s_aim << " m; terrain bumps="
              << sc_.world.terrain.bumps.size() << "\n"
              << "  plan: Dubins " << stats_.plan_word << "  R="
              << guidance.turnRadius() << " m  d_final=" << stats_.d_final
              << " m  path to aim=" << stats_.plan_length
              << " m  (descent at gamma_app="
              << sc_.h0 / std::tan(-sc_.gamma_app) << " m)\n"
              << "  min terrain clearance=";
    if (stats_.min_terrain_clearance < 1e29)
      std::cout << stats_.min_terrain_clearance << " m (CG), "
                << stats_.min_wingtip_clearance << " m (low wingtip, b/2 sin|phi|)";
    else
      std::cout << "n/a (never over terrain)";
    std::cout << "  max |e| on final=" << stats_.max_abs_e_final
              << " m  max cross-track to path=" << stats_.max_xte
              << " m  bank-limit steps=" << stats_.phi_sat_steps << "/"
              << stats_.steps + 1 << "\n";
    if (td.reached) {
      if (td.terrain_strike)
        std::cout << "  TERRAIN STRIKE at N=" << td.x << " E=" << td.y
                  << " (h_T=" << td.hT << " m)\n";
      else
        std::cout << "  touchdown " << (td.in_corridor ? "IN" : "OUTSIDE")
                  << " corridor: s=" << td.s_corr << " m  e=" << td.e_corr
                  << " m  course err=" << td.dchi_corr / kDeg
                  << " deg  heading err=" << td.dpsi_corr / kDeg
                  << " deg  crab=" << td.crab / kDeg << " deg  bank="
                  << td.phi / kDeg << " deg (no flare / decrab / level-off)\n";
      std::cout << "  RESULT: " << (td.success ? "LANDED OK" : "FAILED -- " + td.fail_reasons) << "\n";
      if (td.rollout.ran)
        std::cout << "  rollout: " << (td.rollout.stopped ? "stopped" : "TIMED OUT")
                  << " at t=" << td.rollout.t_stop << " s  s=" << td.rollout.s_stop
                  << " m  e=" << td.rollout.e_stop << " m  run="
                  << td.rollout.distance << " m  inside="
                  << (td.rollout.inside ? "yes" : "NO") << "  stayed inside="
                  << (td.rollout.stayed_inside ? "yes" : "NO") << "\n";
    }
  }
  if (sc_.reach.enabled) {
    const ReachParams& pp = sc_.reach.cbf.prm;
    std::cout << "energy-reach CBF ("
              << (sc_.reach.filter ? "FILTERING" : "monitor only") << ", "
              << (sc_.reach.throttle_idle ? "idle glide" : "nominal throttle")
              << "):\n"
              << "  polar fit @ idle: CLa=" << polar_.CLa
              << "/rad  CD0=" << polar_.CD0 << "  CDa2=" << polar_.CDa2
              << "  alpha0L=" << polar_.alpha0L / kDeg << " deg  (L/D)max="
              << polar_.LD_max << " at at*=" << polar_.alpha_star / kDeg
              << " deg\n"
              << "  target N=" << pp.Nt << " m  E=" << pp.Et << " m  V_min="
              << pp.V_min << " m/s  gamma_app=" << sc_.gamma_app / kDeg
              << " deg\n"
              << "  min h_E=" << stats_.min_hE << " m^2/s^2  active "
              << stats_.cbf_active_steps << "/" << stats_.steps + 1
              << " steps  best-effort=" << stats_.cbf_best_effort_steps << "\n";
    if (td.reached)
      std::cout << "  touchdown: " << td.dist_target
                << " m from target  h_E=" << td.hE << "  V=" << td.V
                << " m/s (V_min=" << pp.V_min << ")\n";
  }
  std::cout << "  trace -> " << csv_path << "\n";
  return td;
}

}  // namespace autoland
