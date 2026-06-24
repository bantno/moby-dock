#pragma once
#include <string>
#include "autoland/types.hpp"

// =============================================================================
// Configuration structures + YAML loaders (aircraft.yaml, scenario.yaml).
//
// Anything marked TODO must be supplied/verified by the user. In particular the
// MASS PROPERTIES and the MIXING MAP are placeholders -- see data/aircraft.yaml.
// =============================================================================
namespace autoland {

struct InertiaParams {
  double mass{0};           // [kg]                 TODO: confirm
  double Ixx{0}, Iyy{0}, Izz{0}, Ixz{0};  // [kg m^2]  TODO: confirm
};

// Offset of the true c.g. from the .stab moment reference (Xcg_,Ycg_,Zcg_),
// expressed in body axes [m]. Zero => c.g. at the .stab reference point.
struct CGOffset {
  double dx{0}, dy{0}, dz{0};  // TODO: confirm if c.g. != .stab reference
};

// Default thrust model: T = throttle * (T_static - k_v * V), clamped >= 0.
// Hook: replace with measured thrust-stand data (see ThrustModel in dynamics).
struct ThrustParams {
  double T_static{0};  // static thrust at throttle=1, V=0     [N]
  double k_v{0};       // linear velocity falloff              [N/(m/s)]
  double zcp{0};       // thrust line z-offset below c.g.      [m] (pitch coupling)
};

struct Environment {
  double rho{1.225};       // air density [kg/m^3]  (NOT taken from .stab Rho_)
  double a_sound{340.3};   // speed of sound [m/s]
  double g{9.80665};       // gravity [m/s^2]
};

// Per-virtual-control deflection and rate limits.
struct SurfaceLimits {
  double de_min{-0.5}, de_max{0.5};   // [rad]
  double da_min{-0.5}, da_max{0.5};
  double dr_min{-0.5}, dr_max{0.5};
  double rate{5.0};                   // max deflection rate [rad/s]
  double dT_min{0.0}, dT_max{1.0};    // throttle bounds [-]
};

struct AircraftConfig {
  InertiaParams inertia;
  CGOffset cg;
  ThrustParams thrust;
  Environment env;
  SurfaceLimits limits;

  // Extra parasite drag (body-x) representing viscous + hull drag that the
  // INVISCID VSPAero panel solution omits. Applied as a constant body-x drag
  // coefficient: CFx -= parasite_CD0. Without it the model has near-zero drag,
  // so shallow descents require ~zero/negative thrust and trim is ill-posed.
  // TODO: calibrate against flight/tow-tank data.
  double parasite_CD0{0.0};

  // Mixing matrix M (numGroups x 3) mapping virtual [delta_e, delta_a, delta_r]
  // -> physical control-group deflections, group order matching the .stab.
  // Empty => build the default identity-style map at load time (see mixing.cpp).
  Mat mixing;                       // TODO: confirm against OpenVSP groups
  std::vector<std::string> mixing_group_order;  // optional name check
};

// Cascaded-PID gains (frontside technique). See controller.hpp.
struct Gains {
  // Airspeed loop (PI on V error -> throttle)
  double Kp_V{0}, Ki_V{0};
  // Glidepath outer loop (altitude/glideslope error -> flight-path/theta cmd)
  double Kp_gs{0}, Kd_gs{0};
  // Pitch inner loop (theta error + pitch-rate damping -> delta_e)
  double Kp_theta{0}, Ki_theta{0}, Kq{0};
  // Cross-track loop (-> phi_cmd, limited). Kd_y damps the bank->heading->y
  // double integrator (cross-track-rate feedback) and is essential for stability.
  double Kp_y{0}, Kd_y{0}, Ki_y{0}, phi_max{0.5};
  // Roll inner loop (phi error + roll-rate damping -> delta_a)
  double Kp_phi{0}, Kp_p{0};
  // Yaw damper (r -> delta_r)
  double Kr{0};
  // Decrab (sideslip -> delta_r near touchdown)
  double Kp_beta{0};
};

struct FlareParams {
  double h_flare{0};   // flare engage height AGL [m]
  double tau{0};       // sink-rate decay time constant [s]
  double w_td{0};      // target touchdown sink rate (positive down) [m/s]
};

struct InitCond {
  double distance{0};  // initial horizontal distance to touchdown [m]
  double dh{0};        // initial altitude error vs glideslope [m]
  double dy{0};        // initial cross-track offset [m]
  double dpsi{0};      // initial heading error [rad]
  double dV{0};        // initial airspeed error [m/s]
};

struct ScenarioConfig {
  double V_app{0};     // target approach airspeed [m/s]
  double gamma_app{0}; // target flight-path angle [rad] (negative = descent)
  double h_decrab{0};  // height below which decrab engages [m]
  Gains gains;
  FlareParams flare;
  InitCond init;
  double dt{0.01};     // integrator step [s]
  double t_max{120.0}; // max sim duration [s]
};

AircraftConfig loadAircraftConfig(const std::string& path);
ScenarioConfig loadScenarioConfig(const std::string& path);

}  // namespace autoland
