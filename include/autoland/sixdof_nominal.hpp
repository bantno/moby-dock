#pragma once
#include <algorithm>
#include <cfloat>
#include <cmath>
#include "autoland/config.hpp"
#include "autoland/px4_tecs.hpp"
#include "autoland/types.hpp"

// =============================================================================
// Nominal controller for the 6-DOF straight-in approach. Two selectable
// LONGITUDINAL outer loops (SixDofNominalConfig::lon_mode) feed one shared
// pitch inner loop; the lateral axes are common.
//
// CASCADE (default) -- successive loop closure (Beard & McLain, "Small
// Unmanned Aircraft", 2012):
//   Airspeed  : PI on V error -> throttle (frontside technique).
//   Glidepath : inertial flight-path angle gamma -> theta_cmd (PI, clamped),
//               mirroring the lon sim's gamma cascade (lon_nominal.hpp).
//
// TECS -- the PX4 Total Energy Control System control law (px4_tecs.hpp, a
// port of PX4-Autopilot TECSControl): total-energy rate -> throttle and
// energy-balance rate -> pitch. It is given PX4's DIRECT height-rate setpoint
// hdot_sp = V_ref sin(gamma_ref) and the TAS setpoint V_ref -- the same
// references the cascade tracks -- so the two nominals differ only in the
// outer energy loops. Its pitch setpoint is relative to the LEVEL-flight trim
// pitch (PX4's FW_PSP_OFF; tecs_pitch_offset), and its throttle feedforward
// is anchored on plant trims (level / full-throttle climb / idle sink at
// V_ref) that the sim solves at construction (sixdof_sim.cpp). The exact
// airspeed rate the sim supplies stands in for PX4's airspeed filter.
//
// Shared:
//   Pitch     : theta_cmd error + pitch-rate damping -> delta_e (PID).
//   Lateral   : cross-track y (+ rate) -> phi_cmd (PD, clamped);
//               roll inner PD -> delta_a; yaw damper r -> delta_r.
//
// No flare/decrab: the aircraft flies the glideslope into the water and crabs
// into a crosswind. Full-state feedback; all feedforwards come from the trim
// solve. The airspeed input is the measured (pitot) airspeed -- the wind
// VECTOR stays unmeasured, per the repo's plant-disturbance philosophy.
// =============================================================================
namespace autoland {

enum class LonMode { Cascade, Tecs };

struct SixDofNominalConfig {
  double V_ref{18.0};        // approach airspeed [m/s]
  double gamma_ref{-0.052};  // inertial flight-path angle [rad] (negative=descent)
  // Trim feedforward (seeded from the Newton trim at (V_ref, gamma_ref)).
  // da/dr matter for the Beaver plant, whose slipstream asymmetry trims with
  // small nonzero aileron/rudder (zero on the symmetric VSPAERO deck).
  double theta_trim{0.0}, de_trim{0.0}, dT_trim{0.0};
  double da_trim{0.0}, dr_trim{0.0};
  // Airspeed loop: throttle = dT_trim + Kp_V eV + Ki_V int(eV).
  double Kp_V{0.05}, Ki_V{0.02};
  // Glidepath: gamma_ref_eff = gamma_ref + clamp(Kv_gamma (V - V_ref)), then
  //            theta_cmd = theta_trim + Kp_gamma e_gam + Ki_gamma int(e_gam).
  // Kv_gamma is the speed <-> path energy trade (a one-line TECS: fast ->
  // command a SHALLOWER path so gravity stops feeding the speed; slow ->
  // steeper). It keeps the speed axis stable whenever the throttle rails at
  // idle (steep approaches / hot entries), where drag alone recovers V only
  // slowly. NOTE it must live in the REFERENCE, not inside theta_cmd: any
  // speed term added to theta_cmd is exactly canceled by the gamma
  // integrator at steady state. As a reference shift the V dynamics gain
  // Vdot ~ -g Kv_gamma (V - V_ref): stable with tau = 1/(g Kv_gamma) even
  // with zero drag slope, and gamma still settles at gamma_ref once V does.
  double Kp_gamma{2.0}, Ki_gamma{0.5};
  double Kv_gamma{0.02};                        // [rad/(m/s)] -> tau ~ 5 s
  double dgamma_V_max{4.0 * M_PI / 180.0};      // cap on the reference shift
  double theta_cmd_max{25.0 * M_PI / 180.0};  // |theta_cmd - theta_trim| limit
  // Pitch inner PID (P gains sized for the real Iyy = 0.052, cf. lon_scenario).
  // Ki_theta closes the DC gap of the fixed de_trim feedforward: the elevator
  // needed to hold a given theta drifts with flight condition, and without an
  // integrator theta lags theta_cmd by (de_needed - de_trim)/Kp_theta -- the
  // gamma loop then never reaches its reference.
  double Kp_theta{1.0}, Ki_theta{0.5}, Kq{0.25};
  // Cross-track -> bank: phi_cmd = -Kp_y y - Kd_y ydot, clamped to phi_max.
  double Kp_y{0.02}, Kd_y{0.05};
  double phi_max{25.0 * M_PI / 180.0};
  // Roll inner PD and yaw damper. Roll gains are ceiling-limited by the
  // dt = 0.01 discrete step, not by performance: the deck's roll authority is
  // ~460 rad/s^2 per rad of aileron, so the rate-loop pole is ~460*Kp_p and
  // Kp_p = 0.5 already limit-cycles at the Nyquist rate (chatter band on
  // delta_a). 0.15 puts the pole at ~69 rad/s (0.69 per step).
  double Kp_phi{1.0}, Kp_p{0.15};
  double Kr{0.2};
  // Control-effectiveness SENSE per axis, multiplying the feedback terms
  // (never the trim feedforward). +1 = the AHAB virtual-control sense
  // (positive de -> nose UP, positive da -> roll RIGHT, positive dr -> yaw
  // damping via -Kr r). The Beaver plant uses the standard Delft/FDC signs
  // (Cm_de < 0, Cl_da < 0, Cn_dr < 0), i.e. -1 on all three.
  double de_sign{1.0}, da_sign{1.0}, dr_sign{1.0};
  SurfaceLimits limits;  // deflection + rate limits (aircraft.yaml)

  // --- Longitudinal outer loop: cascade (above) or the PX4 TECS port. -------
  LonMode lon_mode{LonMode::Cascade};
  // TECS parameters (PX4 FW_T_* defaults in px4_tecs.hpp). The vehicle
  // anchors -- throttle_trim, max_climb_rate, min_sink_rate, max_sink_rate,
  // the pitch/throttle limits and the TAS band -- are filled by the sim from
  // the plant and the nominal limits; the pitch limits are relative to
  // tecs_pitch_offset (PX4 FW_P_LIM_MIN/MAX - FW_PSP_OFF).
  px4::TecsParam tecs;
  double tecs_pitch_offset{0.0};   // FW_PSP_OFF analog: level-flight trim theta at V_ref [rad]
  bool tecs_detect_underspeed{true};
};

class SixDofNominal {
 public:
  explicit SixDofNominal(const SixDofNominalConfig& cfg) : c_(cfg) {}

  // One control step. x is the full (inertial) state; V_air the measured
  // airspeed and Vdot_air its (exact) time derivative -- consumed by the TECS
  // outer loop only (SKE rate = V Vdot). Returns absolute virtual controls
  // [de, da, dr, dT] with the deflection and rate limits applied.
  CtrlVec step(const StateVec& x, double V_air, double Vdot_air, double dt) {
    const double sp = std::sin(x[PHI]), cp = std::cos(x[PHI]);
    const double st = std::sin(x[THETA]), ct = std::cos(x[THETA]);
    const double sy = std::sin(x[PSI]), cy = std::cos(x[PSI]);

    // Inertial climb/cross-track rates from the state kinematics (same rows
    // as Dynamics::xdot), and the inertial flight-path angle.
    const double hdot = x[U] * st - x[V] * sp * ct - x[W] * cp * ct;
    const double ydot = x[U] * ct * sy + x[V] * (sp * st * sy + cp * cy) +
                        x[W] * (cp * st * sy - sp * cy);
    const double Vg = std::max(
        1e-3, std::sqrt(x[U] * x[U] + x[V] * x[V] + x[W] * x[W]));
    const double gamma = std::asin(std::clamp(hdot / Vg, -1.0, 1.0));

    CtrlVec u;
    double dT, theta_cmd;

    if (c_.lon_mode == LonMode::Tecs) {
      // --- PX4 TECS: total energy -> throttle, energy balance -> pitch. -----
      // Direct height-rate setpoint (PX4 altitude_rate_setpoint_direct) and
      // TAS setpoint. The altitude reference is only consulted by
      // initialize(), which PX4 seeds with the CURRENT altitude and climb
      // rate (TECS::initialize -> handle_alt_step semantics) -- mirrored here.
      px4::TecsParam prm = c_.tecs;
      // PX4: _load_factor_from_bank_angle = 1 / max(cos(phi), FLT_EPSILON).
      prm.load_factor = 1.0 / std::max(cp, double(FLT_EPSILON));
      const double hdot_sp = c_.V_ref * std::sin(c_.gamma_ref);
      px4::TecsSetpoint tsp;
      tsp.altitude_reference.alt = x[H];
      tsp.altitude_reference.alt_rate = hdot;
      tsp.altitude_rate_setpoint_direct = hdot_sp;
      tsp.tas_setpoint = c_.V_ref;
      const px4::TecsInput tin{x[H], hdot, V_air, Vdot_air};
      const px4::TecsFlag flag{true, c_.tecs_detect_underspeed};
      if (!tecs_initialized_) {
        // PX4 runs initialize() (not update()) on the first call. Note
        // initialize() has no direct-rate path: its demand is the altitude
        // loop's HRATE_FF * (current climb rate), i.e. half the glideslope
        // energy rate at trim -- the one-sample throttle/pitch kick at t = 0
        // in the traces is upstream behaviour, not a port artefact.
        tecs_.initialize(tsp, tin, prm, flag);
        tecs_initialized_ = true;
      } else {
        tecs_.update(dt, tsp, tin, prm, flag);
      }
      dT = tecs_.getThrottleSetpoint();  // already within [dT_min, dT_max]
      theta_cmd = c_.tecs_pitch_offset + tecs_.getPitchSetpoint();
      hdot_sp_ = hdot_sp;
    } else {
      // --- Airspeed -> throttle (PI, frontside; anti-windup on the clamp). -
      const double eV = c_.V_ref - V_air;
      V_int_ += eV * dt;
      dT = c_.dT_trim + c_.Kp_V * eV + c_.Ki_V * V_int_;
      if (dT > c_.limits.dT_max) { dT = c_.limits.dT_max; V_int_ -= eV * dt; }
      else if (dT < c_.limits.dT_min) { dT = c_.limits.dT_min; V_int_ -= eV * dt; }

      // --- Glidepath: speed-shifted gamma reference -> theta_cmd (PI). ------
      const double gamma_ref_eff =
          c_.gamma_ref + std::clamp(c_.Kv_gamma * (V_air - c_.V_ref),
                                    -c_.dgamma_V_max, c_.dgamma_V_max);
      const double e_gamma = gamma_ref_eff - gamma;
      gamma_int_ += e_gamma * dt;
      theta_cmd = c_.theta_trim + c_.Kp_gamma * e_gamma +
                  c_.Ki_gamma * gamma_int_;
      const double th_lo = c_.theta_trim - c_.theta_cmd_max;
      const double th_hi = c_.theta_trim + c_.theta_cmd_max;
      if (theta_cmd > th_hi) { theta_cmd = th_hi; gamma_int_ -= e_gamma * dt; }
      else if (theta_cmd < th_lo) { theta_cmd = th_lo; gamma_int_ -= e_gamma * dt; }
      hdot_sp_ = c_.V_ref * std::sin(gamma_ref_eff);
    }
    u[DT] = dT;
    theta_cmd_ = theta_cmd;

    // --- Pitch inner PID -> elevator (anti-windup on the deflection clamp). --
    const double e_theta = theta_cmd - x[THETA];
    theta_int_ += e_theta * dt;
    double de = c_.de_trim +
                c_.de_sign * (c_.Kp_theta * e_theta +
                              c_.Ki_theta * theta_int_ - c_.Kq * x[Q]);
    if (de > c_.limits.de_max) { de = c_.limits.de_max; theta_int_ -= e_theta * dt; }
    else if (de < c_.limits.de_min) { de = c_.limits.de_min; theta_int_ -= e_theta * dt; }
    u[DE] = de;

    // --- Cross-track -> bank command (PD, clamped). --------------------------
    double phi_cmd = -c_.Kp_y * x[Y] - c_.Kd_y * ydot;
    phi_cmd = std::clamp(phi_cmd, -c_.phi_max, c_.phi_max);
    phi_cmd_ = phi_cmd;

    // --- Roll inner PD -> aileron; yaw damper -> rudder. ---------------------
    u[DA] = c_.da_trim +
            c_.da_sign * (c_.Kp_phi * (phi_cmd - x[PHI]) - c_.Kp_p * x[P]);
    u[DR] = c_.dr_trim + c_.dr_sign * (-c_.Kr * x[R]);

    return applyLimits(u, dt);
  }

  double thetaCmd() const { return theta_cmd_; }
  double phiCmd() const { return phi_cmd_; }
  // Height-rate reference [m/s]: the TECS direct setpoint, or for the cascade
  // the climb rate at V_ref along its (speed-shifted) gamma reference.
  double hdotSp() const { return hdot_sp_; }
  bool tecsActive() const { return c_.lon_mode == LonMode::Tecs; }
  const px4::TecsDebugOutput& tecsDebug() const { return tecs_.getDebugOutput(); }
  void reset() {
    V_int_ = gamma_int_ = theta_int_ = 0.0;
    have_prev_ = false;
    tecs_ = px4::TecsControl{};
    tecs_initialized_ = false;
  }

 private:
  // Deflection clamps + surface rate limit against the previously APPLIED
  // command (throttle is clamped in the loop; no rate limit on it).
  CtrlVec applyLimits(CtrlVec u, double dt) {
    const SurfaceLimits& L = c_.limits;
    u[DE] = std::clamp(u[DE], L.de_min, L.de_max);
    u[DA] = std::clamp(u[DA], L.da_min, L.da_max);
    u[DR] = std::clamp(u[DR], L.dr_min, L.dr_max);
    if (have_prev_) {
      const double du = L.rate * dt;
      for (int i : {DE, DA, DR})
        u[i] = std::clamp(u[i], u_prev_[i] - du, u_prev_[i] + du);
    }
    u_prev_ = u;
    have_prev_ = true;
    return u;
  }

  SixDofNominalConfig c_;
  double V_int_{0.0};
  double gamma_int_{0.0};
  double theta_int_{0.0};
  double theta_cmd_{0.0};
  double phi_cmd_{0.0};
  double hdot_sp_{0.0};
  px4::TecsControl tecs_;
  bool tecs_initialized_{false};
  CtrlVec u_prev_{CtrlVec::Zero()};
  bool have_prev_{false};
};

}  // namespace autoland
