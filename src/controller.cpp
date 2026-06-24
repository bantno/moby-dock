#include "autoland/controller.hpp"

#include <algorithm>
#include <cmath>

namespace autoland {
namespace {

double clampd(double v, double lo, double hi) {
  return std::max(lo, std::min(hi, v));
}

// Sink rate (positive down) from the full state, same kinematics as Dynamics.
double sinkRate(const StateVec& x) {
  const double ct = std::cos(x[THETA]), st = std::sin(x[THETA]);
  const double cp = std::cos(x[PHI]), sp = std::sin(x[PHI]);
  const double hdot = x[U] * st - x[V] * sp * ct - x[W] * cp * ct;
  return -hdot;
}

}  // namespace

Controller::Controller(const Gains& gains, const FlareParams& flare,
                       const TrimResult& trim, double V_app, double gamma_app,
                       double h_decrab, const Environment& env,
                       const SurfaceLimits& limits)
    : g_(gains),
      flare_(flare),
      V_app_(V_app),
      gamma_app_(gamma_app),
      h_decrab_(h_decrab),
      g_accel_(env.g),
      lim_(limits),
      x_trim_(trim.x),
      u_trim_(trim.u),
      theta_trim_(trim.theta) {}

ControllerOutputs Controller::step(const StateVec& x, double range_to_go,
                                   double dt) {
  ControllerOutputs out;
  t_ += dt;

  // 'Vt' (airspeed) -- not 'V' (that is the State::V enum index).
  const double ub = x[U], vb = x[V], wb = x[W];
  const double Vt = std::max(1e-3, std::sqrt(ub * ub + vb * vb + wb * wb));
  const double beta = std::asin(clampd(vb / Vt, -1.0, 1.0));
  const double h = x[H];
  const double sink = sinkRate(x);                 // positive down
  const double sink_ref = -V_app_ * std::sin(gamma_app_);  // >0 on descent

  // ---- Airspeed loop: PI -> throttle ----------------------------------------
  const double eV = V_app_ - Vt;
  V_int_ += eV * dt;
  V_int_ = clampd(V_int_, -20.0, 20.0);  // anti-windup
  double dT = u_trim_[DT] + g_.Kp_V * eV + g_.Ki_V * V_int_;
  out.V_cmd = V_app_;

  // ---- Glidepath outer loop -> theta_cmd ------------------------------------
  const double h_ref = range_to_go * std::tan(std::abs(gamma_app_));
  out.h_ref = h_ref;
  double theta_cmd;

  if (h <= flare_.h_flare && range_to_go >= 0.0) {
    // Flare: exponentially decay sink rate from entry value toward w_td.
    if (!flaring_) {
      flaring_ = true;
      flare_entry_sink_ = sink;
      t_flare_ = t_;
    }
    const double tau = std::max(1e-3, flare_.tau);
    const double w_cmd =
        flare_.w_td + (flare_entry_sink_ - flare_.w_td) *
                          std::exp(-(t_ - t_flare_) / tau);
    out.w_cmd = w_cmd;
    // Track commanded sink with a pitch command (raise nose if sinking fast).
    // Kd_gs doubles as the flare sink-error gain (sink-rate feedback).
    theta_cmd = theta_trim_ + g_.Kd_gs * (sink - w_cmd);
  } else {
    // Normal glideslope tracking: altitude error + sink-rate damping.
    const double hdot = -sink, hdot_ref = -sink_ref;
    theta_cmd = theta_trim_ - g_.Kp_gs * (h - h_ref) -
                g_.Kd_gs * (hdot - hdot_ref);
    out.w_cmd = sink_ref;
  }
  out.flaring = flaring_;
  out.theta_cmd = theta_cmd;

  // ---- Pitch inner loop -> delta_e ------------------------------------------
  const double eth = theta_cmd - x[THETA];
  theta_int_ += eth * dt;
  theta_int_ = clampd(theta_int_, -0.5, 0.5);
  double de = u_trim_[DE] + g_.Kp_theta * eth + g_.Ki_theta * theta_int_ -
              g_.Kq * x[Q];

  // ---- Lateral: cross-track -> phi_cmd --------------------------------------
  // PD (+ small I): the bank->heading->cross-track path is a double integrator,
  // so cross-track-rate (east velocity) damping Kd_y is required for stability.
  const double ct = std::cos(x[THETA]), st = std::sin(x[THETA]);
  const double cp = std::cos(x[PHI]), sp = std::sin(x[PHI]);
  const double cps = std::cos(x[PSI]), sps = std::sin(x[PSI]);
  const double ydot = ub * ct * sps + vb * (sp * st * sps + cp * cps) +
                      wb * (cp * st * sps - sp * cps);
  y_int_ += x[Y] * dt;
  y_int_ = clampd(y_int_, -50.0, 50.0);
  double phi_cmd = -(g_.Kp_y * x[Y] + g_.Kd_y * ydot + g_.Ki_y * y_int_);
  phi_cmd = clampd(phi_cmd, -g_.phi_max, g_.phi_max);
  // Near touchdown, roll out to wings level for an aligned water contact.
  if (h <= h_decrab_) phi_cmd = 0.0;
  out.phi_cmd = phi_cmd;

  // ---- Roll inner loop -> delta_a -------------------------------------------
  const double eph = phi_cmd - x[PHI];
  double da = u_trim_[DA] + g_.Kp_phi * eph - g_.Kp_p * x[P];

  // ---- Yaw damper + decrab -> delta_r ---------------------------------------
  double dr = u_trim_[DR] - g_.Kr * x[R];
  // Decrab: +delta_r yaws nose-right (+Cn); for +beta (nose left of velocity)
  // we yaw the nose right to align it, hence +Kp_beta*beta.
  if (h <= h_decrab_) dr += g_.Kp_beta * beta;

  CtrlVec u_cmd;
  u_cmd[DE] = de;
  u_cmd[DA] = da;
  u_cmd[DR] = dr;
  u_cmd[DT] = dT;

  out.u = applyLimits(u_cmd, dt);
  u_prev_ = out.u;
  have_prev_ = true;
  return out;
}

CtrlVec Controller::applyLimits(const CtrlVec& u_cmd, double dt) const {
  CtrlVec u = u_cmd;
  // Position limits.
  u[DE] = clampd(u[DE], lim_.de_min, lim_.de_max);
  u[DA] = clampd(u[DA], lim_.da_min, lim_.da_max);
  u[DR] = clampd(u[DR], lim_.dr_min, lim_.dr_max);
  u[DT] = clampd(u[DT], lim_.dT_min, lim_.dT_max);
  // Rate limits on the aerodynamic surfaces.
  if (have_prev_ && dt > 0.0) {
    const double dmax = lim_.rate * dt;
    for (int i = DE; i <= DR; ++i)
      u[i] = u_prev_[i] + clampd(u[i] - u_prev_[i], -dmax, dmax);
  }
  return u;
}

}  // namespace autoland
