#include "autoland/nominal_lon_glideslope.hpp"

#include <algorithm>
#include <cmath>

namespace autoland {
namespace {

double clampd(double v, double lo, double hi) {
  return std::max(lo, std::min(hi, v));
}

// Sink rate (positive down), same kinematics as Dynamics / Controller.
double sinkRate(const StateVec& x) {
  const double ct = std::cos(x[THETA]), st = std::sin(x[THETA]);
  const double cp = std::cos(x[PHI]), sp = std::sin(x[PHI]);
  const double hdot = x[U] * st - x[V] * sp * ct - x[W] * cp * ct;
  return -hdot;
}

}  // namespace

LonGlideslopeController::LonGlideslopeController(const Gains& gains,
                                                const TrimResult& trim,
                                                double V_app, double gamma_app,
                                                const SurfaceLimits& limits)
    : g_(gains),
      V_app_(V_app),
      gamma_app_(gamma_app),
      lim_(limits),
      u_trim_(trim.u),
      theta_trim_(trim.theta) {}

ControllerOutputs LonGlideslopeController::step(const StateVec& x,
                                               double range_to_go, double dt) {
  ControllerOutputs out;

  const double ub = x[U], vb = x[V], wb = x[W];
  const double Vt = std::max(1e-3, std::sqrt(ub * ub + vb * vb + wb * wb));
  const double h = x[H];
  const double sink = sinkRate(x);                          // positive down
  const double sink_ref = -V_app_ * std::sin(gamma_app_);   // >0 on descent

  // ---- Airspeed loop: PI -> throttle ----------------------------------------
  const double eV = V_app_ - Vt;
  V_int_ += eV * dt;
  V_int_ = clampd(V_int_, -20.0, 20.0);
  const double dT = u_trim_[DT] + g_.Kp_V * eV + g_.Ki_V * V_int_;
  out.V_cmd = V_app_;

  // ---- Glideslope outer loop -> theta_cmd (NO flare) ------------------------
  const double h_ref = range_to_go * std::tan(std::abs(gamma_app_));
  out.h_ref = h_ref;
  const double hdot = -sink, hdot_ref = -sink_ref;
  const double theta_cmd =
      theta_trim_ - g_.Kp_gs * (h - h_ref) - g_.Kd_gs * (hdot - hdot_ref);
  out.theta_cmd = theta_cmd;
  out.w_cmd = sink_ref;
  out.flaring = false;

  // ---- Pitch inner loop -> delta_e ------------------------------------------
  const double eth = theta_cmd - x[THETA];
  theta_int_ += eth * dt;
  theta_int_ = clampd(theta_int_, -0.5, 0.5);
  const double de = u_trim_[DE] + g_.Kp_theta * eth + g_.Ki_theta * theta_int_ -
                    g_.Kq * x[Q];

  // ---- Lateral controls held at trim (longitudinal-only nominal) ------------
  CtrlVec u_cmd;
  u_cmd[DE] = de;
  u_cmd[DA] = u_trim_[DA];
  u_cmd[DR] = u_trim_[DR];
  u_cmd[DT] = dT;

  out.u = applyLimits(u_cmd, dt);
  out.phi_cmd = 0.0;
  u_prev_ = out.u;
  have_prev_ = true;
  return out;
}

CtrlVec LonGlideslopeController::applyLimits(const CtrlVec& u_cmd,
                                            double dt) const {
  CtrlVec u = u_cmd;
  u[DE] = clampd(u[DE], lim_.de_min, lim_.de_max);
  u[DA] = clampd(u[DA], lim_.da_min, lim_.da_max);
  u[DR] = clampd(u[DR], lim_.dr_min, lim_.dr_max);
  u[DT] = clampd(u[DT], lim_.dT_min, lim_.dT_max);
  if (have_prev_ && dt > 0.0) {
    const double dmax = lim_.rate * dt;
    for (int i = DE; i <= DR; ++i)
      u[i] = u_prev_[i] + clampd(u[i] - u_prev_[i], -dmax, dmax);
  }
  return u;
}

}  // namespace autoland
