#pragma once
#include "autoland/config.hpp"
#include "autoland/trim.hpp"
#include "autoland/types.hpp"

// =============================================================================
// Cascaded-PID autoland controller, FRONTSIDE technique.
//
//   Airspeed loop : PI on V error -> throttle.            (throttle holds speed)
//   Glidepath     : altitude/glideslope error -> theta_cmd.
//   Pitch inner   : theta_cmd error + pitch-rate damping -> delta_e.
//                                                          (pitch holds attitude)
//   Flare         : below h_flare, the outer reference switches to an
//                   exponentially decaying sink rate blending to w_td.
//   Lateral       : cross-track -> phi_cmd ; roll inner (+ roll-rate damping)
//                   -> delta_a ; yaw damper on r -> delta_r ; decrab drives
//                   sideslip to zero below h_decrab for an aligned water contact.
//
// Full-state feedback from the (linear) plant -- no estimator/sensor model yet.
// All gains come from config. Outputs are ABSOLUTE virtual controls with
// deflection and rate limits applied.
// =============================================================================
namespace autoland {

struct ControllerOutputs {
  CtrlVec u{CtrlVec::Zero()};  // absolute virtual controls [de, da, dr, dT]
  // References / internals, exposed for logging:
  double V_cmd{0};
  double h_ref{0};
  double theta_cmd{0};
  double phi_cmd{0};
  double w_cmd{0};   // commanded sink rate (flare), positive down
  bool flaring{false};
};

class Controller {
 public:
  Controller(const Gains& gains, const FlareParams& flare,
             const TrimResult& trim, double V_app, double gamma_app,
             double h_decrab, const Environment& env,
             const SurfaceLimits& limits);

  // Advance one control step. range_to_go is the horizontal distance to the
  // touchdown point [m]; dt is the step [s]. x is the absolute state.
  ControllerOutputs step(const StateVec& x, double range_to_go, double dt);

 private:
  Gains g_;
  FlareParams flare_;
  double V_app_, gamma_app_, h_decrab_, g_accel_;
  SurfaceLimits lim_;
  StateVec x_trim_;
  CtrlVec u_trim_;
  double theta_trim_;

  // Integrator / mode state.
  double V_int_{0};
  double theta_int_{0};
  double y_int_{0};
  double t_{0};
  bool flaring_{false};
  double flare_entry_sink_{0};
  double t_flare_{0};
  CtrlVec u_prev_;
  bool have_prev_{false};

  CtrlVec applyLimits(const CtrlVec& u_cmd, double dt) const;
};

}  // namespace autoland
