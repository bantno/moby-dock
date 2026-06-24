#pragma once
#include "autoland/config.hpp"
#include "autoland/nominal_controller.hpp"
#include "autoland/trim.hpp"
#include "autoland/types.hpp"

// =============================================================================
// LonGlideslopeController: the minimal, platform-agnostic nominal of doc A.2.
//
// Commands a constant approach flight-path angle (a fixed glide slope) tracked
// by a simple longitudinal law: PI airspeed loop -> throttle, glideslope error
// + sink-rate damping -> theta_cmd, pitch loop -> delta_e. There is NO flare and
// NO lateral control here -- the lateral controls are held at trim and the flare
// is intended to EMERGE from the descent-rate CBF (once promoted to an HOCBF).
//
// It exists to demonstrate the NominalController swap seam: drop it in for the
// cascaded-PID Controller with no change to the CBF filter or the sim loop.
// =============================================================================
namespace autoland {

class LonGlideslopeController : public NominalController {
 public:
  LonGlideslopeController(const Gains& gains, const TrimResult& trim,
                          double V_app, double gamma_app,
                          const SurfaceLimits& limits);

  ControllerOutputs step(const StateVec& x, double range_to_go,
                         double dt) override;

 private:
  Gains g_;
  double V_app_, gamma_app_;
  SurfaceLimits lim_;
  CtrlVec u_trim_;
  double theta_trim_;

  double V_int_{0};
  double theta_int_{0};
  CtrlVec u_prev_;
  bool have_prev_{false};

  CtrlVec applyLimits(const CtrlVec& u_cmd, double dt) const;
};

}  // namespace autoland
