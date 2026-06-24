#pragma once
#include "autoland/types.hpp"

// =============================================================================
// NominalController: the abstract nominal (guidance/inner-loop) controller the
// CBF safety filter sits downstream of.
//
// This is the SWAP SEAM for nominal controllers. The filter (and the sim) depend
// only on this interface and on ControllerOutputs, so any control law -- the
// cascaded-PID autoland (controller.hpp) or the minimal glide-slope law
// (nominal_lon_glideslope.hpp) -- can be substituted. Per the documentation
// (A.2), the nominal layer only has to put the aircraft on a reasonable
// approach; the CBF filter owns the safety-critical shaping, so it is agnostic
// to how the nominal command (delta_e^nom, T^nom) was produced.
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

class NominalController {
 public:
  virtual ~NominalController() = default;

  // Advance one control step. range_to_go is the horizontal distance to the
  // touchdown point [m]; dt is the step [s]. x is the absolute state. Returns
  // absolute virtual controls (with the law's own deflection/rate limits).
  virtual ControllerOutputs step(const StateVec& x, double range_to_go,
                                 double dt) = 0;
};

}  // namespace autoland
