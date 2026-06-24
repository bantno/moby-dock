#pragma once
#include "autoland/dynamics.hpp"
#include "autoland/types.hpp"

// =============================================================================
// Trim solver. Finds (alpha, delta_e, throttle) for steady flight at a target
// airspeed V_app and flight-path angle gamma_app, wings level, zero sideslip.
//
// It is a Newton solve directly on the nonlinear EOM (Dynamics::xdot) -- the
// same single source of truth used by the linear model. The three residuals
// are the body-axis x-force, z-force, and pitching-moment accelerations; the
// lateral-directional accelerations are identically zero by symmetry.
//
// Trim is a steady DESCENT: hdot = V*sin(gamma) != 0, which is expected and is
// not driven to zero.
// =============================================================================
namespace autoland {

struct TrimResult {
  StateVec x{StateVec::Zero()};  // trimmed state (h=y=0 reference)
  CtrlVec u{CtrlVec::Zero()};    // trimmed virtual controls
  double alpha{0};               // [rad]
  double theta{0};               // [rad]
  double V{0};                   // [m/s]
  double gamma{0};               // [rad]
  double residual{0};            // L2 norm of the trim residual
  int iterations{0};
  bool converged{false};
};

// Solve for trim. tol is on the residual L2 norm.
TrimResult trim(const Dynamics& dyn, double V_app, double gamma_app,
                double tol = 1e-9, int max_iter = 100);

}  // namespace autoland
