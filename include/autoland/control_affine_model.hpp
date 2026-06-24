#pragma once
#include "autoland/linear_model.hpp"
#include "autoland/types.hpp"

// =============================================================================
// ControlAffineModel: the prediction model the CBF-QP consumes.
//
// A control-affine model writes the state derivative as
//
//        xdot = f(x) + g(x) u
//
// where f(x) is the drift (NX vector) and g(x) is the control matrix
// (NX x NU). The CBF filter forms each barrier's Lie derivatives from these:
//   Lf h(x) = grad_h(x) . f(x),   Lg h(x) = grad_h(x) . g(x).
//
// This is the SWAP SEAM for dynamics models. The filter depends only on this
// interface, so a reduced longitudinal model, the linearization, or a fully
// nonlinear control-affine model can be substituted with no change to cbf.cpp.
// The first/only concrete model here is LinearModelCAM (the longitudinal-capable
// linearization the project already produces); a nonlinear model can be added
// later as another ControlAffineModel.
// =============================================================================
namespace autoland {

class ControlAffineModel {
 public:
  virtual ~ControlAffineModel() = default;

  // Drift f(x): the state derivative with the control set to zero.
  virtual StateVec drift(const StateVec& x) const = 0;

  // Control matrix g(x): columns are d(xdot)/d(u_j), so xdot = f(x) + g(x) u.
  virtual Mat ctrlMatrix(const StateVec& x) const = 0;
};

// Adapter over the affine LinearModel  xdot = f0 + A (x - x0) + B (u - u0).
// Rewritten in pure control-affine form:
//   f(x) = f0 + A (x - x0) - B u0,   g(x) = B   (so f(x) + B u = the above).
class LinearModelCAM : public ControlAffineModel {
 public:
  explicit LinearModelCAM(const LinearModel& m) : m_(m) {}

  StateVec drift(const StateVec& x) const override {
    return m_.f0 + m_.A * (x - m_.x0) - m_.B * m_.u0;
  }
  Mat ctrlMatrix(const StateVec& /*x*/) const override { return m_.B; }

 private:
  const LinearModel& m_;
};

}  // namespace autoland
