#pragma once
#include <functional>
#include <string>
#include <vector>
#include "autoland/linear_model.hpp"
#include "autoland/types.hpp"

// =============================================================================
// Control Barrier Function (CBF) safety filter -- INTERFACE + PASS-THROUGH STUB.
//
// Designed so an OSQP-backed QP drops in later WITHOUT restructuring callers.
// The filter solves, at each step, the minimal deviation from the nominal
// control that keeps every barrier's forward invariance condition satisfied:
//
//     u* = argmin_u || u - u_nom ||^2_W
//          s.t.  Lf h(x) + Lg h(x) u + alpha(h(x)) >= 0     for each barrier
//                u_min <= u <= u_max
//
// With the control-affine model xdot = f(x) + g(x) u (here f = f0 + A(x-x0),
// g = B), the Lie derivatives are Lf h = grad_h . f,  Lg h = grad_h . g, so each
// barrier contributes one linear inequality in u. That is exactly the LinearModel
// this project already produces.
//
// TODO: back this with OSQP. The stub returns u_nom unchanged so the build
// compiles and the closed-loop sim runs with the filter as a pass-through.
//
// CANDIDATE BARRIERS for autoland (to implement once the QP is wired):
//   * Minimum airspeed above stall:        h = V - V_stall_margin
//   * Bank-angle limit tightening near the surface:
//                                          h = phi_max(h_agl) - |phi|,
//                                          with phi_max shrinking as h_agl -> 0
//   * Sink rate bounded as a function of height (no high sink near the water):
//                                          h = w_max(h_agl) - w   (w = sink rate)
// =============================================================================
namespace autoland {

// A scalar barrier h(x) >= 0 defines the safe set, with its state gradient.
struct Barrier {
  std::string name;
  std::function<double(const StateVec&)> h;            // h(x)
  std::function<Eigen::RowVectorXd(const StateVec&)> grad;  // dh/dx (1 x NX)
  // class-K function alpha(.) applied to h; default linear alpha(h)=k*h.
  std::function<double(double)> alpha = [](double hv) { return 1.0 * hv; };
};

class CBFFilter {
 public:
  CBFFilter() = default;

  // Filter the nominal control. model supplies the control-affine (f0,A,B) used
  // to form the Lie derivatives. Returns the safe control.
  //
  // STUB: currently returns u_nom unchanged (pass-through). The signature is
  // final; only the body changes when OSQP is added.
  CtrlVec filter(const CtrlVec& u_nom, const StateVec& x,
                 const LinearModel& model,
                 const std::vector<Barrier>& barriers) const;

  bool enabled() const { return enabled_; }
  void setEnabled(bool e) { enabled_ = e; }

 private:
  bool enabled_{true};  // when false, hard pass-through
};

}  // namespace autoland
