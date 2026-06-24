#include "autoland/cbf.hpp"

namespace autoland {

CtrlVec CBFFilter::filter(const CtrlVec& u_nom, const StateVec& x,
                          const LinearModel& model,
                          const std::vector<Barrier>& barriers) const {
  (void)x;
  (void)model;
  (void)barriers;

  // ---------------------------------------------------------------------------
  // TODO(OSQP): assemble and solve the safety QP here.
  //
  //   minimize    (u - u_nom)^T W (u - u_nom)
  //   subject to  for each barrier b:
  //                   grad_b . (f0 + A (x - x0))            // Lf h
  //                 + grad_b . B  * u                       // Lg h * u
  //                 + b.alpha(b.h(x))  >= 0
  //               u_min <= u <= u_max
  //
  // Each barrier yields one row  G u <= k  with
  //     G_row = -(grad_b * B),  k_row = grad_b * (f0 + A(x-x0)) + alpha(h).
  // Box constraints come from SurfaceLimits. OSQP solves the resulting QP.
  // Until then we pass the nominal control through unchanged.
  // ---------------------------------------------------------------------------
  if (!enabled_) return u_nom;
  return u_nom;  // pass-through stub
}

}  // namespace autoland
