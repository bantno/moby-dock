#pragma once
#include <array>
#include "autoland/dynamics.hpp"
#include "autoland/types.hpp"

// =============================================================================
// Linear model: a numerical linearization of the nonlinear EOM about a trim
// point, by CENTRAL DIFFERENCES on Dynamics::xdot. The model is the affine
// control-affine form
//
//        xdot ~= f0 + A (x - x0) + B (u - u0)
//
// where f0 = xdot(x0,u0) is the trim drift (nonzero only for the guidance
// states: steady descent hdot and any cross-track rate). Writing it this way
// keeps xdot = A x + B u explicit and control-affine, which is exactly the form
// the CBF-QP filter consumes later (B is the control matrix g(x)).
//
// -----------------------------------------------------------------------------
// SIGN CONVENTIONS (body axes; right-handed, x fwd / y right / z down)
// -----------------------------------------------------------------------------
//   States (deviations from trim):
//     u,v,w   body-axis velocities         +u forward, +v right, +w down
//     p,q,r   body-axis angular rates       +p right-roll, +q nose-up,
//                                            +r nose-right (yaw)
//     phi     roll angle                    +right wing down
//     theta   pitch angle                   +nose up
//     psi     yaw angle                     +nose right of North
//     h       altitude                      +up      (hdot < 0 on descent)
//     y       cross-track position          +east of the approach centerline
//
//   Controls (virtual; deviations from trim). Polarity convention: a POSITIVE
//   command produces a POSITIVE body-axis response (enforced by the mixing map):
//     delta_e elevator   +delta_e -> nose-UP pitch     (Cm wrt delta_e > 0)
//     delta_a aileron    +delta_a -> right roll         (Cl wrt delta_a > 0)
//     delta_r rudder     +delta_r -> nose-right yaw     (Cn wrt delta_r > 0)
//     delta_T throttle   +increases thrust
//   (The .stab's Elevator/Aileron groups give a negative moment per positive
//   group deflection, so the default mixing applies -1 to those; see mixing.cpp.)
//
//   Aero forces/moments are built from BODY-AXIS coefficients
//   (CFx,CFy,CFz,Cl=CMx,Cm=CMy,Cn=CMz), never from wind-axis CL/CD. The .stab
//   reports them in a VSPAero frame with Z UP; dynamics.cpp applies a documented
//   per-coefficient sign transform to standard z-DOWN axes. Rate derivatives use
//   VSPAero's nondimensional rates (phat=pb/2V, qhat=qc/2V, rhat=rb/2V).
//
// -----------------------------------------------------------------------------
// DECOUPLED SUB-MODELS (exact index selections of the full A,B)
//   Longitudinal: state [u, w, q, theta, h]      control [delta_e, delta_T]
//   Lateral-dir.: state [v, p, r, phi, psi, y]    control [delta_a, delta_r]
// Altitude h and cross-track y are integrated guidance states.
// =============================================================================
namespace autoland {

// Index sets for the decoupled sub-models (into the full NX/NU vectors).
inline constexpr std::array<int, 5> kLonStates{U, W, Q, THETA, H};
inline constexpr std::array<int, 2> kLonCtrls{DE, DT};
inline constexpr std::array<int, 6> kLatStates{V, P, R, PHI, PSI, Y};
inline constexpr std::array<int, 2> kLatCtrls{DA, DR};

struct LinearModel {
  StateVec x0{StateVec::Zero()};  // linearization (trim) state
  CtrlVec u0{CtrlVec::Zero()};    // linearization (trim) control
  StateVec f0{StateVec::Zero()};  // xdot(x0,u0): trim drift

  Mat A;  // NX x NX
  Mat B;  // NX x NU

  Mat A_lon, B_lon;  // 5x5, 5x2
  Mat A_lat, B_lat;  // 6x6, 6x2

  // Reconstruct the affine plant derivative at an absolute (x,u).
  StateVec xdot(const StateVec& x, const CtrlVec& u) const {
    return f0 + A * (x - x0) + B * (u - u0);
  }
};

// Linearize the dynamics about (x0,u0) by central differences.
LinearModel linearize(const Dynamics& dyn, const StateVec& x0,
                      const CtrlVec& u0);

}  // namespace autoland
