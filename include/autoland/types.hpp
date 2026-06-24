#pragma once
#include <Eigen/Dense>

// =============================================================================
// Common types and the canonical state / control layout for the whole project.
//
// AXES: body axes throughout (x forward, y right, z down). All angles in
// radians, all lengths in meters, time in seconds, forces in newtons, moments
// in newton-metres. See linear_model.hpp for the full sign-convention block.
// =============================================================================
namespace autoland {

using Vec = Eigen::VectorXd;
using Mat = Eigen::MatrixXd;

// ---- Full nonlinear state (11 states) ---------------------------------------
// Index layout chosen so the longitudinal and lateral-directional sub-models
// (see linear_model) are simple index selections of this one vector.
//
//   u, v, w : body-axis velocities                       [m/s]
//   p, q, r : body-axis angular rates                    [rad/s]
//   phi, theta, psi : Euler angles (roll, pitch, yaw)    [rad]
//   h       : altitude (positive up), integrated guidance state   [m]
//   y       : cross-track position (east of centerline), guidance state [m]
//
// Downrange distance is NOT a state (it does not feed back into the dynamics);
// the sim tracks range-to-touchdown separately for glideslope referencing.
constexpr int NX = 11;
enum State : int {
  U = 0, V = 1, W = 2,
  P = 3, Q = 4, R = 5,
  PHI = 6, THETA = 7, PSI = 8,
  H = 9, Y = 10
};

// ---- Virtual control vector (4 controls) ------------------------------------
//   delta_e : virtual elevator (pitch)     [rad]   -> mixed onto control groups
//   delta_a : virtual aileron  (roll)      [rad]
//   delta_r : virtual rudder   (yaw)       [rad]
//   delta_T : throttle command             [-]  (0..1 nominal)
//
// The mixing matrix (config) maps [delta_e, delta_a, delta_r] onto the physical
// VSPAero control groups; delta_T drives the thrust model.
constexpr int NU = 4;
enum Ctrl : int { DE = 0, DA = 1, DR = 2, DT = 3 };

using StateVec = Eigen::Matrix<double, NX, 1>;
using CtrlVec  = Eigen::Matrix<double, NU, 1>;

// ---- Aerodynamic body-axis coefficient set ----------------------------------
// Order matches the .stab columns we consume:
//   CFx, CFy, CFz : body-axis force coefficients
//   Cl=CMx, Cm=CMy, Cn=CMz : body-axis moment coefficients (roll/pitch/yaw)
enum Coef : int { CFX = 0, CFY = 1, CFZ = 2, CMX = 3, CMY = 4, CMZ = 5 };
constexpr int NCOEF = 6;
using CoefVec = Eigen::Matrix<double, NCOEF, 1>;

}  // namespace autoland
