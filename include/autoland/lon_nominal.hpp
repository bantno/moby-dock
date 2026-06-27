#pragma once
#include <algorithm>
#include "autoland/lon_augmented.hpp"

// =============================================================================
// Nominal controller for the augmented longitudinal system.
//
// Per user direction (replacing the doc's TECS): a CONSTANT THRUST setpoint with
// a PD law on the augmented thrust state, plus a CASCADE elevator law that holds
// a constant negative flight-path angle gamma:
//   * Thrust (doc section 2 PD):  Tddot_nom = Kp_T (T_set - T) - Kd_T Tdot
//   * Elevator (cascade gamma -> theta -> pitch):
//       outer PI:  theta_cmd = theta_trim + Kp_g (gamma_ref - gamma) + Ki_g * int
//       inner PD:  de_nom    = Kp_th (theta_cmd - theta) - Kq q
// The CBF-QP filter then minimally adjusts [de_nom, Tddot_nom].
// =============================================================================
namespace autoland {

struct LonNominalConfig {
  double T_set{6.0};         // constant thrust setpoint [N]
  double gamma_ref{-0.03};   // target flight-path angle [rad] (negative=descent)
  double theta_trim{0.0};    // pitch feedforward [rad]
  double Kp_T{4.0}, Kd_T{4.0};            // thrust PD
  double Kp_gamma{2.0}, Ki_gamma{0.5};    // gamma -> theta_cmd (PI)
  double theta_cmd_max{0.30};             // |theta_cmd - theta_trim| limit [rad]
  double Kp_theta{6.0}, Kq{1.5};          // theta -> de (PD)
};

class LonNominal {
 public:
  explicit LonNominal(const LonNominalConfig& cfg) : c_(cfg) {}

  LonCtrlVec step(const LonStateVec& X, double dt) {
    LonCtrlVec U;

    // Thrust channel: PD tracking of the constant setpoint by the T state.
    U[LTDDOT] = c_.Kp_T * (c_.T_set - X[LT]) - c_.Kd_T * X[LTDOT];

    // Flight-path-angle outer loop -> attitude command (PI, anti-windup clamp).
    const double e_gamma = c_.gamma_ref - X[LGAM];
    gamma_int_ += e_gamma * dt;
    double theta_cmd = c_.theta_trim + c_.Kp_gamma * e_gamma + c_.Ki_gamma * gamma_int_;
    const double lo = c_.theta_trim - c_.theta_cmd_max;
    const double hi = c_.theta_trim + c_.theta_cmd_max;
    if (theta_cmd > hi) { theta_cmd = hi; gamma_int_ -= e_gamma * dt; }  // de-windup
    else if (theta_cmd < lo) { theta_cmd = lo; gamma_int_ -= e_gamma * dt; }

    // Attitude inner loop -> elevator (PD with pitch-rate damping).
    U[LDE] = c_.Kp_theta * (theta_cmd - X[LTH]) - c_.Kq * X[LQ];

    theta_cmd_ = theta_cmd;
    return U;
  }

  double thetaCmd() const { return theta_cmd_; }
  void reset() { gamma_int_ = 0.0; }

 private:
  LonNominalConfig c_;
  double gamma_int_{0.0};
  double theta_cmd_{0.0};
};

}  // namespace autoland
