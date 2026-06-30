#pragma once
#include <algorithm>
#include "autoland/lon_augmented.hpp"

// =============================================================================
// Nominal controller for the augmented longitudinal system.
//
// Per user direction (replacing the doc's TECS): a CONSTANT THRUST setpoint with
// a PD law on the augmented thrust state, plus a CASCADE elevator law that holds
// a PIECEWISE flight-path angle gamma:
//   * Thrust (doc section 2 PD):  Tddot_nom = Kp_T (T_set - T) - Kd_T Tdot
//   * Elevator (cascade gamma -> theta -> pitch):
//       outer PI:  theta_cmd = theta_trim + Kp_g (gamma_ref(h) - gamma) + Ki_g int
//       inner PD:  de_nom    = Kp_th (theta_cmd - theta) - Kq q
// The flight-path reference is a two-segment glideslope: gamma_ref ABOVE the switch
// altitude h_switch, then gamma_ref2 below it to ground contact (e.g. a shallow
// approach that steepens, or vice versa). The integrator carries across the switch
// so the inner loop slews smoothly onto the new segment. The CBF-QP filter then
// minimally adjusts [de_nom, Tddot_nom].
// =============================================================================
namespace autoland {

struct LonNominalConfig {
  double T_set{6.0};         // constant thrust setpoint [N]
  double gamma_ref{-0.03};   // segment-1 flight-path angle [rad] (above h_switch)
  double theta_trim{0.0};    // pitch feedforward [rad]
  double Kp_T{4.0}, Kd_T{4.0};            // thrust PD
  double Kp_gamma{2.0}, Ki_gamma{0.5};    // gamma -> theta_cmd (PI)
  double theta_cmd_max{0.30};             // |theta_cmd - theta_trim| limit [rad]
  double Kp_theta{6.0}, Kq{1.5};          // theta -> de (PD)
  // Piecewise glideslope. Below h_switch the reference switches to gamma_ref2 and
  // holds it to ground contact. gamma_ref2 = gamma_ref (or h_switch <= 0) => single
  // segment (original constant-gamma behavior).
  double gamma_ref2{-0.03};  // segment-2 flight-path angle [rad] (below h_switch)
  double h_switch{0.0};      // glideslope transition altitude [m]
  // Flare / hold-off. Below h_flare the flight-path reference ramps from the
  // glideslope toward gamma_flare (a level / slightly climbing hold-off) as h -> 0,
  // so the cascade pitches the nose up into a high-alpha float that arrests sink AND
  // bleeds airspeed before contact -- using the otherwise-idle elevator authority.
  // h_flare <= 0 => no flare (straight glideslope to contact, original behavior).
  double h_flare{0.0};       // flare engage height AGL [m]
  double gamma_flare{0.0};   // flight-path target at touchdown during flare [rad]
};

class LonNominal {
 public:
  explicit LonNominal(const LonNominalConfig& cfg) : c_(cfg) {}

  LonCtrlVec step(const LonStateVec& X, double dt) {
    LonCtrlVec U;

    // Thrust channel: PD tracking of the constant setpoint by the T state.
    U[LTDDOT] = c_.Kp_T * (c_.T_set - X[LT]) - c_.Kd_T * X[LTDOT];

    // Piecewise glideslope: segment-1 angle aloft, segment-2 below the switch.
    double gamma_ref = (X[LH] > c_.h_switch) ? c_.gamma_ref : c_.gamma_ref2;
    // Flare / hold-off: below h_flare, ramp the reference from the glideslope toward
    // gamma_flare as h -> 0 (s = h/h_flare goes 1 -> 0). The PI then drives gamma up,
    // pitching to a high-alpha float that washes off airspeed before touchdown.
    if (c_.h_flare > 0.0 && X[LH] < c_.h_flare) {
      const double s = std::max(0.0, std::min(1.0, X[LH] / c_.h_flare));
      gamma_ref = c_.gamma_flare + (gamma_ref - c_.gamma_flare) * s;
    }
    gamma_ref_ = gamma_ref;

    // Flight-path-angle outer loop -> attitude command (PI, anti-windup clamp).
    const double e_gamma = gamma_ref - X[LGAM];
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
  double gammaRef() const { return gamma_ref_; }  // active (piecewise) reference
  void reset() { gamma_int_ = 0.0; }

 private:
  LonNominalConfig c_;
  double gamma_int_{0.0};
  double theta_cmd_{0.0};
  double gamma_ref_{0.0};
};

}  // namespace autoland
