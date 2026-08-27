#include "autoland/px4_tecs.hpp"

#include <algorithm>

// Port of PX4-Autopilot src/lib/tecs/TECS.cpp (TECSControl methods only),
// main @ a906b728. Each method mirrors the upstream one of the same name;
// comments marked [PX4] are upstream's own. Where upstream evaluates a
// fast_descend term, it is written out at fast_descend = 0 (mode not ported).
namespace autoland {
namespace px4 {

namespace {
inline double constrain(double v, double lo, double hi) { return std::clamp(v, lo, hi); }
}  // namespace

void TecsControl::initialize(const TecsSetpoint& setpoint, const TecsInput& input,
                             const TecsParam& param, const TecsFlag& flag) {
  resetIntegrals();

  AltitudePitchControl control_setpoint;
  control_setpoint.tas_setpoint = setpoint.tas_setpoint;
  control_setpoint.tas_rate_setpoint = calcAirspeedControlOutput(setpoint, input, param, flag);
  control_setpoint.altitude_rate_setpoint = calcAltitudeControlOutput(setpoint, input, param);

  const SpecificEnergyRates specific_energy_rate{calcSpecificEnergyRates(control_setpoint, input)};

  detectUnderspeed(input, param, flag);

  const SpecificEnergyWeighting weight{updateSpeedAltitudeWeights(param, flag)};
  const ControlValues seb_rate{calcPitchControlSebRate(weight, specific_energy_rate)};

  pitch_setpoint_ = calcPitchControlOutput(input, seb_rate, param, flag);

  const STERateLimit limit{calculateTotalEnergyRateLimit(param)};

  // _ste_rate_estimate_filter.reset(...)
  ste_rate_estimate_filter_ =
      specific_energy_rate.spe_rate.estimate + specific_energy_rate.ske_rate.estimate;

  const ControlValues ste_rate{calcThrottleControlSteRate(limit, specific_energy_rate, param)};

  throttle_setpoint_ = calcThrottleControlOutput(limit, ste_rate, param, flag);

  // Debug output
  debug_output_.total_energy_rate_estimate = ste_rate.estimate;
  debug_output_.total_energy_rate_sp = ste_rate.setpoint;
  debug_output_.throttle_integrator = throttle_integ_state_;
  debug_output_.energy_balance_rate_estimate = seb_rate.estimate;
  debug_output_.energy_balance_rate_sp = seb_rate.setpoint;
  debug_output_.pitch_integrator = pitch_integ_state_;
  debug_output_.altitude_rate_control = control_setpoint.altitude_rate_setpoint;
  debug_output_.true_airspeed_derivative_control = control_setpoint.tas_rate_setpoint;
}

void TecsControl::update(double dt, const TecsSetpoint& setpoint, const TecsInput& input,
                         const TecsParam& param, const TecsFlag& flag) {
  // Input checking (TIMESTAMP_VALID): do not update the states and output.
  if (!(std::isfinite(dt) && dt > FLT_EPSILON)) return;

  AltitudePitchControl control_setpoint;
  control_setpoint.tas_setpoint = setpoint.tas_setpoint;
  control_setpoint.tas_rate_setpoint = calcAirspeedControlOutput(setpoint, input, param, flag);

  if (std::isfinite(setpoint.altitude_rate_setpoint_direct)) {
    // [PX4] direct height rate control
    control_setpoint.altitude_rate_setpoint = setpoint.altitude_rate_setpoint_direct;
  } else {
    // [PX4] altitude is locked, go through altitude outer loop
    control_setpoint.altitude_rate_setpoint = calcAltitudeControlOutput(setpoint, input, param);
  }

  const SpecificEnergyRates specific_energy_rate{calcSpecificEnergyRates(control_setpoint, input)};

  detectUnderspeed(input, param, flag);

  calcPitchControl(dt, input, specific_energy_rate, param, flag);

  calcThrottleControl(dt, specific_energy_rate, param, flag);

  debug_output_.altitude_rate_control = control_setpoint.altitude_rate_setpoint;
  debug_output_.true_airspeed_derivative_control = control_setpoint.tas_rate_setpoint;
  debug_output_.pitch_integrator = pitch_integ_state_;
  debug_output_.throttle_integrator = throttle_integ_state_;
}

TecsControl::STERateLimit TecsControl::calculateTotalEnergyRateLimit(const TecsParam& param) const {
  STERateLimit limit;
  // [PX4] Calculate the specific total energy rate limits from the max throttle limits
  limit.STE_rate_max = std::max(param.max_climb_rate, double(FLT_EPSILON)) * kOneG;
  limit.STE_rate_min = -std::max(param.min_sink_rate, double(FLT_EPSILON)) * kOneG;
  return limit;
}

double TecsControl::calcAirspeedControlOutput(const TecsSetpoint& setpoint, const TecsInput& input,
                                              const TecsParam& param, const TecsFlag& flag) const {
  double airspeed_rate_output{0.0};

  const STERateLimit limit{calculateTotalEnergyRateLimit(param)};

  // [PX4] calculate the demanded true airspeed rate of change based on first order response of
  // true airspeed error. If airspeed measurement is not enabled then always set the rate
  // setpoint to zero in order to avoid constant rate setpoints
  if (flag.airspeed_enabled) {
    // [PX4] Calculate limits for the demanded rate of change of speed based on physical
    // performance limits with a 50% margin to allow the total energy controller to correct
    // for errors. (Upstream: (fast_descend * 0.5 + 0.5) * limit / tas, at fast_descend = 0.)
    const double max_tas_rate_sp = 0.5 * limit.STE_rate_max / std::max(input.tas, double(FLT_EPSILON));
    const double min_tas_rate_sp = 0.5 * limit.STE_rate_min / std::max(input.tas, double(FLT_EPSILON));
    airspeed_rate_output = constrain((setpoint.tas_setpoint - input.tas) * param.airspeed_error_gain,
                                     min_tas_rate_sp, max_tas_rate_sp);
  }

  return airspeed_rate_output;
}

double TecsControl::calcAltitudeControlOutput(const TecsSetpoint& setpoint, const TecsInput& input,
                                              const TecsParam& param) const {
  double altitude_rate_output =
      (setpoint.altitude_reference.alt - input.altitude) * param.altitude_error_gain +
      param.altitude_setpoint_gain_ff * setpoint.altitude_reference.alt_rate;

  altitude_rate_output = constrain(altitude_rate_output, -param.max_sink_rate, param.max_climb_rate);

  return altitude_rate_output;
}

TecsControl::SpecificEnergyRates TecsControl::calcSpecificEnergyRates(
    const AltitudePitchControl& control_setpoint, const TecsInput& input) const {
  SpecificEnergyRates specific_energy_rates;
  // [PX4] Calculate specific energy rate demands in units of (m**2/sec**3)
  specific_energy_rates.spe_rate.setpoint =
      control_setpoint.altitude_rate_setpoint * kOneG;  // potential energy rate of change
  specific_energy_rates.ske_rate.setpoint =
      control_setpoint.tas_setpoint * control_setpoint.tas_rate_setpoint;  // kinetic energy rate of change

  // [PX4] Calculate specific energy rates in units of (m**2/sec**3)
  specific_energy_rates.spe_rate.estimate = input.altitude_rate * kOneG;   // potential energy rate of change
  specific_energy_rates.ske_rate.estimate = input.tas * input.tas_rate;     // kinetic energy rate of change

  return specific_energy_rates;
}

void TecsControl::detectUnderspeed(const TecsInput& input, const TecsParam& param, const TecsFlag& flag) {
  if (!flag.detect_underspeed_enabled || !flag.airspeed_enabled) {
    ratio_undersped_ = 0.0;
    return;
  }

  // [PX4] this is the expected (something like standard) deviation from the airspeed setpoint
  // that we allow the airspeed to vary in before ramping in underspeed mitigation
  const double tas_error_bound = param.tas_error_percentage * param.equivalent_airspeed_trim;

  // [PX4] this is the soft boundary where underspeed mitigation is ramped in
  // NOTE: it's currently the same as the error bound, but separated here to indicate these
  // values do not in general need to be the same
  const double tas_underspeed_soft_bound = param.tas_error_percentage * param.equivalent_airspeed_trim;

  const double tas_fully_undersped =
      std::max(param.tas_min - tas_error_bound - tas_underspeed_soft_bound, 0.0);
  const double tas_starting_to_underspeed = std::max(param.tas_min - tas_error_bound, tas_fully_undersped);

  ratio_undersped_ =
      1.0 - constrain((input.tas - tas_fully_undersped) /
                          std::max(tas_starting_to_underspeed - tas_fully_undersped, double(FLT_EPSILON)),
                      0.0, 1.0);
}

TecsControl::SpecificEnergyWeighting TecsControl::updateSpeedAltitudeWeights(const TecsParam& param,
                                                                             const TecsFlag& flag) const {
  SpecificEnergyWeighting weight;
  // [PX4] Calculate the weight applied to control of specific kinetic energy error
  double pitch_speed_weight = constrain(param.pitch_speed_weight, 0.0, 2.0);

  // [PX4] Underspeed or fast descend: interpolate speed weight from nominal towards max.
  // (fast_descend = 0 here, so max_ratio is the underspeed ratio alone.)
  const double max_ratio = ratio_undersped_;
  pitch_speed_weight = max_ratio * 2.0 + (1.0 - max_ratio) * pitch_speed_weight;

  if (!flag.airspeed_enabled) pitch_speed_weight = 0.0;

  weight.spe_weighting = constrain(2.0 - pitch_speed_weight, 0.0, 2.0);
  weight.ske_weighting = constrain(pitch_speed_weight, 0.0, 2.0);

  return weight;
}

void TecsControl::calcPitchControl(double dt, const TecsInput& input,
                                   const SpecificEnergyRates& specific_energy_rates,
                                   const TecsParam& param, const TecsFlag& flag) {
  const SpecificEnergyWeighting weight{updateSpeedAltitudeWeights(param, flag)};
  const ControlValues seb_rate{calcPitchControlSebRate(weight, specific_energy_rates)};

  calcPitchControlUpdate(dt, input, seb_rate, param);
  const double pitch_setpoint{calcPitchControlOutput(input, seb_rate, param, flag)};

  // [PX4] Comply with the specified vertical acceleration limit by applying a pitch rate limit
  // NOTE: at zero airspeed, the pitch increment is unbounded
  const double pitch_increment = dt * param.vert_accel_limit / std::max(input.tas, double(FLT_EPSILON));
  pitch_setpoint_ = constrain(pitch_setpoint, pitch_setpoint_ - pitch_increment,
                              pitch_setpoint_ + pitch_increment);
  pitch_setpoint_ = constrain(pitch_setpoint_, param.pitch_min, param.pitch_max);

  // Debug Output
  debug_output_.energy_balance_rate_estimate = seb_rate.estimate;
  debug_output_.energy_balance_rate_sp = seb_rate.setpoint;
  debug_output_.pitch_integrator = pitch_integ_state_;
}

TecsControl::ControlValues TecsControl::calcPitchControlSebRate(
    const SpecificEnergyWeighting& weight, const SpecificEnergyRates& specific_energy_rates) const {
  ControlValues seb_rate;
  /* [PX4]
   * The SKE_weighting variable controls how speed and altitude control are prioritized by the
   * pitch demand calculation.
   * A weighting of 1 gives equal speed and altitude priority
   * A weighting of 0 gives 100% priority to altitude control and must be used when no airspeed
   * measurement is available.
   * A weighting of 2 provides 100% priority to speed control and is used when:
   * a) an underspeed condition is detected.
   * b) during climbout where a minimum pitch angle has been set to ensure altitude is gained.
   * The weighting can be adjusted between 0 and 2 depending on speed and altitude accuracy
   * requirements.
   */
  seb_rate.setpoint = specific_energy_rates.spe_rate.setpoint * weight.spe_weighting -
                      specific_energy_rates.ske_rate.setpoint * weight.ske_weighting;

  seb_rate.estimate = (specific_energy_rates.spe_rate.estimate * weight.spe_weighting) -
                      (specific_energy_rates.ske_rate.estimate * weight.ske_weighting);

  return seb_rate;
}

void TecsControl::calcPitchControlUpdate(double dt, const TecsInput& input,
                                         const ControlValues& seb_rate, const TecsParam& param) {
  if (param.integrator_gain_pitch > FLT_EPSILON) {
    // [PX4] Calculate derivative from change in climb angle to rate of change of specific
    // energy balance
    const double climb_angle_to_SEB_rate = std::max(input.tas, param.tas_min) * kOneG;

    // [PX4] Calculate pitch integrator input term
    double pitch_integ_input =
        getControlError(seb_rate) * param.integrator_gain_pitch / climb_angle_to_SEB_rate;

    // [PX4] Prevent the integrator changing in a direction that will increase pitch demand
    // saturation
    if (pitch_setpoint_ >= param.pitch_max) {
      pitch_integ_input = std::min(pitch_integ_input, 0.0);
    } else if (pitch_setpoint_ <= param.pitch_min) {
      pitch_integ_input = std::max(pitch_integ_input, 0.0);
    }

    // [PX4] Update the pitch integrator state
    pitch_integ_state_ = pitch_integ_state_ + pitch_integ_input * dt;
  } else {
    pitch_integ_state_ = 0.0;
  }
}

double TecsControl::calcPitchControlOutput(const TecsInput& input, const ControlValues& seb_rate,
                                           const TecsParam& param, const TecsFlag& flag) const {
  double airspeed_for_seb_rate = param.equivalent_airspeed_trim;

  // [PX4] avoid division by zero by checking if airspeed is finite and greater than zero
  if (flag.airspeed_enabled && std::isfinite(input.tas) && input.tas > FLT_EPSILON) {
    airspeed_for_seb_rate = input.tas;
  }

  // [PX4] Calculate derivative from change in climb angle to rate of change of specific energy
  // balance
  const double climb_angle_to_SEB_rate = airspeed_for_seb_rate * kOneG;

  // [PX4] Calculate a specific energy correction that doesn't include the integrator contribution
  const double SEB_rate_correction =
      getControlError(seb_rate) * param.pitch_damping_gain + param.seb_rate_ff * seb_rate.setpoint;

  // [PX4] Convert the specific energy balance rate correction to a target pitch angle. This
  // calculation assumes:
  // a) The climb angle follows pitch angle with a lag that is small enough not to destabilise
  //    the control loop.
  // b) The offset between climb angle and pitch angle (angle of attack) is constant, excluding
  //    the effect of pitch transients due to control action or turbulence.
  const double pitch_setpoint_unc = SEB_rate_correction / climb_angle_to_SEB_rate + pitch_integ_state_;

  return constrain(pitch_setpoint_unc, param.pitch_min, param.pitch_max);
}

void TecsControl::calcThrottleControl(double dt, const SpecificEnergyRates& specific_energy_rates,
                                      const TecsParam& param, const TecsFlag& flag) {
  const STERateLimit limit{calculateTotalEnergyRateLimit(param)};

  // [PX4] Update STE rate estimate LP filter
  // (AlphaFilter::setParameters(dt, tau): alpha = dt / (tau + dt); update: y += alpha (u - y).)
  const double STE_rate_estimate_raw =
      specific_energy_rates.spe_rate.estimate + specific_energy_rates.ske_rate.estimate;
  const double denominator = param.ste_rate_time_const + dt;
  const double alpha = (denominator > FLT_EPSILON) ? dt / denominator : 0.0;
  ste_rate_estimate_filter_ += alpha * (STE_rate_estimate_raw - ste_rate_estimate_filter_);

  const ControlValues ste_rate{calcThrottleControlSteRate(limit, specific_energy_rates, param)};

  // (Upstream branches on fast_descend == 1 -> throttle_min; not ported.)
  calcThrottleControlUpdate(dt, limit, ste_rate, param, flag);
  double throttle_setpoint = calcThrottleControlOutput(limit, ste_rate, param, flag);

  // [PX4] Rate limit the throttle demand
  if (std::fabs(param.throttle_slewrate) > FLT_EPSILON) {
    const double throttle_increment_limit =
        dt * (param.throttle_max - param.throttle_min) * param.throttle_slewrate;
    throttle_setpoint = constrain(throttle_setpoint, throttle_setpoint_ - throttle_increment_limit,
                                  throttle_setpoint_ + throttle_increment_limit);
  }

  throttle_setpoint_ = constrain(throttle_setpoint, param.throttle_min, param.throttle_max);

  // Debug output
  debug_output_.total_energy_rate_estimate = ste_rate.estimate;
  debug_output_.total_energy_rate_sp = ste_rate.setpoint;
  debug_output_.throttle_integrator = throttle_integ_state_;
}

TecsControl::ControlValues TecsControl::calcThrottleControlSteRate(
    const STERateLimit& limit, const SpecificEnergyRates& specific_energy_rates,
    const TecsParam& param) const {
  ControlValues ste_rate;
  ste_rate.setpoint = specific_energy_rates.spe_rate.setpoint + specific_energy_rates.ske_rate.setpoint;

  // [PX4] Adjust the demanded total energy rate to compensate for induced drag rise in turns.
  // Assume induced drag scales linearly with normal load factor.
  // The additional normal load factor is given by (1/cos(bank angle) - 1)
  ste_rate.setpoint += param.load_factor_correction * (param.load_factor - 1.0);

  ste_rate.setpoint = constrain(ste_rate.setpoint, limit.STE_rate_min, limit.STE_rate_max);
  ste_rate.estimate = ste_rate_estimate_filter_;

  return ste_rate;
}

void TecsControl::calcThrottleControlUpdate(double dt, const STERateLimit& limit,
                                            const ControlValues& ste_rate, const TecsParam& param,
                                            const TecsFlag& flag) {
  // [PX4] Calculate gain scaler from specific energy rate error to throttle
  const double STE_rate_to_throttle = 1.0 / (limit.STE_rate_max - limit.STE_rate_min);

  // [PX4] Integral handling
  if (flag.airspeed_enabled) {
    if (param.integrator_gain_throttle > FLT_EPSILON) {
      // [PX4] underspeed conditions zero out integration
      double throttle_integ_input = (getControlError(ste_rate) * param.integrator_gain_throttle) * dt *
                                    STE_rate_to_throttle * (1.0 - ratio_undersped_);

      // [PX4] only allow integrator propagation into direction which unsaturates throttle
      if (throttle_setpoint_ >= param.throttle_max) {
        throttle_integ_input = std::min(0.0, throttle_integ_input);
      }
      if (throttle_setpoint_ <= param.throttle_min) {
        throttle_integ_input = std::max(0.0, throttle_integ_input);
      }

      // [PX4] Calculate a throttle demand from the integrated total energy rate error
      // This will be added to the total throttle demand to compensate for steady state errors
      throttle_integ_state_ = throttle_integ_state_ + throttle_integ_input;
    } else {
      throttle_integ_state_ = 0.0;
    }
  }
}

double TecsControl::calcThrottleControlOutput(const STERateLimit& limit, const ControlValues& ste_rate,
                                              const TecsParam& param, const TecsFlag& flag) const {
  // [PX4] Calculate gain scaler from specific energy rate error to throttle
  const double STE_rate_to_throttle = 1.0 / (limit.STE_rate_max - limit.STE_rate_min);

  // [PX4] Calculate a predicted throttle from the demanded rate of change of energy, using the
  // cruise throttle as the starting point. Assume:
  // Specific total energy rate = _STE_rate_max is achieved when throttle is set to
  //   _throttle_setpoint_max
  // Specific total energy rate = 0 at cruise throttle
  // Specific total energy rate = _STE_rate_min is achieved when throttle is set to
  //   _throttle_setpoint_min
  // assume airspeed and density-independent delta_throttle to sink/climb rate mapping
  const double throttle_above_trim_per_ste_rate =
      (param.throttle_max - param.throttle_trim) / limit.STE_rate_max;
  const double throttle_below_trim_per_ste_rate =
      (param.throttle_trim - param.throttle_min) / limit.STE_rate_min;

  double throttle_predicted = 0.0;

  if (ste_rate.setpoint >= FLT_EPSILON) {
    // [PX4] throttle is between trim and maximum
    throttle_predicted = param.throttle_trim + ste_rate.setpoint * throttle_above_trim_per_ste_rate;
  } else {
    // [PX4] throttle is between trim and minimum
    throttle_predicted = param.throttle_trim - ste_rate.setpoint * throttle_below_trim_per_ste_rate;
  }

  // [PX4] Add proportional and derivative control feedback to the predicted throttle and
  // constrain to throttle limits
  double throttle_setpoint =
      (getControlError(ste_rate) * param.throttle_damping_gain) * STE_rate_to_throttle + throttle_predicted;

  if (flag.airspeed_enabled) {
    // [PX4] Add the integrator feedback during closed loop operation with an airspeed sensor
    throttle_setpoint += throttle_integ_state_;
  } else {
    // [PX4] We want to avoid reducing the throttle output when switching from airspeed enabled
    // mode into airspeedless mode. Thus, if the throttle integrator has a positive value, add it
    // still to the throttle setpoint.
    throttle_setpoint += std::max(0.0, throttle_integ_state_);
  }

  // [PX4] ramp in max throttle setting with underspeediness value
  throttle_setpoint = ratio_undersped_ * param.throttle_max + (1.0 - ratio_undersped_) * throttle_setpoint;

  return constrain(throttle_setpoint, param.throttle_min, param.throttle_max);
}

void TecsControl::resetIntegrals() {
  pitch_integ_state_ = 0.0;
  throttle_integ_state_ = 0.0;
}

}  // namespace px4
}  // namespace autoland
