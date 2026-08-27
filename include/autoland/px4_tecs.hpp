#pragma once
#include <cfloat>
#include <cmath>

// =============================================================================
// PX4 TECS (Total Energy Control System) -- the CONTROL LAW only.
//
// Port of `TECSControl` from PX4-Autopilot, src/lib/tecs/TECS.hpp + TECS.cpp,
// main branch @ a906b72868d79c89790688e5b41de73597133ce5 (fetched 2026-08-27;
// upstream author Paul Riseborough). Same function decomposition, member
// names, and constants as upstream, in double precision, so the file can be
// read side by side with TECS.cpp. Parameter DEFAULTS are the flown PX4
// parameter defaults (FW_T_* in fw_lateral_longitudinal_control/
// fw_lat_long_params.yaml at the same commit), not the placeholder
// initialisers inside TECS.hpp (which the module overwrites via setters).
//
// PORTED (everything inside TECSControl): the altitude-error and airspeed-
// error outer loops; the specific-energy-rate bookkeeping (SPE_rate = g*hdot,
// SKE_rate = V*Vdot); the energy-BALANCE rate (SEB) -> pitch loop with feed-
// forward, damping, integrator + anti-windup on the pitch limits, and the
// vertical-acceleration pitch-rate limit; the TOTAL energy rate (STE) ->
// throttle loop with the predicted throttle about trim, damping, integrator +
// anti-windup, underspeed ramp-in of full throttle and the slew limit; the
// first-order STE-rate estimate filter (AlphaFilter); underspeed detection;
// the speed/altitude weighting; the bank-angle induced-drag compensation.
//
// NOT PORTED (PX4 overhead outside the control law): the uORB/timestamp
// wrapper (`TECS`), the jerk-limited altitude reference model
// (TECSAltitudeReferenceModel / VelocitySmoothing -- the sim commands a direct
// height-rate setpoint, PX4's `altitude_rate_setpoint_direct` path), the
// airspeed complementary filter (the sim supplies the EXACT airspeed rate),
// the "fast descend" mode (its hooks inside TECSControl are evaluated at
// fast_descend = 0 and dropped), the airspeed-less branches (the flag is kept
// so the structure matches; the sim always has airspeed), and the NaN guards
// (double-precision sim inputs are finite by construction).
//
// CONVENTIONS (upstream): altitude positive UP. The pitch setpoint is RELATIVE
// to the level-flight trim pitch (PX4 subtracts FW_PSP_OFF from the measured
// pitch and the pitch limits before TECS and adds it back to the output; the
// caller here does the same with the level-trim theta). Throttle in [0,1].
// Specific energy rates in m^2/s^3. The direct height-rate setpoint uses NAN
// as the "not set" sentinel exactly as upstream does.
// =============================================================================
namespace autoland {
namespace px4 {

constexpr double kOneG = 9.80665;  // CONSTANTS_ONE_G (lib/geo/geo.h)

// TECSControl::Param. Comments give the PX4 parameter each field is fed from
// and that parameter's default; vehicle-specific fields (rates, trims, limits)
// are filled by the caller from the plant (see sixdof_sim.cpp).
struct TecsParam {
  // Vehicle specific params
  double max_sink_rate{5.0};        // FW_T_SINK_MAX = 5   [m/s] (min throttle, max speed)
  double min_sink_rate{2.0};        // FW_T_SINK_MIN = 2   [m/s] (min throttle, trim speed)
  double max_climb_rate{5.0};       // FW_T_CLMB_MAX = 5   [m/s] (max throttle)
  double vert_accel_limit{7.0};     // FW_T_VERT_ACC = 7   [m/s^2]
  double equivalent_airspeed_trim{15.0};  // FW_AIRSPD_TRIM = 15 [m/s]
  double tas_min{10.0};             // FW_AIRSPD_MIN = 10  [m/s] (demand lower limit)
  double tas_max{20.0};             // FW_AIRSPD_MAX = 20  [m/s] (demand upper limit)
  double pitch_max{0.5};            // above trim [rad] (FW_P_LIM_MAX - FW_PSP_OFF)
  double pitch_min{-0.5};           // below trim [rad] (FW_P_LIM_MIN - FW_PSP_OFF)
  double throttle_trim{0.6};        // FW_THR_TRIM = 0.6 : level flight at the airspeed setpoint
  double throttle_max{1.0};         // FW_THR_MAX = 1
  double throttle_min{0.0};         // FW_THR_MIN = 0

  // Altitude control param
  double altitude_error_gain{1.0 / 5.0};     // 1 / FW_T_ALT_TC (5 s)  [1/s]
  double altitude_setpoint_gain_ff{0.5};     // FW_T_HRATE_FF = 0.5

  // Airspeed control param
  double tas_error_percentage{0.15};         // hard-coded upstream (TECS.hpp)
  double airspeed_error_gain{1.0 / 5.0};     // 1 / FW_T_TAS_TC (5 s)  [1/s]

  // Energy control param
  double ste_rate_time_const{0.4};           // FW_T_STE_R_TC = 0.4 [s]
  double seb_rate_ff{1.0};                   // FW_T_SEB_R_FF = 1.0

  // Pitch control param
  double pitch_speed_weight{1.0};            // FW_T_SPDWEIGHT = 1.0  (0 = height only, 2 = speed only)
  double integrator_gain_pitch{0.1};         // FW_T_I_GAIN_PIT = 0.1
  double pitch_damping_gain{0.1};            // FW_T_PTCH_DAMP = 0.1 [s]

  // Throttle control param
  double integrator_gain_throttle{0.02};     // FW_T_THR_INTEG = 0.02
  double throttle_damping_gain{0.05};        // FW_T_THR_DAMPING = 0.05 [s]
  double throttle_slewrate{0.0};             // FW_THR_SLEW_MAX = 0 (off) [1/s]

  double load_factor_correction{15.0};       // FW_T_RLL2THR = 15 [m^2/s^3]
  double load_factor{1.0};                   // 1 / cos(bank), set per step by the caller
};

// TECSControl::DebugOutput.
struct TecsDebugOutput {
  double altitude_rate_control{0.0};          // altitude-rate setpoint from the altitude loop [m/s]
  double true_airspeed_derivative_control{0.0};  // airspeed-rate setpoint from the airspeed loop [m/s^2]
  double total_energy_rate_estimate{0.0};     // [m^2/s^3]
  double total_energy_rate_sp{0.0};
  double energy_balance_rate_estimate{0.0};
  double energy_balance_rate_sp{0.0};
  double pitch_integrator{0.0};               // [rad]
  double throttle_integrator{0.0};            // [-]
};

// TECSControl::Setpoint.
struct TecsSetpoint {
  struct AltitudeReference {
    double alt{0.0};       // reference altitude [m]
    double alt_rate{0.0};  // reference altitude rate [m/s]
  } altitude_reference;
  double altitude_rate_setpoint_direct{NAN};  // direct height-rate setpoint; NAN => altitude loop
  double tas_setpoint{0.0};                   // true airspeed setpoint [m/s]
};

// TECSControl::Input.
struct TecsInput {
  double altitude{0.0};       // [m]
  double altitude_rate{0.0};  // [m/s], positive up
  double tas{0.0};            // true airspeed [m/s]
  double tas_rate{0.0};       // true airspeed rate [m/s^2]
};

// TECSControl::Flag.
struct TecsFlag {
  bool airspeed_enabled{true};
  bool detect_underspeed_enabled{true};
};

class TecsControl {
 public:
  TecsControl() = default;

  // TECSControl::initialize -- resets the integrators and computes the first
  // pitch/throttle output (PX4 calls this on the first update and on a dt reset).
  void initialize(const TecsSetpoint& setpoint, const TecsInput& input,
                  const TecsParam& param, const TecsFlag& flag);

  // TECSControl::update -- one control step of dt seconds.
  void update(double dt, const TecsSetpoint& setpoint, const TecsInput& input,
              const TecsParam& param, const TecsFlag& flag);

  void resetIntegrals();

  double getRatioUndersped() const { return ratio_undersped_; }
  double getThrottleSetpoint() const { return throttle_setpoint_; }
  double getPitchSetpoint() const { return pitch_setpoint_; }  // above trim [rad]
  const TecsDebugOutput& getDebugOutput() const { return debug_output_; }

 private:
  struct STERateLimit {
    double STE_rate_max;
    double STE_rate_min;
  };
  struct ControlValues {
    double setpoint;
    double estimate;
  };
  struct SpecificEnergyRates {
    ControlValues ske_rate;
    ControlValues spe_rate;
  };
  struct AltitudePitchControl {
    double altitude_rate_setpoint;
    double tas_rate_setpoint;
    double tas_setpoint;
  };
  struct SpecificEnergyWeighting {
    double spe_weighting;
    double ske_weighting;
  };

  static constexpr double getControlError(ControlValues val) {
    return val.setpoint - val.estimate;
  }
  STERateLimit calculateTotalEnergyRateLimit(const TecsParam& param) const;
  double calcAirspeedControlOutput(const TecsSetpoint& setpoint, const TecsInput& input,
                                   const TecsParam& param, const TecsFlag& flag) const;
  double calcAltitudeControlOutput(const TecsSetpoint& setpoint, const TecsInput& input,
                                   const TecsParam& param) const;
  SpecificEnergyRates calcSpecificEnergyRates(const AltitudePitchControl& control_setpoint,
                                              const TecsInput& input) const;
  void detectUnderspeed(const TecsInput& input, const TecsParam& param, const TecsFlag& flag);
  SpecificEnergyWeighting updateSpeedAltitudeWeights(const TecsParam& param,
                                                     const TecsFlag& flag) const;
  void calcPitchControl(double dt, const TecsInput& input,
                        const SpecificEnergyRates& specific_energy_rate,
                        const TecsParam& param, const TecsFlag& flag);
  ControlValues calcPitchControlSebRate(const SpecificEnergyWeighting& weight,
                                        const SpecificEnergyRates& specific_energy_rate) const;
  void calcPitchControlUpdate(double dt, const TecsInput& input, const ControlValues& seb_rate,
                              const TecsParam& param);
  double calcPitchControlOutput(const TecsInput& input, const ControlValues& seb_rate,
                                const TecsParam& param, const TecsFlag& flag) const;
  void calcThrottleControl(double dt, const SpecificEnergyRates& specific_energy_rate,
                           const TecsParam& param, const TecsFlag& flag);
  ControlValues calcThrottleControlSteRate(const STERateLimit& limit,
                                           const SpecificEnergyRates& specific_energy_rate,
                                           const TecsParam& param) const;
  void calcThrottleControlUpdate(double dt, const STERateLimit& limit,
                                 const ControlValues& ste_rate, const TecsParam& param,
                                 const TecsFlag& flag);
  double calcThrottleControlOutput(const STERateLimit& limit, const ControlValues& ste_rate,
                                   const TecsParam& param, const TecsFlag& flag) const;

  // State
  double ste_rate_estimate_filter_{0.0};  // AlphaFilter<float> _ste_rate_estimate_filter
  double pitch_integ_state_{0.0};         // [rad]
  double throttle_integ_state_{0.0};      // [-]

  // Output
  TecsDebugOutput debug_output_;
  double pitch_setpoint_{0.0};     // above trim [rad]
  double throttle_setpoint_{0.0};  // [0,1]
  double ratio_undersped_{0.0};    // [0,1]
};

}  // namespace px4
}  // namespace autoland
