#pragma once
#include <memory>
#include <string>
#include "autoland/aero_table.hpp"
#include "autoland/config.hpp"
#include "autoland/lon_cbf_filter.hpp"
#include "autoland/lon_nominal.hpp"
#include "autoland/mixing.hpp"
#include "autoland/wind_gust.hpp"

// =============================================================================
// Self-consistent augmented-longitudinal landing simulation: the plant IS the
// CBF model (Xdot = f(X) + g(X) U*), so there is zero model mismatch -- a clean
// test of the CBF/QP math. RK4 integration with a zero-order-hold control.
// =============================================================================
namespace autoland {

// Scenario bundle. Defaults are derived from a trim solve at (V_app, gamma_app);
// the YAML may override any field. Held together here to avoid coupling the
// generic config.hpp to the longitudinal control headers.
struct LonScenario {
  LonNominalConfig nominal;
  LonCBFConfig cbf;
  DiscreteGustConfig wind;  // MIL-F-8785C discrete gust, plant-side only
  LonStateVec X0{LonStateVec::Zero()};
  double dt{0.01};
  double t_max{60.0};
  bool cbf_enabled{true};
};

struct LonTouchdown {
  bool reached{false};
  double t{0}, sink{0}, V{0}, theta{0}, gamma{0};
};

class LonSim {
 public:
  LonSim(const std::string& stab_path, const std::string& aircraft_yaml,
         const std::string& scenario_yaml);

  // Run to touchdown (h<=0) or t_max. Writes a CSV trace to csv_path.
  LonTouchdown run(const std::string& csv_path);

  const LonScenario& scenario() const { return sc_; }

 private:
  AeroTable table_;
  AircraftConfig ac_;
  std::unique_ptr<Mixing> mixing_;
  LonScenario sc_;
};

// Full augmented RHS with aero rebuilt from the state (plant + control).
LonStateVec lonXdotFull(const AeroTable& table, const Mixing& mixing,
                        const AircraftConfig& cfg, const LonStateVec& X,
                        const LonCtrlVec& U);

// Wind-perturbed plant RHS. The state keeps its air-relative meaning (V =
// airspeed, gamma = air-path angle, so alpha = theta - gamma is still the true
// aerodynamic AoA); the earth-frame gust W = (u tailwind+, w updraft+) enters
// exactly, via
//   hdot     += W.w                                     (kinematic transport)
//   Vdot     += -( Wdot.u cos(gamma) + Wdot.w sin(gamma) )
//   gammadot += ( Wdot.u sin(gamma) - Wdot.w cos(gamma) ) / V
// i.e. the -m Wdot inertial forcing of m d(v_air)/dt = F - m Wdot projected on
// the air-path axes. A steady wind (Wdot = 0) therefore only transports the
// aircraft; only the gust RAMP forces V/gamma -- the classic airspeed loss in
// a tailwind ramp and the AoA rise in an updraft ramp.
LonStateVec lonXdotFullWind(const AeroTable& table, const Mixing& mixing,
                            const AircraftConfig& cfg, const LonStateVec& X,
                            const LonCtrlVec& U, const GustWind& W,
                            const GustWind& Wdot);

}  // namespace autoland
