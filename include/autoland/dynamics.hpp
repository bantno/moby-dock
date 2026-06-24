#pragma once
#include <functional>
#include "autoland/aero_table.hpp"
#include "autoland/config.hpp"
#include "autoland/mixing.hpp"
#include "autoland/types.hpp"

// =============================================================================
// Dynamics: the ONE nonlinear 6-DOF equation of motion, Dynamics::xdot(x,u).
//
// This is the single source of truth. The trim solver and the linear model are
// both derived from this same function -- the linear model is produced purely
// by central-differencing xdot(), and trim is a Newton solve on xdot(). A
// nonlinear plant can therefore be swapped into the sim later with no rework.
//
// Body axes, SI units. Forces/moments are assembled from the BODY-AXIS aero
// coefficients (CFx,CFy,CFz,Cl=CMx,Cm=CMy,Cn=CMz), never from wind-axis CL/CD.
// Sign conventions are documented in linear_model.hpp.
// =============================================================================
namespace autoland {

// Pluggable thrust model. Default implements T = throttle*(T_static - k_v*V),
// clamped to >= 0. Replace this std::function to drop in thrust-stand data.
using ThrustModel = std::function<double(double throttle, double V)>;

class Dynamics {
 public:
  Dynamics(const AeroTable& table, const Mixing& mixing,
           const AircraftConfig& cfg);

  // Nonlinear state derivative. u = [delta_e, delta_a, delta_r, delta_T].
  StateVec xdot(const StateVec& x, const CtrlVec& u) const;

  // Assemble the 6 body-axis coefficients at a flight condition. Exposed so
  // tests and analysis can inspect the aero buildup directly.
  CoefVec aeroCoeffs(const StateVec& x, const CtrlVec& u) const;

  // Override the thrust model (e.g. measured data).
  void setThrustModel(ThrustModel m) { thrust_ = std::move(m); }

  const AeroTable& table() const { return table_; }
  const AircraftConfig& config() const { return cfg_; }

 private:
  const AeroTable& table_;
  const Mixing& mixing_;
  AircraftConfig cfg_;
  ThrustModel thrust_;
};

// Hand-rolled classic RK4 step for xdot-style dynamics with zero-order-hold u.
StateVec rk4Step(const std::function<StateVec(const StateVec&, const CtrlVec&)>& f,
                 const StateVec& x, const CtrlVec& u, double dt);

}  // namespace autoland
