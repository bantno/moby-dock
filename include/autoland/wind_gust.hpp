#pragma once
#include <cmath>

// =============================================================================
// Discrete wind gust model -- the standard "1 - cosine" gust of U.S. Military
// Specification MIL-F-8785C (5 November 1980), as implemented by the Simulink
// Aerospace Blockset "Discrete Wind Gust Model" block:
//
//                 /  0                              x <  0
//   V_wind(x) = <   (Vm/2) (1 - cos(pi x / dm))     0 <= x <= dm
//                 \  Vm                             x >  dm
//
// where Vm is the gust amplitude [m/s], dm the gust length [m], and x the
// distance penetrated into the gust [m]. Like the Simulink block, x is driven
// by the aircraft airspeed:  xdot = V  for t >= t_start  (the block integrates
// its airspeed input the same way), and each axis has its own (Vm, dm) but all
// axes share the start time and penetration distance -- the block's
// "Gust amplitude [ug vg wg]" / "Gust length [dx dy dz]" parameterization,
// reduced here to the two axes that exist in the longitudinal vertical plane.
//
// MIL-F-8785C ties (Vm, dm) to the turbulence intensities/scale lengths of the
// mission segment; the Simulink block (and this model) leaves both as free,
// signed parameters, so headwind/tailwind and up/downdraft gusts are just sign
// choices on Vm.
//
// FRAME (differs from the spec's body axes, deliberately): this repo's
// longitudinal world is earth-frame with h UP, so the gust vector is given in
// the inertial vertical plane:
//   u  positive = tailwind (air moving along the direction of flight, +x)
//   w  positive = updraft  (air moving up, +h)
// A MIL/body-frame w_g (positive down) is the negation of w here.
//
// The gust is a PLANT disturbance only: the nominal controller and the CBF-QP
// keep the still-air model, so runs with a gust probe the robustness of the
// safety filter to unmeasured wind (TODO.md "gust/shear injection").
// x is a sim state (integrated alongside X in RK4); this header is pure shape
// evaluation and owns no state.
// =============================================================================
namespace autoland {

// One gust axis: the 1-cosine ramp V_wind(x) and its spatial slope. dm <= 0
// degenerates to a step gust (immediate full amplitude), matching the x >= dm
// branch rather than dividing by zero.
struct GustAxis {
  double amp{0.0};  // Vm [m/s], signed
  double len{0.0};  // dm [m]

  double value(double x) const {
    if (x <= 0.0 || amp == 0.0) return 0.0;
    if (len <= 0.0 || x >= len) return amp;
    return 0.5 * amp * (1.0 - std::cos(M_PI * x / len));
  }
  // dV_wind/dx; the time rate is slope(x) * xdot with xdot = V (see gustXdot).
  double slope(double x) const {
    if (x <= 0.0 || amp == 0.0 || len <= 0.0 || x >= len) return 0.0;
    return 0.5 * amp * (M_PI / len) * std::sin(M_PI * x / len);
  }
};

// Block parameters (YAML `wind:` section of the lon scenario). The lateral
// v axis exists for the 6-DOF sim (the longitudinal sims read only u/w); its
// step-gust limit (len <= 0) doubles as a steady crosswind after t_start.
struct DiscreteGustConfig {
  bool enabled{false};
  double t_start{0.0};  // gust onset time [s]
  GustAxis u;           // along-centerline axis, + = tailwind
  GustAxis v;           // lateral axis,          + = wind toward +y (east)
  GustAxis w;           // vertical axis,         + = updraft
};

// Earth-frame gust vector (values [m/s] or rates [m/s^2]).
struct GustWind {
  double u{0.0}, v{0.0}, w{0.0};
};

// Penetration-distance dynamics: xdot = V once the gust has started. x itself
// is integrated by the caller (it is a state of the simulation, like X).
inline double gustXdot(const DiscreteGustConfig& g, double t, double V) {
  return (g.enabled && t >= g.t_start) ? V : 0.0;
}

inline GustWind gustWind(const DiscreteGustConfig& g, double x) {
  if (!g.enabled) return {};
  return {g.u.value(x), g.v.value(x), g.w.value(x)};
}

// Time derivative Wdot = (dW/dx) xdot, needed by the wind-axis EOM forcing.
inline GustWind gustWindRate(const DiscreteGustConfig& g, double x, double xdot) {
  if (!g.enabled) return {};
  return {g.u.slope(x) * xdot, g.v.slope(x) * xdot, g.w.slope(x) * xdot};
}

}  // namespace autoland
