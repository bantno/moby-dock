#pragma once
#include <cmath>
#include <functional>
#include "autoland/world_geometry.hpp"

// =============================================================================
// Post-touchdown planar rollout on the water. Once the keel is on the surface
// the 6-DOF air model stops (no hydrodynamic float physics in this repo); the
// stopping run is a 2-D point mass on the water plane:
//
//   Ndot   = Vg cos(chi)                  chi: ground course [rad]
//   Edot   = Vg sin(chi)                  Vg : ground speed  [m/s]
//   chidot = -k_psi wrap(chi - psi_c)     water-rudder / pilot steering proxy
//   Vgdot  = -(a0 + kq Vg^2)              hydrodynamic decel, Vg >= 0
//
// with psi_c the corridor heading. Closed forms (used by the tests):
//   kq = 0 : d_stop = V0^2 / (2 a0),  t_stop = V0 / a0
//   a0 = 0 : Vg(t) = V0 / (1 + kq V0 t), d(t) = ln(1 + kq V0 t) / kq
// The run ends when Vg <= V_stop or the time cap hits. Integrated with a
// local RK4 at the sim step.
// =============================================================================
namespace autoland {

struct RolloutState {
  double N{0}, E{0}, chi{0}, Vg{0};
};

struct RolloutResult {
  bool ran{false};
  bool stopped{false};      // Vg reached V_stop before t_max
  double t_start{0}, t_stop{0};
  double N_stop{0}, E_stop{0}, chi_stop{0};
  double distance{0};       // path length of the run [m]
  double s_stop{0}, e_stop{0};   // stop point in the corridor frame
  double max_abs_e{0};      // largest lateral excursion during the run
  bool inside{false};       // stop point inside the corridor
  bool stayed_inside{true}; // every sample inside the corridor
};

inline double wrapAngle(double a) { return std::atan2(std::sin(a), std::cos(a)); }

inline RolloutState rolloutRhs(const RolloutState& s, const RolloutConfig& c,
                               double psi_c) {
  RolloutState d;
  d.N = s.Vg * std::cos(s.chi);
  d.E = s.Vg * std::sin(s.chi);
  d.chi = -c.k_psi * wrapAngle(s.chi - psi_c);
  d.Vg = -(c.a0 + c.kq * s.Vg * s.Vg);
  return d;
}

// Integrate the run from the touchdown state. `on_step(t, state)` is called
// at every sample including the first and the last (CSV logging hook).
inline RolloutResult runRollout(
    RolloutState s, double t0, double dt, const RolloutConfig& c,
    const Corridor& cor,
    const std::function<void(double, const RolloutState&)>& on_step = {}) {
  RolloutResult r;
  r.ran = true;
  r.t_start = t0;
  double t = t0;
  auto axpy = [](const RolloutState& a, double h, const RolloutState& d) {
    RolloutState o;
    o.N = a.N + h * d.N; o.E = a.E + h * d.E;
    o.chi = a.chi + h * d.chi; o.Vg = a.Vg + h * d.Vg;
    return o;
  };
  auto sample = [&](double tt, const RolloutState& st) {
    double sc, ec;
    cor.frame(st.N, st.E, sc, ec);
    r.max_abs_e = std::max(r.max_abs_e, std::abs(ec));
    if (!cor.contains(st.N, st.E)) r.stayed_inside = false;
    if (on_step) on_step(tt, st);
  };
  sample(t, s);
  while (t - t0 < c.t_max) {
    if (s.Vg <= c.V_stop) { r.stopped = true; break; }
    const RolloutState k1 = rolloutRhs(s, c, cor.heading);
    const RolloutState k2 = rolloutRhs(axpy(s, 0.5 * dt, k1), c, cor.heading);
    const RolloutState k3 = rolloutRhs(axpy(s, 0.5 * dt, k2), c, cor.heading);
    const RolloutState k4 = rolloutRhs(axpy(s, dt, k3), c, cor.heading);
    RolloutState n;
    n.N = s.N + dt / 6.0 * (k1.N + 2 * k2.N + 2 * k3.N + k4.N);
    n.E = s.E + dt / 6.0 * (k1.E + 2 * k2.E + 2 * k3.E + k4.E);
    n.chi = s.chi + dt / 6.0 * (k1.chi + 2 * k2.chi + 2 * k3.chi + k4.chi);
    n.Vg = std::max(0.0, s.Vg + dt / 6.0 * (k1.Vg + 2 * k2.Vg + 2 * k3.Vg + k4.Vg));
    r.distance += std::hypot(n.N - s.N, n.E - s.E);
    s = n;
    t += dt;
    sample(t, s);
  }
  r.t_stop = t;
  r.N_stop = s.N; r.E_stop = s.E; r.chi_stop = s.chi;
  cor.frame(s.N, s.E, r.s_stop, r.e_stop);
  r.inside = cor.contains(s.N, s.E);
  return r;
}

}  // namespace autoland
