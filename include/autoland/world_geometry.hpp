#pragma once
#include <Eigen/Dense>
#include <autodiff/forward/dual.hpp>
#include <cmath>
#include <string>
#include <vector>

// =============================================================================
// World geometry for the corridor landing problem: an oriented landing
// CORRIDOR on the water and a smooth TERRAIN height field for keep-out.
//
// FRAME: earth-fixed, N north, E east, h up (matches the plant's XN / Y / H
// states). The corridor frame (s, e) is right-handed with s along the landing
// heading psi_c (positive toward the far end) and e to the right of it:
//   s =  (N - N_c) cos psi_c + (E - E_c) sin psi_c
//   e = -(N - N_c) sin psi_c + (E - E_c) cos psi_c
//
// TERRAIN: h_T(N, E) = base + sum_i A_i exp(-q_i^order_i), with q_i the
// squared elliptical distance in the bump's rotated frame. Everything is a
// composition of *, +, exp so the field is C-infinity and evaluates exactly
// under the Taylor-jet / autodiff scalar types of lie_taylor.hpp (the barrier
// h - h_T(N,E) will need L_f^3 of it). Integer `order` (1 = Gaussian,
// 2 = flat-top quartic) is deliberate: no pow(), no log().
//
// The lake is the base plane (0); trees around the shore are a ring of
// small flat-top bumps, a mountain is one large Gaussian, houses are small
// steep quartics.
// =============================================================================
namespace autoland {

struct Corridor {
  double N_c{0.0}, E_c{0.0};   // center [m]
  double length{700.0};        // along-track extent [m]
  double width{80.0};          // across-track extent [m]
  double heading{0.0};         // landing direction psi_c [rad], N -> E positive
  double s_aim{-250.0};        // touchdown aim point, corridor frame [m]

  template <class T>
  void frame(const T& N, const T& E, T& s, T& e) const {
    const double c = std::cos(heading), sn = std::sin(heading);
    const T dN = N - N_c, dE = E - E_c;
    s = dN * c + dE * sn;
    e = dE * c - dN * sn;
  }
  template <class T>
  void inverse(const T& s, const T& e, T& N, T& E) const {
    const double c = std::cos(heading), sn = std::sin(heading);
    N = s * c - e * sn + N_c;
    E = s * sn + e * c + E_c;
  }
  double sNear() const { return -0.5 * length; }
  double sFar() const { return 0.5 * length; }
  bool contains(double N, double E, double margin = 0.0) const {
    double s, e;
    frame(N, E, s, e);
    return std::abs(s) <= 0.5 * length - margin &&
           std::abs(e) <= 0.5 * width - margin;
  }
  Eigen::Vector2d axisPoint(double s) const {
    double N, E, e = 0.0;
    inverse(s, e, N, E);
    return {N, E};
  }
  Eigen::Vector2d aimPoint() const { return axisPoint(s_aim); }
  Eigen::Vector2d nearEdge() const { return axisPoint(sNear()); }
  Eigen::Vector2d farEdge() const { return axisPoint(sFar()); }
};

struct TerrainBump {
  std::string name;
  double N0{0.0}, E0{0.0};     // center [m]
  double height{0.0};          // peak height above base [m]
  double sigma_a{1.0};         // semi-axis scale along the bump's a-axis [m]
  double sigma_b{1.0};         // ... b-axis [m]
  double rot{0.0};             // a-axis heading [rad], N -> E positive
  int order{1};                // 1 Gaussian, 2 flat-top (exp(-q^2)), ...
};

struct TerrainField {
  double base{0.0};
  std::vector<TerrainBump> bumps;

  bool empty() const { return bumps.empty(); }

  // Height above the water datum at (N, E). T = double, autodiff::dual, or
  // a Taylor jet -- only *, +, - and exp are used.
  template <class T>
  T height(const T& N, const T& E) const {
    using std::exp;
    T h = T(base);
    for (const TerrainBump& b : bumps) {
      const double c = std::cos(b.rot), s = std::sin(b.rot);
      const T dN = N - b.N0, dE = E - b.E0;
      const T xa = (dN * c + dE * s) * (1.0 / b.sigma_a);
      const T xb = (dE * c - dN * s) * (1.0 / b.sigma_b);
      const T q = xa * xa + xb * xb;
      T qp = q;
      for (int k = 1; k < b.order; ++k) qp = qp * q;
      h = h + b.height * exp(-qp);
    }
    return h;
  }
  double height(double N, double E) const { return height<double>(N, E); }

  // Closed-form gradient [dh/dN, dh/dE], scalar-generic (the bump sum
  // differentiates termwise: d/dq A exp(-q^n) = -A n q^(n-1) exp(-q^n)).
  template <class T>
  std::array<T, 2> gradientT(const T& N, const T& E) const {
    using std::exp;
    std::array<T, 2> g{T(0.0), T(0.0)};
    for (const TerrainBump& b : bumps) {
      const double c = std::cos(b.rot), s = std::sin(b.rot);
      const T dN = N - b.N0, dE = E - b.E0;
      const T xa = (dN * c + dE * s) * (1.0 / b.sigma_a);
      const T xb = (dE * c - dN * s) * (1.0 / b.sigma_b);
      const T q = xa * xa + xb * xb;
      T qn1 = T(1.0);                         // q^(n-1)
      for (int k = 1; k < b.order; ++k) qn1 = qn1 * q;
      const T dh_dq = (-b.height * b.order) * qn1 * exp(-(qn1 * q));
      const T dq_dN = 2.0 * xa * (c / b.sigma_a) - 2.0 * xb * (s / b.sigma_b);
      const T dq_dE = 2.0 * xa * (s / b.sigma_a) + 2.0 * xb * (c / b.sigma_b);
      g[0] = g[0] + dh_dq * dq_dN;
      g[1] = g[1] + dh_dq * dq_dE;
    }
    return g;
  }

  // Exact gradient [dh/dN, dh/dE] via autodiff::dual.
  Eigen::Vector2d gradient(double N, double E) const {
    using autodiff::dual;
    Eigen::Vector2d g;
    for (int k = 0; k < 2; ++k) {
      dual n = N, e = E;
      if (k == 0) n.grad = 1.0; else e.grad = 1.0;
      g[k] = autodiff::detail::derivative(height<dual>(n, e));
    }
    return g;
  }

  // Ring of `count` bumps on a circle (trees around a lake shore). An
  // optional gap: bumps whose bearing from the center is within
  // `skip_halfangle` of `skip_heading` (either direction) are omitted.
  void addRing(double Nc, double Ec, double radius, int count, double height,
               double sigma, int order, const std::string& name = "ring",
               double skip_heading = NAN, double skip_halfangle = 0.0) {
    for (int i = 0; i < count; ++i) {
      const double th = 2.0 * M_PI * i / count;
      if (std::isfinite(skip_heading)) {
        const double d1 = std::remainder(th - skip_heading, 2.0 * M_PI);
        const double d2 = std::remainder(th - skip_heading - M_PI, 2.0 * M_PI);
        if (std::abs(d1) <= skip_halfangle || std::abs(d2) <= skip_halfangle)
          continue;
      }
      TerrainBump b;
      b.name = name + "_" + std::to_string(i);
      b.N0 = Nc + radius * std::cos(th);
      b.E0 = Ec + radius * std::sin(th);
      b.height = height;
      b.sigma_a = b.sigma_b = sigma;
      b.rot = 0.0;
      b.order = order;
      bumps.push_back(b);
    }
  }
};

// Post-touchdown planar rollout (water run) parameters -- see
// water_rollout.hpp. a0 / kq are PLACEHOLDERS pending POH calibration.
struct RolloutConfig {
  bool enabled{true};
  double a0{0.8};      // constant hydrodynamic deceleration [m/s^2]
  double kq{0.002};    // quadratic drag deceleration coefficient [1/m]
  double k_psi{0.5};   // course steering toward the corridor heading [1/s]
  double V_stop{0.5};  // "stopped" ground-speed threshold [m/s]
  double t_max{90.0};  // rollout time cap [s]
};

// Robustness margin for the Phase 2 terrain-clearance row, from stated
// uncertainties (all in metres / radians):
//   margin = map_err_v + grad_max * pos_err_h + (b/2) sin(phi_max) + extra
// map_err_v: vertical map error; pos_err_h: horizontal position error, priced
// through the largest terrain slope grad_max = sup |grad h_T| (for a bump of
// height A and scale sigma at order 1 it is A/(sigma) * sqrt(2/e), i.e. the
// Gaussian's steepest slope; use the per-bump maximum); the wingtip term
// converts CG clearance into low-wingtip clearance at the bank limit.
inline double terrainRowMargin(double map_err_v, double pos_err_h,
                               double grad_max, double span, double phi_max,
                               double extra = 0.0) {
  return map_err_v + grad_max * pos_err_h + 0.5 * span * std::sin(phi_max) + extra;
}
// Largest slope of the field: per-bump analytic bound, max over bumps.
inline double terrainMaxGradient(const TerrainField& tf) {
  double g = 0.0;
  for (const TerrainBump& b : tf.bumps) {
    const double s = std::min(b.sigma_a, b.sigma_b);
    // d/dx A exp(-(x/s)^(2n)) peaks at A/s * (2n) * q^(1-1/(2n)) e^{-q} with q
    // = (2n-1)/(2n); bounded by A/s * 2n for every order n >= 1.
    g = std::max(g, b.height / s * 2.0 * b.order);
  }
  return g;
}

struct WorldConfig {
  bool enabled{false};
  Corridor corridor;
  TerrainField terrain;
  RolloutConfig rollout;
};

}  // namespace autoland
