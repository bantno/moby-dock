#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include "autoland/impact_clf_table.hpp"
#include "autoland/lon_augmented.hpp"

// =============================================================================
// Hydrodynamic impact-load barrier (NACA TN 1516; see
// documentation/impact_load_barrier_spec.md). Bounds the peak CG load factor a
// water contact at the current state would produce:
//
//   b(X) = (n_limit - n_peak(tau, gamma, V)) + Phi(z),   Phi(z) = Nb(1 - e^{-z/zs})
//   n_peak = K0 * ydot0^2 * Clf(kappa)                   [load factor, g]
//
// with the approach parameter kappa = sin(tau)/sin(gamma0) * cos(tau + gamma0),
// the load-factor coefficient Clf(kappa) (eqs 25/27, precomputed in
// impact_clf_table.hpp), and the hull coefficient K0 = (alpha_hull/(W g^2))^{1/3}.
//
// Relative degree 2 via the elevator (theta->q->de); thrust (T->Tdot->Tddot, deg
// 3) drops out of the degree-2 HOCBF row, so the flare enforces it. Thrust still
// bounds impact load through the descent/sink-rate barrier (spec section 4).
//
// K0 and the local-affine model Clf(kappa) ~ Clf0 + dClf_dk*(kappa - kappa0) are
// FROZEN at the eval point by makeImpactLoadBarrier (mirrors makeDescentBarrier),
// so the templated barrier is smooth and needs no cbrt/pow in the Taylor path,
// while the attitude coupling stays live through kappa(tau, gamma0).
//
// NOTE: Clf is the LOAD-factor coefficient (C_l), NOT the aerodynamic max-lift
// coefficient CLmax (C_L) used by the descent barrier -- distinct quantity.
// =============================================================================
namespace autoland {

// Linear interpolation of the precomputed Clf(kappa) curve. Returns the value
// and the local slope d(Clf)/d(kappa); kappa is clamped to the table's practical
// range [0.2, 10] (value/slope held at the endpoint outside).
struct ClfLocal {
  double value{0};
  double slope{0};
};

inline ClfLocal clfLookup(double kappa) {
  const int n = kImpactClfN;
  const double* K = kImpactClfKappa;
  const double* V = kImpactClfVal;
  if (kappa <= K[0]) return {V[0], (V[1] - V[0]) / (K[1] - K[0])};
  if (kappa >= K[n - 1])
    return {V[n - 1], (V[n - 1] - V[n - 2]) / (K[n - 1] - K[n - 2])};
  int lo = 0, hi = n - 1;
  while (hi - lo > 1) {
    const int mid = (lo + hi) / 2;
    if (K[mid] <= kappa) lo = mid; else hi = mid;
  }
  const double slope = (V[hi] - V[lo]) / (K[hi] - K[lo]);
  return {V[lo] + slope * (kappa - K[lo]), slope};
}

struct ImpactLoadBarrier {
  double n_limit{0};
  double K0{0};                  // (alpha_hull / (W g^2))^(1/3), frozen
  double Clf0{0}, dClf_dk{0};    // local-affine load-factor coef about kappa0
  double kappa0{0};
  double Nb{0}, zs{0};           // Phi(z) = Nb (1 - exp(-z / zs))
  double tau_keel{0};            // tau = theta - tau_keel
  double eps_g0{0.02};           // smooth floor on sin(gamma0) (planing guard)

  template <class T>
  T operator()(const std::array<T, NXA>& X) const {
    using std::cos;
    using std::exp;
    using std::sin;
    using std::sqrt;
    const T tau = X[LTH] - tau_keel;
    const T gamma = X[LGAM];
    const T gamma0 = (-1.0) * gamma;                 // descent positive
    const T ydot0 = (-1.0) * X[LV] * sin(gamma);     // sink speed (>0 in descent)
    // Smooth floor on sin(gamma0): sg0 = sqrt(sin^2 + eps^2) >= eps, which caps
    // the planing singularity gamma0 -> 0 (kappa -> inf) while staying autodiff-
    // smooth. (No 1/0 in the Lie jet near touchdown as the sink rate decays.)
    const T sg0 = sqrt(sin(gamma0) * sin(gamma0) + eps_g0 * eps_g0);
    const T kappa = sin(tau) * cos(tau + gamma0) / sg0;
    const T Clf = Clf0 + dClf_dk * (kappa - kappa0);
    const T n_peak = (K0 * Clf) * (ydot0 * ydot0);
    const T Phi = Nb * (1.0 - exp(X[LH] * (-1.0 / zs)));
    return (n_limit - n_peak) + Phi;
  }
};

// Build the barrier, freezing K0 and the local Clf(kappa) model from the state
// (V, theta, gamma in X). beta is the hull dead-rise in radians.
inline ImpactLoadBarrier makeImpactLoadBarrier(const AeroLocal& a, double n_limit,
                                               double beta, double rho_water,
                                               double Nb, double zs,
                                               double tau_keel, double eps_g0,
                                               const LonStateVec& X) {
  ImpactLoadBarrier b;
  b.n_limit = n_limit;
  b.Nb = Nb;
  b.zs = zs;
  b.tau_keel = tau_keel;
  b.eps_g0 = eps_g0;

  const double tau0_raw = X[LTH] - tau_keel;
  const double gamma0_0 = -X[LGAM];
  // Frozen kappa0 (same smooth gamma0 floor as the templated path), clamped to
  // the table's practical range for the Clf lookup.
  const double sg0 =
      std::sqrt(std::sin(gamma0_0) * std::sin(gamma0_0) + eps_g0 * eps_g0);
  double kappa0 = std::sin(tau0_raw) * std::cos(tau0_raw + gamma0_0) / sg0;
  kappa0 = std::min(10.0, std::max(0.2, kappa0));
  b.kappa0 = kappa0;
  const ClfLocal c = clfLookup(kappa0);
  b.Clf0 = c.value;
  b.dClf_dk = c.slope;

  // Hull coefficient alpha_hull (eq 12a) with dead-rise / end-flow factors
  // (eqs 45/49). Clamp tau0 away from 0 and the phi_A = 0 singularity
  // (tan tau = 2 tan beta) so K0 stays finite and positive in pathological
  // states; landing trims (a few deg) sit well inside this clamp.
  const double deg = M_PI / 180.0;
  const double tau_min = 1.0 * deg;
  const double tau_max = std::atan(2.0 * std::tan(beta)) - 1.0 * deg;
  const double tau0 = std::min(tau_max, std::max(tau_min, tau0_raw));
  const double f_beta = M_PI / (2.0 * beta) - 1.0;                  // eq 45
  const double phi_A = 1.0 - std::tan(tau0) / (2.0 * std::tan(beta));  // eq 49 (>0)
  const double alpha_hull = f_beta * f_beta * phi_A * rho_water * M_PI /
                            (6.0 * std::sin(tau0) * std::cos(tau0) * std::cos(tau0));
  const double W = a.mass * a.g;
  b.K0 = std::cbrt(alpha_hull / (W * a.g * a.g));
  return b;
}

}  // namespace autoland
