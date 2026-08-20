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
// MULTI-FLOAT (n_surfaces): TN 1516 is a SINGLE planing surface carrying the full
// weight (float OR hull). For a twin-float airframe the load is shared, so each
// float takes W/n_surfaces; we evaluate the single-surface theory per float with
// W = mass*g/n_surfaces. Since n_peak ~ (W)^{-1/3}, n_surfaces=2 raises the CG load
// factor by 2^{1/3} ~ 1.26x over full weight on one surface -- the conservative
// (safe) direction for symmetric contact. n_surfaces=1 (default) = single hull.
// See documentation/impact_load_barrier_spec.md sections 7-8.
//
// Relative degree 2 via the elevator (theta->q->de); thrust (T->Tdot->Tddot, deg
// 3) drops out of the degree-2 HOCBF row, so the flare enforces it. Thrust still
// bounds impact load through the descent/sink-rate barrier (spec section 4).
//
// K0(tau) and Clf(kappa) are FROZEN at the eval point by makeImpactLoadBarrier
// as LOCAL AFFINE models (value + first derivative), so the templated barrier
// is smooth and needs no cbrt/pow in the Taylor path while its FIRST-ORDER
// state sensitivities are exact. Freezing K0's VALUE alone is not enough: at
// landing trims d(ln K0)/dtau ~ -(1/3)cot(tau) partially cancels the retained
// kappa channel (elasticity ~0.6 * cot(tau)), so a value-only freeze inflates
// d(n_peak)/d(theta) -- and with it the hard row's elevator coefficient --
// by ~2x. See test "impact barrier frozen model preserves the true attitude
// gradient".
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

// Hull-coefficient model K0(tau) = (alpha_hull(tau)/(W g^2))^(1/3) (eqs 12a/
// 45/49), evaluated at tau clamped away from 0 and the phi_A = 0 singularity
// (tan tau = 2 tan beta) so K0 stays finite and positive in pathological
// states; landing trims (a few deg) sit well inside the clamp. Shared by the
// frozen-model factory below and the exact evaluator impactNPeakExact.
struct ImpactK0Local {
  double K0{0};
  double tau0{0};   // the clamped tau actually used
  double phi_A{0};  // end-flow factor at tau0 (the affine slope needs it)
};

inline ImpactK0Local impactK0(double mass, double g, double beta,
                              double rho_water, double tau_raw) {
  const double deg = M_PI / 180.0;
  const double tau_min = 1.0 * deg;
  const double tau_max = std::atan(2.0 * std::tan(beta)) - 1.0 * deg;
  const double tau0 = std::min(tau_max, std::max(tau_min, tau_raw));
  const double f_beta = M_PI / (2.0 * beta) - 1.0;                     // eq 45
  const double phi_A = 1.0 - std::tan(tau0) / (2.0 * std::tan(beta));  // eq 49 (>0)
  const double alpha_hull = f_beta * f_beta * phi_A * rho_water * M_PI /
                            (6.0 * std::sin(tau0) * std::cos(tau0) * std::cos(tau0));
  const double W = mass * g;
  return {std::cbrt(alpha_hull / (W * g * g)), tau0, phi_A};
}

// Exact (no freeze) TN 1516 peak load n_peak = K0(tau) Clf(kappa) ydot0^2 at
// an arbitrary contact geometry -- plain double math for plant-side truth
// diagnostics (e.g. the wave-slope-referenced touchdown load, where tau and
// gamma0 are tilted by the local surface angle and ydot0 is the closure rate
// onto the moving surface). Mirrors the factory's clamps: kappa to the Clf
// table range, tau via impactK0.
inline double impactNPeakExact(double mass, double g, double beta,
                               double rho_water, double eps_g0, double tau,
                               double gamma0, double ydot0, int n_surfaces = 1) {
  const double sg0 =
      std::sqrt(std::sin(gamma0) * std::sin(gamma0) + eps_g0 * eps_g0);
  double kappa = std::sin(tau) * std::cos(tau + gamma0) / sg0;
  kappa = std::min(10.0, std::max(0.2, kappa));
  const double Clf = clfLookup(kappa).value;
  // W/n_surfaces per float (single-surface theory applied per planing surface).
  return impactK0(mass / n_surfaces, g, beta, rho_water, tau).K0 * Clf * ydot0 *
         ydot0;
}

struct ImpactLoadBarrier {
  double n_limit{0};
  double K00{0}, dK0_dtau{0};    // local-affine hull coef (alpha_hull/(W g^2))^(1/3)
  double tau0{0};                //   about the clamped anchor tau0
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
    const T K0t = K00 + dK0_dtau * (tau - tau0);
    const T n_peak = (K0t * Clf) * (ydot0 * ydot0);
    const T Phi = Nb * (1.0 - exp(X[LH] * (-1.0 / zs)));
    return (n_limit - n_peak) + Phi;
  }
};

// Build the barrier, freezing the local-affine K0(tau) and Clf(kappa) models
// from the state (V, theta, gamma in X). beta is the hull dead-rise in radians.
inline ImpactLoadBarrier makeImpactLoadBarrier(const AeroLocal& a, double n_limit,
                                               double beta, double rho_water,
                                               double Nb, double zs,
                                               double tau_keel, double eps_g0,
                                               const LonStateVec& X,
                                               int n_surfaces = 1) {
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

  // Hull coefficient K0 at the clamped anchor (impactK0: eqs 12a/45/49 with
  // the tau clamps). Effective weight is W/n_surfaces (W/2 per float).
  const ImpactK0Local kl =
      impactK0(a.mass / n_surfaces, a.g, beta, rho_water, tau0_raw);
  const double tau0 = kl.tau0;
  const double phi_A = kl.phi_A;
  b.K00 = kl.K0;
  // Local-affine slope, dK0/dtau = (K0/3) d(ln alpha_hull)/dtau, anchored at the
  // clamped tau0 (bounds |dK0_dtau| via the cot term; mirrors the kappa0 clamp).
  // d(ln alpha_hull)/dtau = -cot(tau) + 2 tan(tau) - sec^2(tau)/(2 tan(beta) phi_A)
  // from -ln sin(tau), -2 ln cos(tau), and ln phi_A. The cbrt/tan stay here in
  // the factory (double math), out of the Taylor path.
  const double ct0 = std::cos(tau0);
  const double dlnA = -1.0 / std::tan(tau0) + 2.0 * std::tan(tau0) -
                      1.0 / (ct0 * ct0 * 2.0 * std::tan(beta) * phi_A);
  b.dK0_dtau = b.K00 * dlnA / 3.0;
  b.tau0 = tau0;
  return b;
}

}  // namespace autoland
