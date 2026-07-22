#pragma once
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

// =============================================================================
// Linear (Airy) surface-wave field for the water-landing sim -- long-crested,
// deep-water, evaluated in the longitudinal vertical plane as eta(x, t).
//
// SPECTRUM (which sea to make): the JONSWAP form (Hasselmann et al. 1973; ITTC
// Recommended Procedure 7.5-02-07-01.1),
//
//   S(w)   = A_gamma S_B(w) gamma^b
//   S_B(w) = (5/16) Hs^2 wp^4 / w^5 exp(-1.25 (wp/w)^4)          [Bretschneider]
//   b      = exp(-(w - wp)^2 / (2 sigma^2 wp^2)),  sigma = 0.07 (w <= wp)
//                                                        = 0.09 (w >  wp)
//   A_gamma = 1 - 0.287 ln(gamma)     [Hs-preserving normalization, DNV-RP-C205]
//
// gamma = 1 recovers EXACTLY the two-parameter Bretschneider spectrum that
// NATO STANAG 4194 prescribes with its open-ocean sea-state table (Annex D
// most-probable pairs: SS2 Hs 0.3 m / Tm 6.3 s, SS3 0.9 / 7.5, SS4 1.9 / 8.8,
// SS5 3.3 / 9.7, SS6 5.0 / 12.4). gamma ~ 3.3 is the standard fetch-limited
// (developing) sea -- the LAKE case, where a developed sea cannot exist. Pick
// lake (Hs, Tp) from wind speed + fetch with the USACE Coastal Engineering
// Manual growth laws (EM 1110-2-1100 Part II-2), e.g. U10 = 8 m/s over 3 km
// of fetch -> Hs ~ 0.22 m, Tp ~ 1.8 s.
//
// REALIZATION (how to make eta(x, t) from S): the standard superposition of N
// harmonics (St. Denis & Pierson 1953; Fossen, Handbook of Marine Craft
// Hydrodynamics and Motion Control, ch. 8):
//
//   eta(x, t) = sum_i A_i cos(k_i x - dir w_i t + phi_i)
//   A_i   = sqrt(2 S(w_i) dw)            [component amplitude from the PSD]
//   k_i   = w_i^2 / g                    [deep-water Airy dispersion]
//   phi_i ~ U[0, 2pi)                    [random phase; seeded, so a scenario
//                                         is ONE locked, reproducible sea]
//
// with each w_i drawn uniformly WITHIN its frequency bin (Fossen's recipe) so
// the realization is non-periodic. dir = -1 makes the waves travel toward -x:
// HEAD seas for an aircraft flying +x, the high-encounter-rate landing case
// (encounter frequency w_e = w (1 + V/c_p)); dir = +1 is following seas. A
// single REGULAR wave -- deterministic demos and sharp unit tests -- is the
// N = 1 degenerate case with A = H/2 (cf. the gust model's step-gust limit).
//
// The field is a PLANT-side truth only (wave-blind filter): touchdown happens
// at h = eta, the radar altimeter measures clearance h - eta, and the touchdown
// diagnostics evaluate the TN 1516 slam load against the local wave slope --
// but the nominal controller and the CBF-QP keep the flat-water model, so wave
// runs probe the safety filter against an unmodeled seaway (TODO.md
// "smooth-water assumption"). This header is pure field evaluation and owns no
// sim state (cf. wind_gust.hpp).
// =============================================================================
namespace autoland {

// Deep-water (Airy) dispersion.
inline double waveNumberDeep(double w, double g) { return w * w / g; }

// One-sided JONSWAP spectral density S(w) [m^2 s]; gamma <= 1 degenerates to
// the two-parameter Bretschneider form (the STANAG 4194 open-ocean spectrum).
inline double waveSpectrum(double w, double Hs, double wp, double gamma) {
  if (w <= 0.0) return 0.0;
  const double r = wp / w;
  const double r4 = (r * r) * (r * r);
  const double SB = (5.0 / 16.0) * Hs * Hs * r4 / w * std::exp(-1.25 * r4);
  if (gamma <= 1.0) return SB;
  const double sigma = (w <= wp) ? 0.07 : 0.09;
  const double dv = (w - wp) / (sigma * wp);
  const double A_gamma = 1.0 - 0.287 * std::log(gamma);
  return A_gamma * SB * std::pow(gamma, std::exp(-0.5 * dv * dv));
}

// Scenario parameters (YAML `waves:` section of the lon scenario).
struct WaveConfig {
  bool enabled{false};
  bool regular{false};    // N = 1 deterministic component; Hs is then the
                          // crest-to-trough height H (A = H/2)
  double Hs{0.22};        // significant wave height [m]
  double Tp{1.8};         // peak (modal) period [s]
  double gamma{3.3};      // JONSWAP peak enhancement; 1 = Bretschneider/STANAG
  int n{128};             // spectral components
  unsigned seed{1};       // phase/frequency seed (one locked sea per scenario)
  double dir{-1.0};       // -1 = head seas (waves travel -x), +1 = following
  double phase_deg{0.0};  // regular wave: phase at (x, t) = (0, 0)
  double g{9.80665};
  double w_lo_fac{0.5};   // realized band [w_lo_fac, w_hi_fac] * wp; the
  double w_hi_fac{5.0};   // Bretschneider tail outside it holds < 0.3% of m0
  double contact_len{0.4};  // hull wetted length [m] for slopeMean: the slam
                            // geometry sees the surface tilt ACROSS the keel,
                            // not the point slope (whose spectrum k^2 S(w) is
                            // dominated by sub-hull-length ripples)
};

struct WaveComponent {
  double A, w, k, phi;
};

// Realized long-crested field. eta / slope / etaDot are the elevation and its
// exact spatial (d/dx) and temporal (d/dt at fixed x) derivatives; a single
// component satisfies the advection identity etaDot = -c_p * slope with phase
// speed c_p = dir w / k.
struct WaveField {
  bool enabled{false};
  double dir{-1.0};
  std::vector<WaveComponent> comps;

  double eta(double x, double t) const {
    if (!enabled) return 0.0;
    double e = 0.0;
    for (const WaveComponent& c : comps)
      e += c.A * std::cos(c.k * x - dir * c.w * t + c.phi);
    return e;
  }
  double slope(double x, double t) const {  // d eta / dx
    if (!enabled) return 0.0;
    double s = 0.0;
    for (const WaveComponent& c : comps)
      s -= c.A * c.k * std::sin(c.k * x - dir * c.w * t + c.phi);
    return s;
  }
  double etaDot(double x, double t) const {  // d eta / dt at fixed x
    if (!enabled) return 0.0;
    double s = 0.0;
    for (const WaveComponent& c : comps)
      s += c.A * dir * c.w * std::sin(c.k * x - dir * c.w * t + c.phi);
    return s;
  }
  // Mean surface slope across a contact of length L centred on x -- the tilt a
  // hull of that wetted length actually meets (chord slope; degenerates to the
  // point slope for L <= 0). Components shorter than L average out, exactly as
  // the keel bridges them.
  double slopeMean(double x, double t, double L) const {
    if (!enabled) return 0.0;
    if (L <= 0.0) return slope(x, t);
    return (eta(x + 0.5 * L, t) - eta(x - 0.5 * L, t)) / L;
  }
};

inline WaveField makeWaveField(const WaveConfig& c) {
  WaveField f;
  f.enabled = c.enabled;
  f.dir = c.dir;
  if (!c.enabled) return f;
  const double wp = 2.0 * M_PI / c.Tp;
  if (c.regular) {
    f.comps.push_back(
        {0.5 * c.Hs, wp, waveNumberDeep(wp, c.g), c.phase_deg * M_PI / 180.0});
    return f;
  }
  std::mt19937 rng(c.seed);
  std::uniform_real_distribution<double> uphi(0.0, 2.0 * M_PI);
  std::uniform_real_distribution<double> ubin(0.02, 0.98);
  const double wlo = c.w_lo_fac * wp;
  const int n = std::max(1, c.n);
  const double dw = (c.w_hi_fac * wp - wlo) / n;
  f.comps.reserve(n);
  for (int i = 0; i < n; ++i) {
    const double w = wlo + (i + ubin(rng)) * dw;
    const double A = std::sqrt(2.0 * waveSpectrum(w, c.Hs, wp, c.gamma) * dw);
    f.comps.push_back({A, w, waveNumberDeep(w, c.g), uphi(rng)});
  }
  return f;
}

}  // namespace autoland
