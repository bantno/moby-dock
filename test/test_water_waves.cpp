#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <string>

#include "autoland/lon_sim.hpp"
#include "autoland/water_waves.hpp"

using namespace autoland;
using Catch::Approx;

namespace {
const std::string kData = AUTOLAND_DATA_DIR;
constexpr double kDeg = M_PI / 180.0;
constexpr double kG = 9.80665;

// Trapezoid integral of S(w) over [wlo, whi] -- the realized variance target.
double spectrumM0(double Hs, double wp, double gamma, double wlo, double whi,
                  int n) {
  const double dw = (whi - wlo) / n;
  double m0 = 0.0;
  for (int i = 0; i <= n; ++i) {
    const double wgt = (i == 0 || i == n) ? 0.5 : 1.0;
    m0 += wgt * waveSpectrum(wlo + i * dw, Hs, wp, gamma) * dw;
  }
  return m0;
}
}  // namespace

// The spectrum layer: the Bretschneider closed form (the STANAG 4194
// two-parameter open-ocean spectrum), the JONSWAP gamma = 1 degeneration, and
// the zeroth-moment identity m0 = integral S dw = Hs^2 / 16 that ties the
// spectral density to the significant wave height it claims to represent.
TEST_CASE("wave spectrum: Bretschneider closed form and Hs variance identity",
          "[waves]") {
  const double Hs = 0.9, Tp = 7.5;  // STANAG 4194 Annex D sea state 3
  const double wp = 2.0 * M_PI / Tp;

  // Closed form at the modal frequency: S(wp) = (5/16) Hs^2 / wp * e^{-5/4}.
  CHECK(waveSpectrum(wp, Hs, wp, 1.0) ==
        Approx((5.0 / 16.0) * Hs * Hs / wp * std::exp(-1.25)).epsilon(1e-12));

  // gamma = 1 IS Bretschneider, everywhere on the band.
  for (double f = 0.4; f <= 4.0; f += 0.15) {
    const double w = f * wp;
    const double r4 = std::pow(wp / w, 4.0);
    const double SB =
        (5.0 / 16.0) * Hs * Hs * r4 / w * std::exp(-1.25 * r4);
    CHECK(waveSpectrum(w, Hs, wp, 1.0) == Approx(SB).epsilon(1e-12));
  }

  // m0 = Hs^2/16 for the Bretschneider form (exact); the JONSWAP
  // A_gamma = 1 - 0.287 ln(gamma) normalization holds it only approximately.
  const double m0_target = Hs * Hs / 16.0;
  CHECK(spectrumM0(Hs, wp, 1.0, 0.05 * wp, 12.0 * wp, 20000) ==
        Approx(m0_target).epsilon(0.01));
  CHECK(spectrumM0(Hs, wp, 3.3, 0.05 * wp, 12.0 * wp, 20000) ==
        Approx(m0_target).epsilon(0.06));

  // The peak-enhancement factor only ADDS energy near wp (before
  // renormalization): far from the peak JONSWAP ~ A_gamma * Bretschneider.
  const double far = 3.0 * wp;
  CHECK(waveSpectrum(far, Hs, wp, 3.3) ==
        Approx((1.0 - 0.287 * std::log(3.3)) * waveSpectrum(far, Hs, wp, 1.0))
            .epsilon(1e-6));
}

// The realization layer: component amplitudes carry the spectrum's variance
// (sum A_i^2/2 = integral S dw over the realized band), every component obeys
// the deep-water dispersion, and the seeded draw is reproducible.
TEST_CASE("wave realization: variance, dispersion, seeded determinism",
          "[waves]") {
  WaveConfig c;
  c.enabled = true;
  c.Hs = 0.22;
  c.Tp = 1.8;
  c.gamma = 3.3;
  c.n = 128;
  c.seed = 1;
  const double wp = 2.0 * M_PI / c.Tp;

  const WaveField f = makeWaveField(c);
  REQUIRE(static_cast<int>(f.comps.size()) == c.n);

  double var = 0.0;
  for (const WaveComponent& wc : f.comps) {
    var += 0.5 * wc.A * wc.A;
    CHECK(wc.k == Approx(wc.w * wc.w / c.g).epsilon(1e-12));  // Airy deep water
    CHECK(wc.w >= c.w_lo_fac * wp);
    CHECK(wc.w <= c.w_hi_fac * wp);
  }
  // Band-limited target (the within-bin frequency jitter adds a few % of
  // sampling scatter on the peaked gamma = 3.3 spectrum).
  const double m0_band =
      spectrumM0(c.Hs, wp, c.gamma, c.w_lo_fac * wp, c.w_hi_fac * wp, 20000);
  CHECK(var == Approx(m0_band).epsilon(0.08));

  // Same seed -> the same locked sea; a different seed -> a different one.
  const WaveField f2 = makeWaveField(c);
  CHECK(f2.comps[17].phi == f.comps[17].phi);
  CHECK(f2.comps[17].w == f.comps[17].w);
  c.seed = 2;
  const WaveField f3 = makeWaveField(c);
  CHECK(f3.comps[17].phi != f.comps[17].phi);

  // Statistical sanity of the field itself: the time variance of eta at a
  // fixed point over many peak periods approaches sum A_i^2/2 (loose bound --
  // one realization, finite window).
  double m = 0.0, m2 = 0.0;
  const int N = 20000;
  const double T = 400.0 * c.Tp;
  for (int i = 0; i < N; ++i) {
    const double e = f.eta(0.0, T * i / N);
    m += e;
    m2 += e * e;
  }
  m /= N;
  m2 = m2 / N - m * m;
  CHECK(m2 == Approx(var).epsilon(0.25));
}

// The regular (N = 1) degenerate case: exact Airy kinematics -- elevation,
// spatial slope, heave rate, the advection identity, and travel direction.
TEST_CASE("regular wave: analytic elevation, slope, advection, direction",
          "[waves]") {
  WaveConfig c;
  c.enabled = true;
  c.regular = true;
  c.Hs = 0.20;  // crest-to-trough H -> amplitude A = 0.10
  c.Tp = 2.0;
  c.phase_deg = 30.0;
  c.dir = -1.0;  // head seas: eta = A cos(k x + w t + phi)
  const WaveField f = makeWaveField(c);
  REQUIRE(f.comps.size() == 1);

  const double A = 0.10, w = 2.0 * M_PI / c.Tp, k = w * w / kG;
  const double phi = 30.0 * kDeg;
  CHECK(f.comps[0].A == Approx(A).epsilon(1e-12));
  CHECK(f.comps[0].k == Approx(k).epsilon(1e-12));

  for (double x : {0.0, 1.7, 4.2}) {
    for (double t : {0.0, 0.3, 1.9}) {
      const double th = k * x + w * t + phi;
      CHECK(f.eta(x, t) == Approx(A * std::cos(th)).margin(1e-12));
      CHECK(f.slope(x, t) == Approx(-A * k * std::sin(th)).margin(1e-12));
      CHECK(f.etaDot(x, t) == Approx(-A * w * std::sin(th)).margin(1e-12));
      // Advection: eta_t = -c_p eta_x with phase speed c_p = dir w / k.
      CHECK(f.etaDot(x, t) ==
            Approx(-(c.dir * w / k) * f.slope(x, t)).margin(1e-12));
    }
  }

  // Direction: riding the crest. Following seas (dir = +1) keep eta constant
  // along x = (w/k) t; head seas along x = -(w/k) t.
  WaveConfig cf = c;
  cf.dir = 1.0;
  const WaveField ff = makeWaveField(cf);
  const double cp = w / k;
  for (double t : {0.0, 0.7, 2.3}) {
    CHECK(ff.eta(cp * t, t) == Approx(ff.eta(0.0, 0.0)).margin(1e-12));
    CHECK(f.eta(-cp * t, t) == Approx(f.eta(0.0, 0.0)).margin(1e-12));
  }

  // Disabled field is identically flat water.
  WaveConfig off;
  const WaveField f0 = makeWaveField(off);
  CHECK(f0.eta(3.0, 5.0) == 0.0);
  CHECK(f0.slope(3.0, 5.0) == 0.0);
  CHECK(f0.etaDot(3.0, 5.0) == 0.0);
}

// End-to-end: the flagship approach landed onto a locked regular wave
// (data/lon_landing_wave_regular.yaml). Touchdown must happen AT the
// instantaneous surface -- not at h = 0 -- with the wave-blind filter flying
// the surface-relative altimeter, and the touchdown record must carry the
// wave-referenced TN 1516 slam-load truth alongside its flat-water twin.
TEST_CASE("lon sim touches down on the instantaneous wave surface", "[waves]") {
  LonSim sim(kData + "/AHAB_combined.stab", kData + "/aircraft.yaml",
             kData + "/lon_landing_wave_regular.yaml");
  REQUIRE(sim.scenario().waves.enabled);
  REQUIRE(sim.scenario().waves.regular);

  const LonTouchdown td = sim.run("test_wave_landing_trace.csv");

  REQUIRE(td.reached);
  // Contact is at the surface: h_td = eta_td within one integration step of
  // closure (NOT at the flat-water h = 0).
  CHECK(std::abs(td.h - td.eta) < 0.02);
  // The locked wave phase at contact is nonzero -- the test would be vacuous
  // if the keel happened to meet a zero crossing.
  CHECK(std::abs(td.eta) > 0.005);
  // Envelope (loose): the surface-relative altimeter still lands it softly.
  CHECK(td.sink < 1.0);
  CHECK(td.sink_rel < 1.5);
  CHECK(td.V <= 14.0 + 1e-6);
  // Slam-load truth populated, and the wave-referenced value differs from the
  // flat-water one (the smooth-water assumption is measurably wrong here).
  CHECK(td.n_peak_flat > 0.0);
  CHECK(td.n_peak_wave > 0.0);
  CHECK(td.n_peak_wave != Approx(td.n_peak_flat).epsilon(0.01));
}
