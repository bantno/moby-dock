#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>

#include "autoland/beaver_lon.hpp"
#include "autoland/impact_barrier.hpp"
#include "autoland/lie_taylor.hpp"
#include "autoland/lon_augmented.hpp"

using namespace autoland;
using Catch::Approx;

namespace {
constexpr double kDeg = M_PI / 180.0;

// A descending approach state [h, V, gamma, theta, q, P, -] for the Beaver.
LonStateVec approachState(double V, double P) {
  LonStateVec X;
  const double gamma = -3.0 * kDeg;
  const double alpha = 5.0 * kDeg;         // representative approach AoA
  X << 5.0, V, gamma, gamma + alpha, 0.0, P, 0.0;
  return X;
}

// Impact-load barrier frozen at X, using the Beaver mass (n_surfaces = 2, W/2 per
// float). makeImpactLoadBarrier only reads a.mass / a.g from the AeroLocal.
ImpactLoadBarrier beaverImpact(const LonStateVec& X) {
  AeroLocal a;
  a.mass = BeaverGeom::mass;
  a.g = 9.81;
  return makeImpactLoadBarrier(a, /*n_limit=*/3.0, /*beta=*/22.5 * kDeg,
                               /*rho_water=*/1000.0, /*Nb=*/10.0, /*zs=*/2.0,
                               /*tau_keel=*/0.0, /*eps_g0=*/0.02, X,
                               /*n_surfaces=*/2);
}
}  // namespace

// THE point of the single-integrator-power model: BOTH the elevator (moment
// channel) and the power rate u_P (force channel) reach the impact barrier at
// relative degree EXACTLY 2 -- a uniform degree-2, two-input HOCBF row. Contrast
// the AHAB double-integrator thrust, which is degree >= 3 and drops out of the
// degree-2 impact row (test "impact barrier is relative degree 2 via the
// elevator"). Derivatives are the exact flow-Taylor Lie jet (no finite diffs).
TEST_CASE("Beaver power reaches the impact barrier at degree 2 (two-actuator row)",
          "[beaver]") {
  const double rho = 1.112, V = 33.0, P = 50.0;
  const LonStateVec X = approachState(V, P);
  const BeaverLonDrift f{rho, 9.81, 0.0, BeaverAeroCoef{}};
  const ImpactLoadBarrier b = beaverImpact(X);
  const auto x0 = toArray(X);
  const auto g_de = beaverGColDe(rho, V);
  const auto g_uP = beaverGColPower();

  // Degree 1: neither control appears (the impact barrier depends on
  // theta/gamma/V/h, not on q or P) -> relative degree is at least 2.
  CHECK(std::abs(lieAlong<1, NXA>(f, b, x0, g_de)) < 1e-9);
  CHECK(std::abs(lieAlong<1, NXA>(f, b, x0, g_uP)) < 1e-9);

  // Degree 2: BOTH controls have authority. The elevator via the pitch moment
  // (Cm_de -> q_dot), the power rate via the force channel (dpt -> V_dot/gamma_dot
  // -> sink). L_g L_f b has two nonzero columns -> the mixed-actuator row.
  const double LgLf_de = lieAlong<2, NXA>(f, b, x0, g_de);
  const double LgLf_uP = lieAlong<2, NXA>(f, b, x0, g_uP);
  CHECK(std::abs(LgLf_de) > 1e-6);  // elevator: degree 2 (unchanged)
  CHECK(std::abs(LgLf_uP) > 1e-6);  // POWER: degree 2 (the new capability)
}

// Complementary authority vs. speed: the elevator's moment authority scales like
// rho V^2 (weak at low speed), while the power's force authority scales like ~1/V
// through dpt ~ 1/V^3 (strong at low speed). So the power-to-elevator authority
// ratio is larger at low speed -- power backs up the failing elevator exactly at
// the low-q touchdown.
TEST_CASE("Beaver power authority on the impact barrier grows relative to the "
          "elevator as speed drops",
          "[beaver]") {
  const double rho = 1.112, P = 50.0;
  auto ratio = [&](double V) {
    const LonStateVec X = approachState(V, P);
    const BeaverLonDrift f{rho, 9.81, 0.0, BeaverAeroCoef{}};
    const ImpactLoadBarrier b = beaverImpact(X);
    const auto x0 = toArray(X);
    const double de = lieAlong<2, NXA>(f, b, x0, beaverGColDe(rho, V));
    const double up = lieAlong<2, NXA>(f, b, x0, beaverGColPower());
    return std::abs(up) / std::abs(de);
  };
  const double r_slow = ratio(25.0);
  const double r_fast = ratio(40.0);
  CHECK(r_slow > r_fast);  // power relatively stronger at low speed
}
