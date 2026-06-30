#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "autoland/aero_table.hpp"
#include "autoland/config.hpp"
#include "autoland/lon_augmented.hpp"
#include "autoland/mixing.hpp"
#include "autoland/stall_model.hpp"

using namespace autoland;
using Catch::Approx;

namespace {
const std::string kData = AUTOLAND_DATA_DIR;
constexpr double kDeg = M_PI / 180.0;
}  // namespace

TEST_CASE("stallLookup: inert below the grid, held flat above", "[stall]") {
  // Far below the table (attached flow): blend weight 0, zero slope -- so cruise
  // is pure VSPAERO and the Lie derivatives pick up no spurious stall gradient.
  const StallBlend lo = stallLookup(-30.0 * kDeg);
  CHECK(lo.w == 0.0);
  CHECK(lo.w_da == 0.0);
  // Far above the table (deep stall): fully blended, held flat.
  const StallBlend hi = stallLookup(120.0 * kDeg);
  CHECK(hi.w == Approx(1.0));
  CHECK(hi.w_da == 0.0);
  CHECK(hi.CLpost == Approx(0.0).margin(1e-6));  // lift craters to ~0 at 90 deg
}

TEST_CASE("stallLookup: attached inert, blend ramps in at stall", "[stall]") {
  CHECK(stallLookup(2.0 * kDeg).w == Approx(0.0).margin(1e-12));   // attached
  CHECK(stallLookup(20.0 * kDeg).w == Approx(1.0));                // fully stalled
  // Post-stall lift declines toward zero with alpha (Viterna flat-plate tail).
  const double cl16 = stallLookup(16.0 * kDeg).CLpost;
  const double cl45 = stallLookup(45.0 * kDeg).CLpost;
  const double cl80 = stallLookup(80.0 * kDeg).CLpost;
  CHECK(cl16 > cl45);
  CHECK(cl45 > cl80);
  CHECK(cl80 < 0.4);
  // Separation drag rises toward CD_max.
  CHECK(stallLookup(60.0 * kDeg).CDpost > stallLookup(20.0 * kDeg).CDpost);
}

TEST_CASE("stallLookup interpolates linearly between grid nodes", "[stall]") {
  // Grid is 0.5 deg through stall; 16.25 deg is the midpoint of [16.0, 16.5].
  const double w1 = stallLookup(16.0 * kDeg).w;
  const double w2 = stallLookup(16.5 * kDeg).w;
  const StallBlend mid = stallLookup(16.25 * kDeg);
  CHECK(mid.w == Approx(0.5 * (w1 + w2)));
  CHECK(mid.w_da == Approx((w2 - w1) / (0.5 * kDeg)));
}

TEST_CASE("stall blend: OFF is a no-op, ON drops lift and raises drag", "[stall]") {
  AeroTable table(AeroTable::fromFile(kData + "/example.stab"));
  AircraftConfig cfg = loadAircraftConfig(kData + "/aircraft.yaml");
  const Mixing mx = Mixing::build(cfg, table);
  const double V = 16.0, alpha = 30.0 * kDeg;  // deep post-stall (w = 1)

  cfg.stall.enabled = false;
  const AeroLocal off = makeAeroLocal(table, mx, cfg, V, alpha);
  CHECK(off.stall_on == false);
  CHECK(off.off_w == 0.0);

  cfg.stall.enabled = true;
  const AeroLocal on = makeAeroLocal(table, mx, cfg, V, alpha);
  CHECK(on.stall_on == true);
  // Fully stalled: the frozen affine evaluates to w = 1 and CL = the absolute
  // post-stall (Viterna) value at this alpha, independent of the (inviscid) deck.
  CHECK(on.off_w + on.dAlpha_w * alpha == Approx(1.0));
  // Fully stalled: CL/CD hand off to the absolute Viterna post-stall values
  // (deck-independent), with substantial separation drag.
  CHECK(on.off_CLp + on.dAlpha_CLp * alpha == Approx(stallLookup(alpha).CLpost));
  CHECK(on.off_CDp + on.dAlpha_CDp * alpha == Approx(stallLookup(alpha).CDpost));
  CHECK(stallLookup(alpha).CDpost > 0.2);

  // The stalled (Viterna) lift is far below the inviscid lift, so flight-path
  // acceleration drops -- this is robust regardless of the deck's extrapolation.
  LonStateVec X;
  X << 30.0, V, 0.0, alpha, 0.0, 5.0, 0.0;  // theta = alpha, gamma = 0
  const auto Xa = toArray(X);
  CHECK(LonDrift(on)(Xa)[LGAM] < LonDrift(off)(Xa)[LGAM]);
}

TEST_CASE("stall blend leaves attached flow unchanged", "[stall]") {
  AeroTable table(AeroTable::fromFile(kData + "/example.stab"));
  AircraftConfig cfg = loadAircraftConfig(kData + "/aircraft.yaml");
  const Mixing mx = Mixing::build(cfg, table);
  const double V = 18.0, alpha = 2.0 * kDeg;  // attached flow (w = 0)

  cfg.stall.enabled = false;
  const AeroLocal off = makeAeroLocal(table, mx, cfg, V, alpha);
  cfg.stall.enabled = true;
  const AeroLocal on = makeAeroLocal(table, mx, cfg, V, alpha);

  LonStateVec X;
  X << 30.0, V, 0.0, alpha, 0.0, 5.0, 0.0;
  const auto Xa = toArray(X);
  const auto dOff = LonDrift(off)(Xa);
  const auto dOn = LonDrift(on)(Xa);
  // w = 0 in attached flow, so the dynamics are bit-identical to stall-off.
  CHECK(dOn[LGAM] == Approx(dOff[LGAM]).margin(1e-12));
  CHECK(dOn[LV] == Approx(dOff[LV]).margin(1e-12));
}
