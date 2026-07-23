#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <string>

#include "autoland/aero_table.hpp"
#include "autoland/config.hpp"
#include "autoland/lon_sim.hpp"
#include "autoland/mixing.hpp"
#include "autoland/wind_gust.hpp"

using namespace autoland;
using Catch::Approx;

namespace {
const std::string kData = AUTOLAND_DATA_DIR;
}

// ---- 1-cosine shape (MIL-F-8785C / Simulink Discrete Wind Gust Model) -------

TEST_CASE("Gust axis reproduces the 1-cosine shape", "[wind_gust]") {
  const double Vm = 4.0, dm = 120.0;
  GustAxis g{Vm, dm};

  // Region boundaries of the piecewise definition.
  CHECK(g.value(-10.0) == 0.0);
  CHECK(g.value(0.0) == 0.0);
  CHECK(g.value(dm) == Approx(Vm));
  CHECK(g.value(10.0 * dm) == Vm);

  // Closed-form interior points: half amplitude at dm/2, quarter-cosine values.
  CHECK(g.value(0.5 * dm) == Approx(0.5 * Vm));
  CHECK(g.value(0.25 * dm) == Approx(0.5 * Vm * (1.0 - std::cos(M_PI * 0.25))));

  // Ramp is monotonic for positive amplitude.
  double prev = 0.0;
  for (int i = 1; i <= 32; ++i) {
    const double v = g.value(dm * i / 32.0);
    CHECK(v >= prev);
    prev = v;
  }

  // C1 join: slope vanishes at both ends, peaks at midpoint with Vm*pi/(2 dm).
  CHECK(g.slope(0.0) == 0.0);
  CHECK(g.slope(dm) == 0.0);
  CHECK(g.slope(0.5 * dm) == Approx(Vm * M_PI / (2.0 * dm)));
}

TEST_CASE("Gust axis handles signed amplitude and degenerate length", "[wind_gust]") {
  // Headwind gust = negative amplitude, same shape mirrored.
  GustAxis head{-3.0, 60.0};
  CHECK(head.value(30.0) == Approx(-1.5));
  CHECK(head.value(60.0) == Approx(-3.0));
  CHECK(head.slope(30.0) < 0.0);

  // dm <= 0 degenerates to a step gust rather than dividing by zero.
  GustAxis step{2.0, 0.0};
  CHECK(step.value(0.0) == 0.0);
  CHECK(step.value(1e-9) == 2.0);
  CHECK(step.slope(1.0) == 0.0);
}

TEST_CASE("Gust model gates on start time and enabled flag", "[wind_gust]") {
  DiscreteGustConfig g;
  g.t_start = 5.0;
  g.u = {4.0, 100.0};
  g.w = {-1.0, 50.0};

  // Disabled: no penetration, no wind, regardless of t and x.
  g.enabled = false;
  CHECK(gustXdot(g, 10.0, 18.0) == 0.0);
  CHECK(gustWind(g, 25.0).u == 0.0);
  CHECK(gustWindRate(g, 25.0, 18.0).w == 0.0);

  // Enabled: penetration integrates the airspeed only after onset.
  g.enabled = true;
  CHECK(gustXdot(g, 4.99, 18.0) == 0.0);
  CHECK(gustXdot(g, 5.0, 18.0) == 18.0);

  // Wind vector and its time rate follow the per-axis shape/slope at x.
  const double x = 25.0, V = 18.0;
  CHECK(gustWind(g, x).u == Approx(g.u.value(x)));
  CHECK(gustWind(g, x).w == Approx(g.w.value(x)));
  CHECK(gustWindRate(g, x, V).u == Approx(g.u.slope(x) * V));
  CHECK(gustWindRate(g, x, V).w == Approx(g.w.slope(x) * V));
}

// ---- Wind injection into the longitudinal EOM -------------------------------

TEST_CASE("Wind-perturbed RHS adds the exact path-axis gust terms", "[wind_gust]") {
  AeroTable table = AeroTable::fromFile(kData + "/example.stab");
  AircraftConfig cfg = loadAircraftConfig(kData + "/aircraft.yaml");
  Mixing mx = Mixing::build(cfg, table);

  LonStateVec X;
  X << 30.0, 18.0, -3.0 * M_PI / 180.0, 2.0 * M_PI / 180.0, 0.01, 6.0, 0.0;
  LonCtrlVec U;
  U << 0.02, 0.0;

  const GustWind W{2.0, -0.8};       // steady part [m/s]
  const GustWind Wdot{0.5, -0.2};    // ramp forcing [m/s^2]
  const LonStateVec base = lonXdotFull(table, mx, cfg, X, U);
  const LonStateVec wind = lonXdotFullWind(table, mx, cfg, X, U, W, Wdot);

  const double sg = std::sin(X[LGAM]), cg = std::cos(X[LGAM]);
  CHECK(wind[LH] == Approx(base[LH] + W.w));
  CHECK(wind[LV] == Approx(base[LV] - (Wdot.u * cg + Wdot.w * sg)));
  CHECK(wind[LGAM] == Approx(base[LGAM] + (Wdot.u * sg - Wdot.w * cg) / X[LV]));
  // Rotational/thrust states see no direct wind forcing (quasi-steady aero on
  // the air-relative state).
  for (int i : {int(LTH), int(LQ), int(LT), int(LTDOT)})
    CHECK(wind[i] == Approx(base[i]));

  // Steady wind (Wdot = 0) only transports the aircraft: V/gamma untouched.
  const LonStateVec steady = lonXdotFullWind(table, mx, cfg, X, U, W, GustWind{});
  CHECK(steady[LV] == Approx(base[LV]));
  CHECK(steady[LGAM] == Approx(base[LGAM]));
  CHECK(steady[LH] == Approx(base[LH] + W.w));

  // Physical sign sanity on the ramp: a tailwind ramp bleeds airspeed, an
  // updraft ramp (Wdot.w > 0) drops gamma_air, i.e. raises alpha = theta-gamma.
  const LonStateVec tail = lonXdotFullWind(table, mx, cfg, X, U, GustWind{}, GustWind{1.0, 0.0});
  CHECK(tail[LV] < base[LV]);
  const LonStateVec updraft = lonXdotFullWind(table, mx, cfg, X, U, GustWind{}, GustWind{0.0, 1.0});
  CHECK(updraft[LGAM] < base[LGAM]);
}
