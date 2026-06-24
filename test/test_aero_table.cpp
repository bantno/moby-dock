#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>

#include "autoland/aero_table.hpp"

using namespace autoland;
using Catch::Approx;

namespace {
constexpr double kDeg2Rad = M_PI / 180.0;
const std::string kStab = std::string(AUTOLAND_DATA_DIR) + "/example.stab";

bool contains(const std::vector<std::string>& v, const std::string& s) {
  return std::find(v.begin(), v.end(), s) != v.end();
}
}  // namespace

TEST_CASE("stab reference quantities are parsed", "[aero_table]") {
  AeroTable t = AeroTable::fromFile(kStab);
  CHECK(t.Sref() == Approx(0.425).epsilon(1e-6));
  CHECK(t.Cref() == Approx(0.25).epsilon(1e-6));
  CHECK(t.Bref() == Approx(1.7).epsilon(1e-6));
  CHECK(t.Xcg() == Approx(0.445).epsilon(1e-6));
}

TEST_CASE("control group names are read from the file header", "[aero_table]") {
  AeroTable t = AeroTable::fromFile(kStab);
  REQUIRE(t.numControlGroups() == 3);
  // Order matters: index 0=Ailerons, 1=Elevator, 2=Rudder (verified vs the
  // perturbation rows). The mixing/dynamics rely on this ordering.
  const auto& n = t.controlGroupNames();
  CHECK(n[0] == "Ailerons");
  CHECK(n[1] == "Elevator");
  CHECK(n[2] == "Rudder");
  CHECK(contains(n, "Elevator"));
}

TEST_CASE("sweep grid is a complete 5x5 (alpha,beta) at single Mach",
          "[aero_table]") {
  AeroTable t = AeroTable::fromFile(kStab);
  CHECK(t.alphaGrid().size() == 5);
  CHECK(t.betaGrid().size() == 5);
  CHECK(t.machGrid().size() == 1);
  CHECK(t.machGrid()[0] == Approx(0.059));
}

TEST_CASE("node lookup returns parsed derivative values", "[aero_table]") {
  // Values taken directly from the (AoA=-10, Beta=-10, Mach=0.059) block of
  // data/example.stab. **USER: these are the derivative values to confirm.**
  AeroTable t = AeroTable::fromFile(kStab);
  AeroLookup L = t.lookup(-10.0 * kDeg2Rad, -10.0 * kDeg2Rad, 0.059);
  CHECK_FALSE(L.clamped);

  // Base coefficients.
  CHECK(L.d.base[CFX] == Approx(-0.0932965).epsilon(1e-5));
  CHECK(L.d.base[CFZ] == Approx(-0.5759537).epsilon(1e-5));
  CHECK(L.d.base[CMY] == Approx(0.0888597).epsilon(1e-5));

  // Stability derivatives (per radian).
  CHECK(L.d.d_alpha[CFZ] == Approx(8.6899099).epsilon(1e-5));
  CHECK(L.d.d_alpha[CMY] == Approx(-5.4178023).epsilon(1e-5));  // Cm_alpha < 0
  CHECK(L.d.d_beta[CMZ] == Approx(-0.7668386).epsilon(1e-5));

  // Rate derivatives (wrt nondimensional rates).
  CHECK(L.d.d_q[CMY] == Approx(-4.9450967).epsilon(1e-5));  // Cm_q < 0 (damping)

  // Control-group derivatives (per radian): 0=Ailerons,1=Elevator,2=Rudder.
  CHECK(L.d.d_ctrl[0][CMX] == Approx(0.3977618).epsilon(1e-5));   // Cl_aileron
  CHECK(L.d.d_ctrl[1][CMY] == Approx(-1.0551404).epsilon(1e-5));  // Cm_elevator
  CHECK(L.d.d_ctrl[2][CMZ] == Approx(-0.1703910).epsilon(1e-5));  // Cn_rudder
}

TEST_CASE("bilinear interpolation hits the analytic midpoint", "[aero_table]") {
  AeroTable t = AeroTable::fromFile(kStab);
  // Midway between (AoA=-10,Beta=-10) CFz=-0.5759537 and
  //                (AoA= -5,Beta=-10) CFz=-0.1192710.
  AeroLookup L = t.lookup(-7.5 * kDeg2Rad, -10.0 * kDeg2Rad, 0.059);
  const double expected = 0.5 * (-0.5759537 + -0.1192710);
  CHECK(L.d.base[CFZ] == Approx(expected).epsilon(1e-5));
  CHECK_FALSE(L.clamped);
}

TEST_CASE("off-grid queries clamp and flag, never extrapolate", "[aero_table]") {
  AeroTable t = AeroTable::fromFile(kStab);
  AeroLookup L = t.lookup(-20.0 * kDeg2Rad, -10.0 * kDeg2Rad, 0.059);
  CHECK(L.clamped);
  // Clamped to the alpha=-10 boundary node.
  CHECK(L.d.base[CFZ] == Approx(-0.5759537).epsilon(1e-5));
}
