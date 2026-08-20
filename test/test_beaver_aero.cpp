#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>

#include "autoland/beaver_aero.hpp"

using namespace autoland;
using Catch::Approx;

namespace {
enum { CX, CY, CZ, CL, CM, CN };  // beaverAeroCoeffs() index order
}  // namespace

// Zero-state: every polynomial term vanishes, so the coefficients reduce to the
// constant (bias) terms of the verified LR-556 Table 3 model.
TEST_CASE("Beaver aero zero-state returns the constant coefficients", "[beaver]") {
  const auto C = beaverAeroCoeffs(0, 0, 0, 0, 0, 50.0, 0, 0, 0, 0, 0);
  CHECK(C[CX] == Approx(-0.03554).margin(1e-9));
  CHECK(C[CY] == Approx(-0.002226).margin(1e-9));
  CHECK(C[CZ] == Approx(-0.05504).margin(1e-9));
  CHECK(C[CL] == Approx(0.000591).margin(1e-9));
  CHECK(C[CM] == Approx(0.09448).margin(1e-9));
  CHECK(C[CN] == Approx(-0.003117).margin(1e-9));
}

// Full-chain oracle: reproduce the model at the LR-556/FDC reference flight
// condition (GitHub Flight_Simulator FLightScript.m IC), exercising the engine
// power -> dpt map and every aero polynomial. Oracle computed from an independent
// Python port of the verified model.
TEST_CASE("Beaver aero + engine reproduce the reference condition", "[beaver]") {
  const double V = 35.0, alpha = 0.21131, beta = -0.020667;
  const double de = -0.093083, da = 0.0096242, dr = -0.049242, df = 0.0;
  const double z = 2000.0 * 0.3048;  // 2000 ft
  const double rho =
      1.225 * std::exp(-9.81 * z / 287.05 / (288.0 - 0.0065 * z));
  CHECK(rho == Approx(1.13837).epsilon(1e-4));

  const double P = beaverEnginePower(20.0, 1800.0, rho);
  const double dpt = beaverDpt(P, rho, V);
  CHECK(P == Approx(88.3875).epsilon(1e-4));
  CHECK(dpt == Approx(0.77939).epsilon(1e-4));

  const auto C =
      beaverAeroCoeffs(alpha, beta, 0, 0, 0, V, de, da, dr, df, dpt);
  CHECK(C[CX] == Approx(0.26758).margin(1e-4));
  CHECK(C[CY] == Approx(0.00221).margin(1e-4));
  CHECK(C[CZ] == Approx(-1.28539).margin(1e-4));
  CHECK(C[CL] == Approx(-0.00009).margin(1e-4));
  CHECK(C[CM] == Approx(-0.01088).margin(1e-4));
  CHECK(C[CN] == Approx(-0.00066).margin(1e-4));
}

// Physical sign conventions (body axes, standard signs) -- the same six
// constraints used to pin the VSPAERO frame in linear_model.hpp.
TEST_CASE("Beaver aero has the correct static-stability signs", "[beaver]") {
  const double V = 40.0, h = 1e-4;
  auto Cm = [&](double a) {
    return beaverAeroCoeffs(a, 0, 0, 0, 0, V, 0, 0, 0, 0, 0)[CM];
  };
  auto Cz = [&](double a) {
    return beaverAeroCoeffs(a, 0, 0, 0, 0, V, 0, 0, 0, 0, 0)[CZ];
  };
  auto Cn = [&](double be) {
    return beaverAeroCoeffs(0, be, 0, 0, 0, V, 0, 0, 0, 0, 0)[CN];
  };
  auto Cl = [&](double be) {
    return beaverAeroCoeffs(0, be, 0, 0, 0, V, 0, 0, 0, 0, 0)[CL];
  };
  CHECK((Cm(0.1 + h) - Cm(0.1 - h)) / (2 * h) < 0.0);  // Cm_alpha < 0 (pitch stable)
  CHECK((Cz(0.1 + h) - Cz(0.1 - h)) / (2 * h) < 0.0);  // dCZ/da < 0 (lift up with a)
  CHECK((Cn(0.05 + h) - Cn(0.05 - h)) / (2 * h) > 0.0);  // Cn_beta > 0 (weathercock)
  CHECK((Cl(0.05 + h) - Cl(0.05 - h)) / (2 * h) < 0.0);  // Cl_beta < 0 (dihedral)

  // Elevator authority is present with the right sign (pitch down for +de).
  const double V2 = 40.0;
  const double Cm_de0 = beaverAeroCoeffs(0, 0, 0, 0, 0, V2, 0, 0, 0, 0, 0)[CM];
  const double Cm_de1 = beaverAeroCoeffs(0, 0, 0, 0, 0, V2, 0.1, 0, 0, 0, 0)[CM];
  CHECK((Cm_de1 - Cm_de0) / 0.1 < 0.0);  // Cm_de < 0
}

// Verified reference geometry / mass / inertia (LR-556 Tables 1-2).
TEST_CASE("Beaver reference geometry and inertia match LR-556", "[beaver]") {
  CHECK(BeaverGeom::S == Approx(23.23));
  CHECK(BeaverGeom::b == Approx(14.63));
  CHECK(BeaverGeom::c == Approx(1.5875));
  CHECK(BeaverGeom::mass == Approx(2288.231));
  CHECK(BeaverGeom::Ix == Approx(5368.39));
  CHECK(BeaverGeom::Iy == Approx(6928.93));
  CHECK(BeaverGeom::Iz == Approx(11158.75));
  CHECK(BeaverGeom::Ixz == Approx(117.64));
}
