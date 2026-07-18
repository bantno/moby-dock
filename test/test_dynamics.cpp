#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <string>

#include "autoland/config.hpp"
#include "autoland/dynamics.hpp"
#include "autoland/mixing.hpp"

using namespace autoland;
using Catch::Approx;

// Rotational dynamics must satisfy the rigid-body angular-momentum balance
//   hdot_body + omega x h = M_total,   h = [Ixx p - Ixz r,  Iyy q,  Izz r - Ixz p]
// independently of how the coupled p/r system is solved. Run with a NONZERO Ixz
// and all rates nonzero so every product-of-inertia coupling term is exercised:
// aircraft.yaml ships Ixz = 0, which hid a sign error on the Ixz*p*q roll term.
// M_total is reconstructed from the same aeroCoeffs()/thrust buildup xdot() uses,
// so any mismatch is in the Euler-equation assembly itself.
TEST_CASE("6-DOF rotational dynamics satisfy the angular-momentum balance",
          "[dynamics]") {
  const std::string data = AUTOLAND_DATA_DIR;
  AeroTable table = AeroTable::fromFile(data + "/example.stab");
  AircraftConfig cfg = loadAircraftConfig(data + "/aircraft.yaml");
  cfg.inertia.Ixz = 0.35 * std::sqrt(cfg.inertia.Ixx * cfg.inertia.Izz);
  Mixing mx = Mixing::build(cfg, table);
  Dynamics dyn(table, mx, cfg);

  StateVec x = StateVec::Zero();
  x[U] = 17.0; x[V] = 1.2; x[W] = 1.5;    // nonzero alpha and beta
  x[P] = 0.4; x[Q] = -0.3; x[R] = 0.25;   // all body rates nonzero
  x[PHI] = 0.1; x[THETA] = 0.05; x[PSI] = -0.2;
  x[H] = 50.0;
  CtrlVec u = CtrlVec::Zero();
  u[DE] = 0.05; u[DA] = -0.03; u[DR] = 0.02; u[DT] = 0.4;

  const StateVec xd = dyn.xdot(x, u);

  // Reconstruct the applied moments from the same coefficient buildup.
  const double Vt = std::sqrt(x[U] * x[U] + x[V] * x[V] + x[W] * x[W]);
  const double qbar = 0.5 * cfg.env.rho * Vt * Vt;
  const CoefVec C = dyn.aeroCoeffs(x, u);
  const double S = table.Sref(), b = table.Bref(), c = table.Cref();
  const double T =
      std::max(0.0, u[DT] * (cfg.thrust.T_static - cfg.thrust.k_v * Vt));
  const double Lm = qbar * S * b * C[CMX];
  const double Mm = qbar * S * c * C[CMY] + T * cfg.thrust.zcp;
  const double Nm = qbar * S * b * C[CMZ];

  const InertiaParams& I = cfg.inertia;
  const double p = x[P], q = x[Q], r = x[R];
  const double hx = I.Ixx * p - I.Ixz * r;
  const double hy = I.Iyy * q;
  const double hz = I.Izz * r - I.Ixz * p;
  const double hxd = I.Ixx * xd[P] - I.Ixz * xd[R];
  const double hyd = I.Iyy * xd[Q];
  const double hzd = I.Izz * xd[R] - I.Ixz * xd[P];

  CHECK(hxd + (q * hz - r * hy) == Approx(Lm).epsilon(1e-9).margin(1e-9));
  CHECK(hyd + (r * hx - p * hz) == Approx(Mm).epsilon(1e-9).margin(1e-9));
  CHECK(hzd + (p * hy - q * hx) == Approx(Nm).epsilon(1e-9).margin(1e-9));
}
