#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "autoland/aero_table.hpp"
#include "autoland/config.hpp"
#include "autoland/dynamics.hpp"
#include "autoland/linear_model.hpp"
#include "autoland/mixing.hpp"
#include "autoland/trim.hpp"

using namespace autoland;
using Catch::Approx;

namespace {
const std::string kData = AUTOLAND_DATA_DIR;
}

TEST_CASE("linear model has correct shape and finite entries", "[linear]") {
  AeroTable t = AeroTable::fromFile(kData + "/example.stab");
  AircraftConfig cfg = loadAircraftConfig(kData + "/aircraft.yaml");
  Mixing mx = Mixing::build(cfg, t);
  Dynamics dyn(t, mx, cfg);
  TrimResult tr = trim(dyn, 18.0, -1.5 * M_PI / 180.0);
  LinearModel m = linearize(dyn, tr.x, tr.u);

  CHECK(m.A.rows() == NX);
  CHECK(m.A.cols() == NX);
  CHECK(m.B.rows() == NX);
  CHECK(m.B.cols() == NU);
  CHECK(m.A_lon.rows() == 5);
  CHECK(m.A_lon.cols() == 5);
  CHECK(m.B_lon.cols() == 2);
  CHECK(m.A_lat.rows() == 6);
  CHECK(m.B_lat.cols() == 2);
  CHECK(m.A.allFinite());
  CHECK(m.B.allFinite());
  CHECK(m.f0.allFinite());
}

TEST_CASE("linear model reproduces the nonlinear EOM to first order",
          "[linear]") {
  AeroTable t = AeroTable::fromFile(kData + "/example.stab");
  AircraftConfig cfg = loadAircraftConfig(kData + "/aircraft.yaml");
  Mixing mx = Mixing::build(cfg, t);
  Dynamics dyn(t, mx, cfg);
  TrimResult tr = trim(dyn, 18.0, -1.5 * M_PI / 180.0);
  LinearModel m = linearize(dyn, tr.x, tr.u);

  // At the trim point the affine model must equal the nonlinear drift exactly.
  StateVec lin0 = m.xdot(tr.x, tr.u);
  StateVec non0 = dyn.xdot(tr.x, tr.u);
  CHECK((lin0 - non0).norm() < 1e-9);

  // For a small state perturbation, the linear prediction matches the nonlinear
  // EOM closely (first-order agreement).
  StateVec dx = StateVec::Zero();
  dx[Q] = 0.01;     // 0.01 rad/s pitch rate
  dx[THETA] = 0.01; // 0.01 rad pitch
  StateVec xp = tr.x + dx;
  StateVec lin = m.xdot(xp, tr.u);
  StateVec non = dyn.xdot(xp, tr.u);
  CHECK((lin - non).norm() < 1e-3);
}

TEST_CASE("longitudinal sub-model shows static + control sign sense",
          "[linear]") {
  AeroTable t = AeroTable::fromFile(kData + "/example.stab");
  AircraftConfig cfg = loadAircraftConfig(kData + "/aircraft.yaml");
  Mixing mx = Mixing::build(cfg, t);
  Dynamics dyn(t, mx, cfg);
  TrimResult tr = trim(dyn, 18.0, -1.5 * M_PI / 180.0);
  LinearModel m = linearize(dyn, tr.x, tr.u);

  // Longitudinal state order [u, w, q, theta, h]; control [delta_e, delta_T].
  // d(qdot)/dw < 0  : positive static pitch stiffness (Cm_alpha < 0).
  CHECK(m.A_lon(2, 1) < 0.0);
  // d(qdot)/d(delta_e) > 0 : virtual elevator polarity is +delta_e -> nose-up.
  CHECK(m.B_lon(2, 0) > 0.0);
  // d(hdot)/d(theta) > 0 : pitching up increases altitude rate.
  CHECK(m.A_lon(4, 3) > 0.0);
}
