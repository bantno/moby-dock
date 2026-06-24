#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "autoland/aero_table.hpp"
#include "autoland/config.hpp"
#include "autoland/dynamics.hpp"
#include "autoland/mixing.hpp"
#include "autoland/trim.hpp"

using namespace autoland;
using Catch::Approx;

namespace {
constexpr double kRad2Deg = 180.0 / M_PI;
const std::string kData = AUTOLAND_DATA_DIR;

Dynamics makeDynamics(AeroTable& t, Mixing& mx, AircraftConfig& cfg) {
  t = AeroTable::fromFile(kData + "/example.stab");
  cfg = loadAircraftConfig(kData + "/aircraft.yaml");
  mx = Mixing::build(cfg, t);
  return Dynamics(t, mx, cfg);
}
}  // namespace

TEST_CASE("trim drives the EOM residual to zero", "[trim]") {
  AeroTable t;
  AircraftConfig cfg;
  Mixing mx;
  Dynamics dyn = makeDynamics(t, mx, cfg);

  const double V = 18.0, gamma = -1.5 * M_PI / 180.0;
  TrimResult tr = trim(dyn, V, gamma);

  INFO("residual=" << tr.residual << " alpha_deg=" << tr.alpha * kRad2Deg
                   << " throttle=" << tr.u[DT]);
  CHECK(tr.residual < 1e-5);

  // Re-evaluate the EOM at the reported trim: the three trimmed accelerations
  // must vanish.
  StateVec xd = dyn.xdot(tr.x, tr.u);
  CHECK(std::abs(xd[U]) < 1e-4);
  CHECK(std::abs(xd[W]) < 1e-4);
  CHECK(std::abs(xd[Q]) < 1e-4);
}

TEST_CASE("trim solution is physically sane", "[trim]") {
  AeroTable t;
  AircraftConfig cfg;
  Mixing mx;
  Dynamics dyn = makeDynamics(t, mx, cfg);

  TrimResult tr = trim(dyn, 18.0, -1.5 * M_PI / 180.0);

  // Sane angle of attack for a normal approach.
  CHECK(tr.alpha * kRad2Deg > -5.0);
  CHECK(tr.alpha * kRad2Deg < 12.0);
  // Throttle within physical bounds (descent => modest positive thrust).
  CHECK(tr.u[DT] >= 0.0);
  CHECK(tr.u[DT] <= 1.0);
  // theta = gamma + alpha for wings-level symmetric flight.
  CHECK(tr.theta == Approx(tr.gamma + tr.alpha).epsilon(1e-9));
  // Steady descent: altitude rate is negative, matching V*sin(gamma).
  StateVec xd = dyn.xdot(tr.x, tr.u);
  CHECK(xd[H] == Approx(18.0 * std::sin(-1.5 * M_PI / 180.0)).epsilon(1e-3));
}
