#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <memory>

#include "autoland/cbf.hpp"
#include "autoland/config.hpp"
#include "autoland/control_affine_model.hpp"
#include "autoland/qp_solver.hpp"

using namespace autoland;
using Catch::Approx;

namespace {
// A constant control-affine model xdot = f + B u, chosen so the CBF QP rows can
// be reasoned about by hand. Doubles as a check that the filter depends only on
// the ControlAffineModel interface (the dynamics-model swap seam).
struct ToyModel : ControlAffineModel {
  StateVec f{StateVec::Zero()};
  Mat B{Mat::Zero(NX, NU)};
  StateVec drift(const StateVec&) const override { return f; }
  Mat ctrlMatrix(const StateVec&) const override { return B; }
};

// Airspeed barrier b = V - V_min with the analytic gradient.
Barrier airspeedBarrier(double V_min, double alpha_gain) {
  Barrier b;
  b.name = "min_airspeed";
  // NB: name the local 'Vt', never 'V' -- 'V' is the State::V enum index, so a
  // local 'V' would make x[V] a self-referential, uninitialized index read.
  b.h = [V_min](const StateVec& x) {
    const double Vt = std::sqrt(x[U] * x[U] + x[V] * x[V] + x[W] * x[W]);
    return Vt - V_min;
  };
  b.grad = [](const StateVec& x) {
    Eigen::RowVectorXd g = Eigen::RowVectorXd::Zero(NX);
    const double Vt = std::max(1e-9, std::sqrt(x[U] * x[U] + x[V] * x[V] +
                                               x[W] * x[W]));
    g[U] = x[U] / Vt; g[V] = x[V] / Vt; g[W] = x[W] / Vt;
    return g;
  };
  b.alpha = [alpha_gain](double hv) { return alpha_gain * hv; };
  return b;
}

CBFConfig lonConfig() {
  CBFConfig c;
  c.ctrl_idx = {DE, DT};
  c.ctrl_weights = {1.0, 1.0};
  c.slack_penalty = 1.0e4;
  return c;
}
}  // namespace

TEST_CASE("OsqpSolver solves a trivial box-constrained QP", "[cbf][qp]") {
  // min (z-3)^2  s.t. 0 <= z <= 1   ->  z* = 1.
  OsqpSolver solver;
  Mat P(1, 1); P(0, 0) = 2.0;
  Vec q(1);    q[0] = -6.0;
  Mat A(1, 1); A(0, 0) = 1.0;
  Vec l(1);    l[0] = 0.0;
  Vec u(1);    u[0] = 1.0;
  QPResult r = solver.solve(P, q, A, l, u);
  CHECK(r.success);
  CHECK(r.z[0] == Approx(1.0).margin(1e-4));
}

TEST_CASE("CBF passes the nominal through when safely inside the set",
          "[cbf]") {
  ToyModel m;
  m.f[U] = 2.0;            // airspeed increasing -> barrier derivative positive
  m.B(U, DT) = 10.0;
  CBFFilter filt(lonConfig(), SurfaceLimits{});
  std::vector<Barrier> bars{airspeedBarrier(15.0, 1.0)};

  StateVec x = StateVec::Zero();
  x[U] = 16.0;             // V = 16 > V_min = 15, comfortable margin
  CtrlVec u_nom = CtrlVec::Zero();

  CtrlVec u = filt.filter(u_nom, x, m, bars);
  CHECK(u[DE] == Approx(0.0).margin(1e-4));
  CHECK(u[DT] == Approx(0.0).margin(1e-4));
}

TEST_CASE("CBF raises throttle to keep airspeed above the margin", "[cbf]") {
  // Lf h = f[U] = -5, Lg h via throttle = B(U,DT) = 10, alpha(h) = 1.
  // Constraint: -5 + 10*dT + 1 >= 0  ->  dT >= 0.4. Nominal dT = 0 violates.
  ToyModel m;
  m.f[U] = -5.0;
  m.B(U, DT) = 10.0;
  CBFFilter filt(lonConfig(), SurfaceLimits{});
  std::vector<Barrier> bars{airspeedBarrier(15.0, 1.0)};

  StateVec x = StateVec::Zero();
  x[U] = 16.0;             // V = 16, h = 1
  CtrlVec u_nom = CtrlVec::Zero();

  CtrlVec u = filt.filter(u_nom, x, m, bars);
  CHECK(u[DT] == Approx(0.4).margin(1e-2));   // minimal correction
  CHECK(u[DE] == Approx(0.0).margin(1e-2));   // elevator unchanged
}

TEST_CASE("CBF stays feasible when a barrier has no control authority",
          "[cbf]") {
  // A high-relative-degree-like barrier: grad points along THETA, but B has no
  // THETA row, so Lg h = 0 (control cannot affect it). With Lf h + alpha < 0 the
  // row is unsatisfiable in u; the slack must absorb it and the QP stay solvable.
  ToyModel m;
  m.f[U] = -5.0;
  m.B(U, DT) = 10.0;
  m.f[THETA] = -10.0;      // uncontrollable, "violated" barrier derivative

  Barrier uncontrollable;
  uncontrollable.name = "no_authority";
  uncontrollable.h = [](const StateVec& x) { return x[THETA]; };
  uncontrollable.grad = [](const StateVec&) {
    Eigen::RowVectorXd g = Eigen::RowVectorXd::Zero(NX);
    g[THETA] = 1.0;
    return g;
  };
  uncontrollable.alpha = [](double hv) { return 1.0 * hv; };

  CBFFilter filt(lonConfig(), SurfaceLimits{});
  std::vector<Barrier> bars{airspeedBarrier(15.0, 1.0), uncontrollable};

  StateVec x = StateVec::Zero();
  x[U] = 16.0;
  CtrlVec u = filt.filter(CtrlVec::Zero(), x, m, bars);

  // Still returns a finite, in-bounds control, and the controllable (airspeed)
  // barrier is still enforced.
  CHECK(std::isfinite(u[DE]));
  CHECK(std::isfinite(u[DT]));
  CHECK(u[DT] >= 0.0);
  CHECK(u[DT] <= 1.0);
  CHECK(u[DT] == Approx(0.4).margin(1e-2));
}

TEST_CASE("Hard barrier with authority is enforced exactly (no slack)",
          "[cbf]") {
  // Same setup as the airspeed-activation case but the barrier is HARD: with
  // control authority (B(U,DT)=10) the constraint -5 + 10 dT + 1 >= 0 must hold
  // exactly -> dT = 0.4, and no feasibility recovery is needed.
  ToyModel m;
  m.f[U] = -5.0;
  m.B(U, DT) = 10.0;
  CBFFilter filt(lonConfig(), SurfaceLimits{});
  Barrier b = airspeedBarrier(15.0, 1.0);
  b.hard = true;
  std::vector<Barrier> bars{b};

  StateVec x = StateVec::Zero();
  x[U] = 16.0;
  CtrlVec u = filt.filter(CtrlVec::Zero(), x, m, bars);
  CHECK(u[DT] == Approx(0.4).margin(1e-2));
  CHECK_FALSE(filt.lastFeasibilityRecovery());
}

TEST_CASE("Hard barrier with no authority triggers feasibility recovery",
          "[cbf]") {
  // A hard barrier whose gradient lies along THETA, with no THETA row in B, so
  // Lg h = 0 and the controls cannot satisfy it while it is violated. The hard
  // solve is infeasible; the filter must recover (soften) and still return a
  // valid, in-box control, with the soft airspeed barrier still respected.
  ToyModel m;
  m.f[U] = -5.0;
  m.B(U, DT) = 10.0;
  m.f[THETA] = -10.0;

  Barrier hard;
  hard.name = "no_authority";
  hard.h = [](const StateVec& x) { return x[THETA]; };
  hard.grad = [](const StateVec&) {
    Eigen::RowVectorXd g = Eigen::RowVectorXd::Zero(NX);
    g[THETA] = 1.0;
    return g;
  };
  hard.hard = true;

  CBFFilter filt(lonConfig(), SurfaceLimits{});
  std::vector<Barrier> bars{airspeedBarrier(15.0, 1.0), hard};

  StateVec x = StateVec::Zero();
  x[U] = 16.0;
  CtrlVec u = filt.filter(CtrlVec::Zero(), x, m, bars);
  CHECK(filt.lastFeasibilityRecovery());     // hard set was infeasible -> softened
  CHECK(std::isfinite(u[DE]));
  CHECK(u[DT] >= 0.0);
  CHECK(u[DT] <= 1.0);
  CHECK(u[DT] == Approx(0.4).margin(1e-2));   // soft airspeed still defended
}

TEST_CASE("Disabled CBF is a hard pass-through", "[cbf]") {
  ToyModel m;
  m.f[U] = -5.0;
  m.B(U, DT) = 10.0;
  CBFFilter filt(lonConfig(), SurfaceLimits{});
  filt.setEnabled(false);
  std::vector<Barrier> bars{airspeedBarrier(15.0, 1.0)};

  StateVec x = StateVec::Zero();
  x[U] = 16.0;
  CtrlVec u_nom = CtrlVec::Zero();
  u_nom[DT] = 0.1;
  CtrlVec u = filt.filter(u_nom, x, m, bars);
  CHECK(u[DT] == Approx(0.1));
  CHECK(u[DE] == Approx(0.0));
}

