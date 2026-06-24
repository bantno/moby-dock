#include "autoland/linear_model.hpp"

#include <cmath>

namespace autoland {
namespace {

// Per-state central-difference steps, scaled to each variable's natural size.
StateVec stateSteps() {
  StateVec s;
  s[U] = s[V] = s[W] = 1e-5;        // m/s
  s[P] = s[Q] = s[R] = 1e-6;        // rad/s
  s[PHI] = s[THETA] = s[PSI] = 1e-6;  // rad
  s[H] = s[Y] = 1e-4;               // m
  return s;
}

CtrlVec ctrlSteps() {
  CtrlVec s;
  s[DE] = s[DA] = s[DR] = 1e-6;  // rad
  s[DT] = 1e-5;                  // throttle
  return s;
}

}  // namespace

LinearModel linearize(const Dynamics& dyn, const StateVec& x0,
                      const CtrlVec& u0) {
  LinearModel m;
  m.x0 = x0;
  m.u0 = u0;
  m.f0 = dyn.xdot(x0, u0);

  m.A = Mat::Zero(NX, NX);
  m.B = Mat::Zero(NX, NU);

  const StateVec sx = stateSteps();
  for (int j = 0; j < NX; ++j) {
    StateVec xp = x0, xm = x0;
    xp[j] += sx[j];
    xm[j] -= sx[j];
    const StateVec fp = dyn.xdot(xp, u0);
    const StateVec fm = dyn.xdot(xm, u0);
    m.A.col(j) = (fp - fm) / (2.0 * sx[j]);
  }

  const CtrlVec su = ctrlSteps();
  for (int j = 0; j < NU; ++j) {
    CtrlVec up = u0, um = u0;
    up[j] += su[j];
    um[j] -= su[j];
    const StateVec fp = dyn.xdot(x0, up);
    const StateVec fm = dyn.xdot(x0, um);
    m.B.col(j) = (fp - fm) / (2.0 * su[j]);
  }

  // Extract the decoupled sub-models as index selections of the full A,B.
  auto sub = [](const Mat& M, const int* rows, int nr, const int* cols, int nc) {
    Mat out = Mat::Zero(nr, nc);
    for (int i = 0; i < nr; ++i)
      for (int j = 0; j < nc; ++j) out(i, j) = M(rows[i], cols[j]);
    return out;
  };

  m.A_lon = sub(m.A, kLonStates.data(), 5, kLonStates.data(), 5);
  m.B_lon = sub(m.B, kLonStates.data(), 5, kLonCtrls.data(), 2);
  m.A_lat = sub(m.A, kLatStates.data(), 6, kLatStates.data(), 6);
  m.B_lat = sub(m.B, kLatStates.data(), 6, kLatCtrls.data(), 2);

  return m;
}

}  // namespace autoland
