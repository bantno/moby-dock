#include "autoland/trim.hpp"

#include <cmath>

namespace autoland {
namespace {

// Build the steady-flight state and controls from the 3 trim unknowns
// v = [alpha, delta_e, throttle], for wings-level zero-sideslip flight at the
// given V and gamma.  theta = gamma + alpha.
void assemble(double V, double gamma, const Eigen::Vector3d& v, StateVec& x,
              CtrlVec& u) {
  const double alpha = v[0];
  x = StateVec::Zero();
  x[U] = V * std::cos(alpha);
  x[W] = V * std::sin(alpha);
  // v_body = 0 (zero sideslip), rates = 0, phi = psi = 0
  x[THETA] = gamma + alpha;

  u = CtrlVec::Zero();
  u[DE] = v[1];
  u[DT] = v[2];
}

// Residual: the three accelerations that must vanish at trim.
Eigen::Vector3d residual(const Dynamics& dyn, double V, double gamma,
                         const Eigen::Vector3d& v) {
  StateVec x;
  CtrlVec u;
  assemble(V, gamma, v, x, u);
  const StateVec xd = dyn.xdot(x, u);
  return Eigen::Vector3d(xd[U], xd[W], xd[Q]);
}

}  // namespace

TrimResult trim(const Dynamics& dyn_in, double V_app, double gamma_app,
                double tol, int max_iter) {
  TrimResult res;
  res.V = V_app;
  res.gamma = gamma_app;

  // Solve against a SMOOTH (unclamped) thrust model so the throttle column of
  // the Jacobian never vanishes (the plant's max(0,.) clamp would make Newton
  // singular wherever throttle goes negative). The returned throttle is then
  // checked against [0,1] by the caller/tests.
  Dynamics dyn = dyn_in;
  const ThrustParams& tp = dyn_in.config().thrust;
  const double Ts = tp.T_static, kv = tp.k_v;
  dyn.setThrustModel(
      [Ts, kv](double throttle, double V) { return throttle * (Ts - kv * V); });

  // Initial guess: small alpha, neutral elevator, mid throttle.
  Eigen::Vector3d v(0.03, 0.0, 0.3);

  const double h = 1e-6;  // central-difference step for the Jacobian
  Eigen::Vector3d R = residual(dyn, V_app, gamma_app, v);

  int it = 0;
  for (; it < max_iter; ++it) {
    if (R.norm() < tol) break;

    // Numerical Jacobian (3x3) via central differences.
    Eigen::Matrix3d J;
    for (int j = 0; j < 3; ++j) {
      Eigen::Vector3d vp = v, vm = v;
      vp[j] += h;
      vm[j] -= h;
      J.col(j) = (residual(dyn, V_app, gamma_app, vp) -
                  residual(dyn, V_app, gamma_app, vm)) /
                 (2.0 * h);
    }

    // Newton step with a simple damping/backtracking guard.
    Eigen::Vector3d dv = J.fullPivLu().solve(-R);
    double step = 1.0;
    Eigen::Vector3d v_new = v + step * dv;
    Eigen::Vector3d R_new = residual(dyn, V_app, gamma_app, v_new);
    for (int bt = 0; bt < 10 && R_new.norm() > R.norm(); ++bt) {
      step *= 0.5;
      v_new = v + step * dv;
      R_new = residual(dyn, V_app, gamma_app, v_new);
    }
    v = v_new;
    R = R_new;
  }

  assemble(V_app, gamma_app, v, res.x, res.u);
  res.alpha = v[0];
  res.theta = gamma_app + v[0];
  res.residual = R.norm();
  res.iterations = it;
  res.converged = (R.norm() < tol);
  return res;
}

}  // namespace autoland
