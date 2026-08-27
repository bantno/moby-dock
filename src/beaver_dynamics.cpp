#include "autoland/beaver_dynamics.hpp"

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

#include <cmath>

// Exact derivatives for the Beaver plant: the state-space linearization and
// the trim Jacobian are computed by seeding autodiff::dual through the SAME
// templated EOM the plant integrates (xdotT) -- no finite differences.
namespace autoland {
namespace {
using autodiff::dual;

// Steady straight wings-level flight assembled from the 6 trim unknowns
// z = [alpha, beta, de, da, dr, dT] at airspeed Vtrim / path angle gamma:
// the six body-acceleration residuals [udot vdot wdot pdot qdot rdot].
template <class T>
Eigen::Matrix<T, 6, 1> trimResidualT(const BeaverDynamics& dyn, double Vtrim,
                                     double gamma,
                                     const Eigen::Matrix<T, 6, 1>& z) {
  using std::asin;
  using std::cos;
  using std::sin;
  const T alpha = z[0], beta = z[1];
  // Wings-level (phi = 0) exact relation: sin(gamma) = cos(beta) sin(theta -
  // alpha)  =>  theta = alpha + asin(sin(gamma)/cos(beta)).
  const T theta = alpha + asin(T(std::sin(gamma)) / cos(beta));

  Eigen::Matrix<T, NX, 1> x;
  x.setZero();
  x[U] = Vtrim * cos(alpha) * cos(beta);
  x[V] = Vtrim * sin(beta);
  x[W] = Vtrim * sin(alpha) * cos(beta);
  x[THETA] = theta;

  Eigen::Matrix<T, NU, 1> u;
  u.setZero();
  u[DE] = z[2];
  u[DA] = z[3];
  u[DR] = z[4];
  u[DT] = z[5];

  const Eigen::Matrix<T, NX, 1> xd =
      dyn.xdotT<T>(x, x[U], x[V], x[W], u);
  // NOTE: the residual local must NOT be named 'R' -- that shadows the
  // State::R enum index in xd[R] (cf. the 'Vt' note in dynamics.cpp).
  Eigen::Matrix<T, 6, 1> res;
  res << xd[U], xd[V], xd[W], xd[P], xd[Q], xd[R];
  return res;
}
}  // namespace

void BeaverDynamics::linearize(const StateVec& x0, const CtrlVec& u0, Mat& A,
                               Mat& B) const {
  A.resize(NX, NX);
  B.resize(NX, NU);
  for (int j = 0; j < NX; ++j) {
    Eigen::Matrix<dual, NX, 1> x = x0.cast<dual>();
    Eigen::Matrix<dual, NU, 1> u = u0.cast<dual>();
    x[j].grad = 1.0;
    const Eigen::Matrix<dual, NX, 1> xd = xdotT<dual>(x, x[U], x[V], x[W], u);
    for (int i = 0; i < NX; ++i) A(i, j) = xd[i].grad;
  }
  for (int j = 0; j < NU; ++j) {
    Eigen::Matrix<dual, NX, 1> x = x0.cast<dual>();
    Eigen::Matrix<dual, NU, 1> u = u0.cast<dual>();
    u[j].grad = 1.0;
    const Eigen::Matrix<dual, NX, 1> xd = xdotT<dual>(x, x[U], x[V], x[W], u);
    for (int i = 0; i < NX; ++i) B(i, j) = xd[i].grad;
  }
}

TrimResult beaverTrim(const BeaverDynamics& dyn, double V_target, double gamma,
                      double tol, int max_iter) {
  TrimResult res;
  res.V = V_target;
  res.gamma = gamma;

  using Vec6 = Eigen::Matrix<double, 6, 1>;
  Vec6 z;
  z << 0.05, 0.0, -0.05, 0.0, 0.0, 0.5;  // [alpha beta de da dr dT]

  Vec6 R = trimResidualT<double>(dyn, V_target, gamma, z);
  int it = 0;
  for (; it < max_iter; ++it) {
    if (R.norm() < tol) break;

    Eigen::Matrix<double, 6, 6> J;
    for (int j = 0; j < 6; ++j) {
      Eigen::Matrix<dual, 6, 1> zd = z.cast<dual>();
      zd[j].grad = 1.0;
      const Eigen::Matrix<dual, 6, 1> Rd =
          trimResidualT<dual>(dyn, V_target, gamma, zd);
      for (int i = 0; i < 6; ++i) J(i, j) = Rd[i].grad;
    }

    const Vec6 dz = J.fullPivLu().solve(-R);
    double step = 1.0;
    Vec6 z_new = z + step * dz;
    Vec6 R_new = trimResidualT<double>(dyn, V_target, gamma, z_new);
    for (int bt = 0; bt < 12 && R_new.norm() > R.norm(); ++bt) {
      step *= 0.5;
      z_new = z + step * dz;
      R_new = trimResidualT<double>(dyn, V_target, gamma, z_new);
    }
    z = z_new;
    R = R_new;
  }

  const double alpha = z[0], beta = z[1];
  res.alpha = alpha;
  res.theta = alpha + std::asin(std::sin(gamma) / std::cos(beta));
  res.x = StateVec::Zero();
  res.x[U] = V_target * std::cos(alpha) * std::cos(beta);
  res.x[V] = V_target * std::sin(beta);
  res.x[W] = V_target * std::sin(alpha) * std::cos(beta);
  res.x[THETA] = res.theta;
  res.u = CtrlVec::Zero();
  res.u[DE] = z[2];
  res.u[DA] = z[3];
  res.u[DR] = z[4];
  res.u[DT] = z[5];
  res.residual = R.norm();
  res.iterations = it;
  res.converged = (R.norm() < tol);
  return res;
}

}  // namespace autoland
