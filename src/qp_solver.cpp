#include "autoland/qp_solver.hpp"

#include <osqp.h>

#include <vector>

namespace autoland {
namespace {

// Convert a dense Eigen matrix to OSQP compressed-sparse-column arrays. When
// upper_only is set (for the cost matrix P) only the upper triangle is emitted,
// as OSQP requires. Structural zeros are dropped.
void denseToCsc(const Mat& M, bool upper_only, std::vector<c_float>& x,
                std::vector<c_int>& i, std::vector<c_int>& p) {
  const int rows = static_cast<int>(M.rows());
  const int cols = static_cast<int>(M.cols());
  x.clear();
  i.clear();
  p.clear();
  p.reserve(cols + 1);
  p.push_back(0);
  for (int c = 0; c < cols; ++c) {
    for (int r = 0; r < rows; ++r) {
      if (upper_only && r > c) continue;
      const double v = M(r, c);
      if (v != 0.0) {
        x.push_back(static_cast<c_float>(v));
        i.push_back(static_cast<c_int>(r));
      }
    }
    p.push_back(static_cast<c_int>(x.size()));
  }
}

}  // namespace

QPResult OsqpSolver::solve(const Mat& P, const Vec& q, const Mat& A,
                           const Vec& l, const Vec& u) const {
  const c_int n = static_cast<c_int>(P.rows());
  const c_int m = static_cast<c_int>(A.rows());

  QPResult r;
  r.z = Vec::Zero(n);
  if (n == 0) {
    r.success = true;
    return r;
  }

  // CSC arrays must outlive the workspace; OSQP copies them during setup but we
  // keep them alive through cleanup to be safe.
  std::vector<c_float> Px, Ax;
  std::vector<c_int> Pi, Pp, Ai, Ap;
  denseToCsc(P, /*upper_only=*/true, Px, Pi, Pp);
  denseToCsc(A, /*upper_only=*/false, Ax, Ai, Ap);

  std::vector<c_float> q_(q.data(), q.data() + n);
  std::vector<c_float> l_(l.data(), l.data() + m);
  std::vector<c_float> u_(u.data(), u.data() + m);

  // Stack-allocated csc descriptors (nz = -1 marks compressed-column form), so
  // there is no per-call heap allocation to free.
  csc Pcsc;
  Pcsc.nzmax = static_cast<c_int>(Px.size());
  Pcsc.m = n;
  Pcsc.n = n;
  Pcsc.p = Pp.data();
  Pcsc.i = Pi.data();
  Pcsc.x = Px.data();
  Pcsc.nz = -1;

  csc Acsc;
  Acsc.nzmax = static_cast<c_int>(Ax.size());
  Acsc.m = m;
  Acsc.n = n;
  Acsc.p = Ap.data();
  Acsc.i = Ai.data();
  Acsc.x = Ax.data();
  Acsc.nz = -1;

  OSQPData data;
  data.n = n;
  data.m = m;
  data.P = &Pcsc;
  data.A = &Acsc;
  data.q = q_.data();
  data.l = l_.data();
  data.u = u_.data();

  OSQPSettings settings;
  osqp_set_default_settings(&settings);
  settings.verbose = 0;
  settings.eps_abs = 1e-8;
  settings.eps_rel = 1e-8;
  settings.polish = 1;
  settings.max_iter = 8000;

  OSQPWorkspace* work = nullptr;
  const c_int flag = osqp_setup(&work, &data, &settings);
  if (flag != 0 || work == nullptr) {
    r.success = false;
    return r;
  }

  osqp_solve(work);
  r.status_val = static_cast<int>(work->info->status_val);
  r.success = (work->info->status_val == OSQP_SOLVED ||
               work->info->status_val == OSQP_SOLVED_INACCURATE);
  if (work->solution && work->solution->x) {
    for (c_int idx = 0; idx < n; ++idx)
      r.z[static_cast<int>(idx)] = work->solution->x[idx];
  }

  osqp_cleanup(work);
  return r;
}

}  // namespace autoland
