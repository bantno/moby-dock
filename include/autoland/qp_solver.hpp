#pragma once
#include "autoland/types.hpp"

// =============================================================================
// QPSolver: a thin interface over a convex quadratic-program solver, in the
// standard OSQP form
//
//        min  1/2 z^T P z + q^T z      s.t.   l <= A z <= u
//
// with P symmetric positive semidefinite. The CBF filter assembles its safety
// QP into (P, q, A, l, u) and calls solve(); keeping the solver behind this
// interface means OSQP can be replaced (e.g. a dense fallback) without touching
// cbf.cpp. One-sided constraints use a large sentinel (>= 1e30) for +/-inf.
// =============================================================================
namespace autoland {

struct QPResult {
  Vec z;              // primal solution (zeros on failure)
  bool success{false};
  int status_val{0};  // solver-specific status code (for diagnostics)
};

class QPSolver {
 public:
  virtual ~QPSolver() = default;
  virtual QPResult solve(const Mat& P, const Vec& q, const Mat& A,
                         const Vec& l, const Vec& u) const = 0;
};

// OSQP-backed solver. Takes dense P/A (P read as upper-triangular) and converts
// to OSQP's compressed-sparse-column form internally.
class OsqpSolver : public QPSolver {
 public:
  QPResult solve(const Mat& P, const Vec& q, const Mat& A, const Vec& l,
                 const Vec& u) const override;
};

}  // namespace autoland
