#pragma once
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include "autoland/config.hpp"
#include "autoland/control_affine_model.hpp"
#include "autoland/qp_solver.hpp"
#include "autoland/types.hpp"

// =============================================================================
// Control Barrier Function (CBF) safety filter -- BASIC CBF, QP-backed (OSQP).
//
// At each step the filter makes the minimal deviation from the nominal control
// that keeps every barrier's forward-invariance condition satisfied:
//
//     u* = argmin_u || u - u_nom ||^2_W
//          s.t.  Lf h(x) + Lg h(x) u + alpha(h(x)) >= 0     for each barrier
//                u_min <= u <= u_max
//
// With a control-affine model xdot = f(x) + g(x) u (provided via
// ControlAffineModel), the Lie derivatives are Lf h = grad_h . f and
// Lg h = grad_h . g, so each barrier is one linear inequality in u. The QP is
// assembled into OSQP standard form and solved through the QPSolver interface.
//
// FIRST PASS SCOPE (longitudinal, actuator-level):
//   * Decision controls are the longitudinal subset {delta_e, delta_T}; the
//     lateral controls pass through from the nominal command.
//   * Each barrier is HARD or SOFT (Barrier::hard). A soft barrier carries a
//     penalized non-negative slack, so its row can never make the QP infeasible
//     (best-effort). A hard barrier has no slack: its row must hold exactly.
//   * If the hard set is infeasible (e.g. a high-relative-degree barrier with
//     Lg h ~= 0 that is already violated -- the descent-rate barrier under
//     direct elevator), the filter re-solves with every barrier softened so it
//     still returns a sane, in-box control; lastFeasibilityRecovery() reports
//     that a guarantee was unavailable on that step. Authoritative hard
//     enforcement of the descent barrier ultimately needs an HOCBF (see
//     documentation/water_landing_cbf.md A.3).
// =============================================================================
namespace autoland {

// A scalar barrier h(x) >= 0 defines the safe set, with its state gradient.
struct Barrier {
  std::string name;
  std::function<double(const StateVec&)> h;            // h(x)
  std::function<Eigen::RowVectorXd(const StateVec&)> grad;  // dh/dx (1 x NX)
  // class-K function alpha(.) applied to h; default linear alpha(h)=k*h.
  std::function<double(double)> alpha = [](double hv) { return 1.0 * hv; };
  // HARD vs SOFT enforcement. A soft barrier carries a penalized slack so its
  // CBF row can never make the QP infeasible (best-effort). A hard barrier has
  // NO slack: its row  Lf h + Lg h u + alpha(h) >= 0  must hold exactly. If the
  // hard set of constraints is infeasible (e.g. a high-relative-degree barrier
  // with Lg h ~= 0 that is already being violated), the filter re-solves with
  // every barrier softened so it still returns a sane, in-box control rather
  // than failing -- see CBFFilter::filter and lastFeasibilityRecovery().
  bool hard{false};
};

// Filter tuning. ctrl_idx selects which virtual controls are QP decision
// variables (the rest pass through); ctrl_weights are the diagonal of W over
// those controls; slack_penalty is rho on each barrier's slack^2.
struct CBFConfig {
  std::vector<int> ctrl_idx{DE, DT};        // longitudinal subset by default
  std::vector<double> ctrl_weights{1.0, 1.0};
  double slack_penalty{1.0e4};
};

class CBFFilter {
 public:
  CBFFilter();
  CBFFilter(CBFConfig cfg, SurfaceLimits limits,
            std::unique_ptr<QPSolver> solver = nullptr);

  // Filter the nominal control. model supplies the control-affine (f,g) used to
  // form the Lie derivatives. Returns the safe control. When disabled or with no
  // barriers, passes u_nom through unchanged.
  CtrlVec filter(const CtrlVec& u_nom, const StateVec& x,
                 const ControlAffineModel& model,
                 const std::vector<Barrier>& barriers) const;

  bool enabled() const { return enabled_; }
  void setEnabled(bool e) { enabled_ = e; }
  const CBFConfig& config() const { return cfg_; }

  // True if the most recent filter() call could not satisfy the HARD barriers
  // and had to re-solve with every barrier softened (best-effort). A guarantee
  // was not available on that step.
  bool lastFeasibilityRecovery() const { return recovered_; }

 private:
  bool enabled_{true};
  CBFConfig cfg_;
  SurfaceLimits lim_;
  std::unique_ptr<QPSolver> solver_;
  mutable bool recovered_{false};

  // [lo, hi] deflection bounds for a virtual control index, from SurfaceLimits.
  std::pair<double, double> ctrlBounds(int ctrl) const;
};

}  // namespace autoland
