#pragma once
#include <string>
#include <vector>
#include "autoland/types.hpp"

// =============================================================================
// VSPAero .stab parser + (alpha, beta, Mach) interpolation.
//
// FILE FORMAT (verified against data/example.stab, AHAB V-tail flying boat):
//   The file is a sequence of blocks separated by lines of '*'. Each block is
//   one operating point of an (alpha, beta) [x Mach] sweep and contains:
//     * a reference-quantity table ("Sref_ ... Vinf_ ...")  -- Name Value Units
//     * a "Case" table: base + finite-difference perturbation rows. The trailing
//       rows (after Mach) are the named control groups, in control-group order.
//     * a "Coef ... Total Alpha Beta p q r Mach U ConGrp_1 .. ConGrp_N" table of
//       analytic derivatives, one row per coefficient.
//   We consume only the body-axis rows CFx,CFy,CFz,CMx(=Cl),CMy(=Cm),CMz(=Cn).
//
// UNITS / CONVENTIONS (read from the file, not assumed):
//   * Lengths (Sref,Cref,Bref,Xcg,...) are model units. Confirmed SI metres for
//     this model (Mach 0.059 == Vinf 20 / a, forcing Vinf in m/s). Air density
//     is NOT taken from the file's Rho_ (a stale imperial default); it comes
//     from config -- see config.hpp / dynamics.hpp.
//   * AoA_ and Beta_ are in degrees in the file; converted to RADIANS on read.
//   * Derivative columns are, per the file's own sub-header:
//        wrt Alpha,Beta : per radian        (stored as-is; angle derivs are
//                                             therefore already per-radian)
//        wrt p,q,r      : per NONDIMENSIONAL rate. VSPAero uses
//                            phat = p*Bref/(2V), qhat = q*Cref/(2V),
//                            rhat = r*Bref/(2V).
//                         (Verified numerically against the perturbation rows.)
//                         The nondimensionalization is applied in Dynamics, not
//                         here -- this table stores the raw d/d(rate-hat) values.
//        wrt Mach       : per unit Mach
//        wrt U          : per unit freestream speed [1/(m/s)] (parsed, unused by
//                         default -- dynamic pressure already carries speed).
//        wrt ConGrp_k   : per radian of control-group deflection.
// =============================================================================
namespace autoland {

// One set of body-axis derivatives at a single sweep node (raw, as parsed).
struct AeroDerivs {
  CoefVec base{CoefVec::Zero()};
  CoefVec d_alpha{CoefVec::Zero()};
  CoefVec d_beta{CoefVec::Zero()};
  CoefVec d_p{CoefVec::Zero()};   // wrt phat = p*Bref/2V
  CoefVec d_q{CoefVec::Zero()};   // wrt qhat = q*Cref/2V
  CoefVec d_r{CoefVec::Zero()};   // wrt rhat = r*Bref/2V
  CoefVec d_mach{CoefVec::Zero()};
  CoefVec d_u{CoefVec::Zero()};
  std::vector<CoefVec> d_ctrl;    // one CoefVec per control group, group order
};

// Result of an interpolated query. alpha_ref/beta_ref/mach_ref are the points
// about which the caller should apply the (angle/Mach) derivative buildup:
//   * swept axis, query inside the grid: ref == query value, so the residual
//     term is zero and the interpolated base already carries the dependence;
//   * degenerate single-node axis (e.g. a single-Mach table): ref == node
//     value, so the derivative provides a local linear correction;
//   * query OFF the grid: ref == the nearest boundary node value, so the
//     residual term becomes a first-order linear extrapolation from that node.
// See Dynamics::aeroCoeffs().
struct AeroLookup {
  AeroDerivs d;
  double alpha_ref{0.0};
  double beta_ref{0.0};
  double mach_ref{0.0};
  bool clamped{false};  // query fell outside the grid (linearly extrapolated)
};

class AeroTable {
 public:
  // Parse a VSPAero .stab file. Throws std::runtime_error on a malformed file.
  static AeroTable fromFile(const std::string& path);

  // Trilinear interpolation at (alpha[rad], beta[rad], Mach). Off-grid queries
  // are LINEARLY EXTRAPOLATED from the nearest boundary node using that node's
  // analytic angle/Mach derivative (realized via alpha_ref/beta_ref/mach_ref in
  // the AeroLookup; see Dynamics::aeroCoeffs), and a rate-limited warning is
  // emitted. The base is still taken at the boundary node (not extrapolated);
  // only the first-order derivative term extends beyond the grid.
  AeroLookup lookup(double alpha, double beta, double mach) const;

  // Reference geometry (model units == metres for this model).
  double Sref() const { return sref_; }
  double Cref() const { return cref_; }
  double Bref() const { return bref_; }
  double Xcg() const { return xcg_; }
  double Ycg() const { return ycg_; }
  double Zcg() const { return zcg_; }

  const std::vector<std::string>& controlGroupNames() const { return cg_names_; }
  int numControlGroups() const { return static_cast<int>(cg_names_.size()); }

  // Grid axes (alpha/beta in radians).
  const std::vector<double>& alphaGrid() const { return alphas_; }
  const std::vector<double>& betaGrid() const { return betas_; }
  const std::vector<double>& machGrid() const { return machs_; }

 private:
  double sref_{0}, cref_{0}, bref_{0}, xcg_{0}, ycg_{0}, zcg_{0};
  std::vector<std::string> cg_names_;
  std::vector<double> alphas_, betas_, machs_;  // sorted unique axis values
  std::vector<AeroDerivs> nodes_;               // size = na*nb*nm, row-major

  int idx(int ia, int ib, int im) const {
    return (ia * static_cast<int>(betas_.size()) + ib) *
               static_cast<int>(machs_.size()) + im;
  }
};

}  // namespace autoland
