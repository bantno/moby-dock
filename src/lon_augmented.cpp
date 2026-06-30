#include "autoland/lon_augmented.hpp"

#include "autoland/stall_model.hpp"

namespace autoland {

AeroLocal makeAeroLocal(const AeroTable& table, const Mixing& mixing,
                        const AircraftConfig& cfg, double V, double alpha) {
  AeroLocal a;
  a.Sref = table.Sref();
  a.cref = table.Cref();
  a.a_sound = cfg.env.a_sound;
  a.mass = cfg.inertia.mass;
  a.Iyy = cfg.inertia.Iyy;
  a.rho = cfg.env.rho;
  a.g = cfg.env.g;
  a.zcp = cfg.thrust.zcp;
  a.parasite_CD0 = cfg.parasite_CD0;

  const double mach = V / cfg.env.a_sound;
  const AeroLookup L = table.lookup(alpha, 0.0, mach);
  const AeroDerivs& d = L.d;

  // off_c folds in the constant (angle/Mach/beta) reference terms so that the
  // EOM can use live alpha/mach/qhat:  C = base + d_alpha*(alpha-aref) + ...
  auto off = [&](int c) {
    return d.base[c] + d.d_alpha[c] * (-L.alpha_ref) +
           d.d_mach[c] * (-L.mach_ref) + d.d_beta[c] * (-L.beta_ref);
  };
  a.off_CFx = off(CFX); a.dAlpha_CFx = d.d_alpha[CFX]; a.dMach_CFx = d.d_mach[CFX]; a.dQ_CFx = d.d_q[CFX];
  a.off_CFz = off(CFZ); a.dAlpha_CFz = d.d_alpha[CFZ]; a.dMach_CFz = d.d_mach[CFZ]; a.dQ_CFz = d.d_q[CFZ];
  a.off_CMy = off(CMY); a.dAlpha_CMy = d.d_alpha[CMY]; a.dMach_CMy = d.d_mach[CMY]; a.dQ_CMy = d.d_q[CMY];

  // Elevator pitch-moment derivative through the mixing matrix: virtual delta_e
  // (column 0) -> physical control-group deflections -> sum of group derivs.
  const Mat& M = mixing.matrix();
  double dDe = 0.0;
  for (int gidx = 0; gidx < static_cast<int>(d.d_ctrl.size()) && gidx < M.rows(); ++gidx)
    dDe += d.d_ctrl[gidx][CMY] * M(gidx, 0);
  a.dDe_CMy = dDe;

  // Viscous-stall overlay: freeze the NACA 4414 blend at the eval point alpha as
  // local affine models {off + slope*alpha}, folding -slope*alpha into the offset
  // (mirrors the off() lambda above). LonDrift forms (1-w)*C_vspaero + w*C_post.
  // severity scales the blend weight w (depth knob: 0 => no stall, 1 => full).
  // Default off => never touched => identical to the inviscid deck.
  a.stall_on = cfg.stall.enabled;
  if (a.stall_on) {
    const StallBlend s = stallLookup(alpha);
    const double sev = cfg.stall.severity;
    a.dAlpha_w = sev * s.w_da;      a.off_w = sev * s.w - a.dAlpha_w * alpha;
    a.dAlpha_CLp = s.CLpost_da;     a.off_CLp = s.CLpost - a.dAlpha_CLp * alpha;
    a.dAlpha_CDp = s.CDpost_da;     a.off_CDp = s.CDpost - a.dAlpha_CDp * alpha;
    a.dAlpha_CMp = s.CMpost_da;     a.off_CMp = s.CMpost - a.dAlpha_CMp * alpha;
  }
  return a;
}

LonGMat gMatrix(const AeroLocal& a, const LonStateVec& X) {
  const double V = X[LV];
  const double qbar = 0.5 * a.rho * V * V;
  LonGMat G;
  G.setZero();
  // delta_e -> pitch moment -> qdot.  (rho V^2 S cbar / 2 Iyy) * C_mde
  G(LQ, LDE) = qbar * a.Sref * a.cref * a.dDe_CMy / a.Iyy;
  // Tddot is the Tdot-integrator input.
  G(LTDOT, LTDDOT) = 1.0;
  return G;
}

}  // namespace autoland
