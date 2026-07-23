#include "autoland/dynamics.hpp"

#include <algorithm>
#include <cmath>

namespace autoland {

Dynamics::Dynamics(const AeroTable& table, const Mixing& mixing,
                   const AircraftConfig& cfg)
    : table_(table), mixing_(mixing), cfg_(cfg) {
  // Default thrust model: T = throttle * (T_static - k_v * V), clamped >= 0.
  const double Tstatic = cfg_.thrust.T_static;
  const double kv = cfg_.thrust.k_v;
  thrust_ = [Tstatic, kv](double throttle, double V) {
    return std::max(0.0, throttle * (Tstatic - kv * V));
  };
}

CoefVec Dynamics::aeroCoeffs(const StateVec& x, const CtrlVec& u) const {
  // NOTE: 'Vt' (true airspeed), NOT 'V' -- 'V' is the State::V enum index and
  // must not be shadowed by a local, or xd[V] would index at the airspeed value.
  const double ub = x[U], vb = x[V], wb = x[W];
  const double Vt = std::max(1e-3, std::sqrt(ub * ub + vb * vb + wb * wb));
  const double alpha = std::atan2(wb, ub);
  const double beta = std::asin(std::clamp(vb / Vt, -1.0, 1.0));
  const double mach = Vt / cfg_.env.a_sound;

  const AeroLookup L = table_.lookup(alpha, beta, mach);
  const AeroDerivs& d = L.d;

  // Nondimensional rates (VSPAero convention):
  //   phat = p*Bref/2V, qhat = q*Cref/2V, rhat = r*Bref/2V.
  const double bref = table_.Bref(), cref = table_.Cref();
  const double phat = x[P] * bref / (2.0 * Vt);
  const double qhat = x[Q] * cref / (2.0 * Vt);
  const double rhat = x[R] * bref / (2.0 * Vt);

  // Mixing: virtual controls -> physical control-group deflections [rad].
  const Vec defl = mixing_.apply(u[DE], u[DA], u[DR]);

  // Coefficient buildup IN THE FILE'S FRAME (VSPAero X-fwd, Y-right, Z-UP). The
  // angle/Mach residual terms vanish on swept axes (alpha_ref==alpha) and give a
  // local linear correction on degenerate single-node axes (single-Mach table).
  //
  // Yaw-rate convention: the file's r-derivatives are in the Z-up frame, so the
  // standard-axes yaw rate enters as (-rhat) here (see the frame note below).
  CoefVec C = d.base;
  C += d.d_alpha * (alpha - L.alpha_ref);
  C += d.d_beta * (beta - L.beta_ref);
  C += d.d_mach * (mach - L.mach_ref);
  C += d.d_p * phat + d.d_q * qhat + d.d_r * (-rhat);
  // d.d_u (wrt freestream speed) intentionally omitted: dynamic pressure
  // already carries the speed dependence. (Parsed and available for analysis.)
  for (int g = 0; g < static_cast<int>(d.d_ctrl.size()) && g < defl.size(); ++g)
    C += d.d_ctrl[g] * defl[g];

  // -------------------------------------------------------------------------
  // FRAME TRANSFORM: file (X-fwd, Y-right, Z-UP) -> standard flight-dynamics
  // body axes (X-fwd, Y-right, Z-DOWN). The .stab reports CFz positive-UP
  // (CFz ~= +CL in the data). The per-coefficient sign vector below is fixed by
  // six independent physical stability constraints (Cm_alpha<0, lift up,
  // Cl_beta<0, Cn_beta>0, Cy_beta<0, Cn_r<0); see linear_model.hpp.
  //   [CFx, CFy, CFz, Cl=CMx, Cm=CMy, Cn=CMz] -> [+1,+1,-1,-1,+1,-1]
  // The yaw-rate variable flip was already applied to its term above.
  // -------------------------------------------------------------------------
  static const CoefVec kFrameSign =
      (CoefVec() << 1.0, 1.0, -1.0, -1.0, 1.0, -1.0).finished();
  C = C.cwiseProduct(kFrameSign);

  // Parasite drag (viscous + hull) absent from the inviscid .stab: a constant
  // body-x drag increment in standard axes. Crude (no alpha rotation) but
  // explicit and configurable; see AircraftConfig::parasite_CD0.
  C[CFX] -= cfg_.parasite_CD0;

  return C;
}

StateVec Dynamics::xdot(const StateVec& x, const CtrlVec& u) const {
  const double ub = x[U], vb = x[V], wb = x[W];
  const double p = x[P], q = x[Q], r = x[R];
  const double phi = x[PHI], theta = x[THETA], psi = x[PSI];

  // 'Vt' (airspeed) -- do NOT name this 'V'; see aeroCoeffs note above.
  const double Vt = std::max(1e-3, std::sqrt(ub * ub + vb * vb + wb * wb));
  const double qbar = 0.5 * cfg_.env.rho * Vt * Vt;
  const double S = table_.Sref(), b = table_.Bref(), c = table_.Cref();

  const CoefVec C = aeroCoeffs(x, u);

  // Body-axis aerodynamic forces and moments (about the .stab moment reference).
  double Fx = qbar * S * C[CFX];
  double Fy = qbar * S * C[CFY];
  double Fz = qbar * S * C[CFZ];
  double Lm = qbar * S * b * C[CMX];  // roll  (Cl = CMx)
  double Mm = qbar * S * c * C[CMY];  // pitch (Cm = CMy)
  double Nm = qbar * S * b * C[CMZ];  // yaw   (Cn = CMz)

  // Thrust along body x. zcp != 0 couples thrust into a pitching moment.
  const double T = thrust_(u[DT], Vt);
  Fx += T;
  Mm += T * cfg_.thrust.zcp;  // +z below c.g. -> pitch coupling (TODO: confirm)

  // TODO: transfer moments to the true c.g. if cg_offset != 0:
  //   M_cg = M_ref + r_cg x F.  Default cg_offset = 0 (c.g. at .stab reference).

  const InertiaParams& I = cfg_.inertia;
  const double m = I.mass, g = cfg_.env.g;

  // Gravity in body axes.
  const double gx = -g * std::sin(theta);
  const double gy = g * std::cos(theta) * std::sin(phi);
  const double gz = g * std::cos(theta) * std::cos(phi);

  // Translational dynamics (body axes).
  const double udot = r * vb - q * wb + Fx / m + gx;
  const double vdot = p * wb - r * ub + Fy / m + gy;
  const double wdot = q * ub - p * vb + Fz / m + gz;

  // Rotational dynamics. Solve the coupled p-r system with the Ixz product term.
  const double Gamma = I.Ixx * I.Izz - I.Ixz * I.Ixz;
  // Roll:  Ixx pdot - Ixz rdot = L + Ixz p q - (Izz - Iyy) q r
  // Yaw:   Izz rdot - Ixz pdot = N - Ixz q r - (Iyy - Ixx) p q
  // (from hdot + omega x h = M with h = [Ixx p - Ixz r, Iyy q, Izz r - Ixz p];
  // expanding the Gamma solve below reproduces Stevens & Lewis.)
  const double rhsL = Lm - (I.Izz - I.Iyy) * q * r + I.Ixz * p * q;
  const double rhsN = Nm - (I.Iyy - I.Ixx) * p * q - I.Ixz * q * r;
  const double pdot = (I.Izz * rhsL + I.Ixz * rhsN) / Gamma;
  const double rdot = (I.Ixx * rhsN + I.Ixz * rhsL) / Gamma;
  const double qdot =
      (Mm - (I.Ixx - I.Izz) * p * r - I.Ixz * (p * p - r * r)) / I.Iyy;

  // Euler kinematics.
  const double ct = std::cos(theta), st = std::sin(theta), tt = std::tan(theta);
  const double cp = std::cos(phi), sp = std::sin(phi);
  const double phidot = p + (q * sp + r * cp) * tt;
  const double thetadot = q * cp - r * sp;
  const double psidot = (q * sp + r * cp) / ct;

  // Guidance states. Approach centerline is along North; psi measured from North.
  const double cps = std::cos(psi), sps = std::sin(psi);
  // Altitude rate (positive up):
  const double hdot = ub * st - vb * sp * ct - wb * cp * ct;
  // Cross-track (east) velocity:
  const double ydot = ub * ct * sps + vb * (sp * st * sps + cp * cps) +
                      wb * (cp * st * sps - sp * cps);

  StateVec xd;
  xd[U] = udot;   xd[V] = vdot;   xd[W] = wdot;
  xd[P] = pdot;   xd[Q] = qdot;   xd[R] = rdot;
  xd[PHI] = phidot; xd[THETA] = thetadot; xd[PSI] = psidot;
  xd[H] = hdot;   xd[Y] = ydot;
  return xd;
}

StateVec rk4Step(const std::function<StateVec(const StateVec&, const CtrlVec&)>& f,
                 const StateVec& x, const CtrlVec& u, double dt) {
  const StateVec k1 = f(x, u);
  const StateVec k2 = f(x + 0.5 * dt * k1, u);
  const StateVec k3 = f(x + 0.5 * dt * k2, u);
  const StateVec k4 = f(x + dt * k3, u);
  return x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

}  // namespace autoland
