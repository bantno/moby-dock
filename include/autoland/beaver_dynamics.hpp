#pragma once
#include <Eigen/Dense>
#include <cmath>

#include "autoland/beaver_aero.hpp"
#include "autoland/config.hpp"
#include "autoland/trim.hpp"  // TrimResult
#include "autoland/types.hpp"

// =============================================================================
// BeaverDynamics: the DHC-2 Beaver 6-DOF nonlinear plant, evaluated DIRECTLY
// from the verified LR-556/FDC polynomials (beaver_aero.hpp) -- no table
// export, no interpolation. Same state/control layout and rigid-body EOM
// structure as Dynamics (body axes, full Ixz coupling, wind through the
// air-relative velocity only), so it drops into SixDofSim as an alternate
// plant.
//
// PROPULSION: the Beaver has no separate thrust force -- power enters the aero
// coefficients through the slipstream coefficient dpt (see beaver_aero.hpp).
// The virtual throttle u[DT] in [0,1] maps AFFINELY onto the manifold pressure
// pz in [pz_idle, pz_max] at a fixed propeller speed n_rpm (constant-speed
// governor on approach); pz -> engine power P (FDC eq. 3.15) -> dpt. Note
// pz below ~13-14 "Hg gives NEGATIVE power = windmilling-propeller drag.
//
// ATMOSPHERE: rho and g are FROZEN at the configured reference altitude h_ref
// (ICAO standard atmosphere + inverse-square gravity, FDC eqs. 3.24-3.31) --
// the same constant-density convention as the VSPAERO plant. x[H] does NOT
// feed back into rho, so linearizations are classic frozen-atmosphere modes.
//
// The scalar type is templated (xdotT) so trim Jacobians and the state-space
// linearization are computed EXACTLY with autodiff::dual -- no finite
// differences (see beaver_dynamics.cpp).
// =============================================================================
namespace autoland {

// --- ICAO standard atmosphere (troposphere) + gravity, FDC eqs. 3.24-3.31 ---
inline double icaoTemperature(double h) { return 288.15 - 0.0065 * h; }
inline double icaoPressure(double h) {
  const double T0 = 288.15, lam = -0.0065, R = 287.05, g0 = 9.80665;
  return 101325.0 * std::pow(1.0 + lam * h / T0, -g0 / (lam * R));
}
inline double icaoDensity(double h) {
  return icaoPressure(h) / (287.05 * icaoTemperature(h));
}
inline double gravityAt(double h) {
  const double Re = 6371020.0;
  const double f = Re / (Re + h);
  return 9.80665 * f * f;
}

// Plant configuration. Defaults: sea level, approach RPM, flaps up, and
// deflection limits representative of the DHC-2 (ASSUMED round numbers, not
// LR-556 data -- the aero model itself carries no limit information).
struct BeaverPlantConfig {
  double n_rpm{1800.0};   // propeller speed [RPM] (constant-speed governor)
  double pz_idle{5.0};    // manifold pressure at u[DT]=0 [inHg] (windmill drag)
  double pz_max{26.0};    // manifold pressure at u[DT]=1 [inHg] (Table C.1 max)
  double flap{0.0};       // fixed flap deflection df [rad]
  double h_ref{0.0};      // reference altitude for rho/g [m]
  BeaverAeroCoef aero{};  // verified LR-556 coefficients
  SurfaceLimits limits{makeDefaultLimits()};

  static SurfaceLimits makeDefaultLimits() {
    SurfaceLimits L;
    const double d = M_PI / 180.0;
    L.de_min = -25.0 * d; L.de_max = 20.0 * d;
    L.da_min = -20.0 * d; L.da_max = 20.0 * d;
    L.dr_min = -25.0 * d; L.dr_max = 25.0 * d;
    L.rate = 1.0;  // [rad/s]
    return L;
  }
};

class BeaverDynamics {
 public:
  explicit BeaverDynamics(const BeaverPlantConfig& cfg = BeaverPlantConfig{})
      : cfg_(cfg), rho_(icaoDensity(cfg.h_ref)), g_(gravityAt(cfg.h_ref)) {}

  const BeaverPlantConfig& config() const { return cfg_; }
  double rho() const { return rho_; }
  double g() const { return g_; }

  // Throttle chain accessors (doubles; the templated versions live in xdotT).
  double pzFromThrottle(double dT) const {
    return cfg_.pz_idle + dT * (cfg_.pz_max - cfg_.pz_idle);
  }
  double power(double dT) const {
    return beaverEnginePower(pzFromThrottle(dT), cfg_.n_rpm, rho_);
  }
  double dpt(double dT, double V) const {
    return beaverDpt(power(dT), rho_, V);
  }

  // Templated EOM core. x is the INERTIAL state; (ua, va, wa) are the
  // AIR-RELATIVE body velocities the aerodynamics sees (equal to x[U..W] in
  // still air). Pass x[U], x[V], x[W] for them when differentiating so the
  // aero terms pick up the state dependence.
  //
  // alpha = atan(wa/ua) (not atan2: ua > 0 over the whole flight envelope and
  // atan is smooth for autodiff); beta = asin(va/V) with |va| <= V by
  // construction.
  template <class T>
  Eigen::Matrix<T, NX, 1> xdotT(const Eigen::Matrix<T, NX, 1>& x, const T& ua,
                                const T& va, const T& wa,
                                const Eigen::Matrix<T, NU, 1>& u) const {
    using std::asin;
    using std::atan;
    using std::cos;
    using std::sin;
    using std::sqrt;
    using std::tan;

    const T Vt = sqrt(ua * ua + va * va + wa * wa + T(1e-12));
    const T alpha = atan(wa / ua);
    const T beta = asin(va / Vt);

    // Propulsion chain: throttle -> pz -> P -> dpt (see class comment).
    const T pz = cfg_.pz_idle + u[DT] * (cfg_.pz_max - cfg_.pz_idle);
    const T P = beaverEnginePowerT<T>(pz, cfg_.n_rpm, rho_);
    const T dpt = beaverDptT<T>(P, rho_, Vt);

    const std::array<T, 6> C = beaverAeroCoeffsT<T>(
        alpha, beta, x[autoland::P], x[Q], x[R], Vt, u[DE], u[DA], u[DR],
        T(cfg_.flap), dpt, cfg_.aero);

    const T qbar = (0.5 * rho_) * Vt * Vt;
    const double S = BeaverGeom::S, b = BeaverGeom::b, c = BeaverGeom::c;
    const T Fx = qbar * S * C[0];
    const T Fy = qbar * S * C[1];
    const T Fz = qbar * S * C[2];
    const T Lm = qbar * S * b * C[3];
    const T Mm = qbar * S * c * C[4];
    const T Nm = qbar * S * b * C[5];

    const T &ub = x[U], &vb = x[V], &wb = x[W];
    const T &p = x[autoland::P], &q = x[Q], &r = x[R];
    const T &phi = x[PHI], &theta = x[THETA], &psi = x[PSI];
    const double m = BeaverGeom::mass;

    // Gravity in body axes.
    const T gx = -g_ * sin(theta);
    const T gy = g_ * cos(theta) * sin(phi);
    const T gz = g_ * cos(theta) * cos(phi);

    // Translational dynamics (inertial body velocities; wind enters only
    // through the aero above -- cf. Dynamics::xdot).
    const T udot = r * vb - q * wb + Fx / T(m) + gx;
    const T vdot = p * wb - r * ub + Fy / T(m) + gy;
    const T wdot = q * ub - p * vb + Fz / T(m) + gz;

    // Rotational dynamics with the Ixz product term (Stevens & Lewis; same
    // Gamma solve as Dynamics::xdot).
    const double Ix = BeaverGeom::Ix, Iy = BeaverGeom::Iy,
                 Iz = BeaverGeom::Iz, Ixz = BeaverGeom::Ixz;
    const double Gamma = Ix * Iz - Ixz * Ixz;
    const T rhsL = Lm - (Iz - Iy) * q * r + Ixz * p * q;
    const T rhsN = Nm - (Iy - Ix) * p * q - Ixz * q * r;
    const T pdot = (Iz * rhsL + Ixz * rhsN) / T(Gamma);
    const T rdot = (Ix * rhsN + Ixz * rhsL) / T(Gamma);
    const T qdot = (Mm - (Ix - Iz) * p * r - Ixz * (p * p - r * r)) / T(Iy);

    // Euler kinematics.
    const T ct = cos(theta), st = sin(theta), tt = tan(theta);
    const T cp = cos(phi), sp = sin(phi);
    const T phidot = p + (q * sp + r * cp) * tt;
    const T thetadot = q * cp - r * sp;
    const T psidot = (q * sp + r * cp) / ct;

    // Guidance rows (inertial kinematics; identical to Dynamics::xdot).
    const T cps = cos(psi), sps = sin(psi);
    const T hdot = ub * st - vb * sp * ct - wb * cp * ct;
    const T ydot = ub * ct * sps + vb * (sp * st * sps + cp * cps) +
                   wb * (cp * st * sps - sp * cps);

    Eigen::Matrix<T, NX, 1> xd;
    xd[U] = udot; xd[V] = vdot; xd[W] = wdot;
    xd[autoland::P] = pdot; xd[Q] = qdot; xd[R] = rdot;
    xd[PHI] = phidot; xd[THETA] = thetadot; xd[PSI] = psidot;
    xd[H] = hdot; xd[Y] = ydot;
    return xd;
  }

  // Still-air state derivative.
  StateVec xdot(const StateVec& x, const CtrlVec& u) const {
    return xdotT<double>(x, x[U], x[V], x[W], u);
  }

  // Wind-perturbed state derivative; W_earth = (north tailwind+, east+,
  // updraft+), identical convention and rotation to Dynamics::xdot.
  StateVec xdot(const StateVec& x, const CtrlVec& u,
                const Eigen::Vector3d& W_earth) const {
    const double ct = std::cos(x[THETA]), st = std::sin(x[THETA]);
    const double cp = std::cos(x[PHI]), sp = std::sin(x[PHI]);
    const double cy = std::cos(x[PSI]), sy = std::sin(x[PSI]);
    const double Wn = W_earth[0], We = W_earth[1], Wd = -W_earth[2];
    const double Wbx = ct * cy * Wn + ct * sy * We - st * Wd;
    const double Wby = (sp * st * cy - cp * sy) * Wn +
                       (sp * st * sy + cp * cy) * We + sp * ct * Wd;
    const double Wbz = (cp * st * cy + sp * sy) * Wn +
                       (cp * st * sy - sp * cy) * We + cp * ct * Wd;
    return xdotT<double>(x, x[U] - Wbx, x[V] - Wby, x[W] - Wbz, u);
  }

  // Exact still-air Jacobians A = d xdot/dx (NX x NX), B = d xdot/du (NX x NU)
  // at (x0, u0), computed with autodiff::dual (beaver_dynamics.cpp).
  void linearize(const StateVec& x0, const CtrlVec& u0, Mat& A, Mat& B) const;

 private:
  BeaverPlantConfig cfg_;
  double rho_, g_;
};

// Steady straight-flight trim at airspeed V and flight-path angle gamma,
// wings level (phi = 0), on the Beaver plant. SIX unknowns
// [alpha, beta, de, da, dr, dT] against the six body accelerations -- the
// Beaver's slipstream asymmetry (Cy0, Cl0, Cn0, and the dpt yaw term) makes
// the true trim carry small nonzero beta/da/dr, exactly as the FDC ACTRIM
// check case shows. theta = alpha + asin(sin(gamma)/cos(beta)) (exact
// wings-level relation). Newton with the EXACT autodiff Jacobian.
// The returned TrimResult carries da/dr in u[DA]/u[DR] and the sideslip in
// x[V] (body v = V sin(beta)).
TrimResult beaverTrim(const BeaverDynamics& dyn, double V, double gamma,
                      double tol = 1e-10, int max_iter = 60);

}  // namespace autoland
