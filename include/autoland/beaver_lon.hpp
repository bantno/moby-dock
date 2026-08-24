#pragma once
#include <array>
#include <cmath>

#include "autoland/beaver_aero.hpp"
#include "autoland/lon_augmented.hpp"  // NXA, LonState enums, LonStateVec

// =============================================================================
// Beaver longitudinal augmented model with SINGLE-INTEGRATOR power.
//
// This is the propulsion architecture that lets ENGINE POWER enter the impact-
// load barrier at RELATIVE DEGREE 2 -- the same degree as the elevator -- so the
// hard impact row is enforced by BOTH actuators (a uniform degree-2, two-input
// HOCBF) instead of the elevator alone. Full justification and the relative-
// degree derivation are in documentation/water_landing_cbf_math.md section 2 and
// water_landing_cbf_design.md section 6.
//
//   State   X = [h, V, gamma, theta, q, P, -]        (reuses NXA=7; index 5 = P)
//   Control U = [delta_e, u_P]                        (u_P = Pdot, the power rate)
//
// vs. the AHAB path (lon_augmented.hpp): thrust there is a body-x FORCE augmented
// with TWO integrators (T, Tdot; control Tddot), which puts thrust at degree 3 to
// the energy barrier and degree 4 (dropped) to the impact barrier. Here power is
// a SINGLE-integrator state P (control u_P = Pdot); because power enters the aero
// through dpt in the force channel (CX_dpt/CZ_dpt -> V_dot/gamma_dot -> sink), one
// integrator places u_P at degree 2 on the impact barrier, uniform with the
// elevator's degree-2 moment channel. See beaver_aero.hpp for the dpt propulsion.
//
// CONVENTIONS: body-axis Beaver coefficients rotated into path axes exactly as
// LonDrift does (verified consistent). The elevator is routed through the pitch
// MOMENT only (Cm_de); the Beaver's direct elevator-lift term Cz_de is dropped
// from the drift/g so the elevator stays cleanly degree 2 (same modeling choice
// as the AHAB path). Rate nondim is the Beaver convention (qhat = q*c/V).
// =============================================================================
namespace autoland {

// Power state / control indices within the reused NXA=7 layout.
enum BeaverLonExtra { BLP = 5 };   // engine-power state P (index 5)
enum BeaverLonCtrl { BLDE = 0, BLUP = 1 };  // U = [delta_e, u_P = Pdot]

// Templated drift f(X) with U = 0 (so delta_e = 0, u_P = 0). Element type is
// double (plant/eval) or a Taylor series (exact Lie derivatives via lie_taylor).
struct BeaverLonDrift {
  double rho{1.112};    // air density [kg/m^3] (constant, config)
  double g{9.81};       // [m/s^2]
  double flap{0.0};     // flap deflection df [rad]
  BeaverAeroCoef k{};   // verified LR-556 coefficients

  template <class T>
  std::array<T, NXA> operator()(const std::array<T, NXA>& X) const {
    using std::cos;
    using std::sin;
    const T V = X[LV], gam = X[LGAM], th = X[LTH], q = X[LQ], P = X[BLP];
    const T alpha = th - gam;

    // Propulsion: dpt (normalized total-pressure rise) is a function of the power
    // STATE P and V (beaverDpt, inlined for the Taylor type). This V-dependence is
    // exactly why the force-channel power authority scales ~1/V (grows at low
    // speed), complementing the elevator's ~rho V^2 (shrinks at low speed).
    // Keep every division Taylor/Taylor (the Taylor type has no Taylor/scalar
    // overload): fold the scalar 2*191.18/rho into a double, then divide by V^3.
    const T dpt = 0.08696 + (191.18 * 2.0 / rho) * P / (V * V * V);
    const T qhat = q * BeaverGeom::c / V;  // Beaver convention: q*c/V
    const T qbar = (0.5 * rho) * V * V;

    // Longitudinal body-axis coefficients (beta = 0, de = 0, da = dr = 0).
    const T CX = k.Cx0 + k.Cx_a * alpha + k.Cx_a2 * alpha * alpha +
                 k.Cx_a3 * alpha * alpha * alpha + k.Cx_q * qhat +
                 k.Cx_df * flap + k.Cx_adf * alpha * flap + k.Cx_dpt * dpt +
                 k.Cx_dpt2a * dpt * dpt * alpha;
    const T CZ = k.Cz0 + k.Cz_a * alpha + k.Cz_a3 * alpha * alpha * alpha +
                 k.Cz_q * qhat + k.Cz_df * flap + k.Cz_adf * alpha * flap +
                 k.Cz_dpt * dpt;
    const T Cm = k.Cm0 + k.Cm_a * alpha + k.Cm_a2 * alpha * alpha +
                 k.Cm_q * qhat + k.Cm_df * flap + k.Cm_dpt * dpt;

    const T sa = sin(alpha), ca = cos(alpha);
    // Body -> path rotation (matches LonDrift's convention; thrust is inside CX).
    const T Fpar = (CX * ca + CZ * sa) * (qbar * BeaverGeom::S);  // along velocity
    const T Lift = (CX * sa - CZ * ca) * (qbar * BeaverGeom::S);  // up
    const T M = Cm * (qbar * BeaverGeom::S * BeaverGeom::c);      // pitch (no zcp)

    const T sg = sin(gam), cg = cos(gam);
    std::array<T, NXA> dX;
    dX[LH] = V * sg;
    dX[LV] = Fpar * (1.0 / BeaverGeom::mass) - g * sg;
    dX[LGAM] = Lift / (BeaverGeom::mass * V) - g * cg / V;
    dX[LTH] = q;
    dX[LQ] = M * (1.0 / BeaverGeom::Iy);
    dX[BLP] = T(0.0);  // Pdot = u_P enters through g (single integrator)
    dX[6] = T(0.0);    // spare slot (unused in this model)
    return dX;
  }
};

// Control-matrix columns for U = [delta_e, u_P].
//   delta_e: pitch MOMENT only (Cm_de) -> the LQ (q_dot) row. Degree 2 to any
//            theta-dependent barrier (theta <- q <- de).
//   u_P    : the power integrator -> the BLP (P_dot) row. Degree 2 to any V/gamma-
//            dependent barrier via the FORCE channel (V_dot/gamma_dot depend on P
//            through dpt).
inline std::array<double, NXA> beaverGColDe(double rho, double V,
                                            const BeaverAeroCoef& k = {}) {
  const double qbar = 0.5 * rho * V * V;
  std::array<double, NXA> col{};
  col[LQ] = qbar * BeaverGeom::S * BeaverGeom::c * k.Cm_de / BeaverGeom::Iy;
  return col;
}
inline std::array<double, NXA> beaverGColPower() {
  std::array<double, NXA> col{};
  col[BLP] = 1.0;
  return col;
}

}  // namespace autoland
