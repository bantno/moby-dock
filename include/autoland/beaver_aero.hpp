#pragma once
#include <array>
#include <cmath>

// =============================================================================
// DHC-2 Beaver nonlinear aerodynamic + engine model.
//
// Source: Tjee & Mulder, "Stability and Control Derivatives of the De Havilland
// DHC-2 'Beaver' Aircraft", TU-Delft Report LR-556 (1988); nonlinear polynomial
// model as implemented in Rauw's FDC toolbox. The coefficient set below was
// VERIFIED against LR-556 Table 3 (aero model) + Table 2 (mass/inertia) on
// 2026-08-19 -- every value matches to the printed digits. See
// documentation/paper_readiness.md section 6.
//
// SCOPE / VALIDITY (LR-556 Table 3 caption): the model is flight-validated only
// within the 30-55 m/s IAS range and over the attached-flow envelope (no stall /
// post-stall aero). See paper_readiness.md for the landing-phase implication (the
// approach is validated; the near-stall flare/touchdown extrapolates below 30 m/s).
//
// CONVENTIONS:
//   * Coefficients are BODY-AXIS with STANDARD flight-dynamics signs (X fwd,
//     Y right, Z down): [CX, CY, CZ, Cl, Cm, Cn]. These map DIRECTLY onto
//     autoland's CoefVec [CFX, CFY, CFZ, CMX, CMY, CMZ] with NO frame flip
//     (unlike the VSPAERO deck, which needs kFrameSign in dynamics.cpp): the
//     Beaver deck already reports Z-down / nose-down-negative-Cm etc.
//   * Rate nondimensionalization is the Delft/Beaver convention and DIFFERS from
//     the VSPAero convention in dynamics.cpp:
//         qhat = q * c / V        (PITCH: full chord, NOT c/2V)
//         phat = p * b / (2V),  rhat = r * b / (2V)   (roll/yaw: half span)
//   * PROPULSION enters ONLY through the slipstream/thrust coefficient dpt (the
//     normalized propeller total-pressure rise). There is NO separate body-x
//     thrust force: net thrust and slipstream effects are folded into the aero
//     coefficients via the Cx_dpt / Cz_dpt / Cm_dpt (etc.) terms. The engine
//     "throttle" is therefore the manifold pressure Pz / rpm n -> power P -> dpt.
//   * Angles/deflections in RADIANS; rates in rad/s; V in m/s.
//
// OMISSION vs LR-556 (standard, matches the FDC implementation): the CY column's
// beta-dot term (Cy_bdot ~ -0.16, ~10% rel.std.dev in LR-556) is dropped -- beta-
// dot is not a state, the term is the least certain in the table, and it is
// negligible for a landing sim. All other ~60 coefficients are reproduced exactly.
// =============================================================================
namespace autoland {

// Reference geometry, mass, and inertia (LR-556 Table 1 / Table 2, verified).
struct BeaverGeom {
  static constexpr double S = 23.23;       // wing area [m^2]
  static constexpr double b = 14.63;       // span [m]
  static constexpr double c = 1.5875;      // mean aerodynamic chord [m]
  static constexpr double mass = 2288.231; // [kg]
  static constexpr double Ix = 5368.39;    // [kg m^2]
  static constexpr double Iy = 6928.93;    // [kg m^2]
  static constexpr double Iz = 11158.75;   // [kg m^2]
  static constexpr double Ixz = 117.64;    // [kg m^2]
  static constexpr double xcg = 0.5996;    // c.g. in F_M [m]
  static constexpr double zcg = -0.8851;   // c.g. in F_M [m]
};

// The verified LR-556 Table 3 polynomial coefficients. Data-driven (a struct of
// defaults) so a future calibration can override without touching the assembly.
struct BeaverAeroCoef {
  // CX
  double Cx0 = -0.03554, Cx_a = 0.00292, Cx_a2 = 5.459, Cx_a3 = -5.162,
         Cx_q = -0.6748, Cx_dr = 0.03412, Cx_df = -0.09447, Cx_adf = 1.106,
         Cx_dpt = 0.1161, Cx_dpt2a = 0.1453;
  // CY
  double Cy0 = -0.002226, Cy_b = -0.7678, Cy_p = -0.124, Cy_r = 0.3666,
         Cy_da = -0.02956, Cy_dr = 0.1158, Cy_dra = 0.5238;
  // CZ
  double Cz0 = -0.05504, Cz_a = -5.578, Cz_a3 = 3.442, Cz_q = -2.988,
         Cz_de = -0.398, Cz_deb2 = -15.93, Cz_df = -1.377, Cz_adf = -1.261,
         Cz_dpt = -0.1563;
  // Cl
  double Cl0 = 0.000591, Cl_b = -0.0618, Cl_p = -0.5045, Cl_r = 0.1695,
         Cl_da = -0.09917, Cl_dr = 0.006934, Cl_daa = -0.08269, Cl_a2dpt = -0.01406;
  // Cm
  double Cm0 = 0.09448, Cm_a = -0.6028, Cm_a2 = -2.14, Cm_q = -15.56,
         Cm_de = -1.921, Cm_b2 = 0.6921, Cm_r = -0.3118, Cm_df = 0.4072,
         Cm_dpt = -0.07895;
  // Cn
  double Cn0 = -0.003117, Cn_b = 0.006719, Cn_p = -0.1585, Cn_r = -0.1112,
         Cn_da = -0.003872, Cn_dr = -0.08265, Cn_q = 0.1595, Cn_b3 = 0.1373,
         Cn_dpt3 = -0.003026;
};

// Engine shaft power [in the model's units] from manifold pressure Pz [inHg],
// engine speed n [rpm], and air density rho [kg/m^3] (LR-556 / FDC engine model).
inline double beaverEnginePower(double Pz, double n, double rho) {
  return 0.7355 * (-326.5 + (0.00412 * (Pz + 7.4) * (n + 2010.0) +
                             (408.0 - 0.0965 * n) * (1.0 - rho / 1.225)));
}

// Normalized propeller total-pressure rise dpt = 0.08696 + 191.18 * 2P/(rho V^3).
// This is the single "throttle" quantity the aero coefficients respond to.
inline double beaverDpt(double P, double rho, double V) {
  const double Vs = (std::abs(V) < 1e-3) ? 1e-3 : V;
  return 0.08696 + 191.18 * (P * 2.0 / rho / (Vs * Vs * Vs));
}

// Body-axis force/moment coefficients [CX, CY, CZ, Cl, Cm, Cn].
//   alpha, beta [rad]; p, q, r [rad/s]; V [m/s]; de, da, dr, df [rad]; dpt [-].
inline std::array<double, 6> beaverAeroCoeffs(
    double alpha, double beta, double p, double q, double r, double V, double de,
    double da, double dr, double df, double dpt,
    const BeaverAeroCoef& k = BeaverAeroCoef{}) {
  const double Vs = (std::abs(V) < 1e-3) ? 1e-3 : V;
  const double half_b_V = BeaverGeom::b / (2.0 * Vs);
  const double phat = p * half_b_V;
  const double rhat = r * half_b_V;
  const double qhat = q * BeaverGeom::c / Vs;  // Beaver convention: q*c/V
  const double a = alpha, be = beta;

  const double CX = k.Cx0 + k.Cx_a * a + k.Cx_a2 * a * a + k.Cx_a3 * a * a * a +
                    k.Cx_q * qhat + k.Cx_dr * dr + k.Cx_df * df +
                    k.Cx_adf * a * df + k.Cx_dpt * dpt +
                    k.Cx_dpt2a * dpt * dpt * a;
  const double CY = k.Cy0 + k.Cy_b * be + k.Cy_p * phat + k.Cy_r * rhat +
                    k.Cy_da * da + k.Cy_dr * dr + k.Cy_dra * a * dr;
  const double CZ = k.Cz0 + k.Cz_a * a + k.Cz_a3 * a * a * a + k.Cz_q * qhat +
                    k.Cz_de * de + k.Cz_deb2 * de * be * be + k.Cz_df * df +
                    k.Cz_adf * a * df + k.Cz_dpt * dpt;
  const double Cl = k.Cl0 + k.Cl_b * be + k.Cl_p * phat + k.Cl_r * rhat +
                    k.Cl_da * da + k.Cl_dr * dr + k.Cl_daa * a * da +
                    k.Cl_a2dpt * a * a * dpt;
  const double Cm = k.Cm0 + k.Cm_a * a + k.Cm_a2 * a * a + k.Cm_q * qhat +
                    k.Cm_de * de + k.Cm_b2 * be * be + k.Cm_r * rhat +
                    k.Cm_df * df + k.Cm_dpt * dpt;
  const double Cn = k.Cn0 + k.Cn_b * be + k.Cn_p * phat + k.Cn_r * rhat +
                    k.Cn_da * da + k.Cn_dr * dr + k.Cn_q * qhat +
                    k.Cn_b3 * be * be * be + k.Cn_dpt3 * dpt * dpt * dpt;
  return {CX, CY, CZ, Cl, Cm, Cn};
}

}  // namespace autoland
