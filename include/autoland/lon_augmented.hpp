#pragma once
#include <array>
#include <cmath>
#include <Eigen/Dense>
#include "autoland/aero_table.hpp"
#include "autoland/config.hpp"
#include "autoland/lie_taylor.hpp"
#include "autoland/mixing.hpp"

// =============================================================================
// Augmented longitudinal system for the water-landing CBF-QP.
//   State   X = [h, V, gamma, theta, q, T, Tdot]   (NXA = 7)
//   Control U = [delta_e, Tddot]                    (NUA = 2)
// Thrust is augmented with two integrators (T, Tdot states; Tddot the control)
// to align the relative degrees of elevator and thrust. See
// documentation/water_landing_cbf_math.md section 1.
//
// AERO FRAME: the EOM are written in wind/path axes, so they need wind-axis lift
// L and drag D. VSPAero's primary data is body-axis (CFx,CFz,...); we ROTATE it
// into the wind frame analytically, using the VSPAERO-documented transform
//   C_L = -CFx*sin a + CFz*cos a,   C_D = CFx*cos a + CFz*sin a   (beta = 0)
// applied to the FILE-frame coefficients (before dynamics' kFrameSign flip).
// Verified against data/example.stab. Pitch moment stays body-axis (CMm == CMy).
//
// LOCAL AFFINE AERO: within an aero-table cell the model is affine in
// (alpha, mach, qhat, de) with constant tabulated slopes. AeroLocal freezes
// those slopes at the evaluation point, so the templated EOM are smooth and the
// autodiff/Taylor Lie derivatives are EXACT (in-cell curvature is zero).
// =============================================================================
namespace autoland {

constexpr int NXA = 7;
enum LonState { LH = 0, LV = 1, LGAM = 2, LTH = 3, LQ = 4, LT = 5, LTDOT = 6 };
constexpr int NUA = 2;
enum LonCtrl { LDE = 0, LTDDOT = 1 };

using LonStateVec = Eigen::Matrix<double, NXA, 1>;
using LonCtrlVec = Eigen::Matrix<double, NUA, 1>;
using LonGMat = Eigen::Matrix<double, NXA, NUA>;

// Frozen local affine aero model + geometry/config scalars. Coefficients are
//   C = off + dAlpha*alpha + dMach*mach + dQ*qhat     (file frame; de=0 drift)
// with mach = V/a_sound, qhat = q*cref/(2V) evaluated live in the EOM.
struct AeroLocal {
  double Sref{0}, cref{0}, a_sound{0};
  double mass{0}, Iyy{0}, rho{0}, g{0}, zcp{0}, parasite_CD0{0};
  double off_CFx{0}, dAlpha_CFx{0}, dMach_CFx{0}, dQ_CFx{0};
  double off_CFz{0}, dAlpha_CFz{0}, dMach_CFz{0}, dQ_CFz{0};
  double off_CMy{0}, dAlpha_CMy{0}, dMach_CMy{0}, dQ_CMy{0};
  double dDe_CMy{0};  // d(C_my)/d(virtual delta_e), through the mixing matrix
};

// Build AeroLocal at flight condition (V, alpha) [SI, rad]. Looks up the aero
// table once (beta = 0), freezes slopes, folds angle/Mach reference offsets in.
// (q is not needed: the table lookup is over (alpha,beta,mach); the qhat rate
// term is reconstructed live in the EOM.)
AeroLocal makeAeroLocal(const AeroTable& table, const Mixing& mixing,
                        const AircraftConfig& cfg, double V, double alpha);

// Templated drift f(X) with U = 0. Element type T is double (plant/eval) or a
// Taylor series (Lie derivatives). All ops go through ADL so Taylor overloads
// are picked up automatically.
struct LonDrift {
  const AeroLocal& a;
  explicit LonDrift(const AeroLocal& aero) : a(aero) {}

  template <class T>
  std::array<T, NXA> operator()(const std::array<T, NXA>& X) const {
    using std::sin;
    using std::cos;
    const T V = X[LV], gam = X[LGAM], th = X[LTH], q = X[LQ];
    const T Tr = X[LT], Td = X[LTDOT];
    const T alpha = th - gam;
    const T mach = V * (1.0 / a.a_sound);
    const T qhat = (q * (a.cref / 2.0)) / V;
    const T qbar = (0.5 * a.rho) * V * V;

    const T CFx = a.off_CFx + a.dAlpha_CFx * alpha + a.dMach_CFx * mach + a.dQ_CFx * qhat;
    const T CFz = a.off_CFz + a.dAlpha_CFz * alpha + a.dMach_CFz * mach + a.dQ_CFz * qhat;
    const T CMy = a.off_CMy + a.dAlpha_CMy * alpha + a.dMach_CMy * mach + a.dQ_CMy * qhat;

    const T sa = sin(alpha), ca = cos(alpha);
    const T CD = CFx * ca + CFz * sa + a.parasite_CD0;  // wind-axis, beta=0
    const T CL = (-1.0) * CFx * sa + CFz * ca;
    const T Lift = qbar * a.Sref * CL;
    const T Drag = qbar * a.Sref * CD;
    const T M = qbar * a.Sref * a.cref * CMy + Tr * a.zcp;

    const T sg = sin(gam), cg = cos(gam);
    std::array<T, NXA> dX;
    dX[LH] = V * sg;
    dX[LV] = (Tr * ca - Drag) * (1.0 / a.mass) - a.g * sg;
    dX[LGAM] = (Tr * sa + Lift) / (a.mass * V) - (a.g * cg) / V;
    dX[LTH] = q;
    dX[LQ] = M * (1.0 / a.Iyy);
    dX[LT] = Td;
    dX[LTDOT] = T(0.0);
    return dX;
  }
};

// Control matrix g(X) (doubles). delta_e acts only through the pitch moment
// (LQ row); Tddot is the LTDOT-row integrator. Sparse by construction -> the
// elevator/thrust relative degrees in the doc follow.
LonGMat gMatrix(const AeroLocal& a, const LonStateVec& X);

// Convenience: convert between Eigen and std::array used by the Lie engine.
inline std::array<double, NXA> toArray(const LonStateVec& x) {
  std::array<double, NXA> a;
  for (int i = 0; i < NXA; ++i) a[i] = x[i];
  return a;
}
inline LonStateVec toVec(const std::array<double, NXA>& a) {
  LonStateVec x;
  for (int i = 0; i < NXA; ++i) x[i] = a[i];
  return x;
}

// Full augmented RHS Xdot = f(X) + g(X) U (doubles), the self-consistent plant.
inline LonStateVec lonXdot(const AeroLocal& a, const LonStateVec& X,
                           const LonCtrlVec& U) {
  std::array<double, NXA> d = LonDrift(a)(toArray(X));
  LonStateVec xd = toVec(d);
  xd += gMatrix(a, X) * U;
  return xd;
}

}  // namespace autoland
