#pragma once
#include "autoland/naca4414_stall_table.hpp"

// =============================================================================
// Viscous-stall plant aero for the longitudinal model.
//
// The base aero (data/AHAB_combined.stab) is inviscid VSPAERO and cannot stall;
// its post-stall values are meaningless. So the plant BLENDS, handing the lift
// off from VSPAERO (attached only) to an absolute Viterna flat-plate post-stall
// curve as the wing stalls:
//
//   C_plant(a) = (1 - w(a)) * C_vspaero(a)  +  w(a) * C_post(a)
//
// stallLookup() returns the blend weight w in [0,1] and the absolute post-stall
// coefficients C_post = {CLpost, CDpost, CMpost} (and their d/dalpha slopes) from
// the precomputed table. LonDrift forms the blend. Mirrors clfLookup() (binary
// search + linear-in-alpha interpolation), returning VALUE AND SLOPE so
// makeAeroLocal can freeze a local affine model and keep the autodiff/Taylor Lie
// derivatives exact.
//
// Outside the grid: held with ZERO slope -- below the table w = 0 (attached flow,
// pure VSPAERO, no spurious gradient), above it the last (deep-stall) values.
// alpha in radians.
// =============================================================================
namespace autoland {

struct StallBlend {
  double w{0}, CLpost{0}, CDpost{0}, CMpost{0};
  double w_da{0}, CLpost_da{0}, CDpost_da{0}, CMpost_da{0};  // d/dalpha [1/rad]
};

inline StallBlend stallLookup(double alpha) {
  const int n = kStallN;
  const double* A = kStallAlpha;
  StallBlend s;
  if (alpha <= A[0]) {            // below grid: attached, w = 0, slopes 0
    s.w = kStallW[0]; s.CLpost = kStallCLpost[0];
    s.CDpost = kStallCDpost[0]; s.CMpost = kStallCMpost[0];
    return s;
  }
  if (alpha >= A[n - 1]) {        // above grid: deep-stall plateau, slopes 0
    s.w = kStallW[n - 1]; s.CLpost = kStallCLpost[n - 1];
    s.CDpost = kStallCDpost[n - 1]; s.CMpost = kStallCMpost[n - 1];
    return s;
  }
  int lo = 0, hi = n - 1;
  while (hi - lo > 1) {
    const int mid = (lo + hi) / 2;
    if (A[mid] <= alpha) lo = mid; else hi = mid;
  }
  const double dA = A[hi] - A[lo];
  auto interp = [&](const double* V, double& val, double& slope) {
    slope = (V[hi] - V[lo]) / dA;
    val = V[lo] + slope * (alpha - A[lo]);
  };
  interp(kStallW, s.w, s.w_da);
  interp(kStallCLpost, s.CLpost, s.CLpost_da);
  interp(kStallCDpost, s.CDpost, s.CDpost_da);
  interp(kStallCMpost, s.CMpost, s.CMpost_da);
  return s;
}

}  // namespace autoland
