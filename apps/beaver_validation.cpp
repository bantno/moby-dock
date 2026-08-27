#include <Eigen/Eigenvalues>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "autoland/beaver_dynamics.hpp"

// =============================================================================
// Beaver flight-dynamics validation harness (see documentation/
// beaver_validation.md):
//
//  1. FDC ACTRIM CHECK CASE (external oracle): the FDC 1.2 manual, figs.
//     10.18/10.19, prints a complete trimmed state/input/xdot for the Beaver at
//     V=35 m/s, h=2000 ft, n=1800 RPM, pz=20 "Hg. We evaluate OUR 6-DOF xdot at
//     exactly that printed (x, u) and compare against the printed xdot -- an
//     end-to-end check of aero + engine + atmosphere + EOM assembly.
//  2. TRIM CURVE: level-flight trim sweep V=30..60 (cf. FDC fig. 10.13
//     trimmed-flight elevator deflection curve).
//  3. LINEARIZATION: exact autodiff A/B at reference conditions; longitudinal /
//     lateral-directional eigenvalues; full A,B dumped to CSV for the
//     independent Python cross-implementation (scripts/validate_beaver_modes.py).
//
// Usage: beaver_validation [out_dir]     (default out_dir = "results/beaver")
// =============================================================================
using namespace autoland;

namespace {
constexpr double kDeg = M_PI / 180.0;

// FDC state order [V alpha beta p q r psi theta phi xe ye H] -> printed values
// (manual fig. 10.19).
struct FdcCheckCase {
  double V = 35.0, alpha = 2.1131e-1, beta = -2.0667e-2;
  double theta = 1.9190e-1;
  double de = -9.3083e-2, da = 9.6242e-3, dr = -4.9506e-2, df = 0.0;
  double n_rpm = 1800.0, pz = 20.0;
  // ACTRIM was PROMPTED with 2000 ft, but the printed trimmed STATE carries
  // H = 0 and the printed accelerations are only reproduced with the
  // sea-level atmosphere (rho = 1.225): at rho(2000 ft) the residuals are
  // ~0.1 m/s^2, at rho(0) they match the printout to ~3e-6. So the check
  // case is evaluated at h = 0 (the altitude prompt seeds the simulation IC,
  // not the trim atmosphere).
  double h = 0.0;
  // Printed xdot (same FDC order).
  double Vdot = -1.8871e-4, alphadot = -1.2348e-5, betadot = 4.6356e-4;
  double pdot = -2.5027e-5, qdot = -2.0660e-5, rdot = -5.2604e-5;
  double xedot = 3.4986e1, yedot = -7.2330e-1, Hdot = -6.7909e-1;
};

// Body-axis (u,v,w) state from the FDC (V, alpha, beta) coordinates.
StateVec fdcStateToBody(const FdcCheckCase& c) {
  StateVec x = StateVec::Zero();
  x[U] = c.V * std::cos(c.alpha) * std::cos(c.beta);
  x[V] = c.V * std::sin(c.beta);
  x[W] = c.V * std::sin(c.alpha) * std::cos(c.beta);
  x[THETA] = c.theta;
  return x;
}

// Convert body-axis (udot,vdot,wdot) to FDC (Vdot, alphadot, betadot).
void bodyRatesToFdc(const StateVec& x, const StateVec& xd, double& Vdot,
                    double& alphadot, double& betadot) {
  const double u = x[U], v = x[V], w = x[W];
  const double Vt = std::sqrt(u * u + v * v + w * w);
  Vdot = (u * xd[U] + v * xd[V] + w * xd[W]) / Vt;
  alphadot = (u * xd[W] - w * xd[U]) / (u * u + w * w);
  betadot = (xd[V] * Vt - v * Vdot) / (Vt * std::sqrt(u * u + w * w));
}

void printRow(const std::string& name, double mine, double fdc) {
  std::cout << "  " << std::left << std::setw(10) << name << std::right
            << std::setw(14) << mine << std::setw(14) << fdc << std::setw(14)
            << mine - fdc << "\n";
}

// Eigenvalues of the 4x4 submatrix of A on the given state indices.
Eigen::Vector4cd subEigs(const Mat& A, const std::array<int, 4>& idx) {
  Eigen::Matrix4d S;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) S(i, j) = A(idx[i], idx[j]);
  return Eigen::EigenSolver<Eigen::Matrix4d>(S).eigenvalues();
}

void printModes(const std::string& label, const Mat& A) {
  const Eigen::Vector4cd lon = subEigs(A, {U, W, Q, THETA});
  const Eigen::Vector4cd lat = subEigs(A, {V, P, R, PHI});
  std::cout << label << "\n  longitudinal [u w q theta]:\n";
  for (int i = 0; i < 4; ++i)
    std::cout << "    " << lon[i].real() << (lon[i].imag() >= 0 ? " + " : " - ")
              << std::abs(lon[i].imag()) << "i\n";
  std::cout << "  lateral-directional [v p r phi]:\n";
  for (int i = 0; i < 4; ++i)
    std::cout << "    " << lat[i].real() << (lat[i].imag() >= 0 ? " + " : " - ")
              << std::abs(lat[i].imag()) << "i\n";
}

void dumpMat(const std::string& path, const Mat& M) {
  std::ofstream f(path);
  f << std::setprecision(17);
  for (int i = 0; i < M.rows(); ++i) {
    for (int j = 0; j < M.cols(); ++j)
      f << M(i, j) << (j + 1 < M.cols() ? ',' : '\n');
  }
}
}  // namespace

int main(int argc, char** argv) {
  const std::string out = argc > 1 ? argv[1] : "results/beaver";
  std::cout << std::setprecision(6);

  // ---------------------------------------------------------------- check case
  FdcCheckCase c;
  BeaverPlantConfig pc;
  pc.n_rpm = c.n_rpm;
  pc.h_ref = c.h;
  BeaverDynamics dyn(pc);

  // Recover the throttle that reproduces pz = 20 "Hg exactly.
  const double dT = (c.pz - pc.pz_idle) / (pc.pz_max - pc.pz_idle);
  CtrlVec u = CtrlVec::Zero();
  u[DE] = c.de;
  u[DA] = c.da;
  u[DR] = c.dr;
  u[DT] = dT;

  const StateVec x = fdcStateToBody(c);
  const StateVec xd = dyn.xdot(x, u);
  double Vdot, alphadot, betadot;
  bodyRatesToFdc(x, xd, Vdot, alphadot, betadot);

  std::cout << "=== 1. FDC ACTRIM check case (V=35, n=1800, pz=20, sea-level "
               "atmosphere) ===\n";
  std::cout << "  atmosphere: rho=" << dyn.rho() << " kg/m^3  g=" << dyn.g()
            << " m/s^2  P=" << dyn.power(dT) << " kW  dpt="
            << dyn.dpt(dT, c.V) << "\n";
  std::cout << "  " << std::left << std::setw(10) << "row" << std::right
            << std::setw(14) << "mine" << std::setw(14) << "FDC printed"
            << std::setw(14) << "diff" << "\n";
  printRow("Vdot", Vdot, c.Vdot);
  printRow("alphadot", alphadot, c.alphadot);
  printRow("betadot", betadot, c.betadot);
  printRow("pdot", xd[P], c.pdot);
  printRow("qdot", xd[Q], c.qdot);
  printRow("rdot", xd[R], c.rdot);
  // North ground speed (psi = phi = 0): xedot = u cos(theta) + w sin(theta).
  printRow("xedot",
           x[U] * std::cos(x[THETA]) + x[W] * std::sin(x[THETA]), c.xedot);
  printRow("yedot", xd[Y], c.yedot);
  printRow("Hdot", xd[H], c.Hdot);

  // ------------------------------------------------------------- trim recovery
  // Reproduce the check case with OUR trim solver: same gamma as the FDC
  // descent (gamma = asin(Hdot/V)), same atmosphere; the solved controls and
  // incidence angles should land on the FDC values.
  const double gamma_fdc = std::asin(c.Hdot / c.V);
  const TrimResult tr = beaverTrim(dyn, c.V, gamma_fdc);
  std::cout << "\n=== 2. Trim recovery at gamma=" << gamma_fdc / kDeg
            << " deg (FDC values in parens) ===\n"
            << "  converged=" << tr.converged << " res=" << tr.residual
            << " iters=" << tr.iterations << "\n";
  std::cout << "  alpha = " << tr.alpha / kDeg << " deg  (" << c.alpha / kDeg
            << ")\n";
  std::cout << "  beta  = " << std::asin(tr.x[V] / c.V) / kDeg << " deg  ("
            << c.beta / kDeg << ")\n";
  std::cout << "  theta = " << tr.theta / kDeg << " deg  (" << c.theta / kDeg
            << ")\n";
  std::cout << "  de    = " << tr.u[DE] / kDeg << " deg  (" << c.de / kDeg
            << ")\n";
  std::cout << "  da    = " << tr.u[DA] / kDeg << " deg  (" << c.da / kDeg
            << ")\n";
  std::cout << "  dr    = " << tr.u[DR] / kDeg << " deg  (" << c.dr / kDeg
            << ")\n";
  std::cout << "  pz    = " << dyn.pzFromThrottle(tr.u[DT]) << " \"Hg  ("
            << c.pz << ")\n";

  // ---------------------------------------------------------------- trim sweep
  // Level-flight elevator curve (cf. FDC fig. 10.13) at n=1800 and n=2300.
  std::ofstream sweep(out + "_trim_sweep.csv");
  sweep << "n_rpm,V,alpha_deg,beta_deg,theta_deg,de_deg,da_deg,dr_deg,dT,pz,"
           "P_kW,converged\n";
  std::cout << "\n=== 3. Level-flight trim sweep (sea level) -> " << out
            << "_trim_sweep.csv ===\n";
  for (double n : {1800.0, 2300.0}) {
    BeaverPlantConfig psw = pc;
    psw.n_rpm = n;
    BeaverDynamics dsw(psw);
    for (double Vs = 30.0; Vs <= 60.0 + 1e-9; Vs += 1.0) {
      const TrimResult t = beaverTrim(dsw, Vs, 0.0);
      sweep << n << ',' << Vs << ',' << t.alpha / kDeg << ','
            << std::asin(t.x[V] / Vs) / kDeg << ',' << t.theta / kDeg << ','
            << t.u[DE] / kDeg << ',' << t.u[DA] / kDeg << ','
            << t.u[DR] / kDeg << ',' << t.u[DT] << ','
            << dsw.pzFromThrottle(t.u[DT]) << ',' << dsw.power(t.u[DT]) << ','
            << t.converged << '\n';
      if (!t.converged)
        std::cout << "  WARNING: trim failed at n=" << n << " V=" << Vs
                  << " (res=" << t.residual << ")\n";
    }
  }
  sweep.close();

  // -------------------------------------------------------------- linearization
  std::cout << "\n=== 4. Linearized modes (exact autodiff Jacobians) ===\n";
  struct Cond {
    const char* name;
    double V, gamma, h, n;
  };
  const std::vector<Cond> conds = {
      {"cruise_45", 45.0, 0.0, 1828.8, 1800.0},
      {"check_35", 35.0, gamma_fdc, 0.0, 1800.0},
      {"approach_35", 35.0, -3.5 * kDeg, 0.0, 1800.0},
  };
  for (const Cond& cd : conds) {
    BeaverPlantConfig pl = pc;
    pl.h_ref = cd.h;
    pl.n_rpm = cd.n;
    BeaverDynamics dl(pl);
    const TrimResult t = beaverTrim(dl, cd.V, cd.gamma);
    if (!t.converged) {
      std::cout << cd.name << ": TRIM FAILED (res=" << t.residual << ")\n";
      continue;
    }
    Mat A, B;
    dl.linearize(t.x, t.u, A, B);
    dumpMat(out + "_A_" + cd.name + ".csv", A);
    dumpMat(out + "_B_" + cd.name + ".csv", B);
    // Trim point dump for the Python cross-check (same file set).
    std::ofstream tp(out + "_trim_" + cd.name + ".csv");
    tp << std::setprecision(17)
       << "V,gamma,h,n_rpm,alpha,beta,theta,de,da,dr,dT,pz,rho,g,flap\n"
       << cd.V << ',' << cd.gamma << ',' << cd.h << ',' << cd.n << ','
       << t.alpha << ',' << std::asin(t.x[V] / cd.V) << ',' << t.theta << ','
       << t.u[DE] << ',' << t.u[DA] << ',' << t.u[DR] << ',' << t.u[DT] << ','
       << dl.pzFromThrottle(t.u[DT]) << ',' << dl.rho() << ',' << dl.g() << ','
       << pl.flap << '\n';
    printModes(std::string(cd.name) + "  (V=" + std::to_string(cd.V) + ")", A);
  }

  std::cout << "\ndone.\n";
  return EXIT_SUCCESS;
}
