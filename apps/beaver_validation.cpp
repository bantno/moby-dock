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
// --- `--xdot in.csv out.csv`: evaluate the plant xdot at arbitrary states. --
// Input rows: h_ref,n_rpm,flap, u,v,w,p,q,r,phi,theta,psi,h,y, de,da,dr,dT
// (SI/rad). Output rows: the 11 xdot components at full double precision.
// Consumed by scripts/validate_beaver_sixdof.py, which compares every row
// against the independent Python implementation over the whole flight
// envelope -- large body rates, bank, sideslip, flaps -- the regime the
// wings-level oracles cannot exercise (quadratic gyroscopic/Coriolis terms).
int dumpXdot(const std::string& in_csv, const std::string& out_csv) {
  std::ifstream in(in_csv);
  if (!in) {
    std::cerr << "cannot open " << in_csv << "\n";
    return EXIT_FAILURE;
  }
  std::ofstream outf(out_csv);
  outf << std::setprecision(17);
  std::string line;
  std::getline(in, line);  // header
  int n = 0;
  while (std::getline(in, line)) {
    std::vector<double> v;
    std::size_t pos = 0;
    while (pos <= line.size()) {
      const std::size_t nxt = line.find(',', pos);
      v.push_back(std::stod(line.substr(pos, nxt - pos)));
      if (nxt == std::string::npos) break;
      pos = nxt + 1;
    }
    if (v.size() != 18) {
      std::cerr << "bad row (" << v.size() << " fields)\n";
      return EXIT_FAILURE;
    }
    BeaverPlantConfig pc;
    pc.h_ref = v[0];
    pc.n_rpm = v[1];
    pc.flap = v[2];
    BeaverDynamics dyn(pc);
    StateVec x;
    for (int i = 0; i < NX; ++i) x[i] = v[3 + i];
    CtrlVec u;
    for (int i = 0; i < NU; ++i) u[i] = v[3 + NX + i];
    const StateVec xd = dyn.xdot(x, u);
    for (int i = 0; i < NX; ++i) outf << xd[i] << (i + 1 < NX ? ',' : '\n');
    ++n;
  }
  std::cout << "evaluated xdot at " << n << " states -> " << out_csv << "\n";
  return EXIT_SUCCESS;
}

// --- `--doublet out.csv`: open-loop 6-DOF maneuver time response. ----------
// From the 40 m/s level trim: elevator doublet (+/-3 deg, t=1..3 s), aileron
// pulse (+5 deg, t=5..5.7 s), rudder pulse (-5 deg, t=8..8.7 s), RK4 at
// dt=0.01 for 14 s. The trace is re-integrated by the independent Python
// implementation from the SAME initial state/schedule and overlaid -- a
// time-domain cross-check of the full coupled dynamics at finite rates.
int dumpDoublet(const std::string& out_csv) {
  constexpr double d = M_PI / 180.0;
  BeaverPlantConfig pc;  // sea level, n=1800, flaps up
  BeaverDynamics dyn(pc);
  const TrimResult t0 = beaverTrim(dyn, 40.0, 0.0);
  if (!t0.converged) {
    std::cerr << "doublet: trim failed\n";
    return EXIT_FAILURE;
  }
  std::ofstream outf(out_csv);
  outf << std::setprecision(17);
  outf << "t,u,v,w,p,q,r,phi,theta,psi,h,y,de,da,dr,dT\n";
  StateVec x = t0.x;
  x[H] = 200.0;
  const double dt = 0.01;
  for (int k = 0; k <= 1400; ++k) {
    const double t = k * dt;
    CtrlVec u = t0.u;
    if (t >= 1.0 && t < 2.0) u[DE] += 3.0 * d;
    else if (t >= 2.0 && t < 3.0) u[DE] -= 3.0 * d;
    if (t >= 5.0 && t < 5.7) u[DA] += 5.0 * d;
    if (t >= 8.0 && t < 8.7) u[DR] -= 5.0 * d;
    outf << t;
    for (int i = 0; i < NX; ++i) outf << ',' << x[i];
    for (int i = 0; i < NU; ++i) outf << ',' << u[i];
    outf << '\n';
    x = rk4Step([&dyn](const StateVec& xs, const CtrlVec& us)
                    { return dyn.xdot(xs, us); },
                x, u, dt);
  }
  std::cout << "doublet response (14 s open loop) -> " << out_csv << "\n";
  return EXIT_SUCCESS;
}
}  // namespace

int main(int argc, char** argv) {
  if (argc > 1 && std::string(argv[1]) == "--xdot") {
    if (argc < 4) {
      std::cerr << "usage: beaver_validation --xdot in.csv out.csv\n";
      return EXIT_FAILURE;
    }
    return dumpXdot(argv[2], argv[3]);
  }
  if (argc > 1 && std::string(argv[1]) == "--doublet") {
    return dumpDoublet(argc > 2 ? argv[2] : "results/beaver_doublet_cpp.csv");
  }

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
