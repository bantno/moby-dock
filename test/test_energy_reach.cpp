#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <fstream>
#include <string>

#include "autoland/beaver_dynamics.hpp"
#include "autoland/energy_reach_cbf.hpp"
#include "autoland/sixdof_sim.hpp"

using namespace autoland;
using Catch::Approx;

namespace {
constexpr double kDeg = M_PI / 180.0;

ReachParams testParams() {
  BeaverDynamics dyn;  // sea level, flaps up
  const ReachPolarFit f = fitBeaverReachPolar(dyn, 35.0, 0.0);
  ReachParams p;
  p.rho = dyn.rho();
  p.g = dyn.g();
  p.CLa = f.CLa;
  p.CD0 = f.CD0;
  p.CDa2 = f.CDa2;
  p.alpha0L = f.alpha0L;
  p.CYb = dyn.config().aero.Cy_b;
  p.Nt = 1500.0;
  p.Et = 100.0;
  p.V_min = 30.0;
  return p;
}

// Independent closed-form gradient of the barrier (chain rules written by
// hand -- a SEPARATE derivation path from the Taylor/autodiff jet engine).
// Valid away from the sinc-series seam (|phi_s| > 0.05) and the d_near switch.
struct HandGrad {
  double h;
  std::array<double, NXER> grad;
};
HandGrad handBarrierGrad(const ReachParams& p, const ReachStateArr& X) {
  const double V = X[ERV], chi = X[ERCHI], alpha = X[ERAL];
  // Smooth-max alpha floor.
  const double s = (alpha - p.alpha0L) - p.alpha_lo;
  const double rt = std::sqrt(s * s + p.alpha_smooth * p.alpha_smooth);
  const double ats = p.alpha_lo + 0.5 * (s + rt);
  const double dats = 0.5 * (1.0 + s / rt);
  // Arc distance.
  const double dN = p.Nt - X[ERN], dE = p.Et - X[ERE];
  const double R = std::sqrt(dN * dN + dE * dE + p.R_eps * p.R_eps);
  const double dpar = dN * std::cos(chi) + dE * std::sin(chi);
  const double dperp = dN * std::sin(chi) - dE * std::cos(chi);
  const double n2 = dpar * dpar + dperp * dperp;
  const double phi = std::atan2(dperp, dpar);
  const double th = std::tanh(phi / p.phi_sat);
  const double phis = p.phi_sat * th;
  const double dphis = 1.0 - th * th;  // d(phis)/d(phi)
  const double sp = std::sin(phis), cp = std::cos(phis);
  const double q = phis / sp;
  const double dq = (sp - phis * cp) / (sp * sp);  // q'(phis)
  const double d = R * q;
  // phi = atan2(dperp, dpar) with dN, dE = target - vehicle, so wrt N/E:
  //   d(dpar)/dN = -cos(chi)   d(dperp)/dN = -sin(chi)
  //   d(dpar)/dE = -sin(chi)   d(dperp)/dE = +cos(chi)
  const double ddpar_dN = -std::cos(chi), ddperp_dN = -std::sin(chi);
  const double ddpar_dE = -std::sin(chi), ddperp_dE = std::cos(chi);
  const double dphi_dN = (dpar * ddperp_dN - dperp * ddpar_dN) / n2;
  const double dphi_dE = (dpar * ddperp_dE - dperp * ddpar_dE) / n2;
  const double dR_dN = -dN / R, dR_dE = -dE / R;
  const double dd_dchi = R * dq * dphis;  // dphi/dchi = 1 exactly
  const double dd_dN = dR_dN * q + R * dq * dphis * dphi_dN;
  const double dd_dE = dR_dE * q + R * dq * dphis * dphi_dE;
  // Glide-ratio term.
  const double Psi = (p.CD0 + p.CDa2 * ats * ats) / (p.CLa * ats);
  const double dPsi = (p.CDa2 * ats * ats - p.CD0) / (p.CLa * ats * ats);

  HandGrad out;
  out.h = 0.5 * V * V + p.g * X[ERZ] - 0.5 * p.V_min * p.V_min - p.g * d * Psi;
  out.grad = {V,
              0.0,
              -p.g * Psi * dd_dchi,
              -p.g * Psi * dd_dN,
              -p.g * Psi * dd_dE,
              p.g,
              -p.g * d * dPsi * dats};
  return out;
}
}  // namespace

TEST_CASE("Beaver idle-glide polar fit is physical", "[energy_reach]") {
  BeaverDynamics dyn;
  const ReachPolarFit f = fitBeaverReachPolar(dyn, 35.0, 0.0);
  INFO("CLa=" << f.CLa << " CD0=" << f.CD0 << " CDa2=" << f.CDa2
              << " alpha0L=" << f.alpha0L / kDeg << " deg  LDmax=" << f.LD_max
              << " at at*=" << f.alpha_star / kDeg << " deg  rmsCL=" << f.rms_CL
              << " rmsCD=" << f.rms_CD);
  CHECK(f.CLa > 4.0);
  CHECK(f.CLa < 7.0);
  CHECK(f.CD0 > 0.0);
  CHECK(f.CDa2 > 0.0);
  // Windmilling-prop glide ratio: well below a clean airframe, above a brick.
  CHECK(f.LD_max > 4.0);
  CHECK(f.LD_max < 14.0);
  // The quadratic polar should track the cubic polynomials tightly in-band.
  CHECK(f.rms_CL < 0.05);
  CHECK(f.rms_CD < 0.01);
  // The windmilling prop inflates CD0 so hard that the polar's best-glide
  // alpha sits ABOVE the trusted AoA band (~17-26 deg vs the 14 deg ceiling):
  // in-band, L/D is monotonically increasing and the CBF rides the ceiling
  // when energy-critical. Only sanity-bound the extrapolated optimum.
  CHECK(f.alpha0L + f.alpha_star > 5.0 * kDeg);
  CHECK(f.alpha0L + f.alpha_star < 30.0 * kDeg);
}

TEST_CASE("Reach Lie derivatives match an independent hand derivation",
          "[energy_reach]") {
  const ReachParams p = testParams();
  // Generic state away from every guard seam: phi ~ 9.3 deg, alpha well above
  // the floor, 1.3 km from the target.
  const ReachStateArr X{34.0, -0.12, 0.3, 200.0, -80.0, 180.0, 6.0 * kDeg};

  const HandGrad hg = handBarrierGrad(p, X);
  const ReachLie lie = reachLie(p, X);

  CHECK(lie.h == Approx(hg.h).epsilon(1e-10));

  // L_f h = grad(h) . f  with f from the (double) drift.
  const std::array<double, NXER> f = ReachDrift(p)(X);
  double Lfh = 0.0;
  for (int i = 0; i < NXER; ++i) Lfh += hg.grad[i] * f[i];
  CHECK(lie.Lfh == Approx(Lfh).epsilon(1e-8));

  // Control row: A_ua = dh/dalpha, A_be = dh/dchi * G(chi, beta).
  const Eigen::Matrix<double, NXER, NUER> G = reachGMatrix(p, X);
  CHECK(lie.A_ua == Approx(hg.grad[ERAL]).epsilon(1e-8));
  CHECK(lie.A_be == Approx(hg.grad[ERCHI] * G(ERCHI, ERBE)).epsilon(1e-8));

  // Signs at this state: alpha (6 deg from zero-lift ~ below best glide) up
  // => more margin; heading error phi > 0 with CYb < 0 => positive beta turns
  // chi DOWN toward the bearing => more margin.
  CHECK(lie.A_ua > 0.0);
  CHECK(lie.A_be > 0.0);
}

TEST_CASE("Filter passes the nominal through when the margin is large",
          "[energy_reach]") {
  EnergyReachCfg cfg;
  cfg.prm = testParams();
  const EnergyReachFilter filt(cfg);
  // 400 m up, 1.3 km out: h_E >> 0.
  const ReachStateArr X{35.0, -0.06, 0.14, 200.0, -80.0, 400.0, 5.0 * kDeg};
  EnergyReachDiag d;
  const Eigen::Vector2d u = filt.filter(0.01, 0.005, X, &d);
  REQUIRE_FALSE(d.best_effort);
  CHECK(d.h > 0.0);
  CHECK_FALSE(d.active);
  CHECK(u[0] == Approx(0.01).margin(1e-6));
  CHECK(u[1] == Approx(0.005).margin(1e-6));
}

TEST_CASE("Binding filter raises alpha toward best glide and steers at the "
          "target", "[energy_reach]") {
  EnergyReachCfg cfg;
  cfg.prm = testParams();
  const EnergyReachFilter filt(cfg);
  // Low on energy (h_E < 0 at this state), heading 20 deg right of the
  // bearing, alpha below best-glide.
  ReachStateArr X{34.0, -0.10, 0.5, 200.0, -80.0, 120.0, 3.0 * kDeg};
  EnergyReachDiag d;
  const Eigen::Vector2d u = filt.filter(-0.02, 0.0, X, &d);
  REQUIRE_FALSE(d.best_effort);
  INFO("h=" << d.h << " Lfh=" << d.Lfh << " A_ua=" << d.A_ua
            << " A_be=" << d.A_be << " ua=" << u[0] << " be=" << u[1]);
  CHECK(d.h < 0.0);
  CHECK(d.active);
  CHECK(u[0] > 0.0);       // pull alpha UP toward best glide, not the dive
  CHECK(u[1] > 0.0);       // skid LEFT (CYb < 0) toward the target bearing
  CHECK(u[0] <= cfg.ua_max + 1e-9);
  CHECK(u[1] <= cfg.beta_max + 1e-9);
}

TEST_CASE("Closed-loop 3-DOF invariance under an adversarial nominal",
          "[energy_reach]") {
  EnergyReachCfg cfg;
  cfg.prm = testParams();
  cfg.prm.Nt = 1600.0;
  cfg.prm.Et = 0.0;
  const EnergyReachFilter filt(cfg);
  const ReachDrift drift(cfg.prm);

  // Start with modest positive margin (alpha near best in-band glide),
  // heading 15 deg off the bearing. The 30 s horizon covers the class-K decay
  // and a long boundary-riding stretch but stops BEFORE the endgame where the
  // adversarial nominal has spent every reserve: a plain (myopic) CBF-QP with
  // BOUNDED inputs only guarantees invariance while the QP stays feasible,
  // and this adversary eventually drives the geometry (target abeam, low,
  // weak skid authority) somewhere no admissible input can hold.
  ReachStateArr X{35.0, -0.09, 15.0 * kDeg, 0.0, 0.0, 260.0, 9.0 * kDeg};
  {
    EnergyReachDiag d0;
    filt.filter(0.0, 0.0, X, &d0);
    REQUIRE(d0.h > 0.0);
  }

  const double dt = 0.02;
  double min_h = 1e30;
  double min_envelope_slip = 1e30;  // h(t) - h0 exp(-ck t), most-negative
  double h0 = 0.0;
  bool first = true;
  for (int k = 0; k < 1500 && X[ERZ] > 1.0; ++k) {
    EnergyReachDiag d;
    // Adversarial nominal: full-rate dive toward min drag, no steering.
    const Eigen::Vector2d u = filt.filter(-cfg.ua_max, 0.0, X, &d);
    // (best-effort steps are allowed near the boundary; the envelope check
    // below is the invariance assertion)
    if (first) { h0 = d.h; first = false; }
    min_h = std::min(min_h, d.h);
    min_envelope_slip =
        std::min(min_envelope_slip, d.h - h0 * std::exp(-cfg.ck * k * dt));
    // RK4 with the filtered control held (ZOH, matching the sim structure).
    auto xdot = [&](const ReachStateArr& x) {
      std::array<double, NXER> dx = drift(x);
      const Eigen::Matrix<double, NXER, NUER> G = reachGMatrix(cfg.prm, x);
      for (int i = 0; i < NXER; ++i) dx[i] += G(i, ERUA) * u[0] + G(i, ERBE) * u[1];
      return dx;
    };
    const auto k1 = xdot(X);
    ReachStateArr t2, t3, t4;
    for (int i = 0; i < NXER; ++i) t2[i] = X[i] + 0.5 * dt * k1[i];
    const auto k2 = xdot(t2);
    for (int i = 0; i < NXER; ++i) t3[i] = X[i] + 0.5 * dt * k2[i];
    const auto k3 = xdot(t3);
    for (int i = 0; i < NXER; ++i) t4[i] = X[i] + dt * k3[i];
    const auto k4 = xdot(t4);
    for (int i = 0; i < NXER; ++i)
      X[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
  }
  INFO("h_E(0)=" << h0 << "  min h_E=" << min_h
                 << "  envelope slip=" << min_envelope_slip);
  // Forward invariance up to the discrete-time (ZOH) slip: the barrier must
  // never dip more than a fraction of a percent of its starting value...
  CHECK(min_h > -0.005 * std::abs(h0));
  // ...and must respect the class-K decay envelope h(t) >= h0 exp(-ck t)
  // throughout (the RD-1 constraint integrated in continuous time), again up
  // to the ZOH slip.
  CHECK(min_envelope_slip > -0.02 * std::abs(h0));
}

TEST_CASE("Positive rudder raises Beaver sideslip (beta-loop sign)",
          "[energy_reach]") {
  // dvdot/ddr > 0 (Cy_dr side force) and the yaw path reinforce: the exact
  // linearization's B column must push v (hence beta) up for +dr.
  BeaverDynamics dyn;
  const TrimResult tr = beaverTrim(dyn, 35.0, -4.0 * kDeg);
  REQUIRE(tr.converged);
  Mat A, B;
  dyn.linearize(tr.x, tr.u, A, B);
  CHECK(B(V, DR) > 0.0);
}

TEST_CASE("6-DOF idle-glide smoke run with the energy-reach layer",
          "[energy_reach][slow]") {
  const std::string path = "test_reach_scenario.yaml";
  {
    std::ofstream f(path);
    f << "plant: beaver\n"
         "V_app: 40.0\n"
         "dt: 0.01\n"
         "t_max: 25.0\n"
         "initial:\n"
         "  h0: 300.0\n"
         "energy_cbf:\n"
         "  enabled: true\n"
         "  throttle: idle\n"
         "  V_min: 30.0\n";
  }
  SixDofSim sim("", "", path);
  REQUIRE(sim.stats().trim_converged);
  // Idle-glide gamma replaced the default approach slope (steeper).
  CHECK(sim.scenario().gamma_app < -4.0 * kDeg);
  CHECK(sim.scenario().reach.enabled);
  const ReachPolarFit& pf = sim.reachPolar();
  CHECK(pf.LD_max > 4.0);
  sim.run("test_reach_smoke.csv");
  const SixDofRunStats& st = sim.stats();
  CHECK(st.steps > 1000);
  CHECK(st.cbf_best_effort_steps == 0);
  CHECK(std::isfinite(st.min_hE));
  // The aim-point target starts exactly on the glide, so the barrier begins
  // near h_E ~ 1/2 (V^2 - V_min^2) > 0 and must stay above a small dip.
  CHECK(st.min_hE > -50.0);
  CHECK(st.min_V > 25.0);
  CHECK(st.max_abs_phi < 10.0 * kDeg);  // wings-level lateral mode
}
