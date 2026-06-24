#include "autoland/sim.hpp"

#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>

#include "autoland/nominal_lon_glideslope.hpp"

namespace autoland {
namespace {

constexpr double kRad2Deg = 180.0 / M_PI;

double airspeed(const StateVec& x) {
  return std::max(1e-3, std::sqrt(x[U] * x[U] + x[V] * x[V] + x[W] * x[W]));
}
double alphaOf(const StateVec& x) { return std::atan2(x[W], x[U]); }
double betaOf(const StateVec& x) {
  return std::asin(std::max(-1.0, std::min(1.0, x[V] / airspeed(x))));
}
double sinkOf(const StateVec& x) {
  const double ct = std::cos(x[THETA]), st = std::sin(x[THETA]);
  const double cp = std::cos(x[PHI]), sp = std::sin(x[PHI]);
  return -(x[U] * st - x[V] * sp * ct - x[W] * cp * ct);  // positive down
}
double northVel(const StateVec& x) {
  const double ct = std::cos(x[THETA]), st = std::sin(x[THETA]);
  const double cp = std::cos(x[PHI]), sp = std::sin(x[PHI]);
  const double cps = std::cos(x[PSI]), sps = std::sin(x[PSI]);
  return x[U] * ct * cps + x[V] * (sp * st * cps - cp * sps) +
         x[W] * (cp * st * cps + sp * sps);
}

// Central-difference state gradient of a scalar barrier h(x). Lets a barrier be
// added with only its h(x) defined (no hand-derived Jacobian), consistent with
// the project's numerical-linearization approach.
Eigen::RowVectorXd numGrad(const std::function<double(const StateVec&)>& h,
                           const StateVec& x) {
  Eigen::RowVectorXd g = Eigen::RowVectorXd::Zero(NX);
  for (int i = 0; i < NX; ++i) {
    const double step = 1e-6 * std::max(1.0, std::abs(x[i]));
    StateVec xp = x, xm = x;
    xp[i] += step;
    xm[i] -= step;
    g[i] = (h(xp) - h(xm)) / (2.0 * step);
  }
  return g;
}

}  // namespace

Sim::Sim(const std::string& stab_path, const std::string& aircraft_yaml,
         const std::string& scenario_yaml)
    : table_(AeroTable::fromFile(stab_path)),
      ac_(loadAircraftConfig(aircraft_yaml)),
      sc_(loadScenarioConfig(scenario_yaml)) {
  mixing_ = std::make_unique<Mixing>(Mixing::build(ac_, table_));
  dyn_ = std::make_unique<Dynamics>(table_, *mixing_, ac_);

  trim_ = trim(*dyn_, sc_.V_app, sc_.gamma_app);
  lin_ = linearize(*dyn_, trim_.x, trim_.u);

  // Configure the CBF safety filter over the longitudinal controls {de, dT}.
  CBFConfig ccfg;
  ccfg.ctrl_idx = {DE, DT};
  ccfg.ctrl_weights = {sc_.cbf.w_de, sc_.cbf.w_dT};
  ccfg.slack_penalty = sc_.cbf.slack_penalty;
  cbf_ = CBFFilter(ccfg, ac_.limits);
  cbf_.setEnabled(sc_.cbf.enabled);

  buildBarriers();
}

void Sim::buildBarriers() {
  barriers_.clear();
  const CBFParams& p = sc_.cbf;

  // Airspeed (stall) barrier  b_V = V - V_min  (relative degree 1, works).
  if (p.airspeed_barrier) {
    const double V_min = (p.V_min > 0.0) ? p.V_min : 0.85 * sc_.V_app;
    const double aV = p.alpha_V;
    Barrier b{
        "min_airspeed",
        [V_min](const StateVec& x) { return airspeed(x) - V_min; },
        [](const StateVec& x) {
          Eigen::RowVectorXd gx = Eigen::RowVectorXd::Zero(NX);
          const double Vt = airspeed(x);  // not 'V' (== State::V enum index)
          gx[U] = x[U] / Vt; gx[V] = x[V] / Vt; gx[W] = x[W] / Vt;
          return gx;
        },
        [aV](double hv) { return aV * hv; }};
    b.hard = p.airspeed_hard;
    barriers_.push_back(std::move(b));
  }

  // Descent-rate (soft-landing) barrier  b = hdot + sqrt(v_safe^2 + 2 a_brk h).
  // Enforced HARD by default (cbf.descent_hard): no slack, the CBF row must hold
  // exactly. It is high relative degree under direct elevator, so when it cannot
  // be met the filter recovers by softening for that step (CBFFilter); a proper
  // HOCBF is still the path to authoritative enforcement + emergent flare
  // (documentation/water_landing_cbf.md A.3).
  if (p.descent_barrier) {
    const double vsafe = p.v_safe, abrk = p.a_brk, ah = p.alpha_h;
    std::function<double(const StateVec&)> hfun =
        [vsafe, abrk](const StateVec& x) {
          const double hgt = std::max(0.0, x[H]);
          const double hdot = -sinkOf(x);  // doc v = hdot < 0 on descent
          return hdot + std::sqrt(vsafe * vsafe + 2.0 * abrk * hgt);
        };
    Barrier b{"descent_rate", hfun,
              [hfun](const StateVec& x) { return numGrad(hfun, x); },
              [ah](double hv) { return ah * hv; }};
    b.hard = p.descent_hard;
    barriers_.push_back(std::move(b));
  }
}

std::unique_ptr<NominalController> Sim::makeNominal() const {
  switch (sc_.nominal) {
    case NominalKind::LonGlideslope:
      return std::make_unique<LonGlideslopeController>(
          sc_.gains, trim_, sc_.V_app, sc_.gamma_app, ac_.limits);
    case NominalKind::CascadedPID:
    default:
      return std::make_unique<Controller>(sc_.gains, sc_.flare, trim_,
                                          sc_.V_app, sc_.gamma_app,
                                          sc_.h_decrab, ac_.env, ac_.limits);
  }
}

TouchdownSummary Sim::run(const std::string& csv_path) {
  // ---- print trim ----------------------------------------------------------
  std::cout << "==== TRIM (V_app=" << sc_.V_app << " m/s, gamma="
            << sc_.gamma_app * kRad2Deg << " deg) ====\n";
  std::cout << std::fixed << std::setprecision(5);
  std::cout << "  converged   : " << (trim_.converged ? "yes" : "NO")
            << "  (residual " << trim_.residual << ", " << trim_.iterations
            << " iters)\n";
  std::cout << "  alpha       : " << trim_.alpha * kRad2Deg << " deg\n";
  std::cout << "  theta       : " << trim_.theta * kRad2Deg << " deg\n";
  std::cout << "  delta_e     : " << trim_.u[DE] * kRad2Deg << " deg\n";
  std::cout << "  throttle    : " << trim_.u[DT] << "\n";
  if (!trim_.converged)
    std::cout << "  WARNING: trim did not converge to tolerance; check mass "
                 "properties / thrust model in aircraft.yaml.\n";

  // ---- initial condition on a stabilized approach ---------------------------
  StateVec x = trim_.x;
  double range = sc_.init.distance;  // horizontal distance to touchdown [m]
  x[H] = range * std::tan(std::abs(sc_.gamma_app)) + sc_.init.dh;
  x[Y] = sc_.init.dy;
  x[PSI] = trim_.x[PSI] + sc_.init.dpsi;
  x[U] = trim_.x[U] + sc_.init.dV;  // perturb airspeed

  std::unique_ptr<NominalController> ctrl = makeNominal();

  // Control-affine view of the linear model, consumed by the CBF filter. Swap
  // this for any other ControlAffineModel (reduced or nonlinear) with no other
  // change to the loop.
  LinearModelCAM cam(lin_);

  // ---- CSV log --------------------------------------------------------------
  std::ofstream csv(csv_path);
  if (!csv) throw std::runtime_error("sim: cannot open CSV for writing: " +
                                     csv_path);
  csv << "t,u,v,w,p,q,r,phi_deg,theta_deg,psi_deg,h,y,"
         "V,alpha_deg,beta_deg,gamma_deg,sink,"
         "de_deg,da_deg,dr_deg,dT,"
         "theta_cmd_deg,phi_cmd_deg,h_ref,w_cmd,V_cmd,flaring,range,"
         "de_nom_deg,dT_nom,b_airspeed,b_descent\n";

  // Airspeed/descent barrier values for logging (independent of enable flags),
  // using the configured CBF parameters. b >= 0 is the safe set.
  const double V_min_log =
      (sc_.cbf.V_min > 0.0) ? sc_.cbf.V_min : 0.85 * sc_.V_app;
  csv << std::fixed << std::setprecision(6);

  auto plant = [&](const StateVec& xs, const CtrlVec& us) {
    return lin_.xdot(xs, us);
  };

  TouchdownSummary td;
  const double dt = sc_.dt;
  double t = 0.0;
  StateVec x_prev = x;

  // Track how many steps a HARD barrier could not be met and the filter had to
  // soften (a guarantee was unavailable on that step).
  int filter_steps = 0, recovery_steps = 0;

  const int nsteps = static_cast<int>(sc_.t_max / dt);
  for (int k = 0; k <= nsteps; ++k) {
    // Nominal controller (full-state feedback) then CBF safety filter.
    ControllerOutputs co = ctrl->step(x, range, dt);
    const CtrlVec u_nom = co.u;  // nominal command, before the safety filter
    CtrlVec u = cbf_.filter(u_nom, x, cam, barriers_);
    ++filter_steps;
    if (cbf_.lastFeasibilityRecovery()) ++recovery_steps;

    const double Vt = airspeed(x);  // not 'V' (== State::V enum index)
    const double gamma = x[THETA] - alphaOf(x);
    // Barrier values at the current state (b >= 0 safe).
    const double b_air = Vt - V_min_log;
    const double b_des = -sinkOf(x) + std::sqrt(sc_.cbf.v_safe * sc_.cbf.v_safe +
                                                2.0 * sc_.cbf.a_brk *
                                                    std::max(0.0, x[H]));
    csv << t << "," << x[U] << "," << x[V] << "," << x[W] << "," << x[P] << ","
        << x[Q] << "," << x[R] << "," << x[PHI] * kRad2Deg << ","
        << x[THETA] * kRad2Deg << "," << x[PSI] * kRad2Deg << "," << x[H] << ","
        << x[Y] << "," << Vt << "," << alphaOf(x) * kRad2Deg << ","
        << betaOf(x) * kRad2Deg << "," << gamma * kRad2Deg << "," << sinkOf(x)
        << "," << u[DE] * kRad2Deg << "," << u[DA] * kRad2Deg << ","
        << u[DR] * kRad2Deg << "," << u[DT] << "," << co.theta_cmd * kRad2Deg
        << "," << co.phi_cmd * kRad2Deg << "," << co.h_ref << "," << co.w_cmd
        << "," << co.V_cmd << "," << (co.flaring ? 1 : 0) << "," << range << ","
        << u_nom[DE] * kRad2Deg << "," << u_nom[DT] << "," << b_air << ","
        << b_des << "\n";

    // Touchdown: altitude crosses zero.
    if (x[H] <= 0.0 && k > 0) {
      // Linear interpolation between the previous and current sample.
      const double h0 = x_prev[H], h1 = x[H];
      const double frac = (h1 != h0) ? h0 / (h0 - h1) : 0.0;
      StateVec xt = x_prev + frac * (x - x_prev);
      td.reached = true;
      td.t = t - dt + frac * dt;
      td.sink_rate = sinkOf(xt);
      td.V = airspeed(xt);
      td.theta = xt[THETA];
      td.phi = xt[PHI];
      td.beta = betaOf(xt);
      td.y = xt[Y];
      break;
    }

    // Advance the plant (RK4, zero-order-hold control).
    x_prev = x;
    x = rk4Step(plant, x, u, dt);
    range -= northVel(x) * dt;  // downrange progress
    t += dt;
  }
  csv.close();

  // ---- touchdown summary ----------------------------------------------------
  std::cout << "\n==== TOUCHDOWN SUMMARY ====\n";
  if (td.reached) {
    std::cout << "  time        : " << td.t << " s\n";
    std::cout << "  sink rate   : " << td.sink_rate << " m/s (positive down)\n";
    std::cout << "  airspeed    : " << td.V << " m/s\n";
    std::cout << "  pitch       : " << td.theta * kRad2Deg << " deg\n";
    std::cout << "  bank        : " << td.phi * kRad2Deg << " deg\n";
    std::cout << "  sideslip    : " << td.beta * kRad2Deg << " deg\n";
    std::cout << "  cross-track : " << td.y << " m\n";
  } else {
    std::cout << "  touchdown NOT reached within t_max=" << sc_.t_max
              << " s.\n";
  }

  // ---- CBF filter health ----------------------------------------------------
  if (cbf_.enabled() && !barriers_.empty()) {
    bool any_hard = false;
    for (const Barrier& b : barriers_) any_hard = any_hard || b.hard;
    if (any_hard) {
      const double pct =
          filter_steps ? 100.0 * recovery_steps / filter_steps : 0.0;
      std::cout << "\n==== CBF (hard barriers) ====\n";
      std::cout << std::setprecision(1)
                << "  feasibility recovery: " << recovery_steps << " / "
                << filter_steps << " steps (" << pct
                << "%) softened a hard barrier\n"
                << std::setprecision(5);
      if (recovery_steps > 0)
        std::cout << "  NOTE: a hard barrier could not be met by the controls "
                     "on those steps\n"
                     "        (high relative degree); an HOCBF is the path to "
                     "authoritative enforcement.\n";
    }
  }

  std::cout << "\nLog written to: " << csv_path << "\n";
  return td;
}

}  // namespace autoland
