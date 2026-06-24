#include "autoland/sim.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

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
  buildCandidateBarriers();
}

void Sim::buildCandidateBarriers() {
  // Candidate barriers (passed to the pass-through filter for now). These show
  // the intended h(x)/grad(x) shape the OSQP-backed filter will consume.
  const double V_stall_margin = 0.85 * sc_.V_app;  // TODO: real V_stall
  barriers_.push_back(Barrier{
      "min_airspeed",
      [V_stall_margin](const StateVec& x) {
        return airspeed(x) - V_stall_margin;
      },
      [](const StateVec& x) {
        Eigen::RowVectorXd gx = Eigen::RowVectorXd::Zero(NX);
        const double Vt = airspeed(x);  // not 'V' (== State::V enum index)
        gx[U] = x[U] / Vt; gx[V] = x[V] / Vt; gx[W] = x[W] / Vt;
        return gx;
      },
      [](double hv) { return 1.0 * hv; }});
  // Further candidates (bank limit vs height, sink-rate vs height) to be added
  // alongside the QP implementation; see cbf.hpp.
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

  Controller ctrl(sc_.gains, sc_.flare, trim_, sc_.V_app, sc_.gamma_app,
                  sc_.h_decrab, ac_.env, ac_.limits);

  // ---- CSV log --------------------------------------------------------------
  std::ofstream csv(csv_path);
  if (!csv) throw std::runtime_error("sim: cannot open CSV for writing: " +
                                     csv_path);
  csv << "t,u,v,w,p,q,r,phi_deg,theta_deg,psi_deg,h,y,"
         "V,alpha_deg,beta_deg,gamma_deg,sink,"
         "de_deg,da_deg,dr_deg,dT,"
         "theta_cmd_deg,phi_cmd_deg,h_ref,w_cmd,V_cmd,flaring,range\n";
  csv << std::fixed << std::setprecision(6);

  auto plant = [&](const StateVec& xs, const CtrlVec& us) {
    return lin_.xdot(xs, us);
  };

  TouchdownSummary td;
  const double dt = sc_.dt;
  double t = 0.0;
  StateVec x_prev = x;

  const int nsteps = static_cast<int>(sc_.t_max / dt);
  for (int k = 0; k <= nsteps; ++k) {
    // Controller (full-state feedback) then CBF filter (pass-through stub).
    ControllerOutputs co = ctrl.step(x, range, dt);
    CtrlVec u = cbf_.filter(co.u, x, lin_, barriers_);

    const double Vt = airspeed(x);  // not 'V' (== State::V enum index)
    const double gamma = x[THETA] - alphaOf(x);
    csv << t << "," << x[U] << "," << x[V] << "," << x[W] << "," << x[P] << ","
        << x[Q] << "," << x[R] << "," << x[PHI] * kRad2Deg << ","
        << x[THETA] * kRad2Deg << "," << x[PSI] * kRad2Deg << "," << x[H] << ","
        << x[Y] << "," << Vt << "," << alphaOf(x) * kRad2Deg << ","
        << betaOf(x) * kRad2Deg << "," << gamma * kRad2Deg << "," << sinkOf(x)
        << "," << u[DE] * kRad2Deg << "," << u[DA] * kRad2Deg << ","
        << u[DR] * kRad2Deg << "," << u[DT] << "," << co.theta_cmd * kRad2Deg
        << "," << co.phi_cmd * kRad2Deg << "," << co.h_ref << "," << co.w_cmd
        << "," << co.V_cmd << "," << (co.flaring ? 1 : 0) << "," << range
        << "\n";

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
  std::cout << "\nLog written to: " << csv_path << "\n";
  return td;
}

}  // namespace autoland
