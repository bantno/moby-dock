#include "autoland/pattern_guidance.hpp"
#include <algorithm>
#include <limits>

// -----------------------------------------------------------------------------
// Geometry conventions (N, E), course chi from N toward E.
//   heading vector           u(chi) = (cos chi, sin chi)
//   right of the heading     chi + pi/2  -> (-sin chi, cos chi)
// A RIGHT turn increases chi; its center is R to the right of the pose and the
// point sits at bearing chi - pi/2 from the center:
//   p = c + R (sin chi, -cos chi)
// A LEFT turn decreases chi; center to the left, point at bearing chi + pi/2:
//   p = c + R (-sin chi, cos chi)
// CSC tangents: the straight leg has course theta; the tangent points are the
// points on each circle where the course equals theta (formulas above). With
// d = c2 - c1 = D (cos phi, sin phi) and n = (sin theta, -cos theta):
//   RSR / LSL : d = L u                 -> theta = phi,             L = D
//   RSL       : d = L u + 2R n          -> sin(theta - phi) =  2R/D, L = sqrt(D^2 - 4R^2)
//   LSR       : d = L u - 2R n          -> sin(theta - phi) = -2R/D, L = sqrt(D^2 - 4R^2)
// -----------------------------------------------------------------------------
namespace autoland {

namespace {
double mod2pi(double a) {
  a = std::fmod(a, 2.0 * M_PI);
  if (a < 0) a += 2.0 * M_PI;
  return a;
}
double wrapPi(double a) { return std::atan2(std::sin(a), std::cos(a)); }

Eigen::Vector2d rightCenter(const Pose2& p, double R) {
  return {p.N - R * std::sin(p.chi), p.E + R * std::cos(p.chi)};
}
Eigen::Vector2d leftCenter(const Pose2& p, double R) {
  return {p.N + R * std::sin(p.chi), p.E - R * std::cos(p.chi)};
}

PathSegment arcSeg(const Pose2& p, double R, bool right, double sweep) {
  PathSegment s;
  s.type = right ? SegType::ArcR : SegType::ArcL;
  s.N0 = p.N; s.E0 = p.E; s.chi0 = p.chi;
  s.R = R;
  s.length = R * sweep;
  const Eigen::Vector2d c = right ? rightCenter(p, R) : leftCenter(p, R);
  s.Nc = c[0]; s.Ec = c[1];
  return s;
}
PathSegment lineSeg(const Pose2& p, double L) {
  PathSegment s;
  s.type = SegType::Line;
  s.N0 = p.N; s.E0 = p.E; s.chi0 = p.chi;
  s.length = L;
  return s;
}

// One CSC candidate. first/last: true = right turn.
DubinsPath csc(const Pose2& a, const Pose2& b, double R, bool first_right,
               bool last_right, const std::string& word) {
  DubinsPath out;
  out.word = word;
  const Eigen::Vector2d c1 = first_right ? rightCenter(a, R) : leftCenter(a, R);
  const Eigen::Vector2d c2 = last_right ? rightCenter(b, R) : leftCenter(b, R);
  const Eigen::Vector2d d = c2 - c1;
  const double D = d.norm();
  const double phi = std::atan2(d[1], d[0]);
  double theta, L;
  if (first_right == last_right) {
    // Coincident circles (D -> 0): no straight leg, keep the start course
    // so the sweeps do not pick up a spurious 2 pi from atan2(0, 0).
    theta = D < 1e-9 ? a.chi : phi;
    L = D;
  } else {
    if (D < 2.0 * R) return out;  // inner tangent does not exist
    const double s = std::clamp(2.0 * R / D, -1.0, 1.0);
    theta = first_right ? phi + std::asin(s) : phi - std::asin(s);
    L = std::sqrt(std::max(0.0, D * D - 4.0 * R * R));
  }
  // Arc sweeps: right turn increases chi, left decreases it.
  const double sw1 = first_right ? mod2pi(theta - a.chi) : mod2pi(a.chi - theta);
  const double sw3 = last_right ? mod2pi(b.chi - theta) : mod2pi(theta - b.chi);
  // Snap near-full circles that are really zero (numerical 2pi wrap).
  auto snap = [](double sw) { return (2.0 * M_PI - sw < 1e-9) ? 0.0 : sw; };

  const PathSegment s1 = arcSeg(a, R, first_right, snap(sw1));
  const Pose2 p1 = s1.poseAt(s1.length);
  const PathSegment s2 = lineSeg(p1, L);
  const Pose2 p2 = s2.poseAt(L);
  const PathSegment s3 = arcSeg(p2, R, last_right, snap(sw3));
  out.segs = {s1, s2, s3};
  out.length = s1.length + s2.length + s3.length;
  out.valid = true;
  return out;
}
}  // namespace

Pose2 PathSegment::poseAt(double s) const {
  s = std::clamp(s, 0.0, length);
  Pose2 p;
  switch (type) {
    case SegType::Line:
      p.N = N0 + s * std::cos(chi0);
      p.E = E0 + s * std::sin(chi0);
      p.chi = chi0;
      break;
    case SegType::ArcR:
      p.chi = chi0 + s / R;
      p.N = Nc + R * std::sin(p.chi);
      p.E = Ec - R * std::cos(p.chi);
      break;
    case SegType::ArcL:
      p.chi = chi0 - s / R;
      p.N = Nc - R * std::sin(p.chi);
      p.E = Ec + R * std::cos(p.chi);
      break;
  }
  return p;
}

void PathSegment::closest(double N, double E, double& s_best,
                          double& dist) const {
  auto dAt = [&](double s) {
    const Pose2 p = poseAt(s);
    return std::hypot(N - p.N, E - p.E);
  };
  double s_proj;
  if (type == SegType::Line) {
    s_proj = (N - N0) * std::cos(chi0) + (E - E0) * std::sin(chi0);
  } else {
    // Bearing of the point from the center -> course on the circle -> s.
    const double bear = std::atan2(E - Ec, N - Nc);
    const double chi_pt = (type == SegType::ArcR) ? bear + M_PI / 2 : bear - M_PI / 2;
    const double dchi = (type == SegType::ArcR) ? mod2pi(chi_pt - chi0)
                                               : mod2pi(chi0 - chi_pt);
    s_proj = R * dchi;
  }
  s_best = 0.0;
  dist = dAt(0.0);
  for (double s : {std::clamp(s_proj, 0.0, length), length}) {
    const double d = dAt(s);
    if (d < dist) { dist = d; s_best = s; }
  }
}

std::vector<DubinsPath> dubinsCSCAll(const Pose2& start, const Pose2& goal,
                                     double R) {
  return {csc(start, goal, R, false, false, "LSL"),
          csc(start, goal, R, true, true, "RSR"),
          csc(start, goal, R, false, true, "LSR"),
          csc(start, goal, R, true, false, "RSL")};
}

DubinsPath dubinsCSC(const Pose2& start, const Pose2& goal, double R) {
  DubinsPath best;
  best.length = std::numeric_limits<double>::infinity();
  for (const auto& cand : dubinsCSCAll(start, goal, R))
    if (cand.valid && cand.length < best.length) best = cand;
  return best;
}

DubinsPath directLine(const Pose2& start, const Pose2& goal) {
  DubinsPath out;
  out.word = "DIRECT";
  const double L = std::hypot(goal.N - start.N, goal.E - start.E);
  Pose2 p = start;
  p.chi = std::atan2(goal.E - start.E, goal.N - start.N);
  out.segs = {lineSeg(p, L)};
  out.length = L;
  out.valid = L > 1e-9;
  return out;
}

void Path::append(const PathSegment& s) {
  cum_.push_back(total_);
  segs_.push_back(s);
  total_ += s.length;
}

Pose2 Path::poseAt(double s_path) const {
  if (segs_.empty()) return {};
  s_path = std::clamp(s_path, 0.0, total_);
  int i = static_cast<int>(segs_.size()) - 1;
  while (i > 0 && s_path < cum_[i]) --i;
  return segs_[i].poseAt(s_path - cum_[i]);
}

double Path::closest(double N, double E, double s_min, double s_max,
                     int* seg, double* dist) const {
  double best_s = 0.0, best_d = std::numeric_limits<double>::infinity();
  int best_i = 0;
  s_min = std::max(0.0, s_min);
  s_max = std::min(total_, std::max(s_max, s_min));
  for (int i = 0; i < static_cast<int>(segs_.size()); ++i) {
    if (cum_[i] + segs_[i].length < s_min || cum_[i] > s_max) continue;
    double s, d;
    segs_[i].closest(N, E, s, d);
    const double sp = std::clamp(cum_[i] + s, s_min, s_max);
    if (sp != cum_[i] + s) {  // clipped by the window: re-measure
      const Pose2 p = poseAt(sp);
      d = std::hypot(N - p.N, E - p.E);
    }
    if (d < best_d) { best_d = d; best_s = sp; best_i = i; }
  }
  if (seg) *seg = best_i;
  if (dist) *dist = best_d;
  return best_s;
}

double PatternGuidance::turnRadius() const {
  return c_.k_R * c_.V_app * c_.V_app / (c_.g * std::tan(c_.phi_plan));
}

bool PatternGuidance::plan(double N, double E, double chi, double h) {
  if (!(c_.gamma_app < 0.0)) return false;  // descent distance undefined
  const double R = turnRadius();
  const Pose2 start{N, E, chi};
  const double descent = h / std::tan(-c_.gamma_app);  // ground distance at gamma_app
  // Is the aircraft already on the axis (direct-line candidate allowed)?
  double s0, e0;
  cor_.frame(N, E, s0, e0);
  const bool on_axis = std::abs(e0) <= c_.direct_e_tol &&
                       std::abs(wrapPi(chi - cor_.heading)) <= c_.direct_chi_tol;
  // Choose d_final so that d_final + L_dubins(d_final) is as close as
  // possible to the descent distance: a GRID SCAN (25 m) over
  // [d_final_min, d_final_max], ties broken toward the longer final leg.
  // Not a bisection: Dubins length is discontinuous in the goal pose (a
  // nearly aligned FAF just ahead of the aircraft needs a full 2 pi R loop
  // when the heading mismatch cannot be absorbed by an inner tangent), so
  // the total is neither continuous nor monotone. The scan simply steps
  // over such loops.
  // Best candidate for a given FAF: the shortest CSC word without a
  // near-full-circle arc, or the direct line when the aircraft is on the
  // axis ahead of the FAF (the direct line wins ties by being shortest).
  auto planTo = [&](double d) {
    const double s_faf = cor_.s_aim - d;
    const Eigen::Vector2d faf = cor_.axisPoint(s_faf);
    const Pose2 goal{faf[0], faf[1], cor_.heading};
    DubinsPath best;
    best.length = std::numeric_limits<double>::infinity();
    for (const auto& cand : dubinsCSCAll(start, goal, R)) {
      if (!cand.valid) continue;
      if (cand.segs[0].length > R * c_.max_arc_sweep ||
          cand.segs[2].length > R * c_.max_arc_sweep) continue;
      if (cand.length < best.length) best = cand;
    }
    if (on_axis && s_faf > s0 + 1.0) {
      const DubinsPath dl = directLine(start, goal);
      if (dl.valid && dl.length < best.length) best = dl;
    }
    return best;
  };
  DubinsPath dp;
  double d_final = c_.d_final_min;
  double best = std::numeric_limits<double>::infinity();
  const double step = 25.0;
  for (double d = c_.d_final_min; d <= c_.d_final_max + 1e-9; d += step) {
    const DubinsPath cand = planTo(d);
    if (!cand.valid) continue;
    const double err = std::abs(d + cand.length - descent);
    // Ties (within the grid resolution) go to the longer final leg.
    if (err < best - 1e-9 || (std::abs(err - best) <= 1e-9 && d > d_final)) {
      best = err;
      d_final = d;
      dp = cand;
    }
  }
  if (!dp.valid) return false;
  path_.clear();
  for (const auto& s : dp.segs) path_.append(s);
  s_faf_ = path_.length();
  const Pose2 faf_pose = path_.poseAt(s_faf_);
  // Final leg: axis line from the FAF through the aim point, past the far edge.
  const double L_axis = d_final + (cor_.sFar() - cor_.s_aim) + c_.axis_overrun;
  path_.append(lineSeg(Pose2{faf_pose.N, faf_pose.E, cor_.heading}, L_axis));
  s_aim_ = s_faf_ + d_final;
  d_final_ = d_final;
  word_ = dp.word;
  s_prog_ = 0.0;
  planned_ = true;
  return true;
}

GuidanceCmd PatternGuidance::step(double N, double E, double h, double xdot_n,
                                  double ydot_e) {
  GuidanceCmd g;
  if (!planned_) return g;
  const double Vg = std::max(1.0, std::hypot(xdot_n, ydot_e));
  const double chi = std::atan2(ydot_e, xdot_n);

  // Progress along the path within the [s_back, s_ahead] window.
  int seg = 0;
  double dist = 0.0;
  const double s_c = path_.closest(N, E, s_prog_ - c_.s_back,
                                   s_prog_ + c_.s_ahead, &seg, &dist);
  s_prog_ = std::max(s_prog_, s_c);
  g.s_path = s_c;
  g.xte = dist;
  g.seg = seg;
  g.on_final = s_c >= s_faf_;
  g.L_to_go = s_aim_ - s_c;

  // --- L1 lateral: lookahead point L1 ahead along the path. ----------------
  const double L1 = std::max(c_.L1_min, c_.L1_scale * c_.L1_period * Vg / M_PI);
  const Pose2 la = path_.poseAt(s_c + L1);
  g.N_ref = la.N; g.E_ref = la.E;
  const double bearing = std::atan2(la.E - E, la.N - N);
  const double eta = wrapPi(bearing - chi);
  const double Ld = std::max(1.0, std::hypot(la.N - N, la.E - E));
  const double a_s = 2.0 * Vg * Vg * std::sin(eta) / Ld;
  g.phi_cmd = std::clamp(std::atan(a_s / c_.g), -c_.phi_max, c_.phi_max);
  g.chi_cmd = bearing;
  g.chi = chi;

  // --- Glideslope cone on distance-to-go. -----------------------------------
  g.h_gs = std::max(0.0, g.L_to_go) * std::tan(-c_.gamma_app);
  g.gamma_ref = std::clamp(c_.gamma_app + c_.K_h * (g.h_gs - h) / Vg,
                           c_.gamma_min, c_.gamma_max);
  return g;
}

}  // namespace autoland
