#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <string>
#include <vector>
#include "autoland/world_geometry.hpp"

// =============================================================================
// Pattern guidance: the target-tracking "pilot" nominal for the corridor
// landing. Plans ONCE at t = 0 a path from the current ground pose to a
// final-approach fix (FAF) on the corridor axis, extended along the axis
// through the aim point and past the far edge, then tracks it:
//
//   lateral      : arclength-carrot pursuit with the Park-Deyst-How gain
//                  (cf. Park, Deyst & How, AIAA GNC 2004; the L1 = ratio *
//                  period * V / pi parameterization follows ArduPilot's
//                  AP_L1_Control). The carrot sits at path ARCLENGTH L1
//                  ahead of the foot point (PDH intersect a circle of radius
//                  L1 with the path); the lateral acceleration
//                  a_s = 2 V^2 sin(eta) / Ld uses the actual chord Ld to the
//                  carrot, which makes a_s = V^2/R EXACT on any circular
//                  segment (Ld = 2R sin(eta)). Linearized about a straight
//                  path: omega_n = sqrt(2) V / L1, zeta = 1/sqrt(2) -- the
//                  "period"/"scale" parameters only set L1. -> phi_cmd for
//                  the existing roll inner loop;
//   longitudinal : glideslope cone on distance-to-go to the aim point,
//                  h_gs = L tan(-gamma_app), gamma_ref = gamma_app +
//                  K_h (h_gs - h) / Vg: the height error obeys
//                  e_dot = -K_h e exactly (tau = 1/K_h), fed as a gamma_ref
//                  into the existing cascade / TECS outer loop.
//
// FRAME: earth N/E, course angle chi from N toward E (right turn = chi
// increasing = positive bank). Ground course, not heading, is tracked so a
// crosswind crab falls out for free (the plan itself is ground-frame with a
// fixed radius, so in wind the bank needed on an arc varies with the ground
// speed -- no wind-aware planning yet).
//
// The FAF distance d_final is chosen by a 25 m GRID SCAN so the pattern
// length matches the descent distance h0 / tan(-gamma_app) as closely as
// possible (whole pattern near gamma_app when the geometry allows; the cone
// absorbs the rest). Not a fixed point / bisection: the path length is
// discontinuous in the goal pose.
//
// Path words: the four Dubins CSC words (LSL, RSR, LSR, RSL) plus a DIRECT
// line when the aircraft is already nearly on the axis ahead of the FAF.
// Shortest OF THESE is chosen, which is the Dubins shortest path only when
// the poses are far enough apart (Shkel & Lumelsky 2001: CCC words can only
// be optimal when the same-sense circle centres are <= 4R apart, which the
// scan's FAF goals often are); minimality is not needed here, the objective
// is length matching. Candidates with a near-full-circle arc (the 2 pi wrap
// of a nearly aligned goal just ahead) are rejected. Tangent geometry is
// derived directly (see pattern_guidance.cpp) and verified against forward
// integration in the tests.
// =============================================================================
namespace autoland {

struct Pose2 {
  double N{0}, E{0}, chi{0};
};

enum class SegType { Line, ArcL, ArcR };

struct PathSegment {
  SegType type{SegType::Line};
  double N0{0}, E0{0}, chi0{0};  // start pose
  double length{0};              // arc length [m]
  double R{0};                   // radius (arcs)
  double Nc{0}, Ec{0};           // center (arcs)
  Pose2 poseAt(double s) const;  // s in [0, length] (clamped)
  // Closest point on the segment to (N, E): returns s and the distance.
  void closest(double N, double E, double& s, double& dist) const;
};

struct DubinsPath {
  bool valid{false};
  std::string word;
  double length{0};
  std::vector<PathSegment> segs;
};

// All four CSC candidates (invalid ones have valid = false) and the
// shortest valid one.
std::vector<DubinsPath> dubinsCSCAll(const Pose2& start, const Pose2& goal,
                                     double R);
DubinsPath dubinsCSC(const Pose2& start, const Pose2& goal, double R);
// Single straight segment from start to goal (course = bearing); the goal
// heading is NOT matched -- valid only when the poses are nearly aligned.
DubinsPath directLine(const Pose2& start, const Pose2& goal);

// Polyline of segments with cumulative arc length.
class Path {
 public:
  void clear() { segs_.clear(); cum_.clear(); total_ = 0.0; }
  void append(const PathSegment& s);
  const std::vector<PathSegment>& segs() const { return segs_; }
  double length() const { return total_; }
  double segStart(int i) const { return cum_[i]; }
  Pose2 poseAt(double s_path) const;
  // Closest point searching only within s_min <= s_path <= s_max (progress
  // window: the pattern may pass over / cross its own final leg, e.g. an
  // overhead entry sits ON the corridor axis at t = 0). Returns the path
  // coordinate; sets seg + dist.
  double closest(double N, double E, double s_min, double s_max, int* seg,
                 double* dist) const;
 private:
  std::vector<PathSegment> segs_;
  std::vector<double> cum_;
  double total_{0.0};
};

struct PatternGuidanceConfig {
  double V_app{40.0};                    // planning speed [m/s]
  double gamma_app{-4.0 * M_PI / 180.0}; // approach glideslope [rad]
  double g{9.80665};
  // Dubins radius: R = k_R V_app^2 / (g tan phi_plan).
  double phi_plan{20.0 * M_PI / 180.0};
  double k_R{1.15};
  // Final leg (FAF to aim point) bounds for the fixed point [m].
  double d_final_min{800.0}, d_final_max{4000.0};
  double axis_overrun{2000.0};           // axis line extends this far past the far edge
  // Glideslope cone: gamma_ref = gamma_app + K_h (h_gs - h)/V, clamped.
  double K_h{0.1};                       // [1/s]
  double gamma_min{-8.0 * M_PI / 180.0};
  double gamma_max{6.0 * M_PI / 180.0};
  // Carrot distance L1 = max(L1_min, L1_scale * L1_period * Vg / pi)
  // (ArduPilot parameterization; NOTE the closed loop has zeta = 1/sqrt(2)
  // and omega_n = sqrt(2) Vg / L1 regardless of these values -- they only
  // size L1. Defaults: L1 = 239 m at 40 m/s, natural period 26.5 s.)
  double L1_period{25.0}, L1_scale{0.75}, L1_min{100.0};
  // Direct-line candidate tolerances: the aircraft is "on the axis" when
  // its cross-axis offset and course error are within these.
  double direct_e_tol{30.0}, direct_chi_tol{5.0 * M_PI / 180.0};
  // Reject CSC candidates containing an arc sweeping more than this (a
  // nearly aligned goal just ahead makes LSL/RSR wrap a full circle).
  double max_arc_sweep{350.0 * M_PI / 180.0};
  double phi_max{25.0 * M_PI / 180.0};
  // Progress window for the closest-point search [m]: the path coordinate
  // may fall back at most `s_back` and advance at most `s_ahead` per step.
  double s_back{50.0}, s_ahead{500.0};
};

struct GuidanceCmd {
  double gamma_ref{0};   // [rad]
  double phi_cmd{0};     // [rad]
  double chi_cmd{0};     // bearing to the lookahead point [rad]
  double chi{0};         // ground course used by the tracker [rad]
  double N_ref{0}, E_ref{0};   // lookahead point
  double h_gs{0};        // glideslope height at the current distance-to-go
  double L_to_go{0};     // path distance to the aim point [m] (negative past it)
  double s_path{0};      // progress along the path [m]
  double xte{0};         // distance to the closest path point [m] (unsigned)
  int seg{0};
  bool on_final{false};  // past the FAF
};

class PatternGuidance {
 public:
  PatternGuidance(const PatternGuidanceConfig& cfg, const Corridor& cor)
      : c_(cfg), cor_(cor) {}

  // Plan from the current ground pose and altitude. Returns false if no
  // valid path (cannot happen for CSC, kept for safety).
  bool plan(double N, double E, double chi, double h);

  GuidanceCmd step(double N, double E, double h, double xdot_n, double ydot_e);

  const Path& path() const { return path_; }
  double dFinal() const { return d_final_; }
  double sFafPath() const { return s_faf_; }
  double sAimPath() const { return s_aim_; }
  // Progress along the path (closest-point coordinate, monotone). Setting it
  // re-centers the search window (tests / re-entry after a go-around).
  double progress() const { return s_prog_; }
  void setProgress(double s) { s_prog_ = s; }
  double turnRadius() const;
  const std::string& word() const { return word_; }

 private:
  PatternGuidanceConfig c_;
  Corridor cor_;
  Path path_;
  std::string word_;
  double d_final_{0}, s_faf_{0}, s_aim_{0};
  double s_prog_{0};
  bool planned_{false};
};

}  // namespace autoland
