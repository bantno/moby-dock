#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <random>
#include "autoland/pattern_guidance.hpp"

using namespace autoland;
using Catch::Approx;

namespace {
constexpr double kDeg = M_PI / 180.0;

// Independent forward integration of the unicycle along the path word:
// N' = cos chi, E' = sin chi, chi' = +1/R (R), -1/R (L), 0 (S). RK4 with
// fine steps -- no use of poseAt.
Pose2 integratePath(const DubinsPath& p, const Pose2& start, double R) {
  Pose2 x = start;
  for (const auto& seg : p.segs) {
    const double kappa = seg.type == SegType::ArcR ? 1.0 / R
                       : seg.type == SegType::ArcL ? -1.0 / R : 0.0;
    const int n = std::max(1, static_cast<int>(seg.length / 0.05));
    const double h = seg.length / n;
    auto f = [&](const Pose2& q) {
      return Pose2{std::cos(q.chi), std::sin(q.chi), kappa};
    };
    auto add = [](const Pose2& a, double s, const Pose2& d) {
      return Pose2{a.N + s * d.N, a.E + s * d.E, a.chi + s * d.chi};
    };
    for (int i = 0; i < n; ++i) {
      const Pose2 k1 = f(x), k2 = f(add(x, 0.5 * h, k1)),
                  k3 = f(add(x, 0.5 * h, k2)), k4 = f(add(x, h, k3));
      x = Pose2{x.N + h / 6 * (k1.N + 2 * k2.N + 2 * k3.N + k4.N),
                x.E + h / 6 * (k1.E + 2 * k2.E + 2 * k3.E + k4.E),
                x.chi + h / 6 * (k1.chi + 2 * k2.chi + 2 * k3.chi + k4.chi)};
    }
  }
  return x;
}
double wrap(double a) { return std::atan2(std::sin(a), std::cos(a)); }
}  // namespace

TEST_CASE("Dubins CSC reaches the goal pose (independent integration)",
          "[guidance]") {
  std::mt19937 rng(11);
  std::uniform_real_distribution<double> un(-1.0, 1.0);
  const double R = 500.0;
  int words[4] = {0, 0, 0, 0};
  for (int k = 0; k < 40; ++k) {
    const Pose2 a{3000.0 * un(rng), 3000.0 * un(rng), M_PI * un(rng)};
    const Pose2 b{3000.0 * un(rng), 3000.0 * un(rng), M_PI * un(rng)};
    const DubinsPath p = dubinsCSC(a, b, R);
    REQUIRE(p.valid);
    REQUIRE(p.segs.size() == 3);
    const Pose2 end = integratePath(p, a, R);
    CHECK(end.N == Approx(b.N).margin(1e-6));
    CHECK(end.E == Approx(b.E).margin(1e-6));
    CHECK(wrap(end.chi - b.chi) == Approx(0.0).margin(1e-9));
    // poseAt agrees with the integration at the segment ends.
    const Pose2 e3 = p.segs[2].poseAt(p.segs[2].length);
    CHECK(e3.N == Approx(b.N).margin(1e-6));
    CHECK(e3.E == Approx(b.E).margin(1e-6));
    // The chosen word is the shortest VALID candidate, every valid
    // candidate reaches the goal, and none is shorter than the chord.
    words[p.word == "LSL" ? 0 : p.word == "RSR" ? 1 : p.word == "LSR" ? 2 : 3]++;
    for (const DubinsPath& c : dubinsCSCAll(a, b, R)) {
      if (!c.valid) continue;
      CHECK(c.length >= p.length - 1e-9);
      CHECK(c.length >= std::hypot(b.N - a.N, b.E - a.E) - 1e-9);
      const Pose2 ce = integratePath(c, a, R);
      CHECK(ce.N == Approx(b.N).margin(1e-6));
      CHECK(ce.E == Approx(b.E).margin(1e-6));
    }
  }
  // All four words appear over the sample.
  for (int w : words) CHECK(w > 0);
}

TEST_CASE("Dubins special cases: on-axis straight line and abeam parallel",
          "[guidance]") {
  const double R = 400.0;
  // Aligned poses on a line: pure straight leg.
  {
    const Pose2 a{0.0, 0.0, 30.0 * kDeg};
    const double d = 1500.0;
    const Pose2 b{d * std::cos(a.chi), d * std::sin(a.chi), a.chi};
    const DubinsPath p = dubinsCSC(a, b, R);
    REQUIRE(p.valid);
    CHECK(p.length == Approx(d).margin(1e-6));
    CHECK(p.segs[0].length == Approx(0.0).margin(1e-6));
    CHECK(p.segs[2].length == Approx(0.0).margin(1e-6));
  }
  // Parallel, opposite heading, abeam by exactly 2R: a half circle + line.
  {
    const Pose2 a{0.0, 0.0, 0.0};
    const double d = 1000.0;
    const Pose2 b{-d, 2.0 * R, M_PI};  // to the right (east), heading south
    const DubinsPath p = dubinsCSC(a, b, R);
    REQUIRE(p.valid);
    CHECK(p.length == Approx(M_PI * R + d).margin(1e-6));
    CHECK(p.word == "RSR");
  }
  // Goal directly behind, same line, opposite heading: LSL/RSR give a
  // 180-deg turn displaced by 2R and back: length pi R + ... > 0, and the
  // path exists (overhead entry case).
  {
    const Pose2 a{0.0, 0.0, 0.0};
    const Pose2 b{-2000.0, 0.0, M_PI};
    const DubinsPath p = dubinsCSC(a, b, R);
    REQUIRE(p.valid);
    const Pose2 end = integratePath(p, a, R);
    CHECK(end.N == Approx(b.N).margin(1e-6));
    CHECK(end.E == Approx(b.E).margin(1e-6));
  }
  // Coincident same-sense circles (goal on the start's own turning circle):
  // a single arc, no spurious 2 pi from the undefined tangent direction.
  {
    const Pose2 a{0.0, 0.0, 0.0};
    const double sw = 60.0 * kDeg;
    // Right circle centre (0, R); point at sweep sw: (R sin sw, R - R cos sw).
    const Pose2 b{R * std::sin(sw), R - R * std::cos(sw), sw};
    const DubinsPath p = dubinsCSC(a, b, R);
    REQUIRE(p.valid);
    CHECK(p.length == Approx(R * sw).margin(1e-6));
  }
  // The nearly-aligned-goal-just-ahead pathology: the same-sense words
  // LSL/RSR wrap a full circle (documented discontinuity); the inner-tangent
  // words exist only on the knife edge D >= 2R (here D - 2R ~ 0.1 m, which
  // is exactly why the planner has a direct-line candidate and rejects
  // near-full-circle arcs instead of relying on them).
  {
    const Pose2 a{0.0, 0.0, -0.8 * kDeg};
    const Pose2 b{25.0, 0.0, 0.0};
    for (const DubinsPath& c : dubinsCSCAll(a, b, R)) {
      if (!c.valid) continue;
      if (c.word == "LSL" || c.word == "RSR") CHECK(c.length > M_PI * R);
    }
    const DubinsPath dl = directLine(a, b);
    CHECK(dl.length == Approx(25.0).margin(1e-9));
  }
}

TEST_CASE("Path closest point and progress guard", "[guidance]") {
  Path path;
  PathSegment l1;
  l1.type = SegType::Line; l1.N0 = 0; l1.E0 = 0; l1.chi0 = 0; l1.length = 1000.0;
  path.append(l1);
  PathSegment arc;  // right turn from north to east, R = 300
  arc.type = SegType::ArcR; arc.N0 = 1000; arc.E0 = 0; arc.chi0 = 0; arc.R = 300;
  arc.Nc = 1000; arc.Ec = 300; arc.length = 300 * M_PI / 2;
  path.append(arc);
  // Point beside the line.
  int seg; double dist;
  double s = path.closest(400.0, 25.0, 0.0, 1e9, &seg, &dist);
  CHECK(s == Approx(400.0).margin(1e-9));
  CHECK(dist == Approx(25.0).margin(1e-9));
  CHECK(seg == 0);
  // Point outside the arc at 45 deg.
  const double ang = 45.0 * kDeg;
  const double Nq = 1000 + 350 * std::sin(ang), Eq = 300 - 350 * std::cos(ang);
  s = path.closest(Nq, Eq, 0.0, 1e9, &seg, &dist);
  CHECK(seg == 1);
  CHECK(s == Approx(1000.0 + 300.0 * ang).margin(1e-9));
  CHECK(dist == Approx(50.0).margin(1e-9));
  // Forward window: with s_max inside the line, the arc is ignored.
  s = path.closest(Nq, Eq, 0.0, 900.0, &seg, &dist);
  CHECK(seg == 0);
  CHECK(s <= 900.0);
  // Progress guard: with s_min past the line, the line is ignored.
  s = path.closest(400.0, 25.0, 1050.0, 1e9, &seg, &dist);
  CHECK(seg == 1);
  CHECK(s >= 1050.0);
  // poseAt continuity at the joint.
  const Pose2 pa = path.poseAt(1000.0 - 1e-9), pb = path.poseAt(1000.0 + 1e-9);
  CHECK(pa.N == Approx(pb.N).margin(1e-6));
  CHECK(pa.E == Approx(pb.E).margin(1e-6));
  // Far from the path (xte >> R): the closest point is well defined and the
  // window still bounds the progress.
  s = path.closest(-5000.0, 4000.0, 300.0, 800.0, &seg, &dist);
  CHECK(s >= 300.0);
  CHECK(s <= 800.0);
  CHECK(dist > 4000.0);
}

// A self-crossing pattern (a right-hand loop that comes back over its own
// entry leg) must not let the foot point jump onto the later leg: the
// window keeps it on the entry leg until progress has been made there.
TEST_CASE("Progress window on a self-crossing path", "[guidance]") {
  const double R = 300.0;
  Path path;
  PathSegment l1;  // north 600 m from the origin
  l1.type = SegType::Line; l1.N0 = 0; l1.E0 = 0; l1.chi0 = 0; l1.length = 600.0;
  path.append(l1);
  PathSegment loop;  // 300 deg right turn, comes back across the entry leg
  loop.type = SegType::ArcR; loop.N0 = 600; loop.E0 = 0; loop.chi0 = 0; loop.R = R;
  loop.Nc = 600; loop.Ec = R; loop.length = R * 300.0 * kDeg;
  path.append(loop);
  const Pose2 end = loop.poseAt(loop.length);
  PathSegment l2;
  l2.type = SegType::Line; l2.N0 = end.N; l2.E0 = end.E; l2.chi0 = end.chi; l2.length = 1500.0;
  path.append(l2);
  // The exit line crosses the entry leg; a point ON the entry leg near the
  // crossing is equidistant to both. With s_prog on the entry leg the window
  // must choose the entry leg.
  double cross_N = -1.0;
  for (double s = 0; s < l2.length; s += 0.5) {
    const Pose2 q = l2.poseAt(s);
    if (std::abs(q.E) < 0.5 && q.N > 0 && q.N < 600) { cross_N = q.N; break; }
  }
  REQUIRE(cross_N > 0.0);
  int seg; double dist;
  const double s = path.closest(cross_N, 0.0, 0.0, 500.0, &seg, &dist);
  CHECK(seg == 0);
  CHECK(s == Approx(cross_N).margin(1e-6));
  // Unwindowed, the exit leg is an equally close candidate (the ambiguity
  // the window exists to break).
  const double s_any = path.closest(cross_N, 0.5, 0.0, 1e9, &seg, &dist);
  CHECK(dist < 1.0);
  (void)s_any;
}

TEST_CASE("Pattern guidance: FAF fixed point, L1 sign, cone", "[guidance]") {
  Corridor cor;
  cor.length = 700.0; cor.width = 80.0; cor.heading = 0.0; cor.s_aim = -250.0;
  PatternGuidanceConfig cfg;
  cfg.V_app = 40.0;
  cfg.gamma_app = -4.0 * kDeg;

  SECTION("FAF placement: bounds, length matching, no bank limit in the plan") {
    for (double h0 : {100.0, 200.0, 300.0, 400.0}) {
      PatternGuidance g(cfg, cor);
      REQUIRE(g.plan(0.0, 0.0, M_PI, h0));  // overhead, opposite heading
      CHECK(g.dFinal() >= cfg.d_final_min);
      CHECK(g.dFinal() <= cfg.d_final_max);
      // The planned length matches the descent distance to within the
      // grid step whenever the scan is not clamped at a bound.
      const double descent = h0 / std::tan(4.0 * kDeg);
      const double planned = g.sAimPath();
      if (g.dFinal() > cfg.d_final_min + 1.0 && g.dFinal() < cfg.d_final_max - 1.0)
        CHECK(std::abs(planned - descent) <= 25.0);
      // No arc in the plan sweeps a near-full circle.
      for (const PathSegment& sg : g.path().segs())
        if (sg.type != SegType::Line) CHECK(sg.length <= sg.R * cfg.max_arc_sweep + 1e-9);
      // The FAF lies on the axis, upwind of the aim point, at corridor heading.
      const Pose2 faf = g.path().poseAt(g.sFafPath());
      double s, e;
      cor.frame(faf.N, faf.E, s, e);
      CHECK(e == Approx(0.0).margin(1e-6));
      CHECK(s == Approx(cor.s_aim - g.dFinal()).margin(1e-6));
      CHECK(wrap(faf.chi - cor.heading) == Approx(0.0).margin(1e-9));
    }
  }
  SECTION("on the final axis: zero bank, glideslope reference") {
    PatternGuidance g(cfg, cor);
    REQUIRE(g.plan(0.0, 0.0, M_PI, 250.0));
    const double s_f = g.sFafPath() + 200.0;
    const Pose2 p = g.path().poseAt(s_f);
    g.setProgress(s_f);  // teleport onto the final: re-center the window
    const double L = g.sAimPath() - s_f;
    const double h_gs = L * std::tan(4.0 * kDeg);
    // Exactly on the path, on speed, on the glideslope.
    GuidanceCmd c = g.step(p.N, p.E, h_gs, 40.0 * std::cos(p.chi), 40.0 * std::sin(p.chi));
    CHECK(c.on_final);
    CHECK(c.phi_cmd == Approx(0.0).margin(1e-9));
    CHECK(c.gamma_ref == Approx(cfg.gamma_app).margin(1e-12));
    CHECK(c.L_to_go == Approx(L).margin(1e-6));
    // Offset to the RIGHT of the axis -> bank LEFT (negative).
    GuidanceCmd r = g.step(p.N, p.E + 30.0, h_gs, 40.0, 0.0);
    CHECK(r.phi_cmd < 0.0);
    GuidanceCmd l = g.step(p.N, p.E - 30.0, h_gs, 40.0, 0.0);
    CHECK(l.phi_cmd > 0.0);
    // High -> steeper; low -> shallower; both clamped.
    GuidanceCmd hi = g.step(p.N, p.E, h_gs + 40.0, 40.0, 0.0);
    CHECK(hi.gamma_ref < cfg.gamma_app);
    CHECK(hi.gamma_ref >= cfg.gamma_min);
    GuidanceCmd lo = g.step(p.N, p.E, std::max(0.0, h_gs - 40.0), 40.0, 0.0);
    CHECK(lo.gamma_ref > cfg.gamma_app);
    GuidanceCmd vhi = g.step(p.N, p.E, h_gs + 5000.0, 40.0, 0.0);
    CHECK(vhi.gamma_ref == Approx(cfg.gamma_min));
  }
  SECTION("on-axis start plans a direct line, robust to small perturbations") {
    // The straight-in anchor: 3 km out on the axis. The trim sideslip and
    // small placement errors must not flip the plan into a looping CSC word
    // (the Dubins discontinuity for a nearly aligned goal just ahead).
    std::mt19937 rng(3);
    std::uniform_real_distribution<double> un(-1.0, 1.0);
    for (int k = 0; k < 25; ++k) {
      const double N0 = -3200.0 + 5.0 * un(rng), E0 = 5.0 * un(rng);
      const double chi0 = 1.0 * kDeg * un(rng);
      const double h0 = 3000.0 * std::tan(4.0 * kDeg);
      PatternGuidance g(cfg, cor);
      REQUIRE(g.plan(N0, E0, chi0, h0));
      CHECK(g.word() == "DIRECT");
      // Path to the aim point = the straight distance to it (aim at s_aim).
      CHECK(g.sAimPath() == Approx(cor.s_aim - N0).margin(30.0));
    }
    // Off the axis by more than the tolerance: a CSC word, still no loop.
    PatternGuidance g(cfg, cor);
    REQUIRE(g.plan(-3200.0, 200.0, 0.0, 3000.0 * std::tan(4.0 * kDeg)));
    CHECK(g.word() != "DIRECT");
    CHECK(g.sAimPath() < 3000.0 + 2.0 * M_PI * g.turnRadius());
    for (const PathSegment& sg : g.path().segs())
      if (sg.type != SegType::Line) CHECK(sg.length <= sg.R * cfg.max_arc_sweep + 1e-9);
  }
  SECTION("non-descending gamma_app is rejected") {
    PatternGuidanceConfig bad = cfg;
    bad.gamma_app = 0.0;
    PatternGuidance g(bad, cor);
    CHECK_FALSE(g.plan(0.0, 0.0, 0.0, 100.0));
  }
  SECTION("overhead start: first command turns, L1 point lies ahead on path") {
    PatternGuidance g(cfg, cor);
    REQUIRE(g.plan(0.0, 0.0, M_PI, 250.0));
    GuidanceCmd c = g.step(0.0, 0.0, 250.0, -40.0, 0.0);
    CHECK_FALSE(c.on_final);
    CHECK(std::abs(c.phi_cmd) > 5.0 * kDeg);
    CHECK(std::abs(c.phi_cmd) <= cfg.phi_max + 1e-12);
    CHECK(c.s_path == Approx(0.0).margin(1e-9));
  }
}
