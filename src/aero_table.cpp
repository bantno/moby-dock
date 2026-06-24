#include "autoland/aero_table.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace autoland {
namespace {

constexpr double kDeg2Rad = M_PI / 180.0;
constexpr double kGridTol = 1e-6;  // tolerance for matching axis grid values

std::vector<std::string> tokenize(const std::string& line) {
  std::vector<std::string> out;
  std::istringstream ss(line);
  std::string tok;
  while (ss >> tok) out.push_back(tok);
  return out;
}

bool startsWith(const std::string& s, const std::string& p) {
  return s.size() >= p.size() && s.compare(0, p.size(), p) == 0;
}

// Map a coefficient row label ("CFx".."CMz") to its Coef index, or -1 if it is
// a row we don't consume (CL, CD, CS, CMl, CMm, CMn, ...).
int coefIndex(const std::string& label) {
  if (label == "CFx") return CFX;
  if (label == "CFy") return CFY;
  if (label == "CFz") return CFZ;
  if (label == "CMx") return CMX;
  if (label == "CMy") return CMY;
  if (label == "CMz") return CMZ;
  return -1;
}

// Insert v into a sorted-unique axis vector (within tolerance); return nothing.
void insertAxisValue(std::vector<double>& axis, double v) {
  for (double a : axis)
    if (std::abs(a - v) < kGridTol) return;
  axis.push_back(v);
}

int axisIndex(const std::vector<double>& axis, double v) {
  for (int i = 0; i < static_cast<int>(axis.size()); ++i)
    if (std::abs(axis[i] - v) < kGridTol) return i;
  return -1;
}

// One parsed block, before it is placed into the grid.
struct Block {
  double alpha{0}, beta{0}, mach{0};  // alpha/beta in radians
  AeroDerivs derivs;
};

}  // namespace

AeroTable AeroTable::fromFile(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("aero_table: cannot open file: " + path);

  std::vector<std::string> lines;
  std::string line;
  while (std::getline(in, line)) lines.push_back(line);

  AeroTable table;
  std::vector<Block> blocks;
  std::vector<std::string> cg_names;  // captured once from the first block

  // Known finite-difference perturbation row labels in the "Case" section; any
  // trailing rows beyond these are the named control groups (in group order).
  const std::vector<std::string> known_case_rows = {
      "Base_Aero", "Alpha", "Beta", "Roll__Rate",
      "Pitch_Rate", "Yaw___Rate", "Mach"};
  auto isKnownCaseRow = [&](const std::string& l) {
    return std::find(known_case_rows.begin(), known_case_rows.end(), l) !=
           known_case_rows.end();
  };

  // State machine over the file.
  bool in_block = false;
  Block cur;
  bool have_ref = false;
  enum Section { NONE, CASE, DERIV } section = NONE;
  int n_ctrl = -1;  // number of control-group derivative columns

  auto flushBlock = [&]() {
    if (in_block && have_ref) blocks.push_back(cur);
    cur = Block{};
    have_ref = false;
    section = NONE;
  };

  for (const std::string& raw : lines) {
    // Block separator.
    if (raw.find("****") != std::string::npos) {
      flushBlock();
      in_block = true;
      continue;
    }
    if (!in_block) continue;

    std::string trimmed = raw;
    // strip leading spaces for prefix checks
    size_t first = trimmed.find_first_not_of(" \t");
    std::string body = (first == std::string::npos) ? "" : trimmed.substr(first);

    // ---- reference quantities: "Name   value   units" --------------------
    auto refVal = [&](const char* name, double& dst, double scale = 1.0) -> bool {
      if (startsWith(body, name)) {
        auto t = tokenize(body);
        if (t.size() >= 2) {
          dst = std::stod(t[1]) * scale;
          have_ref = true;
          return true;
        }
      }
      return false;
    };
    if (refVal("Sref_", table.sref_)) continue;
    if (refVal("Cref_", table.cref_)) continue;
    if (refVal("Bref_", table.bref_)) continue;
    if (refVal("Xcg_", table.xcg_)) continue;
    if (refVal("Ycg_", table.ycg_)) continue;
    if (refVal("Zcg_", table.zcg_)) continue;
    if (refVal("Mach_", cur.mach)) continue;
    if (refVal("AoA_", cur.alpha, kDeg2Rad)) continue;   // deg -> rad
    if (refVal("Beta_", cur.beta, kDeg2Rad)) continue;   // deg -> rad

    // ---- section headers ---------------------------------------------------
    if (startsWith(body, "Case")) {
      section = CASE;
      continue;
    }
    if (startsWith(body, "Coef")) {
      // Header: Coef Total Alpha Beta p q r Mach U ConGrp_1 ... ConGrp_N
      auto t = tokenize(body);
      // columns after the label: Total + 7 std derivs + control groups
      n_ctrl = static_cast<int>(t.size()) - 1 - 8;
      if (n_ctrl < 0) n_ctrl = 0;
      section = DERIV;
      continue;
    }

    if (body.empty() || body[0] == '#') continue;

    auto t = tokenize(body);
    if (t.empty()) continue;

    if (section == CASE) {
      // Capture control-group names (rows after the known perturbation rows).
      if (!isKnownCaseRow(t[0]) && cg_names.size() < known_case_rows.size() + 16) {
        // Only the leading alpha token must be a label; guard against stray rows.
        if (coefIndex(t[0]) == -1)  // not a coefficient row
          cg_names.push_back(t[0]);
      }
      continue;
    }

    if (section == DERIV) {
      int ci = coefIndex(t[0]);
      if (ci < 0) continue;  // skip CL/CD/CS/CMl/CMm/CMn rows
      // Expect: label Total Alpha Beta p q r Mach U [ctrl...]
      if (static_cast<int>(t.size()) < 9)
        throw std::runtime_error("aero_table: short derivative row: " + body);
      AeroDerivs& d = cur.derivs;
      d.base[ci]    = std::stod(t[1]);
      d.d_alpha[ci] = std::stod(t[2]);
      d.d_beta[ci]  = std::stod(t[3]);
      d.d_p[ci]     = std::stod(t[4]);
      d.d_q[ci]     = std::stod(t[5]);
      d.d_r[ci]     = std::stod(t[6]);
      d.d_mach[ci]  = std::stod(t[7]);
      d.d_u[ci]     = std::stod(t[8]);
      const int avail = static_cast<int>(t.size()) - 9;
      const int ncg = (n_ctrl >= 0) ? std::min(n_ctrl, avail) : avail;
      if (d.d_ctrl.empty()) d.d_ctrl.assign(ncg, CoefVec::Zero());
      for (int g = 0; g < ncg; ++g) d.d_ctrl[g][ci] = std::stod(t[9 + g]);
      continue;
    }
  }
  flushBlock();

  if (blocks.empty())
    throw std::runtime_error("aero_table: no data blocks parsed from " + path);

  // Control-group names: the trailing case-row labels captured per block repeat
  // every block; keep the first full set.
  if (cg_names.size() > 0) {
    // de-dup while preserving order, keeping only the distinct leading run
    std::vector<std::string> uniq;
    for (const auto& n : cg_names) {
      if (std::find(uniq.begin(), uniq.end(), n) == uniq.end())
        uniq.push_back(n);
      else
        break;  // names start repeating on the next block
    }
    table.cg_names_ = uniq;
  }
  // Fall back to generic names if the Case section was not present.
  if (table.cg_names_.empty() && !blocks.empty()) {
    for (size_t g = 0; g < blocks.front().derivs.d_ctrl.size(); ++g)
      table.cg_names_.push_back("ConGrp_" + std::to_string(g + 1));
  }

  // ---- build axes -----------------------------------------------------------
  for (const auto& b : blocks) {
    insertAxisValue(table.alphas_, b.alpha);
    insertAxisValue(table.betas_, b.beta);
    insertAxisValue(table.machs_, b.mach);
  }
  std::sort(table.alphas_.begin(), table.alphas_.end());
  std::sort(table.betas_.begin(), table.betas_.end());
  std::sort(table.machs_.begin(), table.machs_.end());

  const int na = static_cast<int>(table.alphas_.size());
  const int nb = static_cast<int>(table.betas_.size());
  const int nm = static_cast<int>(table.machs_.size());
  table.nodes_.assign(static_cast<size_t>(na) * nb * nm, AeroDerivs{});
  std::vector<char> filled(table.nodes_.size(), 0);

  for (const auto& b : blocks) {
    int ia = axisIndex(table.alphas_, b.alpha);
    int ib = axisIndex(table.betas_, b.beta);
    int im = axisIndex(table.machs_, b.mach);
    int k = table.idx(ia, ib, im);
    table.nodes_[k] = b.derivs;
    filled[k] = 1;
  }
  for (size_t k = 0; k < filled.size(); ++k)
    if (!filled[k])
      throw std::runtime_error(
          "aero_table: sweep grid is not complete (missing nodes); the .stab "
          "must be a full alpha x beta x Mach grid");

  return table;
}

AeroLookup AeroTable::lookup(double alpha, double beta, double mach) const {
  AeroLookup out;

  // Per-axis: clamp, find bracketing indices i0/i1 and blend fraction f.
  auto bracket = [](const std::vector<double>& axis, double v, int& i0, int& i1,
                    double& f, bool& clamped) {
    const int n = static_cast<int>(axis.size());
    if (n == 1) { i0 = i1 = 0; f = 0.0; return; }
    if (v <= axis.front()) {
      if (v < axis.front() - kGridTol) clamped = true;
      i0 = i1 = 0; f = 0.0; return;
    }
    if (v >= axis.back()) {
      if (v > axis.back() + kGridTol) clamped = true;
      i0 = i1 = n - 1; f = 0.0; return;
    }
    for (int i = 0; i < n - 1; ++i) {
      if (v >= axis[i] && v <= axis[i + 1]) {
        i0 = i; i1 = i + 1;
        f = (v - axis[i]) / (axis[i + 1] - axis[i]);
        return;
      }
    }
    i0 = i1 = n - 1; f = 0.0;
  };

  int a0, a1, b0, b1, m0, m1;
  double fa, fb, fm;
  bool clamped = false;
  bracket(alphas_, alpha, a0, a1, fa, clamped);
  bracket(betas_, beta, b0, b1, fb, clamped);
  bracket(machs_, mach, m0, m1, fm, clamped);
  out.clamped = clamped;

  if (clamped) {
    static int warn_count = 0;
    if (warn_count < 20) {
      std::cerr << "[aero_table] WARNING: query (alpha=" << alpha
                << " rad, beta=" << beta << " rad, Mach=" << mach
                << ") off-grid; clamped (no extrapolation across swept axes)."
                << std::endl;
      if (++warn_count == 20)
        std::cerr << "[aero_table] (further off-grid warnings suppressed)"
                  << std::endl;
    }
  }

  // Trilinear blend over the (up to) 8 corners. Lambda blends one field.
  auto blendField = [&](CoefVec AeroDerivs::*field) -> CoefVec {
    auto N = [&](int ia, int ib, int im) -> const CoefVec& {
      return nodes_[idx(ia, ib, im)].*field;
    };
    CoefVec c00 = (1 - fm) * N(a0, b0, m0) + fm * N(a0, b0, m1);
    CoefVec c01 = (1 - fm) * N(a0, b1, m0) + fm * N(a0, b1, m1);
    CoefVec c10 = (1 - fm) * N(a1, b0, m0) + fm * N(a1, b0, m1);
    CoefVec c11 = (1 - fm) * N(a1, b1, m0) + fm * N(a1, b1, m1);
    CoefVec c0 = (1 - fb) * c00 + fb * c01;
    CoefVec c1 = (1 - fb) * c10 + fb * c11;
    return (1 - fa) * c0 + fa * c1;
  };

  out.d.base    = blendField(&AeroDerivs::base);
  out.d.d_alpha = blendField(&AeroDerivs::d_alpha);
  out.d.d_beta  = blendField(&AeroDerivs::d_beta);
  out.d.d_p     = blendField(&AeroDerivs::d_p);
  out.d.d_q     = blendField(&AeroDerivs::d_q);
  out.d.d_r     = blendField(&AeroDerivs::d_r);
  out.d.d_mach  = blendField(&AeroDerivs::d_mach);
  out.d.d_u     = blendField(&AeroDerivs::d_u);

  const int ncg = numControlGroups();
  out.d.d_ctrl.assign(ncg, CoefVec::Zero());
  for (int g = 0; g < ncg; ++g) {
    auto Ng = [&](int ia, int ib, int im) -> const CoefVec& {
      return nodes_[idx(ia, ib, im)].d_ctrl[g];
    };
    CoefVec c00 = (1 - fm) * Ng(a0, b0, m0) + fm * Ng(a0, b0, m1);
    CoefVec c01 = (1 - fm) * Ng(a0, b1, m0) + fm * Ng(a0, b1, m1);
    CoefVec c10 = (1 - fm) * Ng(a1, b0, m0) + fm * Ng(a1, b0, m1);
    CoefVec c11 = (1 - fm) * Ng(a1, b1, m0) + fm * Ng(a1, b1, m1);
    CoefVec c0 = (1 - fb) * c00 + fb * c01;
    CoefVec c1 = (1 - fb) * c10 + fb * c11;
    out.d.d_ctrl[g] = (1 - fa) * c0 + fa * c1;
  }

  // Reference point for the derivative buildup (see header): for a swept axis
  // use the query value so the residual term vanishes (base interpolation
  // already carries the dependence); for a degenerate single-node axis use the
  // node value so the derivative gives a local linear correction.
  out.alpha_ref = (alphas_.size() > 1) ? alpha : alphas_[0];
  out.beta_ref  = (betas_.size() > 1) ? beta : betas_[0];
  out.mach_ref  = (machs_.size() > 1) ? mach : machs_[0];
  return out;
}

}  // namespace autoland
