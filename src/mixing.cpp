#include "autoland/mixing.hpp"

#include <algorithm>
#include <cctype>
#include <iostream>

namespace autoland {
namespace {

std::string lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return s;
}

}  // namespace

Mixing Mixing::build(const AircraftConfig& cfg, const AeroTable& table) {
  Mixing m;
  const int ng = table.numControlGroups();

  // 1) Explicit matrix from config wins.
  if (cfg.mixing.size() > 0) {
    if (cfg.mixing.rows() != ng || cfg.mixing.cols() != 3)
      throw std::runtime_error(
          "mixing: config matrix is " + std::to_string(cfg.mixing.rows()) +
          "x" + std::to_string(cfg.mixing.cols()) + " but the .stab has " +
          std::to_string(ng) + " control groups (need ng x 3).");
    m.M_ = cfg.mixing;
    return m;
  }

  // 2) Name-based default. Columns are [delta_e, delta_a, delta_r].
  //
  // POLARITY: the virtual controls are defined so a POSITIVE command produces a
  // POSITIVE body-axis response in standard (z-down) axes:
  //     +delta_e -> nose-UP pitch (Cm > 0)
  //     +delta_a -> right roll    (Cl > 0)
  //     +delta_r -> nose-right yaw (Cn > 0)
  // For THIS .stab (after the z-up -> z-down frame transform in dynamics) the
  // Elevator and Aileron groups produce a NEGATIVE moment per positive group
  // deflection, while Rudder produces a positive one. The default map therefore
  // uses -1 for elevator/aileron and +1 for rudder so the controller's
  // natural-sign PID loops are stable. TODO: confirm against your OpenVSP setup.
  m.M_ = Mat::Zero(ng, 3);
  const auto& names = table.controlGroupNames();
  bool all_mapped = true;
  for (int i = 0; i < ng; ++i) {
    const std::string n = lower(names[i]);
    if (n.find("eleva") != std::string::npos) {
      m.M_(i, 0) = -1.0;  // delta_e column  (+delta_e -> nose up)
    } else if (n.find("aile") != std::string::npos) {
      m.M_(i, 1) = -1.0;  // delta_a column  (+delta_a -> right roll)
    } else if (n.find("rudder") != std::string::npos) {
      m.M_(i, 2) = 1.0;  // delta_r column  (+delta_r -> nose right)
    } else {
      all_mapped = false;
      std::cerr << "[mixing] WARNING: control group '" << names[i]
                << "' not recognized by the default name map; its row is zero. "
                << "Set 'mixing.matrix' in aircraft.yaml (TODO)." << std::endl;
    }
  }
  if (!all_mapped)
    std::cerr << "[mixing] NOTE: default mixing is incomplete; results are only "
                 "valid once the mixing map is confirmed in config." << std::endl;
  return m;
}

Vec Mixing::apply(double delta_e, double delta_a, double delta_r) const {
  Eigen::Vector3d v(delta_e, delta_a, delta_r);
  return M_ * v;
}

}  // namespace autoland
