#pragma once
#include "autoland/aero_table.hpp"
#include "autoland/config.hpp"
#include "autoland/types.hpp"

// =============================================================================
// Control mixing: virtual controls [delta_e, delta_a, delta_r] -> physical
// VSPAero control-group deflections (group order matches the .stab).
//
// For the AHAB model OpenVSP already exposes independent "Elevator" and
// "Rudder" groups (the V-tail ruddervator mixing was defined inside OpenVSP),
// plus "Ailerons", so the default map is effectively identity:
//     Ailerons <- delta_a,  Elevator <- delta_e,  Rudder <- delta_r.
//
// If a future .stab instead exposes raw ruddervator halves, set the matrix in
// aircraft.yaml to the standard V-tail map (symmetric = elevator, anti = rudder):
//     ruddervator_L <-  delta_e + delta_r
//     ruddervator_R <-  delta_e - delta_r
//
// The matrix lives in config because it depends entirely on how the control
// groups were defined in OpenVSP. TODO: confirm for your model.
// =============================================================================
namespace autoland {

class Mixing {
 public:
  // Build from config; if config.mixing is empty, construct a name-based
  // default from the .stab control-group names.
  static Mixing build(const AircraftConfig& cfg, const AeroTable& table);

  // Map virtual controls (radians) to a deflection per control group (radians).
  Vec apply(double delta_e, double delta_a, double delta_r) const;

  const Mat& matrix() const { return M_; }
  int numGroups() const { return static_cast<int>(M_.rows()); }

 private:
  Mat M_;  // numGroups x 3, columns ordered [delta_e, delta_a, delta_r]
};

}  // namespace autoland
