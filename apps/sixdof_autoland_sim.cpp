#include <cstdlib>
#include <iostream>
#include <string>

#include "autoland/sixdof_sim.hpp"

// 6-DOF straight-in water-landing simulation (nonlinear plant, cascaded-PID
// nominal, plant-side wind gust + surface waves; no CBF filter).
//
// Usage:
//   sixdof_autoland_sim [stab] [aircraft.yaml] [scenario.yaml] [out.csv]
//
// All arguments are optional. The default scenario flies the DHC-2 Beaver
// plant (data/beaver_landing_calm.yaml); scenarios with `plant: vspaero`
// (data/sixdof_*.yaml) use the AHAB deck instead, for which the stab and
// aircraft.yaml paths apply (the defaults bundle the beta-symmetric sweep so
// crosswind sideslip never leaves the table).
int main(int argc, char** argv) {
  const std::string data = AUTOLAND_DATA_DIR;
  std::string stab = data + "/AHAB_combined_betasym.stab";
  std::string aircraft = data + "/aircraft.yaml";
  std::string scenario = data + "/beaver_landing_calm.yaml";
  std::string out = "sixdof_autoland_log.csv";

  if (argc > 1) stab = argv[1];
  if (argc > 2) aircraft = argv[2];
  if (argc > 3) scenario = argv[3];
  if (argc > 4) out = argv[4];

  try {
    autoland::SixDofSim sim(stab, aircraft, scenario);
    autoland::SixDofTouchdown td = sim.run(out);
    return td.reached ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
}
