#include <cstdlib>
#include <iostream>
#include <string>

#include "autoland/lon_sim.hpp"

// Augmented-longitudinal CBF-QP water-landing simulation.
//
// Usage:
//   lon_autoland_sim [stab] [aircraft.yaml] [lon_scenario.yaml] [out.csv]
//
// All arguments are optional and default to the bundled data/ files.
int main(int argc, char** argv) {
  const std::string data = AUTOLAND_DATA_DIR;
  std::string stab = data + "/AHAB_combined.stab";  // real vehicle aero deck
  std::string aircraft = data + "/aircraft.yaml";
  std::string scenario = data + "/lon_scenario.yaml";
  std::string out = "lon_autoland_log.csv";

  if (argc > 1) stab = argv[1];
  if (argc > 2) aircraft = argv[2];
  if (argc > 3) scenario = argv[3];
  if (argc > 4) out = argv[4];

  try {
    autoland::LonSim sim(stab, aircraft, scenario);
    autoland::LonTouchdown td = sim.run(out);
    return td.reached ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
}
