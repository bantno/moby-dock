#include <cstdlib>
#include <iostream>
#include <string>

#include "autoland/sim.hpp"

// Standalone autoland simulation.
//
// Usage:
//   autoland_sim [stab] [aircraft.yaml] [scenario.yaml] [out.csv]
//
// All arguments are optional and default to the bundled data/ files and
// ./autoland_log.csv. Plot the result with scripts/plot_results.py.
int main(int argc, char** argv) {
  const std::string data = AUTOLAND_DATA_DIR;
  std::string stab = data + "/AHAB_sweep.stab";
  std::string aircraft = data + "/aircraft.yaml";
  std::string scenario = data + "/scenario.yaml";
  std::string out = "autoland_log.csv";

  if (argc > 1) stab = argv[1];
  if (argc > 2) aircraft = argv[2];
  if (argc > 3) scenario = argv[3];
  if (argc > 4) out = argv[4];

  try {
    autoland::Sim sim(stab, aircraft, scenario);
    autoland::TouchdownSummary td = sim.run(out);
    return td.reached ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
}
