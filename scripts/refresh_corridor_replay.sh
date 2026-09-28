#!/usr/bin/env bash
# Rebuild figures/corridor_replay.html from the corridor scenarios: runs every
# data/beaver_corridor_*.yaml through the 6-DOF sim into runs/ and bundles the
# traces into one self-contained HTML page (scripts/build_corridor_replay.py).
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p runs
args=()
for y in data/beaver_corridor_straight.yaml data/beaver_corridor_overhead.yaml \
         data/beaver_corridor_crosswind.yaml data/beaver_corridor_fail_*.yaml data/beaver_corridor_cbf_*.yaml; do
  name=$(basename "$y" .yaml); name=${name#beaver_corridor_}; name=${name/fail_/fail:}; name=${name/cbf_/cbf:}
  csv="runs/corridor_${name/:/_}.csv"
  ./build/sixdof_autoland_sim "" "" "$y" "$csv" > /dev/null
  args+=("$name=$csv:$y")
done
python3 scripts/build_corridor_replay.py figures/corridor_replay.html "${args[@]}"
