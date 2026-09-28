#!/usr/bin/env python3
"""Build a self-contained interactive replay page for corridor landing runs.

Usage: build_corridor_replay.py OUT.html NAME=RUN.csv:SCENARIO.yaml [...]
Embeds downsampled traces (every `--stride` samples), the terrain bump list
and the corridor from each scenario YAML into ONE plain HTML file (no server,
no build step, no external assets except the optional Google Fonts link,
which falls back to the system fonts offline): a top-down map over the
terrain with the planned/flown/rollout tracks and a scrubbable aircraft, plus
time-synced strips (altitude vs cone, bank, cross-track, speed). Open the
output in any browser (file:// is fine). The committed instance is
figures/corridor_replay.html; regenerate it with scripts/refresh_corridor_replay.sh.
"""
import argparse
import json
import sys

import numpy as np
import yaml

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from plot_corridor_landing import terrain_from_yaml  # noqa: E402

DEG = np.pi / 180.0
COLS = ["t", "x", "y", "h", "psi_deg", "phi_deg", "phi_cmd_deg", "beta_deg", "V_air",
        "gamma_deg", "chi_cmd_deg", "N_ref", "E_ref", "h_gs", "L_to_go", "s_corr",
        "e_corr", "hT", "Vg", "xte", "phase", "seg", "sink", "theta_deg", "alpha_deg",
        "cbf6_active", "cbf6_be", "h6_alpha", "h6_phi", "h6_beta", "de_nom", "de"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("runs", nargs="+", help="NAME=RUN.csv:SCENARIO.yaml")
    ap.add_argument("--stride", type=int, default=10)
    ap.add_argument("--template", default=__file__.rsplit("/", 1)[0] + "/corridor_replay_template.html")
    a = ap.parse_args()
    runs = []
    for spec in a.runs:
        name, rest = spec.split("=", 1)
        csv_path, yaml_path = rest.split(":", 1)
        d = np.genfromtxt(csv_path, delimiter=",", names=True)
        with open(yaml_path) as f:
            sc = yaml.safe_load(f)
        w = sc["world"]
        bumps, base = terrain_from_yaml(w)
        # keep every stride-th airborne row, every stride-th rollout row, the
        # touchdown row and the last row.
        idx = list(range(0, len(d), a.stride))
        td_i = int(np.where(d["phase"] == 0)[0][-1])
        for k in (td_i, len(d) - 1):
            if k not in idx:
                idx.append(k)
        idx.sort()
        cols = {c: [round(float(v), 3) for v in d[c][idx]] for c in COLS}
        c = w["corridor"]
        runs.append(dict(
            name=name, csv=csv_path, scenario=yaml_path,
            corridor=dict(N=float(c.get("center_N", 0)), E=float(c.get("center_E", 0)),
                          length=float(c["length"]), width=float(c["width"]),
                          heading=float(c.get("heading_deg", 0)) * DEG,
                          aim_s=float(c.get("aim_s", -float(c["length"]) / 2 + 100))),
            terrain=dict(base=base, bumps=[[b["N"], b["E"], b["A"], b["sa"], b["sb"], b["rot"], b["order"]] for b in bumps]),
            gamma_app=float(sc.get("gamma_app_deg", -3.5)), V_app=float(sc.get("V_app", 40.0)),
            h0=float(sc["initial"]["h0"]), wind=sc.get("wind", {}) or {},
            cbf=sc.get("sixdof_cbf", {}) or {},
            data=cols, td_index=idx.index(td_i)))
    tpl = open(a.template).read()
    body = tpl.replace("/*__RUNS__*/", json.dumps(runs, separators=(",", ":")))
    html = ("<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
            "<meta name=\"color-scheme\" content=\"light dark\">\n" + body +
            "\n</body>\n</html>\n")
    open(a.out, "w").write(html)
    print(f"wrote {a.out} ({len(html) / 1e6:.2f} MB, {len(runs)} runs)")


if __name__ == "__main__":
    main()
