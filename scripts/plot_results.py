#!/usr/bin/env python3
"""Plot an autoland_sim CSV log.

Standalone helper -- intentionally NOT part of the CMake build. Requires only
matplotlib + numpy:

    python3 scripts/plot_results.py autoland_log.csv
    python3 scripts/plot_results.py autoland_log.csv --save run.png
"""
import argparse
import csv
import sys

import matplotlib.pyplot as plt
import numpy as np


def load_csv(path):
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    if not rows:
        sys.exit(f"no data in {path}")
    cols = {k: np.array([float(r[k]) for r in rows]) for k in rows[0].keys()}
    return cols


def main():
    ap = argparse.ArgumentParser(description="Plot autoland sim log")
    ap.add_argument("csv", help="path to autoland_sim CSV log")
    ap.add_argument("--save", help="save figure to this path instead of showing")
    args = ap.parse_args()

    d = load_csv(args.csv)
    t = d["t"]

    fig, ax = plt.subplots(3, 2, figsize=(13, 10))
    fig.suptitle("Autoland closed-loop response", fontsize=14)

    # Altitude + glideslope reference
    ax[0, 0].plot(t, d["h"], label="h")
    ax[0, 0].plot(t, d["h_ref"], "--", label="h_ref (glideslope)")
    ax[0, 0].axhline(0.0, color="k", lw=0.6)
    ax[0, 0].set_ylabel("altitude [m]")
    ax[0, 0].legend(); ax[0, 0].grid(True)

    # Airspeed
    ax[0, 1].plot(t, d["V"], label="V")
    ax[0, 1].plot(t, d["V_cmd"], "--", label="V_cmd")
    ax[0, 1].set_ylabel("airspeed [m/s]")
    ax[0, 1].legend(); ax[0, 1].grid(True)

    # Pitch attitude vs command
    ax[1, 0].plot(t, d["theta_deg"], label="theta")
    ax[1, 0].plot(t, d["theta_cmd_deg"], "--", label="theta_cmd")
    ax[1, 0].plot(t, d["alpha_deg"], ":", label="alpha")
    ax[1, 0].set_ylabel("pitch [deg]")
    ax[1, 0].legend(); ax[1, 0].grid(True)

    # Sink rate + flare command
    ax[1, 1].plot(t, d["sink"], label="sink rate")
    ax[1, 1].plot(t, d["w_cmd"], "--", label="w_cmd")
    if "flaring" in d:
        flare_on = t[d["flaring"] > 0.5]
        if flare_on.size:
            ax[1, 1].axvspan(flare_on.min(), flare_on.max(), color="orange",
                             alpha=0.15, label="flare")
    ax[1, 1].set_ylabel("sink [m/s, +down]")
    ax[1, 1].legend(); ax[1, 1].grid(True)

    # Lateral: cross-track + bank command
    ax[2, 0].plot(t, d["y"], label="cross-track y")
    ax[2, 0].plot(t, d["phi_deg"], "--", label="phi [deg]")
    ax[2, 0].plot(t, d["phi_cmd_deg"], ":", label="phi_cmd [deg]")
    ax[2, 0].plot(t, d["beta_deg"], "-.", label="beta [deg]")
    ax[2, 0].axhline(0.0, color="k", lw=0.6)
    ax[2, 0].set_xlabel("time [s]"); ax[2, 0].set_ylabel("lateral")
    ax[2, 0].legend(); ax[2, 0].grid(True)

    # Controls
    ax[2, 1].plot(t, d["de_deg"], label="delta_e [deg]")
    ax[2, 1].plot(t, d["da_deg"], label="delta_a [deg]")
    ax[2, 1].plot(t, d["dr_deg"], label="delta_r [deg]")
    ax[2, 1].plot(t, d["dT"] * 10.0, "--", label="throttle x10")
    ax[2, 1].set_xlabel("time [s]"); ax[2, 1].set_ylabel("controls")
    ax[2, 1].legend(); ax[2, 1].grid(True)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    if args.save:
        fig.savefig(args.save, dpi=120)
        print(f"saved {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
