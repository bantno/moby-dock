#!/usr/bin/env python3
"""Plot the CBF internals from one autoland_sim CSV log.

Shows what the safety filter is doing: descent rate, pitch, the nominal vs the
CBF-filtered (safe) control on the longitudinal channels {delta_e, delta_T}, and
the two barrier function values (airspeed and descent-rate) with the b=0 safety
boundary.

    python3 scripts/plot_cbf_detail.py autoland_log.csv --save figures/cbf_detail.png
"""
import argparse
import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        sys.exit(f"no data in {path}")
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0].keys()}


def main():
    ap = argparse.ArgumentParser(description="Plot CBF internals from a sim log")
    ap.add_argument("csv", help="path to autoland_sim CSV log")
    ap.add_argument("--vsafe", type=float, default=0.6,
                    help="hull-safe touchdown sink rate [m/s] (reference line)")
    ap.add_argument("--save", help="save figure to this path instead of showing")
    args = ap.parse_args()

    d = load(args.csv)
    need = ["de_nom_deg", "dT_nom", "b_airspeed", "b_descent"]
    missing = [c for c in need if c not in d]
    if missing:
        sys.exit(f"{args.csv} is missing columns {missing}; re-run autoland_sim "
                 "after the logging update.")
    t = d["t"]

    fig, ax = plt.subplots(3, 2, figsize=(13, 11))
    fig.suptitle("CBF safety filter — internals", fontsize=14)

    # Descent rate (sink, +down).
    ax[0, 0].plot(t, d["sink"], color="C0", label="sink rate (+down)")
    ax[0, 0].axhline(args.vsafe, color="green", ls="--", lw=1,
                     label=f"v_safe = {args.vsafe}")
    ax[0, 0].axhline(0.0, color="k", lw=0.6)
    ax[0, 0].set_ylabel("descent rate [m/s]"); ax[0, 0].set_title("Descent rate")
    ax[0, 0].legend(); ax[0, 0].grid(True)

    # Pitch (and alpha for context).
    ax[0, 1].plot(t, d["theta_deg"], color="C1", label="theta (pitch)")
    if "alpha_deg" in d:
        ax[0, 1].plot(t, d["alpha_deg"], color="C2", ls=":", label="alpha")
    ax[0, 1].axhline(0.0, color="k", lw=0.6)
    ax[0, 1].set_ylabel("angle [deg]"); ax[0, 1].set_title("Pitch attitude")
    ax[0, 1].legend(); ax[0, 1].grid(True)

    # Elevator: nominal vs safe.
    ax[1, 0].plot(t, d["de_nom_deg"], color="C3", ls="--", label="delta_e nominal")
    ax[1, 0].plot(t, d["de_deg"], color="C0", label="delta_e safe (CBF)")
    ax[1, 0].set_ylabel("elevator [deg]")
    ax[1, 0].set_title("Elevator: nominal vs safe")
    ax[1, 0].legend(); ax[1, 0].grid(True)

    # Throttle: nominal vs safe.
    ax[1, 1].plot(t, d["dT_nom"], color="C3", ls="--", label="throttle nominal")
    ax[1, 1].plot(t, d["dT"], color="C0", label="throttle safe (CBF)")
    ax[1, 1].set_ylabel("throttle [-]")
    ax[1, 1].set_title("Throttle: nominal vs safe")
    ax[1, 1].legend(); ax[1, 1].grid(True)

    # Airspeed barrier value.
    _barrier_panel(ax[2, 0], t, d["b_airspeed"], "C4",
                   "Airspeed barrier  b_V = V - V_min")
    ax[2, 0].set_xlabel("time [s]")

    # Descent-rate barrier value.
    _barrier_panel(ax[2, 1], t, d["b_descent"], "C5",
                   "Descent barrier  b = hdot + sqrt(v_safe^2 + 2 a_brk h)")
    ax[2, 1].set_xlabel("time [s]")

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    if args.save:
        fig.savefig(args.save, dpi=120)
        print(f"saved {args.save}")
    else:
        plt.show()


def _barrier_panel(axis, t, b, color, title):
    axis.plot(t, b, color=color, label="b(x)")
    axis.axhline(0.0, color="red", ls="--", lw=1.2, label="b = 0 (safety boundary)")
    lo = min(b.min(), 0.0)
    axis.fill_between(t, lo - 0.1 * abs(lo) - 1e-3, 0.0, color="red", alpha=0.06)
    axis.set_ylabel("barrier value"); axis.set_title(title)
    axis.legend(); axis.grid(True)


if __name__ == "__main__":
    main()
