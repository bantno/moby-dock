#!/usr/bin/env python3
"""Plot a CBF-filter ON vs OFF comparison from two autoland_sim CSV logs.

Companion to plot_results.py. Pass the two logs (same scenario, CBF enabled in
one and disabled in the other) and the airspeed-barrier margin:

    python3 scripts/cbf_compare.py cbf_on.csv cbf_off.csv --vmin 15 --save fig.png
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
    ap = argparse.ArgumentParser(description="CBF ON vs OFF comparison plot")
    ap.add_argument("on_csv", help="log with the CBF filter ENABLED")
    ap.add_argument("off_csv", help="log with the CBF filter DISABLED")
    ap.add_argument("--vmin", type=float, default=15.0,
                    help="airspeed-barrier margin V_min [m/s]")
    ap.add_argument("--save", help="save figure to this path instead of showing")
    args = ap.parse_args()

    on, off = load(args.on_csv), load(args.off_csv)
    vmin = args.vmin

    fig, ax = plt.subplots(2, 2, figsize=(13, 8))
    fig.suptitle("CBF safety filter ON vs OFF — approach pushed to the stall margin",
                 fontsize=13)

    # Airspeed vs the barrier.
    ax[0, 0].plot(on["t"], on["V"], color="C0", label="V  (CBF ON)")
    ax[0, 0].plot(off["t"], off["V"], color="C3", label="V  (CBF OFF)")
    ax[0, 0].axhline(vmin, color="k", ls="--", lw=1.2, label=f"V_min barrier = {vmin}")
    ax[0, 0].fill_between(on["t"], 0, vmin, color="red", alpha=0.06)
    lo = min(on["V"].min(), off["V"].min(), vmin) - 1.0
    hi = max(on["V"].max(), off["V"].max()) + 0.5
    ax[0, 0].set_ylim(lo, hi)
    ax[0, 0].set_ylabel("airspeed [m/s]"); ax[0, 0].set_title("Airspeed defends the barrier")
    ax[0, 0].legend(); ax[0, 0].grid(True)

    # Throttle.
    ax[0, 1].plot(on["t"], on["dT"], color="C0", label="throttle (ON)")
    ax[0, 1].plot(off["t"], off["dT"], color="C3", label="throttle (OFF)")
    ax[0, 1].set_ylabel("throttle [-]"); ax[0, 1].set_title("Filter trims throttle/elevator to hold V")
    ax[0, 1].legend(); ax[0, 1].grid(True)

    # Elevator.
    ax[1, 0].plot(on["t"], on["de_deg"], color="C0", label="delta_e (ON)")
    ax[1, 0].plot(off["t"], off["de_deg"], color="C3", label="delta_e (OFF)")
    ax[1, 0].set_xlabel("time [s]"); ax[1, 0].set_ylabel("elevator [deg]")
    ax[1, 0].legend(); ax[1, 0].grid(True)

    # Altitude.
    ax[1, 1].plot(on["t"], on["h"], color="C0", label="h (ON)")
    ax[1, 1].plot(off["t"], off["h"], color="C3", label="h (OFF)")
    ax[1, 1].axhline(0, color="k", lw=0.6)
    ax[1, 1].set_xlabel("time [s]"); ax[1, 1].set_ylabel("altitude [m]")
    ax[1, 1].set_title("Glidepath to touchdown")
    ax[1, 1].legend(); ax[1, 1].grid(True)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    if args.save:
        fig.savefig(args.save, dpi=120)
        print(f"saved {args.save}  (min V: ON={on['V'].min():.3f}, OFF={off['V'].min():.3f})")
    else:
        plt.show()


if __name__ == "__main__":
    main()
