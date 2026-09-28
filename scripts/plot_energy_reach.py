#!/usr/bin/env python3
"""Plot the reach-the-target energy-CBF runs (sixdof_autoland_sim CSVs).

Two modes:
  * A/B comparison (filter ON vs OFF, same scenario):
      python3 scripts/plot_energy_reach.py runs/beaver_reach_glide.csv \
          --off runs/beaver_reach_glide_off.csv --save figures/reach_glide.png
  * single run:
      python3 scripts/plot_energy_reach.py runs/beaver_reach_powered.csv \
          --save figures/reach_powered.png

Panels: barrier h_E, ground track vs the target, airspeed + AoA (command vs
flown), heading error phi with the sideslip command, arc distance to target,
and the CBF's alpha-rate channel.
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


def target_of(d):
    """Recover the target from the trace: d_tgt is the ARC distance, which
    equals the straight-line distance when phi ~ 0; take it at the last sample
    with the smallest |phi| unavailable -- instead reconstruct from the final
    point plus its distance along the bearing. Simpler: the run prints the
    target; here we only need it for the track panel, so reconstruct from the
    closest-approach sample assuming the target sits on d_tgt along the
    line-of-sight: N_t = x + d*cos(chi - phi), E_t = y + d*sin(chi - phi)."""
    chi = np.arctan2(np.gradient(d["y"]), np.gradient(d["x"]) + 1e-12)
    phi = np.radians(d["phi_arc_deg"])
    los = chi - phi
    k = len(d["t"]) // 2  # mid-flight sample: phi small, geometry clean
    return (d["x"][k] + d["d_tgt"][k] * np.cos(los[k]),
            d["y"][k] + d["d_tgt"][k] * np.sin(los[k]))


def main():
    ap = argparse.ArgumentParser(description="energy-reach CBF plots")
    ap.add_argument("on_csv", help="log with the reach layer (filter on)")
    ap.add_argument("--off", help="A/B log with filter: false (optional)")
    ap.add_argument("--save", help="save figure here instead of showing")
    ap.add_argument("--title", default="Reach-the-target energy CBF (3-DOF "
                    "formulation on the 6-DOF Beaver)")
    args = ap.parse_args()

    on = load(args.on_csv)
    off = load(args.off) if args.off else None

    fig, ax = plt.subplots(3, 2, figsize=(13, 11))
    fig.suptitle(args.title)

    a = ax[0][0]
    a.plot(on["t"], on["hE"], "C0", label="CBF on")
    if off is not None:
        a.plot(off["t"], off["hE"], "C3--", label="filter off")
    a.axhline(0.0, color="k", lw=0.8)
    a.set_ylabel(r"$h_E$ [m$^2$/s$^2$]")
    a.set_xlabel("t [s]")
    a.set_title("reach-energy barrier")
    a.legend()

    a = ax[0][1]
    a.plot(on["x"], on["y"], "C0", label="CBF on")
    if off is not None:
        a.plot(off["x"], off["y"], "C3--", label="filter off")
    tN, tE = target_of(on)
    a.plot([tN], [tE], "k*", ms=14, label="target")
    a.set_xlabel("north [m]")
    a.set_ylabel("east [m]")
    a.set_title("ground track")
    a.axis("equal")
    a.legend()

    a = ax[1][0]
    a.plot(on["t"], on["alpha_deg"], "C0", label=r"$\alpha$ (CBF on)")
    a.plot(on["t"], on["alpha_cmd_deg"], "C1", lw=0.9, label=r"$\alpha_{cmd}$")
    if off is not None:
        a.plot(off["t"], off["alpha_deg"], "C3--", label=r"$\alpha$ (off)")
    a.set_ylabel(r"$\alpha$ [deg]")
    a.set_xlabel("t [s]")
    a.set_title("angle of attack: filtered command vs flown")
    a.legend()

    a = ax[1][1]
    a.plot(on["t"], on["V_air"], "C0", label="V (CBF on)")
    if off is not None:
        a.plot(off["t"], off["V_air"], "C3--", label="V (off)")
    a.set_ylabel("V [m/s]")
    a.set_xlabel("t [s]")
    a.set_title("airspeed")
    a.legend()

    a = ax[2][0]
    a.plot(on["t"], on["phi_arc_deg"], "C0", label=r"$\phi$ (CBF on)")
    a.plot(on["t"], on["beta_cmd_deg"], "C2", lw=0.9, label=r"$\beta_{cmd}$")
    a.plot(on["t"], on["beta_deg"], "C1", lw=0.9, label=r"$\beta$ flown")
    if off is not None:
        a.plot(off["t"], off["phi_arc_deg"], "C3--", label=r"$\phi$ (off)")
    a.set_ylabel("[deg]")
    a.set_xlabel("t [s]")
    a.set_title("heading error to bearing + sideslip channel")
    a.legend()

    a = ax[2][1]
    a.plot(on["t"], on["d_tgt"], "C0", label="arc d (CBF on)")
    if off is not None:
        a.plot(off["t"], off["d_tgt"], "C3--", label="arc d (off)")
    a2 = a.twinx()
    a2.plot(on["t"], on["LD_inst"], "C2", lw=0.9, label="L/D inst")
    a2.set_ylabel("instantaneous L/D", color="C2")
    a.set_ylabel("d [m]")
    a.set_xlabel("t [s]")
    a.set_title("arc distance to target / instantaneous L/D")
    a.legend(loc="upper right")

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    if args.save:
        fig.savefig(args.save, dpi=130)
        print(f"saved {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
