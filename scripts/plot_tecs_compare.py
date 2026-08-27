#!/usr/bin/env python3
"""Overlay a cascaded-PID and a PX4-TECS 6-DOF landing trace (sixdof_autoland_sim CSVs).

Six panels: altitude, airspeed, sink rate vs reference, pitch / pitch command,
throttle, elevator. Both runs must share the scenario (only nominal.type differs)
for the overlay to mean anything.

Usage: plot_tecs_compare.py cascade.csv tecs.csv [out.png] [title]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

if len(sys.argv) < 3:
    sys.exit(__doc__)
a = np.genfromtxt(sys.argv[1], delimiter=",", names=True)
b = np.genfromtxt(sys.argv[2], delimiter=",", names=True)
out = sys.argv[3] if len(sys.argv) > 3 else "figures/sixdof_tecs_compare.png"
title = sys.argv[4] if len(sys.argv) > 4 else "Beaver straight-in: cascaded PID vs PX4 TECS nominal"
deg = 180.0 / np.pi

fig, ax = plt.subplots(2, 3, figsize=(16, 8))
fig.suptitle(title, fontsize=14, weight="bold")
runs = [(a, "cascaded PID", "tab:blue"), (b, "PX4 TECS", "tab:orange")]

for d, lab, c in runs:
    ax[0, 0].plot(d["t"], d["h"], color=c, label=lab)
    ax[0, 1].plot(d["t"], d["V_air"], color=c, label=lab)
    ax[0, 2].plot(d["t"], d["sink"], color=c, label=f"{lab} sink")
    ax[0, 2].plot(d["t"], -d["hdot_sp"], color=c, ls="--", lw=0.8, label=f"{lab} ref")
    ax[1, 0].plot(d["t"], d["theta_deg"], color=c, label=f"{lab} " + r"$\theta$")
    ax[1, 0].plot(d["t"], d["theta_cmd_deg"], color=c, ls="--", lw=0.8,
                  label=f"{lab} " + r"$\theta_{cmd}$")
    ax[1, 1].plot(d["t"], d["dT"], color=c, label=lab)
    ax[1, 2].plot(d["t"], d["de"] * deg, color=c, label=lab)

ax[0, 0].set_title("Altitude"); ax[0, 0].set_ylabel("h [m]")
ax[0, 1].set_title("Airspeed"); ax[0, 1].set_ylabel("V_air [m/s]")
ax[0, 2].set_title("Sink rate (down +) vs reference"); ax[0, 2].set_ylabel("[m/s]")
ax[1, 0].set_title("Pitch / pitch command"); ax[1, 0].set_ylabel("[deg]")
ax[1, 1].set_title("Throttle"); ax[1, 1].set_ylabel(r"$\delta_T$ [-]")
ax[1, 2].set_title("Elevator"); ax[1, 2].set_ylabel(r"$\delta_e$ [deg]")
for row in ax:
    for a_ in row:
        a_.grid(alpha=0.3)
        a_.legend(fontsize=7)
        a_.set_xlabel("t [s]")
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(out, dpi=130)
print("wrote", out)
