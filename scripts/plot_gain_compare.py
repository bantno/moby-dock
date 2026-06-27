#!/usr/bin/env python3
"""Compare two lon_autoland_sim runs (e.g. different class-K gains).

Each physical quantity gets its OWN axes (no mixing of different units/scales):
altitude, descent rate, elevator (nominal vs QP), Tddot (nominal vs QP), and the
four barrier functions b (descent, airspeed, min-thrust, max-thrust). b >= 0 is
the safety criterion -- if the b panels stay positive, the run was safe.

Usage:
  plot_gain_compare.py run_a.csv labelA run_b.csv labelB out.png [Tmax]

Generate the two CSVs first, e.g.:
  ./lon_autoland_sim data/AHAB_combined.stab data/aircraft.yaml scn_a.yaml run_a.csv
  ./lon_autoland_sim data/AHAB_combined.stab data/aircraft.yaml scn_b.yaml run_b.csv
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

a_csv, a_lbl, b_csv, b_lbl, out = sys.argv[1:6]
Tmax = float(sys.argv[6]) if len(sys.argv) > 6 else 50.0

A = np.genfromtxt(a_csv, delimiter=",", names=True)
B = np.genfromtxt(b_csv, delimiter=",", names=True)
runs = [(A, a_lbl, "tab:blue"), (B, b_lbl, "tab:red")]

fig, ax = plt.subplots(4, 2, figsize=(14, 15))
fig.suptitle(f"Longitudinal autoland: {a_lbl} vs {b_lbl}", fontsize=14)

def plot(axis, getter, title, ylabel, zero=False):
    for D, lbl, col in runs:
        axis.plot(D["t"], getter(D), color=col, label=lbl, lw=1.4)
    if zero:
        axis.axhline(0.0, color="k", lw=0.8, ls=":")
    axis.set_title(title)
    axis.set_xlabel("t [s]")
    axis.set_ylabel(ylabel)
    axis.grid(alpha=0.3)
    axis.legend(fontsize=8)

# Row 1: trajectory
plot(ax[0, 0], lambda D: D["h"], "Altitude", "h [m]", zero=True)
plot(ax[0, 1], lambda D: D["sink"], "Descent rate (sink, +down)", "sink [m/s]", zero=True)

# Row 2: controls -- nominal (dashed) vs QP output (solid), each control own axes
for D, lbl, col in runs:
    ax[1, 0].plot(D["t"], D["de_nom"], color=col, ls="--", lw=1.0, label=f"{lbl} nominal")
    ax[1, 0].plot(D["t"], D["de"], color=col, ls="-", lw=1.5, label=f"{lbl} QP")
ax[1, 0].set_title("Elevator command de (nominal vs QP)")
ax[1, 0].set_xlabel("t [s]"); ax[1, 0].set_ylabel("de [rad]")
ax[1, 0].grid(alpha=0.3); ax[1, 0].legend(fontsize=8)

for D, lbl, col in runs:
    ax[1, 1].plot(D["t"], D["Tddot_nom"], color=col, ls="--", lw=1.0, label=f"{lbl} nominal")
    ax[1, 1].plot(D["t"], D["Tddot"], color=col, ls="-", lw=1.5, label=f"{lbl} QP")
ax[1, 1].set_title("Thrust 2nd-deriv command Tddot (nominal vs QP)")
ax[1, 1].set_xlabel("t [s]"); ax[1, 1].set_ylabel("Tddot [N/s^2]")
ax[1, 1].grid(alpha=0.3); ax[1, 1].legend(fontsize=8)

# Row 3-4: the four barrier functions (each own axes). b >= 0 is the safety
# criterion -- if these stay positive the run was safe.
plot(ax[2, 0], lambda D: D["b_descent"], "Barrier: descent-rate b", "b_descent", zero=True)
plot(ax[2, 1], lambda D: D["b_airspeed"], "Barrier: airspeed b_V", "b_airspeed", zero=True)
plot(ax[3, 0], lambda D: D["T"], "Barrier: min-thrust (b = T)", "T [N]", zero=True)
plot(ax[3, 1], lambda D: Tmax - D["T"], "Barrier: max-thrust (b = Tmax - T)", "Tmax - T [N]", zero=True)

fig.tight_layout(rect=[0, 0, 1, 0.985])
fig.savefig(out, dpi=110)
print("wrote", out)
