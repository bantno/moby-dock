#!/usr/bin/env python3
"""Plot a 6-DOF straight-in landing trace (sixdof_autoland_sim CSV).

Nine panels: descent (h + wave surface), airspeed, longitudinal angles,
lateral attitude/sideslip, heading, body rates, surface deflections,
throttle, and cross-track (with the lateral gust overlay).

Usage: plot_sixdof_results.py [sixdof_autoland_log.csv] [out.png]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

csv = sys.argv[1] if len(sys.argv) > 1 else "sixdof_autoland_log.csv"
out = sys.argv[2] if len(sys.argv) > 2 else "sixdof_trajectory.png"

d = np.genfromtxt(csv, delimiter=",", names=True)
t = d["t"]

fig, ax = plt.subplots(3, 3, figsize=(16, 10))
fig.suptitle("6-DOF straight-in water landing (cascaded-PID nominal)",
             fontsize=14, weight="bold")

# (1) Descent profile with the instantaneous wave surface under the keel.
ax[0, 0].plot(t, d["h"], color="tab:blue", label="h")
if np.any(d["eta"] != 0.0):
    ax[0, 0].plot(t, d["eta"], color="tab:cyan", lw=0.8, label=r"$\eta(x,t)$")
ax[0, 0].axhline(0, color="k", lw=0.8, ls=":")
ax[0, 0].set_ylabel("altitude h [m]")
ax[0, 0].set_title("Descent profile")
ax[0, 0].legend(fontsize=8)

# (2) Airspeed + sink rate.
ax[0, 1].plot(t, d["V_air"], color="tab:green", label="V_air")
ax2 = ax[0, 1].twinx()
ax2.plot(t, d["sink"], color="tab:red", lw=0.9, label="sink")
ax2.set_ylabel("sink [m/s] (down +)", color="tab:red")
ax[0, 1].set_ylabel("V_air [m/s]", color="tab:green")
ax[0, 1].set_title("Airspeed / sink rate")

# (3) Longitudinal angles.
ax[0, 2].plot(t, d["theta_deg"], label=r"$\theta$")
ax[0, 2].plot(t, d["theta_cmd_deg"], ls="--", label=r"$\theta_{cmd}$")
ax[0, 2].plot(t, d["gamma_deg"], label=r"$\gamma$")
ax[0, 2].plot(t, d["alpha_deg"], label=r"$\alpha$")
ax[0, 2].set_ylabel("[deg]")
ax[0, 2].set_title("Longitudinal angles")
ax[0, 2].legend(fontsize=8)

# (4) Lateral attitude + sideslip.
ax[1, 0].plot(t, d["phi_deg"], label=r"$\phi$")
ax[1, 0].plot(t, d["phi_cmd_deg"], ls="--", label=r"$\phi_{cmd}$")
ax[1, 0].plot(t, d["beta_deg"], label=r"$\beta$")
ax[1, 0].axhline(0, color="k", lw=0.8, ls=":")
ax[1, 0].set_ylabel("[deg]")
ax[1, 0].set_title("Bank / sideslip")
ax[1, 0].legend(fontsize=8)

# (5) Heading (crab develops under crosswind).
ax[1, 1].plot(t, d["psi_deg"], color="tab:purple")
ax[1, 1].axhline(0, color="k", lw=0.8, ls=":")
ax[1, 1].set_ylabel(r"$\psi$ [deg]")
ax[1, 1].set_title("Heading (crab)")

# (6) Body rates.
ax[1, 2].plot(t, d["p"], label="p")
ax[1, 2].plot(t, d["q"], label="q")
ax[1, 2].plot(t, d["r"], label="r")
ax[1, 2].set_ylabel("[rad/s]")
ax[1, 2].set_title("Body rates")
ax[1, 2].legend(fontsize=8)

# (7) Surface deflections.
deg = 180.0 / np.pi
ax[2, 0].plot(t, d["de"] * deg, label=r"$\delta_e$")
ax[2, 0].plot(t, d["da"] * deg, label=r"$\delta_a$")
ax[2, 0].plot(t, d["dr"] * deg, label=r"$\delta_r$")
ax[2, 0].set_ylabel("[deg]")
ax[2, 0].set_xlabel("t [s]")
ax[2, 0].set_title("Surface deflections")
ax[2, 0].legend(fontsize=8)

# (8) Throttle.
ax[2, 1].plot(t, d["dT"], color="tab:orange")
ax[2, 1].set_ylim(-0.02, 1.02)
ax[2, 1].set_ylabel(r"$\delta_T$ [-]")
ax[2, 1].set_xlabel("t [s]")
ax[2, 1].set_title("Throttle")

# (9) Cross-track, with the (earth-frame) gust components if any.
ax[2, 2].plot(t, d["y"], color="tab:blue", label="y")
ax[2, 2].axhline(0, color="k", lw=0.8, ls=":")
if np.any(d["W_u"] != 0) or np.any(d["W_v"] != 0) or np.any(d["W_h"] != 0):
    ax3 = ax[2, 2].twinx()
    ax3.plot(t, d["W_u"], color="tab:gray", lw=0.8, label="W_u (tail+)")
    ax3.plot(t, d["W_v"], color="tab:olive", lw=0.8, label="W_v (east+)")
    ax3.plot(t, d["W_h"], color="tab:cyan", lw=0.8, label="W_h (up+)")
    ax3.set_ylabel("gust [m/s]")
    ax3.legend(fontsize=7, loc="lower right")
ax[2, 2].set_ylabel("cross-track y [m]")
ax[2, 2].set_xlabel("t [s]")
ax[2, 2].set_title("Cross-track / gust")
ax[2, 2].legend(fontsize=8, loc="upper left")

for a in ax.flat:
    a.grid(alpha=0.3)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(out, dpi=140)
print(f"wrote {out}")
