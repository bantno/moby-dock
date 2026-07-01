#!/usr/bin/env python3
"""Plot an augmented-longitudinal CBF-QP landing trace (lon_autoland_sim CSV).

Recovery barrier set: impact-load (hard) + stall/AoA + nose-up + total-energy
ceiling, plus the thrust actuator guards.

Usage: plot_lon_results.py [lon_autoland_log.csv] [out.png] [V_td_max] [g_eff] [Tmax]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

csv = sys.argv[1] if len(sys.argv) > 1 else "lon_autoland_log.csv"
out = sys.argv[2] if len(sys.argv) > 2 else "lon_trajectory.png"
V_td_max = float(sys.argv[3]) if len(sys.argv) > 3 else 14.0
g_eff = float(sys.argv[4]) if len(sys.argv) > 4 else 16.0
Tmax = float(sys.argv[5]) if len(sys.argv) > 5 else 50.0
G = 9.80665
ALPHA_LIM_DEG = 9.0   # alpha_stall - margin (reference line only)
THETA_MIN_DEG = 3.0

d = np.genfromtxt(csv, delimiter=",", names=True)
t = d["t"]

fig, ax = plt.subplots(3, 2, figsize=(13, 10))
fig.suptitle("Augmented-longitudinal CBF-QP water landing", fontsize=14, weight="bold")

# (1) Altitude
ax[0, 0].plot(t, d["h"], color="tab:blue")
ax[0, 0].axhline(0, color="k", lw=0.8, ls=":")
ax[0, 0].set_ylabel("altitude h [m]")
ax[0, 0].set_title("Descent profile")

# (2) Sink rate (the impact barrier regulates the touchdown sink)
ax[0, 1].plot(t, d["sink"], color="tab:red", label="sink rate")
ax[0, 1].axhline(0, color="k", ls=":", lw=0.9)
ax[0, 1].set_ylabel("sink [m/s] (down +)")
ax[0, 1].set_title("Sink rate (impact barrier regulates touchdown)")
ax[0, 1].legend(fontsize=8)

# (3) Airspeed vs the energy-ceiling cap V <= sqrt(V_td_max^2 + 2(g_eff-g)h)
vcap = np.sqrt(V_td_max**2 + 2 * (g_eff - G) * np.clip(d["h"], 0, None))
ax[1, 0].plot(t, d["V"], color="tab:green", label="V")
ax[1, 0].plot(t, vcap, color="tab:gray", ls="--",
              label=r"energy cap $\sqrt{V_{td}^2+2(g_{eff}-g)h}$")
ax[1, 0].axhline(V_td_max, color="k", ls=":", lw=0.9, label=f"V_td_max={V_td_max}")
ax[1, 0].set_ylabel("V [m/s]")
ax[1, 0].set_title("Airspeed vs energy-ceiling cap")
ax[1, 0].legend(fontsize=8)

# (4) Attitude / path angles, with the stall + nose-up reference lines
ax[1, 1].plot(t, d["theta_deg"], label=r"$\theta$")
ax[1, 1].plot(t, d["theta_cmd_deg"], ls="--", label=r"$\theta_{cmd}$")
ax[1, 1].plot(t, d["gamma_deg"], label=r"$\gamma$")
ax[1, 1].plot(t, d["alpha_deg"], label=r"$\alpha$")
ax[1, 1].axhline(ALPHA_LIM_DEG, color="tab:red", ls=":", lw=0.9,
                 label=f"alpha_lim={ALPHA_LIM_DEG:.0f}")
ax[1, 1].axhline(THETA_MIN_DEG, color="tab:brown", ls=":", lw=0.9,
                 label=f"theta_min={THETA_MIN_DEG:.0f}")
ax[1, 1].set_ylabel("angle [deg]")
ax[1, 1].set_title("Attitude / flight-path angles")
ax[1, 1].legend(fontsize=7, ncol=2)

# (5) Controls: nominal vs filtered
ax[2, 0].plot(t, np.rad2deg(d["de"]), color="tab:purple", label=r"$\delta_e$ (filtered)")
ax[2, 0].plot(t, np.rad2deg(d["de_nom"]), color="tab:purple", ls=":", alpha=0.7,
              label=r"$\delta_e$ (nominal)")
ax[2, 0].set_ylabel("elevator [deg]")
ax[2, 0].set_xlabel("t [s]")
ax[2, 0].set_title("Elevator: nominal vs CBF-filtered")
ax[2, 0].legend(fontsize=8)
axT = ax[2, 0].twinx()
axT.plot(t, d["T"], color="tab:orange", alpha=0.6, label="T [N]")
axT.set_ylabel("thrust T [N]", color="tab:orange")

# (6) All barrier values over time. Impact is the only HARD safety row; stall /
#     nose-up / energy are soft (may dip). Actuator guards: b=T and b=Tmax-T.
bars = {
    "b_stall": d["b_stall"],
    "b_noseup": d["b_noseup"],
    "b_energy": d["b_energy"],
    "b_impact (HARD)": d["b_impact"],
    "b_thrust_min (T)": d["T"],
    "b_thrust_max (Tmax-T)": Tmax - d["T"],
}
for name, vals in bars.items():
    ax[2, 1].plot(t, vals, label=f"{name}  (min {np.nanmin(vals):.3f})")
    print(f"  min {name:24s} = {np.nanmin(vals):+.4f}")
ax[2, 1].axhline(0, color="k", lw=1.0, ls="--")
ax[2, 1].set_ylabel("barrier value")
ax[2, 1].set_xlabel("t [s]")
ax[2, 1].set_title("Barriers (impact hard >= 0; soft rows may dip)")
ax[2, 1].legend(fontsize=7)

for a in ax.flat:
    a.grid(alpha=0.3)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(out, dpi=110)
print("wrote", out)
