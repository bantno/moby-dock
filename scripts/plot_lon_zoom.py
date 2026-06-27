#!/usr/bin/env python3
"""Zoomed view of a lon_autoland_sim trace: the last WINDOW seconds before
touchdown. Same 6-panel layout as plot_lon_results.py, restricted to the final
approach/flare so the touchdown behaviour is legible.

Usage: plot_lon_zoom.py [csv] [out.png] [v_safe] [a_brk] [Vmin] [Tmax] [window_s]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

csv = sys.argv[1] if len(sys.argv) > 1 else "lon_autoland_log.csv"
out = sys.argv[2] if len(sys.argv) > 2 else "lon_trajectory_zoom.png"
v_safe = float(sys.argv[3]) if len(sys.argv) > 3 else 0.6
a_brk = float(sys.argv[4]) if len(sys.argv) > 4 else 3.0
Vmin = float(sys.argv[5]) if len(sys.argv) > 5 else 13.5
Tmax = float(sys.argv[6]) if len(sys.argv) > 6 else 6.0
window = float(sys.argv[7]) if len(sys.argv) > 7 else 5.0

d = np.genfromtxt(csv, delimiter=",", names=True)
t_all = d["t"]
t_td = t_all.max()                       # CSV ends at touchdown
t0 = t_td - window
m = t_all >= t0                          # last `window` seconds
t = t_all[m]

fig, ax = plt.subplots(3, 2, figsize=(13, 10))
fig.suptitle(f"Augmented-longitudinal CBF-QP water landing "
             f"— final {window:g} s before touchdown (t={t_td:.2f} s)",
             fontsize=14, weight="bold")

# (1) Altitude
ax[0, 0].plot(t, d["h"][m], color="tab:blue")
ax[0, 0].axhline(0, color="k", lw=0.8, ls=":")
ax[0, 0].set_ylabel("altitude h [m]")
ax[0, 0].set_title("Descent profile")

# (2) Sink rate vs the descent-barrier envelope
env = np.sqrt(v_safe**2 + 2 * a_brk * np.clip(d["h"][m], 0, None))
ax[0, 1].plot(t, d["sink"][m], color="tab:red", label="sink rate")
ax[0, 1].plot(t, env, color="tab:gray", ls="--", label=r"envelope $\sqrt{v_{safe}^2+2a_{brk}h}$")
ax[0, 1].axhline(v_safe, color="k", ls=":", lw=0.9, label=f"v_safe={v_safe}")
ax[0, 1].set_ylabel("sink [m/s] (down +)")
ax[0, 1].set_title("Sink rate vs descent envelope (flare emerges from CBF)")
ax[0, 1].legend(fontsize=8)

# (3) Airspeed vs floor
ax[1, 0].plot(t, d["V"][m], color="tab:green")
ax[1, 0].axhline(Vmin, color="k", ls=":", lw=0.9, label=f"Vmin={Vmin}")
ax[1, 0].set_ylabel("V [m/s]")
ax[1, 0].set_title("Airspeed vs stall floor")
ax[1, 0].legend(fontsize=8)

# (4) Attitude / path angles
ax[1, 1].plot(t, d["theta_deg"][m], label=r"$\theta$")
ax[1, 1].plot(t, d["theta_cmd_deg"][m], ls="--", label=r"$\theta_{cmd}$")
ax[1, 1].plot(t, d["gamma_deg"][m], label=r"$\gamma$")
ax[1, 1].plot(t, d["alpha_deg"][m], label=r"$\alpha$")
ax[1, 1].set_ylabel("angle [deg]")
ax[1, 1].set_title("Attitude / flight-path angles")
ax[1, 1].legend(fontsize=8, ncol=2)

# (5) Controls: nominal vs filtered
ax[2, 0].plot(t, np.rad2deg(d["de"][m]), color="tab:purple", label=r"$\delta_e$ (filtered)")
ax[2, 0].plot(t, np.rad2deg(d["de_nom"][m]), color="tab:purple", ls=":", alpha=0.7, label=r"$\delta_e$ (nominal)")
ax[2, 0].set_ylabel("elevator [deg]")
ax[2, 0].set_xlabel("t [s]")
ax[2, 0].set_title("Elevator: nominal vs CBF-filtered")
ax[2, 0].legend(fontsize=8)
axT = ax[2, 0].twinx()
axT.plot(t, d["T"][m], color="tab:orange", alpha=0.6, label="T [N]")
axT.set_ylabel("thrust T [N]", color="tab:orange")

# (6) All barrier values over time (must stay >= 0 under the QP).
bars = {
    "b_descent (sink)": d["b_descent"][m],
    "b_airspeed (V-Vmin)": d["b_airspeed"][m],
    "b_thrust_min (T)": d["T"][m],
    "b_thrust_max (Tmax-T)": Tmax - d["T"][m],
}
for name, vals in bars.items():
    ax[2, 1].plot(t, vals, label=f"{name}  (min {vals.min():.3f})")
    print(f"  [last {window:g}s] min {name:24s} = {vals.min():+.4f}  -> {'OK' if vals.min() >= -1e-6 else 'VIOLATION'}")
ax[2, 1].axhline(0, color="k", lw=1.0, ls="--")
ax[2, 1].set_ylabel("barrier value")
ax[2, 1].set_xlabel("t [s]")
ax[2, 1].set_title("All barriers (>= 0 = safe) under the CBF-QP")
ax[2, 1].legend(fontsize=7)

for a in ax.flat:
    a.grid(alpha=0.3)
    a.set_xlim(t0, t_td)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(out, dpi=110)
print("wrote", out)
