#!/usr/bin/env python3
"""Overlay CBF-on vs CBF-off (nominal-only) longitudinal landing traces.

Usage: plot_lon_compare.py cbf_on.csv cbf_off.csv out.png [v_safe] [a_brk] [Vmin]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

on = np.genfromtxt(sys.argv[1], delimiter=",", names=True)
off = np.genfromtxt(sys.argv[2], delimiter=",", names=True)
out = sys.argv[3] if len(sys.argv) > 3 else "lon_compare.png"
v_safe = float(sys.argv[4]) if len(sys.argv) > 4 else 0.6
a_brk = float(sys.argv[5]) if len(sys.argv) > 5 else 3.0
Vmin = float(sys.argv[6]) if len(sys.argv) > 6 else 13.5

C_ON, C_OFF = "tab:blue", "tab:red"
fig, ax = plt.subplots(2, 2, figsize=(13, 9))
fig.suptitle("CBF-QP filter vs nominal-only (no CBF)", fontsize=14, weight="bold")


def last(d):  # touchdown index = first h<=0 (or end)
    idx = np.where(d["h"] <= 0)[0]
    return idx[0] if len(idx) else len(d["h"]) - 1


i_on, i_off = last(on), last(off)

# (1) Altitude vs time
ax[0, 0].plot(on["t"], on["h"], C_ON, label="CBF on")
ax[0, 0].plot(off["t"], off["h"], C_OFF, label="nominal only")
ax[0, 0].axhline(0, color="k", lw=0.8, ls=":")
ax[0, 0].set(xlabel="t [s]", ylabel="altitude h [m]", title="Descent profile")
ax[0, 0].legend()

# (2) Sink rate vs altitude (descent envelope clamp is clearest here).
# Altitude on x (high -> touchdown at right), sink on y.
hgrid = np.linspace(0, max(on["h"].max(), off["h"].max()), 200)
env = np.sqrt(v_safe**2 + 2 * a_brk * hgrid)
ax[0, 1].plot(hgrid, env, color="tab:gray", ls="--", label=r"CBF envelope $\sqrt{v_{safe}^2+2a_{brk}h}$")
ax[0, 1].plot(on["h"], on["sink"], C_ON, label="CBF on")
ax[0, 1].plot(off["h"], off["sink"], C_OFF, label="nominal only")
ax[0, 1].axhline(v_safe, color="k", ls=":", lw=0.9, label=f"v_safe={v_safe}")
ax[0, 1].set(xlabel="altitude h [m]", ylabel="sink [m/s] (down +)",
             title="Sink vs altitude (CBF clamps to envelope; nominal sails through)")
ax[0, 1].set_xlim(12, 0)  # descend left -> right (touchdown at h=0 on the right)
ax[0, 1].set_ylim(-0.1, max(on["sink"].max(), off["sink"].max()) * 1.1)
ax[0, 1].legend(fontsize=8)

# (3) Airspeed vs floor
ax[1, 0].plot(on["t"], on["V"], C_ON, label="CBF on")
ax[1, 0].plot(off["t"], off["V"], C_OFF, label="nominal only")
ax[1, 0].axhline(Vmin, color="k", ls=":", lw=0.9, label=f"Vmin={Vmin}")
ax[1, 0].set(xlabel="t [s]", ylabel="V [m/s]", title="Airspeed")
ax[1, 0].legend(fontsize=8)

# (4) Touchdown sink bar + flight-path angle
ax[1, 1].plot(on["t"], on["gamma_deg"], C_ON, label=r"$\gamma$ CBF on")
ax[1, 1].plot(off["t"], off["gamma_deg"], C_OFF, label=r"$\gamma$ nominal only")
ax[1, 1].plot(on["t"], on["theta_deg"], C_ON, ls="--", alpha=0.7, label=r"$\theta$ CBF on")
ax[1, 1].plot(off["t"], off["theta_deg"], C_OFF, ls="--", alpha=0.7, label=r"$\theta$ nominal only")
ax[1, 1].set(xlabel="t [s]", ylabel="angle [deg]",
             title=f"Path/attitude  (touchdown sink: on={on['sink'][i_on]:.3f}, "
                   f"off={off['sink'][i_off]:.3f} m/s)")
ax[1, 1].legend(fontsize=8, ncol=2)

for a in ax.flat:
    a.grid(alpha=0.3)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(out, dpi=120)
print("wrote", out)
