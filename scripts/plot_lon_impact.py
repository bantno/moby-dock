#!/usr/bin/env python3
"""Plot the hydrodynamic impact-load barrier over a lon_autoland_sim trace.

Shows the trajectory plus the impact-barrier diagnostics logged by the sim
(b_impact, n_peak, kappa_imp, psi1/2_imp, res_imp) and the elevator nominal-vs-
filtered command, so you can see whether/where the barrier engages.

Usage: plot_lon_impact.py [lon_autoland_log.csv] [out.png] [n_limit] [z_gate]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

csv = sys.argv[1] if len(sys.argv) > 1 else "lon_autoland_log.csv"
out = sys.argv[2] if len(sys.argv) > 2 else "lon_impact.png"
n_limit = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0
z_gate = float(sys.argv[4]) if len(sys.argv) > 4 else 10.0

d = np.genfromtxt(csv, delimiter=",", names=True)
t, h = d["t"], d["h"]
# active, model-valid window: below z_gate, descending, positive trim
active = (h < z_gate) & (-d["gamma_deg"] > 0) & (d["theta_deg"] > 0)

fig, ax = plt.subplots(4, 2, figsize=(13, 13))
fig.suptitle("Impact-load barrier (NACA TN 1516) over the landing", fontsize=14, weight="bold")


def shade_active(a):
    a.fill_between(t, *a.get_ylim(), where=active, color="tab:orange", alpha=0.08,
                   label="impact row active")


# (1) Altitude with z_gate
ax[0, 0].plot(t, h, color="tab:blue")
ax[0, 0].axhline(z_gate, color="tab:orange", ls=":", lw=1, label=f"z_gate={z_gate} m")
ax[0, 0].axhline(0, color="k", lw=0.8, ls=":")
ax[0, 0].set_ylabel("altitude h [m]"); ax[0, 0].set_title("Altitude (descent profile)")
ax[0, 0].legend(fontsize=8)

# (2) Sink rate (positive = descending)
ax[0, 1].plot(t, d["sink"], color="tab:red")
ax[0, 1].axhline(0, color="k", lw=0.8, ls=":")
ax[0, 1].set_ylabel("sink [m/s] (down +)")
ax[0, 1].set_title(f"Sink rate (touchdown {d['sink'][-1]:.4f} m/s)")

# (3) Aircraft pitch attitude theta (with gamma for context)
ax[1, 0].plot(t, d["theta_deg"], color="tab:blue", label=r"$\theta$ (pitch)")
ax[1, 0].plot(t, d["gamma_deg"], color="tab:green", alpha=0.7, label=r"$\gamma$ (flight path)")
ax[1, 0].axhline(0, color="k", lw=0.8, ls=":")
ax[1, 0].set_ylabel("angle [deg]")
ax[1, 0].set_title(r"Aircraft pitch $\theta$ (impact row needs $\theta>0$)")
ax[1, 0].legend(fontsize=8)

# (4) Predicted peak load factor vs the limit
ax[1, 1].plot(t, d["n_peak"], color="tab:red", label="n_peak [g]")
ax[1, 1].axhline(n_limit, color="k", ls="--", lw=1, label=f"n_limit={n_limit} g")
ax[1, 1].set_ylabel("load factor [g]")
ax[1, 1].set_title(f"Predicted peak CG load (max active = {d['n_peak'][active].max():.3f} g)" if active.any()
                   else "Predicted peak CG load")
ax[1, 1].legend(fontsize=8)

# (5) Approach parameter kappa with the practical band
ax[2, 0].plot(t, d["kappa_imp"], color="tab:purple")
ax[2, 0].axhspan(0.2, 10.0, color="tab:green", alpha=0.08, label="practical [0.2,10]")
ax[2, 0].set_ylabel(r"$\kappa$"); ax[2, 0].set_title("Approach parameter")
ax[2, 0].legend(fontsize=8)

# (6) Barrier value + nested psi (must be >= 0 where active)
ax[2, 1].plot(t, d["b_impact"], label=f"b_impact (min act {d['b_impact'][active].min():.3f})" if active.any() else "b_impact")
ax[2, 1].plot(t, d["psi1_imp"], label="psi1_imp", alpha=0.8)
ax[2, 1].plot(t, d["psi2_imp"], label="psi2_imp", alpha=0.8)
ax[2, 1].axhline(0, color="k", lw=1, ls="--")
ax[2, 1].set_ylabel("value"); ax[2, 1].set_title("Barrier b and nested psi (>=0 = safe)")
ax[2, 1].legend(fontsize=8)

# (7) Elevator: nominal vs filtered (does the barrier shape it?)
ax[3, 0].plot(t, np.rad2deg(d["de_nom"]), color="tab:gray", ls=":", label=r"$\delta_e$ nominal")
ax[3, 0].plot(t, np.rad2deg(d["de"]), color="tab:purple", label=r"$\delta_e$ filtered")
ax[3, 0].set_ylabel("elevator [deg]"); ax[3, 0].set_xlabel("t [s]")
ax[3, 0].set_title("Elevator: nominal vs CBF-filtered"); ax[3, 0].legend(fontsize=8)

# (8) Impact-row residual rhs - a.U (~0 => binding/active-set)
ax[3, 1].plot(t, d["res_imp"], color="tab:brown")
ax[3, 1].axhline(0, color="k", lw=1, ls="--", label="binding (res=0)")
ax[3, 1].set_ylabel("rhs - a.U"); ax[3, 1].set_xlabel("t [s]")
ax[3, 1].set_title("Impact-row residual (slack stays 0 if never binds)")
ax[3, 1].legend(fontsize=8)

for a in ax.flat:
    a.grid(alpha=0.3)
    if active.any():
        shade_active(a)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(out, dpi=110)

# console summary
na = active.sum()
print(f"active window rows: {na}")
if na:
    print(f"  n_peak max  = {d['n_peak'][active].max():.4f} g  (limit {n_limit})")
    print(f"  kappa range = [{d['kappa_imp'][active].min():.3f}, {d['kappa_imp'][active].max():.3f}]")
    print(f"  b_impact min= {d['b_impact'][active].min():.4f}")
    print(f"  res_imp min = {d['res_imp'][active].min():.4f}  ({'BINDS' if d['res_imp'][active].min() < 1e-3 else 'never binds'})")
print(f"  elevator |de - de_nom| max = {np.rad2deg(np.abs(d['de'] - d['de_nom'])).max():.4f} deg")
print("wrote", out)
