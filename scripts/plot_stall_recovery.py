#!/usr/bin/env python3
"""Overlay CBF-on (recovery) vs CBF-off (departure) stall-demo traces.

Companion to data/lon_stall_recovery_cbf.yaml / data/lon_stall_recovery.yaml:
the same aggressive pitch-up on the NACA 4414 viscous-stall plant, with and
without the safety filter. Shows the stall/AoA HOCBF overriding the nominal's
nose-up command at alpha_lim = alpha_stall - margin (the pilot recovery) while
the unfiltered run departs.

Usage: plot_stall_recovery.py cbf_on.csv cbf_off.csv [out.png] [alpha_stall_deg] [margin_deg]
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

on = np.genfromtxt(sys.argv[1], delimiter=",", names=True)
off = np.genfromtxt(sys.argv[2], delimiter=",", names=True)
out = sys.argv[3] if len(sys.argv) > 3 else "figures/stall_recovery_compare.png"
a_stall = float(sys.argv[4]) if len(sys.argv) > 4 else 11.0
margin = float(sys.argv[5]) if len(sys.argv) > 5 else 2.0
a_lim = a_stall - margin

C_ON, C_OFF = "tab:blue", "tab:red"  # fixed per entity, same on every panel
t_end = max(on["t"][-1], off["t"][-1])
t_dep = off["t"][-1]  # unfiltered log ends here (departure past model validity)

fig, ax = plt.subplots(2, 2, figsize=(13, 8.5))
fig.suptitle("Stall recovery via the CBF-QP — aggressive pitch-up on the "
             "viscous-stall plant", fontsize=14, weight="bold")

# (1) Angle of attack vs the barrier: the recovery itself.
a = ax[0, 0]
a.axhspan(a_stall, 30.0, color="tab:red", alpha=0.06)
a.axhline(a_stall, color="k", lw=1.0, ls="-",
          label=rf"$\alpha_{{stall}}$ = {a_stall:.0f}°")
a.axhline(a_lim, color="k", lw=1.0, ls="--",
          label=rf"$\alpha_{{lim}}$ = {a_lim:.0f}° (barrier)")
a.plot(on["t"], on["alpha_deg"], C_ON, label="CBF on")
a.plot(off["t"], off["alpha_deg"], C_OFF, label="CBF off")
i_ex = np.argmax(off["alpha_deg"] > 28.0)
a.annotate(rf"departs $\to$ {off['alpha_deg'][-1]:.0f}° deep stall",
           xy=(off["t"][i_ex], 27.0), xytext=(off["t"][i_ex] + 2.5, 22.0),
           color=C_OFF, fontsize=9,
           arrowprops=dict(arrowstyle="->", color=C_OFF, lw=1))
a.annotate("rides the limit", xy=(20.0, 10.1), color=C_ON, fontsize=9)
a.set_ylim(-1, 30)
a.set(xlabel="t [s]", ylabel=r"$\alpha$ [deg]",
      title="AoA: the barrier caps the pitch-up at the stall margin")
a.legend(fontsize=8, loc="upper right")

# (2) Altitude: controlled mush vs falling departure.
a = ax[0, 1]
a.plot(on["t"], on["h"], C_ON, label="CBF on")
a.plot(off["t"], off["h"], C_OFF, label="CBF off")
a.plot(off["t"][-1], off["h"][-1], "o", color=C_OFF, ms=5)
a.annotate(f"log ends mid-fall\n(sink {off['sink'][-1]:.1f} m/s)",
           xy=(off["t"][-1], off["h"][-1]),
           xytext=(off["t"][-1] + 2.0, off["h"][-1] + 6.0),
           color=C_OFF, fontsize=9,
           arrowprops=dict(arrowstyle="->", color=C_OFF, lw=1))
a.set(xlabel="t [s]", ylabel="altitude h [m]",
      title="Altitude: recovered flight vs departure")
a.legend(fontsize=8, loc="lower right")

# (3) Elevator: the QP shaves the nominal's nose-up command.
a = ax[1, 0]
de_box = 28.6
a.axhline(de_box, color="k", lw=0.9, ls=":", label=rf"$\pm${de_box}° box (QP)")
a.plot(on["t"], np.degrees(on["de_nom"]), C_ON, ls="--", lw=1.2, alpha=0.4,
       label=r"$\delta_e$ nominal wants (on)")
a.plot(on["t"], np.degrees(on["de"]), C_ON, label=r"$\delta_e$ applied (on)")
a.plot(off["t"], np.degrees(off["de"]), C_OFF,
       label=r"$\delta_e$ applied (off, unfiltered)")
a.set_ylim(0, 32)
a.annotate("nominal demand exceeds panel (clipped)", xy=(t_end - 0.8, 30.6),
           color=C_ON, alpha=0.7, fontsize=8, ha="right")
a.set(xlabel="t [s]", ylabel=r"$\delta_e$ [deg]  (+ = nose-up)",
      title="Elevator: stall row overrides the nose-up demand")
a.legend(fontsize=8, loc="lower right")

# (4) Lift coefficient: attached flow held vs post-stall collapse.
a = ax[1, 1]
a.plot(on["t"], on["CL"], C_ON, label="CBF on")
a.plot(off["t"], off["CL"], C_OFF, label="CBF off")
a.annotate("lift collapse", xy=(off["t"][-1] - 0.3, off["CL"][-1] + 0.06),
           color=C_OFF, fontsize=9, ha="right")
a.set(xlabel="t [s]", ylabel=r"$C_L$",
      title="Realized lift: the recovery keeps the wing attached")
a.legend(fontsize=8, loc="center right")

for a in ax.flat:
    a.grid(alpha=0.3)
    a.set_xlim(0, t_end)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig(out, dpi=120)
i_bind = np.argmax(on["alpha_deg"] > a_lim - 0.5)
print(f"wrote {out}")
print(f"  ON : alpha max {on['alpha_deg'].max():.2f} deg (limit {a_lim}), "
      f"binds from t~{on['t'][i_bind]:.1f} s, min CL {on['CL'].min():.2f}")
print(f"  OFF: alpha max {off['alpha_deg'].max():.1f} deg, min CL {off['CL'].min():.2f}, "
      f"max sink {off['sink'].max():.1f} m/s, log ends t={t_dep:.1f} s")
