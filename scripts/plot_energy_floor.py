#!/usr/bin/env python3
"""Overlay the powered stall recovery (energy floor ON) vs the unpowered one.

Companion to data/lon_stall_recovery_cbf_efloor.yaml: the same energy-deficient
pull-up at near-idle thrust, with the opt-in energy-floor row enabled vs
disabled (both runs have the full CBF-QP on, so the stall/AoA row caps alpha in
both). The floor adds the THROTTLE half of the pilot recovery: thrust ramps up
to hold V >= V_floor and the aircraft climbs away instead of mushing down.

Usage: plot_energy_floor.py floor_on.csv floor_off.csv [out.png] [V_floor] [V_stall_speed]
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

on = np.genfromtxt(sys.argv[1], delimiter=",", names=True)
off = np.genfromtxt(sys.argv[2], delimiter=",", names=True)
out = sys.argv[3] if len(sys.argv) > 3 else "figures/stall_recovery_efloor.png"
v_floor = float(sys.argv[4]) if len(sys.argv) > 4 else 12.0
v_stall = float(sys.argv[5]) if len(sys.argv) > 5 else 9.7  # sqrt(2W/(rho S CLmax))

C_ON, C_OFF = "tab:blue", "tab:orange"  # fixed per entity on every panel
t_end = max(on["t"][-1], off["t"][-1])

fig, ax = plt.subplots(2, 2, figsize=(13, 8.5))
fig.suptitle("Powered stall recovery — opt-in energy floor adds the throttle-up "
             "(idle pull-up, CBF on in both runs)", fontsize=13.5, weight="bold")

# (1) Airspeed vs the floor: the sag the floor exists to catch.
a = ax[0, 0]
a.axhline(v_floor, color="k", lw=1.0, ls="--", label=f"V_floor = {v_floor:.0f} m/s")
a.axhline(v_stall, color="k", lw=1.0, ls=":",
          label=f"V_stall(CLmax) ≈ {v_stall:.1f} m/s")
a.plot(on["t"], on["V"], C_ON, label="floor on")
a.plot(off["t"], off["V"], C_OFF, label="floor off")
a.annotate("held at the floor", xy=(22.0, 12.35), color=C_ON, fontsize=9)
a.annotate("sags below stall speed\n(unpowered mush)", xy=(24.0, 8.6),
           color=C_OFF, fontsize=9)
a.set_ylim(7, 15)
a.set(xlabel="t [s]", ylabel="V [m/s]", title="Airspeed: the floor catches the energy sag")
a.legend(fontsize=8, loc="upper right")

# (2) Thrust: the throttle-up itself (the channel the stall row cannot command).
a = ax[0, 1]
a.plot(on["t"], on["T"], C_ON, label="floor on")
a.plot(off["t"], off["T"], C_OFF, label="floor off")
a.axhline(0.3, color="k", lw=0.9, ls=":", label="nominal T_set = 0.3 N (idle)")
a.annotate("QP throttles up", xy=(12.0, 12.5), color=C_ON, fontsize=9)
a.annotate("thrust stays at idle", xy=(20.0, 2.1), color=C_OFF, fontsize=9)
a.set(xlabel="t [s]", ylabel="thrust T [N]",
      title="Thrust: energy-floor row commands the power")
a.legend(fontsize=8, loc="center right")

# (3) AoA: the pitch half is the stall row's, identical in both runs.
a = ax[1, 0]
a.axhspan(11.0, 14.0, color="tab:red", alpha=0.06)
a.axhline(11.0, color="k", lw=1.0, ls="-", label=r"$\alpha_{stall}$ = 11°")
a.axhline(9.0, color="k", lw=1.0, ls="--", label=r"$\alpha_{lim}$ = 9° (barrier)")
a.plot(on["t"], on["alpha_deg"], C_ON, label="floor on")
a.plot(off["t"], off["alpha_deg"], C_OFF, label="floor off")
a.set_ylim(-1, 14)
a.set(xlabel="t [s]", ylabel=r"$\alpha$ [deg]",
      title="AoA: the stall row still caps the pitch-up in both")
a.legend(fontsize=8, loc="lower right", ncol=2)

# (4) Altitude: climb-away vs mush-down.
a = ax[1, 1]
a.plot(on["t"], on["h"], C_ON, label="floor on")
a.plot(off["t"], off["h"], C_OFF, label="floor off")
a.annotate("climbs away under power", xy=(20.0, 120.0), color=C_ON, fontsize=9,
           rotation=32)
a.annotate("mushes down", xy=(26.0, 24.0), color=C_OFF, fontsize=9)
a.set(xlabel="t [s]", ylabel="altitude h [m]",
      title="Altitude: powered recovery vs unpowered mush")
a.legend(fontsize=8, loc="upper left")

for a in ax.flat:
    a.grid(alpha=0.3)
    a.set_xlim(0, t_end)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig(out, dpi=120)
print(f"wrote {out}")
print(f"  floor ON : V min {on['V'].min():.1f}, T max {on['T'].max():.1f} N, "
      f"alpha max {on['alpha_deg'].max():.2f} deg, final h {on['h'][-1]:.0f} m")
print(f"  floor OFF: V min {off['V'].min():.1f}, T max {off['T'].max():.1f} N, "
      f"alpha max {off['alpha_deg'].max():.2f} deg, final h {off['h'][-1]:.0f} m")
