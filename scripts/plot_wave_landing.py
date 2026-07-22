#!/usr/bin/env python3
"""Wave-landing trace: the wave-blind CBF landing onto a JONSWAP/Airy seaway.

Companion to data/lon_landing_waves_lake.yaml (and the regular-wave variant):
the flagship -3 deg approach with the surface-wave field enabled as plant-side
truth. Shows the aircraft flying the surface-relative radar altimeter down onto
the moving surface, and the gap the smooth-water assumption leaves: the
flat-referenced sink meets the spec while the SURFACE-relative closure (and the
wave-referenced TN 1516 slam load printed by the sim) explodes on a rising face.

Usage: plot_wave_landing.py waves.csv [flat.csv] [out.png]
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

wv = np.genfromtxt(sys.argv[1], delimiter=",", names=True)
flat = None
out = "figures/lon_landing_waves.png"
args = sys.argv[2:]
if args and args[0].endswith(".csv"):
    flat = np.genfromtxt(args[0], delimiter=",", names=True)
    args = args[1:]
if args:
    out = args[0]

t, h, eta = wv["t"], wv["h"], wv["eta"]
x, slope, hf = wv["x"], wv["eta_slope"], wv["h_filt"]
clr = h - eta                      # true clearance over the instantaneous surface
t_td = t[-1]
# Surface-relative closure -d/dt(h - eta), numerically from the trace.
dt = np.diff(t)
closure = np.empty_like(t)
closure[1:] = -np.diff(clr) / dt
closure[0] = closure[1]

C_W, C_S, C_F = "tab:blue", "tab:cyan", "tab:red"
fig, ax = plt.subplots(2, 2, figsize=(13, 8.5))
fig.suptitle("Landing on waves — wave-blind CBF, plant-side JONSWAP/Airy sea "
             "(surface-relative altimeter)", fontsize=14, weight="bold")

# (1) Landing profile over ground: the corrugated surface the keel actually met
# (last ~100 m of track, where the wavelength-scale structure is visible).
a = ax[0, 0]
m = x > x[-1] - 100.0
a.plot(x[m], h[m], C_W, label="keel altitude h")
a.plot(x[m], eta[m], C_S, lw=1.2, label=r"surface $\eta(x)$ under the aircraft")
a.fill_between(x[m], eta[m], eta[m].min() - 0.1, color=C_S, alpha=0.15)
a.plot(x[-1], eta[-1], "v", color="k", ms=7, zorder=5, label="touchdown")
a.annotate(rf"contact on a {np.degrees(np.arctan(slope[-1])):.1f}° face"
           f"\n$\\eta$ = {eta[-1]:+.2f} m",
           xy=(x[-1], eta[-1]), xytext=(x[-1] - 55, eta[-1] + 0.8),
           fontsize=9, arrowprops=dict(arrowstyle="->", lw=1))
a.set(xlabel="ground track x [m]", ylabel="altitude [m]",
      title="Final approach profile: flying down onto the seaway")
a.legend(fontsize=8, loc="upper left")

# (2) What the filter flew: measured clearance vs the true one.
a = ax[0, 1]
m = t > t_td - 12.0
a.plot(t[m], clr[m], C_W, label=r"true clearance $h-\eta$")
a.plot(t[m], hf[m], "k", lw=0.9, alpha=0.6, label=r"altimeter $h_{filt}$ (noisy, LPF)")
a.axhline(0.0, color=C_S, lw=1.0)
a.set(xlabel="t [s]", ylabel="clearance [m]",
      title="Radar-altimeter chain: the filter sees the surface-relative altitude")
a.legend(fontsize=8, loc="upper right")

# (3) The smooth-water gap: flat-referenced sink vs surface-relative closure.
# Only the last couple of encounter periods matter -- aloft, the "closure"
# is just the surface heaving far below the keel.
a = ax[1, 0]
m = t > t_td - 2.5
a.plot(t[m], wv["sink"][m], C_W, label="inertial sink $-\\dot h$ (flat-water reference)")
a.plot(t[m], closure[m], C_F, lw=1.2,
       label=r"surface closure $-\frac{d}{dt}(h-\eta)$ (what the keel feels)")
if flat is not None:
    mf = flat["t"] > flat["t"][-1] - 2.5
    a.plot(flat["t"][mf] - flat["t"][-1] + t_td, flat["sink"][mf], C_W, ls="--",
           lw=1.0, alpha=0.5, label="sink, flat-water run")
a.plot(t_td, closure[-2], "o", color=C_F, ms=5)
a.annotate(f"closure at contact ~{closure[-2]:.1f} m/s\n"
           f"(flat-referenced sink {wv['sink'][-1]:.2f} m/s)",
           xy=(t_td, closure[-2]), xytext=(t_td - 1.1, closure[-2] + 3.2),
           color=C_F, fontsize=9,
           arrowprops=dict(arrowstyle="->", color=C_F, lw=1))
a.set(xlabel="t [s]", ylabel="[m/s]",
      title="Sink: the flat-water spec is met while the face closes 10x faster")
a.legend(fontsize=8, loc="upper left")

# (4) Attitude vs the local surface: tau is referenced to the WRONG plane.
# (Hull-length mean slope; final seconds, one dot per surface the keel passed.)
a = ax[1, 1]
m = t > t_td - 2.5
surf_deg = np.degrees(np.arctan(slope[m]))
a.plot(t[m], wv["theta_deg"][m], C_W, label=r"$\theta$ (attitude, flat-water frame)")
a.plot(t[m], surf_deg, C_S, lw=1.0,
       label=r"surface tilt over the hull length")
a.plot(t[m], wv["theta_deg"][m] - surf_deg, C_F, lw=1.0, alpha=0.8,
       label=r"$\theta$ relative to the surface")
a.axhline(0.0, color="k", lw=0.7, alpha=0.5)
th_rel_td = wv["theta_deg"][-1] - np.degrees(np.arctan(slope[-1]))
a.annotate(f"nose-up on flat water ({wv['theta_deg'][-1]:.1f}°),\n"
           f"nose-INTO the wave face ({th_rel_td:.1f}°)",
           xy=(t_td, th_rel_td), xytext=(t_td - 1.35, 14.5),
           color=C_F, fontsize=9,
           arrowprops=dict(arrowstyle="->", color=C_F, lw=1))
a.set(xlabel="t [s]", ylabel="[deg]",
      title="Attitude: the nose-up touchdown is nose-down relative to the face")
a.legend(fontsize=8, loc="lower left")

for a in ax.flat:
    a.grid(alpha=0.3)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig(out, dpi=120)
print(f"wrote {out}")
print(f"  touchdown t={t_td:.2f} s  eta={eta[-1]:+.3f} m  "
      f"face={np.degrees(np.arctan(slope[-1])):.1f} deg  "
      f"sink(flat)={wv['sink'][-1]:.2f}  closure={closure[-1]:.2f} m/s")
