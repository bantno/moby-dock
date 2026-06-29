#!/usr/bin/env python3
"""Altitude-sensor study figure: CBF ideal vs CBF filtered (noise+LPF) vs no-CBF
(nominal only) longitudinal landing, as a 2x4 panel grid.

Panels (row 0): altitude | descent rate vs max-allowable envelope | airspeed vs
  velocity limits | thrust vs thrust limits.
Panels (row 1): elevator (nominal vs CBF) | thrust jerk (nominal vs CBF) |
  altitude sensor chain (true/raw/filtered) | all barriers over time.

The three input CSVs are lon_autoland_sim traces for the same scenario with:
  ideal = perfect altitude (h_meas_stddev=0, h_lpf_tau=0),
  filt  = noisy + low-passed altitude (the scenario's h_meas_stddev / h_lpf_tau),
  nocbf = CBF disabled (cbf_enabled=false), nominal-only on ideal altitude.

All panels drop the non-physical sub-surface touchdown frame (first h<=0). The
elevator and thrust-jerk panels ADDITIONALLY clip the near-surface band (h<=clip_h)
to hide the steep-barrier flare chatter; every other panel shows the full descent.

Usage: plot_lon_sensor_panels.py ideal.csv filt.csv nocbf.csv out.png
                                 [sigma] [tau] [Vmin] [Vmax] [Tmax] [clip_h]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ideal = np.genfromtxt(sys.argv[1], delimiter=",", names=True)
filt = np.genfromtxt(sys.argv[2], delimiter=",", names=True)
nocbf = np.genfromtxt(sys.argv[3], delimiter=",", names=True)


# Always drop the non-physical sub-surface touchdown frame (first h<=0; the
# control logged there is never integrated). The ELEVATOR panel additionally
# clips the near-surface band (h<=clip_h) to hide the steep-barrier saturation
# chatter; every other panel shows the full physical descent. (clip_h = 10th arg.)
clip_h = float(sys.argv[10]) if len(sys.argv) > 10 else 0.0


def trim_h(d, thresh):
    idx = np.where(d["h"] <= thresh)[0]
    return d[: idx[0]] if (len(idx) and idx[0] > 0) else d


ideal, filt, nocbf = trim_h(ideal, 0.0), trim_h(filt, 0.0), trim_h(nocbf, 0.0)
ideal_e, filt_e, nocbf_e = trim_h(ideal, clip_h), trim_h(filt, clip_h), trim_h(nocbf, clip_h)
out = sys.argv[4] if len(sys.argv) > 4 else "lon_lpf_panels.png"
sigma = float(sys.argv[5]) if len(sys.argv) > 5 else 0.5
tau = float(sys.argv[6]) if len(sys.argv) > 6 else 0.15
Vmin = float(sys.argv[7]) if len(sys.argv) > 7 else 13.5
Vmax = float(sys.argv[8]) if len(sys.argv) > 8 else 20.0
Tmax = float(sys.argv[9]) if len(sys.argv) > 9 else 50.0

C_ID, C_FI, C_NO = "tab:blue", "tab:green", "tab:red"
cases = [("ideal CBF", ideal, C_ID), ("filtered CBF", filt, C_FI), ("no CBF", nocbf, C_NO)]

fig, ax = plt.subplots(2, 4, figsize=(25, 10))
fig.suptitle(rf"CBF (ideal / filtered) vs no-CBF longitudinal landing  "
             rf"($\sigma$={sigma:g} m, $\tau$={tau:g} s)", fontsize=15, weight="bold")

# (1) Altitude.
for lab, d, c in cases:
    ax[0, 0].plot(d["t"], d["h"], c, label=lab)
ax[0, 0].axhline(0, color="k", lw=0.8, ls=":")
ax[0, 0].set(xlabel="t [s]", ylabel="altitude h [m]", title="Altitude")
ax[0, 0].legend(fontsize=9)

# (2) Descent rate + max allowable (descent-barrier envelope sqrt(v_safe^2 +
#     2 a_brk(V,gamma) h) = b_descent + sink; the CBF keeps sink below it). The
#     envelope is huge at altitude, so cap the y-axis to keep the detail readable.
smax = max(np.nanmax(d["sink"]) for _, d, _ in cases)
for lab, d, c in cases:
    ax[0, 1].plot(d["t"], d["sink"], c, label=lab)
    ax[0, 1].plot(d["t"], d["b_descent"] + d["sink"], c, ls="--", lw=1.0, alpha=0.6)
ax[0, 1].plot([], [], color="0.4", ls="--", lw=1.0, label="max allowable (envelope)")
ax[0, 1].axhline(0, color="k", lw=0.8, ls=":")
ax[0, 1].set(xlabel="t [s]", ylabel="descent rate [m/s] (down +)",
             title="Descent rate vs max allowable", ylim=(-0.6, 1.25 * smax))
ax[0, 1].legend(fontsize=8)

# (3) Airspeed vs lower/upper velocity limits. If the trace carries the altitude-
#     SCHEDULED limits (Vmin_sched/Vmax_sched columns), draw those time-varying
#     envelopes (from the filtered case) instead of constant lines; else fall back
#     to the constant Vmin/Vmax passed on the command line.
for lab, d, c in cases:
    ax[0, 2].plot(d["t"], d["V"], c, label=lab)
if filt.dtype.names is not None and "Vmin_sched" in filt.dtype.names:
    ax[0, 2].plot(filt["t"], filt["Vmin_sched"], color="k", ls=":", lw=1.2,
                  label=r"Vmin(h) sched")
    ax[0, 2].plot(filt["t"], filt["Vmax_sched"], color="k", ls="--", lw=1.2,
                  label=r"Vmax(h) sched")
else:
    ax[0, 2].axhline(Vmin, color="k", ls=":", lw=1.0, label=f"Vmin={Vmin:g}")
    ax[0, 2].axhline(Vmax, color="k", ls="--", lw=1.0, label=f"Vmax={Vmax:g}")
ax[0, 2].set(xlabel="t [s]", ylabel="airspeed V [m/s]", title="Airspeed vs velocity limits")
ax[0, 2].legend(fontsize=8)

# (4) Thrust vs thrust limits [0, Tmax].
for lab, d, c in cases:
    ax[0, 3].plot(d["t"], d["T"], c, label=lab)
ax[0, 3].axhline(0, color="k", ls=":", lw=1.0, label="Tmin=0")
ax[0, 3].axhline(Tmax, color="k", ls="--", lw=1.0, label=f"Tmax={Tmax:g}")
ax[0, 3].set(xlabel="t [s]", ylabel="thrust T [N]", title="Thrust vs thrust limits")
ax[0, 3].legend(fontsize=8)

# (5) Elevator: nominal vs CBF-filtered (CBF cases) + no-CBF applied. ONLY this
#     panel uses the near-surface-clipped data (h>clip_h) to hide flare chatter.
ce = [("ideal", ideal_e, C_ID), ("filtered", filt_e, C_FI)]
for lab, d, c in ce:
    ax[1, 0].plot(d["t"], np.rad2deg(d["de_nom"]), c, ls="--", alpha=0.55, label=f"{lab}: nominal")
    ax[1, 0].plot(d["t"], np.rad2deg(d["de"]), c, label=f"{lab}: CBF")
ax[1, 0].plot(nocbf_e["t"], np.rad2deg(nocbf_e["de"]), C_NO, label="no CBF (=nominal)")
ax[1, 0].set(xlabel="t [s]", ylabel=r"elevator $\delta_e$ [deg]",
             title=f"Elevator: nominal vs CBF (h<{clip_h:g} m clipped)")
ax[1, 0].legend(fontsize=8, ncol=2)

# (6) Thrust jerk (throttle control) Tddot: nominal vs CBF-filtered + no-CBF.
#     Also near-surface-clipped (h>clip_h), like the elevator.
for lab, d, c in ce:
    ax[1, 1].plot(d["t"], d["Tddot_nom"], c, ls="--", alpha=0.55, label=f"{lab}: nominal")
    ax[1, 1].plot(d["t"], d["Tddot"], c, label=f"{lab}: CBF")
ax[1, 1].plot(nocbf_e["t"], nocbf_e["Tddot"], C_NO, label="no CBF (=nominal)")
ax[1, 1].set(xlabel="t [s]", ylabel=r"thrust jerk $\ddot{T}$ [N/s$^2$]",
             title=f"Thrust jerk: nominal vs CBF (h<{clip_h:g} m clipped)")
ax[1, 1].legend(fontsize=8, ncol=2)

# (7) Altitude sensor chain (filtered case): true / raw / filtered.
ax[1, 2].plot(filt["t"], filt["h_meas"], color="0.6", lw=0.7, alpha=0.8, label=r"raw $h_{meas}=h+\mathcal{N}$")
ax[1, 2].plot(filt["t"], filt["h_filt"], C_FI, lw=1.6, label=r"filtered $h_{filt}$ (CBF input)")
ax[1, 2].plot(filt["t"], filt["h"], "k--", lw=1.0, label=r"true $h$")
ax[1, 2].axhline(0, color="k", lw=0.8, ls=":")
ax[1, 2].set(xlabel="t [s]", ylabel="altitude [m]", title="Altitude sensor chain (filtered case)")
ax[1, 2].legend(fontsize=9)

# (8) All barriers over time (filtered case).
t = filt["t"]
bars = {
    "descent": filt["b_descent"],
    "airspeed (V-Vmin)": filt["b_airspeed"],
    "airspeed upper (Vmax-V)": filt["b_airspeed_upper"],
    "thrust min (T)": filt["T"],
    "thrust max (Tmax-T)": Tmax - filt["T"],
    "impact": filt["b_impact"],
}
for name, vals in bars.items():
    ax[1, 3].plot(t, vals, lw=1.2, label=f"{name}  (min {np.nanmin(vals):.2f})")
ax[1, 3].axhline(0, color="k", lw=1.0, ls="--")
ax[1, 3].set(xlabel="t [s]", ylabel="barrier value", title="All barriers (filtered case; >=0 = safe)")
ax[1, 3].legend(fontsize=7)

for a in ax.flat:
    a.grid(alpha=0.3)
fig.tight_layout(rect=[0, 0, 1, 0.98])
fig.savefig(out, dpi=120)
print("wrote", out)
for lab, d, c in cases:
    i = np.where(d["h"] <= 0)[0]
    i = i[0] if len(i) else len(d["h"]) - 1
    print(f"  {lab:14s} touchdown t={d['t'][i]:.2f}s  descent rate={d['sink'][i]:.3f} m/s  "
          f"V={d['V'][i]:.2f}  T={d['T'][i]:.2f}")
