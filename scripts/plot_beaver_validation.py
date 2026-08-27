#!/usr/bin/env python3
"""Beaver flight-dynamics validation summary figure.

Panels:
  (1) trimmed-flight elevator curve -de(V) vs the FDC 1.2 fig. 10.13
      reference (solid curve digitized by eye, +/-0.5 deg) and the exact
      ACTRIM check point;
  (2) trim incidence + engine power over the sweep;
  (3) s-plane: linearized modes from the C++ autodiff linearization (x) vs
      the independent Python cross-implementation (o) at three conditions.

Inputs: the CSVs written by `beaver_validation [prefix]`.
Usage: plot_beaver_validation.py [prefix] [out.png]
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import validate_beaver_modes as vbm  # noqa: E402  (the independent EOM)

prefix = sys.argv[1] if len(sys.argv) > 1 else "results/beaver"
out = sys.argv[2] if len(sys.argv) > 2 else "results/beaver_validation.png"

sweep = np.genfromtxt(f"{prefix}_trim_sweep.csv", delimiter=",", names=True)

fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
fig.suptitle("DHC-2 Beaver plant validation vs FDC 1.2 / LR-556 references",
             fontsize=13, weight="bold")

# (1) Elevator trim curve vs FDC fig 10.13. Digitized solid ("low power")
# curve, read off the manual figure to ~+/-0.5 deg.
fdc_fig = np.array([[30, 10.2], [32.5, 8.0], [35, 6.2], [40, 3.6], [45, 1.9],
                    [50, 0.7], [55, -0.3], [60, -1.1]])
for n, style in [(1800.0, "tab:blue"), (2300.0, "tab:orange")]:
    m = sweep["n_rpm"] == n
    ax[0].plot(sweep["V"][m], -sweep["de_deg"][m], color=style,
               label=f"this plant, level trim, n={int(n)}")
ax[0].errorbar(fdc_fig[:, 0], fdc_fig[:, 1], yerr=0.5, fmt="s", ms=4,
               color="k", alpha=0.6, capsize=3,
               label="FDC fig. 10.13 solid curve (digitized)")
ax[0].plot([35], [5.3333], "r*", ms=14,
           label="FDC ACTRIM check case (exact)")
ax[0].set_xlabel("V [m/s]")
ax[0].set_ylabel(r"$-\delta_e$ [deg]")
ax[0].set_title("Trimmed-flight elevator curve")
ax[0].grid(alpha=0.3)
ax[0].legend(fontsize=8)

# (2) Trim incidence and power over the sweep (n=1800).
m = sweep["n_rpm"] == 1800.0
ax[1].plot(sweep["V"][m], sweep["alpha_deg"][m], color="tab:green",
           label=r"$\alpha$ [deg]")
ax[1].set_xlabel("V [m/s]")
ax[1].set_ylabel(r"$\alpha$ [deg]", color="tab:green")
ax[1].grid(alpha=0.3)
ax1b = ax[1].twinx()
ax1b.plot(sweep["V"][m], sweep["P_kW"][m], color="tab:red", label="P [kW]")
ax1b.set_ylabel("engine power [kW]", color="tab:red")
ax[1].set_title("Level trim: incidence + power (n=1800)")
ax[1].axvspan(35, 55, color="tab:blue", alpha=0.06)
ax[1].text(45, ax[1].get_ylim()[1] * 0.92, "LR-556 validity band",
           ha="center", fontsize=8, color="tab:blue")

# (3) s-plane: C++ eigenvalues (from the dumped A) vs the independent Python
# implementation, all three conditions.
IDX8 = [0, 1, 2, 3, 4, 5, 7, 6]
colors = {"cruise_45": "tab:blue", "check_35": "tab:orange",
          "approach_35": "tab:green"}
for cond, col in colors.items():
    import csv as csvmod
    with open(f"{prefix}_trim_{cond}.csv") as f:
        t = {k: float(v) for k, v in next(iter(csvmod.DictReader(f))).items()}
    prm = dict(pz=t["pz"], n_rpm=t["n_rpm"], rho=t["rho"], g=t["g"],
               flap=t["flap"])
    s0 = np.array([t["V"], t["alpha"], t["beta"], 0, 0, 0, t["theta"], 0])
    ctl = [t["de"], t["da"], t["dr"]]
    A_py = vbm.complex_step_jacobian(lambda s: vbm.xdot(s, ctl, prm), s0)
    e_py = np.linalg.eigvals(A_py)
    A_cpp = vbm.load_matrix(f"{prefix}_A_{cond}.csv")[np.ix_(IDX8, IDX8)]
    e_cpp = np.linalg.eigvals(A_cpp)
    ax[2].plot(e_cpp.real, e_cpp.imag, "x", ms=9, color=col,
               label=f"{cond} (C++ autodiff)")
    ax[2].plot(e_py.real, e_py.imag, "o", ms=5, mfc="none", color=col,
               label=f"{cond} (independent Python)")
ax[2].axhline(0, color="k", lw=0.6)
ax[2].axvline(0, color="k", lw=0.6)
ax[2].set_xlabel("Re [1/s]")
ax[2].set_ylabel("Im [rad/s]")
ax[2].set_title("Linearized modes: cross-implementation")
ax[2].grid(alpha=0.3)
ax[2].legend(fontsize=7)

fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(out, dpi=140)
print("wrote", out)
