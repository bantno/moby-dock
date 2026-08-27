#!/usr/bin/env python3
"""Beaver flight-dynamics validation figures.

Figure 1 (<out>): three panels --
  (1) trimmed-flight elevator curve: our level trim, our fixed-pz=20 "Hg
      descending trim, the pixel-digitized FDC fig. 10.13 solid (low-power)
      curve (data/fdc_fig1013_solid_digitized.csv), and the exact ACTRIM
      check point;
  (2) trim incidence + engine power over the sweep;
  (3) s-plane: C++ autodiff modes (x) vs the independent Python
      cross-implementation (o) at three conditions.

Figure 2 (<out base>_checkcase.png): the FDC ACTRIM check-case oracle --
  per-row |xdot difference| vs the printout (and vs the WRONG 2000-ft
  atmosphere, showing the check's discrimination), plus the trim-recovery
  comparison.

Figure 3 (<out base>_fig1013_sidebyside.png): the scanned manual figure next
  to our model curves (needs figures/fdc_fig1013_scan.png; skipped if absent).

Inputs: the CSVs written by `beaver_validation [prefix]`.
Usage: plot_beaver_validation.py [prefix] [out.png]
"""
import csv as csvmod
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import validate_beaver_modes as vbm  # noqa: E402  (the independent EOM)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
prefix = sys.argv[1] if len(sys.argv) > 1 else "results/beaver"
out = sys.argv[2] if len(sys.argv) > 2 else "results/beaver_validation.png"
base = out[:-4] if out.endswith(".png") else out

kDeg = np.pi / 180.0


def trim_fixed_pz(V, pz, n=1800.0, rho=1.225, g=9.80665, flap=0.0):
    """Wings-level trim at FIXED manifold pressure (gamma free): unknowns
    [alpha, beta, de, da, dr, gamma]; Newton with complex-step Jacobian on the
    independent-implementation EOM."""
    prm = dict(pz=pz, n_rpm=n, rho=rho, g=g, flap=flap)

    def resid(z):
        al, be, de, da, dr, gam = z
        th = al + np.arcsin(np.sin(gam) / np.cos(be))
        s = np.array([V + 0j, al, be, 0, 0, 0, th, 0], dtype=complex)
        return vbm.xdot(s, [de, da, dr], prm)[:6]

    z = np.array([0.05, 0.0, -0.05, 0.0, 0.0, -0.02], dtype=complex)
    for _ in range(40):
        R = resid(z).real
        if np.max(np.abs(R)) < 1e-11:
            break
        J = np.zeros((6, 6))
        for j in range(6):
            zp = z.copy()
            zp[j] += 1e-30j
            J[:, j] = np.imag(resid(zp)) / 1e-30
        z = z - np.linalg.solve(J, R)
    return z.real


def load_trim(cond):
    with open(f"{prefix}_trim_{cond}.csv") as f:
        return {k: float(v) for k, v in next(iter(csvmod.DictReader(f))).items()}


sweep = np.genfromtxt(f"{prefix}_trim_sweep.csv", delimiter=",", names=True)
dig = np.genfromtxt(os.path.join(REPO, "data/fdc_fig1013_solid_digitized.csv"),
                    delimiter=",", names=True, skip_header=4)

# ============================== Figure 1 ====================================
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
fig.suptitle("DHC-2 Beaver plant validation vs FDC 1.2 / LR-556 references",
             fontsize=13, weight="bold")

m = sweep["n_rpm"] == 1800.0
ax[0].plot(sweep["V"][m], -sweep["de_deg"][m], color="tab:blue",
           label="this plant: level trim (n=1800)")
Vfix = np.arange(31.0, 60.1, 1.0)
defix = [-np.degrees(trim_fixed_pz(v, 20.0)[2]) for v in Vfix]
ax[0].plot(Vfix, defix, color="tab:green", ls="--",
           label='this plant: fixed pz=20 "Hg trim')
ax[0].plot(dig["V"], dig["minus_de_deg"], color="k", lw=1, alpha=0.65,
           label="FDC fig. 10.13 solid curve (digitized scan)")
ax[0].plot([35], [5.3333], "r*", ms=14, zorder=5,
           label="FDC ACTRIM check case (exact, matched to 0.0004$^\\circ$)")
ax[0].set_xlabel("V [m/s]")
ax[0].set_ylabel(r"$-\delta_e$ [deg]")
ax[0].set_title("Trimmed-flight elevator curve")
ax[0].grid(alpha=0.3)
ax[0].legend(fontsize=7.5)

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

IDX8 = [0, 1, 2, 3, 4, 5, 7, 6]
colors = {"cruise_45": "tab:blue", "check_35": "tab:orange",
          "approach_35": "tab:green"}
for cond, col in colors.items():
    t = load_trim(cond)
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

# ============================== Figure 2 ====================================
# The check-case oracle, straight from the FDC printout (figs. 10.18/10.19).
FDC = dict(V=35.0, alpha=2.1131e-1, beta=-2.0667e-2, theta=1.9190e-1,
           de=-9.3083e-2, da=9.6242e-3, dr=-4.9506e-2, pz=20.0, n=1800.0)
XDOT_FDC = np.array([-1.8871e-4, -1.2348e-5, 4.6356e-4,
                     -2.5027e-5, -2.0660e-5, -5.2604e-5])
ROWS = [r"$\dot V$", r"$\dot\alpha$", r"$\dot\beta$",
        r"$\dot p$", r"$\dot q$", r"$\dot r$"]


def eval_case(rho):
    prm = dict(pz=FDC["pz"], n_rpm=FDC["n"], rho=rho, g=9.80665, flap=0.0)
    s = np.array([FDC["V"], FDC["alpha"], FDC["beta"], 0, 0, 0,
                  FDC["theta"], 0], dtype=complex)
    return vbm.xdot(s, [FDC["de"], FDC["da"], FDC["dr"]], prm)[:6].real


fig2, bx = plt.subplots(1, 2, figsize=(13, 4.6))
fig2.suptitle("FDC ACTRIM check case: external oracle for the nonlinear model",
              fontsize=13, weight="bold")

d_sl = np.abs(eval_case(1.225) - XDOT_FDC)
d_2k = np.abs(eval_case(1.1549) - XDOT_FDC)
xpos = np.arange(6)
bx[0].bar(xpos - 0.18, d_2k, 0.36, color="lightgray", edgecolor="gray",
          label=r"wrong atmosphere: $\rho$(2000 ft)")
bx[0].bar(xpos + 0.18, d_sl, 0.36, color="tab:blue",
          label=r"sea-level $\rho$ = 1.225 (ours)")
bx[0].axhline(1e-5, color="tab:red", ls="--", lw=1)
bx[0].text(2.5, 1.35e-5, "FDC printout precision (5 sig. digits)",
           color="tab:red", fontsize=8, ha="center")
bx[0].set_yscale("log")
bx[0].set_ylim(1e-8, 1.0)
bx[0].set_xticks(xpos, ROWS)
bx[0].set_ylabel(r"| ours $-$ FDC printed |  [SI/s]")
bx[0].set_title("State-derivative match at the printed (x, u)")
bx[0].grid(alpha=0.3, axis="y")
bx[0].legend(fontsize=8)

# Trim recovery: FDC printed trim vs our converged 6-axis Newton trim
# (values from beaver_validation / test_beaver_dynamics.cpp).
names = [r"$\alpha$", r"$\theta$", r"$\delta_e$", r"$\delta_a$",
         r"$\delta_r$", r"$\beta$"]
fdc_v = np.array([12.1072, 10.9951, -5.3333, 0.5514, -2.8365, -1.1841])
ours_v = np.array([12.1050, 10.9930, -5.3337, 0.4630, -2.8194, -1.0155])
ypos = np.arange(len(names))[::-1]
bx[1].plot(fdc_v, ypos, "s", ms=9, mfc="none", color="k",
           label="FDC printed trim (fmins, $\\dot\\beta$ residual 4.6e-4)")
bx[1].plot(ours_v, ypos, "x", ms=9, color="tab:blue",
           label="our Newton trim (residual ~1e-14)")
for y, a, b in zip(ypos, fdc_v, ours_v):
    bx[1].plot([a, b], [y, y], color="tab:blue", lw=0.8, alpha=0.5)
    bx[1].annotate(f"$\\Delta$={abs(a-b):.3f}$^\\circ$", (max(a, b) + 0.4, y),
                   fontsize=7.5, va="center", color="gray")
bx[1].set_yticks(ypos, names)
bx[1].set_xlabel("[deg]")
bx[1].set_xlim(-7, 16)
bx[1].set_title("Trim recovery at the check-case condition")
bx[1].grid(alpha=0.3, axis="x")
bx[1].legend(fontsize=8, loc="lower right")

fig2.tight_layout(rect=(0, 0, 1, 0.93))
fig2.savefig(base + "_checkcase.png", dpi=140)
print("wrote", base + "_checkcase.png")

# ============================== Figure 3 ====================================
scan = os.path.join(REPO, "figures/fdc_fig1013_scan.png")
if os.path.exists(scan):
    fig3, cx = plt.subplots(1, 2, figsize=(14, 5.2))
    fig3.suptitle("Trim curve vs the FDC 1.2 manual, side by side",
                  fontsize=13, weight="bold")
    cx[0].imshow(plt.imread(scan), cmap="gray")
    cx[0].set_axis_off()
    cx[0].set_title("FDC 1.2 manual fig. 10.13 (scan)\n"
                    "solid = low power, dotted = high power, o/x = flight test",
                    fontsize=9)
    cx[1].plot(sweep["V"][m], -sweep["de_deg"][m], color="tab:blue",
               label="this plant: level trim (n=1800)")
    cx[1].plot(Vfix, defix, color="tab:green", ls="--",
               label='this plant: fixed pz=20 "Hg trim')
    cx[1].plot(dig["V"], dig["minus_de_deg"], color="k", lw=1, alpha=0.65,
               label="scan solid curve (digitized)")
    cx[1].plot([35], [5.3333], "r*", ms=14, zorder=5,
               label="ACTRIM check case (exact)")
    cx[1].set_xlim(30, 60)
    cx[1].set_ylim(-2, 12)
    cx[1].set_xlabel("V [m/s]")
    cx[1].set_ylabel(r"$-\delta_e$ [deg]")
    cx[1].set_title("This plant, same axes", fontsize=9)
    cx[1].grid(alpha=0.3)
    cx[1].legend(fontsize=8)
    fig3.tight_layout(rect=(0, 0, 1, 0.92))
    fig3.savefig(base + "_fig1013_sidebyside.png", dpi=140)
    print("wrote", base + "_fig1013_sidebyside.png")
