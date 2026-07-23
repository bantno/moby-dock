#!/usr/bin/env python3
"""Generate a NACA 4414 viscous post-stall polar (CL/CD/CM vs alpha) across the
landing Reynolds range, using NeuralFoil (a trained Xfoil surrogate -> smooth,
differentiable, includes post-stall). This is the 2D section stall physics that
VSPAERO (inviscid) cannot produce; it gets spliced onto the VSPAERO 3D deck.

Vehicle: NACA 4414, c=0.25 m, V_landing ~ 13-18 m/s  ->  Re ~ 2.2e5 - 3.1e5.
"""
import numpy as np, aerosandbox as asb, neuralfoil as nf

af = asb.Airfoil("naca4414")                       # 4-digit generator
alpha = np.arange(-8.0, 24.01, 0.5)                # resolve through stall
Res   = [150_000, 200_000, 250_000, 300_000, 400_000]

rows = []
for Re in Res:
    a = nf.get_aero_from_airfoil(af, alpha=alpha, Re=Re, model_size="xxlarge")
    CL, CD, CM = np.array(a["CL"]), np.array(a["CD"]), np.array(a["CM"])
    conf = np.array(a["analysis_confidence"])
    imax = int(np.argmax(CL))
    print(f"Re={Re:>7}:  CLmax={CL[imax]:.3f} @ a={alpha[imax]:+.1f} deg   "
          f"CD@CLmax={CD[imax]:.4f}   CM@CLmax={CM[imax]:+.3f}   conf~{conf.mean():.2f}")
    for i,al in enumerate(alpha):
        rows.append((Re, al, CL[i], CD[i], CM[i], conf[i]))

import csv
out = "data/naca4414_polar.csv"
with open(out,"w",newline="") as f:
    w = csv.writer(f); w.writerow(["Re","alpha_deg","CL","CD","CM","confidence"])
    w.writerows([(r[0], f"{r[1]:.2f}", f"{r[2]:.5f}", f"{r[3]:.6f}",
                  f"{r[4]:.5f}", f"{r[5]:.3f}") for r in rows])
print("wrote", out, f"({len(rows)} rows, {len(Res)} Re x {len(alpha)} alpha)")
