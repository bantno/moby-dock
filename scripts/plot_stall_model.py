#!/usr/bin/env python3
"""Validation overlay for the NACA 4414 stall model: plot the spliced 3D PLANT
CL/CD/CM(alpha) against its two inputs -- the 2D 4414 viscous polar and the
inviscid VSPAERO 3D deck -- so the splice can be eyeballed.

    python3 scripts/plot_stall_model.py        # -> figures/stall_model_check.png
"""
import importlib.util
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(__file__)
spec = importlib.util.spec_from_file_location(
    "st", os.path.join(HERE, "precompute_stall_table.py"))
st = importlib.util.module_from_spec(spec)
spec.loader.exec_module(st)

d = st.build_table()
ag = d["ag_deg"]
w = d["w"]

# Plant curves (what LonDrift computes): (1-w)*VSPAERO + w*Viterna_post.
cl_pl = (1 - w) * d["cl_vsp"] + w * d["CLpost"]
cd_pl = (1 - w) * d["cd_vsp"] + w * d["CDpost"]
cm_pl = (1 - w) * d["cm_vsp"] + w * d["CMpost"]

im = int(np.argmax(cl_pl))
fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
for a, (yp, yv, ttl) in zip(ax, [
        (cl_pl, d["cl_vsp"], "CL"),
        (cd_pl, d["cd_vsp"], "CD"),
        (cm_pl, d["cm_vsp"], "CM (pitch)")]):
    a.plot(ag, yv, "k:", lw=1, label="VSPAERO (inviscid, thrown out post-stall)")
    a.plot(ag, yp, "C3", lw=2, label="PLANT (VSPAERO -> Viterna blend)")
    a.set_xlabel("alpha [deg]"); a.set_ylabel(ttl); a.grid(alpha=.3)
    a.set_xlim(-5, 90); a.legend(fontsize=8)
ax[0].axvline(ag[im], color="C3", ls=":", lw=.8)
ax[0].set_ylim(-0.2, 2.4)
ax[0].set_title(f"Lift -- CLmax={cl_pl[im]:.2f} @ {ag[im]:.0f} deg, craters to 0 at 90 deg")
ax[1].set_title("Drag -- rises to flat-plate CD_max"); ax[2].set_title("Pitch -- stall break")
fig.suptitle("NACA 4414 viscous stall: VSPAERO (attached) handed off to a Viterna "
             f"flat-plate post-stall curve (Re={st.RE_ANCHOR:.0f}, A_stall={st.A_STALL_DEG} deg)")
fig.tight_layout()
out = os.path.join(HERE, os.pardir, "figures", "stall_model_check.png")
fig.savefig(out, dpi=110)
print("wrote", os.path.relpath(out, os.path.join(HERE, os.pardir)))
