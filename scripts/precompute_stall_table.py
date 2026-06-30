#!/usr/bin/env python3
"""Precompute the viscous-stall plant aero (CL/CD/CM vs alpha) for the
longitudinal plant, and emit it as a constexpr C++ table.

WHY: the plant aero comes from data/AHAB_combined.stab, a VSPAERO (inviscid
vortex-lattice) deck whose lift curve is dead-linear to +-20deg -- it physically
cannot stall (no viscous separation), and its post-stall values are MEANINGLESS.
So above stall onset we throw VSPAERO out entirely and hand the lift over to a
viscous + flat-plate post-stall model.

MODEL (blend, not additive):
    C_plant(a) = (1 - w(a)) * C_vspaero(a)  +  w(a) * C_post(a)
  * w(a) is a smoothstep blend weight: 0 below the wing stall onset (attached
    flow == the trustworthy VSPAERO deck), ramping to 1 across the onset, then 1
    (VSPAERO discarded). The C++ plant computes C_vspaero itself and applies this
    blend, so above onset its inviscid values are weighted to zero.
  * C_post(a) is the ABSOLUTE post-stall curve from the **Viterna flat-plate
    extrapolation** (Viterna & Corrigan; the standard wind-turbine/BEM method),
    anchored at the wing stall point so it is continuous with the attached curve
    and declines to CL=0 at 90deg like a real fully-separated wing. The 4414 2D
    viscous polar (NeuralFoil) sets where/how high the wing stalls (CLmax,
    gentleness); Viterna provides the physical deep-stall tail.

So: VSPAERO below onset, Viterna flat-plate above -- lift craters off CLmax then
follows sin(2a)/sin^2(a) to zero at 90deg; drag rises toward CD_max ~ 1.2.

The stall point is a TUNABLE engineering value (NeuralFoil inherits Xfoil's
low-Re/post-stall uncertainty); knobs below. See documentation/stall_model_spec.md.

Usage:
    python3 scripts/precompute_stall_table.py            # regenerate the header
    python3 scripts/precompute_stall_table.py --check     # verify + diagnostics, no write
"""
import csv
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(__file__)
ROOT = os.path.join(HERE, os.pardir)
POLAR_CSV = os.path.join(ROOT, "data", "naca4414_polar.csv")
STAB = os.path.join(ROOT, "data", "AHAB_combined.stab")
HEADER = os.path.join(ROOT, "include", "autoland", "naca4414_stall_table.hpp")

# --- Tunable knobs ------------------------------------------------------------
RE_ANCHOR = 250_000.0       # 2D polar block (landing Re ~ 2.2-3.1e5); sets CLmax
AR = 1.7 ** 2 / 0.425       # aspect ratio b^2/S ~= 6.8 (sets Viterna CD_max)
PARASITE_CD0 = 0.030        # hull/viscous parasite drag (matches aircraft.yaml)
A_STALL_DEG = 11.0          # WING stall angle [geometric deg]: where the lift
                            # leaves the VSPAERO line. CLmax = CL_vspaero(A_STALL)
                            # (kept <= the 2D section CLmax; a tailed vehicle's
                            # CLmax is typically below the wing-section value).
BLEND_HALF_DEG = 2.5        # half-width of the VSPAERO->Viterna handoff (rounds CLmax)
CM_DEEP = -0.45             # flat-plate pitching moment in deep stall (nose-down, cp aft)
SEVERITY = 1.0              # depth multiplier (also a runtime scale; 0 => no stall)

# --- Output geometric-alpha grid (radians, ascending), out to 90deg -----------
DEG = np.pi / 180.0
GRID_DEG = np.concatenate([np.arange(-5.0, 30.0, 0.5),   # fine through stall
                           np.arange(30.0, 90.0001, 2.0)])  # coarse deep-stall tail


def load_stab_attached():
    """(alpha_deg, CL, CD, CM) base-aero from the AHAB .stab, lowest Mach, beta=0
    -- the VSPAERO attached reference (valid only below stall)."""
    txt = open(STAB).read()
    out = []
    for b in re.split(r"\*{5,}", txt):
        am = re.search(r"AoA_\s+(-?\d+\.\d+)", b)
        bm = re.search(r"Beta_\s+(-?\d+\.\d+)", b)
        mm = re.search(r"Mach_\s+(-?\d+\.\d+)", b)
        base = re.search(r"Base_Aero\s+\S+\s+\S+\s+(.+)", b)
        if not (am and bm and mm and base) or abs(float(bm.group(1))) > 1e-6:
            continue
        v = base.group(1).split()  # CFx CFy CFz CMx CMy CMz CL CD CS CMl CMm CMn
        out.append((float(mm.group(1)), float(am.group(1)),
                    float(v[6]), float(v[7]), float(v[10])))
    m0 = min(o[0] for o in out)
    sel = sorted([o for o in out if abs(o[0] - m0) < 1e-9], key=lambda x: x[1])
    return (np.array([s[1] for s in sel]), np.array([s[2] for s in sel]),
            np.array([s[3] for s in sel]), np.array([s[4] for s in sel]))


def polar_clmax_2d():
    """2D NACA 4414 CLmax at the anchor Re (informs the wing CLmax sanity check)."""
    cl = [float(r["CL"]) for r in csv.DictReader(open(POLAR_CSV))
          if abs(float(r["Re"]) - RE_ANCHOR) < 1.0]
    return max(cl) if cl else float("nan")


def viterna(a_stall, cl_stall, cd_stall, ar):
    """Viterna-Corrigan flat-plate post-stall extrapolation. Returns cl(a), cd(a)
    valid for a >= a_stall (radians), continuous with (cl_stall, cd_stall) at
    a_stall and -> (0, CD_max) at 90deg. Standard AeroDyn/BEM method."""
    CDmax = 1.11 + 0.018 * ar
    sa, ca = np.sin(a_stall), np.cos(a_stall)
    B2 = (cd_stall - CDmax * sa * sa) / ca
    A2 = (cl_stall - CDmax * sa * ca) * sa / (ca * ca)

    def cl(a):
        s = np.maximum(np.sin(a), 1e-6)
        return 0.5 * CDmax * np.sin(2 * a) + A2 * np.cos(a) ** 2 / s

    def cd(a):
        return CDmax * np.sin(a) ** 2 + B2 * np.cos(a)
    return cl, cd, CDmax


def build_table():
    sa_deg, scl, scd, scm = load_stab_attached()
    a_stall = A_STALL_DEG * DEG
    ag = GRID_DEG * DEG

    # VSPAERO attached curves at every grid alpha (valid below onset; the C++
    # plant recomputes these itself -- here only for the anchor + the check/plot).
    cl_vsp = np.interp(GRID_DEG, sa_deg, scl)
    cd_vsp = np.interp(GRID_DEG, sa_deg, scd) + PARASITE_CD0
    cm_vsp = np.interp(GRID_DEG, sa_deg, scm)

    # Wing stall anchor = the VSPAERO value at A_STALL (last trustworthy point).
    cl_s = float(np.interp(A_STALL_DEG, sa_deg, scl))
    cd_s = float(np.interp(A_STALL_DEG, sa_deg, scd)) + PARASITE_CD0
    cm_s = float(np.interp(A_STALL_DEG, sa_deg, scm))
    cl_v, cd_v, CDmax = viterna(a_stall, cl_s, cd_s, AR)

    # Absolute post-stall curves (VSPAERO-independent). Below a_stall the Viterna
    # form is only used inside the blend, where w is small; clamp its sin
    # singularity by evaluating at >= a_stall.
    a_eval = np.maximum(ag, a_stall)
    CLpost = cl_v(a_eval)
    CDpost = cd_v(a_eval)
    # Pitching moment: hold the attached value at the stall, move smoothly to a
    # nose-down flat-plate value (cp aft) by ~45deg. Simple, documented approx.
    t = np.clip((GRID_DEG - A_STALL_DEG) / (45.0 - A_STALL_DEG), 0.0, 1.0)
    CMpost = cm_s + (CM_DEEP - cm_s) * (t * t * (3 - 2 * t))

    # Blend weight w: smoothstep 0->1 across [A_STALL -/+ BLEND_HALF]; 0 below
    # (pure VSPAERO), 1 above (pure Viterna -- VSPAERO discarded).
    w = np.clip((GRID_DEG - (A_STALL_DEG - BLEND_HALF_DEG)) / (2 * BLEND_HALF_DEG),
                0.0, 1.0)
    w = w * w * (3 - 2 * w)

    cl_plant = (1 - w) * cl_vsp + w * CLpost      # for the check/plot
    return dict(ag_deg=GRID_DEG, ag_rad=ag, w=w, CLpost=CLpost, CDpost=CDpost,
                CMpost=CMpost, cl_vsp=cl_vsp, cd_vsp=cd_vsp, cm_vsp=cm_vsp,
                cl_plant=cl_plant, CDmax=CDmax, cl_s=cl_s, a_stall=A_STALL_DEG)


def check(d):
    ag, clp = d["ag_deg"], d["cl_plant"]
    im = int(np.argmax(clp))
    i90 = int(np.argmin(np.abs(ag - 90)))
    i45 = int(np.argmin(np.abs(ag - 45)))
    pre = ag < (d["a_stall"] - BLEND_HALF_DEG)
    max_pre_w = float(np.max(d["w"][pre])) if pre.any() else 0.0
    print(f"Re_anchor={RE_ANCHOR:.0f}  AR={AR:.2f}  A_STALL={d['a_stall']}deg  "
          f"blend=+/-{BLEND_HALF_DEG}deg  severity={SEVERITY}")
    print(f"2D 4414 CLmax={polar_clmax_2d():.3f} (sets the realistic wing CLmax level)")
    print(f"WING stall: CLmax={clp[im]:.3f} at geometric alpha={ag[im]:.1f}deg "
          f"(VSPAERO would give {d['cl_vsp'][im]:.2f} unstalled)")
    print(f"deep stall: CL(45deg)={clp[i45]:.3f}  CL(90deg)={clp[i90]:.3f}  "
          f"CD_max(Viterna)={d['CDmax']:.2f}  CD(90deg)={d['CDpost'][i90]:.2f}")
    print(f"blend w below onset (must be 0): max={max_pre_w:.2e}")
    print(f"grid: {len(ag)} pts over [{ag[0]:.0f}, {ag[-1]:.0f}]deg")
    ok = (max_pre_w == 0.0 and abs(clp[i90]) < 0.05 and clp[im] > 1.0
          and clp[i45] < clp[im])
    print(f"sanity PASS: {ok}")
    return ok


def fmt(xs):
    return ",\n    ".join(
        ", ".join(f"{x:.8f}" for x in xs[i:i + 6]) for i in range(0, len(xs), 6))


def emit(d):
    body = f"""#pragma once
// AUTO-GENERATED by scripts/precompute_stall_table.py -- DO NOT EDIT BY HAND.
//
// Viscous-stall plant aero for the longitudinal model, as a BLEND that hands the
// lift off from the (inviscid, stall-free) VSPAERO deck to an absolute Viterna
// flat-plate post-stall curve:
//
//   C_plant(a) = (1 - w(a)) * C_vspaero(a)  +  w(a) * C_post(a)
//
// w(a) is 0 below the wing stall onset (attached flow == the VSPAERO deck) and 1
// above it -- so VSPAERO's meaningless post-stall values are weighted to zero and
// the lift comes entirely from C_post. C_post is the Viterna & Corrigan flat-plate
// extrapolation anchored at the wing stall point (CLmax = CL_vspaero(A_STALL) =
// {d['cl_s']:.3f} at {d['a_stall']:.1f} deg geometric), declining to CL=0 at 90 deg
// with CD -> CD_max = {d['CDmax']:.2f}. The NACA 4414 2D viscous polar (NeuralFoil,
// Re = {RE_ANCHOR:.0f}) sets the realistic CLmax level; Viterna gives the deep-stall tail.
//
// NOTE: the stall point is a TUNABLE engineering value, NOT ground truth -- see
// documentation/stall_model_spec.md. alpha grid ascending in RADIANS over the
// geometric range [{d['ag_deg'][0]:.0f}, {d['ag_deg'][-1]:.0f}] deg; runtime
// stallLookup() does linear-in-alpha interpolation, holding w=0 below the grid
// (attached) and the last (deep-stall) values above it.
namespace autoland {{

inline constexpr int kStallN = {len(d['ag_rad'])};

inline constexpr double kStallAlpha[kStallN] = {{
    {fmt(d['ag_rad'])}
}};

// Blend weight w(alpha) in [0,1]: 0 = pure VSPAERO (attached), 1 = pure post-stall.
inline constexpr double kStallW[kStallN] = {{
    {fmt(d['w'])}
}};

// Absolute post-stall coefficients (Viterna flat-plate), used where w > 0.
inline constexpr double kStallCLpost[kStallN] = {{
    {fmt(d['CLpost'])}
}};

inline constexpr double kStallCDpost[kStallN] = {{
    {fmt(d['CDpost'])}
}};

inline constexpr double kStallCMpost[kStallN] = {{
    {fmt(d['CMpost'])}
}};

}}  // namespace autoland
"""
    with open(HEADER, "w") as f:
        f.write(body)
    print(f"wrote {os.path.relpath(HEADER, ROOT)}")


if __name__ == "__main__":
    d = build_table()
    ok = check(d)
    if not ok:
        print("ERROR: stall-table sanity check failed", file=sys.stderr)
        sys.exit(1)
    if "--check" not in sys.argv:
        emit(d)
