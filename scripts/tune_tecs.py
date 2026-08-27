#!/usr/bin/env python3
"""Tune the PX4 TECS nominal on the Beaver plant -- with a hold-out set, so the result is a
setting that generalises rather than one fitted to a single run.

Protocol
--------
* TRAINING cases (used to pick gains): longitudinal disturbances at the design condition
  (40 m/s clean, 1800 RPM, h0 = 60 m): calm regulation, hot entry (+3 m/s, +2 deg), cold entry
  (-3 m/s, -2 deg), tailwind shear (+5 m/s over 120 m at t = 5 s).
* VALIDATION cases (held out -- never used to pick gains): the POH flaps-35 float
  configuration at 33.5 m/s / 2000 RPM (calm and hot; a different trim, throttle map and pitch
  offset), a downdraft gust, a 5 m/s crosswind (lateral coupling through the load factor), a
  120 m approach (slow drift), and 45 m/s / -3 deg (different airspeed and glideslope).
* Gains searched: ptch_damp, i_gain_pit, thr_damping, thr_integ on coarse grids inside PX4's
  parameter ranges (FW_T_PTCH_DAMP 0..2, FW_T_I_GAIN_PIT 0..2, FW_T_THR_DAMPING 0..1,
  FW_T_THR_INTEG 0..1). Everything else stays at the PX4 default; the vehicle anchors are the
  plant-derived ones. A one-at-a-time sensitivity pass over tas_tc / vert_acc / spdweight /
  ste_r_tc / seb_r_ff follows at the chosen point.
* Score per run: J = IAE(sink - sink_ref) + 0.5 IAE(V - V_ref) [m] + 0.02 TV(delta_e) [deg]
  (init sample excluded). For the grid each case's J is normalised by the PX4-DEFAULT TECS on
  that case (Jrel = 1 is the starting point; the cascaded PID is shown for reference only --
  it has an exact approach-trim feedforward, so its calm-case J is ~0 and would swamp a
  ratio), and the training score is the mean Jrel over the training cases. Hard checks are
  reported alongside: max sink, min airspeed, max alpha, max per-step elevator increment
  (limit-cycle detector), touchdown sink / V / gamma, settling time.
* Selection: the ROBUST pick minimises the worst score over the combo and its grid neighbours
  (+-1 index in every gain) -- a flat good region, not a lucky spike. The raw best is reported
  next to it so the two can be compared.

Usage:
  tune_tecs.py grid  [--jobs N] [--out runs/tecs_tune]      # grid on the training cases
  tune_tecs.py eval  [--gains "ptch_damp=1.0,i_gain_pit=0.2,..."] [more --gains ...]
                                                            # train + validation tables for the
                                                            # PX4 defaults, the tuned set, and any
                                                            # --gains sets
  tune_tecs.py sens  --gains "..."                          # one-at-a-time sensitivity
"""
import argparse
import concurrent.futures as cf
import itertools
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "build" / "sixdof_autoland_sim"
DEG = np.pi / 180.0

# ----------------------------------------------------------------------------- cases
BASE = {
    "plant": "beaver", "V_app": 40.0, "gamma_app_deg": -3.5,
    "beaver": {"n_rpm": 1800.0, "pz_idle": 5.0, "pz_max": 26.0, "flap_deg": 0.0},
    "dt": 0.01, "t_max": 120.0, "initial": {"h0": 60.0}, "hull": {"n_surfaces": 2},
}

def case(name, **over):
    d = yaml.safe_load(yaml.safe_dump(BASE))  # deep copy
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(d.get(k), dict):
            d[k].update(v)
        else:
            d[k] = v
    d["_name"] = name
    return d

TRAIN = [
    case("calm"),
    case("hot_entry", initial={"h0": 60.0, "dV": 3.0, "dtheta_deg": 2.0}),
    case("cold_entry", initial={"h0": 60.0, "dV": -3.0, "dtheta_deg": -2.0}),
    case("tailwind_shear", wind={"enabled": True, "t_start": 5.0, "u_amp": 5.0, "u_len": 120.0}),
]
VALID = [
    case("poh_calm", V_app=33.5, beaver={"n_rpm": 2000.0, "flap_deg": 35.0}),
    case("poh_hot", V_app=33.5, beaver={"n_rpm": 2000.0, "flap_deg": 35.0},
         initial={"h0": 60.0, "dV": 3.0, "dtheta_deg": 2.0}),
    case("downdraft", wind={"enabled": True, "t_start": 5.0, "w_amp": -3.0, "w_len": 60.0}),
    case("crosswind", wind={"enabled": True, "t_start": 5.0, "v_amp": 5.0, "v_len": 0.0}),
    case("long_approach", initial={"h0": 120.0}),
    case("fast_shallow", V_app=45.0, gamma_app_deg=-3.0),
    case("steep_idle", gamma_app_deg=-6.0),  # near the idle sink (4.3 m/s): throttle on its rail
]

# ----------------------------------------------------------------------------- running
def run_case(c, nominal, gains, workdir, tag):
    """Run one scenario; return the metrics dict (None if the sim failed)."""
    d = {k: v for k, v in c.items() if not k.startswith("_")}
    d["nominal"] = {"type": nominal}
    if nominal == "tecs" and gains:
        d["tecs"] = dict(gains)
    workdir.mkdir(parents=True, exist_ok=True)
    scen = workdir / f"{tag}.yaml"
    csv = workdir / f"{tag}.csv"
    scen.write_text(yaml.safe_dump(d))
    p = subprocess.run([str(BIN), "", "", str(scen), str(csv)], capture_output=True, text=True)
    if not csv.exists():
        return None
    m = metrics(csv, c)
    m["reached"] = "TOUCHDOWN" in p.stdout
    m["case"] = c["_name"]
    m["nominal"] = nominal
    return m

def metrics(csv, c):
    d = np.genfromtxt(csv, delimiter=",", names=True)
    d = d[1:]  # drop the t = 0 sample (TECS initialize() kick, upstream behaviour)
    t = d["t"]
    V_ref = c["V_app"]
    sink_ref = -V_ref * np.sin(c["gamma_app_deg"] * DEG)
    e_h = d["sink"] - sink_ref
    e_V = d["V_air"] - V_ref
    dt = np.diff(t, prepend=t[0] - (t[1] - t[0]))
    iae_h = float(np.sum(np.abs(e_h) * dt))
    iae_V = float(np.sum(np.abs(e_V) * dt))
    tv_de = float(np.sum(np.abs(np.diff(d["de"]))) / DEG)
    tv_dT = float(np.sum(np.abs(np.diff(d["dT"]))))
    de_step = float(np.max(np.abs(np.diff(d["de"]))) / DEG)  # deg per 10 ms step
    bad = (np.abs(e_h) > 0.15) | (np.abs(e_V) > 0.3)
    t_settle = float(t[bad].max()) if bad.any() else 0.0
    return dict(
        J=iae_h + 0.5 * iae_V + 0.02 * tv_de,
        iae_h=iae_h, iae_V=iae_V, tv_de=tv_de, tv_dT=tv_dT, de_step=de_step, t_settle=t_settle,
        sink_max=float(d["sink"].max()), sink_min=float(d["sink"].min()),
        V_min=float(d["V_air"].min()), V_max=float(d["V_air"].max()),
        alpha_max=float(d["alpha_deg"].max()),
        td_sink=float(d["sink"][-1]), td_V=float(d["V_air"][-1]),
        td_gamma=float(d["gamma_deg"][-1]), t_td=float(t[-1]),
        e_td_sink=float(d["sink"][-1] - sink_ref), e_td_V=float(d["V_air"][-1] - V_ref),
    )

def run_many(jobs, tasks):
    out = [None] * len(tasks)
    with cf.ThreadPoolExecutor(max_workers=jobs) as ex:
        futs = {ex.submit(run_case, *args): i for i, args in enumerate(tasks)}
        for f in cf.as_completed(futs):
            out[futs[f]] = f.result()
    return out

def parse_gains(s):
    g = {}
    for kv in s.split(","):
        k, v = kv.split("=")
        g[k.strip()] = float(v)
    return g

def fmt_gains(g):
    return ", ".join(f"{k}={v:g}" for k, v in g.items())

# ----------------------------------------------------------------------------- tables
COLS = ["J", "Jrel", "t_settle", "sink_max", "sink_min", "V_min", "V_max", "alpha_max",
        "td_sink", "td_V", "td_gamma", "tv_de", "de_step", "tv_dT"]

def table(rows, title):
    print(f"\n### {title}\n")
    hdr = ["case", "nominal"] + COLS
    print("| " + " | ".join(hdr) + " |")
    print("|" + "---|" * len(hdr))
    for r in rows:
        vals = [r["case"], r["nominal"]] + [f"{r.get(k, float('nan')):.3g}" for k in COLS]
        print("| " + " | ".join(vals) + " |")

def with_rel(rows, ref):
    for r in rows:
        r["Jrel"] = r["J"] / ref[r["case"]]["J"]
    return rows

# ----------------------------------------------------------------------------- modes
PX4DEF = dict(ptch_damp=0.1, i_gain_pit=0.1, thr_damping=0.05, thr_integ=0.02)
# Result of this protocol (2026-08-27) -- the Beaver defaults in sixdof_sim.cpp: three gains
# changed, thr_damping left at the PX4 default (insensitive: 0.2..1.0 gains ~1% on validation).
TUNED = dict(ptch_damp=1.0, i_gain_pit=0.4, thr_damping=0.05, thr_integ=0.3)
# Pass 1 (0.1..2 / 0.1..0.8 / 0.05..1 / 0.05..0.5) put the optimum on the i_gain_pit and
# thr_integ edges; pass 2 extends both toward the PX4 maxima (2.0 / 1.0).
GRID = {
    "ptch_damp": [0.5, 1.0, 1.5, 2.0],
    "i_gain_pit": [0.4, 0.8, 1.2, 1.6, 2.0],
    "thr_damping": [0.2, 0.5, 1.0],
    "thr_integ": [0.3, 0.5, 0.75, 1.0],
}

def cascade_reference(cases, jobs, workdir):
    res = run_many(jobs, [(c, "cascade", None, workdir, f"cascade_{c['_name']}") for c in cases])
    return {r["case"]: r for r in res}

def px4def_reference(cases, jobs, workdir):
    res = run_many(jobs, [(c, "tecs", PX4DEF, workdir, f"px4def_{c['_name']}") for c in cases])
    return {r["case"]: r for r in res}

def mode_grid(args):
    workdir = Path(args.out)
    ref = px4def_reference(TRAIN, args.jobs, workdir / "ref")
    keys = list(GRID)
    combos = list(itertools.product(*[GRID[k] for k in keys]))
    tasks = []
    for ci, combo in enumerate(combos):
        g = dict(zip(keys, combo))
        for c in TRAIN:
            tasks.append((c, "tecs", g, workdir / "grid", f"c{ci}_{c['_name']}"))
    print(f"grid: {len(combos)} combos x {len(TRAIN)} training cases = {len(tasks)} runs", flush=True)
    res = run_many(args.jobs, tasks)
    # score per combo = mean relative J over the training cases; inf if any run failed
    score = np.full(len(combos), np.inf)
    detail = {}
    for ci in range(len(combos)):
        rs = res[ci * len(TRAIN):(ci + 1) * len(TRAIN)]
        if any(r is None or not r["reached"] for r in rs):
            continue
        score[ci] = float(np.mean([r["J"] / ref[r["case"]]["J"] for r in rs]))
        detail[ci] = rs
    # robust score: worst over the combo and its +-1 grid neighbours in every gain
    shape = tuple(len(GRID[k]) for k in keys)
    S = score.reshape(shape)
    robust = np.full(shape, np.inf)
    for idx in itertools.product(*[range(n) for n in shape]):
        vals = [S[idx]]
        for dim in range(len(shape)):
            for dlt in (-1, 1):
                j = list(idx); j[dim] += dlt
                if 0 <= j[dim] < shape[dim]:
                    vals.append(S[tuple(j)])
        robust[idx] = max(vals)
    best = int(np.argmin(score))
    rob = int(np.argmin(robust.reshape(-1)))
    # write the whole grid
    with open(workdir / "grid_results.csv", "w") as f:
        f.write(",".join(keys + ["score", "robust"] + [f"{c['_name']}_J" for c in TRAIN]) + "\n")
        for ci, combo in enumerate(combos):
            js = [detail[ci][k]["J"] for k in range(len(TRAIN))] if ci in detail else [np.nan] * len(TRAIN)
            f.write(",".join(f"{v:g}" for v in list(combo) + [score[ci], robust.reshape(-1)[ci]] + js) + "\n")
    print(f"\nPX4 defaults (score 1 by construction): {fmt_gains(PX4DEF)}")
    print(f"raw best   (score {score[best]:.3f}, robust {robust.reshape(-1)[best]:.3f}): {fmt_gains(dict(zip(keys, combos[best])))}")
    print(f"robust pick (score {score[rob]:.3f}, robust {robust.reshape(-1)[rob]:.3f}): {fmt_gains(dict(zip(keys, combos[rob])))}")
    order = np.argsort(score)
    print("\ntop 15 by training score:")
    for ci in order[:15]:
        print(f"  score {score[ci]:.3f}  robust {robust.reshape(-1)[ci]:.3f}  {fmt_gains(dict(zip(keys, combos[ci])))}")
    print("\ntop 10 by robust score:")
    for ci in np.argsort(robust.reshape(-1))[:10]:
        print(f"  score {score[ci]:.3f}  robust {robust.reshape(-1)[ci]:.3f}  {fmt_gains(dict(zip(keys, combos[ci])))}")
    # marginals: mean score over the rest of the grid for each value of each gain
    print("\nmarginal mean score per gain value (finite combos only):")
    for k in keys:
        vals = []
        for i, v in enumerate(GRID[k]):
            sl = [slice(None)] * len(shape); sl[keys.index(k)] = i
            sub = S[tuple(sl)]
            fin = sub[np.isfinite(sub)]
            vals.append(f"{v:g}: {fin.mean():.3f} (min {fin.min():.3f})" if fin.size else f"{v:g}: -")
        print(f"  {k:12s} " + " | ".join(vals))

def evaluate(gainsets, cases, jobs, workdir, title):
    ref = px4def_reference(cases, jobs, workdir / "ref")
    casc = cascade_reference(cases, jobs, workdir / "ref")
    rows = list(casc.values())
    for label, g in gainsets:
        res = run_many(jobs, [(c, "tecs", g, workdir / label, f"{label}_{c['_name']}") for c in cases])
        for r in res:
            r["nominal"] = f"tecs:{label}"
        rows += res
    with_rel(rows, ref)
    rows.sort(key=lambda r: ([c["_name"] for c in cases].index(r["case"]), r["nominal"]))
    table(rows, title)
    print("\nper setting (Jrel = J / J of the PX4-default TECS on the same case):")
    for label in ["cascade"] + [f"tecs:{l}" for l, _ in gainsets]:
        sel = [r for r in rows if r["nominal"] == label]
        print(f"  {label:22s} mean Jrel {np.mean([r['Jrel'] for r in sel]):.3f}  worst Jrel {max(r['Jrel'] for r in sel):.3f}"
              f"  sum J {sum(r['J'] for r in sel):6.1f} m  max sink {max(r['sink_max'] for r in sel):.2f}"
              f"  min V {min(r['V_min'] for r in sel):.1f}  max de_step {max(r['de_step'] for r in sel):.2f} deg")

def mode_eval(args):
    gainsets = [(f"g{i}", parse_gains(s)) for i, s in enumerate(args.gains)]
    gainsets = [("px4def", PX4DEF),
                ("integ_only", dict(ptch_damp=0.1, i_gain_pit=0.4, thr_damping=0.05, thr_integ=0.3)),
                ("tuned", TUNED)] + gainsets
    for label, g in gainsets:
        print(f"{label}: {fmt_gains(g)}")
    evaluate(gainsets, TRAIN, args.jobs, Path(args.out) / "eval_train", "TRAINING cases")
    evaluate(gainsets, VALID, args.jobs, Path(args.out) / "eval_valid", "VALIDATION cases (held out)")

def mode_sens(args):
    g0 = parse_gains(args.gains[0])
    sens = {"ptch_damp": [0.75, 1.25], "i_gain_pit": [0.2, 0.3, 0.6],
            "tas_tc": [3.0, 5.0, 8.0], "vert_acc": [3.0, 7.0, 10.0], "spdweight": [0.5, 1.0, 1.5],
            "ste_r_tc": [0.1, 0.4, 1.0], "seb_r_ff": [0.5, 1.0, 1.5]}
    workdir = Path(args.out) / "sens"
    ref = px4def_reference(TRAIN, args.jobs, workdir / "ref")
    base = run_many(args.jobs, [(c, "tecs", g0, workdir / "base", f"base_{c['_name']}") for c in TRAIN])
    print(f"base {fmt_gains(g0)}: mean Jrel {np.mean([r['J']/ref[r['case']]['J'] for r in base]):.3f}")
    for k, vals in sens.items():
        for v in vals:
            g = dict(g0); g[k] = v
            res = run_many(args.jobs, [(c, "tecs", g, workdir / f"{k}_{v}", f"{k}{v}_{c['_name']}") for c in TRAIN])
            ok = all(r is not None and r["reached"] for r in res)
            jr = np.mean([r["J"] / ref[r["case"]]["J"] for r in res]) if ok else float("inf")
            worst = max(r["J"] / ref[r["case"]]["J"] for r in res) if ok else float("inf")
            print(f"  {k}={v:g}: mean Jrel {jr:.3f}  worst {worst:.3f}  sink_max {max(r['sink_max'] for r in res):.2f}  V_min {min(r['V_min'] for r in res):.2f}")

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["grid", "eval", "sens"])
    ap.add_argument("--gains", action="append", default=[], help="k=v,k=v (tecs: keys)")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--out", default=str(ROOT / "runs" / "tecs_tune"))
    args = ap.parse_args()
    if not BIN.exists():
        sys.exit(f"build the sim first: {BIN}")
    {"grid": mode_grid, "eval": mode_eval, "sens": mode_sens}[args.mode](args)

if __name__ == "__main__":
    main()
