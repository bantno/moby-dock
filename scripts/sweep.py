#!/usr/bin/env python3
"""Parameter sweep + scorer for the longitudinal CBF-QP water-landing sim.

Runs the sim across the CARTESIAN PRODUCT of a set of parameter grids, scores
each run on the objectives that matter for a good landing, and writes a ranked
summary so you can pick (or further tune toward) the best configuration.

  # sweep V_land x descent gains x energy-barrier slack, on the IDEAL sensor:
  python scripts/sweep.py \
      --set cbf.h_meas_stddev=0 --set cbf.h_lpf_tau=0 \
      --grid cbf.V_land=10;13 \
      --grid "cbf.c_descent=[1,4,38];[2,2,2]" \
      --grid cbf.w_slack_energy=1e2;1e3

`--grid DOTTED.KEY=v1;v2;...` (repeatable) is one swept axis; values are
';'-separated so list values like [1,4,38] survive. `--set DOTTED.KEY=VALUE`
(repeatable) is a FIXED override applied to every run (same syntax as harness.py).
Each value is parsed as YAML, so numbers, bools and lists all work.

Objectives recorded per run (all "lower is better" after the sign in []):
  td_sink   [+] touchdown sink rate [m/s]              -- the soft-landing metric
  td_V      [+] touchdown airspeed [m/s]               -- the anti-bounce metric
  worst_b   [-] min over the flight of every barrier   -- < 0 is a SAFETY VIOLATION
  rec_frac  [+] fraction of steps the QP fell back to best-effort
  chat_de   [+] RMS step-to-step elevator change [rad] -- control chatter
  chat_Td   [+] RMS step-to-step thrust-jerk change    -- control chatter

A composite SCORE (lower = better) min-max-normalizes the soft objectives across
the sweep and adds hard penalties for safety violations / QP recoveries; tune the
weights with --w-* or just re-sort the table on any raw column with --rank.
Artifacts land in runs/sweep/<name>/ (per-run CSV + resolved YAML, summary.csv,
summary.md); the best run is also plotted (harness single-run figure).
"""
import argparse
import itertools
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness  # reuse scenario resolution / trace loading / annotations

REPO = harness.REPO
DATA = harness.DATA


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=os.path.join(DATA, "lon_scenario.yaml"))
    ap.add_argument("--stab", default=os.path.join(DATA, "AHAB_combined.stab"))
    ap.add_argument("--aircraft", default=os.path.join(DATA, "aircraft.yaml"))
    ap.add_argument("--binary",
                    default=os.path.join(REPO, "build", "Release", "lon_autoland_sim.exe"))
    ap.add_argument("--name", default="sweep")
    ap.add_argument("--grid", action="append", default=[], metavar="K=v1;v2;...",
                    help="cartesian-product axis (repeatable); ';'-separated values")
    # Random-sampling mode (for high-dim "all CBFs + nominal" sweeps):
    ap.add_argument("--sample", type=int, default=0,
                    help="draw N random configs from --range/--choice instead of --grid")
    ap.add_argument("--range", action="append", default=[], dest="ranges",
                    metavar="K=lo,hi[,log]", help="continuous scalar axis, uniform "
                    "(or log-uniform with ',log'); only with --sample (repeatable)")
    ap.add_argument("--choice", action="append", default=[], dest="choices",
                    metavar="K=v1;v2;...", help="discrete axis (e.g. class-K vectors); "
                    "only with --sample (repeatable)")
    ap.add_argument("--seed", type=int, default=0, help="RNG seed for --sample")
    ap.add_argument("--set", action="append", default=[], dest="fixed",
                    metavar="K=V", help="fixed override applied to every run (repeatable)")
    ap.add_argument("--top", type=int, default=8, help="rows to print / plot")
    ap.add_argument("--rank", default="score", help="column to sort by (asc)")
    # composite-score weights (soft objectives are min-max normalized first).
    ap.add_argument("--w-sink", type=float, default=1.0)
    ap.add_argument("--w-vel", type=float, default=1.0)
    ap.add_argument("--w-chat", type=float, default=1.0)
    ap.add_argument("--w-viol", type=float, default=10.0, help="penalty per m/s of barrier violation")
    ap.add_argument("--w-rec", type=float, default=5.0, help="penalty per unit recovery fraction")
    return ap.parse_args()


def parse_axis(spec):
    """'cbf.V_land=10;13' -> ('cbf.V_land', [10.0, 13.0]) (values parsed as YAML)."""
    if "=" not in spec:
        sys.exit(f"bad --grid '{spec}' (expected DOTTED.KEY=v1;v2;...)")
    key, raw = spec.split("=", 1)
    vals = [harness.parse_value(v.strip()) for v in raw.split(";") if v.strip()]
    if not vals:
        sys.exit(f"axis '{spec}' has no values")
    return key.strip(), vals


def parse_range(spec):
    """'nominal.Kq=0.1,0.5' or 'cbf.w_slack_energy=1e1,1e4,log' -> (key, lo, hi, log)."""
    if "=" not in spec:
        sys.exit(f"bad --range '{spec}' (expected K=lo,hi[,log])")
    key, raw = spec.split("=", 1)
    parts = [p.strip() for p in raw.split(",")]
    if len(parts) < 2:
        sys.exit(f"bad --range '{spec}' (need lo,hi)")
    log = len(parts) > 2 and parts[2].lower() == "log"
    return key.strip(), float(parts[0]), float(parts[1]), log


def gen_runs(a):
    """Return (param_keys, [param_dict, ...]) for either sampling or grid mode."""
    if a.sample > 0:
        ranges = [parse_range(s) for s in a.ranges]
        choices = [parse_axis(s) for s in a.choices]
        keys = [k for k, *_ in ranges] + [k for k, _ in choices]
        if not keys:
            sys.exit("[sweep] --sample needs at least one --range or --choice")
        rng = np.random.default_rng(a.seed)
        runs = []
        for _ in range(a.sample):
            p = {}
            for k, lo, hi, log in ranges:
                u = rng.uniform(np.log(lo), np.log(hi)) if log else rng.uniform(lo, hi)
                p[k] = round(float(np.exp(u) if log else u), 6)
            for k, vals in choices:
                p[k] = vals[int(rng.integers(len(vals)))]
            runs.append(p)
        return keys, runs
    axes = [parse_axis(g) for g in a.grid]
    if not axes:
        sys.exit("[sweep] need --grid (cartesian) or --sample with --range/--choice")
    keys = [k for k, _ in axes]
    runs = [dict(zip(keys, combo))
            for combo in itertools.product(*[v for _, v in axes])]
    return keys, runs


def run_sim(args, scenario, csv_path, scn_path):
    """Write the resolved scenario, invoke the exe directly, return True on a CSV."""
    import yaml
    with open(scn_path, "w") as f:
        yaml.safe_dump(scenario, f, sort_keys=False)
    proc = subprocess.run([args.binary, args.stab, args.aircraft, scn_path, csv_path],
                          capture_output=True, text=True)
    return os.path.exists(csv_path), proc.stdout + proc.stderr


def chatter(x):
    """RMS of the step-to-step change -- a simple, robust control-smoothness proxy."""
    if len(x) < 2:
        return float("nan")
    dx = np.diff(np.asarray(x, float))
    return float(np.sqrt(np.nanmean(dx * dx)))


def score_run(d, Tmax):
    """Compute the per-run objective metrics from a loaded trace `d`.

    Barriers are split: HARD ones (descent sink-rate + thrust limits) carry the
    real forward-invariance guarantee, so worst_hard < 0 is a genuine safety
    violation; the SOFT ones (airspeed, energy, impact) are slack-penalized comfort
    constraints, so worst_soft < 0 just means the QP traded them off (reported, but
    only lightly penalized).
    """
    i_td = harness.touchdown_index(d)        # first h<=0 (the touchdown frame)
    phys = slice(0, max(i_td, 1))            # physical descent frames (h>0)
    hard = [d["b_descent"], d["T"], Tmax - d["T"]]
    soft = [d["b_airspeed"], d["b_airspeed_upper"], d["b_energy"], d["b_impact"]]
    rec = d["recovered"][phys]
    return dict(
        td_sink=float(d["sink"][i_td]),
        td_V=float(d["V"][i_td]),
        worst_hard=min(float(np.nanmin(b[phys])) for b in hard),
        worst_soft=min(float(np.nanmin(b[phys])) for b in soft),
        rec_frac=float(np.nanmean(rec)) if len(rec) else 0.0,
        chat_de=chatter(d["de"][phys]),
        chat_Td=chatter(d["Tddot"][phys]),
        reached=bool(d["h"][i_td] <= 0.0 and i_td > 0),
    )


def normalize(rows, key):
    vals = np.array([r[key] for r in rows], float)
    lo, hi = np.nanmin(vals), np.nanmax(vals)
    rng = hi - lo
    return {id(r): (0.0 if rng == 0 else (r[key] - lo) / rng) for r in rows}


def composite(rows, a):
    """Lower = better. Min-max-normalize soft objectives; add hard penalties."""
    n_sink = normalize(rows, "td_sink")
    n_vel = normalize(rows, "td_V")
    n_cde = normalize(rows, "chat_de")
    n_cTd = normalize(rows, "chat_Td")
    for r in rows:
        viol_h = max(0.0, -r["worst_hard"])   # genuine safety violation
        viol_s = max(0.0, -r["worst_soft"])   # soft-constraint slack (light penalty)
        r["score"] = (a.w_sink * n_sink[id(r)] + a.w_vel * n_vel[id(r)]
                      + a.w_chat * 0.5 * (n_cde[id(r)] + n_cTd[id(r)])
                      + a.w_viol * viol_h + 0.05 * a.w_viol * viol_s
                      + a.w_rec * r["rec_frac"]
                      + (0.0 if r["reached"] else 100.0))  # no-touchdown = unusable


def fmt_table(rows, keys, cols, top):
    head = keys + cols
    widths = {h: max(len(h), 9) for h in head}
    lines = ["  ".join(h.rjust(widths[h]) for h in head)]
    for r in rows[:top]:
        cells = []
        for k in keys:
            cells.append(str(r["params"][k]).rjust(widths[k]))
        for c in cols:
            v = r[c]
            s = f"{v:.4g}" if isinstance(v, float) else str(v)
            cells.append(s.rjust(widths[c]))
        lines.append("  ".join(cells))
    return "\n".join(lines)


def main():
    a = parse_args()
    a.binary = os.path.abspath(a.binary)
    if not os.path.exists(a.binary):
        sys.exit(f"[sweep] no binary at {a.binary} -- build lon_autoland_sim first")
    keys, runs = gen_runs(a)
    outdir = os.path.join(REPO, "runs", "sweep", a.name)
    os.makedirs(outdir, exist_ok=True)
    mode = f"random sample N={a.sample}" if a.sample > 0 else "cartesian grid"
    print(f"[sweep] {len(runs)} runs ({mode}) over: {', '.join(keys)}")

    import yaml as _yaml

    def ov_for(params):
        # fixed overrides + this run's key=value (each re-serialized as a YAML scalar)
        return list(a.fixed) + [f"{k}={_yaml.safe_dump(v).strip()}"
                                for k, v in params.items()]

    rows = []
    for i, params in enumerate(runs):
        scenario = harness.resolve_scenario(a.base, ov_for(params))
        ann = harness.ann(scenario)
        csv_path = os.path.join(outdir, f"run_{i:03d}.csv")
        scn_path = os.path.join(outdir, f"run_{i:03d}.resolved.yaml")
        ok, log = run_sim(a, scenario, csv_path, scn_path)
        if not ok:
            print(f"  run {i:03d}: FAILED (no CSV) -- {log.strip()[:120]}")
            continue
        d = harness.load_trace(csv_path)
        m = score_run(d, ann["Tmax"])
        m["params"] = dict(params)
        m["idx"] = i
        rows.append(m)
        flag = "" if m["reached"] else "  [NO TOUCHDOWN]"
        print(f"  run {i:03d}: sink={m['td_sink']:.3f} V={m['td_V']:.2f} "
              f"hard_b={m['worst_hard']:+.3f} soft_b={m['worst_soft']:+.1f} "
              f"rec={m['rec_frac']:.2%} chat_de={m['chat_de']:.4f}{flag}")

    if not rows:
        sys.exit("[sweep] no successful runs")
    composite(rows, a)
    cols = ["score", "td_sink", "td_V", "worst_hard", "worst_soft", "rec_frac", "chat_de", "chat_Td"]
    if a.rank not in cols and a.rank not in keys:
        sys.exit(f"[sweep] --rank '{a.rank}' not a known column {cols}")
    rows.sort(key=lambda r: (r[a.rank] if a.rank in r else r["params"][a.rank]))

    table = fmt_table(rows, keys, cols, a.top)
    print(f"\n[sweep] top {min(a.top, len(rows))} by {a.rank} (lower=better):\n{table}")

    # summary.csv (all runs) + summary.md (top-N).
    import csv as _csv
    with open(os.path.join(outdir, "summary.csv"), "w", newline="") as f:
        w = _csv.writer(f)  # proper quoting: list values like "[2, 2, 2]" hold commas
        w.writerow(["idx"] + keys + cols + ["reached"])
        for r in rows:
            w.writerow([r["idx"]] + [r["params"][k] for k in keys]
                       + [f"{r[c]:.6g}" for c in cols] + [r["reached"]])
    with open(os.path.join(outdir, "summary.md"), "w") as f:
        f.write(f"# Sweep `{a.name}` -- top {a.top} by `{a.rank}`\n\n")
        f.write("| " + " | ".join(keys + cols) + " |\n")
        f.write("|" + "---|" * len(keys + cols) + "\n")
        for r in rows[:a.top]:
            cells = [str(r["params"][k]) for k in keys] \
                + [(f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c])) for c in cols]
            f.write("| " + " | ".join(cells) + " |\n")
    print(f"[sweep] wrote {outdir}/summary.{{csv,md}}")

    # Plot the best run with the harness single-run figure.
    best = rows[0]
    d = harness.load_trace(os.path.join(outdir, f"run_{best['idx']:03d}.csv"))
    ann = harness.ann(harness.resolve_scenario(a.base, ov_for(best["params"])))
    harness.plot_single(d, ann, os.path.join(outdir, "best.png"),
                        f"BEST ({a.rank}) -- " + ", ".join(f"{k}={best['params'][k]}" for k in keys))


if __name__ == "__main__":
    main()
