#!/usr/bin/env python3
"""Run-and-plot harness for the longitudinal CBF-QP water-landing sim.

One command builds (if stale), runs `lon_autoland_sim` on a scenario, and plots
the trace. Every plot annotation (V_td_max, alpha_stall, theta_min, Tmax, ...) is
read from the *resolved* scenario, so the figure can never drift out of sync
with the run the way the hand-passed `--vmin 15` plot-script args could.

  # single run, plot to runs/<name>.png
  scripts/harness.py run
  scripts/harness.py run --name dive --set gamma_app_deg=-60 --set cbf.V_td_max=12
  scripts/harness.py run --set cbf.enabled=false -o figures/nominal_only.png

  # overlay several cases (subsumes plot_lon_compare / plot_gain_compare)
  scripts/harness.py compare \
      --case "CBF on:cbf.enabled=true" \
      --case "CBF off:cbf.enabled=false"
  scripts/harness.py compare --name gains \
      --case "cd2:cbf.c_descent=[2,2,2]" \
      --case "cd10:cbf.c_descent=[10,10,10]"

`--set DOTTED.KEY=VALUE` (repeatable) overrides any scenario field; the value is
parsed as YAML, so numbers, bools and lists (`[2,2,2]`) all work. Each `--case`
is `LABEL:k=v;k=v;...` -- the same overrides applied on top of the base scenario,
';'-separated so list values like `[2,2,2]` survive (a ',' would be eaten).
"""
import argparse
import os
import subprocess
import sys

import numpy as np
import yaml

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, "data")
DEFAULTS = dict(
    scenario=os.path.join(DATA, "lon_scenario.yaml"),
    stab=os.path.join(DATA, "AHAB_combined.stab"),
    aircraft=os.path.join(DATA, "aircraft.yaml"),
    binary=os.path.join(REPO, "build", "lon_autoland_sim"),
    build_dir=os.path.join(REPO, "build"),
    outdir=os.path.join(REPO, "runs"),
)


# --------------------------------------------------------------------------- #
# Scenario resolution: load YAML, apply dotted-key overrides.
# --------------------------------------------------------------------------- #
def parse_value(raw):
    """Parse an override RHS as a YAML scalar (int/float/bool/list/str)."""
    return yaml.safe_load(raw)


def set_dotted(d, dotted, value):
    """Set d['a']['b'] = value for dotted key 'a.b', creating maps as needed."""
    keys = dotted.split(".")
    node = d
    for k in keys[:-1]:
        if not isinstance(node.get(k), dict):
            node[k] = {}
        node = node[k]
    node[keys[-1]] = value


def apply_overrides(scenario, overrides):
    """overrides: list of 'dotted.key=value' strings, applied in order."""
    for ov in overrides:
        if "=" not in ov:
            sys.exit(f"bad override '{ov}' (expected dotted.key=value)")
        key, raw = ov.split("=", 1)
        set_dotted(scenario, key.strip(), parse_value(raw.strip()))
    return scenario


def resolve_scenario(base_path, overrides):
    with open(base_path) as f:
        scenario = yaml.safe_load(f) or {}
    return apply_overrides(scenario, overrides)


def ann(scenario):
    """Pull the handful of plot-annotation knobs out of a resolved scenario,
    falling back to the same defaults lon_sim.cpp uses."""
    V_app = scenario.get("V_app", 18.0)
    cbf = scenario.get("cbf", {}) or {}
    return dict(
        V_app=V_app,
        V_td_max=cbf.get("V_td_max", 14.0),
        g_eff=cbf.get("g_eff", 16.0),
        alpha_stall_deg=cbf.get("alpha_stall_deg", 11.0),
        stall_margin_deg=cbf.get("stall_margin_deg", 2.0),
        theta_min_deg=cbf.get("theta_min_deg", 3.0),
        h_noseup=cbf.get("h_noseup", 3.0),
        Tmax=cbf.get("Tmax", 50.0),
        n_limit=cbf.get("n_limit", 3.0),
        cbf_on=bool(cbf.get("enabled", True)) and bool(scenario.get("cbf_enabled", True)),
    )


# --------------------------------------------------------------------------- #
# Build + run
# --------------------------------------------------------------------------- #
def newest_mtime(*dirs):
    newest = 0.0
    for d in dirs:
        for root, _, files in os.walk(d):
            for fn in files:
                if fn.endswith((".cpp", ".hpp", ".h")):
                    newest = max(newest, os.path.getmtime(os.path.join(root, fn)))
    return newest


def ensure_binary(args):
    """Build the sim if the binary is missing, stale vs sources, or --build."""
    need = args.build or not os.path.exists(args.binary)
    if not need and not args.no_build:
        src_t = newest_mtime(os.path.join(REPO, "src"), os.path.join(REPO, "include"),
                             os.path.join(REPO, "apps"))
        if src_t > os.path.getmtime(args.binary):
            print("[harness] sources newer than binary -> rebuilding")
            need = True
    if args.no_build:
        if not os.path.exists(args.binary):
            sys.exit(f"[harness] no binary at {args.binary} and --no-build given")
        return
    if not need:
        return
    target = os.path.basename(args.binary)
    print(f"[harness] building {target} ...")
    r = subprocess.run(["cmake", "--build", args.build_dir, "-j", "--target", target])
    if r.returncode != 0:
        sys.exit(f"[harness] build failed ({r.returncode})")


def run_sim(args, scenario, name):
    """Write the resolved scenario, run the sim, return (csv_path, stdout)."""
    os.makedirs(args.outdir, exist_ok=True)
    scn_path = os.path.join(args.outdir, f"{name}.resolved.yaml")
    csv_path = os.path.join(args.outdir, f"{name}.csv")
    with open(scn_path, "w") as f:
        yaml.safe_dump(scenario, f, sort_keys=False)
    print(f"[harness] run '{name}': {args.stab} + {scn_path}")
    proc = subprocess.run(
        [args.binary, args.stab, args.aircraft, scn_path, csv_path],
        capture_output=True, text=True)
    sys.stdout.write(proc.stdout)
    if proc.returncode != 0 and not os.path.exists(csv_path):
        sys.stderr.write(proc.stderr)
        sys.exit(f"[harness] sim failed for '{name}' ({proc.returncode})")
    with open(os.path.join(args.outdir, f"{name}.log"), "w") as f:
        f.write(proc.stdout + proc.stderr)
    return csv_path


def load_trace(csv_path):
    return np.genfromtxt(csv_path, delimiter=",", names=True)


def touchdown_index(d):
    idx = np.where(d["h"] <= 0)[0]
    return int(idx[0]) if len(idx) else len(d["h"]) - 1


def barrier_report(d, a):
    """Min of each barrier over the trace. Only the HARD rows (impact + thrust)
    fail the verdict; the soft rows (stall/nose-up/energy) report their min and
    are allowed to dip transiently. Nose-up is judged only over its active window
    (h < h_noseup), where the row is actually assembled."""
    h = d["h"]
    # name -> (value, active-mask or None, hard?)
    bars = {
        "b_stall":               (d["b_stall"],       None,               False),
        "b_noseup":              (d["b_noseup"],      h < a["h_noseup"],  False),
        "b_energy":              (d["b_energy"],      None,               False),
        "b_impact":              (d["b_impact"],      None,               True),
        "b_thrust_min (T)":      (d["T"],             None,               True),
        "b_thrust_max (Tmax-T)": (a["Tmax"] - d["T"], None,               True),
    }
    print("[harness] barrier minima (hard rows must be >= 0; soft rows may dip):")
    worst_ok = True
    for nm, (v, mask, hard) in bars.items():
        vv = v[mask] if mask is not None else v
        if len(vv) == 0:
            print(f"    {nm:24s} (never active)")
            continue
        nnan = int(np.isnan(vv).sum())
        mn = float(np.nanmin(vv)) if nnan < len(vv) else float("nan")
        ok = nnan == 0 and mn >= -1e-6
        if hard:
            worst_ok &= ok
            flag = "OK" if ok else ("NaN" if nnan else "VIOLATION")
        else:
            flag = "soft" if mn >= -1e-6 else "soft (dips)"
        tail = f"  ({nnan} NaN)" if nnan else ""
        print(f"    {nm:24s} min={mn:+.4f}  {flag}{tail}")
    i = touchdown_index(d)
    print(f"[harness] touchdown: t={d['t'][i]:.2f}s  sink={d['sink'][i]:.3f} m/s  "
          f"V={d['V'][i]:.2f} m/s  vs V_td_max={a['V_td_max']}")
    return worst_ok


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #
def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_single(d, a, out, title):
    plt = _mpl()
    t = d["t"]
    fig, ax = plt.subplots(3, 2, figsize=(13, 10))
    fig.suptitle(title, fontsize=14, weight="bold")

    ax[0, 0].plot(t, d["h"], color="tab:blue")
    ax[0, 0].axhline(0, color="k", lw=0.8, ls=":")
    ax[0, 0].set(ylabel="altitude h [m]", title="Descent profile")

    ax[0, 1].plot(t, d["sink"], color="tab:red", label="sink rate")
    ax[0, 1].axhline(0, color="k", ls=":", lw=0.9)
    ax[0, 1].set(ylabel="sink [m/s] (down +)",
                 title="Sink rate (impact barrier regulates touchdown)")
    ax[0, 1].legend(fontsize=8)

    g = 9.80665
    vcap = np.sqrt(a["V_td_max"] ** 2 + 2 * (a["g_eff"] - g) * np.clip(d["h"], 0, None))
    ax[1, 0].plot(t, d["V"], color="tab:green", label="V")
    ax[1, 0].plot(t, vcap, color="tab:gray", ls="--",
                  label=r"energy cap $\sqrt{V_{td}^2+2(g_{eff}-g)h}$")
    ax[1, 0].axhline(a["V_td_max"], color="k", ls=":", lw=0.9, label=f"V_td_max={a['V_td_max']}")
    ax[1, 0].set(ylabel="V [m/s]", title="Airspeed vs energy-ceiling cap")
    ax[1, 0].legend(fontsize=8)

    alpha_lim = a["alpha_stall_deg"] - a["stall_margin_deg"]
    ax[1, 1].plot(t, d["theta_deg"], label=r"$\theta$")
    ax[1, 1].plot(t, d["theta_cmd_deg"], ls="--", label=r"$\theta_{cmd}$")
    ax[1, 1].plot(t, d["gamma_deg"], label=r"$\gamma$")
    ax[1, 1].plot(t, d["alpha_deg"], label=r"$\alpha$")
    ax[1, 1].axhline(alpha_lim, color="tab:red", ls=":", lw=0.9,
                     label=f"alpha_lim={alpha_lim:.0f}")
    ax[1, 1].axhline(a["theta_min_deg"], color="tab:brown", ls=":", lw=0.9,
                     label=f"theta_min={a['theta_min_deg']:.0f}")
    ax[1, 1].set(ylabel="angle [deg]", title="Attitude / flight-path angles")
    ax[1, 1].legend(fontsize=7, ncol=2)

    ax[2, 0].plot(t, np.rad2deg(d["de"]), color="tab:purple", label=r"$\delta_e$ filtered")
    ax[2, 0].plot(t, np.rad2deg(d["de_nom"]), color="tab:purple", ls=":", alpha=0.7,
                  label=r"$\delta_e$ nominal")
    ax[2, 0].set(ylabel="elevator [deg]", xlabel="t [s]",
                 title="Elevator: nominal vs CBF-filtered")
    ax[2, 0].legend(fontsize=8, loc="upper left")
    axT = ax[2, 0].twinx()
    axT.plot(t, d["T"], color="tab:orange", alpha=0.6)
    axT.set_ylabel("thrust T [N]", color="tab:orange")

    bars = {
        "b_stall": d["b_stall"],
        "b_noseup": d["b_noseup"],
        "b_energy": d["b_energy"],
        "b_impact": d["b_impact"],
        "b_thrust_min": d["T"],
        "b_thrust_max": a["Tmax"] - d["T"],
    }
    for nm, v in bars.items():
        ax[2, 1].plot(t, v, lw=1.1, label=f"{nm} (min {np.min(v):.2f})")
    ax[2, 1].axhline(0, color="k", lw=1.0, ls="--")
    ax[2, 1].set(ylabel="barrier value", xlabel="t [s]",
                 title="All barriers (>= 0 = safe)")
    ax[2, 1].legend(fontsize=7)

    for x in ax.flat:
        x.grid(alpha=0.3)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out, dpi=110)
    print(f"[harness] wrote {out}")


def plot_compare(cases, a, out, title):
    """cases: list of (label, trace dict). Annotations a come from the base."""
    plt = _mpl()
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    fig, ax = plt.subplots(2, 3, figsize=(17, 9))
    fig.suptitle(title, fontsize=14, weight="bold")

    hmax = max(c["h"].max() for _, c in cases)

    tds = []
    for k, (lbl, d) in enumerate(cases):
        c = colors[k % len(colors)]
        i = touchdown_index(d)
        tds.append((lbl, d["sink"][i]))
        ax[0, 0].plot(d["t"], d["h"], color=c, label=lbl)
        ax[0, 1].plot(d["h"], d["sink"], color=c, label=lbl)
        ax[0, 2].plot(d["t"], d["V"], color=c, label=lbl)
        ax[1, 0].plot(d["t"], d["gamma_deg"], color=c, label=fr"{lbl} $\gamma$")
        ax[1, 0].plot(d["t"], d["theta_deg"], color=c, ls="--", alpha=0.7)
        ax[1, 1].plot(d["t"], np.rad2deg(d["de"]), color=c, label=f"{lbl} filtered")
        ax[1, 1].plot(d["t"], np.rad2deg(d["de_nom"]), color=c, ls=":", alpha=0.6)
        ax[1, 2].plot(d["t"], d["b_impact"], color=c, label=lbl)

    ax[0, 0].axhline(0, color="k", lw=0.8, ls=":")
    ax[0, 0].set(xlabel="t [s]", ylabel="h [m]", title="Descent profile")
    ax[0, 1].set(xlabel="altitude h [m]", ylabel="sink [m/s]",
                 title="Sink vs altitude")
    ax[0, 1].set_xlim(min(12, hmax), 0)
    ax[0, 2].axhline(a["V_td_max"], color="k", ls=":", lw=0.9, label=f"V_td_max={a['V_td_max']}")
    ax[0, 2].set(xlabel="t [s]", ylabel="V [m/s]", title="Airspeed")
    ax[1, 0].set(xlabel="t [s]", ylabel="angle [deg]",
                 title=r"$\gamma$ (solid) / $\theta$ (dashed)")
    ax[1, 1].set(xlabel="t [s]", ylabel="elevator [deg]",
                 title="Elevator (filtered solid, nominal dotted)")
    td_str = "  ".join(f"{l}: {s:.3f}" for l, s in tds)
    ax[1, 2].axhline(0, color="k", lw=1.0, ls="--")
    ax[1, 2].set(xlabel="t [s]", ylabel="b_impact",
                 title=f"Impact barrier  (td sink m/s -- {td_str})")

    for x in ax.flat:
        x.grid(alpha=0.3)
        x.legend(fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out, dpi=120)
    print(f"[harness] wrote {out}")


# --------------------------------------------------------------------------- #
# Subcommands
# --------------------------------------------------------------------------- #
def cmd_run(args):
    ensure_binary(args)
    scenario = resolve_scenario(args.scenario, args.set or [])
    a = ann(scenario)
    csv_path = run_sim(args, scenario, args.name)
    d = load_trace(csv_path)
    barrier_report(d, a)
    if args.no_plot:
        return
    out = args.fig or os.path.join(args.outdir, f"{args.name}.png")
    cbf = "CBF on" if a["cbf_on"] else "CBF OFF"
    title = (f"Lon water landing -- {args.name}  "
             f"(V_app={a['V_app']}, {cbf}, V_td_max={a['V_td_max']})")
    plot_single(d, a, out, title)


def cmd_compare(args):
    ensure_binary(args)
    if not args.case:
        sys.exit("[harness] compare needs at least one --case LABEL:k=v,...")
    cases = []
    base_ann = None
    for ci, spec in enumerate(args.case):
        if ":" not in spec:
            sys.exit(f"bad --case '{spec}' (expected LABEL:k=v,k=v)")
        label, ovr = spec.split(":", 1)
        label = label.strip()
        # ';' (not ',') separates overrides so list values like [2,2,2] survive.
        overrides = [o.strip() for o in ovr.split(";") if o.strip()]
        scenario = resolve_scenario(args.scenario, overrides)
        if base_ann is None:
            base_ann = ann(scenario)
        name = f"{args.name}_{ci}_{label.replace(' ', '_')}"
        csv_path = run_sim(args, scenario, name)
        d = load_trace(csv_path)
        print(f"--- case '{label}' ---")
        barrier_report(d, ann(scenario))
        cases.append((label, d))
    if args.no_plot:
        return
    out = args.fig or os.path.join(args.outdir, f"{args.name}.png")
    plot_compare(cases, base_ann, out, f"Compare -- {args.name}")


# --------------------------------------------------------------------------- #
def add_common(p):
    p.add_argument("--scenario", default=DEFAULTS["scenario"], help="base scenario YAML")
    p.add_argument("--stab", default=DEFAULTS["stab"], help="aero .stab deck")
    p.add_argument("--aircraft", default=DEFAULTS["aircraft"], help="aircraft.yaml")
    p.add_argument("--binary", default=DEFAULTS["binary"], help="sim executable")
    p.add_argument("--build-dir", default=DEFAULTS["build_dir"])
    p.add_argument("--outdir", default=DEFAULTS["outdir"], help="run artifact dir")
    p.add_argument("--build", action="store_true", help="force rebuild before run")
    p.add_argument("--no-build", action="store_true", help="never build")
    p.add_argument("--no-plot", action="store_true", help="run only, skip the figure")
    p.add_argument("-o", "--fig", help="figure path (default <outdir>/<name>.png)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")

    pr = sub.add_parser("run", help="run one scenario and plot it")
    add_common(pr)
    pr.add_argument("--name", default="run", help="run name (artifact prefix)")
    pr.add_argument("--set", action="append", metavar="K=V",
                    help="override scenario field (dotted key), repeatable")
    pr.set_defaults(func=cmd_run)

    pc = sub.add_parser("compare", help="run N labeled cases and overlay them")
    add_common(pc)
    pc.add_argument("--name", default="compare", help="figure name (artifact prefix)")
    pc.add_argument("--case", action="append", metavar='"LABEL:K=V;K=V"',
                    help="a labeled case = base scenario + these overrides "
                         "(';'-separated so list values like [2,2,2] work)")
    pc.set_defaults(func=cmd_compare)

    args = ap.parse_args()
    if not args.cmd:
        ap.print_help()
        sys.exit(0)
    args.func(args)


if __name__ == "__main__":
    main()
