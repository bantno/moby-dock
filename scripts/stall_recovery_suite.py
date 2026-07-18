#!/usr/bin/env python3
"""Emergent stall-recovery test suite for the CBF-QP water-landing filter.

The nominal controller is a constant-gamma / constant-T_set cascade with no
flare and no stall-recovery logic. Every scenario here drives a slow approach
into (or past) the stall near the water and checks that the SAFETY FILTER
recovers and lands the aircraft anyway -- the recovery must be emergent CBF
behavior, not commanded. Each fixed case runs CBF-on (scored PASS/FAIL) and
CBF-off (informational A/B evidence: the same approach without the filter is
expected to depart or land unsafely). A recovery-floor sweep maps the minimum
altitude from which a stalled entry is still recoverable.

  scripts/stall_recovery_suite.py                 # everything -> runs/stall_suite/
  scripts/stall_recovery_suite.py --only stalled_entry --skip-sweep
  scripts/stall_recovery_suite.py --replot        # re-verdict/re-plot, no sim runs
  scripts/stall_recovery_suite.py --set cbf.stall=false   # regression injection

Exit code 0 iff every CBF-on run (all seeds) PASSes. CBF-off runs never gate,
except with --strict, where a CBF-off run that unexpectedly meets every
threshold is flagged UNEXPECTED-PASS (the case stopped stressing the filter).

Verdict philosophy: hard-row checks (impact psi1/b/residual over the active
window, no dropped HARD rows, thrust bounds) are tuning-independent safety
invariants. Sink/theta/alpha bounds are per-case envelopes frozen with >= 1.5x
margin over calibrated values. psi2 of the impact row is drift-only in the
trace and legitimately negative when the QP spends control authority -- it is
deliberately NOT checked; the enforced row is checked via its residual.
"""
import argparse
import csv as csv_mod
import os
import sys
from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np
import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from harness import DEFAULTS, apply_overrides, ann, ensure_binary, load_trace, \
    resolve_scenario, run_sim, touchdown_index

DATA = os.path.join(harness.REPO, "data")
SUITE_OUTDIR = os.path.join(harness.REPO, "runs", "stall_suite")
ALPHA_MODEL_LIMIT_DEG = 85.0   # matches kAlphaModelLimit (lon_sim.hpp)
QP_RES_TOL = 1e-3              # OSQP residual tolerance for the enforced row
TRIM_WARNING = "trim did not converge"


# --------------------------------------------------------------------------- #
# Scenario table
# --------------------------------------------------------------------------- #
@dataclass
class Thresholds:
    sink_max: float                    # touchdown sink bound [m/s], per-case
    theta_td_min_deg: float = 0.0      # touchdown pitch floor [deg]
    max_alpha_deg: float = 12.0        # departure guard (soft row may pass 9)
    t_above_alpha_lim_max: float = 2.0  # [s] alpha allowed past alpha_lim
    V_td_tol: float = 0.5              # tolerance on the V_td_max cap [m/s]
    require_touchdown: bool = True


@dataclass
class Case:
    name: str
    desc: str
    thresholds: Thresholds
    base: str = DEFAULTS["scenario"]
    overrides: tuple = ()              # dotted-key=value strings (A and B runs)
    locked: bool = False               # base yaml is self-contained: skip the
                                       # generated-case hygiene overrides
    expect_cbf_off_fail: bool = True   # B run expected to miss >= 1 threshold
    seeds: tuple = ()                  # cbf.h_meas_seed values (noisy cases)
    gust_hit_h: float = None           # 2-pass wind.t_start calibration target


# Hygiene for generated (non-locked) cases: the plant must be able to stall,
# and clean cases are deterministic (the base lon_scenario.yaml ships noise).
GENERATED_HYGIENE = (
    "stall.enabled=true",
    "stall.severity=1.0",
    "cbf.h_meas_stddev=0",
    "cbf.h_lpf_tau=0",
)
NOISY = ("cbf.h_meas_stddev=0.5", "cbf.h_lpf_tau=0.15")
NOISE_SEEDS = (0xA17B0A11, 12345, 67890)

CASES = [
    Case(
        name="slow_bleed",
        desc="Approach at ~1.1 Vs with T_set below the trim need: V bleeds "
             "along the glide and alpha creeps to the 9-deg cap near the "
             "surface; the stall row must take the elevator.",
        overrides=("V_app=10.8", "gamma_app_deg=-3.0", "initial.h0=30.0",
                   "nominal.T_set=1.2", "t_max=90.0"),
        thresholds=Thresholds(sink_max=1.5, theta_td_min_deg=0.0),
    ),
    Case(
        name="pull_into_stall",
        desc="Nominal commands an unsustainable +22 deg climb at 18 m (the "
             "CBF-off twin is the departure demo); the filter caps alpha and "
             "must land the resulting mush descent.",
        base=os.path.join(DATA, "lon_stall_pull_cbf.yaml"),
        locked=True,
        thresholds=Thresholds(sink_max=2.0, theta_td_min_deg=-5.0,
                              max_alpha_deg=14.0,
                              t_above_alpha_lim_max=5.0),
    ),
    Case(
        name="stalled_entry",
        desc="+11 deg pitch kick at 12 m puts alpha0 ~ 19.8 deg (deep past "
             "stall) with the water right below; the filter must arrest the "
             "departure and land.",
        base=os.path.join(DATA, "lon_stall_entry_cbf.yaml"),
        locked=True,
        thresholds=Thresholds(sink_max=3.5, theta_td_min_deg=-8.0,
                              max_alpha_deg=45.0,
                              t_above_alpha_lim_max=8.0),
    ),
    Case(
        name="tailwind_shear",
        desc="MIL-F-8785C tailwind gust (+4.5 m/s over 60 m) hitting at "
             "~13 m AGL robs airspeed on short final; alpha rises toward the "
             "cap with no room to dive for speed.",
        overrides=("V_app=11.5", "gamma_app_deg=-3.5", "initial.h0=25.0",
                   "nominal.T_set=1.8", "t_max=90.0", "wind.enabled=true",
                   "wind.u_amp=4.5", "wind.u_len=60.0", "wind.w_amp=0.0"),
        thresholds=Thresholds(sink_max=2.0),
        gust_hit_h=13.0,
    ),
    Case(
        name="downdraft_flare",
        desc="MIL-F-8785C downdraft (-1.5 m/s over 40 m) hitting at ~8 m AGL, "
             "inside the impact row's z_gate: forces sink exactly where the "
             "flare must happen. KNOWN GAP (encoded in V_td_tol=2.0): the "
             "energy ceiling's hdot model is air-relative, so an unmeasured "
             "downdraft makes the cap overshoot by ~1.3-1.6 m/s (V_td 15.3 at "
             "-1.5, 15.6 at -2.0, cap 14). The tolerance freezes today's gap "
             "as the regression baseline; wind-aware margins are future work.",
        overrides=("V_app=11.5", "gamma_app_deg=-3.5", "initial.h0=20.0",
                   "nominal.T_set=1.8", "t_max=90.0", "wind.enabled=true",
                   "wind.u_amp=0.0", "wind.w_amp=-1.5", "wind.w_len=40.0"),
        thresholds=Thresholds(sink_max=2.5, V_td_tol=2.0),
        gust_hit_h=8.0,
    ),
    Case(
        name="slow_bleed_noisy",
        desc="slow_bleed with the noisy AGL sensor (sigma 0.5 m, LPF 0.15 s): "
             "the recovery must not depend on perfect altitude.",
        overrides=("V_app=10.8", "gamma_app_deg=-3.0", "initial.h0=30.0",
                   "nominal.T_set=1.2", "t_max=90.0") + NOISY,
        thresholds=Thresholds(sink_max=2.0, theta_td_min_deg=-2.0),
        seeds=NOISE_SEEDS,
    ),
    Case(
        name="tailwind_shear_noisy",
        desc="tailwind_shear with the noisy AGL sensor.",
        overrides=("V_app=11.5", "gamma_app_deg=-3.5", "initial.h0=25.0",
                   "nominal.T_set=1.8", "t_max=90.0", "wind.enabled=true",
                   "wind.u_amp=4.5", "wind.u_len=60.0", "wind.w_amp=0.0") + NOISY,
        thresholds=Thresholds(sink_max=2.5, theta_td_min_deg=-2.0),
        seeds=NOISE_SEEDS,
        gust_hit_h=13.0,
    ),
]

# Recovery-floor sweep: stalled-entry template over entry altitude x pitch kick.
SWEEP_H0 = (5.0, 8.0, 12.0, 16.0, 20.0)
SWEEP_KICK = (8.0, 12.0, 16.0)   # deg; 12/16 exceed the +-20 deg VSPAERO grid
                                 # at entry (extrapolated + Viterna-blended)
SWEEP_SAFE_SINK = 1.5            # recovered-safe touchdown sink bound [m/s]
SWEEP_CRASH_SINK = 4.0           # above this = crash even if rows held


# --------------------------------------------------------------------------- #
# Metrics + verdicts
# --------------------------------------------------------------------------- #
def scenario_params(scenario):
    """Verdict-relevant knobs with the lon_sim.cpp defaults."""
    cbf = scenario.get("cbf", {}) or {}
    return dict(
        dt=scenario.get("dt", 0.01),
        t_max=scenario.get("t_max", 60.0),
        z_gate=cbf.get("z_gate", 10.0),
        tau_keel_deg=cbf.get("tau_keel_deg", 0.0),
        Tmax=cbf.get("Tmax", 50.0),
        V_td_max=cbf.get("V_td_max", 14.0),
        alpha_lim_deg=(cbf.get("alpha_stall_deg", 11.0)
                       - cbf.get("stall_margin_deg", 2.0)),
        h_noseup=cbf.get("h_noseup", 3.0),
        cbf_on=bool(cbf.get("enabled", True))
               and bool(scenario.get("cbf_enabled", True)),
    )


def compute_metrics(d, scenario, log_text):
    """All verdict inputs from one trace (columns are the current lon_sim CSV
    schema). The trace ends at touchdown, so full-trace minima are
    up-to-touchdown minima."""
    p = scenario_params(scenario)
    i = touchdown_index(d)
    reached = bool(d["h"][i] <= 0.0)
    alpha = d["alpha_deg"]

    # Impact HARD row over its ACTIVE window only (h < z_gate, descending,
    # positive trim) -- b_impact goes deeply negative aloft where the row is
    # not assembled, so a full-trace min would be meaningless.
    mask = (d["h"] < p["z_gate"]) & (d["sink"] > 0.0) & \
           (d["theta_deg"] > p["tau_keel_deg"])
    def masked_min(col):
        v = d[col][mask]
        return float(np.nanmin(v)) if len(v) else float("inf")

    # "Stuck beyond the limit" excursion: riding the boundary at alpha_lim
    # (b_stall ~ 0) is healthy soft-row behavior, so the excursion clock only
    # runs past alpha_lim + 0.5 deg.
    alpha_over = p["alpha_lim_deg"] + 0.5
    j_up = np.where(alpha > p["alpha_lim_deg"])[0]
    arrest_h = float("nan")
    if len(j_up):
        j0 = j_up[0]
        j_ar = np.where(d["sink"][j0:] <= 0.05)[0]
        if len(j_ar):
            arrest_h = float(d["h"][j0 + j_ar[0]])
        elif reached:
            arrest_h = 0.0     # never arrested before the water

    return dict(
        reached=reached,
        t_td=float(d["t"][i]), sink_td=float(d["sink"][i]),
        V_td=float(d["V"][i]), theta_td=float(d["theta_deg"][i]),
        max_alpha=float(np.max(alpha)),
        max_abs_alpha=float(np.max(np.abs(alpha))),
        t_above_alpha_lim=float(p["dt"] * np.sum(alpha > alpha_over)),
        arrest_h=arrest_h,
        min_V=float(np.min(d["V"])),
        steps_hard_dropped=int(np.sum(d["hard_dropped"])),
        steps_rows_dropped=int(np.sum(d["n_rows_dropped"] > 0)),
        recoveries=int(np.sum(d["recovered"])),
        min_psi1_imp_active=masked_min("psi1_imp"),
        min_b_impact_active=masked_min("b_impact"),
        min_res_imp_active=masked_min("res_imp"),
        min_T=float(np.min(d["T"])),
        max_T=float(np.max(d["T"])),
        min_b_stall=float(np.min(d["b_stall"])),
        min_dCL_stall=float(np.min(d["dCL_stall"])),
        out_of_model=bool(np.max(np.abs(alpha)) >= ALPHA_MODEL_LIMIT_DEG),
        trim_failed=TRIM_WARNING in log_text,
        noisy=(scenario.get("cbf", {}) or {}).get("h_meas_stddev", 0.0) > 0.0,
        V_td_max=p["V_td_max"], Tmax=p["Tmax"],
        alpha_lim_deg=p["alpha_lim_deg"],
    )


def judge(m, th):
    """Return (verdict, [reasons]). INVALID trumps FAIL trumps PASS."""
    if m["trim_failed"]:
        return "INVALID", ["initial-condition trim did not converge"]
    if m["out_of_model"]:
        return "INVALID", [f"|alpha| reached {m['max_abs_alpha']:.0f} deg "
                           f">= {ALPHA_MODEL_LIMIT_DEG:.0f} (plant invalid)"]
    fails = []
    def need(ok, msg):
        if not ok:
            fails.append(msg)
    # Tuning-independent safety invariants.
    need(m["steps_hard_dropped"] == 0,
         f"{m['steps_hard_dropped']} step(s) dropped a HARD row")
    need(m["min_psi1_imp_active"] >= -1e-6,
         f"impact psi1 min {m['min_psi1_imp_active']:.4f} < 0 (active window)")
    need(m["min_b_impact_active"] >= -1e-6,
         f"b_impact min {m['min_b_impact_active']:.4f} < 0 (active window)")
    # The enforced-row residual is a MEASURED-state guarantee: with a noisy
    # altitude sensor the row the filter enforced (on h_filt) differs from the
    # true-state row recomputed here, so the residual check applies only to
    # noise-free runs; noisy runs are still held to psi1/b on the true state
    # and the touchdown envelope.
    if not m["noisy"]:
        need(m["min_res_imp_active"] >= -QP_RES_TOL,
             f"impact row residual min {m['min_res_imp_active']:.4f} < 0")
    need(m["min_T"] >= -1e-6, f"thrust min {m['min_T']:.3f} < 0")
    need(m["max_T"] <= m["Tmax"] + 1e-6,
         f"thrust max {m['max_T']:.2f} > Tmax {m['Tmax']}")
    # Per-case touchdown envelope.
    if th.require_touchdown:
        need(m["reached"], "no touchdown within t_max")
    if m["reached"]:
        need(m["sink_td"] <= th.sink_max,
             f"sink {m['sink_td']:.2f} > {th.sink_max} m/s")
        need(m["V_td"] <= m["V_td_max"] + th.V_td_tol,
             f"V_td {m['V_td']:.2f} > cap {m['V_td_max']}+{th.V_td_tol}")
        need(m["theta_td"] >= th.theta_td_min_deg,
             f"theta_td {m['theta_td']:.1f} < {th.theta_td_min_deg} deg")
    need(m["max_alpha"] <= th.max_alpha_deg,
         f"max alpha {m['max_alpha']:.1f} > {th.max_alpha_deg} deg")
    need(m["t_above_alpha_lim"] <= th.t_above_alpha_lim_max,
         f"alpha above alpha_lim for {m['t_above_alpha_lim']:.2f} s "
         f"> {th.t_above_alpha_lim_max}")
    return ("PASS" if not fails else "FAIL"), fails


def classify_sweep(m):
    """Recovery-floor outcome classes, by mechanism:
    - crash: physically hard impact, or the safe set itself was left
      (psi1/b_impact negative over the active window / a HARD row dropped).
    - recovered-degraded: landed inside the sink envelope and the safe set
      was never left, but the QP went through best-effort feasibility
      recovery (enforced-row residual < 0): the hard GUARANTEE broke
      transiently even though the OUTCOME stayed safe.
    - recovered-safe / hard-landing: clean enforcement, split on touchdown
      sink / V_td cap.
    - aloft: no touchdown within t_max (mushing above the water)."""
    if m["trim_failed"]:
        return "invalid"
    if m["out_of_model"]:
        return "out-of-model"
    if not m["reached"]:
        return "aloft"
    invariance_ok = (m["steps_hard_dropped"] == 0
                     and m["min_psi1_imp_active"] >= -1e-6
                     and m["min_b_impact_active"] >= -1e-6)
    if not invariance_ok or m["sink_td"] > SWEEP_CRASH_SINK:
        return "crash"
    if m["min_res_imp_active"] < -QP_RES_TOL:
        return "recovered-degraded"
    if m["sink_td"] <= SWEEP_SAFE_SINK and m["V_td"] <= m["V_td_max"] + 0.5:
        return "recovered-safe"
    return "hard-landing"


# --------------------------------------------------------------------------- #
# Running
# --------------------------------------------------------------------------- #
def sim_args(g, outdir):
    return SimpleNamespace(binary=g.binary, stab=g.stab, aircraft=g.aircraft,
                           outdir=outdir, build=False, no_build=True,
                           build_dir=g.build_dir)


def build_scenario(case, g, cbf_on, seed=None):
    ovr = list(case.overrides)
    if not case.locked:
        ovr = list(GENERATED_HYGIENE) + ovr
    if seed is not None:
        ovr.append(f"cbf.h_meas_seed={seed}")
    if not cbf_on:
        ovr.append("cbf_enabled=false")
    ovr += list(g.set or [])       # global regression-injection overrides last
    return resolve_scenario(case.base, ovr)


def read_log(outdir, name):
    path = os.path.join(outdir, f"{name}.log")
    return open(path).read() if os.path.exists(path) else ""


def run_or_load(g, case_dir, scenario, name):
    """Run the sim (or, with --replot, reuse existing artifacts)."""
    csv_path = os.path.join(case_dir, f"{name}.csv")
    scn_path = os.path.join(case_dir, f"{name}.resolved.yaml")
    if g.replot:
        if not os.path.exists(csv_path):
            sys.exit(f"[suite] --replot but no artifact {csv_path}")
        with open(scn_path) as f:
            scenario = yaml.safe_load(f)
    else:
        run_sim(sim_args(g, case_dir), scenario, name)
    d = load_trace(csv_path)
    return d, scenario, read_log(case_dir, name)


def calibrate_gust_t_start(g, case, case_dir, scenario):
    """2-pass gust timing: run the identical scenario with the gust disabled,
    read the time the aircraft first descends through gust_hit_h, and start
    the gust there (pre-onset trajectories are identical, so this is exact)."""
    cal = dict(scenario)
    cal = yaml.safe_load(yaml.safe_dump(cal))       # deep copy
    apply_overrides(cal, ["wind.enabled=false"])
    d, _, _ = run_or_load(g, case_dir, cal, "cal_no_wind")
    below = np.where(d["h"] <= case.gust_hit_h)[0]
    if not len(below):
        print(f"[suite] WARNING {case.name}: never reached gust_hit_h="
              f"{case.gust_hit_h} m; gust starts at t=0")
        return 0.0
    return float(d["t"][below[0]])


def run_case(g, case):
    """Run a fixed case (A = CBF-on scored, per seed; B = CBF-off informational).
    Returns a list of per-run record dicts."""
    case_dir = os.path.join(g.outdir, case.name)
    os.makedirs(case_dir, exist_ok=True)
    records = []

    scenario_a = build_scenario(case, g, cbf_on=True)
    t_start = None
    if case.gust_hit_h is not None:
        t_start = calibrate_gust_t_start(g, case, case_dir, scenario_a)
        print(f"[suite] {case.name}: gust t_start={t_start:.2f} s "
              f"(hits h~{case.gust_hit_h} m)")

    def one(name, cbf_on, seed=None):
        scn = build_scenario(case, g, cbf_on=cbf_on, seed=seed)
        if t_start is not None:
            apply_overrides(scn, [f"wind.t_start={t_start}"])
        d, scn, log = run_or_load(g, case_dir, scn, name)
        m = compute_metrics(d, scn, log)
        rec = dict(case=case.name, run=name, cbf_on=cbf_on, seed=seed, **m)
        if cbf_on:
            rec["verdict"], rec["fail_reasons"] = judge(m, case.thresholds)
        else:
            v, reasons = judge(m, case.thresholds)
            rec["verdict"] = "INFO-UNSAFE" if v != "PASS" else "UNEXPECTED-PASS"
            rec["fail_reasons"] = reasons
        records.append(rec)
        return d

    seeds = case.seeds or (None,)
    traces_a = []
    for si, seed in enumerate(seeds):
        name = "on" if seed is None else f"on_s{si}"
        traces_a.append((name, one(name, cbf_on=True, seed=seed)))
    d_off = one("off", cbf_on=False, seed=seeds[0])

    if not g.no_plot:
        plot_case_ab(case, traces_a, d_off,
                     ann(build_scenario(case, g, cbf_on=True)),
                     os.path.join(case_dir, f"{case.name}_ab.png"))
    return records


def run_sweep(g):
    """Recovery-floor grid from the locked stalled-entry template."""
    sweep_dir = os.path.join(g.outdir, "sweep")
    os.makedirs(sweep_dir, exist_ok=True)
    base = os.path.join(DATA, "lon_stall_entry_cbf.yaml")
    rows = []
    for h0 in SWEEP_H0:
        for kick in SWEEP_KICK:
            name = f"h{int(h0):02d}_k{int(kick):02d}"
            ovr = [f"initial.h0={h0}", f"initial.dtheta_deg={kick}",
                   "t_max=90.0"] + list(g.set or [])
            scn = resolve_scenario(base, ovr)
            d, scn, log = run_or_load(g, sweep_dir, scn, name)
            m = compute_metrics(d, scn, log)
            rows.append(dict(name=name, h0=h0, kick=kick,
                             outcome=classify_sweep(m), **m))
    if not g.no_plot:
        plot_sweep_heatmap(rows, os.path.join(sweep_dir, "sweep_heatmap.png"))
    write_csv(os.path.join(sweep_dir, "sweep.csv"), rows)
    return rows


# --------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------- #
def plot_case_ab(case, traces_a, d_off, a, out):
    """CBF-on (all seeds) vs CBF-off overlay; the alpha panel is the
    load-bearing trace for this suite."""
    plt = harness._mpl()
    fig, ax = plt.subplots(2, 3, figsize=(17, 9))
    fig.suptitle(f"stall-recovery: {case.name} -- CBF-on (scored) vs CBF-off "
                 f"(informational)", fontsize=13, weight="bold")
    alpha_lim = a["alpha_stall_deg"] - a["stall_margin_deg"]

    def draw(d, label, color, ls="-"):
        i = touchdown_index(d)
        ax[0, 0].plot(d["t"], d["h"], color=color, ls=ls, label=label)
        ax[0, 1].plot(d["t"], d["alpha_deg"], color=color, ls=ls, label=label)
        ax[0, 2].plot(d["h"], d["sink"], color=color, ls=ls, label=label)
        ax[1, 0].plot(d["t"], d["V"], color=color, ls=ls, label=label)
        ax[1, 1].plot(d["t"], d["theta_deg"], color=color, ls=ls, label=label)
        ax[1, 1].plot(d["t"], d["gamma_deg"], color=color, ls=ls, alpha=0.45)
        ax[1, 2].plot(d["t"], d["b_impact"], color=color, ls=ls,
                      label=f"{label} (td sink {d['sink'][i]:.2f})")

    colors = ["tab:blue", "tab:cyan", "tab:green"]
    for k, (name, d) in enumerate(traces_a):
        draw(d, f"CBF on ({name})", colors[k % len(colors)])
    draw(d_off, "CBF off", "tab:red", ls="--")

    ax[0, 0].axhline(0, color="k", lw=0.8, ls=":")
    ax[0, 0].set(xlabel="t [s]", ylabel="h [m]", title="Altitude")
    ax[0, 1].axhline(alpha_lim, color="tab:red", ls=":", lw=1.0,
                     label=f"alpha_lim={alpha_lim:.0f}")
    ax[0, 1].axhline(a["alpha_stall_deg"], color="k", ls=":", lw=0.8,
                     label=f"wing stall={a['alpha_stall_deg']:.0f}")
    ax[0, 1].set(xlabel="t [s]", ylabel="alpha [deg]",
                 title="AoA vs stall barrier (the emergent recovery)")
    ax[0, 2].set(xlabel="h [m]", ylabel="sink [m/s]", title="Sink vs altitude")
    ax[0, 2].set_xlim(12, 0)   # the near-water endgame is the story
    ax[1, 0].axhline(a["V_td_max"], color="k", ls=":", lw=0.9,
                     label=f"V_td_max={a['V_td_max']}")
    ax[1, 0].set(xlabel="t [s]", ylabel="V [m/s]", title="Airspeed")
    ax[1, 1].axhline(a["theta_min_deg"], color="tab:brown", ls=":", lw=0.9,
                     label=f"theta_min={a['theta_min_deg']:.0f}")
    ax[1, 1].set(xlabel="t [s]", ylabel="angle [deg]",
                 title=r"$\theta$ (solid) / $\gamma$ (faint)")
    ax[1, 2].axhline(0, color="k", lw=1.0, ls="--")
    ax[1, 2].set(xlabel="t [s]", ylabel="b_impact",
                 title="Impact barrier (HARD)")
    for x in ax.flat:
        x.grid(alpha=0.3)
        x.legend(fontsize=7)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out, dpi=120)
    print(f"[suite] wrote {out}")


SWEEP_ORDER = ["recovered-safe", "recovered-degraded", "hard-landing",
               "crash", "out-of-model", "aloft", "invalid"]
SWEEP_COLORS = ["#2e9e4f", "#e07b28", "#e0b13d", "#c8402f", "#8a8a8a",
                "#4a7fd4", "#000000"]


def plot_sweep_heatmap(rows, out):
    plt = harness._mpl()
    from matplotlib.colors import BoundaryNorm, ListedColormap
    h0s, kicks = sorted({r["h0"] for r in rows}), sorted({r["kick"] for r in rows})
    grid = np.zeros((len(h0s), len(kicks)))
    for r in rows:
        grid[h0s.index(r["h0"]), kicks.index(r["kick"])] = \
            SWEEP_ORDER.index(r["outcome"])
    fig, ax = plt.subplots(figsize=(7, 5.5))
    cmap = ListedColormap(SWEEP_COLORS)
    norm = BoundaryNorm(np.arange(-0.5, len(SWEEP_ORDER)), cmap.N)
    ax.pcolormesh(np.arange(len(kicks) + 1), np.arange(len(h0s) + 1),
                  grid, cmap=cmap, norm=norm, edgecolors="w", lw=1.5)
    letters = {"recovered-safe": "S", "recovered-degraded": "D",
               "hard-landing": "H", "crash": "C", "out-of-model": "O",
               "aloft": "A", "invalid": "!"}
    for r in rows:
        i, j = h0s.index(r["h0"]), kicks.index(r["kick"])
        txt = letters[r["outcome"]]
        if r["reached"]:
            txt += f"\n{r['sink_td']:.1f}"
        ax.text(j + 0.5, i + 0.5, txt, ha="center", va="center",
                color="w", fontsize=9, weight="bold")
    ax.set_xticks(np.arange(len(kicks)) + 0.5, [f"+{k:.0f}" for k in kicks])
    ax.set_yticks(np.arange(len(h0s)) + 0.5, [f"{h:.0f}" for h in h0s])
    ax.set(xlabel="pitch kick dtheta [deg]", ylabel="entry altitude h0 [m]",
           title="Recovery floor: stalled entry, CBF on\n"
                 "(letter = outcome class, number = touchdown sink m/s)")
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c) for c in SWEEP_COLORS]
    ax.legend(handles, SWEEP_ORDER, fontsize=7, loc="center left",
              bbox_to_anchor=(1.02, 0.5))
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    print(f"[suite] wrote {out}")


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #
CSV_FIELDS = ["case", "run", "cbf_on", "seed", "verdict", "fail_reasons",
              "reached", "t_td", "sink_td", "V_td", "theta_td", "max_alpha",
              "t_above_alpha_lim", "arrest_h", "min_V", "recoveries",
              "steps_rows_dropped", "steps_hard_dropped",
              "min_psi1_imp_active", "min_b_impact_active",
              "min_res_imp_active", "min_T", "max_T", "min_b_stall",
              "min_dCL_stall", "out_of_model", "trim_failed"]


def write_csv(path, rows):
    keys = list(rows[0].keys()) if rows else []
    fields = [k for k in CSV_FIELDS if k in keys] + \
             [k for k in keys if k not in CSV_FIELDS]
    with open(path, "w", newline="") as f:
        w = csv_mod.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            r = dict(r)
            if isinstance(r.get("fail_reasons"), list):
                r["fail_reasons"] = " | ".join(r["fail_reasons"])
            w.writerow(r)
    print(f"[suite] wrote {path}")


def write_summary_md(path, records, sweep_rows, cases_run):
    by_case = {}
    for r in records:
        by_case.setdefault(r["case"], []).append(r)
    lines = ["# Stall-recovery suite report", ""]
    lines += ["The nominal controller has no flare and no stall-recovery "
              "logic; every safe landing below is emergent CBF behavior. "
              "CBF-off rows are the informational A/B evidence.", ""]
    lines += ["| case | run | verdict | td t [s] | sink [m/s] | V_td | "
              "theta_td [deg] | max alpha | recoveries | notes |",
              "|---|---|---|---|---|---|---|---|---|---|"]
    for case in cases_run:
        for r in by_case.get(case.name, []):
            td = (f"{r['t_td']:.1f}" if r["reached"] else "--")
            sink = (f"{r['sink_td']:.2f}" if r["reached"] else "--")
            vtd = (f"{r['V_td']:.2f}" if r["reached"] else "--")
            th = (f"{r['theta_td']:.1f}" if r["reached"] else "--")
            notes = "; ".join(r["fail_reasons"]) or "-"
            lines.append(f"| {r['case']} | {r['run']} | **{r['verdict']}** | "
                         f"{td} | {sink} | {vtd} | {th} | "
                         f"{r['max_alpha']:.1f} | {r['recoveries']} | {notes} |")
    lines += ["", "## Case descriptions", ""]
    for case in cases_run:
        lines.append(f"- **{case.name}**: {case.desc}")
    if sweep_rows:
        lines += ["", "## Recovery floor (stalled entry, CBF on)", "",
                  "Rows = entry altitude h0, columns = pitch kick. "
                  "Cell = outcome (touchdown sink m/s).", ""]
        kicks = sorted({r["kick"] for r in sweep_rows})
        lines.append("| h0 \\ kick | " +
                     " | ".join(f"+{k:.0f} deg" for k in kicks) + " |")
        lines.append("|---|" + "---|" * len(kicks))
        for h0 in sorted({r["h0"] for r in sweep_rows}):
            cells = []
            for k in kicks:
                r = next(x for x in sweep_rows
                         if x["h0"] == h0 and x["kick"] == k)
                c = r["outcome"]
                if r["reached"]:
                    c += f" ({r['sink_td']:.2f})"
                cells.append(c)
            lines.append(f"| {h0:.0f} m | " + " | ".join(cells) + " |")
        lines.append("")
        for k in kicks:
            col = [r for r in sweep_rows if r["kick"] == k]
            safe = sorted(r["h0"] for r in col
                          if r["outcome"] == "recovered-safe")
            degr = sorted(r["h0"] for r in col
                          if r["outcome"] == "recovered-degraded")
            floor = f"{safe[0]:.0f} m" if safe else "none"
            dtag = (f"; degraded (outcome safe, guarantee broke) down to "
                    f"{degr[0]:.0f} m" if degr else "")
            mono = all(r["outcome"] == "recovered-safe"
                       for r in col if safe and r["h0"] >= safe[0])
            tag = "" if mono else "  (NON-MONOTONE above the floor -- inspect)"
            lines.append(f"- kick +{k:.0f} deg: clean recovery floor = "
                         f"{floor}{dtag}{tag}")
        lines += ["", "See `sweep/sweep_heatmap.png` and `sweep/sweep.csv`."]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"[suite] wrote {path}")


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", help="comma-separated case names (default all)")
    ap.add_argument("--skip-sweep", action="store_true")
    ap.add_argument("--replot", action="store_true",
                    help="recompute verdicts/plots from existing artifacts")
    ap.add_argument("--strict", action="store_true",
                    help="also fail on UNEXPECTED-PASS of a CBF-off run")
    ap.add_argument("--set", action="append", metavar="K=V", default=[],
                    help="global scenario override (regression injection)")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--outdir", default=SUITE_OUTDIR)
    ap.add_argument("--stab", default=DEFAULTS["stab"])
    ap.add_argument("--aircraft", default=DEFAULTS["aircraft"])
    ap.add_argument("--binary", default=DEFAULTS["binary"])
    ap.add_argument("--build-dir", dest="build_dir", default=DEFAULTS["build_dir"])
    ap.add_argument("--build", action="store_true", help="force rebuild")
    ap.add_argument("--no-build", action="store_true")
    g = ap.parse_args()

    if not g.replot:
        ensure_binary(g)
    os.makedirs(g.outdir, exist_ok=True)

    names = [n.strip() for n in g.only.split(",")] if g.only else None
    cases_run = [c for c in CASES if names is None or c.name in names]
    if names:
        unknown = set(names) - {c.name for c in CASES}
        if unknown:
            sys.exit(f"[suite] unknown case(s): {sorted(unknown)} "
                     f"(have: {[c.name for c in CASES]})")

    records = []
    for case in cases_run:
        print(f"\n[suite] ===== case {case.name} =====")
        records += run_case(g, case)

    sweep_rows = []
    if not g.skip_sweep:
        print("\n[suite] ===== recovery-floor sweep =====")
        sweep_rows = run_sweep(g)

    write_csv(os.path.join(g.outdir, "summary.csv"), records)
    write_summary_md(os.path.join(g.outdir, "summary.md"),
                     records, sweep_rows, cases_run)

    scored = [r for r in records if r["cbf_on"]]
    failed = [r for r in scored if r["verdict"] != "PASS"]
    surprises = [r for r in records
                 if not r["cbf_on"] and r["verdict"] == "UNEXPECTED-PASS"
                 and next(c for c in cases_run
                          if c.name == r["case"]).expect_cbf_off_fail]
    print("\n[suite] ================= VERDICT =================")
    for r in scored:
        tag = "" if r["verdict"] == "PASS" else \
            "  <-- " + "; ".join(r["fail_reasons"])
        print(f"  {r['case']:22s} {r['run']:6s} {r['verdict']}{tag}")
    for r in surprises:
        print(f"  {r['case']:22s} {r['run']:6s} UNEXPECTED-PASS "
              f"(CBF-off met every threshold -- case no longer stresses "
              f"the filter)")
    ok = not failed and not (g.strict and surprises)
    print(f"[suite] {'ALL PASS' if ok else 'FAILURES PRESENT'} "
          f"({len(scored) - len(failed)}/{len(scored)} scored runs pass)")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
