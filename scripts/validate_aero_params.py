#!/usr/bin/env python3
"""Validate the placeholder aero/thrust parameters against a PX4 flight log.

Targets the TODO-flagged placeholders in data/aircraft.yaml:
    parasite_CD0 (0.030)     thrust.T_static (50 N)     thrust.k_v (0.05 N/(m/s))

METHOD. No pitot and no airborne idle-glide exist in the log, so drag cannot
be read off directly; and the flown throttle band in the clean window
(delta in [0.23, 0.38]) is too narrow to identify a throttle->thrust map
(delta and delta*V are collinear there: naive joint fits return negative
T_static). The primary identification therefore uses BATTERY POWER as the
thrust proxy, which sidesteps the throttle-map form entirely:

  air-relative specific-energy rate (exact for steady wind, zero vertical wind):
      m*(V*Vdot + g*hdot)/V = T - D
  thrust proxy:   T = eta_p * P_battery / V     (one lumped propulsive
                  efficiency: prop x motor x ESC; fitted, must land ~0.3-0.6)
  drag model:     D = qbar*S*(CD0 + k_ind*CL^2)

  CL is measured PER SAMPLE from the accelerometer load factor
      n = -f_z/g  (vehicle_acceleration specific force),  CL = n*m*g/(qbar*S),
  and k_ind comes from the VSPAERO deck's own polar -- so the drag model needs
  NO angle-of-attack estimate. (alpha = theta - gamma_air is UNUSABLE here:
  with no vane its scatter (~3-4 deg std) is as large as its signal; the
  CL-vs-alpha panel documents this.)

  Everything is LINEAR in p = [eta_p, CD0]:
      y := m*(V*Vdot + g*hdot)/V + qbar*S*k_ind*CL^2
         = [P/V, -qbar*S] . p

THRUST-MAP CHARACTERIZATION (secondary): with eta_p fitted, T_est =
eta_p*P/V is plotted against throttle across the WHOLE flight -- including
the high-throttle takeoff, where P is wind-immune (airspeed there uses
ground speed +- the ~2 m/s wind the estimator later converged to). This
tests the aircraft.yaml form T = delta*(T_static - k_v*V) directly; a
delta-squared curve is fitted alongside, since cruise (~6 N at delta 0.3)
vs takeoff (~tens of N at delta 0.7) is far off the linear map.

DATA WINDOW for the energy fit: t in [45, 133] s -- the EKF wind estimate
(sideslip-fused, no pitot) starts at zero and only converges by ~40 s, which
silently corrupts airspeed for the entire takeoff.

CAVEATS: mass from aircraft.yaml (3.6139 kg) -- not identifiable here, errors
scale T_static/CD0 one-for-one; rho from the log (1.193), not config 1.225.

Usage:
  scripts/validate_aero_params.py [data/Takeoff_2_18_07_24.ulg] [figures/flightlog_param_validation.png]
"""
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import pyulog
except ImportError:
    sys.exit("pyulog is required: pip install pyulog")

ROOT = os.path.join(os.path.dirname(__file__), os.pardir)

# --- config values under test (aircraft.yaml) + geometry (.stab) --------------
MASS = 3.6139           # kg (aircraft.yaml; NOT identifiable from this log)
G = 9.80665
SREF = 0.425            # m^2
CFG_CD0 = 0.030
CFG_TSTATIC = 50.0
CFG_KV = 0.05

# --- sample-selection gates ----------------------------------------------------
T_WIN = (45.0, 133.0)   # wind-estimator-converged airborne window [s]
H_MIN = 3.0             # m AGL
PHI_MAX = np.radians(15.0)
Q_MAX = 0.3             # rad/s pitch-rate quasi-steady gate
V_MIN = 8.0             # m/s
DT = 0.05               # uniform grid [s]
SMOOTH_S = 1.0          # smoothing window for rates [s]


def smooth(x, n):
    if n <= 1:
        return x
    return np.convolve(x, np.ones(n) / n, mode="same")


class Topics:
    def __init__(self, ulog):
        self.ulog = ulog
        self.t0 = self.raw("vehicle_local_position").data["timestamp"][0]

    def raw(self, name):
        for d in self.ulog.data_list:
            if d.name == name and d.multi_id == 0:
                return d
        return None

    def series(self, name, field):
        d = self.raw(name)
        if d is None or field not in d.data:
            return None, None
        t = (d.data["timestamp"] - self.t0) / 1e6
        v = np.asarray(d.data[field], float)
        o = np.argsort(t, kind="stable")
        t, v = t[o], v[o]
        keep = np.concatenate(([True], np.diff(t) > 0))
        return t[keep], v[keep]


def resample(t, v, grid):
    if t is None:
        return np.full_like(grid, np.nan)
    out = np.interp(grid, t, v)
    out[(grid < t[0]) | (grid > t[-1])] = np.nan
    return out


def load_states(path, t_win=T_WIN):
    """Uniform-grid state table (plot_log_states.py recipes + accel/power)."""
    T = Topics(pyulog.ULog(path))
    grid = np.arange(t_win[0], t_win[1], DT)
    nsm = max(1, int(round(SMOOTH_S / DT)))

    tx, x = T.series("vehicle_local_position", "x")
    _, y = T.series("vehicle_local_position", "y")
    _, z = T.series("vehicle_local_position", "z")
    win = max(1, int(round(0.4 / np.median(np.diff(tx)))))
    vN = resample(tx, smooth(np.gradient(x, tx), win), grid)
    vE = resample(tx, smooth(np.gradient(y, tx), win), grid)
    vD = resample(tx, smooth(np.gradient(z, tx), win), grid)
    h = resample(tx, -z, grid)

    twn, wn = T.series("estimator_wind", "windspeed_north")
    twe, we = T.series("estimator_wind", "windspeed_east")
    wind_n, wind_e = resample(twn, wn, grid), resample(twe, we, grid)

    vNa, vEa = vN - wind_n, vE - wind_e
    V_air = smooth(np.sqrt(vNa**2 + vEa**2 + vD**2), nsm)
    gamma = np.arctan2(-vD, np.hypot(vNa, vEa))
    Vdot = np.gradient(V_air, grid)
    hdot = -vD

    tq, q0 = T.series("vehicle_attitude", "q[0]")
    _, q1 = T.series("vehicle_attitude", "q[1]")
    _, q2 = T.series("vehicle_attitude", "q[2]")
    _, q3 = T.series("vehicle_attitude", "q[3]")
    theta = resample(tq, np.arcsin(np.clip(2 * (q0 * q2 - q3 * q1), -1, 1)), grid)
    phi = resample(tq, np.arctan2(2 * (q0 * q1 + q2 * q3),
                                  1 - 2 * (q1 * q1 + q2 * q2)), grid)
    tw, wy = T.series("vehicle_angular_velocity", "xyz[1]")
    qrate = resample(tw, wy, grid) if tw is not None else np.full_like(grid, np.nan)

    ta, az = T.series("vehicle_acceleration", "xyz[2]")   # body-z specific force
    n_load = resample(ta, smooth(-az / G, max(1, int(round(0.5 / np.median(np.diff(ta))))))
                      if ta is not None else None, grid)

    tm, c0 = T.series("actuator_motors", "control[0]")
    delta = resample(tm, c0, grid)

    _, rho = T.series("vehicle_air_data", "rho")
    rho_mean = float(np.nanmean(rho)) if rho is not None else 1.225

    tb, Ib = T.series("battery_status", "current_a")
    _, Vb = T.series("battery_status", "voltage_v")
    P_elec = resample(tb, Ib * Vb, grid)

    return dict(t=grid, h=h, V=V_air, Vdot=Vdot, hdot=hdot, gamma=gamma,
                theta=theta, phi=phi, qrate=qrate, n_load=n_load, delta=delta,
                rho=rho_mean, P_elec=P_elec)


def load_deck(stab_path):
    """(alpha [rad], CL_deck, CD_deck) at beta=0, lowest Mach: Base_Aero
    CFx/CFz rotated analytically to wind axes (never the pre-rotated cols)."""
    txt = open(stab_path).read()
    rows = []
    for b in re.split(r"\*{5,}", txt):
        am = re.search(r"AoA_\s+(-?\d+\.?\d*)", b)
        bm = re.search(r"Beta_\s+(-?\d+\.?\d*)", b)
        mm = re.search(r"Mach_\s+(-?\d+\.?\d*)", b)
        base = re.search(r"Base_Aero\s+\S+\s+\S+\s+(.+)", b)
        if not (am and bm and mm and base) or abs(float(bm.group(1))) > 1e-9:
            continue
        v = base.group(1).split()  # CFx CFy CFz ...
        rows.append((float(mm.group(1)), float(am.group(1)),
                     float(v[0]), float(v[2])))
    m0 = min(r[0] for r in rows)
    rows = sorted([r for r in rows if abs(r[0] - m0) < 1e-9], key=lambda r: r[1])
    a = np.radians([r[1] for r in rows])
    cfx = np.array([r[2] for r in rows])
    cfz = np.array([r[3] for r in rows])
    return a, -cfx * np.sin(a) + cfz * np.cos(a), cfx * np.cos(a) + cfz * np.sin(a)


def deck_induced_factor(cl_deck, cd_deck):
    """k_ind of the deck's own polar, CD_deck ~ c0 + k_ind*CL^2 (CL in [0,1.2])."""
    m = (cl_deck > 0.0) & (cl_deck < 1.2)
    k, c0 = np.polyfit(cl_deck[m] ** 2, cd_deck[m], 1)
    return float(k), float(c0)


def fit(X, y, clip=None):
    """LS fit; clip=k applies two rounds of k-sigma residual clipping (the
    energy-rate signal has sparse maneuver-transient outliers)."""
    keep = np.ones(len(y), bool)
    for _ in range(3 if clip else 1):
        p, *_ = np.linalg.lstsq(X[keep], y[keep], rcond=None)
        r = y[keep] - X[keep] @ p
        if clip:
            newkeep = np.ones(len(y), bool)
            newkeep[keep] = np.abs(r) < clip * np.std(r)
            newkeep &= keep
            if newkeep.sum() == keep.sum():
                break
            keep = newkeep
    dof = max(1, keep.sum() - X.shape[1])
    cov = np.linalg.inv(X[keep].T @ X[keep]) * (r @ r / dof)
    se = np.sqrt(np.diag(cov))
    return p, se, cov / np.outer(se, se), r, keep


def main():
    ulg = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "data", "Takeoff_2_18_07_24.ulg")
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "figures", "flightlog_param_validation.png")
    S = load_states(ulg)
    a_deck, cl_deck, cd_deck = load_deck(os.path.join(ROOT, "data", "AHAB_combined.stab"))
    k_ind, cd_min_deck = deck_induced_factor(cl_deck, cd_deck)

    qbarS = 0.5 * S["rho"] * S["V"] ** 2 * SREF
    cl = S["n_load"] * MASS * G / qbarS          # accelerometer load factor -> CL
    sel = (np.isfinite(S["V"]) & np.isfinite(S["delta"]) & np.isfinite(cl) &
           (S["h"] > H_MIN) & (np.abs(S["phi"]) < PHI_MAX) &
           (np.abs(S["qrate"]) < Q_MAX) & (S["V"] > V_MIN) &
           (cl > 0.0) & (cl < 1.3))
    n_air = np.sum(np.isfinite(S["V"]) & (S["h"] > H_MIN))
    print(f"window t=[{T_WIN[0]:.0f},{T_WIN[1]:.0f}] s (wind estimator converged); "
          f"{sel.sum()} of {n_air} airborne samples pass the gates")
    print(f"rho (log) = {S['rho']:.3f} kg/m^3;  V = [{S['V'][sel].min():.1f}, "
          f"{S['V'][sel].max():.1f}] m/s;  throttle = [{S['delta'][sel].min():.2f}, "
          f"{S['delta'][sel].max():.2f}];  CL = [{cl[sel].min():.2f}, {cl[sel].max():.2f}]")
    print(f"deck polar: k_ind = {k_ind:.4f} (CD ~ {cd_min_deck:.4f} + {k_ind:.3f} CL^2)")

    sel &= np.isfinite(S["P_elec"]) & (S["P_elec"] > 30.0)
    y = (MASS * (S["V"][sel] * S["Vdot"][sel] + G * S["hdot"][sel]) / S["V"][sel]
         + qbarS[sel] * k_ind * cl[sel] ** 2)
    d, V, P = S["delta"][sel], S["V"][sel], S["P_elec"][sel]

    # Primary: power-proxy energy fit, linear in [eta_p, CD0].
    Xp = np.column_stack([P / V, -qbarS[sel]])
    pp, sep, corrp, rp, keep = fit(Xp, y, clip=3.0)
    eta_p, cd0_fit = pp
    print(f"\nconfig (aircraft.yaml):  T_static = {CFG_TSTATIC:.1f} N   "
          f"k_v = {CFG_KV:.3f}   parasite_CD0 = {CFG_CD0:.3f}")
    print(f"power-proxy fit:  eta_p = {eta_p:.3f} +- {sep[0]:.3f}   "
          f"CD0 = {cd0_fit:.4f} +- {sep[1]:.4f}   "
          f"(corr {corrp[0,1]:+.2f}, rms {np.std(rp):.2f} N, "
          f"{keep.sum()}/{len(y)} kept after 3-sigma clip)")
    print(f"  -> total zero-lift CD = deck {cd_min_deck:.4f} + parasite: fitted "
          f"parasite ~ {cd0_fit:.4f} vs config {CFG_CD0:.3f}")

    # Secondary: thrust-map characterization across the WHOLE flight using
    # T_est = eta_p*P/V (power is wind-immune; airspeed pre-45 s uses ground
    # speed, +-2 m/s wind caveat).
    Sf = load_states(sys.argv[1] if len(sys.argv) > 1 else
                     os.path.join(ROOT, "data", "Takeoff_2_18_07_24.ulg"),
                     t_win=(12.0, 133.0))
    mf = (np.isfinite(Sf["V"]) & np.isfinite(Sf["delta"]) &
          np.isfinite(Sf["P_elec"]) & (Sf["h"] > H_MIN) & (Sf["V"] > V_MIN) &
          (Sf["P_elec"] > 30.0))
    dT, VT, PT = Sf["delta"][mf], Sf["V"][mf], Sf["P_elec"][mf]
    T_est = eta_p * PT / VT
    # Fit both map forms to the (delta, V, T_est) cloud.
    lin, *_ = np.linalg.lstsq(np.column_stack([dT, -dT * VT]), T_est, rcond=None)
    quad, *_ = np.linalg.lstsq(np.column_stack([dT * dT, -dT * dT * VT]), T_est,
                               rcond=None)
    rms_lin = np.std(T_est - dT * (lin[0] - lin[1] * VT))
    rms_quad = np.std(T_est - dT * dT * (quad[0] - quad[1] * VT))
    print(f"\nthrust map vs T_est = eta_p*P/V over the whole flight "
          f"(delta in [{dT.min():.2f}, {dT.max():.2f}]):")
    print(f"  linear    T = d (Ts - kv V):  Ts = {lin[0]:6.1f} N  kv = {lin[1]:6.3f}"
          f"   rms {rms_lin:.2f} N")
    print(f"  quadratic T = d^2 (Ts - kv V):  Ts = {quad[0]:6.1f} N  kv = {quad[1]:6.3f}"
          f"   rms {rms_quad:.2f} N")
    print(f"  (config: linear with Ts = {CFG_TSTATIC:.0f}, kv = {CFG_KV}; the "
          f"better-fitting FORM matters more than the numbers here)")

    # alpha diagnostic (why the drag model is CL-, not alpha-, parameterized).
    alpha = S["theta"] - S["gamma"]
    print(f"alpha = theta - gamma_air over the window: "
          f"mean {np.degrees(np.nanmean(alpha[sel])):.1f} deg, "
          f"std {np.degrees(np.nanstd(alpha[sel])):.1f} deg (no vane -> "
          f"scatter ~ signal; unusable as a regressor)")

    # ---- figure ---------------------------------------------------------------
    fig, ax = plt.subplots(2, 2, figsize=(13, 8.5))
    fig.suptitle("Flight-log validation of aircraft.yaml placeholders "
                 "(energy-rate regression, wind-converged window)",
                 fontsize=13, weight="bold")

    a = ax[0, 0]
    a.plot(S["t"], S["V"], color="tab:blue", lw=1.0, label="V_air")
    a.plot(S["t"][sel], S["V"][sel], ".", color="tab:blue", ms=2.5, label="selected")
    a.plot(S["t"], 20 * S["delta"], color="tab:gray", lw=1.0, alpha=0.8,
           label="throttle x 20")
    a.plot(S["t"], np.degrees(np.abs(S["phi"])) / 3, color="tab:orange", lw=0.8,
           alpha=0.6, label="|bank| / 3 [deg]")
    a.set(xlabel="t [s]", ylabel="V [m/s] (and scaled inputs)",
          title=f"Coverage: {sel.sum()} quasi-steady wings-level samples")
    a.legend(fontsize=8)

    a = ax[0, 1]
    yhat = Xp @ pp
    a.plot(yhat[keep], y[keep], ".", ms=3, color="tab:blue", alpha=0.5,
           label="kept")
    a.plot(yhat[~keep], y[~keep], "x", ms=4, color="tab:red", alpha=0.6,
           label="clipped")
    lo = min(np.percentile(yhat[keep], 1), np.percentile(y[keep], 1)) - 2
    hi = max(np.percentile(yhat[keep], 99), np.percentile(y[keep], 99)) + 2
    a.plot([lo, hi], [lo, hi], "k--", lw=1)
    a.set_xlim(lo, hi); a.set_ylim(lo, hi)
    a.legend(fontsize=8)
    a.set(xlabel=r"model: $\eta_p P/V - \bar{q}S\,CD_0$ [N]",
          ylabel="measured energy rate [N]",
          title=f"Power-proxy energy fit (rms {np.std(rp):.2f} N)")

    a = ax[1, 0]
    a.plot(np.degrees(alpha[sel]), cl[sel], ".", ms=3, color="tab:blue", alpha=0.4,
           label="flight (accel-based CL)")
    aa = np.linspace(-np.radians(6), np.radians(11), 80)
    a.plot(np.degrees(aa), np.interp(aa, a_deck, cl_deck), "k-", lw=1.5,
           label="VSPAERO deck")
    a.set(xlabel=r"$\alpha=\theta-\gamma_{air}$ [deg] (NO VANE — diagnostic only)",
          ylabel=r"$C_L$",
          title="Why alpha is not used: scatter ~ signal (no vane)")
    a.legend(fontsize=8)

    a = ax[1, 1]
    sc = a.scatter(dT, T_est, c=VT, s=6, cmap="viridis", alpha=0.6)
    fig.colorbar(sc, ax=a, label="V [m/s]")
    dd = np.linspace(0.0, max(0.75, dT.max()), 60)
    Vref = float(np.median(VT))
    a.plot(dd, dd * (CFG_TSTATIC - CFG_KV * Vref), color="tab:gray", ls="-",
           lw=1.5, label=f"config: d(50-0.05V), V={Vref:.0f}")
    a.plot(dd, dd * (lin[0] - lin[1] * Vref), color="tab:blue", ls="--", lw=1.2,
           label=f"fit linear (rms {rms_lin:.1f} N)")
    a.plot(dd, dd * dd * (quad[0] - quad[1] * Vref), color="tab:red", ls="-",
           lw=1.2, label=f"fit quadratic (rms {rms_quad:.1f} N)")
    a.set(xlabel="throttle delta (actuator_motors)", ylabel=r"$T_{est}=\eta_p P/V$ [N]",
          title="Thrust map: config linear form vs flight (power-implied)")
    a.legend(fontsize=7)

    for a in ax.flat:
        a.grid(alpha=0.3)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    fig.savefig(out, dpi=120)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
