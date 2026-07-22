#!/usr/bin/env python3
"""Recorded-state landing figure from a PX4 ULog, laid out like lon_alt_noise_lpf.

Reads the real water-landing in data/Takeoff_2_18_07_24.ulg and plots the
longitudinal states the CBF-QP filter consumes (altitude, descent rate,
airspeed, thrust) plus the signals needed to derive them (wind, vertical
acceleration / impact, flight-path & pitch angle, the filtered altitude the
CBF would actually see).

Time zero is the first ``vehicle_local_position`` sample (when the EKF/flight
data begins); the landing window the operator identified is t in [100, 142] s,
with water impact at the vertical-acceleration spike (~141 s).

This module is the front-half of the pipeline: ``load_landing_states`` returns a
clean, uniformly-resampled state table (dict of numpy arrays on a common time
grid) so the *next* step -- running the CBF barriers on the recorded states --
can consume the same arrays without re-parsing the log:

    from plot_log_states import load_landing_states
    S = load_landing_states("data/Takeoff_2_18_07_24.ulg")
    # S["t"], S["h"], S["V_air"], S["gamma"], S["theta"], S["q"], S["T_N"], ...

Signals (per the operator):
  altitude NED   : vehicle_local_position/z         -> h = -z
  descent rate   : vehicle_local_position/vz        (NED, down +)
  ground vel     : d/dt of vehicle_local_position x,y,z
  airspeed       : ground velocity - wind estimate  (no sensor airspeed in log)
  wind           : estimator_wind/windspeed_{north,east}
  norm thrust    : actuator_motors/control[0]       in [0,1]  (T_N = T_norm*Tmax)
  water impact   : vehicle_acceleration/xyz[2]      (sharp spike at touchdown)

Usage:
  scripts/plot_log_states.py [data/Takeoff_2_18_07_24.ulg] [out.png] [t0] [t1]
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import pyulog
except ImportError:
    sys.exit("pyulog is required: pip install pyulog")

# --- limits / constants (match data/lon_scenario.yaml; the CBF uses these) ----
DEFAULTS = dict(
    t0=100.0,            # landing-window start [s] (since first local_position sample)
    t1=142.0,            # landing-window end   [s] (just past water impact)
    dt_grid=0.02,        # uniform resample step [s] (50 Hz)
    h_lpf_tau=0.15,      # altitude low-pass time constant [s] (scenario h_lpf_tau)
    Tmax_N=50.0,         # static thrust [N] -> T_N = control[0] * Tmax_N
    V_app=18.0,          # nominal approach airspeed [m/s]
    Vmin=13.5,           # stall-margin floor [m/s]
    V_max=25.0,          # never-exceed airspeed [m/s]
    v_safe=0.1,          # hull-safe touchdown sink rate [m/s]
)


# -----------------------------------------------------------------------------
# helpers
# -----------------------------------------------------------------------------
def _lpf1(x, dt, tau):
    """Causal first-order low-pass, matching the CBF's altitude prefilter."""
    if tau <= 0.0:
        return np.asarray(x, float)
    a = dt / (tau + dt)
    y = np.empty_like(x, dtype=float)
    y[0] = x[0]
    for k in range(1, len(x)):
        y[k] = y[k - 1] + a * (x[k] - y[k - 1])
    return y


def _smooth(x, n):
    """Centered moving average (edge effects are cropped by the pad window)."""
    if n <= 1:
        return x
    k = np.ones(n) / n
    return np.convolve(x, k, mode="same")


class _Topics:
    """Lazy ULog topic accessor with time relative to the first local_position."""

    def __init__(self, ulog):
        self.ulog = ulog
        lp = self._raw("vehicle_local_position")
        self.t0_us = lp.data["timestamp"][0]

    def _raw(self, name, multi_id=0):
        for d in self.ulog.data_list:
            if d.name == name and d.multi_id == multi_id:
                return d
        return None

    def has(self, name):
        return self._raw(name) is not None

    def series(self, name, field):
        """(t_rel [s], values) for one field, sorted & de-duplicated in time."""
        d = self._raw(name)
        if d is None or field not in d.data:
            return None, None
        t = (d.data["timestamp"] - self.t0_us) / 1e6
        v = np.asarray(d.data[field], float)
        order = np.argsort(t, kind="stable")
        t, v = t[order], v[order]
        keep = np.concatenate(([True], np.diff(t) > 0))
        return t[keep], v[keep]


def _resample(tp, t, v, grid):
    """Interpolate (t, v) onto grid; NaN outside the source span."""
    if t is None:
        return np.full_like(grid, np.nan)
    out = np.interp(grid, t, v)
    out[(grid < t[0]) | (grid > t[-1])] = np.nan
    return out


# -----------------------------------------------------------------------------
# state loader (front-half of the pipeline; reused by the CBF step)
# -----------------------------------------------------------------------------
def load_landing_states(ulg_path, t0=None, t1=None, dt_grid=None,
                        h_lpf_tau=None, Tmax_N=None):
    """Parse the ULog and return CBF-ready landing states on a uniform grid.

    Returns a dict of numpy arrays all sharing ``t`` (the uniform time grid over
    [t0, t1]) plus native-rate acceleration (for the sharp impact spike) and
    scalar metadata. NaNs mark grid points outside a source topic's span.
    """
    t0 = DEFAULTS["t0"] if t0 is None else t0
    t1 = DEFAULTS["t1"] if t1 is None else t1
    dt_grid = DEFAULTS["dt_grid"] if dt_grid is None else dt_grid
    h_lpf_tau = DEFAULTS["h_lpf_tau"] if h_lpf_tau is None else h_lpf_tau
    Tmax_N = DEFAULTS["Tmax_N"] if Tmax_N is None else Tmax_N

    ulog = pyulog.ULog(ulg_path)
    T = _Topics(ulog)
    grid = np.arange(t0, t1 + 0.5 * dt_grid, dt_grid)

    # --- position & EKF velocity (pad the window so derivatives have no edge bias)
    pad = 3.0
    tx, x = T.series("vehicle_local_position", "x")
    _, y = T.series("vehicle_local_position", "y")
    _, z = T.series("vehicle_local_position", "z")
    _, vz_ekf = T.series("vehicle_local_position", "vz")
    sel = (tx >= t0 - pad) & (tx <= t1 + pad)
    txp = tx[sel]
    # ground velocity from position derivative (operator's recipe), lightly
    # smoothed to tame the 10 Hz differentiation noise, then put on the grid.
    win = max(1, int(round(0.4 / np.median(np.diff(txp)))))  # ~0.4 s window
    vN = _resample(T, txp, _smooth(np.gradient(x[sel], txp), win), grid)
    vE = _resample(T, txp, _smooth(np.gradient(y[sel], txp), win), grid)
    vD = _resample(T, txp, _smooth(np.gradient(z[sel], txp), win), grid)

    h = _resample(T, tx, -z, grid)            # altitude AGL [m] (= -z)
    z_g = _resample(T, tx, z, grid)           # NED z [m]
    vz = _resample(T, tx, vz_ekf, grid)       # descent rate [m/s] (down +, EKF field)

    # --- wind estimate (slow; linear interp is fine) --------------------------
    twn, wn = T.series("estimator_wind", "windspeed_north")
    twe, we = T.series("estimator_wind", "windspeed_east")
    wind_n = _resample(T, twn, wn, grid)
    wind_e = _resample(T, twe, we, grid)

    # --- airspeed = ground velocity - wind (no valid sensor airspeed in log) ---
    V_gnd = np.hypot(vN, vE)                       # horizontal ground speed
    V_gnd3d = np.sqrt(vN**2 + vE**2 + vD**2)
    V_air = np.sqrt((vN - wind_n)**2 + (vE - wind_e)**2 + vD**2)
    wind_mag = np.hypot(wind_n, wind_e)

    # --- flight-path angle & pitch (CBF states gamma, theta) ------------------
    gamma = np.arctan2(-vD, V_gnd)                 # +climb / -descent [rad]
    tq, q0 = T.series("vehicle_attitude", "q[0]")
    _, q1 = T.series("vehicle_attitude", "q[1]")
    _, q2 = T.series("vehicle_attitude", "q[2]")
    _, q3 = T.series("vehicle_attitude", "q[3]")
    if tq is not None:
        pitch = np.arcsin(np.clip(2.0 * (q0 * q2 - q3 * q1), -1.0, 1.0))
        theta = _resample(T, tq, pitch, grid)
    else:
        theta = np.full_like(grid, np.nan)
    # pitch rate q (CBF state); body-y angular velocity if available
    tw, wy = T.series("vehicle_angular_velocity", "xyz[1]")
    q_pitch = _resample(T, tw, wy, grid) if tw is not None else np.full_like(grid, np.nan)

    # --- normalized thrust -> Newtons -----------------------------------------
    tt, c0 = T.series("actuator_motors", "control[0]")
    T_norm = _resample(T, tt, c0, grid)
    T_N = T_norm * Tmax_N

    # --- commanded pitch: pilot stick (flight was fully manual) ---------------
    # manual_control_setpoint/pitch is the normalized [-1,1] pitch command; in
    # manual mode it passes through to vehicle_torque_setpoint/xyz[1] (sign
    # flipped). This is the recorded analog of the sim's commanded pitch.
    tmc, pcmd = T.series("manual_control_setpoint", "pitch")
    pitch_cmd = _resample(T, tmc, pcmd, grid)
    # elevon servo positions kept for reference / the CBF step (normalized [-1,1])
    ts2, e0 = T.series("actuator_servos", "control[0]")
    _, e1 = T.series("actuator_servos", "control[1]")
    elevon_l = _resample(T, ts2, e0, grid)
    elevon_r = _resample(T, ts2, e1, grid)
    elevon_pitch = 0.5 * (elevon_l + elevon_r)   # symmetric (elevator-equiv) part

    # --- vertical acceleration: keep native rate so the impact spike survives --
    tacc, az = T.series("vehicle_acceleration", "xyz[2]")
    am = (tacc >= t0) & (tacc <= t1)
    t_acc, az_w = tacc[am], az[am]
    j = int(np.argmax(np.abs(az_w)))               # impact = largest |specific force|
    t_impact, az_impact = t_acc[j], az_w[j]

    # --- altitude the CBF would see: same low-pass prefilter as the sim --------
    h_lpf = _lpf1(h, dt_grid, h_lpf_tau)

    return dict(
        t=grid, h=h, z=z_g, vz=vz,
        vN=vN, vE=vE, vD=vD,
        V_gnd=V_gnd, V_gnd3d=V_gnd3d, V_air=V_air,
        wind_n=wind_n, wind_e=wind_e, wind_mag=wind_mag,
        gamma=gamma, gamma_deg=np.rad2deg(gamma),
        theta=theta, theta_deg=np.rad2deg(theta),
        q=q_pitch,
        T_norm=T_norm, T_N=T_N,
        pitch_cmd=pitch_cmd,
        elevon_l=elevon_l, elevon_r=elevon_r, elevon_pitch=elevon_pitch,
        h_lpf=h_lpf,
        t_acc=t_acc, az=az_w, t_impact=t_impact, az_impact=az_impact,
        meta=dict(ulg=os.path.basename(ulg_path), t0=t0, t1=t1, dt_grid=dt_grid,
                  h_lpf_tau=h_lpf_tau, Tmax_N=Tmax_N),
    )


# -----------------------------------------------------------------------------
# figure
# -----------------------------------------------------------------------------
def plot_states(S, out_path, lim=None):
    """Recorded-state figure: the reference figure's directly-recorded quantities.

    Panels mirror lon_alt_noise_lpf one-for-one for the quantities the log
    supplies on its own (altitude, descent, airspeed, thrust, commanded pitch,
    altitude sensor chain). The CBF-output panels of the reference (nominal-vs-
    filtered elevator / thrust-jerk, the barrier stack) are intentionally left
    to the CBF step. Water impact is marked on every panel.
    """
    lim = lim or DEFAULTS
    t = S["t"]
    ti = S["t_impact"]
    m = S["meta"]
    fig, ax = plt.subplots(2, 3, figsize=(17, 8))
    fig.suptitle(
        f"Recorded landing states  —  {m['ulg']}  (t = {m['t0']:.0f}–{m['t1']:.0f} s, "
        f"impact {ti:.1f} s)",
        fontsize=14, weight="bold")

    def mark_impact(a, label=False):
        a.axvline(ti, color="tab:red", ls="--", lw=1.0,
                  label=("water impact" if label else None))

    # (1) Altitude h = -z -------------------------------------------------------
    ax[0, 0].plot(t, S["h"], color="tab:blue")
    ax[0, 0].axhline(0, color="k", lw=0.8, ls=":", label="water (h=0)")
    mark_impact(ax[0, 0], label=True)
    ax[0, 0].set_ylabel("altitude h [m]"); ax[0, 0].set_title("Altitude (h = -z, NED)")
    ax[0, 0].legend(fontsize=8)

    # (2) Descent rate (vz, down +) --------------------------------------------
    ax[0, 1].plot(t, S["vz"], color="tab:red")
    ax[0, 1].axhline(0, color="k", lw=0.8, ls=":")
    ax[0, 1].axhline(lim["v_safe"], color="tab:green", ls=":", lw=1,
                     label=f"v_safe={lim['v_safe']} m/s")
    mark_impact(ax[0, 1])
    ax[0, 1].set_ylabel("descent rate [m/s] (down +)")
    ax[0, 1].set_title("Descent rate (vehicle_local_position/vz)")
    ax[0, 1].legend(fontsize=8)

    # (3) Airspeed vs velocity limits ------------------------------------------
    ax[0, 2].plot(t, S["V_air"], color="tab:green", label="airspeed V (gnd vel - wind)")
    ax[0, 2].plot(t, S["V_gnd"], color="tab:gray", alpha=0.6, lw=0.9,
                  label="ground speed (horiz)")
    ax[0, 2].axhline(lim["Vmin"], color="k", ls=":", lw=0.9, label=f"Vmin={lim['Vmin']}")
    ax[0, 2].axhline(lim["V_max"], color="k", ls="--", lw=0.9, label=f"Vmax={lim['V_max']}")
    mark_impact(ax[0, 2])
    ax[0, 2].set_ylabel("V [m/s]"); ax[0, 2].set_title("Airspeed vs velocity limits")
    ax[0, 2].legend(fontsize=8)

    # (4) Thrust vs thrust limits [N] ------------------------------------------
    ax[1, 0].plot(t, S["T_N"], color="tab:orange")
    ax[1, 0].axhline(0, color="k", lw=0.9, ls=":", label="Tmin=0")
    ax[1, 0].axhline(lim["Tmax_N"], color="k", lw=0.9, ls="--", label=f"Tmax={lim['Tmax_N']:.0f}")
    mark_impact(ax[1, 0])
    ax[1, 0].set_ylabel("thrust T [N]"); ax[1, 0].set_xlabel("t [s]")
    ax[1, 0].set_title(f"Thrust vs thrust limits (control[0]·{lim['Tmax_N']:.0f} N)")
    ax[1, 0].legend(fontsize=8)

    # (5) Commanded pitch (manual stick) ---------------------------------------
    ax[1, 1].plot(t, S["pitch_cmd"], color="tab:purple", label="commanded pitch (stick)")
    ax[1, 1].axhline(0, color="k", lw=0.8, ls=":")
    mark_impact(ax[1, 1])
    ax[1, 1].set_ylabel("pitch command [-1,1]"); ax[1, 1].set_xlabel("t [s]")
    ax[1, 1].set_title("Commanded pitch (manual_control_setpoint/pitch)")
    ax[1, 1].set_ylim(-1.05, 1.05)
    ax[1, 1].legend(fontsize=8)

    # (6) Altitude sensor chain: EKF h vs the LPF the CBF will see -------------
    ax[1, 2].plot(t, S["h"], color="tab:gray", lw=0.9, label="EKF h (CBF measurement)")
    ax[1, 2].plot(t, S["h_lpf"], color="tab:green",
                  label=fr"LPF h ($\tau$={m['h_lpf_tau']} s, CBF input)")
    ax[1, 2].axhline(0, color="k", lw=0.8, ls=":")
    mark_impact(ax[1, 2])
    ax[1, 2].set_ylabel("altitude [m]"); ax[1, 2].set_xlabel("t [s]")
    ax[1, 2].set_title("Altitude sensor chain (CBF altitude prefilter)")
    ax[1, 2].legend(fontsize=8)

    for a in ax.flat:
        a.grid(alpha=0.3)
        a.set_xlim(m["t0"], m["t1"])
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_path, dpi=110)
    print("wrote", out_path)


# -----------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("ulg", nargs="?", default="data/Takeoff_2_18_07_24.ulg")
    p.add_argument("out", nargs="?", default="figures/log_landing_states.png")
    p.add_argument("t0", nargs="?", type=float, default=DEFAULTS["t0"])
    p.add_argument("t1", nargs="?", type=float, default=DEFAULTS["t1"])
    args = p.parse_args()

    S = load_landing_states(args.ulg, t0=args.t0, t1=args.t1)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    plot_states(S, args.out)

    # console summary
    def rng(k):
        v = S[k]; v = v[~np.isnan(v)]
        return f"[{v.min():.2f}, {v.max():.2f}]" if v.size else "n/a"
    print(f"  window t=[{args.t0}, {args.t1}] s, {len(S['t'])} grid pts @ "
          f"{S['meta']['dt_grid']*1e3:.0f} ms")
    print(f"  h        {rng('h')} m")
    print(f"  descent  {rng('vz')} m/s")
    print(f"  V_air    {rng('V_air')} m/s   V_gnd {rng('V_gnd')} m/s")
    print(f"  thrust   {rng('T_norm')} (norm)")
    print(f"  impact   t={S['t_impact']:.2f} s  a_z={S['az_impact']:.1f} m/s^2 "
          f"({abs(S['az_impact'])/9.81:.1f} g)")


if __name__ == "__main__":
    main()
