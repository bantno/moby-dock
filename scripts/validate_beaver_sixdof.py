#!/usr/bin/env python3
"""Full-envelope 6-DOF validation of the Beaver plant.

The wings-level oracles (FDC check case, trim curves, trim-point Jacobians)
all sit at p = q = r = 0, where the quadratic gyroscopic/Coriolis terms and
the large-attitude kinematics contribute nothing. This script exercises
exactly that remaining regime, against the same independent Python
implementation used by validate_beaver_modes.py (FDC state coordinates,
freshly-typed coefficients):

 1. RANDOM-STATE SWEEP: N states spanning the whole envelope -- V 28-58 m/s,
    alpha/beta to +/-15 deg, body rates to +/-1.2 rad/s, bank to +/-60 deg,
    full control deflections, flaps 0/20/35, three altitudes/RPMs. Every row
    of the C++ xdot must match the independent implementation pointwise.
 2. STEADY COORDINATED TURNS: level-turn equilibria (phi up to +/-45 deg,
    nonzero p,q,r) solved on the INDEPENDENT implementation, checked for the
    physical invariants (load factor ~ 1/cos(phi), turn rate g tan(phi)/V),
    then handed to the C++ plant, which must agree they are equilibria --
    the classic test of the gyroscopic terms.
 3. OPEN-LOOP DOUBLET: a 14 s elevator-doublet + aileron + rudder maneuver
    integrated by BOTH implementations from the same state/controls; the
    trajectories must overlay.

Runs the C++ side via `<build>/beaver_validation --xdot / --doublet`.
Writes results/beaver_sixdof_validation.png. Exits nonzero on failure.

Usage: validate_beaver_sixdof.py [build_dir] [out_dir]
       (defaults: build, results)
"""
import os
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import validate_beaver_modes as vbm  # noqa: E402

build = sys.argv[1] if len(sys.argv) > 1 else "build"
outdir = sys.argv[2] if len(sys.argv) > 2 else "results"
kDeg = np.pi / 180.0
PZ_IDLE, PZ_MAX = 5.0, 26.0  # the plant's throttle->pz map


def icao_rho(h):
    T = 288.15 - 0.0065 * h
    p = 101325.0 * (T / 288.15) ** (9.80665 / (0.0065 * 287.05))
    return p / (287.05 * T)


def grav(h):
    return 9.80665 * (6371020.0 / (6371020.0 + h)) ** 2


def full_xdot(s9, ctl, prm):
    """Independent-implementation xdot with the full kinematics.

    s9 = [V, alpha, beta, p, q, r, psi, theta, phi] -> returns
    [Vdot, alphadot, betadot, pdot, qdot, rdot, psidot, thetadot, phidot,
     Hdot, yedot].
    """
    V, al, be, p, q, r, psi, th, ph = s9
    d8 = vbm.xdot(np.array([V, al, be, p, q, r, th, ph]), ctl, prm)
    psidot = (q * np.sin(ph) + r * np.cos(ph)) / np.cos(th)
    u = V * np.cos(al) * np.cos(be)
    v = V * np.sin(be)
    w = V * np.sin(al) * np.cos(be)
    st, ct = np.sin(th), np.cos(th)
    sp, cp = np.sin(ph), np.cos(ph)
    sy, cy = np.sin(psi), np.cos(psi)
    Hdot = u * st - v * sp * ct - w * cp * ct
    yedot = (u * ct * sy + v * (sp * st * sy + cp * cy)
             + w * (cp * st * sy - sp * cy))
    return np.concatenate([d8[:6], [psidot], d8[6:8], [Hdot, yedot]])


# ============ 1. random-state sweep ==========================================
rng = np.random.default_rng(2026)
N = 250
states = []
for _ in range(N):
    states.append(dict(
        V=rng.uniform(28, 58), al=rng.uniform(-12, 16) * kDeg,
        be=rng.uniform(-15, 15) * kDeg,
        p=rng.uniform(-1.2, 1.2), q=rng.uniform(-0.8, 0.8),
        r=rng.uniform(-0.8, 0.8),
        psi=rng.uniform(-np.pi, np.pi), th=rng.uniform(-25, 35) * kDeg,
        ph=rng.uniform(-60, 60) * kDeg, h=rng.uniform(0, 100),
        de=rng.uniform(-25, 20) * kDeg, da=rng.uniform(-20, 20) * kDeg,
        dr=rng.uniform(-25, 25) * kDeg, dT=rng.uniform(0, 1),
        flap=rng.choice([0.0, 20.0, 35.0]) * kDeg,
        h_ref=rng.choice([0.0, 609.6, 1828.8]),
        n_rpm=rng.choice([1650.0, 1800.0, 2300.0])))

# ============ 2. steady coordinated level turns ==============================
def turn_trim(V, phi, n_rpm=1800.0, h_ref=0.0, flap=0.0):
    """Level coordinated-turn equilibrium on the INDEPENDENT implementation.
    Unknowns z = [alpha, beta, theta, de, da, dr, pz]; Omega = g tan(phi)/V;
    residuals = 6 accelerations + Hdot. Complex-step Newton."""
    rho, g = icao_rho(h_ref), grav(h_ref)
    Om = g * np.tan(phi) / V

    def rates(th):
        return (-Om * np.sin(th), Om * np.sin(phi) * np.cos(th),
                Om * np.cos(phi) * np.cos(th))

    def resid(z):
        al, be, th, de, da, dr, pz = z
        p, q, r = rates(th)
        prm = dict(pz=pz, n_rpm=n_rpm, rho=rho, g=g, flap=flap)
        s9 = np.array([V + 0j, al, be, p, q, r, 0.0, th, phi + 0j])
        d = full_xdot(s9, [de, da, dr], prm)
        return np.array([d[0], d[1], d[2], d[3], d[4], d[5], d[9]])

    z = np.array([0.06, 0.0, 0.06, -0.06, 0.0, 0.0, 18.0], dtype=complex)
    for _ in range(60):
        R = resid(z).real
        if np.max(np.abs(R)) < 1e-11:
            break
        J = np.zeros((7, 7))
        for j in range(7):
            zp = z.copy()
            zp[j] += 1e-30j
            J[:, j] = np.imag(resid(zp)) / 1e-30
        z = z - np.linalg.solve(J, R)
    al, be, th, de, da, dr, pz = z.real
    p, q, r = rates(th)
    return dict(V=V, al=al, be=be, th=th, ph=phi, p=p, q=q, r=r, psi=0.0,
                de=de, da=da, dr=dr, pz=pz, h=50.0, flap=flap,
                h_ref=h_ref, n_rpm=n_rpm, Om=Om,
                res=np.max(np.abs(resid(z).real)))


turns = [turn_trim(40.0, s * b * kDeg) for b in (15.0, 30.0, 45.0)
         for s in (+1, -1)]
print("=== steady coordinated level turns (independent implementation) ===")
turn_rows = []
for t in turns:
    rho = icao_rho(t["h_ref"])
    prm = dict(pz=t["pz"], n_rpm=t["n_rpm"], rho=rho, g=grav(t["h_ref"]),
               flap=t["flap"])
    P = vbm.engine_power(t["pz"], t["n_rpm"], rho)
    dpt = 0.08696 + 191.18 * 2 * P / (rho * t["V"] ** 3)
    cz = vbm.coeffs(t["al"], t["be"], t["p"], t["q"], t["r"], t["V"],
                    t["de"], t["da"], t["dr"], t["flap"], dpt)[2]
    nz = -cz * 0.5 * rho * t["V"] ** 2 * vbm.S / (vbm.MASS * 9.80665)
    n_exp = 1.0 / np.cos(t["ph"])
    R = t["V"] ** 2 / (grav(0) * np.tan(t["ph"])) if t["ph"] != 0 else np.inf
    print(f"  phi={np.degrees(t['ph']):+5.1f} deg: alpha={np.degrees(t['al']):5.2f}"
          f"  de={np.degrees(t['de']):+6.2f}  pz={t['pz']:5.2f}"
          f"  n_z={nz:5.3f} (1/cos_phi={n_exp:5.3f})"
          f"  R={R:+7.1f} m  psidot={np.degrees(t['Om']):+6.2f} deg/s"
          f"  res={t['res']:.1e}")
    turn_rows.append((t, nz, n_exp))
    assert abs(nz - n_exp) < 0.06 * n_exp, "load factor off"

# append turn states to the C++ sweep with dT from pz
for t in turns:
    states.append(dict(
        V=t["V"], al=t["al"], be=t["be"], p=t["p"], q=t["q"], r=t["r"],
        psi=0.0, th=t["th"], ph=t["ph"], h=t["h"], de=t["de"], da=t["da"],
        dr=t["dr"], dT=(t["pz"] - PZ_IDLE) / (PZ_MAX - PZ_IDLE),
        flap=t["flap"], h_ref=t["h_ref"], n_rpm=t["n_rpm"]))

# ============ write states, run the C++ side ================================
os.makedirs(outdir, exist_ok=True)
in_csv = os.path.join(outdir, "beaver_sixdof_states.csv")
out_csv = os.path.join(outdir, "beaver_sixdof_xdot_cpp.csv")
with open(in_csv, "w") as f:
    f.write("h_ref,n_rpm,flap,u,v,w,p,q,r,phi,theta,psi,h,y,de,da,dr,dT\n")
    for s in states:
        u = s["V"] * np.cos(s["al"]) * np.cos(s["be"])
        v = s["V"] * np.sin(s["be"])
        w = s["V"] * np.sin(s["al"]) * np.cos(s["be"])
        f.write(",".join(f"{x:.17g}" for x in [
            s["h_ref"], s["n_rpm"], s["flap"], u, v, w, s["p"], s["q"],
            s["r"], s["ph"], s["th"], s["psi"], s["h"], 0.0,
            s["de"], s["da"], s["dr"], s["dT"]]) + "\n")

exe = os.path.join(build, "beaver_validation")
subprocess.run([exe, "--xdot", in_csv, out_csv], check=True)
subprocess.run([exe, "--doublet", os.path.join(outdir, "beaver_doublet_cpp.csv")],
               check=True)

# ============ compare the sweep pointwise ===================================
cpp = np.loadtxt(out_csv, delimiter=",")
ROWS = ["Vdot", "alphadot", "betadot", "pdot", "qdot", "rdot",
        "psidot", "thetadot", "phidot", "Hdot", "yedot"]
errs = np.zeros((len(states), 11))
for i, s in enumerate(states):
    prm = dict(pz=PZ_IDLE + s["dT"] * (PZ_MAX - PZ_IDLE), n_rpm=s["n_rpm"],
               rho=icao_rho(s["h_ref"]), g=grav(s["h_ref"]), flap=s["flap"])
    s9 = np.array([s["V"], s["al"], s["be"], s["p"], s["q"], s["r"],
                   s["psi"], s["th"], s["ph"]], dtype=complex)
    py = full_xdot(s9, [s["de"], s["da"], s["dr"]], prm).real
    # C++ order: udot vdot wdot pdot qdot rdot phidot thetadot psidot Hdot ydot
    u = s["V"] * np.cos(s["al"]) * np.cos(s["be"])
    v = s["V"] * np.sin(s["be"])
    w = s["V"] * np.sin(s["al"]) * np.cos(s["be"])
    ud, vd, wd = cpp[i, 0], cpp[i, 1], cpp[i, 2]
    Vt = s["V"]
    Vd = (u * ud + v * vd + w * wd) / Vt
    ald = (u * wd - w * ud) / (u * u + w * w)
    bed = (vd * Vt - v * Vd) / (Vt * np.sqrt(u * u + w * w))
    mine = np.array([Vd, ald, bed, cpp[i, 3], cpp[i, 4], cpp[i, 5],
                     cpp[i, 8], cpp[i, 7], cpp[i, 6], cpp[i, 9], cpp[i, 10]])
    errs[i] = np.abs(mine - py) / np.maximum(1.0, np.abs(py))

max_err = errs.max()
print(f"\n=== random-state sweep: {len(states)} states "
      f"({N} random + {len(turns)} turn equilibria) ===")
print(f"  max normalized |cpp - python| over all 11 rows: {max_err:.2e}")
worst = np.unravel_index(errs.argmax(), errs.shape)
print(f"  worst: state {worst[0]}, row {ROWS[worst[1]]}")
sweep_ok = max_err < 1e-8

# turn equilibria: C++ must also see xdot ~ 0 on the dynamic rows
turn_res_cpp = np.abs(cpp[N:, :6]).max()
print(f"  C++ dynamic-row residual at the {len(turns)} turn equilibria: "
      f"{turn_res_cpp:.2e}  (equilibria solved on the independent model)")
turns_ok = turn_res_cpp < 1e-8

# ============ doublet: re-integrate on the independent implementation ========
cppd = np.genfromtxt(os.path.join(outdir, "beaver_doublet_cpp.csv"),
                     delimiter=",", names=True)
dt = cppd["t"][1] - cppd["t"][0]
u0, v0, w0 = cppd["u"][0], cppd["v"][0], cppd["w"][0]
V0 = np.sqrt(u0**2 + v0**2 + w0**2)
s = np.array([V0, np.arctan(w0 / u0), np.arcsin(v0 / V0),
              cppd["p"][0], cppd["q"][0], cppd["r"][0],
              cppd["psi"][0], cppd["theta"][0], cppd["phi"][0]])
prm0 = dict(n_rpm=1800.0, rho=icao_rho(0.0), g=grav(0.0), flap=0.0, pz=0.0)
traj = np.zeros((len(cppd["t"]), 9))
for k in range(len(cppd["t"])):
    traj[k] = s
    prm0["pz"] = PZ_IDLE + cppd["dT"][k] * (PZ_MAX - PZ_IDLE)
    ctl = [cppd["de"][k], cppd["da"][k], cppd["dr"][k]]
    f = lambda ss: full_xdot(np.array(ss, dtype=complex), ctl, prm0).real[:9]
    k1 = f(s); k2 = f(s + 0.5 * dt * k1); k3 = f(s + 0.5 * dt * k2)
    k4 = f(s + dt * k3)
    s = s + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)

Vc = np.sqrt(cppd["u"]**2 + cppd["v"]**2 + cppd["w"]**2)
alc = np.arctan(cppd["w"] / cppd["u"])
dV = np.max(np.abs(Vc - traj[:, 0]))
dang = np.degrees(max(np.max(np.abs(cppd["phi"] - traj[:, 8])),
                      np.max(np.abs(cppd["theta"] - traj[:, 7])),
                      np.max(np.abs(cppd["psi"] - traj[:, 6]))))
drate = np.max(np.abs(np.column_stack([cppd["p"], cppd["q"], cppd["r"]])
                      - traj[:, 3:6]))
print(f"\n=== 14 s open-loop doublet, C++ vs independent integration ===")
print(f"  max |dV| = {dV:.2e} m/s   max attitude diff = {dang:.2e} deg   "
      f"max rate diff = {drate:.2e} rad/s")
doublet_ok = (dV < 1e-3) and (dang < 1e-2) and (drate < 1e-4)

# ============ figure ========================================================
fig = plt.figure(figsize=(16, 8.6))
fig.suptitle("Beaver 6-DOF full-envelope validation: C++ plant vs independent "
             "implementation", fontsize=13, weight="bold")
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15])

axA = fig.add_subplot(gs[0, 0])
axA.boxplot([np.log10(np.maximum(errs[:, j], 1e-18)) for j in range(11)],
            tick_labels=ROWS, showfliers=True, flierprops=dict(ms=2))
axA.axhline(-8, color="tab:red", ls="--", lw=1)
axA.text(6, -7.6, "pass threshold 1e-8", color="tab:red", fontsize=8,
         ha="center")
axA.set_ylabel(r"$\log_{10}$ normalized $|\Delta \dot x|$")
axA.set_title(f"Random-state sweep: {len(states)} states\n"
              "(rates to 1.2 rad/s, bank to 60$^\\circ$, flaps, full controls)",
              fontsize=9)
axA.tick_params(axis="x", rotation=60, labelsize=7)
axA.grid(alpha=0.3, axis="y")

axB = fig.add_subplot(gs[0, 1])
phis = [np.degrees(t["ph"]) for t, _, _ in turn_rows]
nzs = [nz for _, nz, _ in turn_rows]
nexp = [ne for _, _, ne in turn_rows]
axB.plot(phis, nzs, "x", ms=10, color="tab:blue",
         label="model load factor at turn trim")
ph_line = np.linspace(-50, 50, 100)
axB.plot(ph_line, 1 / np.cos(ph_line * kDeg), color="k", lw=1, alpha=0.6,
         label=r"$1/\cos\phi$ (physics)")
axB.set_xlabel(r"bank $\phi$ [deg]")
axB.set_ylabel(r"$n_z$")
axB.set_title("Steady coordinated turns: load factor\n"
              "(equilibria with nonzero p, q, r)", fontsize=9)
axB.grid(alpha=0.3)
axB.legend(fontsize=8)

axC = fig.add_subplot(gs[0, 2])
axC.semilogy(np.arange(len(turns)),
             np.abs(cpp[N:, :6]).max(axis=1), "s", color="tab:blue")
axC.axhline(1e-8, color="tab:red", ls="--", lw=1)
axC.set_xticks(np.arange(len(turns)),
               [f"{np.degrees(t['ph']):+.0f}$^\\circ$" for t in turns])
axC.set_xlabel("turn bank")
axC.set_ylabel(r"max $|\dot x|$ dynamic rows [SI/s]")
axC.set_title("C++ residual at turn equilibria solved\non the independent "
              "model (gyroscopic terms)", fontsize=9)
axC.grid(alpha=0.3, axis="y")

t = cppd["t"]
axD = fig.add_subplot(gs[1, 0])
axD.plot(t, np.degrees(cppd["phi"]), color="tab:blue", label=r"$\phi$ C++")
axD.plot(t, np.degrees(traj[:, 8]), "--", color="tab:orange", label=r"$\phi$ py")
axD.plot(t, np.degrees(cppd["theta"]), color="tab:green", label=r"$\theta$ C++")
axD.plot(t, np.degrees(traj[:, 7]), "--", color="tab:red", label=r"$\theta$ py")
axD.plot(t, np.degrees(cppd["psi"]), color="tab:purple", label=r"$\psi$ C++")
axD.plot(t, np.degrees(traj[:, 6]), "--", color="tab:brown", label=r"$\psi$ py")
axD.set_xlabel("t [s]"); axD.set_ylabel("[deg]")
axD.set_title("Doublet response: attitude", fontsize=9)
axD.grid(alpha=0.3); axD.legend(fontsize=7, ncol=3)

axE = fig.add_subplot(gs[1, 1])
for nm, cpp_v, py_v, col in [("p", cppd["p"], traj[:, 3], "tab:blue"),
                             ("q", cppd["q"], traj[:, 4], "tab:green"),
                             ("r", cppd["r"], traj[:, 5], "tab:purple")]:
    axE.plot(t, cpp_v, color=col, label=f"{nm} C++")
    axE.plot(t, py_v, "--", color="k", lw=0.8)
axE.set_xlabel("t [s]"); axE.set_ylabel("[rad/s]")
axE.set_title("Doublet response: body rates (dashed = independent py)",
              fontsize=9)
axE.grid(alpha=0.3); axE.legend(fontsize=7)

axF = fig.add_subplot(gs[1, 2])
axF.plot(t, Vc, color="tab:blue", label="V C++")
axF.plot(t, traj[:, 0], "--", color="tab:orange", label="V py")
axFb = axF.twinx()
axFb.plot(t, np.degrees(alc), color="tab:green", label=r"$\alpha$ C++")
axFb.plot(t, np.degrees(traj[:, 1]), "--", color="tab:red", label=r"$\alpha$ py")
axFb.set_ylabel(r"$\alpha$ [deg]", color="tab:green")
axF.set_xlabel("t [s]"); axF.set_ylabel("V [m/s]", color="tab:blue")
axF.set_title(f"Doublet: V, alpha  (max dV {dV:.1e} m/s, "
              f"max att {dang:.1e} deg)", fontsize=9)
axF.grid(alpha=0.3)
axF.legend(fontsize=7, loc="lower left"); axFb.legend(fontsize=7)

fig.tight_layout(rect=(0, 0, 1, 0.95))
out_png = os.path.join(outdir, "beaver_sixdof_validation.png")
fig.savefig(out_png, dpi=140)
print("wrote", out_png)

ok = sweep_ok and turns_ok and doublet_ok
print("\nRESULT:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
