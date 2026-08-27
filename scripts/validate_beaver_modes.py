#!/usr/bin/env python3
"""Independent cross-implementation check of the Beaver 6-DOF linearization.

The C++ plant (beaver_dynamics.hpp) integrates body-axis (u,v,w) equations and
linearizes with autodiff. This script implements the SAME published model
(FDC 1.2 manual, appendix C) a second time, independently:

  * FDC state coordinates [V, alpha, beta, p, q, r, theta, phi] (not body u,v,w),
  * coefficients typed in fresh from the manual's Table C.3/C.4 (not imported
    from the C++ header),
  * complex-step differentiation for the Jacobian (exact to machine precision,
    no finite-difference truncation).

The full-state Jacobians of the two implementations are similar matrices
(smooth invertible change of state coordinates), so their eigenvalues must
agree to numerical precision. Agreement validates the EOM assembly and every
rate/attitude derivative -- the parts a single-point trim oracle cannot see.

Reads the trim/state-matrix dumps written by `beaver_validation` and compares
mode-by-mode. Exits nonzero on mismatch.

Usage: validate_beaver_modes.py [results_prefix]   (default: results/beaver)
"""
import csv
import sys

import numpy as np

# --- Aircraft data, typed from FDC 1.2 manual Table C.1/C.2 ------------------
S, B, CBAR = 23.23, 14.63, 1.5875
MASS = 2288.231
IX, IY, IZ, IXZ = 5368.39, 6928.93, 11158.75, 117.64

# --- Aero polynomial, typed from Table C.3 / C.4 (independent transcription) --
def coeffs(al, be, p, q, r, V, de, da, dr, df, dpt):
    ph = p * B / (2 * V)
    rh = r * B / (2 * V)
    qh = q * CBAR / V
    cx = (-0.03554 + 0.002920 * al + 5.459 * al**2 - 5.162 * al**3
          - 0.6748 * qh + 0.03412 * dr - 0.09447 * df + 1.106 * al * df
          + 0.1161 * dpt + 0.1453 * al * dpt**2)
    cy = (-0.002226 - 0.7678 * be - 0.1240 * ph + 0.3666 * rh
          - 0.02956 * da + 0.1158 * dr + 0.5238 * al * dr)
    cz = (-0.05504 - 5.578 * al + 3.442 * al**3 - 2.988 * qh - 0.3980 * de
          - 15.93 * de * be**2 - 1.377 * df - 1.261 * al * df - 0.1563 * dpt)
    cl = (0.0005910 - 0.06180 * be - 0.5045 * ph + 0.1695 * rh
          - 0.09917 * da + 0.006934 * dr - 0.08269 * al * da
          - 0.01406 * al**2 * dpt)
    cm = (0.09448 - 0.6028 * al - 2.140 * al**2 - 15.56 * qh - 1.921 * de
          + 0.6921 * be**2 - 0.3118 * rh + 0.4072 * df - 0.07895 * dpt)
    cn = (-0.003117 + 0.006719 * be - 0.1585 * ph - 0.1112 * rh
          - 0.003872 * da - 0.08265 * dr + 0.1595 * qh + 0.1373 * be**3
          - 0.003026 * dpt**3)
    return cx, cy, cz, cl, cm, cn


def engine_power(pz, n, rho):
    """FDC eq. 3.15 [kW]; altitude correction on the (408 - 0.0965 n) term."""
    return 0.7355 * (-326.5 + (0.00412 * (pz + 7.4) * (n + 2010.0)
                               + (408.0 - 0.0965 * n) * (1.0 - rho / 1.225)))


def xdot(s, ctl, prm):
    """FDC-coordinate state derivative.

    s   = [V, alpha, beta, p, q, r, theta, phi]
    ctl = [de, da, dr]
    prm = dict(pz, n_rpm, rho, g, flap)

    Written for complex-step: all operations are analytic in the state.
    """
    V, al, be, p, q, r, th, ph = s
    de, da, dr = ctl
    rho, g = prm["rho"], prm["g"]

    u = V * np.cos(al) * np.cos(be)
    v = V * np.sin(be)
    w = V * np.sin(al) * np.cos(be)

    P = engine_power(prm["pz"], prm["n_rpm"], rho)
    dpt = 0.08696 + 191.18 * (2.0 * P / (rho * V**3))

    cx, cy, cz, cl, cm, cn = coeffs(al, be, p, q, r, V, de, da, dr,
                                    prm["flap"], dpt)
    qbar = 0.5 * rho * V * V
    X, Yf, Z = qbar * S * cx, qbar * S * cy, qbar * S * cz
    L, M, N = qbar * S * B * cl, qbar * S * CBAR * cm, qbar * S * B * cn

    udot = r * v - q * w + X / MASS - g * np.sin(th)
    vdot = p * w - r * u + Yf / MASS + g * np.cos(th) * np.sin(ph)
    wdot = q * u - p * v + Z / MASS + g * np.cos(th) * np.cos(ph)

    Vdot = (u * udot + v * vdot + w * wdot) / V
    aldot = (u * wdot - w * udot) / (u * u + w * w)
    bedot = (vdot * V - v * Vdot) / (V * np.sqrt(u * u + w * w))

    gam = IX * IZ - IXZ * IXZ
    rhs_l = L - (IZ - IY) * q * r + IXZ * p * q
    rhs_n = N - (IY - IX) * p * q - IXZ * q * r
    pdot = (IZ * rhs_l + IXZ * rhs_n) / gam
    rdot = (IX * rhs_n + IXZ * rhs_l) / gam
    qdot = (M - (IX - IZ) * p * r - IXZ * (p * p - r * r)) / IY

    thdot = q * np.cos(ph) - r * np.sin(ph)
    phdot = p + (q * np.sin(ph) + r * np.cos(ph)) * np.tan(th)

    return np.array([Vdot, aldot, bedot, pdot, qdot, rdot, thdot, phdot])


def complex_step_jacobian(f, s0):
    n = len(s0)
    A = np.zeros((n, n))
    h = 1e-30
    for j in range(n):
        sp = np.array(s0, dtype=complex)
        sp[j] += 1j * h
        A[:, j] = np.imag(f(sp)) / h
    return A


def load_matrix(path):
    with open(path) as f:
        return np.array([[float(x) for x in row] for row in csv.reader(f)])


def match_eigs(a, b):
    """Greedy nearest-pair matching; returns max pair distance."""
    b = list(b)
    worst = 0.0
    for lam in a:
        j = int(np.argmin([abs(lam - mu) for mu in b]))
        worst = max(worst, abs(lam - b.pop(j)))
    return worst


def main():
    prefix = sys.argv[1] if len(sys.argv) > 1 else "results/beaver"
    # C++ state order (types.hpp): u v w p q r phi theta psi h y
    IDX8 = [0, 1, 2, 3, 4, 5, 7, 6]  # -> [u v w p q r theta phi]
    ok = True
    for cond in ["cruise_45", "check_35", "approach_35"]:
        with open(f"{prefix}_trim_{cond}.csv") as f:
            rows = list(csv.DictReader(f))
        t = {k: float(v) for k, v in rows[0].items()}
        prm = dict(pz=t["pz"], n_rpm=t["n_rpm"], rho=t["rho"], g=t["g"],
                   flap=t["flap"])
        s0 = np.array([t["V"], t["alpha"], t["beta"], 0.0, 0.0, 0.0,
                       t["theta"], 0.0])
        ctl = [t["de"], t["da"], t["dr"]]

        # Sanity: the C++ trim must also be an equilibrium of THIS model.
        res = np.max(np.abs(xdot(s0, ctl, prm)))
        A_py = complex_step_jacobian(lambda s: xdot(s, ctl, prm), s0)
        eig_py = np.linalg.eigvals(A_py)

        A_cpp_full = load_matrix(f"{prefix}_A_{cond}.csv")
        A_cpp = A_cpp_full[np.ix_(IDX8, IDX8)]
        eig_cpp = np.linalg.eigvals(A_cpp)

        worst = match_eigs(eig_cpp, eig_py)
        scale = max(abs(eig_cpp))
        rel = worst / scale
        status = "OK" if (rel < 1e-8 and res < 1e-8) else "MISMATCH"
        if status != "OK":
            ok = False
        print(f"[{status}] {cond}: trim residual {res:.2e}, "
              f"max eigenvalue diff {worst:.2e} (rel {rel:.2e})")
        srt = sorted(eig_py, key=lambda z: (round(z.real, 6), abs(z.imag)))
        print("        modes:", "  ".join(
            f"{z.real:+.4f}{z.imag:+.4f}j" for z in srt if z.imag >= 0))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
