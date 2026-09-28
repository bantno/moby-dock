#!/usr/bin/env python3
"""Calibrate the water-rollout deceleration (world.rollout.a0, kq) against
published landing-run data.

The rollout model (include/autoland/water_rollout.hpp) is
    Vg' = -(a0 + kq Vg^2)
whose stopping distance from touchdown speed V0 is
    d(V0) = ln((a0 + kq V0^2) / a0) / (2 kq).

Usage:
  calibrate_rollout.py --point V0,d [--point V0,d ...] [--kq KQ]
Examples:
  # one POH point (touchdown speed [m/s], water run [m]) with kq fixed:
  calibrate_rollout.py --point 33.5,320 --kq 0.002
  # two points: solve a0 and kq jointly
  calibrate_rollout.py --point 33.5,320 --point 40.0,450

NOTE: the values shipped in data/beaver_corridor_*.yaml (a0 = 0.8 m/s^2,
kq = 0.002 /m) are PLACEHOLDERS. Supply the DHC-2 float POH landing-run
distance (calm water, sea level, gross weight) at its recommended touchdown
speed to replace them; the script prints the YAML lines.
"""
import argparse
import math

import numpy as np
from scipy.optimize import least_squares


def d_stop(V0, a0, kq):
    if kq <= 0:
        return V0 * V0 / (2 * a0)
    return math.log((a0 + kq * V0 * V0) / a0) / (2 * kq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--point", action="append", required=True,
                    help="V0[m/s],distance[m] (repeatable)")
    ap.add_argument("--kq", type=float, default=None,
                    help="fix kq [1/m] and solve a0 only")
    a = ap.parse_args()
    pts = [tuple(float(v) for v in p.split(",")) for p in a.point]
    if a.kq is not None or len(pts) == 1:
        kq = a.kq if a.kq is not None else 0.002
        V0, d = pts[0]
        # d = ln((a0 + kq V0^2)/a0)/(2kq)  ->  a0 = kq V0^2 / (exp(2 kq d) - 1)
        a0 = kq * V0 * V0 / (math.exp(2 * kq * d) - 1.0)
        print(f"kq fixed at {kq:.5f} /m -> a0 = {a0:.4f} m/s^2   (check d({V0}) = {d_stop(V0, a0, kq):.1f} m)")
    else:
        def res(p):
            a0, kq = p
            return [d_stop(V0, a0, kq) - d for V0, d in pts]
        sol = least_squares(res, x0=[0.8, 0.002], bounds=([1e-3, 1e-5], [10.0, 0.1]))
        a0, kq = sol.x
        print(f"a0 = {a0:.4f} m/s^2   kq = {kq:.5f} /m   residuals: "
              + ", ".join(f"{r:+.1f} m" for r in res(sol.x)))
    print("\nrollout:\n  enabled: true\n  a0: %.4f\n  kq: %.5f\n  k_psi: 0.5\n  V_stop: 0.5\n  t_max: 90.0" % (a0, kq))


if __name__ == "__main__":
    main()
