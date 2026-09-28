#!/usr/bin/env python3
"""Plot a corridor landing run (sixdof_autoland_sim CSV with a world: block).

Panels: ground track over the terrain height field with the corridor
rectangle, the planned path (lookahead trace), touchdown and stop markers and
the rollout track; altitude vs distance-to-go against the glideslope cone and
the terrain under the track; corridor-frame cross-track e(s); bank (flown vs
commanded); course (flown vs commanded) ; airspeed + flight-path angle;
rollout ground speed.

Usage: plot_corridor_landing.py RUN.csv SCENARIO.yaml [out.png]
The terrain and corridor are re-evaluated from the scenario YAML (same
super-Gaussian bump formula as world_geometry.hpp).
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

DEG = np.pi / 180.0
# Fixed categorical roles (never cycled): flown / commanded / reference / marker.
C_FLOWN, C_CMD, C_REF, C_TD, C_STOP = "#1f5fa8", "#d1691f", "#5c6b73", "#b0234f", "#2a8a5f"


def load_csv(path):
    d = np.genfromtxt(path, delimiter=",", names=True)
    if d.size == 0:
        sys.exit(f"no data in {path}")
    return d


def terrain_from_yaml(w):
    """Return (bumps, base) with bumps as dicts, rings expanded."""
    t = w.get("terrain", {}) or {}
    base = float(t.get("base", 0.0))
    bumps = []
    for b in t.get("bumps", []) or []:
        bumps.append(dict(N=float(b["N"]), E=float(b["E"]), A=float(b["height"]),
                          sa=float(b.get("sigma_a", 10.0)),
                          sb=float(b.get("sigma_b", b.get("sigma_a", 10.0))),
                          rot=float(b.get("rot_deg", 0.0)) * DEG,
                          order=int(b.get("order", 1))))
    for r in t.get("rings", []) or []:
        n = int(r.get("count", 32))
        gap_h = r.get("gap_heading_deg", None)
        gap_w = float(r.get("gap_halfangle_deg", 0.0)) * DEG
        for i in range(n):
            th = 2 * np.pi * i / n
            if gap_h is not None:
                gh = float(gap_h) * DEG
                d1 = np.remainder(th - gh + np.pi, 2 * np.pi) - np.pi
                d2 = np.remainder(th - gh, 2 * np.pi) - np.pi
                if abs(d1) <= gap_w or abs(d2) <= gap_w:
                    continue
            bumps.append(dict(N=float(r.get("center_N", 0)) + float(r["radius"]) * np.cos(th),
                              E=float(r.get("center_E", 0)) + float(r["radius"]) * np.sin(th),
                              A=float(r.get("height", 15.0)), sa=float(r.get("sigma", 25.0)),
                              sb=float(r.get("sigma", 25.0)), rot=0.0,
                              order=int(r.get("order", 2))))
    return bumps, base


def terrain_height(bumps, base, N, E):
    h = np.full_like(np.asarray(N, dtype=float), base)
    for b in bumps:
        c, s = np.cos(b["rot"]), np.sin(b["rot"])
        dN, dE = N - b["N"], E - b["E"]
        xa = (dN * c + dE * s) / b["sa"]
        xb = (dE * c - dN * s) / b["sb"]
        q = xa * xa + xb * xb
        h = h + b["A"] * np.exp(-(q ** b["order"]))
    return h


def corridor_poly(c):
    Nc, Ec = float(c.get("center_N", 0)), float(c.get("center_E", 0))
    L, W = float(c["length"]), float(c["width"])
    psi = float(c.get("heading_deg", 0.0)) * DEG
    cs, sn = np.cos(psi), np.sin(psi)
    pts = []
    for s, e in [(-L / 2, -W / 2), (L / 2, -W / 2), (L / 2, W / 2), (-L / 2, W / 2), (-L / 2, -W / 2)]:
        pts.append((Nc + s * cs - e * sn, Ec + s * sn + e * cs))
    return np.array(pts), (Nc, Ec, psi)


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    csv_path, yaml_path = sys.argv[1], sys.argv[2]
    out = sys.argv[3] if len(sys.argv) > 3 else "corridor_landing.png"
    d = load_csv(csv_path)
    with open(yaml_path) as f:
        sc = yaml.safe_load(f)
    w = sc["world"]
    bumps, base = terrain_from_yaml(w)
    poly, (Nc, Ec, psi_c) = corridor_poly(w["corridor"])

    air = d[d["phase"] == 0]
    roll = d[d["phase"] == 1]
    td = air[-1]

    fig, ax = plt.subplots(2, 3, figsize=(17, 10))
    fig.suptitle(f"Corridor landing: {sc.get('plant', 'beaver')} plant, pattern nominal "
                 f"(no CBF) -- {csv_path}", fontsize=13, weight="bold")

    # (1) Ground track over terrain.
    a = ax[0, 0]
    # Extent: the flown track + corridor (distant terrain is clipped, not
    # squeezed in -- the equal-aspect panel would collapse otherwise).
    allN = np.concatenate([air["x"], roll["x"], poly[:, 0]])
    allE = np.concatenate([air["y"], roll["y"], poly[:, 1]])
    pad = 250.0
    Nlo, Nhi = allN.min() - pad, allN.max() + pad
    Elo, Ehi = allE.min() - pad, allE.max() + pad
    NN, EE = np.meshgrid(np.linspace(Nlo, Nhi, 300), np.linspace(Elo, Ehi, 300), indexing="ij")
    HT = terrain_height(bumps, base, NN, EE)
    cf = a.contourf(EE, NN, HT, levels=[0.5, 5, 15, 50, 100, 200, 400], cmap="YlOrBr", alpha=0.6)
    fig.colorbar(cf, ax=a, label="terrain height [m]", shrink=0.8)
    a.plot(air["E_ref"], air["N_ref"], color=C_REF, lw=1.0, ls="--", label="planned path (lookahead)")
    a.plot(air["y"], air["x"], color=C_FLOWN, lw=1.6, label="flown")
    if roll.size:
        a.plot(roll["y"], roll["x"], color=C_STOP, lw=2.0, label="rollout")
        a.plot(roll["y"][-1], roll["x"][-1], "s", color=C_STOP, ms=8, label="stop")
    a.plot(poly[:, 1], poly[:, 0], color="k", lw=1.2, label="corridor")
    a.plot(td["y"], td["x"], "v", color=C_TD, ms=9, label="touchdown")
    a.plot(air["y"][0], air["x"][0], "o", color=C_FLOWN, ms=7, mfc="none", label="start")
    a.set_xlabel("E [m]")
    a.set_ylabel("N [m]")
    a.set_aspect("equal")
    a.set_xlim(Elo, Ehi)
    a.set_ylim(Nlo, Nhi)
    a.set_title("Ground track over terrain")
    a.legend(fontsize=7, loc="upper left", bbox_to_anchor=(0.0, -0.12), ncol=3, frameon=False)

    # (2) Altitude vs distance-to-go with the cone and terrain under track.
    a = ax[0, 1]
    L = air["L_to_go"]
    a.plot(L, air["h"], color=C_FLOWN, lw=1.6, label="h flown")
    a.plot(L, air["h_gs"], color=C_CMD, lw=1.2, ls="--", label="glideslope cone h_gs")
    a.fill_between(L, 0, air["hT"], color="#a56d2b", alpha=0.5, label="terrain under track")
    a.axvline(0, color="k", lw=0.8, ls=":")
    a.invert_xaxis()
    a.set_xlabel("path distance to aim point [m]")
    a.set_ylabel("altitude [m]")
    a.set_title("Descent vs the glideslope cone")
    a.legend(fontsize=8)

    # (3) Corridor frame: e(s) on final + rollout.
    a = ax[0, 2]
    fin = air[air["seg"] == air["seg"].max()]
    a.plot(fin["s_corr"], fin["e_corr"], color=C_FLOWN, lw=1.6, label="final (airborne)")
    if roll.size:
        a.plot(roll["s_corr"], roll["e_corr"], color=C_STOP, lw=2.0, label="rollout")
    Wc, Lc = float(w["corridor"]["width"]), float(w["corridor"]["length"])
    for e in (-Wc / 2, Wc / 2):
        a.axhline(e, color="k", lw=1.0)
    for s in (-Lc / 2, Lc / 2):
        a.axvline(s, color="k", lw=1.0)
    a.axvline(float(w["corridor"].get("aim_s", -Lc / 2 + 100)), color=C_TD, lw=0.9, ls=":", label="aim s")
    a.plot(td["s_corr"], td["e_corr"], "v", color=C_TD, ms=9)
    a.set_xlabel("along-corridor s [m]")
    a.set_ylabel("cross-corridor e [m]")
    a.set_title("Corridor frame (final + rollout)")
    a.legend(fontsize=8)

    # (4) Bank flown vs commanded.
    a = ax[1, 0]
    a.plot(air["t"], air["phi_deg"], color=C_FLOWN, lw=1.4, label=r"$\phi$")
    a.plot(air["t"], air["phi_cmd_deg"], color=C_CMD, lw=1.0, ls="--", label=r"$\phi_{cmd}$ (L1)")
    a.plot(air["t"], air["beta_deg"], color=C_REF, lw=0.9, label=r"$\beta$")
    a.set_xlabel("t [s]")
    a.set_ylabel("[deg]")
    a.set_title("Bank command tracking, sideslip")
    a.legend(fontsize=8)

    # (5) Course flown vs commanded, cross-track to path.
    a = ax[1, 1]
    chi = np.degrees(np.unwrap(np.arctan2(np.gradient(air["y"]), np.gradient(air["x"]) + 1e-12)))
    a.plot(air["t"], chi, color=C_FLOWN, lw=1.4, label=r"course $\chi$")
    a.plot(air["t"], np.degrees(np.unwrap(np.radians(air["chi_cmd_deg"]))), color=C_CMD, lw=1.0, ls="--",
           label="bearing to L1 point")
    a.set_xlabel("t [s]")
    a.set_ylabel("[deg]")
    a.set_title(f"Course; cross-track to path max {air['xte'].max():.0f} m")
    a.legend(fontsize=8)

    # (6) Airspeed, gamma, and rollout ground speed on a common time axis.
    a = ax[1, 2]
    a.plot(air["t"], air["V_air"], color=C_FLOWN, lw=1.4, label="V_air [m/s]")
    a.plot(air["t"], air["gamma_deg"], color=C_CMD, lw=1.0, label=r"$\gamma$ [deg]")
    if roll.size:
        a.plot(roll["t"], roll["Vg"], color=C_STOP, lw=2.0, label="rollout Vg [m/s]")
    a.axvline(td["t"], color=C_TD, lw=0.9, ls=":")
    a.set_xlabel("t [s]")
    a.set_title("Speed and flight-path angle; rollout")
    a.legend(fontsize=8)

    for row in ax:
        for a in row:
            a.grid(alpha=0.25)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
