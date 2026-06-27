#!/usr/bin/env python3
"""Splice control-surface derivative columns from a small donor .stab into a
full (alpha,beta)-sweep .stab that was exported WITHOUT control groups.

Why this exists
---------------
A VSPAERO stability .stab carries the control-surface effectiveness as extra
"ConGrp_k" columns in the derivative table (and Ailerons/Elevator/Rudder rows in
each Case block). If the sweep was run without the control-deflection cases,
those columns are missing and the aircraft has ZERO control authority in the sim
(elevator does nothing -> the CBF-QP can't flare). Control effectiveness is
nearly constant in (alpha,beta) for attached flow, so it does NOT need the full
grid: run ONE (or a few-alpha) VSPAERO stability solve WITH the control groups
defined, and broadcast those columns across every block of the big sweep.

Usage
-----
    splice_control_derivs.py DONOR.stab TARGET.stab OUTPUT.stab

  DONOR   small .stab WITH control groups (>=1 block; ideally a beta=0 alpha
          sub-sweep so control power can be alpha-interpolated).
  TARGET  the full sweep missing control columns (e.g. "data/AHAB 3.stab").
  OUTPUT  written .stab with ConGrp columns + Ailerons/Elevator/Rudder case rows
          spliced into every block.

The donor's control columns are linearly interpolated in alpha (clamped at the
ends; constant if the donor has a single alpha) and held constant in beta. The
injected Case rows exist only to give the control groups their NAMES so the
name-based mixing map (mixing.cpp) resolves "Elevator"->delta_e; their numeric
content is copied verbatim from the donor and is otherwise unused by the parser.
"""
import sys

KNOWN_CASE_ROWS = {"Base_Aero", "Alpha", "Beta", "Roll__Rate",
                   "Pitch_Rate", "Yaw___Rate", "Mach"}
# Derivative-table coefficient row labels (same set aero_table.cpp recognises).
COEF_LABELS = {"CFx", "CFy", "CFz", "CMx", "CMy", "CMz",
               "CL", "CD", "CS", "CMl", "CMm", "CMn"}


def parse_donor(path):
    """Return (alpha_deg -> {coef_label: [ctrl numbers]}, header_tail_names,
    [verbatim control Case-row lines]). Only beta~=0 blocks are used."""
    with open(path, encoding="latin-1") as f:
        lines = f.readlines()

    by_alpha = {}          # alpha_deg -> {label: [floats]}
    header_tail = None     # e.g. ["ConGrp_1","ConGrp_2","ConGrp_3"]
    ctrl_case_lines = []   # verbatim Ailerons/Elevator/Rudder Case rows
    n_std = 8              # Total + 7 std derivs before control columns

    alpha = beta = None
    section = None         # "case" | "deriv"
    cur = None
    for raw in lines:
        if "****" in raw:
            if cur is not None and beta is not None and abs(beta) < 1e-6:
                by_alpha.setdefault(round(alpha, 4), cur)
            alpha = beta = None
            section = None
            cur = None
            continue
        s = raw.split()
        if not s:
            continue
        if raw.lstrip().startswith("AoA_"):
            alpha = float(s[1])
        elif raw.lstrip().startswith("Beta_"):
            beta = float(s[1])
        elif s[0] == "Case":
            section = "case"
        elif s[0] == "Coef":
            section = "deriv"
            if header_tail is None and len(s) > 1 + n_std:
                header_tail = s[1 + n_std:]   # ConGrp_* names
            if cur is None:
                cur = {}
        elif section == "case":
            # Control case rows exist only to supply the group NAMES, so keep one
            # representative row per name (the first beta=0 occurrence).
            if s[0] not in KNOWN_CASE_ROWS and s[0] not in COEF_LABELS:
                seen = {ln.split()[0] for ln in ctrl_case_lines}
                if abs((beta if beta is not None else 1)) < 1e-6 \
                        and s[0] not in seen:
                    ctrl_case_lines.append(raw.rstrip("\n"))
        elif section == "deriv" and s[0] in COEF_LABELS:
            if len(s) > 1 + n_std:
                if cur is None:
                    cur = {}
                cur[s[0]] = [float(x) for x in s[1 + n_std:]]

    # Flush the final block (no trailing "****" separator after the last one).
    if cur is not None and beta is not None and abs(beta) < 1e-6:
        by_alpha.setdefault(round(alpha, 4), cur)

    if not by_alpha:
        sys.exit("donor: no beta=0 blocks with control columns found")
    if header_tail is None:
        sys.exit("donor: no ConGrp columns in the derivative-table header")
    return by_alpha, header_tail, ctrl_case_lines


def interp(by_alpha, label, alpha_deg):
    """Linear interpolation in alpha (clamped), constant if single point."""
    pts = sorted(a for a in by_alpha if label in by_alpha[a])
    if not pts:
        return None
    if alpha_deg <= pts[0]:
        return by_alpha[pts[0]][label]
    if alpha_deg >= pts[-1]:
        return by_alpha[pts[-1]][label]
    for i in range(len(pts) - 1):
        a0, a1 = pts[i], pts[i + 1]
        if a0 <= alpha_deg <= a1:
            w = (alpha_deg - a0) / (a1 - a0)
            v0, v1 = by_alpha[a0][label], by_alpha[a1][label]
            return [x0 + w * (x1 - x0) for x0, x1 in zip(v0, v1)]
    return by_alpha[pts[-1]][label]


def fmt(vals):
    return "".join(f" {v:>12.7f}" for v in vals)


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    donor, target, out = sys.argv[1:4]
    by_alpha, header_tail, ctrl_lines = parse_donor(donor)

    with open(target, encoding="latin-1") as f:
        lines = f.readlines()

    out_lines = []
    alpha = None
    section = None
    injected = False        # control Case rows added in this block yet?
    for raw in lines:
        if "****" in raw:
            alpha = None
            section = None
            injected = False
            out_lines.append(raw)
            continue
        s = raw.split()
        if s and raw.lstrip().startswith("AoA_"):
            alpha = float(s[1])
        elif s and s[0] == "Case":
            section = "case"
        elif s and s[0] == "Coef":
            section = "deriv"
            out_lines.append(raw.rstrip("\n") + "   " + "   ".join(header_tail) + "\n")
            continue
        elif section == "case" and s and s[0] == "Mach" and not injected:
            # Emit the Mach case row, then inject the control case rows (names).
            out_lines.append(raw)
            for cl in ctrl_lines:
                out_lines.append(cl + "\n")
            injected = True
            continue
        elif section == "deriv" and s and s[0] in COEF_LABELS:
            vals = interp(by_alpha, s[0], alpha if alpha is not None else 0.0)
            if vals is not None:
                out_lines.append(raw.rstrip("\n") + fmt(vals) + "\n")
                continue
        out_lines.append(raw)

    with open(out, "w", encoding="latin-1") as f:
        f.writelines(out_lines)
    print(f"wrote {out}: control groups {header_tail}, "
          f"{len(by_alpha)} donor alpha point(s), "
          f"{'alpha-interpolated' if len(by_alpha) > 1 else 'broadcast constant'}")


if __name__ == "__main__":
    main()
