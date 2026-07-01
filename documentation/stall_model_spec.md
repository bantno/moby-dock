# NACA 4414 Viscous-Stall Plant Model for AHAB Autoland

Implementation reference for the viscous-stall overlay added to the longitudinal
plant. The base aerodynamic deck (`data/AHAB_combined.stab`) is an inviscid
VSPAERO vortex-lattice solution and **cannot stall** — its lift curve is linear to
±20° and its post-stall values are meaningless. This model **hands the lift off**
from VSPAERO (attached flow only) to a viscous + flat-plate post-stall curve so the
plant can stall, depart, and be recovered. Self-contained: it describes the source
data, the blend, and the design decisions; it does not prescribe a code structure.

---

## 1. Purpose

Give the longitudinal plant a physically-grounded stall so that stall-recovery
control work has something to recover *from*. The overlay must:

- reproduce a real CLmax, a post-stall **lift crater** (lift declines toward zero,
  not a plateau), the separation **drag rise**, and the **pitch break**;
- be **inert in attached flow**, so the nominal approach/trim and all existing
  behaviour are unchanged when enabled (and bit-identical when disabled);
- **not use VSPAERO above stall onset** — its inviscid post-stall lift is fiction;
- stay compatible with the exact-Lie-derivative CBF machinery (smooth, frozen
  local-affine, no state-dependent branch in the templated EOM).

It is a **plant** model only. Designing a stall-recovery *barrier* (e.g. an
angle-of-attack-margin HOCBF) is separate, future work.

---

## 2. Source data: NACA 4414 2D viscous polar (NeuralFoil)

The stall physics come from the **2D section** polar of the NACA 4414 (the AHAB
wing section), generated with **NeuralFoil** — a neural-network surrogate of Xfoil
(pure-Python, smooth, differentiable, includes the post-stall branch). Generated
into `data/naca4414_polar.csv` by `scripts/gen_4414_polar.py`.

- Reynolds number: `Re = V·c/ν`. With `c = 0.25 m`, landing `V ≈ 13–18 m/s`,
  `Re ≈ 2.2–3.1×10⁵`. The table is anchored at **Re = 2.5×10⁵** (`CLmax₂d ≈ 1.46`).
- The 2D polar sets the realistic **CLmax level**; the Viterna extrapolation (§3.2)
  provides the deep-stall tail.

### 2.1 Validity and trust (read before calibrating)

NeuralFoil is validated **against Xfoil** (mean ~few-% error), *not* against
experiment. The trust chain is **NeuralFoil → Xfoil → reality**, weakest exactly
where we operate: **low Re and post-stall**, where Xfoil over-predicts CLmax/α_stall
and the deep-stall branch is uncertain. Therefore the stall *onset* and *level* are
**tunable engineering values, not ground truth** (knobs in §3.3); calibrate against
experimental low-Re 4414 data or flight ID when available.

---

## 3. The model: a VSPAERO → Viterna blend (not additive deltas)

The plant blends each wind-axis coefficient `C ∈ {CL, CD, CM}` from the VSPAERO
attached value to an absolute post-stall value:

    C_plant(α) = (1 − w(α))·C_vspaero(α)  +  w(α)·C_post(α)

- **`w(α)`** is a smoothstep blend weight: `0` below the wing stall onset (attached
  flow ≡ the trustworthy VSPAERO deck), ramping to `1` across the onset, then `1`.
  Because the C++ plant computes `C_vspaero` itself, weighting it by `(1−w)→0`
  **discards VSPAERO's meaningless post-stall values** above the onset — VSPAERO is
  thrown out, exactly where it has no validity.
- **`C_post(α)`** is the absolute post-stall curve (§3.2), VSPAERO-independent.

### 3.1 Wing stall onset

The wing leaves the VSPAERO line at `A_STALL` (geometric α), with
`CLmax = CL_vspaero(A_STALL)`. Default `A_STALL = 11°` → **CLmax ≈ 1.44** (kept ≤
the 2D section CLmax 1.46, since a tailed vehicle's CLmax is typically below the
wing-section value). The wing stalls earlier in geometric α than the 2D *section*
because the high VSPAERO lift slope (≈5.6/rad, α₀=−4°) reaches a realistic CLmax by
~11°; the section's own CLmax sits higher (~17–20°) but the vehicle never gets
there.

### 3.2 Viterna flat-plate post-stall tail

`C_post` is the **Viterna & Corrigan** flat-plate extrapolation — the standard
wind-turbine/BEM method for post-stall polars — anchored at the wing stall point
`(A_STALL, CLmax, CD_stall)` so it is continuous with the attached curve and
declines to **CL = 0 at 90°** like a real fully-separated wing:

    CD_max = 1.11 + 0.018·AR           (≈ 1.23 for AR ≈ 6.8)
    CL_post(α) = (CD_max/2)·sin 2α  +  A₂·cos²α / sin α
    CD_post(α) = CD_max·sin²α       +  B₂·cos α

with `A₂, B₂` fixed by continuity at the stall point. So past onset the lift
craters off CLmax and follows `sin 2α` to zero at 90°, while drag rises toward
`CD_max`. `CM_post` moves smoothly from the stall value to a nose-down deep-stall
value (`CM_DEEP ≈ −0.45`, flat-plate centre-of-pressure aft) — a simple, documented
approximation (the section moment reference ≠ the vehicle reference).

### 3.3 Tunable knobs (in `scripts/precompute_stall_table.py`)

| Knob | Default | Effect |
|---|---|---|
| `RE_ANCHOR` | 2.5×10⁵ | which 2D polar Re sets CLmax |
| `A_STALL_DEG` | 11 | wing stall angle / CLmax (= CL_vspaero there) |
| `BLEND_HALF_DEG` | 2.5 | handoff half-width (rounds CLmax; bigger ⇒ softer stall) |
| `AR` | 6.8 | Viterna `CD_max` |
| `CM_DEEP` | −0.45 | deep-stall pitching moment |
| `SEVERITY` | 1.0 | scales the blend weight `w` (0 ⇒ no stall; runtime knob too) |

**Resulting plant stall (defaults):** CLmax ≈ **1.44 at 11°**, cratering to ≈0.95
by 30°, ≈0.80 by 45°, **0 at 90°**; CD → 1.23; nose-down CM break.

---

## 4. Implementation (frozen-affine, shared plant + CBF)

The overlay reuses the existing frozen-local-affine aero path, so the exact-Lie
machinery is preserved:

1. **Generated table** `include/autoland/naca4414_stall_table.hpp` — α-indexed
   `w / CLpost / CDpost / CMpost` arrays out to 90° geometric, emitted by
   `scripts/precompute_stall_table.py` (`--check` prints diagnostics; mirrors
   `precompute_impact_clf.py`). Committed; regenerate by hand.
2. **Lookup** `include/autoland/stall_model.hpp` — `stallLookup(α)` binary-search +
   linear interpolation returning `w` and the post-stall coefficients **with their
   d/dα slopes** (mirrors `clfLookup`). Held with zero slope below the grid (`w=0`,
   attached) and at the last deep-stall values above it.
3. **Freeze** `makeAeroLocal` (`src/lon_augmented.cpp`) — when `cfg.stall.enabled`,
   look up at the eval α and fold value+slope into local-affine `{off + slope·α}`
   models for `w` (× `severity`) and each post coefficient.
4. **Apply** `LonDrift` (`include/autoland/lon_augmented.hpp`) — when `stall_on`,
   form `C = (1−w)·C_vspaero + w·C_post` for CL/CD/CMy. The disabled path runs
   **zero extra ops** (bit-identical). `w` is a per-evaluation constant freeze, not
   a state branch, so the Taylor/autodiff jet stays valid and the affine blend keeps
   the in-cell Lie derivatives exact.

**Shared with the CBF:** because it lives in `makeAeroLocal`, both the integrated
plant and the barrier Lie derivatives see the (locally-linearized) blend. In
attached flight `w=0`, so the CBF is unaffected; past stall its local model matches
the plant — the right foundation for a future α-margin barrier.

### 4.1 Configuration

`AircraftConfig::stall` (`StallParams{enabled=false, severity=1.0}`), parsed from a
`stall:` block in `data/aircraft.yaml` (default **OFF**). A scenario YAML may
override it (see `data/lon_stall_recovery.yaml`). `severity` scales `w` at runtime
(0 ⇒ off, 1 ⇒ as-generated).

### 4.2 Diagnostics

`lon_sim` logs `CL, CD, dCL_stall` (the blended wind-axis coefficients + the lift
lost to stall = blended − attached) per step. `scripts/plot_stall_model.py` overlays
the plant CL/CD/CM(α) against the VSPAERO deck → `figures/stall_model_check.png`.

---

## 5. Verification

- **Default-off / attached-flow invariance:** the nominal approach is **bit-identical**
  with stall on vs off (`max|Δh| = 0`, `max|Δδe| = 0`; `dCL_stall ≡ 0`, max α ≈ 2°).
- **Unit tests** (`test/test_stall_model.cpp`): `stallLookup` inert-below / held-above
  / linear-interp / post-stall lift-decline-to-zero + drag-rise; and the blend at a
  deep-stall α hands CL/CD to the absolute Viterna values and cuts the lift.
- **End-to-end departure** (`data/lon_stall_recovery.yaml`, CBF off): an aggressive
  low-altitude pitch-up climbs to CLmax ≈ 1.42 at α ≈ 11°, then the lift **craters**
  (0.90 at 28°, 0.45 at 69°, 0 at 97°), CD → 1.23, and the sink builds to ~10 m/s — a
  real stall departure (the un-recovered fall tumbles past 90°; see §6).

---

## 6. Caveats & follow-ups

- Stall onset/level are tunable, not experiment-validated (§2.1).
- **Valid to ~90° geometric α.** Above 90° the table clamps (CL held at 0, CD at
  CD_max); a fully-departed tumble past 90° is outside the model (a flat plate
  reverses lift there). Fine for stall *onset/recovery* studies, which live on the
  front side of the curve.
- Single anchor Reynolds number; a 2nd table axis over Re is a follow-up.
- Static (memoryless) model. A Goman–Krabrov separation **state** would add
  dynamic-stall lag/hysteresis (and would double as a CBF state) — future work.
- The descent barrier's separate `CL_max` knob (`LonCBFConfig`) is *not* auto-linked
  to this model's CLmax (~1.44); update it if you want the barrier's max-lift
  assumption to match the plant.
