# Water-Landing CBF — Design Document

**Status:** living document — the canonical, maintained statement of *what we are building and
why*. Update it when the design, approach, or assumptions change.

## Purpose & doc map

This doc captures the **current** theory, approach, assumptions, and key decisions for the
augmented-longitudinal CBF-QP water-landing safety filter. It deliberately does **not**
reproduce the derivations — those live in the math spec.

| Doc | Job |
|---|---|
| **`water_landing_cbf_design.md`** (this) | Current design: theory, approach, assumptions, decisions, open questions. |
| `water_landing_cbf_math.md` | Formal derivations: EOM, HOCBF relative-degree analysis, control-affine forms, QP. |
| `CHANGELOG.md` | Chronological "what changed, when, who, why." |
| `archive/` | Superseded historical docs (original pitch; the 2026-06-25 implementation notes). |

---

## 1. What this is

A **safety filter** for autonomous fixed-wing **water landing**. A trivial nominal controller
flies a constant-glide-slope approach; a CBF-QP minimally adjusts its commands to guarantee the
aircraft stays in a provably safe set — hull-safe touchdown sink rate and stall margin — passing
the nominal through untouched until a constraint is about to bind. The flare is **not** designed
into the nominal; it **emerges** from the descent-rate barrier near the surface.

Scope is **longitudinal (pitch-plane) only**. Lateral/heading dynamics, wave-phase exploitation,
and the post-touchdown on-water regime are out of scope (see `archive/water_landing_cbf.md`).

This is a separate, self-contained `lon_*` module living alongside the original 11-state
body-axis code (`dynamics.cpp`, `cbf.cpp`, `sim.cpp`), which it leaves untouched.

## 2. Augmented longitudinal model

- **State** `X = [h, V, γ, θ, q, T, Ṫ]` (NXA = 7); **control** `U = [δe, T̈]` (NUA = 2).
- Thrust is **integrator-augmented**: `T` and `Ṫ` become states, `T̈` the control. This aligns
  the relative degree of thrust with that of elevator (both **degree 3** for the descent and
  airspeed barriers), so a single well-posed QP can share authority between them.
- Code: `include/autoland/lon_augmented.hpp` (`LonDrift`, `gMatrix`, `AeroLocal`).
- **Aero frame:** the EOM are in wind/path axes (need `L`, `D`), while VSPAero data is body-axis.
  We rotate body-axis coefficients into the wind frame analytically
  (`C_L = −CFx·sinα + CFz·cosα`, `C_D = CFx·cosα + CFz·sinα + CD0`, β = 0), verified against the
  `.stab` data. Pitch moment stays body-axis.

Full EOM and symbols: `water_landing_cbf_math.md` §1.

## 3. Nominal controller (current implementation)

A minimal, platform-agnostic approach controller — **constant-thrust powered approach** plus a
cascade holding a constant flight-path angle:

- **Thrust:** PD on the augmented thrust state toward a constant setpoint,
  `T̈_nom = Kp_T (T_set − T) − Kd_T Ṫ`.
- **Elevator:** outer PI `γ → θ_cmd` (limited), inner PD `θ → δe` with pitch-rate damping.

This replaces the math spec's TECS placeholder (§2 there) by directive. The point is that the
nominal only has to put the aircraft on a reasonable glide slope; the CBF owns the
safety-critical shaping. Code: `include/autoland/lon_nominal.hpp`.

## 4. Barriers

| Barrier | Definition | Rel. degree | Enforcement |
|---|---|---|---|
| Descent-rate `b` | `V·sinγ + √(v_safe² + 2·a_brk·h)` | 3 (δe & T̈) | **Hard** |
| Airspeed `b_V` | `V − V_min` | 3 (δe & T̈) | Soft (slack) |
| Min thrust | `T ≥ 0` | 2 (T̈) | Hard (closed form) |
| Max thrust | `T_max − T ≥ 0` | 2 (T̈) | Hard (closed form) |
| Contact-force `b_comp` | `V·sinγ + √(v_safe²(θ) + 2·a_brk·h)` | 2 (δe), 3 (T̈) | **Deferred — not implemented** |

The descent and airspeed HOCBF rows are built generically from the exact Lie stack via
`hocbfRow()`; the thrust actuator barriers use closed forms. Code: `include/autoland/hocbf.hpp`,
`src/lon_cbf_filter.cpp`. Derivations and control-affine authority terms: math spec §3.

**The flare emerges from the degree-3 descent barrier.** `b, ḃ, b̈` are control-free; the
elevator first appears in `b⃛`. The class-K cascade `ψ₃ = ψ̇₂ + c₃ψ₂ ≥ 0` is *predictive*, so the
QP must begin pitching up early — that early pitch-up *is* the flare. A relative-degree-1 view
of `b` cannot flare (`L_g b = 0`); the HOCBF is what makes the constraint enforceable.

## 5. Method: exact Lie derivatives (no finite differences)

The HOCBF rows need `L_f^k b` and the control row `L_g L_f^{r-1} b` over the augmented dynamics.
These are computed **exactly**, not by central differences (contrast: `main`'s `numGrad` and
central-difference linearization):

- **Drift stack** `{b, …, L_f^R b}` via the **flow Taylor-jet** — Picard iteration in truncated-
  Taylor arithmetic, reading `L_f^k b = k!·[t^k] b(X(t))`.
- **Control row** `L_g L_f^{r-1} b` via forward-mode autodiff (`autodiff::dual`) seeded along a
  g-column.
- Code: `include/autoland/lie_taylor.hpp`. Dependency: **autodiff** (header-only, FetchContent).

**"Exact" has a boundary:** it is exact differentiation of a **frozen local-affine aero** model
(`AeroLocal` freezes the table slopes at the evaluation point), so 2nd/3rd Lie derivatives zero
out true aero curvature. Exact math ≠ exact physics.

## 6. Key design decisions

- **`a_brk` is a constant set-shaping parameter, not a measured acceleration.** The barrier
  `√(v_safe² + 2·a_brk·h)` is the kinematic braking envelope (boundary of the maximal
  control-invariant set of a double integrator with bounded deceleration `a_brk`). It defines
  *which set we stay in*. Currently `a_brk = 3.0 m/s²` (config).
- **Drag and thrust ARE accounted for — in the exact derivatives, not in `a_brk`.** The HOCBF
  differentiates through the full wind-axis EOM (lift, drag, thrust), so the *enforcement* side
  carries the real forces even though the *set-shaping* constant does not. This is the cleaner
  split: kinematic envelope shapes the set; exact dynamics realize the braking, with the
  degree-3 HOCBF handling the lift-build-up lag (the relative-degree problem).
- **Open:** make `a_brk` itself force-based and speed-dependent,
  `a_brk(V) ≈ ρV²S·C_L,max/(2m)·cosγ − g` (+ drag/thrust), the "honest form" the original doc
  calls for. This couples the descent and airspeed barriers and still requires a stall ceiling
  `C_L,max` (the inviscid aero has no lift cap). It addresses magnitude, not the lag.
- **Doc-vs-autodiff divergence (recorded as a test):** the math spec's closed-form elevator
  authority omits the pitch-rate aero terms (`∂C_L/∂q̂`, `∂C_D/∂q̂`); the autodiff captures them.
  They differ by a few percent. The autodiff quantity is the more complete one.
- **Touchdown sink is a tuning consequence, not a barrier property.** A safety *filter* yields
  `sink ≤ v_safe` with slack; the exact value follows from the class-K gains. A *defined*
  touchdown sink belongs in a nominal flare reference, with the CBF as backstop.
- **No `h ≥ 0` CBF.** Any class-K has `α(0)=0`, so an `h ≥ 0` barrier forces `sink → 0` at the
  surface — it would forbid touchdown. The descent-rate barrier (value `v_safe > 0` at `h=0`) is
  precisely the ground-proximity condition in the only form compatible with landing.

## 7. Vehicle & scenario data (current)

As of `2b610a3`, the model uses **real AHAB vehicle data**, not placeholders:

- `data/aircraft.yaml`: mass 3.6139 kg; Iyy 0.0521 kg·m²; real Ixx/Izz; thrust `T_static = 50 N`,
  `zcp = 0.15 m`; `parasite_CD0 = 0.030`. (Thrust `k_v` and `parasite_CD0` still uncalibrated.)
- `data/AHAB_combined.stab`: real combined VSPAero aero deck.
- `data/lon_scenario.yaml`: powered approach `T_set = 2.0 N`, `V_app = 18`, `γ = −3°`,
  `v_safe = 0.1`, `a_brk = 3.0`, `Vmin = 13.5`, `Tmax = 50`, descent class-K `[2,2,2]`; gains
  retuned for the real mass.

## 8. Assumptions & limitations

Folded from the 2026-06-25 implementation notes (`archive/`), updated for current state:

1. **Self-consistent plant.** The sim plant *is* the CBF model (zero model mismatch), so clean
   barrier-invariance partly reflects structure; it validates the math, not robustness to model
   error.
2. **`L_f³b` drift not independently verified.** Tests check the *control* authorities vs the
   spec; the spec never writes `L_f³b`, so a drift error would pass. (Follow-up: a test-only
   finite-difference cross-check as an oracle.)
3. **Frozen local-affine aero** — autodiff exactness is relative to a C0 piecewise-linear table
   surrogate; aero curvature is lost in 2nd/3rd derivatives.
4. **Initial condition is not a true longitudinal-model trim** (seeded from the 11-state
   body-axis trim); airspeed drifts before the flare. Follow-up: Newton trim on the lon model.
5. **No airspeed control in the nominal** (constant-thrust directive) — airspeed leans on the
   airspeed CBF; a high `Vmin` parks the energy-limited glider in a safe hover instead of landing.
6. **Soft airspeed barrier** (slack-penalized) — not a hard guarantee. Descent + thrust are hard.
7. **Discrete-time / ZOH** at dt = 0.01 approximates the continuous-time CBF guarantee.
8. **No controllability guard** on `L_g L_f²b` (elevator authority → 0 at very low dynamic
   pressure would make the descent barrier unenforceable by elevator; not yet observed).
9. **Idealized scope:** longitudinal-only, perfect full-state feedback, no wind/gust/sensor
   noise, instantaneous elevator (no actuator lag).
10. **§3.3 deferred ⇒ "hull-safe" is incomplete.** `v_safe` is a flat constant and touchdown
    attitude is unconstrained; real slamming load is attitude-dependent.

## 9. Open questions / follow-ups

- **Force-based `a_brk(V)`** (§6) — magnitude honesty + descent/airspeed coupling.
- **§3.3 contact-force barrier** — attitude-coupled `v_safe(θ)` from von Kármán/Wagner slamming
  theory; the contribution most likely to make this a paper.
- **Independent `L_f³b` check** — finite-difference oracle test.
- **True longitudinal trim** — Newton solve on the lon model for a consistent initial condition.
- **Hardware path** — PX4 insertion point (actuator vs rate/attitude setpoint sets the relative
  degree).

---

*See `CHANGELOG.md` for the running log of changes.*
