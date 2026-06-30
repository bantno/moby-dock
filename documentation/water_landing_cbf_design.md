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
| Airspeed (stall) `b_V` | `V − V_min` | 3 (δe & T̈) | Soft (slack) |
| Airspeed (over-speed) `b_V,max` | `V_max − V` | 3 (δe & T̈) | Soft (slack) |
| Min thrust | `T ≥ 0` | 2 (T̈) | Hard (closed form) |
| Max thrust | `T_max − T ≥ 0` | 2 (T̈) | Hard (closed form) |
| Impact-load `b_imp` | `(n_limit − n_peak(θ,γ,V)) + Φ(z)` | **2 (δe)** | Soft (height-sched. slack) |

The descent and airspeed HOCBF rows are built generically from the exact Lie stack via
`hocbfRow()`; the thrust actuator barriers use closed forms. Code: `include/autoland/hocbf.hpp`,
`src/lon_cbf_filter.cpp`. Derivations and control-affine authority terms: math spec §3.

**The flare emerges from the degree-3 descent barrier.** `b, ḃ, b̈` are control-free; the
elevator first appears in `b⃛`. The class-K cascade `ψ₃ = ψ̇₂ + c₃ψ₂ ≥ 0` is *predictive*, so the
QP must begin pitching up early — that early pitch-up *is* the flare. A relative-degree-1 view
of `b` cannot flare (`L_g b = 0`); the HOCBF is what makes the constraint enforceable.

**Impact-load barrier (NACA TN 1516).** This is the attitude-coupled contact-load constraint
(the §3.3 idea, realized rigorously rather than via a hand-shaped `v_safe(θ)`). It bounds the
peak CG load factor a water contact at the current state would produce:
`n_peak = K0·ẏ₀²·Clf(κ)`, `κ = sinτ/sinγ₀·cos(τ+γ₀)`, with the load-factor coefficient `Clf(κ)`
precomputed offline (eqs 25/27, anchor `Clf(0)=0.6123`) and the hull coefficient
`K0 = (α_hull/(W g²))^{1/3}` (dead rise via `f(β)`, `φ(A)`). A height-relaxed term
`Φ(z)=Nb(1−e^{−z/zs})` makes it touchdown-only (Option A), backed by a height-scheduled slack
(Option C). It is **relative degree 2 via the elevator** (θ→q→δe) and 3 via thrust; a single
affine QP row can't mix the two (the elevator's degree-2 entry pushed to degree 3 gives
non-affine u²/u̇ terms), so it's a clean degree-2 row — thrust's column in `L_gL_f b` is exactly
zero and drops out, and the flare enforces it. Thrust still bounds impact load through the
descent/sink-rate barrier (they overlap by design). `K0` and a local-affine `Clf(κ)` are frozen
at the eval point (like `makeDescentBarrier`) so the templated barrier is smooth; the row is
assembled only in the model-valid window (below `z_gate`, descending, positive trim). Soft +
enabled + non-binding by default. Code: `include/autoland/impact_barrier.hpp`,
`scripts/precompute_impact_clf.py`; reference: `documentation/impact_load_barrier_spec.md` and
NACA TN 1516 (`documentation/19930082553.pdf`). Open: full mixed-degree (2/3) / predictive CBF
(see §9 and `TODO.md`).

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

- **`a_brk(V,γ)` is force-based and speed/path-angle-dependent** (since this branch):
  `a_brk(V,γ) = (ρV²S/2m)[C_L,max cosγ − C_D,maxlift sinγ] − g` — the max upward (lift − gravity,
  with the small drag correction) acceleration available at the current airspeed — substituted
  into the kinematic envelope `√(v_safe² + 2·a_brk·h)`. The envelope now tightens at low speed
  and loosens at high speed, as physics demands. Built by `makeDescentBarrier()` and
  differentiated **exactly** by the autodiff engine (the gradient picks up the new
  `∂b/∂V`, `∂b/∂γ` terms for free).
- **Thrust is deliberately omitted** (lift + drag + gravity only). Putting the thrust state `T`
  or pitch `θ` into `a_brk` would make `b` depend on them directly and drop the relative degree
  below 3, breaking the degree-3 HOCBF alignment (that θ-coupling is the §3.3 mixed-degree
  problem). Lift/drag depend only on `(V,γ)` — already in `b` — so degree 3 is preserved. The
  thrust contribution is a backlog item, coupled to the §3.3 work.
- **`C_L,max` is a config placeholder (1.2), `C_D,maxlift` is computed.** The inviscid aero has
  no stall ceiling, so `C_L,max` is an input (TODO: calibrate). `C_D,maxlift` is the rotated drag
  at the linearly-extrapolated max-lift α (a small, ∝sinγ correction).
- **Conservatism / positivity.** `a_brk(V,γ)` overestimates *instantly* recoverable braking
  (it ignores the lift-build-up lag — the relative-degree caveat), and goes negative below
  ~stall (lift < weight), which would make the radicand ill-posed. The airspeed barrier keeps
  `V` above stall so `a_brk > 0` in the operating envelope; the filter warns once if not. A
  lag-aware / committed-flare backup set is the proper conservative treatment (backlog).
- **Result on the real deck.** With `AHAB_combined.stab`, `C_L,max = 1.2`, `c_descent = [2,2,2]`:
  touchdown sink **0.011 m/s** (within the 0.1 budget), **0 QP recoveries**, and the descent
  barrier `b ≥ 0` for the whole flight (min ≈ 0.10). The small residual descent `ψ < 0` is **not**
  the initial condition (the IC is now a true level-flight trim, §8.4) and **not** control
  saturation (0 recoveries; smaller class-K gains make it *worse*): it is a **terminal** effect,
  confined to the last ~0.4 s as `h → 0` at touchdown (the discrete-time crossing). The realized
  trajectory stays safe — `b`, the actual criterion, never goes negative.
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
- `data/AHAB_combined.stab`: real combined VSPAero aero deck. This is now the **default** deck
  for `lon_autoland_sim` (the old `example.stab` placeholder is kept only for the unit tests;
  removing it is a backlog item). The deck is **inviscid** (lift linear to ±20°, no stall); an
  optional NACA 4414 viscous-stall overlay can be enabled on top of it — see below and
  `stall_model_spec.md`.
- `data/lon_scenario.yaml`: powered approach `T_set = 2.0 N`, `V_app = 18`, `γ = −3°`,
  `v_safe = 0.1`, `CL_max = 1.2` (placeholder, drives `a_brk(V,γ)`), `Vmin = 13.5`,
  `V_max = 27` (placeholder over-speed ceiling, ~1.5·V_app; non-binding on this approach), `Tmax = 50`,
  descent class-K `[2,2,2]`; gains retuned for the real mass. (`a_brk = 3.0` remains in the file
  but is **deprecated/unused** — kept for logging/plot-script compatibility.)

## 8. Assumptions & limitations

Folded from the 2026-06-25 implementation notes (`archive/`), updated for current state:

1. **Self-consistent plant.** The sim plant *is* the CBF model (zero model mismatch), so clean
   barrier-invariance partly reflects structure; it validates the math, not robustness to model
   error.
2. **`L_f³b` drift cross-checked (2026-06-27).** A test-only finite-difference *flow oracle*
   now verifies the full drift stack `{b, L_f b, L_f²b, L_f³b}` for the descent and airspeed
   barriers against the Taylor-jet engine (agreement ~2.5e-4 on `L_f³b`, finite-difference-
   limited; a 1% engine error is caught). This closes the prior gap where only the control
   authorities were checked. See `test/test_lon_cbf.cpp` ("drift Lie stack matches a
   finite-difference flow oracle").
3. **Frozen local-affine aero** — autodiff exactness is relative to a C0 piecewise-linear table
   surrogate; aero curvature is lost in 2nd/3rd derivatives.
4. **Initial condition is a steady level-flight trim of the lon model** (`lonTrim()` in
   `src/lon_sim.cpp` — a Newton solve on the lon EOM for `γ=0`, `V̇=γ̇=q̇=0`). The run starts at a
   true equilibrium (cruising level), then the nominal commands `γ_ref` and the aircraft pushes
   over into the approach — no more body-axis-seed startup transient / airspeed drift. (The
   nominal's *feedforward* `θ_trim`/`T_set` still come from the body-axis approach trim; that
   only feeds a controller with feedback, so the small model mismatch there is harmless.)
5. **No airspeed control in the nominal** (constant-thrust directive) — airspeed leans on the
   airspeed CBF; a high `Vmin` parks the energy-limited glider in a safe hover instead of landing.
6. **Soft airspeed barriers** (slack-penalized) — not a hard guarantee. Both the stall
   (`V ≥ V_min`) and over-speed (`V ≤ V_max`) barriers are soft; descent + thrust are hard. The
   over-speed barrier reuses the stall barrier's degree-3 machinery with flipped authority signs;
   `V_max` is a placeholder set generously (non-binding on the nominal approach) — protective
   only until calibrated to a real never-exceed / structural speed.
7. **Discrete-time / ZOH** at dt = 0.01 approximates the continuous-time CBF guarantee.
8. **No controllability guard** on `L_g L_f²b` (elevator authority → 0 at very low dynamic
   pressure would make the descent barrier unenforceable by elevator; not yet observed).
9. **Idealized scope:** longitudinal-only, perfect full-state feedback, no wind/gust/sensor
   noise, instantaneous elevator (no actuator lag).
10. **§3.3 deferred ⇒ "hull-safe" is incomplete.** `v_safe` is a flat constant and touchdown
    attitude is unconstrained; real slamming load is attitude-dependent.
11. **Optional viscous stall (plant).** An optional NACA 4414 stall overlay (NeuralFoil-derived,
    `stall_model_spec.md`) gives the otherwise-stall-free inviscid deck a real post-stall lift
    drop / drag rise / pitch break. OFF by default (nominal behaviour bit-identical). Its onset
    and depth are **tunable engineering values, not experiment-validated** (NeuralFoil→Xfoil→
    reality is weakest at low-Re post-stall), and it is frozen local-affine like the rest of the
    aero (item 3). It is shared into the CBF path, so the self-consistent-plant property (item 1)
    still holds when it is enabled. Designing a stall-recovery barrier on top is future work.

## 9. Open questions / follow-ups

See `TODO.md` for the full backlog. The major threads:

- **Thrust term in `a_brk`** — needs the mixed-relative-degree treatment (§6); couples with §3.3.
- **Lag-aware / conservative `a_brk`** — a committed-flare backup set so the set is truly
  invariant despite the lift-build-up lag (§6).
- **Calibrate `C_L,max`** — currently a placeholder 1.2.
- **§3.3 contact-force barrier** — attitude-coupled `v_safe(θ)` from von Kármán/Wagner slamming
  theory; the contribution most likely to make this a paper.
- **True longitudinal trim** — Newton solve on the lon model for a consistent initial condition
  (also clears the residual descent `ψ < 0` start transient).
- **Hardware path** — PX4 insertion point (actuator vs rate/attitude setpoint sets the relative
  degree).

---

*See `CHANGELOG.md` for the running log of changes.*
