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
flies a constant-glide-slope approach; a CBF-QP minimally adjusts its commands to keep the
aircraft in a provably safe set, passing the nominal through untouched until a constraint is about
to bind. The current **"recovery" barrier set** guarantees a bounded hydrodynamic touchdown load
(**hard**) and softly enforces angle-of-attack (stall) margin, a nose-up touchdown attitude, and a
bounded touchdown energy/speed. The AoA guard reproduces the pilot **low-altitude stall recovery** —
pitch down + full elevator as α approaches stall. The flare is **not** designed into the nominal;
near the surface it emerges from the impact-load barrier (and, once the nominal moves to an
exponential-altitude flare, from the nominal itself).

Scope is **longitudinal (pitch-plane) only**. Lateral/heading dynamics, wave-phase exploitation,
and the post-touchdown on-water regime are out of scope (see `archive/water_landing_cbf.md`).

This is a separate, self-contained `lon_*` module living alongside the original 11-state
body-axis code (`dynamics.cpp`, `cbf.cpp`, `sim.cpp`), which it leaves untouched.

## 2. Augmented longitudinal model

- **State** `X = [h, V, γ, θ, q, T, Ṫ]` (NXA = 7); **control** `U = [δe, T̈]` (NUA = 2).
- Thrust is **integrator-augmented**: `T` and `Ṫ` become states, `T̈` the control. This aligns
  the relative degree of thrust with that of elevator for barriers that involve airspeed (both
  **degree 3**, e.g. the total-energy barrier), so a single well-posed QP can share authority
  between them. The elevator-only barriers (stall, nose-up, impact) are **degree 2**.
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

This replaces the math spec's TECS placeholder by directive. The point is that the nominal only has
to put the aircraft on a reasonable glide slope; the CBF owns the safety-critical shaping. **The
nominal is slated to move to an exponential-altitude flare `h(t) = h₀ e^{−t/τ}`** (sink ∝ height);
because every barrier uses only the *instantaneous* state, that switch needs **no barrier changes**.
Code: `include/autoland/lon_nominal.hpp`.

## 4. Barriers

| Barrier | Definition | Rel. degree | Enforcement |
|---|---|---|---|
| Impact-load `b_imp` | `(n_limit − n_peak(θ,γ,V)) + Φ(z)` | **2 (δe)** | **Hard** — the only hard *safety* row |
| Stall / AoA `b_stall` | `α_lim − (θ − γ)`, `α_lim = α_stall − margin` | 2 (δe) | Soft, weight 1e5 |
| Total-energy `b_E` | `½(V_td_max² − V²) + (g_eff − g)·h` | 3 (δe & T̈) | Soft, weight 1e4 |
| Nose-up `b_nose` | `θ − θ_min` (active `h < h_noseup`) | 2 (δe) | Soft, weight 1e3 |
| Min thrust | `T ≥ 0` | 2 (T̈) | Hard (closed form) |
| Max thrust | `T_max − T ≥ 0` | 2 (T̈) | Hard (closed form) |

Priority is set by the slack weights: **impact (hard) > stall (1e5) > energy (1e4) > nose-up
(1e3)**, so under a conflict on the shared elevator the QP sacrifices nose-up first, then energy,
and protects the stall guard. There is deliberately **no airspeed-floor barrier**. The HOCBF rows
are built generically from the exact Lie stack via `hocbfRow()`; the thrust actuator barriers use
closed forms. Code: `include/autoland/hocbf.hpp`, `src/lon_cbf_filter.cpp`. Derivations and
control-affine authority terms: math spec §3.

**The recovery emerges from the stall barrier.** `b_stall = α_lim − (θ−γ)` is degree 2 via the
elevator; as α → α_lim the class-K cascade drives δe **nose-down**, saturating at the boundary —
i.e. *pitch down + full elevator*, exactly the pilot low-altitude stall recovery. This is emergent,
not scripted. `α_stall` comes from the NACA 4414 plant stall model (11°).

**The nose-up barrier keeps the impact model valid.** `b_nose = θ − θ_min` (θ_min ≥ θ_keel) holds a
nose-up attitude in the final metres so the impact model's gate `τ = θ − θ_keel > 0` stays valid and
the touchdown is nose-up. It is **θ-based, not AoA-based** (a θ floor directly guarantees the gate;
an α floor would not, since a steep γ lets θ = α + γ go negative). As the lowest-weight row it
**yields to the stall guard**: when a steep γ makes the demanded nose-up conflict with the stall
limit, stall wins and the nose drops — which is the recovery.

**The energy barrier is a loose never-exceed `(V,h)` ceiling.** `b_E` bounds touchdown speed via a
height-scheduled total-specific-energy cap `E ≤ ½V_td_max² + g_eff·h`, which reduces to a descending
airspeed cap `V ≤ √(V_td_max² + 2(g_eff−g)h)`. A *constant* energy cap is useless (potential energy
converts to kinetic on descent, so a constant-E glider arrives fast); the height schedule fixes
that. `V_td_max` is the hull/structural touchdown limit; `g_eff` is sized so the ceiling clears the
whole approach and starts satisfied (`g_eff ≥ g + (V₀²−V_td_max²)/(2h₀)`). `b_E` is a polynomial in
`(V,h)` — no √, no floors — and touches only `V,h`, preserving the degree-3 alignment.

**Impact-load barrier (NACA TN 1516).** The attitude-coupled contact-load constraint, realized
rigorously rather than via a hand-shaped `v_safe(θ)`. It bounds the peak CG load factor a water
contact at the current state would produce: `n_peak = K0·ẏ₀²·Clf(κ)`,
`κ = sinτ/sinγ₀·cos(τ+γ₀)`, with the load-factor coefficient `Clf(κ)` precomputed offline
(eqs 25/27, anchor `Clf(0)=0.6123`) and the hull coefficient `K0 = (α_hull/(W g²))^{1/3}` (dead rise
via `f(β)`, `φ(A)`). A height-relaxed term `Φ(z)=Nb(1−e^{−z/zs})` makes it touchdown-only (Option A).
It is **relative degree 2 via the elevator** (θ→q→δe) and 3 via thrust; a single affine QP row can't
mix the two (the elevator's degree-2 entry pushed to degree 3 gives non-affine u²/u̇ terms), so it's
a clean degree-2 row — thrust's column in `L_gL_f b` is exactly zero and drops out, and the flare
enforces it. It is now the **only hard safety row** (`impact_hard`); thrust bounds the *approach*
energy/speed through the energy barrier instead. `K0` and a local-affine `Clf(κ)` are frozen at the
eval point so the templated barrier is smooth; the row is assembled only in the model-valid window
(below `z_gate`, descending, positive trim). Code: `include/autoland/impact_barrier.hpp`,
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
out true aero curvature. Exact math ≠ exact physics. (The stall, nose-up, and energy barriers are
polynomial/linear in the state, so they are exact regardless.)

## 6. Key design decisions

- **Impact is the only hard *safety* barrier; the thrust CBFs are the actuator-effectiveness /
  HOCBF-validity guards.** Every elevator-driven barrier's control coefficient scales with dynamic
  pressure (`L_gL_f^{r-1}b ∝ ρV²`); keeping the augmented thrust chain bounded (and the aircraft
  flying) preserves the relative-degree structure the whole formulation rests on. This is why there
  is **no separate airspeed-floor CBF** — stall is handled directly in AoA.
- **Soft-row priority via slack weights** `w_stall (1e5) > w_energy (1e4) > w_noseup (1e3)`. This
  encodes impact > stall > energy > nose-up so the recovery is protected while nose-up yields first.
- **Stall is AoA-direct, not an airspeed proxy.** `b_stall = α_lim − (θ−γ)` with `α_lim = α_stall −
  margin`. `α_stall` (11°) comes from the NACA 4414 plant stall model, so the barrier and the plant
  agree on where stall is (the old airspeed floor conflated speed with AoA). The pitch-down recovery
  is the emergent behavior when α → α_lim.
- **Nose-up is θ-based and gated to the final metres.** It exists to keep the impact model's gate
  (`τ > 0`) valid and give a nose-up touchdown, and is the lowest-priority row so it yields to stall
  (the recovery). `θ_min ≥ θ_keel`.
- **Energy is a loose never-exceed ceiling, not a tracked schedule.** It bounds touchdown speed but
  is sized to clear a good approach; it binds only on a hot/steep arrival. **Operating rule:**
  keep `V_td_max ≳ 1.1·V_stall` (~11 m/s) for a soft landing — demanding a sub-stall touchdown speed
  is physically un-flareable and forces a hard nose-down landing (the honest insufficient-margin
  case; see below).
- **Glide-path independence.** All barriers are functions of the instantaneous state only, so the
  set does not assume the constant-γ nominal and needs no changes when the exponential-altitude
  flare lands (§3).
- **Doc-vs-autodiff divergence (recorded as a test):** the math spec's closed-form elevator
  authority omits the pitch-rate aero terms (`∂C_L/∂q̂`, `∂C_D/∂q̂`); the autodiff captures them.
  The autodiff quantity is the more complete one.
- **Touchdown sink/speed are tuning consequences, not barrier properties.** A safety *filter* yields
  the constraint with slack; the exact touchdown values follow from the class-K gains and (once it
  lands) the nominal flare, with the CBFs as backstop.
- **No `h ≥ 0` CBF.** Any class-K has `α(0)=0`, so an `h ≥ 0` barrier forces `sink → 0` at the
  surface — it would forbid touchdown. Vertical rate at contact is bounded by the impact-load
  barrier (load factor), which is the ground-proximity condition in the only form compatible with
  landing.
- **Result on the real deck** (`AHAB_combined.stab`, α_lim = 9°): the default approach lands
  **soft nose-up** (sink ≈ 0.21 m/s, θ ≈ +3°, V ≈ 13.4 m/s), impact (hard) satisfied throughout,
  0 hard rows dropped. A tight-but-feasible touchdown (`V_td_max = 10`) lands even softer/slower
  (≈ 0.03 m/s at V ≈ 11, θ ≈ +7.3°) with the stall guard protective (α peaks ~7.4°). Only when a
  **sub-stall** touchdown speed is demanded (`V_td_max = 9 < V_stall ≈ 9.7`) does the guard bind at
  9° and take a hard nose-down landing (~2.7 m/s) rather than stall — the honest recovery tradeoff.

## 7. Vehicle & scenario data (current)

The model uses **real AHAB vehicle data**, not placeholders:

- `data/aircraft.yaml`: mass 3.6139 kg; Iyy 0.0521 kg·m²; real Ixx/Izz; thrust `T_static = 50 N`,
  `zcp = 0.15 m`; `parasite_CD0 = 0.030`. (Thrust `k_v` and `parasite_CD0` still uncalibrated.)
- `data/AHAB_combined.stab`: real combined VSPAero aero deck, the **default** for `lon_autoland_sim`
  (`example.stab` is kept only for the unit tests). The deck is **inviscid** (lift linear to ±20°,
  no stall); an optional NACA 4414 viscous-stall overlay can be enabled on top of it — see
  `stall_model_spec.md`.
- `data/lon_scenario.yaml`: powered approach `T_set = 2.0 N`, `V_app = 18`, `γ = −3°`, and the
  recovery-barrier knobs: `alpha_stall_deg = 11` / `stall_margin_deg = 2` (⇒ α_lim = 9°),
  `theta_min_deg = 3` active below `h_noseup = 3 m`, `V_td_max = 14` and `g_eff = 16` (the energy
  ceiling), `impact_hard: true`. Class-K: stall `[4,4]` / nose-up `[2,2]` (deg 2), energy `[2,2,2]`
  (deg 3), impact `[2,2]`. Slack weights `1e5 / 1e4 / 1e3` (stall / energy / nose-up). Impact
  hydrodynamic knobs (`n_limit`, `beta_deg`, `Nb`, `zs`, …) remain placeholders pending calibration.

## 8. Assumptions & limitations

Folded from the 2026-06-25 implementation notes (`archive/`), updated for current state:

1. **Self-consistent plant.** The sim plant *is* the CBF model (zero model mismatch), so clean
   barrier-invariance partly reflects structure; it validates the math, not robustness to model
   error.
2. **Drift Lie stacks cross-checked.** A test-only finite-difference *flow oracle* verifies the
   drift stack `{b, …, L_f^R b}` for the stall/nose-up (deg 2) and energy (deg 3) barriers against
   the Taylor-jet engine, alongside closed-form and hard-enforcement tests for all three new
   barriers and the impact barrier. See `test/test_lon_cbf.cpp`.
3. **Frozen local-affine aero** — autodiff exactness is relative to a C0 piecewise-linear table
   surrogate; aero curvature is lost in 2nd/3rd derivatives. (The stall/nose-up/energy barriers are
   polynomial in the state, so they are unaffected.)
4. **Initial condition is a steady level-flight trim of the lon model** (`lonTrim()` in
   `src/lon_sim.cpp`). The run starts at a true equilibrium (cruising level), then the nominal
   commands `γ_ref` and the aircraft pushes over into the approach. (The nominal's *feedforward*
   `θ_trim`/`T_set` still come from the body-axis approach trim; harmless, it only feeds feedback.)
5. **No airspeed control in the nominal, and no airspeed-floor CBF.** Airspeed is bounded *above*
   by the energy ceiling; stall is guarded directly in AoA. There is no lower-airspeed floor, so on
   an energy-limited approach the aircraft can get slow — the **stall guard** (not an airspeed floor)
   is what prevents an actual stall, taking a hard nose-down landing if the touchdown demand is
   un-flareable (§6).
6. **Only impact is a hard safety guarantee.** Stall, nose-up, and energy are soft (slack-penalized)
   with priority weights `1e5 > 1e4 > 1e3`; the thrust actuator barriers are hard. So the soft rows
   can transiently dip (e.g. the energy ceiling on a hot segment) without violating the hard set.
7. **Discrete-time / ZOH** at dt = 0.01 approximates the continuous-time CBF guarantee.
8. **Elevator authority ∝ ρV².** The thrust CBFs keep the aircraft flying (dynamic pressure up), so
   the elevator-driven rows' control coefficient stays away from zero — the actuator-effectiveness
   role (§6). There is still no *explicit* guard that trips if elevator authority collapses; not yet
   observed in the operating envelope.
9. **Idealized scope:** longitudinal-only, perfect full-state feedback, no wind/gust; an
   altitude-sensor noise + low-pass model is available (`h_meas_stddev`, `h_lpf_tau`) but off by
   default; instantaneous elevator (no actuator lag).
10. **Hull-safe = impact-load + nose-up attitude.** The impact-load barrier (hard) bounds the
    attitude-dependent slamming load at contact, and the nose-up barrier keeps a valid nose-up
    touchdown attitude. Vertical rate at contact is bounded only through the load factor (no separate
    `v_safe` row); horizontal/total speed is bounded by the energy ceiling (soft).
11. **Optional viscous stall (plant).** The NACA 4414 stall overlay (NeuralFoil-derived,
    `stall_model_spec.md`) gives the inviscid deck a real post-stall lift drop / drag rise / pitch
    break. OFF by default. Its onset/depth are **tunable engineering values, not experiment-
    validated**, and it is frozen local-affine like the rest of the aero (item 3). The stall
    barrier's `α_stall` is set to the overlay's onset so barrier and plant agree.

## 9. Open questions / follow-ups

See `TODO.md` for the full backlog. The major threads:

- **Exponential-altitude nominal flare** — swap the constant-γ cascade for `h(t)=h₀e^{−t/τ}`; the
  barriers are already glide-path agnostic, so this is a nominal-only change (retune
  `V_td_max`/`g_eff`/`theta_min` after).
- **Calibrate the impact hydrodynamic knobs** — `n_limit`, `beta_deg`, `Nb`, `zs`, `tau_keel` are
  placeholders; and validate the stall overlay's onset/depth.
- **Harder energy/speed guarantee** — the energy ceiling is soft; if a firm touchdown-speed bound is
  needed, consider a hard energy row (feasibility permitting) or fold it into the nominal.
- **Full mixed-degree (2/3) impact CBF / predictive form** — the impact row is a clean degree-2
  approximation; the full 2/3 treatment is backlog.
- **Update `water_landing_cbf_design.md`'s sibling docs** as the nominal changes; keep the math spec
  (§3/§4) in sync.
- **Hardware path** — PX4 insertion point (actuator vs rate/attitude setpoint sets the relative
  degree).

---

*See `CHANGELOG.md` for the running log of changes.*
