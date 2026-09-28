# 6-DOF Corridor Water Landing with Terrain Keep-Out: Roadmap

*Implementation roadmap agreed with the author on 2026-09-26 (plan-mode interview). **Status
(2026-09-26): Phase 0 DONE, Phase 1 DONE, Phase 2 IN PROGRESS with a pinned negative result**
(see the "as built" sections below); Phases 3-5 not started. This doc is the living roadmap for the corridor / keep-out / go-around effort; each
phase updates its section as it lands and logs the session in `CHANGELOG.md`.*

---

## Phase 0 as built (2026-09-26, revised after controls review the same day)

Everything in the Phase 0 plan below landed, with these deviations / lessons:

* **State:** `XN = 11` appended (`types.hpp`), north row in both plants; the sim's
  trapezoid `x_pos` is gone. Legacy scenarios are identical to the committed baseline in
  every printed column but `x` (RK4 vs trapezoid, <= 4e-7 m), checked by diffing the CSVs
  of a HEAD build; the outer-loop refactor (`stepOuter` = `stepOuterRef` at the config
  references) is pinned bit-identical by a 500-step random-sequence test for both
  longitudinal modes. `validate_beaver_sixdof.py` checks 12 rows.
* **Geometry:** `world_geometry.hpp` — `Corridor` (frame/inverse/contains/axisPoint),
  `TerrainField` (super-Gaussian bumps + `addRing` with an optional gap; templated
  `height<T>`, exact `gradient` via dual), `RolloutConfig`, `WorldConfig`. Verified
  against a hand chain-rule gradient and under the Taylor jet (`lieDrift<2>` of the
  field along a constant flow = grad·v, v'Hv).
* **Rollout:** `water_rollout.hpp` (`runRollout`, double-only RK4; the templated RHS was
  not needed by anything). Result nested as `SixDofTouchdown::rollout`; rollout rows
  continue the CSV with `phase = 1`. Tests pin the kq = 0, a0 = 0 and combined
  (`d = ln((a0 + kq V0²)/a0)/(2 kq)`) closed forms. `a0 = 0.8`, `kq = 0.002` are
  placeholders. The rollout knows nothing about the crab at contact (float side-loads,
  water-looping) — "stayed inside" is a planar statement only.
* **Guidance:** `pattern_guidance.{hpp,cpp}`. What the plan got wrong, and what the
  review caught:
  1. **FAF placement is a 25 m grid scan, not a fixed point / bisection.** Dubins length
     is discontinuous in the goal pose (a nearly aligned FAF just ahead needs a 2πR loop:
     the trim sideslip alone puts the initial ground course 0.8° off the axis), so
     `d + L(d)` is neither continuous nor monotone.
  2. **The straight-in anchor needed a DIRECT-line candidate.** The first version got a
     near-straight plan only because the inner-tangent word existed on a knife edge
     (centre distance minus 2R = 0.11 m). Now: when the aircraft is within 30 m / 5° of
     the axis ahead of the FAF the direct line competes (and wins by length); CSC
     candidates with a near-full-circle arc (> 350°) are rejected. A test perturbs the
     start by ±5 m / ±1° and requires the DIRECT plan every time.
  3. **The closest-point search needs a forward window too** (`s_ahead = 500 m`, not only
     `s_back`): an overhead entry sits ON the corridor axis at t = 0 and the search
     snapped to the final leg, skipping the pattern. A self-crossing-path test pins it.
  4. **The lateral law is arclength-carrot pursuit with the Park–Deyst–How gain**, not
     PDH's circle-intersection reference point; the `L1 = scale·period·V/π`
     parameterization is ArduPilot's. Linearized about a straight path
     ω_n = √2·V/L1, ζ = 1/√2 regardless of the parameters (they only size L1: 239 m at
     40 m/s, natural period 26.5 s). Using the actual chord to the carrot makes the
     commanded acceleration exactly V²/R on any circular segment. Only CSC words are
     searched: "shortest of these", never "Dubins shortest path" (CCC can win when the
     same-sense circle centres are within 4R, which the scan's FAF goals often are;
     minimality is not the objective).
  5. **Glideslope cone:** with the horizontal ground speed in the denominator the height
     error obeys ė = −K_h e exactly (τ = 10 s). The cone's ±[−8°, +6°] clamp is followed
     by the cascade's own speed shift `Kv_gamma` (τ ≈ 5 s, ±4°), so `gamma_ref_eff` can
     reach −12° and the two loops sit within a factor of two in bandwidth; the runs
     stay within γ ∈ [−4.5°, −3.0°], a hot-and-high entry would exercise it.
* **Turn coordination was missing — twice.** The legacy yaw damper `-Kr r` opposes every
  steady turn: at 22° bank the Beaver skidded round at 3°/s (coordinated 5.8°/s) with
  560 m cross-track. A turn-rate feedforward (`r_coord = g sin φ cos θ / V`, the exact
  steady-turn body yaw rate under μ ≈ φ, β = 0) fixed the rate but still left |β| = 11°
  at the roll-in (adverse yaw, rudder parked at trim). `stepInnerPhi` now also closes a
  sideslip PI on the rudder (same loop as the reach layer's `stepInnerBeta`,
  β_cmd = trim sideslip). `stepInner` (legacy) is untouched. Whether LR-556 is valid at
  the |β| ≈ 6° still reached in the roll-in is **not stated in the source** — see TODO.
* **Results** (all: touchdown in corridor, rollout stopped inside, never left it; no
  flare / decrab / level-off exists, so bank at contact is simply the converged tracker):

  | Scenario | s_td [m] (aim −200) | e_td [m] | course err [deg] | crab [deg] | max xte to path [m] | max \|e\| on final [m] | max \|φ\| / limit-steps | max \|β\| [deg] | CG clearance over trees [m] |
  |---|---|---|---|---|---|---|---|---|---|
  | `beaver_corridor_straight` | −199.6 | 0.00 | 0.00 | +0.79 (trim) | 1.3 | 1.3 | 1.1° / 0 | 0.9 | 9.1 |
  | `beaver_corridor_overhead` | −198.9 | 0.02 | −0.01 | +0.79 | 9.4 | 4.8 | 24.0° / 0 | 4.5 | 9.6 |
  | `beaver_corridor_crosswind` (5 m/s east) | −190.0 | 0.07 | −0.02 | −6.4 | 12.2 | 5.5 | 24.6° / 0 | 6.3 | 10.4 |

  Before the sideslip loop the overhead case had 35 m max cross-track, e_td = −1.5 m,
  |β| = 11°, and both pattern cases spent 10–20 % of the flight at the 25° bank limit
  (crosswind: 100 m cross-track). Figures: `figures/corridor_{straight,overhead,crosswind}.png`
  (`scripts/plot_corridor_landing.py`); interactive replay: `scripts/build_corridor_replay.py`.
  Test bounds are stated requirements (course within 3°, |e| within a quarter half-width,
  wings within 3° at contact, bank command never at the limit, LR-556 speed band, stop
  point 50 m short of the far edge), not multiples of the observed values.
* **Not done / deferred (TODO):** replan-on-large-error, terrain-aware planner, wind-aware
  (air-mass) Dubins or a continuous-curvature planner, LR-556 sideslip validity band,
  any flare / decrab / wings-level gate, TECS + corridor (untested).

---

## Phase 1 as built (2026-09-26): 6-DOF surfaces-only CBF engine + envelope rows

Built to the amended plan (no engine or servo lag, no turbulence or sensor models by
decision; those are stated limitations). `include/autoland/sixdof_cbf.hpp`,
`src/sixdof_cbf.cpp`, `test/test_sixdof_cbf.cpp`, scenario block `sixdof_cbf:`.

**Model.** Decision variables are the three surfaces; the throttle stays with the nominal and
enters the barrier model as a known constant, so every row has uniform relative degree 2 in
the surfaces with no added state. Split at each step: `xdot = f + g U + F (U − u_prev)` with
`g` the moment columns (exact: all moment terms are affine in the surfaces), `f` the plant
with forces frozen at the last applied surfaces, `F` the direct control-force columns left
out of the control map. Lie derivatives are Taylor jets of the plant polynomials themselves
(`BeaverDynamics::xdotCoreT` is scalar-generic; `lie_taylor.hpp` gained `tan`, `atan`, `asin`
and scalar division); the jets are pinned to the plant derivative and to the autodiff
linearization, and `f + g U` to the plant with wind, to 1e-9.

**Wind.** The aerodynamics are air-relative: the model takes an earth-frame wind estimate and
rotates it into the body frame inside the jet. Without it a crosswind crab read as sideslip and
the soft sideslip row flew the crosswind case into the trees. In the sim the estimate is the
true wind (idealization: a perfect GPS-minus-air-data wind estimate); the nominal still sees
only pitot air data.

**Robustness margin (the reviewer's model-mismatch point).** The neglected force term enters
the degree-2 condition at BOTH orders: `hdot_true = L_f h + (L_F h) dU`,
`hddot_true = … + (L_F L_f h + L_f L_F h) dU`. Each row subtracts
`Σ_k (|D2_k| + (c1 + c2)|D1_k|)·min(rate·dt, range_k)` with `D1 = L_F h`,
`D2 = L_F L_f h + L_f L_F h`, both exact: the first-order-in-ε coefficients of the Lie stack of
the perturbed drift `f + ε F_k`, ε seeded as an autodiff dual through the jet. (The first draft
used `|∇(L_f h)·F + c1 ∇h·F|` — wrong by the flow derivative of the first-order leak and by
`c2`; the hand-derivation test caught it.) A per-row constant `margin_extra` covers the
sampled-data gap. The claim everywhere: invariance of the modelled system with margin.

**Anti-windup through the filter.** `SixDofNominal::commitApplied(u)`: any axis the filter
moved has its integrator increment undone that step (the nominal's own conditional
integration), and the rate limiter's reference becomes the applied command.

**QP.** `[de, da, dr, slacks]`; rows normalized by `max|a|`; input box = deflection AND
one-step rate limit about the applied command; hard rows without slack; soft rows with
quadratic slack; best-effort minimum-violation fallback when the hard set is infeasible.
Lesson: the best-effort hard-slack penalty must be O(1e2), not 1e4-1e6 — OSQP (eps 1e-8) hit
its iteration limit on the violated-row problem and the filter silently returned the nominal
(the slow case showed "active 1/12722 steps"). The solver status is now in the diagnostics.

**Envelope limits sit above the nominal's clamps** (bank 30 vs the tracker's 25; AoA 13 vs the
band; V 30-55) so the rows are rare interventions by construction; a case that wants a row to
bind raises the tracker's clamp instead (tight pattern: tracker 35, barrier 25).

**Measured (tests):**
* Oracles: bank row vs an independent kinematic derivation; AoA first order by hand, second
  order vs a flow difference; margin structure vs hand `D1`, `D2`.
* Adversarial invariance, hostile nominal at full deflection against each hard row, 15 s at
  10 ms AND 2.5 ms: bank ≥ −0.1°, AoA ≥ −0.05°, sideslip ≥ −0.05°, no best-effort.
* Feasibility sweep over 200 random in-set states (both HOCBF levels): each family alone with
  the deflection box 200/200; **all four hard at once 200/200 — the rows are compatible**;
  with the one-step rate box: V 200, bank 199, sideslip 198, **AoA 155/200** — a rate-limited
  surface cannot always meet a hard row within one 10 ms sample from an arbitrary in-set
  state (it can over a few); hence one hard row per scenario and best-effort counted.

**Results:**

| Case | Filter | Outcome |
|---|---|---|
| overhead / crosswind / straight | soft rows, limits above clamps | LANDED OK, filter active ≤ 1 step, best-effort 0 |
| `cbf_tight_pattern` (tracker allowed 35°) | bank row HARD at 25° | flown bank 25.00° (min h −0.00006°), active 2742/10577 steps, best-effort 0, LANDED OK; nominal alone: 26.4° |
| `cbf_slow` (36 → 31 m/s deceleration) | AoA row HARD at 13° | α max 12.94°, stall never entered, V min 34.8; the aircraft then runs out of altitude mid-pattern (an approach at 31 m/s cannot be flown inside the band) — envelope kept, corridor not claimed; nominal alone: 18.4°, stall region |
| `fail_tailwind` (new) | none | airspeed on target, 8 m/s tailwind, water run 477 m, stops 109 m past the far edge — Phase 3 |

Compute: mean 0.47 ms per step for 7 rows (3 jets of order 2 over 12 states, 3 g-columns,
3 ε-perturbed jets, one OSQP setup+solve), worst step ~10 ms (OSQP setup outliers); a
warm-started solver would remove the outliers if this ever runs at rate.

**Stall region = failed run.** `envelope.alpha_stall_deg` (16): any step above it flags the run
FAILED whatever else happens; the console, the verdict (`SixDofTouchdown::success` /
`fail_reasons`) and the replay chips say so. Also added: wingtip clearance
(CG clearance − (b/2)|sin φ|), `terrainRowMargin()` / `terrainMaxGradient()` for the Phase 2
row (map error + position error through the slope bound + wingtip at the bank limit), and
`scripts/calibrate_rollout.py` (fits `a0`, `kq` to POH landing-run points — the POH numbers
still have to be supplied; the shipped values remain placeholders).

**Not done (by decision or deferred):** engine/servo lag, Dryden turbulence, sensor models;
cross-plant (AHAB) check waits for the full row set; TECS + corridor untested.

## Phase 1 addendum: the first case a barrier row actually flips

`fail_gust_turn` / `cbf_gust_turn(_off)`: an 18 m/s crosswind gust arriving in the first turn
of the overhead pattern (t = 22 s, 25° bank). The sideslip transient feeds the Beaver's
β² pitch-up term (`Cm_b2 = +0.69`) and α spikes to 17.4° — the nominal is fine before and
after, but the run enters the stall region. With the hard AoA row at 13.5°: α peaks at
13.4°, 635 active steps, no best-effort, **LANDED OK** with the touchdown unchanged. Negative
results from the same sweep: vertical gusts do not produce the excursion (the pitch loop
follows a 2 s ramp; a 0.75 s ramp changes the wind faster than any surface can answer, and
the discrete-gust model is a persistent step, so both nominal and filter fail); no lateral
gust produced a bank excursion the roll loop could not handle (max 26.9°). Earlier "the CBF
improved nothing" tally (tight pattern: envelope only; slow: not a landing) stands.

## Phase 2 as built (2026-09-26, in progress): terrain row — works straight ahead, does NOT steer through a gap

Built: the terrain keep-out row on the 6-DOF filter, `TerrainBrakingRow` in
`src/sixdof_cbf.cpp`, YAML `row_terrain` + `terrain_*` keys, degree-3 generalization of the
row builder (`pushR<R>`, margin `Σ coeff_j |D_j|·dU` — the pure clearance row is kept as
`terrain_degree3`), closed-form scalar-generic terrain gradient (`TerrainField::gradientT`).

Four formulations were tried, in order, each fixing the previous one's failure:
1. **Pure clearance, degree 3, linear class-K** — far too permissive far from the ground (at
   57 m it tolerates a 57 m/s sink), pull-up comes too late for a 1–2 g airplane.
2. **Braking form** `ż + sqrt(v² + 2 a z⁺) + k z⁻ ≥ 0` (the lon filter's descent envelope with the
   terrain slope; the `k z⁻` term is needed because the lon form is built to ARRIVE at z = 0
   and otherwise keeps sinking through the margin) — myopic: it lets the aircraft settle at the
   margin on the flat and then meet a slope no climb can fly (closure from a slope is not
   arrested by a pull-up).
3. **Climb-capability lookahead**: h_T replaced by max over points ahead of
   `h_T(s) − s·climb_grad`, the height needed now to clear the terrain at s on a sustainable
   climb — with the lookahead direction frozen along the course to keep degree 2. Works for
   terrain straight ahead: against a full nose-down hostile nominal at climb power, CG
   clearance never drops below 2.8 m over a 3 m margin at 10 ms and 2.5 ms
   (`[terrain][invariance]`). But the barrier quantity itself goes to −43 and the row runs
   best-effort (the capability form is conservative); the test claims the clearance, not b ≥ 0.
4. **Fan of directions** about the course (a single frozen ray swept across a peak under yaw
   made b jump by hundreds of metres between steps; the hard row flung the aircraft into a
   stall). Continuous, but:

**The gap (negative result, pinned by `[baseline]`).** `fail_gap_final`: two 350 m mountains
1 km short of the aim point, gap 150 m east of the axis, floor ~43 m, flank on the axis
224 m; the nominal strikes at 70 m. `cbf_gap_final` and every variant (terrain hard/soft,
fan 0/6/25/30°, lookahead 0.8–2.5 km, envelope soft/hard, elevator priced 4× the aileron):
terrain strike, usually with a stall. Structural reason: within a fan wide enough to keep
the barrier continuous, some direction always sees a mountain, so the safe set along the
final is EMPTY and the QP is in best-effort for the whole approach; a ray narrow enough to fit
a 300 m gap is discontinuous under yaw. A one-step QP cannot plan a 150 m lateral jog 1 km
ahead from a clearance barrier. **Conclusion: the gap is a corridor problem, not a
terrain-clearance problem** — the Phase 3 cross-track funnel, with the gap as an intermediate
gate on the corridor axis (a piecewise axis), gives the QP a lateral barrier that IS
satisfiable along the final; the terrain row then guards. That is the plan for Phase 3.

Compute with the fan row: ~1.3 ms per step mean (51 bumps × 21 lookahead points in
order-2 jets).

## Rework proposal (2026-09-27): two-layer CBFs — steer at the guidance level, guard at the actuator level

What the evidence says. The actuator-level (6-DOF, surfaces-only) rows work for what they
are: envelope invariance holds against hostile nominals and they clip a real transient
(gust in the turn). They do NOT steer: a position barrier seen through five dynamic stages
(surface → rate → attitude → velocity → position) with a one-step QP either acts too late
(clearance forms), or must freeze its geometry per step and becomes discontinuous (rays), or
must widen its geometry to stay continuous and then has an empty safe set (fans). Steering
needs a barrier that is a smooth function of the COURSE, with the course a state and the
turn command close to it in relative degree.

Proposed architecture (hierarchical CBFs):

1. **Guidance-level CBF-QP** on the kinematic airplane, the layer the energy-reach filter
   already uses (bank-to-turn instead of skid-to-turn, as `energy_ceiling_touchdown_notes.md`
   §2 asked for): states (N, E, h, χ, γ, μ, V), drift Ṅ = V cosγ cosχ, Ė = V cosγ sinχ,
   ḣ = V sinγ, χ̇ = g tanμ / V, V̇ = (V_c − V)/τ_V (speed stays with the nominal for now),
   inputs (μ̇, γ̇) so every geometric barrier has relative degree 2 in both: terrain clearance
   with a SINGLE lookahead ray h − h_T(N + s cosχ, E + s sinχ) − m (smooth in χ, so no frozen
   direction and no fan: the gap is a ray through the floor), the corridor cross-track funnel
   and course alignment (Phase 3), the energy floor/ceiling (Phase 3, `ReachEnergyBarrier`
   already lives here), state boxes on μ and γ (RD 1), rate boxes on the inputs. Exact Lie
   derivatives from `lie_taylor.hpp` on a 7-state model: cheap. Output: (γ_c, μ_c) to the
   nominal's existing `stepOuterRef` / `stepInnerPhi`; the nominal (pattern guidance) supplies
   the nominal rates ((γ_ref − γ_c)/τ, (φ_cmd − μ_c)/τ). Anti-windup through the guidance
   filter as in Phase 1.
2. **Actuator-level CBF-QP** (Phase 1 as built) stays downstream as the ENVELOPE guard: α, φ,
   β, V hard, limits consistent with the guidance model's assumed capability (guidance |μ| ≤ 30°
   ⇒ bank row at 35°), plus the Phase 4 impact/flare row, which is a pitch-axis problem and
   belongs at this level. The terrain climb-capability row can stay as a last-resort guard.
3. **Interconnection argument** (what the theorist will ask): the inner loops track (γ_c, μ_c)
   with a bounded lag; measure the closed-loop step responses of the 6-DOF loops (τ_γ, τ_μ,
   overshoot), fold them into the guidance rows as first-order lag states or as a margin
   (ISSf), and state the guarantee as: kinematic invariance with margin, envelope invariance on
   the full model.

Sequence: (a) guidance model + terrain ray row + boxes, oracle and adversarial tests on the
kinematic model; (b) wiring and the A/B on `fail_gap_final`, `fail_ridge`, `fail_mountain`
(these should flip; if the gap still does not, the gap gate of Phase 3 is next); (c) the
landing cases untouched; (d) Phase 3 rows at this level.

## Nominal-failure cases (2026-09-26): where a barrier row is actually needed

`data/beaver_corridor_fail_*.yaml`, all derived from the overhead entry. Each is built so
the pattern nominal alone does NOT land safely; the failure signature is pinned by
`test_sixdof_sim.cpp` `[corridor][baseline]` so the phase that adds the row has a ground
truth to flip (then the case moves into a `filter: on/off` A/B pair). Nominal-only outcomes:

| Case | What is wrong | Nominal-only outcome | Row that must flip it |
|---|---|---|---|
| `fail_ridge` | 45 m tree line 450 m short of the aim point (glideslope 27.5 m there) | **terrain strike** on final, h_T = 30.5 m | Phase 2 terrain clearance (RD 3, hard) → climb over, steeper bounded final; or Phase 5 go-around |
| `fail_mountain` | 300 m mountain under the LSR pattern loop | **terrain strike** mid-pattern at 164 m | Phase 2 terrain row + terrain-aware candidate selection; Phase 5 go-around if infeasible |
| `fail_short_lake` | 500 m corridor, 46 m/s approach (water run ≈ 460 m, far edge 400 m past the aim point) | touchdown in corridor at 46 m/s, **rollout stops 63 m past the far edge** (into the far-shore trees) | Phase 3 energy CEILING to far edge − stopping run; Phase 4 impact row for the resulting flare |
| `fail_engine_out` | throttle authority idle from t = 0 | follows the plan on an idle glide, **lands 1.4 km short**, α → 26° (past stall, outside LR-556), V → 28 m/s | Phase 3 energy FLOOR to the near edge (reach barrier: best-glide α, direct arc, abandon the pattern) |
| `fail_gust_final` | 40 m corridor, 12 m/s crosswind gust on final (~70 m up) | cross-track held to < 1 m (gusts to 15 m/s were tried — the ground-course tracker never loses the axis), but **floats contact with 17° of crab** | Phase 4 decrab row (ψ − χ, RD 2 in the rudder, gated below h_decrab) with the Phase 1 bank / sideslip rows |
| `fail_slow` | starts at 36 m/s, nominal decelerates to 31 m/s | **enters the stall region** (α 18.4°, FAILED by definition; the plant has no post-stall aero, numbers past 16° are not physical) | Phase 1 AoA ceiling — DONE, see `cbf_slow` |
| `fail_tailwind` | 500 m corridor, 8 m/s tailwind at the normal 40 m/s airspeed | touchdown in corridor on airspeed, **rollout stops 109 m past the far edge** | Phase 3 energy ceiling to far edge − stopping run (ground speed) |
| `fail_gust_turn` | 18 m/s crosswind gust in the first turn | β² pitch-up, **α 17.4°, stall region** | Phase 1 AoA ceiling — **DONE, flips to LANDED OK** (`cbf_gust_turn`) |
| `fail_gap_final` | two mountains on the final, gap 150 m off the axis | **terrain strike** on the flank | Phase 2 terrain row — tried, **does not** (negative result above); Phase 3 corridor gate |
| `fail_tight_pattern` | planner radius 307 m (needs 29° bank) vs the 25° guidance limit | arrives, but the **bank command is saturated 23 % of the flight**, flown bank 26.4°, 47 m off path | Phase 1 bank row — DONE, see `cbf_tight_pattern` (flown bank held at 25.00°) |

Two negative results worth keeping: (a) lateral gusts on final do not defeat the tracker
(course tracking absorbs crab), so a cross-track-funnel demonstration needs a different
forcing (e.g. a late, large lateral offset at low altitude, or the tight-pattern
overshoot); (b) the cone handles "high" entries (up to the γ_min clamp) without any
barrier — the energy ceiling is needed for "hot / short lake", not for "high".
Figures: `figures/corridor_fail_*.png`; all ten runs are tabs in the committed replay page
`figures/corridor_replay.html` (self-contained; `scripts/refresh_corridor_replay.sh`).

## Context

Goal: set a **landing corridor** (oriented rectangle long enough for touchdown plus the stopping run) and **keep-out terrain** (trees ringing the lake, a mountain, houses) and have CBFs on the full 6-DOF Beaver plant guard a pattern-flying nominal so the aircraft enters from anywhere, flies the right approach, touches down gently inside the corridor, and stops inside it. Powered go-around when the corridor becomes infeasible.

What exists today (branch `6dof`, energy-reach layer uncommitted): validated Beaver 6-DOF plant with autodiff-templated `xdotT`; cascade/TECS longitudinal nominal + cross-track-PD lateral; RK4 sim that **terminates at touchdown**; longitudinal HOCBF stack (`hocbf.hpp`, `lon_cbf_filter.cpp`) on a reduced model; a point-mass energy-reach CBF that knows one aim point. **Missing**: north position in the state, any ground/rollout dynamics, any heading/course guidance, any corridor/terrain/obstacle geometry, and any CBF at the actuator level on the 6-DOF plant.

User decisions (fixed): water landing with a simple planar rollout model; **full 6-DOF actuator-level CBFs** over (de, da, dr, dT); keep-out as a **smooth terrain height field**; **target-tracking pattern nominal** that CBFs guard; oriented rectangle corridor, **one landing direction**; start **anywhere** (overhead, opposite heading); **go-around** in the target solution; **this session implements Phase 0 (foundations) only, no CBF rows**.

## Phase overview

| Phase | Deliverable |
|---|---|
| **0 Foundations (now)** | North in the state; corridor + terrain geometry + `world:` YAML; post-touchdown rollout; Dubins-pattern + glideslope-cone + L1 nominal; tests; plot; docs. Demo: pattern nominal alone lands in the corridor on a flat lake with terrain defined, from overhead/opposite heading. |
| 1 6-DOF CBF engine | Engine-lag augmented CBF model, moment-only g, Taylor-jet completeness for `xdotT`, generic multi-row QP filter over (de,da,dr,dT), envelope/actuator rows, A/B twin YAMLs. |
| 2 Terrain keep-out | RD-3 terrain-clearance hard row; terrain-aware Dubins candidate selection. |
| 3 Corridor rows | Cross-track funnel, course alignment, energy floor to near edge / ceiling to far edge minus stopping run. |
| 4 Touchdown rows | TN 1516 impact row on the 6-DOF plant (flare emerges), nose-up floor, decrab, waves on. |
| 5 Go-around + suite | Mode machine, infeasibility triggers, start-anywhere suite, Monte Carlo, model-mismatch study. |

---

## Phase 0 — Foundations (implement now)

### 0.1 Promote north position into the state (NX 11 → 12)

Append `XN = 11` at the end of `enum State` in `include/autoland/types.hpp` (name `XN`, not `N`: `lie_taylor.hpp` uses `template <int N,...>`); `NX = 12`. Existing indices 0..10, `head<3>()`, and `kLonStates/kLatStates` (`linear_model.hpp:55-58`) stay unchanged; `XN` belongs to neither sub-model (comment it).

North row (same formula as `groundKinematics().xdot_n`, `src/sixdof_sim.cpp:64-65`):
`ndot = u·cθ·cψ + v·(sφ·sθ·cψ − cφ·sψ) + w·(cφ·sθ·cψ + sφ·sψ)`

Consumers:
- `include/autoland/beaver_dynamics.hpp` `xdotT` (lines 169-179): add `xd[XN] = ndot`.
- `src/dynamics.cpp` `Dynamics::xdot` (VSPAERO plant, ~177-186): same row.
- `src/linear_model.cpp:14`: add `s[XN]` to the FD step vector (legacy path only).
- `src/beaver_dynamics.cpp`, `src/trim.cpp`: nothing (loop over NX / `Zero()`).
- `src/sixdof_sim.cpp`: delete `x_pos` (569, 706); use `x[XN]` for the wave field (600-601), CSV `x` column (649), `td.x` (678), `dist_target` (685), `sink_rel` (688).
- `apps/beaver_validation.cpp`: `--xdot` keeps the 18-field input, sets `x[XN]=0`, emits 12 columns (149); `--doublet` header (173) adds `n`.
- `scripts/validate_beaver_sixdof.py`: 11→12 rows (`errs`, `range(11)`, ROWS add `Ndot`), add `ndot` to `full_xdot` (~line 58). `validate_beaver_modes.py` slices 8 states, unaffected.
- `test/test_sixdof_sim.cpp`: `sampleState` sets `x[XN]`; wind-invariant kinematic-row list (~98) adds `XN`.
- `test/test_beaver_dynamics.cpp`: level-turn case adds a check on `xd[XN]`; new test "north row matches ground-kinematics oracle" on random states.
- Docs: `types.hpp` header comment (drop "downrange is NOT a state"), `beaver_validation.md` 12-state note.

Bit-identity: rows 0..10 unchanged; the CSV `x` column moves from trapezoid to RK4 (~1e-6 m), which only feeds the wave field and diagnostics. Do this step first and get all existing tests green before anything else.

### 0.2 World geometry — new `include/autoland/world_geometry.hpp` (header-only, templated on scalar)

```cpp
struct Corridor {            // oriented rectangle, ONE landing direction
  double N_c, E_c;           // center
  double length, width;      // along / across [m]
  double heading;            // landing direction psi_c [rad], N->E positive
  double s_aim;              // aim point in corridor frame (near edge = -length/2)
  template<class T> void frame(const T& N,const T& E,T& s,T& e) const;  // s along (+ toward heading), e right+
  template<class T> void inverse(const T& s,const T& e,T& N,T& E) const;
  bool contains(double N,double E,double margin=0) const;
  double sNear() const; double sFar() const;
  Eigen::Vector2d aimPoint() const; Eigen::Vector2d axisPoint(double s) const;
};
struct TerrainBump { double N0,E0,height,sigma_a,sigma_b,rot; int order{1}; }; // h = A·exp(−q^order), q=(xa/σa)²+(xb/σb)²
struct TerrainField {
  double base{0}; std::vector<TerrainBump> bumps;
  template<class T> T height(const T& N,const T& E) const;   // C∞; only * + exp so Taylor jets evaluate it exactly in Phase 2
  Eigen::Vector2d gradient(double N,double E) const;         // autodiff::dual
  void addRing(double Nc,double Ec,double radius,int count,double height,double sigma,int order,
               double skip_heading=NAN,double skip_halfangle=0);  // trees around the shore, optional gap
};
struct RolloutConfig { bool enabled; double a0,kq,k_psi,V_stop,t_max; };
struct WorldConfig { bool enabled{false}; Corridor corridor; TerrainField terrain; RolloutConfig rollout; };
```

Integer `order` (1 Gaussian, 2 flat-top) keeps the field to `*`, `+`, `exp`, which `Taylor` in `lie_taylor.hpp` already supports (no `pow`/`log`).

`world:` YAML block (parsed in the `SixDofSim` constructor with `getOr/getOrB`, `src/sixdof_sim.cpp:19-24`; degrees→radians at the boundary):

```yaml
world:
  enabled: true              # absent/false -> legacy straight-in, new CSV columns log 0
  corridor: {center_N: 0.0, center_E: 0.0, length: 700.0, width: 80.0, heading_deg: 0.0, aim_s: -250.0}
  terrain:
    base: 0.0
    bumps:
      - {name: mountain, N: 3000.0, E: 900.0, height: 250.0, sigma_a: 400.0, sigma_b: 250.0, rot_deg: 30.0, order: 1}
      - {name: house_a,  N: -800.0, E: 250.0, height: 8.0,   sigma_a: 6.0,   sigma_b: 6.0,   order: 2}
    rings:
      - {name: shore_trees, center_N: 0.0, center_E: 0.0, radius: 600.0, count: 48, height: 15.0, sigma: 25.0, order: 2}
  rollout: {enabled: true, a0: 0.8, kq: 0.002, k_psi: 0.5, V_stop: 0.5, t_max: 90.0}
```

Demo geometry is chosen so the terrain-blind Phase-0 nominal clears the trees by construction (aim point ~350 m inside the ring; at γ = −4° the glideslope is ~24 m over 15 m trees). Phases 1+ tighten it so rows bind.

### 0.3 Post-touchdown rollout — new `include/autoland/water_rollout.hpp` (templated RHS, local RK4)

State `[N, E, chi, Vg]`; corridor heading `psi_c`:
```
Ndot = Vg cos chi;  Edot = Vg sin chi
chidot = -k_psi * wrap(chi - psi_c)      (water-rudder/pilot steering proxy)
Vgdot  = -(a0 + kq Vg^2), Vg >= 0
```
Init from touchdown: `N,E = x[XN], x[Y]`; `Vg = hypot(xdot_n, ydot_e)`; `chi = atan2(ydot_e, xdot_n)`; record crab `psi_td − chi_td`. Stop when `Vg <= V_stop` or `t − t_td > t_max`. Closed forms for tests: `kq=0` → `d = V0²/(2a0)`; `a0=0` → `d = ln(1+kq V0 t)/kq`. `a0`, `kq` are labelled placeholders (POH-calibrate later; add to TODO Modeling/data).

Result nested in `SixDofTouchdown` (signature of `run()` unchanged):
```cpp
struct SixDofRollout { bool ran; double t_start,t_stop,N_stop,E_stop,chi_stop,distance,s_stop,e_stop,max_abs_e; bool inside; };
```
Logging: same CSV, new `phase` column (0 airborne, 1 rollout); rollout rows carry `t,x,y,psi_deg,Vg,s_corr,e_corr`, everything else 0. All new columns are **appended** to the header (`src/sixdof_sim.cpp:578-585`) so `plot_sixdof_results.py` / `plot_energy_reach.py` (both read by name) keep working.

### 0.4 Pattern / target-tracking nominal — new `include/autoland/pattern_guidance.hpp` + `src/pattern_guidance.cpp`

Plan once at t=0, track a fixed path (replanning every step chatters on-path; leave "replan if cross-track > threshold" as a TODO).

1. **Dubins CSC planner** `dubinsCSC(start_pose, goal_pose, R)`: best of LSL/RSR/LSR/RSL (LSL/RSR always exist → a path exists for any start incl. overhead/opposite heading; CCC omitted, documented). Shkel & Lumelsky (2001) formulas, self-implemented. `R = k_R·V_app²/(g·tan φ_plan)`, defaults `φ_plan = 20°`, `k_R = 1.15` (~500 m at 40 m/s).
2. **Path** = `{Arc | Line}` segments; final segment is the corridor axis line from the FAF through the aim point and past the far edge (tracker never runs out). API: `closest(N,E) → {seg, s_path}`, `pointAt(s)`, `L_to_aim(s)`.
3. **FAF placement**: on-axis `d_final` before the aim point, `h_faf = d_final·tan(−γ_app)`. Fixed-point (5 iters) on `d_final = clamp(h0/tan(−γ_app) − L_dubins(d_final), d_final_min, d_final_max)` so the whole pattern flies near `γ_app`; excess altitude handled by the cone below.
4. **Longitudinal reference = glideslope cone on distance-to-go**: `h_gs = L_to_aim·tan(−γ_app)`; `gamma_ref = clamp(γ_app + K_h·(h_gs − h)/V, γ_min, γ_max)` (defaults −8°/+6°). Add `SixDofNominal::stepOuterRef(x, V_air, Vdot_air, gamma_ref, V_ref, dt)` in `sixdof_nominal.hpp`; existing `stepOuter` (line 136) becomes a one-line call with `c_.gamma_ref, c_.V_ref` (bit-identical). TECS works unchanged (`hdot_sp = V_ref·sin γ_ref`, line 159).
5. **Lateral = L1 course tracking → phi_cmd** (Park, Deyst, How, AIAA GNC 2004): lookahead `L1 = max(L1_min, L1_period·0.75·Vg/π)` along the path from the closest point; `η` = angle from ground-velocity course `chi = atan2(ydot_e, xdot_n)` (crab in wind automatic) to the lookahead; `a_s = 2Vg²·sin η/L1`; `phi_cmd = clamp(atan(a_s/g), ±phi_max)`. Add `SixDofNominal::stepInnerPhi(x, theta_cmd, dT, phi_cmd, dt)` mirroring `stepInnerBeta` (line 239); `stepInner` stays bit-identical. Roll PD / yaw damper reused (Kp_p Nyquist-limited, untouched); `phi_max_deg: 25` only in corridor YAMLs.
6. No flare (as today). Flare comes from the Phase 4 impact row.

Per-step outputs for CSV: `chi_cmd_deg, N_ref, E_ref, h_gs, L_to_go, seg, s_corr, e_corr, hT, phase`.

Sim wiring (`SixDofSim::run` control block, `src/sixdof_sim.cpp:609-640`): third branch `else if (world.enabled)`: `guid = guidance.step(x, gk)`; `oc = nominal.stepOuterRef(..., guid.gamma_ref, V_app, dt)`; `u = nominal.stepInnerPhi(x, oc.theta_cmd, oc.dT, guid.phi_cmd, dt)`. `reach.enabled && world.enabled` → throw.

### 0.5 Initial condition, touchdown semantics, stats, console

- `initial:` gains absolute `N0, E0, psi_deg` when `world.enabled` (plus existing `h0, dV`): `x0_[XN]=N0; x0_[Y]=E0; x0_[PSI]=psi0`.
- Touchdown test (line 666): `h <= max(eta, hT(N,E))`. `hT > 0.05` → `td.terrain_strike = true` (no rollout). Water touchdown → `td.in_corridor`, `td.s_corr, td.e_corr, td.dpsi_corr = wrap(psi − psi_c)`, `td.crab`, then rollout.
- `SixDofRunStats` add `min_terrain_clearance, max_abs_e_final, plan_length, d_final`.
- Console: corridor summary block (touchdown s/e/heading error, rollout stop s/e, inside yes/no).

### 0.6 Scenarios (`data/`)

- `beaver_corridor_straight.yaml` — 3 km out on the extended centerline at glideslope altitude (regression baseline; should reproduce today's calm touchdown numbers).
- `beaver_corridor_overhead.yaml` — the demo: over the lake center at 250 m, heading opposite the corridor, terrain block above, rollout on.
- `beaver_corridor_crosswind.yaml` — overhead + 5 m/s step crosswind (`wind:` block).
These become the `filter: false` baselines Phase 1 pairs with.

### 0.7 Tests (Catch2; add files to `unit_tests` in `CMakeLists.txt`, sources to `add_library(autoland_core)`)

- `test/test_world_geometry.cpp`: rotated-corridor frame hand values + round-trip; `contains` at corners ±ε; terrain height at bump center/far field; **gradient oracle** (hand chain rule for a rotated Gaussian and an order-2 bump vs `gradient()`, 1e-12); **Lie oracle** `lieDrift<1>` of `height()` along a constant-velocity flow = `grad·v`; ring count/radius.
- `test/test_water_rollout.cpp`: both closed forms (1e-6), heading converges to `psi_c`, 33 m/s touchdown at `s_aim` stops inside a 700 m corridor, 45 m/s at the far third stops outside.
- `test/test_pattern_guidance.cpp`: Dubins endpoint reached by independent forward integration for 20 random pose pairs (1e-9); on-axis start → straight line; abeam-2R parallel case → `πR + d`; `d_final` fixed point monotone in `h0`; L1 sign on lateral offset, zero on-path; cone clamps.
- `test/test_sixdof_sim.cpp` additions: NX edits; end-to-end `corridor_straight` (inside, `|e|<3 m`, `|dpsi|<2°`, rollout inside) and `corridor_overhead` (`|e|<8 m`, `|dpsi|<4°`, rollout inside, `min_terrain_clearance > 5 m`, tag `[slow]`). Existing calm/hot/xwind/wave tests must pass unchanged.

### 0.8 Plot + docs

- `scripts/plot_corridor_landing.py CSV SCENARIO.yaml OUT.png`: ground track with terrain contours (bumps re-evaluated from YAML), corridor rectangle, planned path, touchdown + stop markers, rollout track; altitude vs distance-to-go with cone and `hT` under track; corridor-frame (s,e); φ/φ_cmd, χ/χ_cmd; V, γ; rollout Vg. Figures → `figures/corridor_*.png`.
- Docs: `documentation/CHANGELOG.md` entry (newest on top, template); `TODO.md` (check the TECS "downrange origin" item, add Phase 1+ items and the rollout-coefficient placeholders); new `documentation/corridor_landing_design.md` = Phase 0 model + this roadmap; cross-link from `energy_ceiling_touchdown_notes.md` §3 and `water_landing_cbf_design.md` doc map; `beaver_validation.md` 12-state note.

### 0.9 Session order / minimum cut

1. NX promotion, build, all existing tests green. 2. `world_geometry.hpp` + tests. 3. rollout + tests. 4. Dubins/path/L1/cone + tests. 5. sim wiring + YAMLs + end-to-end tests. 6. plot + docs.
Minimum viable cut: CSC-only Dubins, plan once, cascade nominal only, straight-in test with tight tolerances; overhead scenario as a demo with looser tolerances.

---

## Phases 1–5 (roadmap; not implemented this session)

### Phase 1 — 6-DOF CBF engine
- **CBF model** `include/autoland/sixdof_cbf_model.hpp`: `X = [x(12); dT_s]` (NXC=13), `U = [de, da, dr, dT_cmd]`, engine lag `dT_s_dot = (dT_cmd − dT_s)/τ_eng` (τ ≈ 0.5–1 s). Puts every non-affine throttle term (`dpt²`, `dpt³`, `dpt ∝ P/V³`) into the drift where the Taylor jet is exact, and raises throttle RD by one so it aligns with the elevator on every row. 6-DOF analogue of the lon model's single-integrator power (`water_landing_cbf_math.md` §5).
- **Direct control-force terms** (`Cz_de`, `Cy_da`, `Cy_dr`, `Cx_dr`): drop from `g`, keep in the drift frozen at the last applied `u` (same RD outcome as `lon_augmented.hpp`). `g(X)` is moment-only: `de→qdot` via `Cm_de`; `da/dr→pdot,rdot` via the Γ solve; `dT_cmd→dT_s_dot`. Remaining control terms (`Cz_deb2`, `Cl_daa`, `Cy_dra`) are affine.
- **Taylor plumbing**: refactor `xdotT` into a container-agnostic `xdotCoreT<T,VecX,VecU>` (Eigen wrapper kept, bit-identical); add `operator/(double,Taylor)`, `atan`, `asin`, `tan` to `lie_taylor.hpp`; identity tests: order-0 jet == `xdot`, order-1 jet == `A·f` from `linearize`.
- **Filter** `sixdof_cbf_filter.{hpp,cpp}`: generalize `src/lon_cbf_filter.cpp:36-220` (hard/soft stacking, quadratic slacks, best-effort with 1e6 hard penalty) to NU=4, with the reach filter's row normalization by `max|A|` and closed-form feasibility check (`src/energy_reach_cbf.cpp:66-80`); weights `1/u_max²`; box = deflection **and** rate limits vs `u_prev`. Wiring: `u = filter(u_nom, X)` between `stepInnerPhi` and RK4, gated by `sixdof_cbf: {enabled, filter}`; `filter: false` = monitor-only twin.
- **Rows**: AoA ceiling (RD 2 all, soft 1e5); bank ±φ_max (RD 2 da/dr); sideslip ±β_max (RD 2); airspeed floor/ceiling (RD 2 all); throttle-state box (RD 1 dT, hard). Oracle tests: hand-derived `L_f h, L_f² h` + control row for the bank and airspeed barriers.

### Phase 2 — Terrain keep-out
`h_T = h − TerrainField::height(N,E) − m_T` (optional `+k_v·V`). Uniform **RD 3** in all four controls (kinematic → forces with `dT_s` → `pdot,qdot,rdot,dT_s_dot`). Hard row, class-K sized to roll/path time constants, gated by clearance (< ~150 m) for cost. Planner becomes terrain-aware: evaluate the CSC candidates (and `d_final` extensions) for clearance along the path, choose the shortest clear one, add a downwind extension if none.

### Phase 3 — Corridor rows
| Row | h | RD | policy |
|---|---|---|---|
| Cross-track funnel | `w(s) ∓ e`, `w(s) = w_td + k_w·smoothmax(s_aim − s, 0)` (sqrt smooth-max as in `energy_reach_cbf.hpp:166-171`) | 3 all | gated past FAF; soft 1e3 |
| Course alignment | `dχ_max(s) ∓ wrap(χ − ψ_c)`, `χ = atan2(Ė, Ṅ)` (Taylor `atan2` exists) | 2 | gated on final |
| Energy floor to **near edge** | `E − (½V_min² + g·d_near·Ψ(α))`, `d_near` via `ReachEnergyBarrier::arcDistance` to the near-edge axis point | 2 all | small ck |
| Energy ceiling to **far edge − stopping run** | `(½V_td_max² + g·(d_far − d_stop(V))·Ψ_cap) − E`, `d_stop = V²/(2a0)` | 2 all | capability-priced Ψ; thrust in drift via `dT_s` |

### Phase 4 — Touchdown rows
TN 1516 impact barrier on the 6-DOF state (`τ = θ − τ_keel`, `γ0` from `hdot/V`, sink `−hdot`), RD 2 (de, dT via force channel), hard, gated `h < z_gate` & descending & `τ > 0` as in `lon_cbf_filter.cpp:104-116`. Nose-up floor (RD 2 de, soft); decrab `ψ − χ` (RD 2 dr, gated `h < h_decrab`); waves on. Demo: gentle touchdown inside the corridor, rollout inside.

### Phase 5 — Go-around + suite
Mode machine `{PATTERN, FINAL, GO_AROUND, ROLLOUT}` (CSV `mode` column). Triggers above `h_decision` (~15 m; below it, commit): corridor soft-row slack above threshold for `T_hold`; empty floor/ceiling band at current distance; predicted touchdown `s + h/tan(−γ) > s_far − d_stop(V)`; hard-row best-effort steps > N. Action: `γ_ref = +γ_ga`, `V_ref = V_climb`, straight along the axis (corridor rows off, terrain/envelope rows on) to `h_pattern`, replan Dubins to the FAF, back to PATTERN. Suite: 8 start poses × wind × waves; Monte Carlo; mismatch study (frozen direct-force terms, engine-lag) for the paper.

### Risks and mitigations
- **RD-3 HOCBF at dt=0.01**: degree 3 already exercised by the lon energy rows; keep g moment-only, do **not** add surface-lag states (would make position rows RD 4). Log predicted vs realised `ḧ` per row.
- **Direct-lift mismatch at the flare** (`Cz_de` ~1.2 m/s² per 0.3 rad at 40 m/s): rate limit in the QP box bounds step-to-step change; quantify in Phase 5.
- **Corridor squeeze infeasibility**: small class-K gains, capability pricing, corridor rows soft, go-around trigger.
- **Terrain local minima**: keep-out is a guard, not a planner; terrain-aware Dubins in Phase 2; report `best_effort_steps` honestly.
- **Envelope honesty**: Beaver aero valid 30–55 m/s, no post-stall; AoA/airspeed rows keep the QP in-band; state in docs.
- **Roll authority**: Kp_p Nyquist-limited; φ_cmd clamp 25°; planner radius margin `k_R`.

---

## Verification (Phase 0)

1. `cmake -S . -B build && cmake --build build -j && ctest --test-dir build --output-on-failure` — all existing tests green after the NX promotion before any new feature lands; new tests green at the end.
2. `python scripts/validate_beaver_sixdof.py` exits 0 with 12 rows (independent Python `ndot`).
3. `./build/sixdof_autoland_sim data/example.stab data/aircraft.yaml data/beaver_corridor_straight.yaml runs/corridor_straight.csv` reproduces today's calm touchdown sink/V/y within tolerance; console reports in-corridor + rollout inside.
4. Same for `beaver_corridor_overhead.yaml` and `_crosswind.yaml`; `python scripts/plot_corridor_landing.py` renders `figures/corridor_overhead.png` showing the Dubins pattern, terrain contours, corridor, touchdown and stop markers.
5. Legacy scenarios (`beaver_landing_calm.yaml`, `beaver_reach_glide.yaml`) produce identical results for columns other than `x` (which changes at ~1e-6 m).
