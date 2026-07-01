# TODO / Backlog

Granular, forward-looking capture of things we want to do but haven't yet — the low-friction
place to jot tasks, cleanups, and ideas so a collaborator can see them. Keep entries short.

**How this differs from the other docs:** `CHANGELOG.md` records what's *been done*
(chronological); `documentation/water_landing_cbf_design.md` §9 holds the *rationale* for the
big research threads. This file is the *actionable checklist*. When an item is done, check it
off and summarize it in `CHANGELOG.md` (then it can be removed here).

Optional per item: `(owner)` and a one-line note.

---

## CBF / theory

- [ ] **Thrust term in `a_brk`.** Add the thrust contribution to the braking acceleration. Needs
      the mixed-relative-degree treatment — putting the thrust state `T`/pitch `θ` into the
      descent barrier drops it below degree 3 — so it couples with the §3.3 work below.
- [x] **§3.3 contact-force barrier (partial).** Attitude-coupled contact-load constraint
      now exists: the hydrodynamic **impact-load barrier** (NACA TN 1516) bounds the peak CG
      load factor at touchdown, built as a **degree-2 HOCBF enforced by the elevator/flare**
      (`impact_barrier.hpp`). A single affine QP row can't mix the elevator's degree 2 with
      thrust's degree 3 (it yields non-affine u²/u̇ terms), so thrust is left to bound impact
      via the sink-rate barrier. *Still open:* (a) a **full mixed-degree (2/3)** construction
      so thrust jointly enforces the impact barrier; (b) a **predictive / backup-set CBF**
      (forward-integrate the planned flare; spec §5.1) that would use both controls and drop
      the height term. Highest remaining research value.
- [ ] **Independent control-row oracle.** The finite-difference oracle covers the drift stack
      `L_f^k b`; add an independent check of the control row `L_g L_f^{r-1} b`.
- [x] **Upper airspeed barrier** `b = V_max − V` (over-speed / high-energy impact). Wired into
      the filter (`airspeed_upper` flag + `Vmax_air`), degree 3, soft by default; reuses the
      lower barrier's machinery with flipped authority signs. `V_max` is a placeholder (see
      Modeling/data). Like the lower barrier it is currently **soft** (see "Harden the airspeed
      barrier" below — applies to both).
- [ ] **Harden the airspeed barrier** — currently soft (slack-penalized), not a hard guarantee.
- [ ] **Tune the impact-load barrier.** `Nb` (Φ budget) from the worst nominal-descent
      excess load, `zs` (Φ altitude scale) so `Φ'(0)=Nb/zs` is matched by the flare
      authority, the `c_impact` class-K gains against the flare timescale, and the Option C
      slack schedule (`impact_slack_lo/hi`). Currently placeholders (`Nb=10`, `zs=2`,
      `c_impact=[2,2]`) — soft + non-binding on the nominal approach.
- [ ] **Impact-barrier altitude estimation (deferred, spec §5.4).** `Φ(z)` currently uses
      the true keel height. With a real estimator, overestimating `z` is the dangerous
      direction — consume a lower-bound `z`, desensitize `Φ'` with a larger `zs` + light
      smoothing, and convert any rangefinder beam range to vertical height via attitude.
- [ ] **Controllability guard** on the elevator authority `L_g L_f^{r-1} b` (→ 0 at very low
      dynamic pressure). Affects the descent barrier and **especially the impact-load barrier**,
      which is elevator-only and must flare at the lowest-q point of the flight — its soft slack
      would silently *absorb* an authority shortfall rather than enforce the load limit. Thrust
      (the actuator that could add flare energy) is excluded by the degree-2 construction.

## Modeling / data

- [x] **True longitudinal trim (for the IC).** `lonTrim()` Newton-solves the lon EOM for a steady
      **level-flight** initial condition (`γ=0`, `V̇=γ̇=q̇=0`) — `src/lon_sim.cpp`. Note: this
      revealed the residual descent-`psi < 0` is **not** the IC (it persists, unchanged) — it is a
      terminal effect at touchdown (`h→0`); the actual barrier `b` stays ≥ 0 in flight. *Optional
      follow-up:* also trim the nominal feedforward on the lon model (still uses the body-axis
      approach trim — harmless, since the controller has feedback).
- [ ] **Compute `C_L,max`** for the airframe. `a_brk(V,gamma)` uses a placeholder `CL_max = 1.2`
      (config). Needs a real value (CFD/wind-tunnel/flight); also the max-lift drag it implies.
      *Note:* the new viscous-stall plant model gives a vehicle CLmax ≈ **1.44** (NeuralFoil 4414,
      `documentation/stall_model_spec.md`) — a better default than 1.2; the descent barrier's
      `CL_max` knob is **not** auto-linked to it, so update it if you want them consistent.
- [ ] **Calibrate / extend the viscous-stall plant model** (`stall-model`). Onset/level are
      tunable knobs, not experiment (NeuralFoil→Xfoil→reality is weakest at low-Re post-stall).
      Calibrate `A_STALL_DEG`/`BLEND_HALF_DEG`/`SEVERITY` in `scripts/precompute_stall_table.py`
      against experimental low-Re 4414 data or flight ID. Follow-ups: a 2nd table axis over Re
      (deck spans 1.5–4×10⁵ in `naca4414_polar.csv`); extend the flat-plate tail past 90°
      (currently clamped) for full tumble; an optional **Goman–Krabrov** separation state for
      dynamic-stall lag/hysteresis (doubles as a CBF state); and an α-margin **stall-recovery
      barrier** (the motivating downstream work).
- [ ] **Calibrate `V_max`** (over-speed barrier ceiling). Currently a placeholder `27` m/s
      (~1.5·V_app) in `lon_scenario.yaml`. Set to the real never-exceed / structural / hull-slam
      speed limit.
- [x] **Retire the deprecated aero decks.** Removed `AHAB_sweep.stab` and `AHAB 2.stab`; the
      body-axis `autoland_sim` default and `test_cbf` now use `AHAB_combined.stab`. **Decision:
      `example.stab` is kept** as the small synthetic fixture for the deck-agnostic unit tests
      (`test_aero_table`/`test_linear_model`/`test_trim`/`test_lon_cbf`).
- [ ] **(Optional) Migrate the unit tests off `example.stab` onto `AHAB_combined.stab`.** Lower
      priority since the tests are deck-agnostic. Non-trivial: the combined deck's β grid is not
      symmetric (β ∈ [0,20], not [−10,10]), so `test_aero_table`'s node-value / clamp / midpoint
      assertions need new query points and re-derived ground-truth values; and `test_linear_model`'s
      static sign-sense checks (e.g. `Cm_α < 0`) may flip on the real near-neutral airframe — a real
      finding to surface, not silently patch. (Was also Brian's "remove `example.stab`" item; folded
      in here under the keep-decision.)
- [ ] **Calibrate the impact-load barrier hull/limit params.** `n_limit` (3 g placeholder)
      from the structural allowable *with knockdowns* (it's a CG load factor normal to the
      water, not local hull pressure — keep that check separate). Confirm hull dead-rise
      `beta` (22.5° placeholder), `tau_keel` (keel incidence, 0), and `rho_water`
      (1000 fresh / 1025 sea). Absolute load level is dominated by dead rise and contact
      sink rate, so placeholders can be off by large factors.
- [ ] **Calibrate `parasite_CD0` and thrust `k_v`** against flight / tow-tank / thrust-stand data.
- [ ] **Confirm the mixing map** — OpenVSP control-group definitions vs the default
      (Elevator←δe, Ailerons←δa, Rudder←δr), or set `mixing.matrix`.
- [ ] **c.g. offset** vs the `.stab` moment reference (default 0) — confirm if it differs.

## Validation / realism gaps

Soft spots in the current impact-load barrier worth being honest about — most are
"optimistic in the dangerous direction" and none have been exercised in the regime the
barrier exists for.

- [ ] **Exercise the impact barrier off-nominal — it is currently unverified.** On the nominal
      approach `n_peak` peaks ~0.17 g vs the 3 g limit, and the sink-rate barrier flares to a
      ~0.014 m/s touchdown, so the impact row never takes command (slack ≡ 0) — even when forced
      *hard* with `n_limit=0.05 g` the trajectory was unchanged. Run a sweep where the kinematic
      flare *can't* save it (steeper γ, higher `V_app`, gust/shear injection, degraded/disabled
      descent barrier) and confirm the impact row actually becomes the active constraint and
      holds `n_peak ≤ n_limit`. Until then "it works" is indistinguishable from "it's off."
- [ ] **Validity gate leaves the worst case unprotected.** The row is assembled only while
      descending with positive trim (NACA TN 1516's valid domain), so a nose-down / high-sink
      botched approach — plausibly the highest-load case — gets *no* protection. The hard on/off
      of the gate is also a discontinuity (1 QP feasibility-recovery seen in the hard binding
      test). Consider a smooth fallback / out-of-envelope load guard.
- [ ] **Idealized estimation & actuators.** Sim uses perfect θ/γ/V/z and no actuator lag; the
      barrier is the most sensitive part (Lie-derivative gains are large — elevator coeff ~3e4 —
      so it amplifies state noise, and the flare assumes an instantaneous elevator). Add sensor
      noise + a servo-lag model and re-check. (Altitude specifically: see the Φ(z) estimation
      item under CBF/theory.)
- [ ] **Smooth-water assumption (no wave model).** `τ`/`γ₀` are referenced to flat water; in a
      seaway the trim/path *relative to the local wave slope* dominate the slam load (a contact
      on a wave face spikes κ). Reference the contact state to the wave surface. Probably the
      biggest single physical unrealism for real ops.
- [ ] **Lift = weight during contact (inherited TN 1516 assumption).** At the decelerating,
      low-speed touchdown the wing may carry < weight, dumping more onto the hull → true load
      *higher* than predicted (non-conservative). Check/bound this.
- [ ] **CG load factor ≠ local hull pressure.** `n_peak` is the rigid-body CG load normal to the
      water, not the local panel pressure that usually drives structural failure. A separate
      local-pressure check should exist (spec §7); `n_limit` alone can't stand in for it.
      *(Calibration of `n_limit` is tracked under Modeling/data.)*
- [ ] **Longitudinal / symmetric only.** No lateral DOF — a wing-down or crabbed touchdown slams
      one chine first and produces the worst *local* loads, entirely outside this barrier.

## Infra / build

- [ ] (none currently)

## Hardware / future

- [ ] **PX4 path.** Insertion point for the QP filter (actuator-level vs rate/attitude setpoint
      sets the relative degree) — a design point worth reporting.

## Ideas / maybe

- [ ] (jot loose ideas here)
