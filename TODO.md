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
- [ ] **§3.3 contact-force barrier.** Attitude-coupled `v_safe(θ)` from von Kármán/Wagner
      slamming theory (mixed degree: 2 for elevator, 3 for thrust). Highest research value.
- [ ] **Independent control-row oracle.** The finite-difference oracle covers the drift stack
      `L_f^k b`; add an independent check of the control row `L_g L_f^{r-1} b`.
- [x] **Upper airspeed barrier** `b = V_max − V` (over-speed / high-energy impact). Wired into
      the filter (`airspeed_upper` flag + `Vmax_air`), degree 3, soft by default; reuses the
      lower barrier's machinery with flipped authority signs. `V_max` is a placeholder (see
      Modeling/data). Like the lower barrier it is currently **soft** (see "Harden the airspeed
      barrier" below — applies to both).
- [ ] **Harden the airspeed barrier** — currently soft (slack-penalized), not a hard guarantee.
- [ ] **Controllability guard** on `L_g L_f²b` (elevator authority → 0 at very low dynamic
      pressure would make the descent barrier unenforceable by elevator).

## Modeling / data

- [x] **True longitudinal trim (for the IC).** `lonTrim()` Newton-solves the lon EOM for a steady
      **level-flight** initial condition (`γ=0`, `V̇=γ̇=q̇=0`) — `src/lon_sim.cpp`. Note: this
      revealed the residual descent-`psi < 0` is **not** the IC (it persists, unchanged) — it is a
      terminal effect at touchdown (`h→0`); the actual barrier `b` stays ≥ 0 in flight. *Optional
      follow-up:* also trim the nominal feedforward on the lon model (still uses the body-axis
      approach trim — harmless, since the controller has feedback).
- [ ] **Compute `C_L,max`** for the airframe. `a_brk(V,gamma)` uses a placeholder `CL_max = 1.2`
      (config). Needs a real value (CFD/wind-tunnel/flight); also the max-lift drag it implies.
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
      finding to surface, not silently patch.
- [ ] **Calibrate `parasite_CD0` and thrust `k_v`** against flight / tow-tank / thrust-stand data.
- [ ] **Confirm the mixing map** — OpenVSP control-group definitions vs the default
      (Elevator←δe, Ailerons←δa, Rudder←δr), or set `mixing.matrix`.
- [ ] **c.g. offset** vs the `.stab` moment reference (default 0) — confirm if it differs.

## Infra / build

- [ ] (none currently)

## Hardware / future

- [ ] **PX4 path.** Insertion point for the QP filter (actuator-level vs rate/attitude setpoint
      sets the relative degree) — a design point worth reporting.

## Ideas / maybe

- [ ] (jot loose ideas here)
