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
- [ ] **Calibrate the impact-load barrier hull/limit params.** `n_limit` (3 g placeholder)
      from the structural allowable *with knockdowns* (it's a CG load factor normal to the
      water, not local hull pressure — keep that check separate). Confirm hull dead-rise
      `beta` (22.5° placeholder), `tau_keel` (keel incidence, 0), and `rho_water`
      (1000 fresh / 1025 sea). Absolute load level is dominated by dead rise and contact
      sink rate, so placeholders can be off by large factors.
- [ ] **Remove the old `example.stab` placeholder aero deck.** Migrate the unit tests
      (`test/*` use `example.stab` via `Setup`) and the body-axis `autoland_sim` default onto
      `AHAB_combined.stab` first, then delete `example.stab`.
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
