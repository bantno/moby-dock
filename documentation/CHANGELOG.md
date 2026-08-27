# Changelog

Chronological log of meaningful changes to this project. **Read the top entry at the start
of a session** to catch up on what your collaborator did since you last worked.

This doc's job — and how it differs from the others:

| Doc | Job |
|---|---|
| **`CHANGELOG.md`** (this file) | *When / who / what / why* — chronological, append-only. |
| `../TODO.md` | *Backlog* — granular open tasks/ideas to pick up. |
| `water_landing_cbf_design.md` | *Current state* — the living theory/approach/assumptions/decisions. |
| `water_landing_cbf_math.md` | *The math* — formal derivations. |
| `archive/` | Superseded historical docs. |

## Conventions

- **Newest entry on top.** Append above the existing entries; never rewrite past ones.
- One entry per working session (or per logically-grouped change). Keep it skimmable.
- Use the template below. Date is `YYYY-MM-DD`; author is a first name so collaborators can
  tell each other apart.

```
## YYYY-MM-DD — <author> — <short title>
**Branch/commit:** <branch> (<short-sha> if committed)
**What changed:** <the concrete change>
**Why:** <motivation / what problem it solves>
**Follow-ups / notes for collaborator:** <anything they should know or pick up>
**Files touched:** <key paths>
```

### Collaborator setup (one-time, per developer)

Each developer keeps their own local `CLAUDE.md` (gitignored, not shared via git). Add this
line to it so your Claude session reads this changelog at startup:

> At the start of each session, read `documentation/CHANGELOG.md` to see what changed since
> the last session before doing any work.

---

## 2026-08-26 — Brian — Beaver 6-DOF plant wired + validated vs FDC references; Beaver landing-case suite
**Branch/commit:** 6dof
**What changed:** The DHC-2 Beaver is now a real 6-DOF plant (`beaver_dynamics.hpp/.cpp`) —
direct evaluation of the verified LR-556/FDC polynomials (no table export), templated scalar so
trim Jacobians and the state-space linearization are EXACT autodiff (no finite differences).
Propulsion is throttle→manifold-pressure→power→dpt at fixed RPM; atmosphere is the exact ICAO
troposphere + inverse-square gravity (both confirmed as FDC's own models). `plant: beaver` is the
`sixdof_autoland_sim` DEFAULT; the AHAB deck stays selectable (`plant: vspaero`, existing
scenarios/tests tagged and green). **Validation (documentation/beaver_validation.md):** (1) the
FDC 1.2 manual's printed ACTRIM check case is reproduced to print precision (all six acceleration
rows ≤3.3e-6); (2) our new 6-axis Newton trim recovers the FDC trim to 0.002° in α/θ and 0.0004°
in δe (the lateral deltas are FDC's own fmins tolerance — its printed β̇ residual is 4.6e-4);
(3) the level-flight trim sweep tracks fig. 10.13 through every digitized point; (4) an
INDEPENDENT Python re-implementation (FDC state coords, freshly-typed coefficients, complex-step)
agrees with the C++ linearization to 4e-10 relative on all eigenvalues — classic Beaver modes
(SP 3.2 rad/s ζ0.67, phugoid 0.26 ζ0.06, roll −5.1, DR 1.1 ζ0.44, spiral −0.05). **Landing
cases** (`data/beaver_landing_{calm,crosswind,gust,waves_lake,poh}.yaml`, plots in results/):
all land on target; gains re-sized for 2288 kg; the nominal gained per-axis control-SENSE signs
(`de_sign/da_sign/dr_sign`) — the Beaver's standard Delft signs (Cm_de, Cl_da, Cn_dr < 0) are
OPPOSITE the AHAB virtual controls on all three axes (with AHAB signs the roll loop is positive
feedback → inverted in ~6 s) — plus δa/δr trim feedforwards (slipstream asymmetry trims at
β≈−0.9°, δr≈−2.3°). Standard cases fly 40 m/s clean (35 clean trims at α≈12.9°, ~3° from stall);
the POH float configuration (flaps 35°, 33.5 m/s, θ_trim 0.2°) is its own case. 76 tests green
(+5: check-case oracle, trim recovery, cruise modes, wind coupling, closed-loop calm landing).
**Why:** paper_readiness §6 — validate the flight-validated plant is implemented correctly before
moving the CBF machinery onto it.
**Addendum (same session):** fig. 10.13's solid curve was PIXEL-DIGITIZED from a 400-dpi render
(`data/fdc_fig1013_solid_digitized.csv`) and identified as a fixed-pz≈20 "Hg trim (γ free): our
fixed-pz curve reproduces it with a uniform ~0.3° offset while the exact-condition ACTRIM point
matches to 0.0004°. New verification figures: `figures/beaver_validation_checkcase.png` (oracle
match + wrong-atmosphere discrimination + trim recovery) and `_fig1013_sidebyside.png` (scan vs
model). Sim CSV precision raised to 10 sig. digits — the default 6 quantized a 40 m/s airspeed
to 0.1 mm/s steps, rendering the tightly-held V channel as a staircase on autoscaled plots.
**Follow-ups / notes for collaborator:** two reference-condition traps documented in
beaver_validation.md: the ACTRIM check case is at the SEA-LEVEL atmosphere (the 2000 ft prompt
seeds only the sim IC), and FDC eq. 3.15's altitude correction multiplies only the (408−0.0965n)
term. Remaining: CBF re-tune on the Beaver, flare/decrab, u_P spool-rate bound.
**Files touched:** `include/autoland/{beaver_dynamics,beaver_aero,sixdof_sim,sixdof_nominal}.hpp`,
`src/{beaver_dynamics,sixdof_sim}.cpp`, `apps/{beaver_validation,sixdof_autoland_sim}.cpp`,
`scripts/{validate_beaver_modes,plot_beaver_validation}.py`, `test/test_beaver_dynamics.cpp`,
`data/beaver_landing_*.yaml`, `data/sixdof_*.yaml` (plant tag), `documentation/beaver_validation.md`,
`CMakeLists.txt`, `results/beaver_*`.

## 2026-08-20 — Brian — Beaver plant (LR-556, verified) + W/2 float impact split + single-integrator-power degree-2 impact barrier
**Branch/commit:** 6dof
**What changed:** Groundwork to swap the self-authored VSPAERO plant for the **flight-validated
DHC-2 Beaver** (Tjee & Mulder LR-556, TU-Delft 1988), to isolate the CBF machinery from the aero
model for the paper (fixes the "self-consistent plant" reviewer gap). Three pieces landed, all
parallel to the AHAB path (69→71 tests, still green):
1. **`beaver_aero.hpp`** — the nonlinear Beaver aero + engine model. The ~60 coefficients were
   VERIFIED against LR-556 Table 3 / Table 2 to every printed digit (cross-checked the scanned
   `documentation/LR-556.pdf` against the GitHub `ftmeeet/Flight_Simulator` `dhc2_vars.mat`; they
   agree, one standard omission — the CY β̇ term). Body-axis, standard signs (no frame flip), Beaver
   rate convention (`qhat=q·c/V`), propulsion via the slipstream coefficient `dpt` (no separate
   thrust force). `test_beaver_aero.cpp` reproduces the LR-556 reference condition end-to-end.
2. **W/2-per-float impact split** — `impact_barrier.hpp` gains `n_surfaces` (default 1 = single
   hull, unchanged). For a twin-float floatplane the single-surface TN 1516 theory is applied per
   float at W/2, raising `n_peak` by 2^(1/3) (conservative for symmetric contact). Threaded through
   the lon filter, lon sim, and 6-DOF sim; documented in `impact_load_barrier_spec.md` §7-8.
3. **`beaver_lon.hpp`** — the Beaver longitudinal augmented model with **single-integrator power**
   (state P, control u_P=Ṗ), which puts engine power on the impact barrier at **relative degree 2**
   (force channel: dpt→V̇/γ̇→sink), uniform with the elevator's degree-2 moment channel → a hard
   impact row enforced by BOTH actuators. `test_beaver_lon.cpp` proves it with the exact Lie jet
   (degree-1 both zero, degree-2 both nonzero) and the complementary power/elevator authority scaling
   (power ~1/V grows as the elevator's ~ρV² fades near touchdown).
**Why:** validated plant + a mixed-actuator hard-safety row that survives the low-q touchdown where
the elevator-only impact barrier is weakest (the TODO "controllability guard" gap).
**Follow-ups / notes for collaborator:** Beaver is a 2288 kg aircraft (vs AHAB 3.6 kg) — wiring it
in as the default plant is a re-parameterization (trim ~33-35 m/s inside LR-556's 30-55 m/s validity
band; the near-stall flare/touchdown extrapolates below 30 m/s — a stated limitation), plus
re-tuning the CBF and bounding `u_P` to the real engine spool rate. Justification:
`water_landing_cbf_math.md` §5, `water_landing_cbf_design.md` §6, `paper_readiness.md` §6.
**Files touched:** `include/autoland/{beaver_aero,beaver_lon,impact_barrier,lon_cbf_filter,sixdof_sim}.hpp`,
`src/{lon_cbf_filter,lon_sim,sixdof_sim}.cpp`, `test/{test_beaver_aero,test_beaver_lon,test_lon_cbf}.cpp`,
`CMakeLists.txt`, `documentation/{paper_readiness,impact_load_barrier_spec,water_landing_cbf_math,water_landing_cbf_design,CHANGELOG}.md`.

## 2026-07-22 — Brian — 6-DOF straight-in landing sim (nonlinear plant + cascaded-PID nominal)
**Branch/commit:** 6dof
**What changed:** New 6-DOF closed-loop water-landing sim (`sixdof_autoland_sim`), replacing the
legacy linear-plant `autoland_sim` path as the 6-DOF workflow. The plant is the existing nonlinear
body-axis EOM extended with an exact wind-aware overload `Dynamics::xdot(x, u, W_earth)`: the
state velocities keep their inertial meaning and the wind enters only through the aero/thrust
(air-relative velocity), so zero wind is bit-identical to the old path and no `Wdot` forcing is
needed. The MIL-F-8785C gust model gained a lateral `v` axis (`v_amp`/`v_len`, + = toward +y east;
step limit = steady crosswind); the Airy/JONSWAP wave field is consumed unchanged (touchdown at
`h = eta(x, t)`, TN 1516 flat- vs wave-referenced slam-load truth in the touchdown record). The
nominal is a cascaded-PID successive-loop-closure design (Beard & McLain 2012): airspeed→throttle
PI, gamma→theta→elevator cascade (pitch inner is PID — `Ki_theta` closes the DC gap of the fixed
`de_trim` feedforward), a `Kv_gamma` speed→path reference shift (one-line TECS; see below), 
cross-track→bank→aileron, yaw damper. No CBF filter and no flare/decrab yet — straight-in only;
the aircraft crabs into crosswind (~9.5 deg at 3 m/s) and holds centerline within ~0.4 m.
**Bug found & fixed on the way: the body-axis EOM's CFx frame sign was inverted.** `kFrameSign`
in `dynamics.cpp` carried `+1` for CFx — the ONE coefficient the six physical sign constraints
(Cm_alpha<0, lift up, Cl_beta<0, Cn_beta>0, Cy_beta<0, Cn_r<0) never pin — so the deck's axial
force entered the standard axes backwards: NEGATIVE total drag (~−0.03 net "suction thrust"
canceling `parasite_CD0` exactly around approach alphas) and an absurdly over-curved polar above
(implied Oswald e ≈ 0.15). It hid because near-zero drag only made trim thrust ≈ 0. Diagnosed by
comparing both sign conventions against the raw deck: flipped, the polar is textbook (min CD_deck
≈ 0.016 near zero lift) and matches the lon path's VSPAERO wind-axis transform (which was always
correct — the lon/CBF results are untouched). Fix: CFx → −1 in `kFrameSign`; total CD is now
~0.045–0.07 across the envelope, idle glide ≈ −6°, and the −3° approach is properly powered.
The bug surfaced as a speed runaway: a hot entry (+2 m/s) accelerated to a ~29 m/s touchdown
because there was no drag path back to V_ref and any speed term inside theta_cmd is exactly
canceled by the gamma integrator at steady state. `Kv_gamma` (the gamma-REFERENCE shift, giving
Vdot ≈ −g·Kv_gamma·(V−V_ref) independent of drag slope) is kept as idle-rail robustness, plus a
hot-entry regression test.
**Why:** The CBF work to date is longitudinal-only; a working 6-DOF dynamic simulation against the
same wind/wave models is the base for lateral-directional landing work (crabbed/wing-down
touchdowns are called out in TODO as entirely outside the lon barriers).
**Follow-ups / notes for collaborator:** Wind and waves each toggle with one `enabled:` line in
the scenario YAML. Roll gains are ceiling-limited by dt=0.01 (roll authority ~460 rad/s² per rad
of aileron → `Kp_p ≥ ~0.4` limit-cycles); pitch gains stay Iyy-scaled per the lon tuning. The
cross-track loop is PD-only, so a small steady y-offset (~0.4 m, from the deck's tiny lateral
asymmetry at beta=0) remains — add `Ki_y` if it matters. Demo scenarios:
`data/sixdof_scenario.yaml` (calm), `data/sixdof_crosswind.yaml`, and
`data/sixdof_landing_waves_lake.yaml`; plot with `scripts/plot_sixdof_results.py`. Default deck
is `AHAB_combined_betasym.stab` (full ±20° beta grid — the plain combined deck clamps at beta<0).
**Files touched:** `include/autoland/{dynamics,wind_gust,sixdof_nominal,sixdof_sim}.hpp`,
`src/{dynamics,sixdof_sim}.cpp`, `apps/sixdof_autoland_sim.cpp`, `data/sixdof_*.yaml`,
`scripts/plot_sixdof_results.py`, `test/test_sixdof_sim.cpp`, `CMakeLists.txt`, `README.md`.

## 2026-07-18 — Brian — Surface-wave model (STANAG 4194 / JONSWAP), plant-side truth
**Branch/commit:** stall-recovery-suite
**What changed:** Linear Airy surface-wave field (`include/autoland/water_waves.hpp`, header-only,
mirrors `wind_gust.hpp`): one JONSWAP spectrum implementation `S(w) = A_γ S_B(w) γ^b` whose
`γ = 1` limit IS the two-parameter Bretschneider spectrum NATO STANAG 4194 prescribes with its
Annex D open-ocean sea states, and whose `γ ≈ 3.3` is the fetch-limited (developing/LAKE) sea —
one code path, the sea is a scenario parameter. Realization per the textbook chain (St. Denis &
Pierson 1953; Fossen ch. 8; ITTC 7.5-02-07-01.1): N seeded components, `A_i = sqrt(2 S(w_i) dw)`,
random phase + within-bin frequency jitter (non-periodic), deep-water dispersion `k = w²/g`;
`regular: true` degenerates to one deterministic component. `eta/slope/etaDot` are exact closed
forms; `slopeMean(x,t,L)` is the tilt across the hull contact length (`contact_len`, default
0.4 m) — the point slope's spectrum `k²S` is ripple-dominated and would overstate the contact
geometry a keel actually bridges.
**Wiring (wave-BLIND filter, per design decision):** plant truth only — touchdown test is now
`h <= eta(x,t)`; the altimeter measures RADAR clearance `h - eta` (so the filter flies the
surface-relative altitude but its barrier model stays flat-water); new earth-frame ground track
`x` (trapezoid outside RK4 — waves never force the airborne dynamics). `LonTouchdown` gains the
contact record: `eta/slope/sink_rel` (closure `-d/dt(h-eta)`) and the exact TN 1516 load truth
`n_peak_flat` vs `n_peak_wave` (tau/gamma0 tilted by the hull-mean surface angle, closure onto
the moving face; via the new unfrozen `impactNPeakExact`, factored out of the barrier factory).
CSV gains `x,eta,eta_slope`. Scenarios: `data/lon_landing_waves_lake.yaml` (JONSWAP γ3.3,
Hs 0.22 / Tp 1.8 from the USACE CEM fetch laws at U10=8 m/s, F=3 km; head seas; STANAG SS3
ocean alternative in comments), `data/lon_landing_wave_regular.yaml` (locked ctest twin).
`scripts/plot_wave_landing.py` → `figures/lon_landing_waves.png`.
**Why:** TODO's top physical-realism gap ("smooth-water assumption"). Source preference was
mil-spec > textbook > paper: STANAG 4194 + USACE CEM parameterize the sea, the textbook layer
realizes it, Hasselmann/Pierson-Moskowitz/Bretschneider are the roots.
**Results:** 57/57 ctest (4 new `[waves]` cases: spectrum closed form + `m0 = Hs²/16` identity,
realization variance/dispersion/determinism, regular-wave analytics incl. the advection identity,
and the end-to-end surface touchdown). Lake demo: the wave-blind filter still meets its
flat-water spec at contact (sink 0.209 m/s, V 13.37 — near-identical to the flat baseline) but
the TRUTH is contact on a 7.7° rising face at 2.36 m/s surface closure: wave-referenced
`n_peak` = 18.8 g vs 0.32 g flat-referenced (~59×, ~6× over `n_limit=3`). The smooth-water gap
is now a measured number, and it motivates the stage-2 wave-aware barrier (TODO updated).
**Caveats:** one seeded sea = one draw (sweep seeds/encounter phase before quoting statistics);
`n_peak_wave` uses the clamped TN 1516 model well outside its validated attitude range when the
face is steep (τ_w clamps at 1°) — read it as "far outside the flat-water design point", not as
a calibrated load; deep-water dispersion assumed (fine for lake chop; shallow-water k would need
the full dispersion relation).
**Files touched:** `include/autoland/water_waves.hpp` (new), `include/autoland/impact_barrier.hpp`
(`impactK0`/`impactNPeakExact` factored out), `include/autoland/lon_sim.hpp`, `src/lon_sim.cpp`,
`test/test_water_waves.cpp` (new), `CMakeLists.txt`, `data/lon_scenario.yaml` (waves block, off),
`data/lon_landing_waves_lake.yaml` (new), `data/lon_landing_wave_regular.yaml` (new),
`scripts/plot_wave_landing.py` (new), `README.md`, `TODO.md`.

## 2026-07-01 — Brian — Emergent stall-recovery test suite (Python driver + ctest layer)
**Branch/commit:** stall-recovery-cbf
**What changed:** End-to-end suite testing the core claim of the recovery-CBF design: on slow
approaches that stall near the water, the CBF-QP recovers and lands safely even though the nominal
(constant-γ + constant-T_set, no flare/recovery logic) never commands it. Three layers:
- **`LonRunStats`** (`lon_sim.hpp/.cpp`, additive): psi minima per barrier over their ACTIVE
  windows, `min_b_impact_active`, `min_res_impact_active` (enforced impact row with the APPLIED
  control), envelope extrema (max α, min V, T range), `out_of_model` (|α| > 85° = stall-table
  validity ceiling), `ic_trim_converged`; exposed via `LonSim::stats()` so tests assert without
  CSV parsing. NOTE: drift-only psi2 goes legitimately negative when the QP spends authority —
  the state-only invariance witness is psi1; enforcement is checked via the row residual.
- **Locked scenarios** `data/lon_stall_pull_cbf.yaml` (nominal pulls into stall at 18 m, CBF-on
  counterpart of the departure demo; T_set 1.2 — at 2.5 the α-capped mush CLIMBS) and
  `data/lon_stall_entry_cbf.yaml` (+11° kick at 12 m → α₀≈19.8°, in-grid) +
  **`test/test_stall_recovery.cpp`** (3 ctest cases: pull recovers @ sink 0.23 / θ +7.7°; entry
  recovers @ sink 0.21 / θ +3.0°; CBF-off demo departs — guards the "emergent" premise).
- **`scripts/stall_recovery_suite.py`** (reuses harness.py): 7 fixed cases × CBF-on(scored)/
  CBF-off(informational A/B) — slow_bleed, pull_into_stall, stalled_entry, tailwind_shear (+4.5
  m/s gust hitting 13 m, 2-pass t_start calibration), downdraft_flare, 2 noisy-sensor variants ×
  3 seeds — plus a recovery-floor sweep (h0 {5..20} × kick {8,12,16}°). PASS/FAIL verdicts, exit
  1 on failure; `--only/--skip-sweep/--replot/--strict/--set` (regression injection). Artifacts →
  `runs/stall_suite/` (summary.md/.csv, per-case A/B PNGs, sweep heatmap).
**Why:** Nothing exercised the filter's recovery behavior end-to-end (TODO: impact barrier
"unverified off-nominal"); the suite makes the emergent-safety claim testable and rerunnable.
**Results (current tuning):** 11/11 scored runs PASS; 53/53 ctest. Recovery floor: clean recovery
from stalled entry down to **12 m** for all kicks; at 8 m and 5 m the outcome stays soft
(sink ≤ 0.25 m/s) but the QP passes through best-effort feasibility recovery (enforced-row
guarantee breaks transiently) → classified "recovered-degraded". Known gaps encoded in
thresholds: (1) unmeasured downdraft breaks the still-air V_td cap by ~1.3–1.6 m/s (energy row's
ḣ model is air-relative) — downdraft_flare carries V_td_tol=2.0 as the frozen baseline;
wind-aware margins are future work. (2) With sensor noise the enforced-row residual on the TRUE
state can dip at the touchdown step (filter still measures +6 cm) — residual check applies to
noise-free runs only. (3) stalled_entry survives `cbf.stall=false` (nominal θ-loop damps the
kick; ensemble lands it) — the stall row's regression detector is pull_into_stall, which departs
to α=133° → INVALID without it.
**Follow-ups / notes for collaborator:** impact row still never binds (see TODO note); sweep
kicks 12/16° extrapolate the VSPAERO deck above +20° α at entry (Viterna blend dominates);
`--strict` flags CBF-off runs that stop failing (scenario drift).
**Files touched:** `include/autoland/lon_sim.hpp`, `src/lon_sim.cpp`, `test/test_stall_recovery.cpp`
(new), `CMakeLists.txt`, `data/lon_stall_pull_cbf.yaml` (new), `data/lon_stall_entry_cbf.yaml`
(new), `scripts/stall_recovery_suite.py` (new), `TODO.md`.

## 2026-07-01 — Brian — MIL-F-8785C discrete wind gust model ("1-cosine")
**Branch/commit:** stall-recovery-cbf
**What changed:** Added the MIL-F-8785C discrete gust — the Simulink Aerospace Blockset
"Discrete Wind Gust Model" block — as a plant-side disturbance for the longitudinal sim
(`include/autoland/wind_gust.hpp`, header-only). Per axis V_wind = Vm/2·(1−cos(πx/dm)) over
penetration distance x, with ẋ = V (airspeed) gated on `t_start`, exactly like the block's
airspeed-integrator input; dm ≤ 0 degenerates to a step gust. Two axes in the vertical plane,
**earth frame, h-up**: `u_amp` + = tailwind, `w_amp` + = updraft (a MIL body-frame w_g, positive
down, is the negation). Plumbed into `lon_sim`: new `lonXdotFullWind` keeps the state air-relative
(V, γ air-path ⇒ α = θ−γ stays the true AoA) and adds the exact wind terms — ḣ += W_h (kinematic
transport) and the −Ẇ inertial forcing projected on path axes (steady wind only transports; only
the gust *ramp* forces V/γ). x_gust is integrated jointly with X in RK4. **The nominal controller
and CBF-QP keep the still-air model** — a gust run probes the filter against unmeasured wind (the
TODO's "gust/shear injection"). YAML `wind:` block in `lon_scenario.yaml` (disabled by default;
wind-off runs keep the old rk4 path, bit-identical). CSV gains `W_u,W_h,x_gust` columns; the
`sink` column is now the *true inertial* −ḣ (= air-relative sink + W_h; unchanged wind-off).
**Why:** Robustness testing of the recovery-CBF set demands a standard, citable gust disturbance.
**Follow-ups / notes for collaborator:** Demo (tailwind 3 m/s + downdraft 1.5 m/s at t=20 s,
dm=120 m, `figures/gust_demo_compare.png`): touchdown 31.3 s @ sink 1.27 m/s / V 15.45 m/s vs
still-air 68.5 s @ 0.21 m/s / 13.37 m/s — V_td cap EXCEEDED because the filter can't see the wind.
Wind-aware barrier margins (ISSf/tightened class-K) are the natural next step. MIL-8785C ties
(Vm, dm) to turbulence intensity/scale; here they're free signed parameters like the block's.
Removed the deprecated `test_cbf.cpp` case "CBF runs on the real linearized longitudinal model"
(it exercised the old body-axis LinearModelCAM seam against `data/AHAB_sweep.stab`, a file no
longer in the repo, so it could never run); full suite now 46/46 green.
**Files touched:** `include/autoland/wind_gust.hpp` (new), `include/autoland/lon_sim.hpp`,
`src/lon_sim.cpp`, `test/test_wind_gust.cpp` (new), `test/test_cbf.cpp`, `CMakeLists.txt`,
`data/lon_scenario.yaml`, `figures/gust_demo_compare.png`.

## 2026-07-01 — Brian — Recovery-set CBFs: impact + stall(AoA) + nose-up + total-energy
**Branch/commit:** stall-recovery-cbf
**What changed:** Replaced the CBF-QP barrier set. **Removed** the descent-rate, lower-airspeed,
and upper-airspeed barriers (structs, config, YAML, wiring, CSV diagnostics, tests). **Kept** the
hydrodynamic impact-load barrier — now the *only* hard safety row (`impact_hard` default flipped to
true) — and the min/max-thrust actuator barriers (kept hard, reframed as the actuator-effectiveness
/ HOCBF-validity guards: every elevator-driven barrier's control coefficient scales with dynamic
pressure ∝ρV², so the thrust guards preserve the relative-degree structure). There is deliberately
**no airspeed-floor CBF** — the AoA barrier handles stall directly. **Added** three barriers, all
soft (slacked), with priority impact(hard) > stall > energy > nose-up encoded by slack weights
1e5 > 1e4 > 1e3:
- `StallBarrier` b = α_lim − (θ−γ), α_lim = α_stall − margin (rel. deg 2, elevator). The stall
  guard; as α→α_lim the HOCBF drives the elevator nose-down — the pilot low-altitude stall recovery.
- `NoseUpBarrier` b = θ − θ_min, gated to the final `h_noseup` metres (rel. deg 2). Holds θ above
  the keel so the impact model's τ=θ−τ_keel>0 gate stays valid and the touchdown is nose-up; lowest
  weight → yields to the stall guard on the shared elevator.
- `EnergyBarrier` b = ½(V_td_max²−V²) + (g_eff−g)h (rel. deg 3, both controls). A **loose
  never-exceed** total-specific-energy ceiling ⇔ descending airspeed cap
  V ≤ √(V_td_max²+2(g_eff−g)h); bounds touchdown speed. `V_td_max` is the hull/structural touchdown
  limit; `g_eff` sized so the ceiling starts satisfied. A *constant* energy cap is useless here (PE→KE
  makes a glider arrive fast), hence the height schedule.
All three are functions of *instantaneous* state only, so the set is **glide-path agnostic** — no
barrier code change when the nominal later switches to the exponential-altitude flare.
**Why:** Recover the pilot behavior at a low-altitude stall (pitch down + full elevator): the impact
CBF regulates vertical touchdown rate, the energy CBF bounds horizontal/kinetic energy at touchdown,
and the AoA barrier prevents stall while keeping the impact model valid.
**Verification:** clean approach lands soft nose-up (sink 0.21 m/s, θ=+3°, V=13.4 < cap), impact
(hard) satisfied throughout, 0 hard rows dropped. Feasible tight touchdown (`V_td_max=10`) → 0.03 m/s
at V=11, θ=+7.3° (guard protective, α peaks 7.4°). The guard only *binds* and forces a hard nose-down
landing when a **sub-stall** touchdown speed is demanded (`V_td_max=9` < V_stall≈9.7 → 3.18 m/s, α
pinned <9°, i.e. no actual stall) — the honest insufficient-margin tradeoff. Operating rule:
`V_td_max ≳ 1.1·V_stall` (~11) for a soft landing. `ctest` 42/43 (the 1 failure is **pre-existing**:
legacy `test_cbf.cpp` needs missing `data/AHAB_sweep.stab`). Rewrote `test/test_lon_cbf.cpp`: closed
forms, relative degrees, finite-difference flow-oracle drift stacks, and hard-enforcement for all
three new barriers.
**Follow-ups / notes for collaborator:** the nominal is still the constant-γ cascade (no flare yet),
so the current soft touchdown comes from the impact barrier flaring near the water; when the
exponential-altitude nominal lands, retune only `V_td_max`/`g_eff`/`theta_min`. `water_landing_cbf_math.md`
is updated for the new set; `water_landing_cbf_design.md` still describes the old barriers — not yet updated.
**Files touched:** `include/autoland/{hocbf,lon_cbf_filter,impact_barrier}.hpp`,
`src/{lon_cbf_filter,lon_sim}.cpp`, `data/lon_scenario.yaml`, `test/test_lon_cbf.cpp`,
`scripts/{harness,plot_lon_results,plot_lon_zoom,plot_gain_compare,plot_lon_sensor_panels}.py`

## 2026-06-30 — Brian — Add NACA 4414 viscous-stall plant model
**Branch/commit:** stall-model
**What changed:** The longitudinal plant can now **stall**. The base deck
(`AHAB_combined.stab`) is inviscid VSPAERO (lift linear to ±20°, no stall, and its
post-stall values are meaningless), so above stall onset the plant **hands the lift off**
from VSPAERO to a viscous + flat-plate post-stall curve via a blend:
`C_plant = (1−w)·C_vspaero + w·C_post`. The blend weight `w(α)` smoothsteps 0→1 across the
wing stall, so VSPAERO is **discarded** (weighted to zero) above onset. `C_post` is a
**Viterna & Corrigan flat-plate extrapolation** (standard BEM post-stall method) anchored at
the wing stall point, so the lift **craters off CLmax and declines to 0 at 90°** while drag
rises to `CD_max ≈ 1.23`. The 4414 **2D viscous polar** (NeuralFoil, an Xfoil surrogate, Re
~2.5×10⁵) sets the realistic CLmax level. New files: generators
`scripts/{gen_4414_polar,precompute_stall_table}.py` → committed table
`include/autoland/naca4414_stall_table.hpp` (`w/CLpost/CDpost/CMpost(α)` to 90°); lookup
`include/autoland/stall_model.hpp` (`stallLookup`, value+slope, mirrors `clfLookup`); blend
frozen in `makeAeroLocal` and applied in `LonDrift` (`src/lon_augmented.cpp`,
`lon_augmented.hpp`). Config: `AircraftConfig::stall {enabled,severity}` (`config.{hpp,cpp}`,
`aircraft.yaml`, default **OFF**), with an optional scenario-level `stall:` override
(`lon_sim.cpp`). Diagnostics: CSV gains `CL,CD,dCL_stall`; `scripts/plot_stall_model.py` →
`figures/stall_model_check.png`. Demo scenario `data/lon_stall_recovery.yaml`.
**Why:** Stall-recovery control work needs a plant that actually stalls. The inviscid deck
cannot — it lifts linearly forever — so there is nothing to recover from. VSPAERO's only
"stall" option is a CLmax clamp (a plateau, not a departure), and its inviscid lift past stall
is fiction; we need the real lift crater + drag rise + pitch break, with VSPAERO thrown out
above onset.
**Design choices:**
- **Blend, not additive deltas.** `(1−w)·C_vspaero + w·C_post` discards VSPAERO above onset
  (the user's point: VSPAERO doesn't model stall at all, so blending *toward* its post-stall
  values is wrong). An earlier additive `ΔCL = CL_vspaero·(f−1)` formulation was scrapped: it
  kept multiplying the (fictional) rising inviscid lift, and the deep-stall tail plateaued
  instead of cratering.
- **Viterna flat-plate tail.** Standard post-stall extrapolation: lift follows `sin 2α` to 0 at
  90°, drag → `CD_max = 1.11+0.018·AR`. Replaces the non-physical plateau; deep stall now
  craters.
- **Frozen-affine, shared plant+CBF.** Lives in `makeAeroLocal`, frozen as local-affine
  `{off + slope·α}` for `w` and each post-coefficient, so it reuses the exact-Lie path; the
  plant re-freezes each RK4 substep. `w=0` in attached flow ⇒ the CBF is untouched.
- **Stall point is a tunable knob, not truth.** NeuralFoil→Xfoil→reality is weakest at low-Re
  post-stall; onset/level are generator knobs (`A_STALL_DEG`, `BLEND_HALF_DEG`, `SEVERITY`).
  Defaults: CLmax ≈ **1.44 @ 11°** (≤ the 2D section CLmax 1.46), → 0.80 by 45°, 0 at 90°.
- **OFF by default**, the disabled `LonDrift` path runs zero extra ops. Valid to ~90° α (table
  clamps above; a tumble past 90° is outside the model).
**Verified:** Nominal approach **bit-identical** stall on vs off (`max|Δh|=max|Δδe|=0`,
`dCL_stall≡0`, max α≈2°). 5 new `[stall]` tests pass (lookup inert/held/interp; post-stall
lift-decline-to-zero + drag-rise; blend hands CL/CD to the Viterna values and cuts lift).
End-to-end `lon_stall_recovery` departs: climbs to CLmax≈1.42 @ α≈11°, then lift craters
(0.90 @ 28°, 0.45 @ 69°, 0 @ 97°), CD→1.23, sink builds to ~10 m/s. Full suite: only the
**pre-existing** `#19` fails (missing `AHAB_sweep.stab`, unrelated).
**Follow-ups / notes for collaborator:** calibrate onset/level vs experimental low-Re 4414 or
flight ID; optional 2nd table axis over Re; optional Goman–Krabrov separation **state** for
dynamic-stall lag (doubles as a CBF state); wire an α-margin recovery barrier; optionally link
the descent barrier's `CL_max` knob to the plant CLmax (~1.44). CSV gained 3 columns —
name-keyed plot scripts unaffected.
**Files touched:** `include/autoland/{naca4414_stall_table,stall_model,lon_augmented,config}.hpp`,
`src/{lon_augmented,config,lon_sim}.cpp`, `scripts/{gen_4414_polar,precompute_stall_table,
plot_stall_model}.py`, `data/{aircraft.yaml,lon_stall_recovery.yaml,naca4414_polar.csv}`,
`test/test_stall_model.cpp`, `CMakeLists.txt`,
`documentation/{stall_model_spec,water_landing_cbf_design,CHANGELOG}.md`, `README.md`, `TODO.md`.

## 2026-06-28 — Brian — Add hydrodynamic impact-load HOCBF barrier (NACA TN 1516)
**Branch/commit:** water-impact-cbf
**What changed:** New independent CBF-QP row that bounds the **peak CG load factor at
water touchdown**, derived from NACA TN 1516 (Milwitzky 1948). The barrier is
`b = (n_limit − n_peak(θ,γ,V)) + Φ(z)` with `n_peak = K0·ẏ₀²·Clf(κ)` (peak load factor,
g), the approach parameter `κ = sin τ/sin γ₀ · cos(τ+γ₀)`, and a height-relaxed term
`Φ(z) = Nb(1−e^{−z/zs})` that makes it touchdown-only. New files: an offline precompute
`scripts/precompute_impact_clf.py` (solves the transcendental eq 27 for the velocity
ratio, evaluates eq 25, verifies the `Clf(0)=0.6123` anchor) → generated table
`include/autoland/impact_clf_table.hpp`; the barrier itself in
`include/autoland/impact_barrier.hpp` (`clfLookup`, `ImpactLoadBarrier`,
`makeImpactLoadBarrier`). Wired into `LonCBFFilter` behind an `impact`/`impact_hard`
flag set with new config (`n_limit`, `beta`, `rho_water`, `Nb`, `zs`, `tau_keel`,
`z_gate`, `eps_g0`, `impact_slack_lo/hi`, `c_impact`) in `LonCBFConfig` +
`lon_scenario.yaml` + the `lon_sim` loader. Added an `exp()` overload to the Taylor-jet
Lie engine (`lie_taylor.hpp`) for `Φ`. The sim logs `b_impact,n_peak,kappa_imp` and the
degree-2 `psi1/psi2_imp` + `res_imp`, and prints the impact ψ-minima.
**Why:** Bound the structural *slam load* at hull contact — which depends on dead rise,
trim, and sink rate, not just sink rate (the descent barrier's kinematic cap). This is
the attitude-coupled contact-load constraint flagged as the §3.3 research item.
**Design choices:**
- **Degree-2 HOCBF, elevator-enforced.** The barrier is relative degree 2 via the
  elevator (θ→q→δe) and 3 via thrust (T→Tdot→Tddot). A single affine row can't mix the
  two (the elevator's degree-2 entry pushed to degree 3 yields non-affine u²/u̇ terms),
  so it's built as a clean degree-2 row (`barrierLie<2>` + 2 class-K gains): thrust's
  column in `L_gL_f b` is exactly 0 and drops out; the flare enforces it. Thrust still
  bounds impact load via the existing descent/sink-rate barrier (spec §4).
- **Frozen K0 + local-affine `Clf(κ)`** at the eval point (mirrors `makeDescentBarrier`),
  so the templated barrier is smooth with no cbrt/pow in the Taylor path while the
  attitude coupling stays live through `κ`.
- **Activation:** Option A height-relaxed `Φ(z)` (primary) + Option C height-scheduled
  slack (`impact_slack_lo→hi`) backstop. Plus a **model-validity gate** — the row is only
  assembled below `z_gate` AND while descending with positive trim (NACA TN 1516's
  validity domain); outside it the prediction is meaningless (κ<0) and would feed a
  huge-coefficient row to the QP.
- **Placeholders:** `n_limit=3 g`, `beta=22.5°`, `rho_water=1000` (fresh), `Nb=10`,
  `zs=2 m`. Soft + enabled by default and sized **non-binding** (additive safety, zero
  trajectory change), mirroring how the upper-airspeed barrier was introduced.
**Verified:** **34/35 tests pass** — incl. 6 new (clfLookup monotone/clamp, closed-form
value, relative-degree thrust-drops-out, drift flow-oracle for the degree-2 stack +
`exp()`, Φ height-relaxation, hard-enforcement under a violating nominal). The 1
remaining failure (`#19`, `test_cbf.cpp`) is **pre-existing** (old body-axis path, throws
unrelated to this work). Default sim run **unchanged**: touchdown 0.0137 m/s, 0
recoveries; `impact: min psi1=5.84, psi2=14.11 ≥ 0`; over the active window
`n_peak ≤ 0.17 g « n_limit=3`, `κ∈[0.001,1.149]`, slack ≡ 0 (row never binds).
**Follow-ups / notes for collaborator:** `n_limit`/`beta`/`Nb`/`zs`/`c_impact` need
calibration/tuning (new `TODO.md` items). The full mixed-degree (2/3) treatment and a
predictive/backup-set CBF remain open (§3.3). CSV gained 6 columns — name-keyed plot
scripts are unaffected. Source paper added at `documentation/19930082553.pdf`.
**Files touched:** `include/autoland/{lie_taylor,impact_barrier,impact_clf_table,
lon_cbf_filter}.hpp`, `src/{lon_cbf_filter,lon_sim}.cpp`, `scripts/precompute_impact_clf.py`,
`data/lon_scenario.yaml`, `test/test_lon_cbf.cpp`,
`documentation/{impact_load_barrier_spec,water_landing_cbf_design,CHANGELOG}.md`, `TODO.md`.

## 2026-06-27 — Jack — Add maximum-airspeed (over-speed) HOCBF barrier

**Branch/commit:** corbin-dev
**What changed:** Wired the previously-stubbed `AirspeedUpperBarrier` (`b = V_max − V`) into the
CBF-QP filter — a symmetric over-speed / high-energy-impact guard alongside the stall barrier.
New config: `airspeed_upper` / `airspeed_upper_hard` flags, `Vmax_air`, `c_airspeed_upper`
(`LonCBFConfig` + `lon_scenario.yaml` + loader). The sim now logs `b_airspeed_upper`, its ψ-minima
(console), and a `res_airup` residual column. Because the barrier is a function of `V` only, it
reuses the lower barrier's degree-3 machinery unchanged: its drift Lie stack and control row are
the **sign-flip** of the stall barrier's (`L_f^k b_up = −L_f^k b_lo` for `k≥1`).
**Why:** Bound high-energy water impact / structural limits — slam loads grow with impact speed
on a flying-boat hull. Closes the `TODO.md` "Upper airspeed barrier" item.
**Design choices:** `V_max` is a placeholder (`27` m/s ≈ 1.5·V_app) set generously so it is
protective but **non-binding** on the current approach; soft by default (mirrors the stall
barrier); enabled by default (additive safety, zero trajectory change).
**Verified:** all **29 tests** pass — incl. 2 new (oracle SECTION revalidating the upper
barrier's drift stack; a sign-flip identity test vs the lower barrier; a hard-enforcement test
under an over-speed nominal). Default sim run **unchanged** (touchdown 0.0137 m/s, 0 recoveries);
`b_airspeed_upper` stays in [8.67, 10.48] (V_max=27), upper ψ-minima ≈ 17/33 — confirmed
non-binding. Spot-check: lowering `V_max` / setting it hard makes the row actively cap `V`.
**Follow-ups / notes for collaborator:** `V_max` needs calibration (new `TODO.md` item, under
Modeling/data); both airspeed barriers remain soft (the "harden the airspeed barrier" item now
applies to the pair). CSV gained a column — name-keyed plot scripts are unaffected.
**Files touched:** `include/autoland/{hocbf,lon_cbf_filter}.hpp`, `src/{lon_cbf_filter,lon_sim}.cpp`,
`data/lon_scenario.yaml`, `test/test_lon_cbf.cpp`,
`documentation/{water_landing_cbf_math,water_landing_cbf_design,CHANGELOG}.md`, `TODO.md`.

## 2026-06-27 — Jack — Steady level-flight trim for the initial condition

**Branch/commit:** corbin-dev
**What changed:** Added `lonTrim()` (`src/lon_sim.cpp`) — a Newton solve on the lon EOM for a
steady **level-flight** equilibrium (`γ=0`, `V̇=γ̇=q̇=0`, mirroring the body-axis `trim()`). The
sim IC now uses it, so the run starts cruising level and the nominal pushes over into the −3°
approach (more realistic; no body-axis-seed startup transient).
**Why:** Asked for a realistic IC; also a step on the "true longitudinal trim" backlog item.
**Finding (corrects an earlier note):** the level IC did **not** change the residual descent
`ψ < 0` — so that dip is **not** the initial condition. It is a **terminal** effect in the last
~0.4 s as `h → 0` at touchdown; the actual barrier `b` stays ≥ 0 in flight (the realized
trajectory is safe). Design doc §6/§8.4 updated accordingly.
**Verified:** all 27 tests pass; default run lands at 0.014 m/s within budget, 0 recoveries,
`b_descent ≥ 0.10` airborne. (Also generated a class-K gain comparison figure showing the four
`b` barriers; stiff `[100,100,100]` gains under-flare and violate `b_descent` at touchdown,
gentle `[1,1,1]` stay safe.)
**Files touched:** `src/lon_sim.cpp`, `documentation/water_landing_cbf_design.md`, `TODO.md`,
`documentation/CHANGELOG.md`.

## 2026-06-27 — Jack — Force-based `a_brk(V,γ)` for the descent barrier; real deck by default

**Branch/commit:** corbin-dev
**What changed:**
- Replaced the constant `a_brk = 3.0` with a speed/path-angle-dependent
  `a_brk(V,γ) = (ρV²S/2m)[CL_max cosγ − CD_maxlift sinγ] − g` (lift + drag + gravity) in
  `DescentBarrier`, built by `makeDescentBarrier()` and differentiated exactly by autodiff.
  **Thrust deliberately omitted** — it would drop the barrier's relative degree below 3 and
  break the augmented HOCBF alignment (backlog, tied to §3.3).
- `C_L,max` is a new config field (`CL_max`, placeholder 1.2); `C_D,maxlift` is computed from the
  frozen aero at the extrapolated max-lift α. Added a positivity warning in the filter.
- Pointed `lon_autoland_sim`'s default aero deck at the **real `AHAB_combined.stab`** (was the
  `example.stab` placeholder — the source of earlier "violations").
- Tests: oracle test now auto-revalidates the new barrier's drift stack; replaced the descent
  closed-form authority check (invalid for state-dependent `a_brk`) with an `a_brk(V,γ)` formula
  value test. **All 27 pass.**
**Why:** Make the soft-landing envelope honestly track available lift (tighten at low speed,
loosen at high) instead of a hand-picked constant.
**Verified (real AHAB deck):** touchdown sink **0.011 m/s** (within the 0.1 budget), **0 QP
recoveries** — a big improvement over the constant-`a_brk` baseline (sink 2.18 m/s, 11
recoveries on the same deck). The residual descent `ψ < 0` at the start is the non-trim initial
condition, **not** control saturation (smaller class-K gains made it *worse*).
**Follow-ups / notes for collaborator:** `CL_max` needs calibration; thrust term + lag-aware
`a_brk` + a true lon trim are in `TODO.md`. The unit tests still use `example.stab` (deck-
agnostic math); removing `example.stab` is a backlog item.
**Files touched:** `include/autoland/{hocbf,lon_cbf_filter}.hpp`,
`src/{lon_cbf_filter,lon_sim}.cpp`, `apps/lon_autoland_sim.cpp`, `data/lon_scenario.yaml`,
`test/test_lon_cbf.cpp`, `documentation/{water_landing_cbf_math,water_landing_cbf_design}.md`,
`TODO.md`, `documentation/CHANGELOG.md`.

## 2026-06-27 — Jack — Add `TODO.md` backlog; retire stale README TODO

**Branch/commit:** corbin-dev
**What changed:** Added a root-level `TODO.md` (granular backlog, grouped by area) seeded from
the design doc's open items. Replaced the stale README "⚠️ TODO" section (placeholder-data
items resolved by the AHAB refactor) with a pointer to it. Wired `TODO.md` into the
collaboration convention (CLAUDE.md startup note + changelog doc-map).
**Why:** Give a low-friction, collaborator-visible place to capture tasks/ideas, distinct from
the changelog (done) and the design doc §9 (rationale for big threads).
**Files touched:** `TODO.md`, `README.md`, `CLAUDE.md` (gitignored), `documentation/CHANGELOG.md`.

## 2026-06-27 — Jack — Verify drift Lie stack (`L_f³b`) + make build MSVC-portable

**Branch/commit:** corbin-dev
**What changed:**
- Added a test-only **finite-difference flow oracle** that independently cross-checks the
  drift Lie stack `{b, L_f b, L_f²b, L_f³b}` (descent + airspeed barriers) against the
  bespoke Taylor-jet engine in `lie_taylor.hpp`. Engine vs oracle agree to ~2.5e-4 on
  `L_f³b`; a seeded 1% engine error is caught. Closes the prior "L_f³b not independently
  verified" gap (only the control authorities were checked before).
- Made the build work under **MSVC** (the only compiler on this Windows machine): guarded the
  GCC/Clang `-Wall -Wextra` flags behind `if(NOT MSVC)` and defined `_USE_MATH_DEFINES` for
  `M_PI`. Strictly additive — does not change the GCC/Clang build.
**Why:** Cheap insurance that the high-order Lie machinery is correct (the autodiff *primitive*
was never in doubt; the custom Picard/Taylor *iteration* was the untested part). Build fix was
needed to compile/run the suite locally on Windows.
**Verified:** all 26 tests pass (`ctest -C Release`); confirmed the new test fails on a 1%
perturbation, then reverted.
**Follow-ups / notes for collaborator:** to build on Windows/MSVC, configure with the Visual
Studio generator (`cmake -S . -B build -G "Visual Studio 17 2022" -A x64`) — Ninja hits an
OSQP `osqp.lib` multiple-rules error.
**Files touched:** `test/test_lon_cbf.cpp`, `CMakeLists.txt`,
`documentation/water_landing_cbf_design.md`, `documentation/CHANGELOG.md`.

## 2026-06-27 — Jack — Set up `corbin-dev` branch + collaboration docs

**Branch/commit:** corbin-dev (off `lon-cbf-water-landing` @ `2b610a3`)
**What changed:**
- Created `corbin-dev` from the latest `lon-cbf-water-landing` (the augmented-longitudinal
  HOCBF water-landing module with exact autodiff/Taylor-jet Lie derivatives).
- Gitignored personal Claude context (`CLAUDE.md`, `PROJECT_CONTEXT.md`).
- Added this changelog and a living design doc (`water_landing_cbf_design.md`) that supersedes
  the dated implementation notes.
- Moved historical docs (`water_landing_cbf.md`, `water_landing_cbf_implementation_notes.md`)
  into `documentation/archive/`.
**Why:** Establish a shared working branch and a lightweight way for two parallel
Claude-using collaborators to stay in sync, plus a single maintained source of truth for the
design instead of a point-in-time session snapshot.
**Follow-ups / notes for collaborator:**
- Add the one-line changelog pointer to your own local `CLAUDE.md` (see Collaborator setup).
- Open design questions captured in `water_landing_cbf_design.md`: force-based `a_brk(V)`,
  the §3.3 contact-force barrier, and an independent `L_f³b` drift check.
**Files touched:** `.gitignore`, `documentation/CHANGELOG.md`,
`documentation/water_landing_cbf_design.md`, `documentation/archive/`.
