# PX4 TECS port — what was ported, how it is wired, what it does on the Beaver

**Source:** PX4-Autopilot `src/lib/tecs/TECS.hpp` + `TECS.cpp`, `main` @
`a906b72868d79c89790688e5b41de73597133ce5` (fetched 2026-08-27; upstream author Paul
Riseborough). Parameter defaults from
`src/modules/fw_lateral_longitudinal_control/fw_lat_long_params.yaml` at the same commit; the
call-site conventions (pitch offset, load factor, zeroed airspeed-rate input) from
`FwLateralLongitudinalControl.cpp`; the test fixture from `TECSTest.cpp`.

**Port:** `include/autoland/px4_tecs.hpp`, `src/px4_tecs.cpp` (`autoland::px4::TecsControl`).
Double precision, upstream's function decomposition and member names, `[PX4]` comments are
upstream's own. Diff it against `TECS.cpp` method by method.

## Scope

| Upstream piece | Status | Note |
|---|---|---|
| `TECSControl` (altitude/airspeed outer loops, SPE/SKE rates, SEB → pitch, STE → throttle, integrators + anti-windup, underspeed detection, speed/altitude weighting, pitch-rate limit from `vert_accel_limit`, throttle slew, bank-angle drag compensation, `AlphaFilter` on the STE-rate estimate) | **ported** | the control law |
| `TECS` wrapper (uORB, `hrt` timestamps, `DT_MIN/DT_MAX` resets) | dropped | the sim calls `initialize()` once, then `update(dt)` |
| `TECSAltitudeReferenceModel` (jerk-limited `VelocitySmoothing` trajectory on the altitude setpoint) | dropped | the sim commands the **direct height-rate** setpoint; the altitude loop is ported but bypassed |
| `TECSAirspeedFilter` (2-state fixed-gain complementary filter; upstream feeds it a zero acceleration — `HOTFIX` in the module) | dropped | the sim supplies the **exact** airspeed rate (`airspeedRate()` in `sixdof_sim.cpp`) |
| fast-descend mode | dropped | its hooks inside `TECSControl` are written out at `fast_descend = 0` |
| airspeed-less branches | kept structurally | `Flag::airspeed_enabled` is always true here |
| NaN guards / `PX4_ISFINITE` | dropped | double-precision sim inputs are finite by construction |

## Wiring (`sixdof_nominal.hpp`, `sixdof_sim.cpp`)

* `nominal.type: cascade | tecs` selects the longitudinal outer loops. The pitch PID
  (`Kp_theta / Ki_theta / Kq`) and the lateral axes are shared, so the two nominals differ only in
  the outer energy loops. Default `cascade`; every pre-existing CSV column of the cascade runs is
  bit-identical to before (checked on calm / crosswind / POH).
* Setpoints: `altitude_rate_setpoint_direct = V_ref sin γ_ref`, `tas_setpoint = V_ref`.
* Inputs: `altitude = h`, `altitude_rate = ḣ` (state kinematics), `tas = V_air` (pitot),
  `tas_rate = V̇_air` exact — `v_a · (v̇_b − Ẇ_b)/|v_a|` with `v̇_b` the plant `xdot` under the
  control held over the last step (ZOH; no algebraic loop) and `Ẇ_b = Rᵀ Ẇ_e − ω × W_b`.
* Pitch: TECS returns pitch **above the level-flight trim pitch** (PX4 subtracts `FW_PSP_OFF`);
  `theta_cmd = tecs_pitch_offset + pitch_sp`, limits `±theta_cmd_max` about the offset.
* Load factor per step: `1 / max(cos φ, FLT_EPSILON)` (upstream expression).
* First call runs `initialize()` (as PX4 does): it has no direct-rate path, its demand is the
  altitude loop's `HRATE_FF × ḣ`, so there is a one-sample throttle/pitch kick at `t = 0`.

## Vehicle anchors — solved from the plant, not parameters

PX4 takes these as parameters scaled by `PerformanceModel`; here they are the plant's own
equilibria (`deriveTecsAnchors()`), so the throttle feedforward's three anchors and the pitch
offset are exact for the airframe:

| TECS param | PX4 source | Here | Beaver, 40 m/s, 1800 RPM |
|---|---|---|---|
| `throttle_trim` | `FW_THR_TRIM` (0.6) | level trim at `V_ref` | 0.832 |
| pitch offset | `FW_PSP_OFF` | level-trim θ at `V_ref` | 8.82° (approach trim θ 6.18°) |
| `max_climb_rate` | `FW_T_CLMB_MAX` (5) | steady climb at `dT_max`, `V_ref` | 1.085 m/s |
| `min_sink_rate` | `FW_T_SINK_MIN` (2) | steady sink at `dT_min`, `V_ref` | 4.30 m/s |
| `max_sink_rate` | `FW_T_SINK_MAX` (5) | steady sink at `dT_min`, `tas_max` | 5.84 m/s (55 m/s) |
| `tas_min / tas_max` | `FW_AIRSPD_MIN/MAX` | LR-556 validity band | 30 / 55 m/s |
| `pitch_min / max` | `FW_P_LIM_MIN/MAX − FW_PSP_OFF` | `∓theta_cmd_max` | ±10° |
| `throttle_min / max` | `FW_THR_MIN/MAX` | `SurfaceLimits.dT_min/max` | 0 / 1 |

The climb/sink solves march γ 1° at a time from the approach γ until the 6-axis trim's throttle
crosses the target, then bisect — trim evaluations only, no derivatives. The small
full-throttle climb rate is the scenario's engine setting (1800 RPM / 26 "Hg), not a TECS
matter.

## Parameter mapping (`tecs:` YAML block)

| YAML key | PX4 parameter | PX4 default | Beaver default |
|---|---|---|---|
| `alt_tc` | `FW_T_ALT_TC` | 5 s | 5 s (altitude loop bypassed) |
| `hrate_ff` | `FW_T_HRATE_FF` | 0.5 | 0.5 |
| `tas_tc` | `FW_T_TAS_TC` | 5 s | 5 s |
| `spdweight` | `FW_T_SPDWEIGHT` | 1.0 | 1.0 |
| `ptch_damp` | `FW_T_PTCH_DAMP` | 0.1 | **1.0** |
| `i_gain_pit` | `FW_T_I_GAIN_PIT` | 0.1 | **0.4** |
| `thr_damping` | `FW_T_THR_DAMPING` | 0.05 | 0.05 |
| `thr_integ` | `FW_T_THR_INTEG` | 0.02 | **0.3** |
| `thr_slew_max` | `FW_THR_SLEW_MAX` | 0 (off) | 0 |
| `vert_acc` | `FW_T_VERT_ACC` | 7 m/s² | 7 |
| `rll2thr` | `FW_T_RLL2THR` | 15 | 15 |
| `ste_r_tc` | `FW_T_STE_R_TC` | 0.4 s | 0.4 |
| `seb_r_ff` | `FW_T_SEB_R_FF` | 1.0 | 1.0 |
| `tas_error_percentage` | (hard-coded) | 0.15 | 0.15 |
| `detect_underspeed` | (module flag) | on | on |
| `thr_trim`, `pitch_offset_deg`, `clmb_max`, `sink_min`, `sink_max`, `tas_min`, `tas_max`, `pitch_min_deg`, `pitch_max_deg` | anchors above | — | plant-derived |

## Tuning (`scripts/tune_tecs.py`) — three gains, with a hold-out set

**Why tune at all.** At PX4's flown defaults the Beaver touches down 5 % steep (sink 2.57 vs
2.44 m/s, γ −3.68°) and rings for ~15 s after a hot entry. Two structural reasons: the predicted
throttle is linear between (idle, −4.30 m/s) and (level, 0) and gives 0.36 where the true −3.5°
approach throttle is 0.409, and the pitch feedforward `γ_sp` about the *level*-trim pitch lands
0.9° below the descent trim pitch — both gaps close on the integrators, which at `THR_INTEG
0.02` / `I_GAIN_PIT 0.1` take ~100 s. The ringing is the SEB loop being integrator-dominated:
`PTCH_DAMP 0.1` is a proportional gain of only 0.1 rad pitch per rad of flight-path error
(`Δθ = damp·Δḣ/V`), so the PI zero sits at 4 rad/s, far above the ~0.5 rad/s loop crossover.

**Protocol (anti-overfitting).**
* Training cases — four longitudinal disturbances at the design condition (40 m/s clean,
  1800 RPM, h0 60 m): calm, hot entry (+3 m/s, +2°), cold entry (−3 m/s, −2°), tailwind shear
  (+5 m/s over 120 m at t = 5 s).
* Validation cases — held out, never used to pick gains: POH flaps-35 float configuration at
  33.5 m/s / 2000 RPM (calm and hot — a different trim, throttle map and pitch offset),
  downdraft (−3 m/s), 5 m/s crosswind (load-factor coupling), 120 m approach, 45 m/s / −3°,
  and a −6° approach with the throttle on its idle rail (anti-windup).
* Score per run `J = IAE(sink − sink_ref) + 0.5·IAE(V − V_ref) [m] + 0.02·TV(δe) [deg]`
  (t = 0 sample excluded), normalised by the PX4-default TECS on the same case (`Jrel`, 1 = the
  starting point); the cascade is shown for reference only (its exact approach-trim feedforward
  makes its calm-case J ≈ 0, which would swamp a ratio). Hard checks alongside: max sink, min
  V, max α, max per-step elevator increment (limit-cycle detector), touchdown errors.
* Search: coarse grids inside PX4's parameter ranges over `ptch_damp`, `i_gain_pit`,
  `thr_damping`, `thr_integ` only; everything else at the PX4 default. Pass 1 (0.1–2 / 0.1–0.8
  / 0.05–1 / 0.05–0.5, 384 combos) put the optimum on the `i_gain_pit` and `thr_integ` edges;
  pass 2 (0.5–2 / 0.4–2 / 0.2–1 / 0.3–1, 240 combos) extended them to the PX4 maxima. The pick
  is from the flat region, judged by a *robust* score (worst over the combo and its ±1 grid
  neighbours), not the single best; a gain whose effect is small stays at the PX4 default.

**Grid marginals (pass 2, mean `Jrel` over the training cases at each gain value).**

| gain | values → mean Jrel |
|---|---|
| `ptch_damp` | 0.5: 1.21 · 1.0: 0.85 · 1.5: 0.79 · 2.0: 0.78 — steep below 1, flat above |
| `i_gain_pit` | 0.4: 0.67 · 0.8: 0.70 · 1.2: 0.78 · 1.6: 0.98 · 2.0: 1.42 — interior optimum |
| `thr_damping` | 0.2: 0.91 · 0.5: 0.90 · 1.0: 0.92 — insensitive |
| `thr_integ` | 0.3: 0.91 · 0.5: 0.91 · 0.75: 0.91 · 1.0: 0.92 — insensitive above 0.3 |

Raw best 0.624 (`1 / 0.4 / 1 / 1`), robust pick 0.635 (`1.5 / 0.4 / 1 / 1`), conservative
corner of the same region 0.654 (`1 / 0.4 / 0.5 / 0.3`) — a 5 % spread. One-at-a-time
sensitivity at the conservative corner: `ptch_damp` 0.75 / 1.25 → 0.663 / 0.656; `i_gain_pit`
0.2 / 0.3 / 0.6 → 0.666 / 0.656 / 0.655; `thr_damping` 0.05 → 0.685 (+5 % train, +1 %
validation); `tas_tc` 3 / 8 → 0.612 / 0.712; `vert_acc` 3 / 10 → 0.805 / 0.640; `spdweight`
0.5 / 1.5 → 0.720 / 0.686 (0.5 cuts the shear sink excursion 5.6 → 4.6 m/s at the cost of
airspeed tracking — a landing-design choice, left at 1); `ste_r_tc` and `seb_r_ff`: flat / worse.

**Decision — three gains change, the rest stay at PX4 defaults:**

| gain | PX4 | Beaver | reading |
|---|---|---|---|
| `FW_T_PTCH_DAMP` | 0.1 | **1.0** | pitch per unit flight-path error; the cascade's hand-tuned `Kp_gamma` is 1.5 |
| `FW_T_I_GAIN_PIT` | 0.1 | **0.4** | pitch rate per unit flight-path error; the cascade's `Ki_gamma` is 0.3 |
| `FW_T_THR_INTEG` | 0.02 | **0.3** | trims the throttle-map mismatch in ~10 s instead of ~100 s |
| `FW_T_THR_DAMPING` | 0.05 | 0.05 | insensitive (0.2–1.0: +1 % on validation) → unchanged |

The pitch gains are the cascade's outer-loop gains expressed in energy form — two independent
routes (hand tuning of the cascade in July, this grid) landing on the same numbers is the best
evidence the result is a property of the airframe rather than of the training set.

**Validation (held out; `Jrel` vs the PX4-default TECS on each case).**

| setting | train mean / worst | validation mean / worst | validation max sink | validation max δe step |
|---|---|---|---|---|
| PX4 defaults | 1.00 / 1.00 | 1.00 / 1.00 | 5.48 m/s | 0.24° |
| integrators only (0.1 / 0.4 / 0.05 / 0.3) | 0.99 / 1.36 | 0.57 / 0.98 | 4.78 | 0.24° |
| **tuned (1.0 / 0.4 / 0.05 / 0.3)** | **0.69 / 1.02** | **0.40 / 0.73** | 4.63 | 0.27° |
| raw best (1.0 / 0.4 / 1.0 / 1.0) | 0.62 / 0.98 | 0.36 / 0.73 | 4.64 | 0.28° |
| cascaded PID (reference) | 0.54 / 0.80 | 0.24 / 0.80 | 4.19 | 0.44° |

Every held-out case improves over the PX4 defaults, including the flaps-35 configuration
(different anchors), the throttle-railed −6° approach (no windup: touchdown sink 4.20 vs the
4.19 reference) and the downdraft, where TECS beats the cascade outright (touchdown sink 2.57 vs
3.11 m/s — the energy loop sees the loss that an airspeed PI cannot). The elevator carries no
high-frequency content (dominant 0.04–0.17 Hz; the largest per-step increment is the t = 0
init kick). The cascade still wins on the pure-longitudinal cases in absolute terms because its
feedforward is the exact approach trim; TECS's is anchored on level flight by construction.

## Results (Beaver, calm and hot entry; `figures/sixdof_tecs_compare_{calm,hot}.png`)

| Case | Nominal | sink [m/s] | γ [°] | V [m/s] | settle (\|sink−ref\| < 0.15) | max sink |
|---|---|---|---|---|---|---|
| calm (nominal 2.442 / −3.5 / 40) | cascade | 2.443 | −3.50 | 40.00 | 0 s | 2.45 |
| calm | TECS, PX4 defaults | 2.567 | −3.68 | 39.97 | 17.8 s | 2.97 |
| calm | TECS, tuned | 2.443 | −3.50 | 40.01 | 4.4 s | 2.76 |
| hot entry (+3 m/s, +2°) | cascade | 2.497 | −3.59 | 39.89 | 18.7 s | 2.70 |
| hot entry | TECS, PX4 defaults | 2.646 | −3.78 | 40.13 | 23.2 s | 4.13 |
| hot entry | TECS, tuned | 2.457 | −3.51 | 40.08 | 9.3 s | 3.21 |

The lateral axes are untouched (φ, y, ψ within the cascade's numbers in every case).

## Tests (`test/test_tecs.cpp`)

Unit checks run on PX4's own fixture numbers (`TECSTest.cpp makeParam()`): setpoint hold /
no integrator drift; energy signs (total energy → throttle, balance → pitch); altitude-loop
and SPE/SKE formulas; `vert_accel_limit` pitch-rate limit and anti-windup at the pitch limit;
underspeed ramp (ratio, full throttle, speed priority, frozen throttle integrator); throttle
slew; direct height-rate bypass; load-factor correction. Closed loop: Beaver calm to the
cascade's tolerances, and the hot entry settling on `V_app`.
