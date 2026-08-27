# DHC-2 Beaver 6-DOF plant — implementation validation

**Status (2026-08-27): validated, including full-envelope 6-DOF.** The Beaver 6-DOF
plant (`beaver_dynamics.hpp`, direct evaluation of the verified LR-556/FDC polynomials —
no table export) reproduces the published FDC 1.2 references at every level we can
check, matches an independent implementation over the whole flight envelope (large
rates, bank, turning flight — §5), and flies the closed-loop 6-DOF landing cases
inside the model's validity band.

Reproduce everything with:

```
build/beaver_validation results/beaver          # checks 1-3 + linearization dumps
python3 scripts/validate_beaver_modes.py results/beaver   # check 4 (exits nonzero on fail)
python3 scripts/validate_beaver_sixdof.py build results    # check 5 (exits nonzero on fail)
python3 scripts/plot_beaver_validation.py results/beaver results/beaver_validation.png
build/sixdof_autoland_sim "" "" data/beaver_landing_<case>.yaml results/beaver_<case>.csv
```

Unit-test versions of checks 1, 2, 3, the 30° turn equilibrium (check 5), and the
calm landing live in `test/test_beaver_dynamics.cpp` (part of the 77-test suite).

## The four validation levels

### 1. Coefficient level (pre-existing)

`test_beaver_aero.cpp`: the ~60 polynomial coefficients were verified against LR-556
Table 3 / FDC manual Table C.3–C.4 to every printed digit, and the aero+engine chain
reproduces an independent Python port's oracle values.

### 2. Full nonlinear ẋ — the FDC ACTRIM check case (external oracle)

The FDC 1.2 manual (figs. 10.18–10.19) prints a complete trimmed state, input vector,
and state derivative for **V = 35 m/s, n = 1800 RPM, pz = 20 ″Hg**. Evaluating our
6-DOF `xdot` at exactly that printed (x, u):

| row | ours | FDC printed | diff |
|---|---|---|---|
| V̇ [m/s²] | −1.898e−4 | −1.887e−4 | 1.1e−6 |
| α̇ [rad/s] | −1.514e−5 | −1.235e−5 | 2.8e−6 |
| β̇ [rad/s] | 4.669e−4 | 4.636e−4 | 3.3e−6 |
| ṗ [rad/s²] | −2.497e−5 | −2.503e−5 | 5.8e−8 |
| q̇ [rad/s²] | −2.152e−5 | −2.066e−5 | 8.6e−7 |
| ṙ [rad/s²] | −5.009e−5 | −5.260e−5 | 2.5e−6 |
| ẏe [m/s] | −0.72329 | −0.72330 | 6.5e−6 |
| Ḣ [m/s] | −0.67916 | −0.67909 | 7.2e−5 |

Every row matches to the printout's own 5-digit precision (forces are O(g), so this is
~10⁻⁶ relative). This validates aero + engine + atmosphere + gravity + EOM assembly +
kinematics end-to-end against an independent published implementation. The check also
*discriminates*: with the wrong (2000 ft) atmosphere the same rows miss by 3–5 orders
of magnitude — see `figures/beaver_validation_checkcase.png`.

**Two reference-condition discoveries** (both matter if anyone re-derives oracles):

* **The check case is at the sea-level atmosphere.** ACTRIM was *prompted* with
  2000 ft, but the printed trim state carries H = 0 and the printed accelerations are
  only reproduced with ρ = 1.225 (at ρ(2000 ft) the residuals are ~0.1 m/s², at ρ(0)
  they are ~3e−6): the altitude prompt seeds the simulation IC, not the trim
  atmosphere. (The GitHub `Flight_Simulator` port that supplied the older
  `test_beaver_aero.cpp` oracle used its own ρ(2000 ft) ≈ 1.138 — self-consistent
  for the coefficient chain, but not the FDC trim condition.)
* **The engine altitude correction multiplies only the (408 − 0.0965 n) term** in
  FDC eq. 3.15 — the manual's typesetting is ambiguous, but the other reading gives a
  constant −240 W at sea level for any manifold pressure, which is absurd (and the
  check case confirms our reading numerically).

Also confirmed en route: FDC's atmosphere is the exact ICAO troposphere
(T₀ = 288.15 K, λ = −0.0065, exponent g₀/λR) and gravity is g₀·(Rₑ/(Rₑ+h))² —
both now implemented in `beaver_dynamics.hpp`.

### 3. Trim solve — recovery of the FDC trim point + the fig. 10.13 elevator curve

Our 6-axis Newton trim (unknowns α, β, δe, δa, δr, throttle; exact autodiff Jacobian)
at the check case's γ = −1.112°:

| quantity | ours | FDC | |
|---|---|---|---|
| α | 12.105° | 12.107° | 0.002° |
| θ | 10.993° | 10.995° | 0.002° |
| δe | −5.3337° | −5.3333° | 0.0004° |
| pz | 19.998 ″Hg | 20 ″Hg | 0.002 |
| β | −1.016° | −1.184° | 0.17° |
| δa | 0.463° | 0.551° | 0.09° |
| δr | −2.819° | −2.836° | 0.02° |

The longitudinal trim recovers FDC to print precision. The lateral deltas are FDC's
own convergence tolerance, not model error: FDC's fmins stopped with a β̇ residual of
4.6e−4 rad/s (its own printout), i.e. its printed lateral trim is not a converged
zero of its own equations; our Newton drives all six residuals to ~1e−14 on the same
model.

The manual's fig. 10.13 elevator curve was **pixel-digitized** from a 400-dpi render
(frame-calibrated axes; `data/fdc_fig1013_solid_digitized.csv`) rather than read by
eye. The figure's exact conditions are unstated ("low power"), but a **fixed
pz = 20 ″Hg** wings-level trim (γ free — apparently TRIMDEMO's constraint style, and
the ACTRIM default power) reproduces the digitized solid curve's shape with a uniform
~0.3° offset over the whole 32–60 m/s range — inside the scan/digitization/unknown-
condition uncertainty — while the exact-condition ACTRIM point is matched to 0.0004°.
Our level-flight trim curve lies just above (more power → more nose-up trim via
Cm_dpt), the published high-power dotted curve above that: the ordering and spacing
are exactly the model's power sensitivity. See `figures/beaver_validation.png` and
the side-by-side against the actual scan, `figures/beaver_validation_fig1013_sidebyside.png`.
The trim power bucket (~99 kW minimum near 33 m/s, rising both ways) is consistent
with the real aircraft.

### 4. Derivative level — independent cross-implementation of the modes

A single-point oracle cannot see rate-derivative errors (p = q = r = 0 there), so
`scripts/validate_beaver_modes.py` implements the same published model **a second
time, independently**: FDC state coordinates (V, α, β, …) instead of body (u, v, w),
coefficients typed in fresh from the manual's tables, complex-step differentiation
instead of autodiff. The full-state Jacobians of the two implementations are similar
matrices, so their eigenvalues must agree:

```
[OK] cruise_45:   trim residual 4.8e-10, max eigenvalue diff 1.2e-09 (rel 2.4e-10)
[OK] check_35:    trim residual 6.2e-10, max eigenvalue diff 1.9e-09 (rel 4.0e-10)
[OK] approach_35: trim residual 1.3e-10, max eigenvalue diff 1.9e-09 (rel 4.0e-10)
```

The modes at 45 m/s cruise are the classic Beaver set: short period −2.16 ± 2.41i
(ωn 3.2 rad/s, ζ 0.67), phugoid −0.016 ± 0.263i (ζ 0.06), roll −5.13, Dutch roll
−0.48 ± 0.97i (ωn 1.1, ζ 0.44), spiral −0.046 (weakly stable).

### 5. Full-envelope 6-DOF — large rates, bank, turning flight (added 2026-08-26)

The four checks above all sit at p = q = r = 0 wings-level, where the quadratic
gyroscopic/Coriolis terms ((Izz−Iyy)qr, Ixz(p²−r²), ω×v) and the large-attitude
kinematics contribute nothing — even the trim-point Jacobians can't see them.
`scripts/validate_beaver_sixdof.py` (figure: `figures/beaver_sixdof_validation.png`)
closes that regime against the same independent implementation:

* **Random-state sweep** — 256 states spanning the envelope (V 28–58 m/s, α to 16°,
  β to ±15°, body rates to ±1.2 rad/s, bank to ±60°, full control deflections,
  flaps 0/20/35°, three altitudes/RPMs). All 11 ẋ rows match the independent
  implementation pointwise to **3.3e−9** normalized (the floor is the C++ smooth-V
  guard, not a modeling difference; pass threshold 1e−8).
* **Steady coordinated turns** — level-turn equilibria at φ = ±15/30/45° solved on
  the *independent* implementation (nonzero p, q, r), with the physics anchors
  checked: n_z tracks 1/cos φ, ψ̇ = g tan φ/V, R = V²/(g tan φ). The C++ plant
  agrees they are equilibria to **1.4e−9** — the gyroscopic terms validated by
  cross-implementation *and* physics. The 30° case is a permanent unit test
  (`test_beaver_dynamics.cpp`, 77 tests total).
* **Open-loop doublet** — a 14 s elevator-doublet + aileron + rudder maneuver
  (q to 0.15 rad/s, φ to 10°, full short-period/Dutch-roll/spiral content)
  integrated by both implementations from the same state/controls: trajectories
  agree to **8e−9 m/s** in V and **3e−8 deg** in attitude over the full maneuver.

With this, every term of the 6-DOF EOM has been exercised by a check that could
fail: the static rows by the external FDC oracle, the first-order dynamics by the
mode cross-check, and the rate-quadratic/attitude terms by the sweep, the turns,
and the doublet. Remaining caveats are *model-validity* (attached flow, 35–55 m/s
band, landplane config), not implementation.

## Closed-loop 6-DOF landing cases

`plant: beaver` is now the `sixdof_autoland_sim` default (the AHAB deck remains
selectable with `plant: vspaero`; its scenarios/tests are unchanged and green). The
cascaded-PID nominal gains were re-sized for the 2288 kg airframe, and the controller
gained per-axis **control-sense signs**: the Beaver uses the standard Delft/FDC
deflection conventions (Cm_δe, Cl_δa, Cn_δr all < 0), opposite to the AHAB virtual
controls on all three axes — with the AHAB signs the aileron loop is positive
feedback and the aircraft rolls inverted in ~6 s. Feedforwards now include the
δa/δr trim (the Beaver trims with β ≈ −0.9°, δr ≈ −2.3° from slipstream asymmetry).

| case (`data/beaver_landing_*.yaml`) | result |
|---|---|
| **calm** (40 m/s, −3.5°, clean) | touchdown exactly on slope: sink 2.443 m/s = V·sin γ, V 40.00, y 0.33 m, φ −0.16° |
| **crosswind** (5 m/s step) | lands y 0.64 m, crab ψ = −7.9° (wind triangle: atan(5/40) ≈ 7.1° + trim offset), β ≈ 0 in the crab |
| **gust** (4.5 m/s tail shear + 2 m/s downdraft on short final) | V sags to 37, throttle briefly firewalls, sink spike arrested; touchdown sink 3.04 m/s, γ −3.94° |
| **waves_lake** (JONSWAP Hs 0.22 m head sea) | lands on η(x,t); wave-referenced TN 1516 per-float load 9.6 g vs 4.2 g flat-referenced — the wave-blind-controller truth |
| **poh** (flaps 35°, 33.5 m/s — the POH float approach) | trims at θ 0.18° (realistic flat float attitude), δe +8.1° against the flap's nose-up Cm_δf, lands y 0.09 m |

Control activity is quiet everywhere (calm-case elevator increments ~2 µrad/step, V
held to ±3 mm/s); α stays ≥5° below stall and V inside the validated band in all
non-POH cases. Per-case plots: `figures/beaver_<case>.png`.

**Logging note:** the sim CSV is written at 10 significant digits. At the C++ default
6, a 40 m/s airspeed quantizes to 0.1 mm/s steps, which rendered the (genuinely
±3 mm/s-tight) airspeed channel as a staircase on autoscaled plots — a logging
artifact, not dynamics.

**Scenario-design note:** clean-configuration 35 m/s trims at α ≈ 12.9°, only ~3°
below stall — that is why the standard cases fly 40 m/s clean (α ≈ 9.7°) and the
POH configuration gets its own flapped scenario at 33.5 m/s (α ≈ 3.7°).

## Stated limitations (unchanged from paper_readiness §6)

* LR-556 Table 3 is flight-validated over **35–55 m/s TAS** (manual Table C.3
  caption), attached flow only; a future flare/touchdown decelerating below ~30 m/s
  extrapolates. The landing cases here touch down at V_app (no flare yet).
* Landplane aero as-is (float drag not modeled) — decided 2026-08-26.
* ρ and g are frozen at the scenario reference altitude (classic constant-density
  convention, same as the AHAB plant); fine for a 60 m approach.
* Deflection limits in `BeaverPlantConfig` are representative round numbers, not
  LR-556 data.

## Follow-ups

* Re-tune / re-derive the lon CBF stack on the Beaver (beaver_lon.hpp exists; the
  6-DOF path runs nominal-only, like the AHAB 6-DOF sim).
* Flare + decrab for realistic touchdown speeds (couples with the CBF work).
* Bound the power-rate control u_P to the real engine spool rate when the CBF
  filter moves onto this plant.
