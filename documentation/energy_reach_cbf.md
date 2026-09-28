# Reach-the-Target Energy CBF (3-DOF formulation, evaluated on the 6-DOF Beaver)

*Implements the reformulated energy CBF from the EnergyCBF slide deck (2026-09): the
alpha-dependent glide ratio replaces the constant (L/D)max, and a heading-dependent
circular-arc distance replaces the straight NE line. Written 2026-09-06.*

Design decisions fixed with the author up front: the sideslip command is realized with a
**rudder sideslip loop** (wings level — faithful to the model's skid-to-turn mechanism);
evaluation covers **both idle-glide and powered** scenarios; the target is a **YAML point
defaulting to the glideslope aim point**; only the **new formulation** is implemented (no
constant-L/D baseline).

## 1. Formulation

$$h_E = \left(\tfrac12 V^2 + g z\right) - \left(\tfrac12 V_{min}^2 + \frac{g\,d}{L/D}\right)$$

* **Glide ratio** $L/D = \dfrac{C_{L\alpha}\,\tilde\alpha}{C_{D0} + C_{D\alpha 2}\tilde\alpha^2}$
  with $\tilde\alpha = \alpha - \alpha_{0L}$ (alpha from the zero-lift line — the Beaver is
  cambered, so the literal slide form $C_L = C_{L\alpha}\alpha$ needs the shift).
* **Arc distance** $d = R\,\varphi/\sin\varphi$, $\varphi = \chi - \mathrm{atan2}(\Delta_E,
  \Delta_N)$ wrapped via $\varphi = \mathrm{atan2}(d_\perp, d_\parallel)$ with
  $d_\parallel = \Delta_N\cos\chi + \Delta_E\sin\chi$, $d_\perp = \Delta_N\sin\chi -
  \Delta_E\cos\chi$ (slide 7's $d_\parallel, d_\perp$ appear naturally).
* **Model** (3-DOF point mass, unpowered, skid-to-turn): state
  $X = [V, \gamma, \chi, N, E, z, \alpha]$, controls $u = (\dot\alpha, \beta)$,
  $\dot\chi = \bar q S C_{Y\beta}\beta/(mV\cos\gamma)$. Uniform relative degree 1;
  constraint $L_f h + L_g h\,u \ge -c_k h$.

Smooth (autodiff-safe) guards: sqrt-form smooth-max floor on $\tilde\alpha$; tanh
saturation of $\varphi$ at `phi_sat` (< π, caps the arc at ~5R for a target behind);
series form of $\varphi/\sin\varphi$ below 0.05 rad; straight-line fallback inside
`d_near` of the target (the bearing is ill-defined there).

## 2. Implementation map

| Piece | Where |
|---|---|
| Taylor-jet `atan2` / `tanh` + `scalarValue` branch helper | `include/autoland/lie_taylor.hpp` |
| Point-mass drift/g-matrix, barrier, exact RD-1 Lie bundle | `include/autoland/energy_reach_cbf.hpp` |
| CBF-QP filter + idle-glide polar fit | `src/energy_reach_cbf.cpp` |
| Nominal outer/inner split + rudder sideslip inner loop | `include/autoland/sixdof_nominal.hpp` |
| Sim wiring (`energy_cbf:` block, idle-gamma solve, CSV, stats) | `src/sixdof_sim.cpp` |
| Tests (fit, Lie cross-check, QP behavior, invariance, 6-DOF smoke) | `test/test_energy_reach.cpp` |
| Scenarios | `data/beaver_reach_{glide,offset,powered}[_off].yaml` |
| Plots | `scripts/plot_energy_reach.py` → `figures/reach_*.png` |

Architecture in the 6-DOF loop: `stepOuter()` (cascade/TECS) → map `theta_cmd` to
`alpha_des = theta_cmd − gamma` → rate-shape to $\dot\alpha_{nom}$ through `tau_alpha`,
$\beta_{nom} = 0$ → `EnergyReachFilter` QP at the measured point-mass state → integrate
`alpha_cmd` → `stepInnerBeta(x, alpha_cmd + gamma, dT, beta_meas, beta_trim + beta*, dt)`.
Throttle forced to idle in glide mode (`gamma_app` is replaced by the plant's own
idle-glide slope at `V_app`, solved by trim bisection, so trim/reference/model agree).

The Lie derivatives come from the existing flow-Taylor-jet engine (`lieDrift<1>` /
`lieAlong<1>`) — the slide-7 hand algebra was **not** transcribed; the unit test instead
cross-checks the jet against an independently hand-derived closed-form gradient (1e-8).

**QP** (z = [\dot\alpha, β]): range-scaled tracking weights (1/u_max²); AoA band folded
into the $\dot\alpha$ bounds (each bound is itself an RD-1 CBF row); the energy row is
**hard**, normalized by max(|A|) (raw $A_{u_\alpha}\sim g\,d\,\Psi'\sim 10^4$ makes OSQP
emit spurious infeasibility certificates otherwise); infeasibility against the input box
is detected in closed form first and falls back to the **max-ḣ rail** (min-violation
action). A soft slack was tried and rejected — pricing the *true* violation needs
w·s² ≈ 1e12 in P (kills conditioning); pricing the *normalized* violation makes the row
soft as mush (h_E sailed 300+ below zero at "slack ≈ 0.001").

**Polar fit**: least squares of the wind-axis rotation of the Beaver body polynomials over
α ∈ [−3°, 14°] at zero rates/deflections and the **idle** dpt, so `C_D0` absorbs the
windmilling-propeller drag (P ≈ −97 kW at pz = 5″Hg!). At 40 m/s:
CLa = 5.45/rad, CD0 = 0.078, CDa2 = 0.85, (L/D)max = 10.6 at α* ≈ 17.3° — i.e. **best
glide sits above the 14° trusted-AoA ceiling**; in-band, L/D increases monotonically with
α and the CBF rides the ceiling when energy-critical. Note the Beaver glides fast and
nose-high at idle: flaps-up trim α is ~10.5° at 40 m/s.

## 3. Results (all in `runs/beaver_reach_*.csv`, figures in `figures/`)

| Scenario | min h_E | miss | verdict |
|---|---|---|---|
| glide, diving nominal, filter ON | **+20.9** | 247 m | dive refused, phugoid suppressed, invariant held |
| glide, diving nominal, filter OFF | −270 | 272 m | dive + phugoid, barrier violated mid-flight |
| 10° offset, filter ON (ck = 0.01) | **+0.02** | 249 m | early skid turn, rides h_E = 0 exactly |
| 10° offset, filter OFF | −460 | 453 m | never turns, passes abeam |
| powered −5.8°, filter ON | +56 → +265 | 42 m | conservative under thrust, quiet |

## 4. What the exercise taught (read before the 6-DOF port)

*(Forward-looking design discussion for the next step — flipping this barrier into an
energy CEILING and pairing it with the impact CBF for a gentle touchdown within a
region — lives in `energy_ceiling_touchdown_notes.md`.)*

1. **h_E is conserved along any on-polar glide toward the target.** The barrier prices the
   remaining flight at the *current* α's L/D — the same L/D the drift flies — so steady
   gliding (at ANY α) and lossless zoom/dive trades cost nothing. Only genuine losses
   drive it down: flying at worse L/D during transients, heading error, unmodeled drag,
   wind. Consequences:
   * The class-K gain is the **"how early" dial**: h_E ~ O(10³) m²/s² while its natural
     loss rates are O(1–10) m²/s³, so ck = 0.3 lets the barrier sleep until the endgame
     (it then rails everything and still loses). ck = 0.01 binds early and wins. The
     scenario YAMLs pin ck = 0.01; the header default stays 0.3 for the unit tests.
   * The barrier is **necessary, not sufficient**, for arrival: landing short with
     V > V_min is "safe" per h_E ≥ 0 (the kinetic reserve counts as spendable but the
     ground interrupts the trade). A terminal z–V coupling (or a proper target-tracking
     nominal) is future work.
2. **Skid-to-turn authority is marginal at km scale**: χ̇_max = q̄SC_{Yβ}β_max/(mV) ≈ 2°/s
   at β = 10° → a ~1.2 km minimum "turn radius". Late binding + weak authority is why
   large offsets (≥ 15°) are unrecoverable; 10° with early binding works. The 6-DOF port
   should either go bank-to-turn or encode the turn-radius limit in d (the deck's Dubins
   alternative does exactly this).
3. **The plant pays sideslip drag the model omits** (wind-frame drag from C_{Yβ}β·β ≈ 25%
   extra at β = 10°); the barrier does not price it, so hard skid phases bleed real h_E.
4. The QP uses the **α channel first, β last** (range-scaled weights + |A_ua| ≫ |A_be|);
   β engages when the AoA ceiling saturates. Visible in `figures/reach_offset.png`.
5. **The rudder sideslip loop works**: flown β tracks an 11.5° railed command within ~1°
   (`Kp_beta = 2, Ki_beta = 0.5`, yaw damper retained, wings level via the roll loop).
6. Best-effort (max-ḣ rail) engages exactly where the input-bounded QP is genuinely
   infeasible — a plain RD-1 CBF-QP is myopic and input bounds void the invariance
   guarantee there; report `cbf_best_effort_steps` honestly.

## 5. Verification

* `Beaver idle-glide polar fit is physical` — fit sanity + rms bounds.
* `Reach Lie derivatives match an independent hand derivation` — jet engine vs closed
  form (h to 1e-10; L_f h, A_ua, A_be to 1e-8), plus sign semantics.
* `Filter passes the nominal through when the margin is large` / `Binding filter raises
  alpha toward best glide and steers at the target` — QP behavior at both regimes.
* `Closed-loop 3-DOF invariance under an adversarial nominal` — full-rate dive nominal;
  h_E must respect the class-K decay envelope h₀e^{−ck t} and never dip below −0.5% h₀
  over the feasible stretch.
* `Positive rudder raises Beaver sideslip (beta-loop sign)` — exact linearization.
* `6-DOF idle-glide smoke run with the energy-reach layer` — end-to-end wiring.
