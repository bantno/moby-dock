# Hydrodynamic Impact-Load Barrier for AHAB Autoland

Implementation reference for a separate control barrier function that bounds the
peak water-impact CG load factor at touchdown, derived from NACA TN 1516. This
document is self-contained and describes the model, the barrier, and the design
decisions. It does not prescribe a particular code structure.

---

## 1. Purpose

Add an independent barrier to the existing CBF-QP landing filter, alongside stall
margin, sink rate, glide slope, and flare, that bounds the peak CG load factor the
hull experiences at water contact. The barrier must enforce the load limit at
touchdown and stay inactive earlier in the flight, where contact acceleration is
irrelevant. Activation is achieved with a height-relaxed barrier (Option A) so the
flare engages from a shrinking altitude budget, backed by a height-scheduled slack
(Option C) for feasibility.

---

## 2. Source paper: NACA TN 1516 (Milwitzky, 1948)

### 2.1 What it models

A theoretical and experimental study of the motion and hydrodynamic loads of a
prismatic V-bottom seaplane during a step-landing impact, from the instant of keel
contact until rebound. The analysis treats the flow about the immersed hull in
transverse planes fixed in space and applies the virtual-mass concept (Wagner) to
each plane, with an end-flow correction for the finite three-dimensional case.
Buoyancy and viscosity are neglected, which holds when inertia forces dominate.

### 2.2 The key result

Nondimensionalizing the equations of motion collapses the entire impact to a single
dimensionless group, the approach parameter

    kappa = sin(tau) / sin(gamma0) * cos(tau + gamma0)            (eq 20)

which depends only on trim `tau` and flight-path angle at contact `gamma0`. For a
given `kappa`, the dimensionless time histories of load-factor coefficient, draft
coefficient, vertical-velocity ratio, and time coefficient are each a single curve,
independent of weight, dead rise, attitude, or initial velocity. The absolute load
factor is recovered through scale factors that carry the hull properties and the
contact velocity.

### 2.3 The quantity we use

The peak CG load factor is set at the instant of maximum acceleration. Combining the
load-factor coefficient (eq 14/15) with the max-acceleration relations (eqs 23a, 25,
27), the peak load reduces to a near-closed form in the contact state:

    n_peak = Cl_max(kappa) * ydot0^2 * (alpha / (W * g^2))^(1/3)

with the hull coefficient (the constant in the transformed equation of motion,
eq 12a)

    alpha  = f_beta^2 * phi_A * rho * pi / (6 * sin(tau) * cos(tau)^2)
    f_beta = pi / (2*beta) - 1                                    (eq 45, beta in rad)
    phi_A  = 1 - tan(tau) / (2*tan(beta))                         (eq 49)

`Cl_max(kappa)` is the peak load-factor coefficient. The paper does not give it in
closed form, since eqs 26 and 27 are transcendental. It is obtained by precomputing
a curve and interpolating (Section 3.1). The closed-form anchor
`Cl_max(0) = 0.6123` (eq 40, with velocity ratio 7/9 at max acceleration, eq 39) is
used to verify the precomputation.

### 2.4 Interpretation that makes this a barrier

`n_peak` is the peak load that would result if contact occurred at the current
`(tau, gamma0, ydot0)`. The impact lasts roughly 0.1 to 0.3 s, and during it the
aerodynamic surfaces have no authority, so the peak is locked in at the moment the
keel touches. The barrier therefore acts before contact and shapes the contact
state. This is why an anticipatory CBF is the right tool and an in-impact CBF is
pointless.

### 2.5 Validity and assumptions

| Assumption | Consequence for the barrier |
|---|---|
| Prismatic float, constant cross section | Valid forward of the step for conventional hulls. Pulled-up bow neglected. |
| Constant trim through the impact | Reasonable for the short contact duration. |
| Wing lift equal to weight during contact | If lift is dumped (spoilers), the load differs. |
| No chine immersion | Conservative. Chine immersion only reduces the peak (eq 42 region), so the prediction is an upper bound on the true load. |
| Smooth water | For waves, define `tau` and `gamma0` relative to the local wave slope (wave-axis approximation). |
| Buoyancy neglected | Small except at very steep flight-path angles beyond normal landings. |
| Dead rise 22.5 to 40 deg validated, ~15 deg plausible | Use AHAB's actual hull dead rise. |
| Practical `kappa` range 0.2 to 10 | Clamp near the planing singularity `gamma0 -> 0`. |

Two further points. The prediction is the rigid-body CG load factor normal to the
water, not the local hull bottom pressure, so `n_limit` must come from the
structural allowable with knockdowns, and the local-pressure check stays separate.
The model singularities are real: `gamma0 -> 0` sends `kappa -> infinity`,
`tau -> 0` sends `alpha -> infinity`, and `phi_A` crosses zero when
`tan(tau) = 2*tan(beta)` (near 40 deg for beta = 22.5 deg, outside landing trims).

---

## 3. The barrier model

### 3.1 Predicted peak load

Evaluating the barrier requires `Cl_max(kappa)`. Precompute it offline. Over a grid
of `kappa` across the practical range, solve the max-acceleration relation (eq 27)
for the velocity ratio at maximum acceleration, evaluate the peak coefficient
(eq 25), and store a monotone interpolant. Verify against the anchor
`Cl_max(0) = 0.6123`. At run time, `n_peak` follows from the contact state
`(tau, gamma0, ydot0)` and the hull coefficients for AHAB's dead rise, mass, and
water density.

### 3.2 Contact-state mapping

In smooth water, `tau = theta` (pitch attitude, minus any keel incidence),
`gamma0 = atan2(ydot, xdot)`, and `ydot0 = ydot`, with `ydot > 0` sink and
`gamma0 > 0` descent. In waves, reference `tau` and `gamma0` to the local wave
slope.

### 3.3 Relative degree

Trim is the pitch attitude, relative degree 2 from elevator. The velocity channel
is relative degree 1 from thrust. The barrier is therefore a high-order CBF with
two layers, and the control enters at the second layer. The height term in
Section 5 keeps it second order, since altitude enters through `zdot = -ydot`. The
constraint row is assembled from the control-affine longitudinal dynamics of the
landing model.

---

## 4. Separate-barrier integration

Add the impact barrier as an independent inequality row in the QP, not by modifying
the sink-rate barrier. The QP is

    min_u  || u - u_des ||_R^2
    s.t.   stall, sink-rate, glide-slope, flare rows
           impact row          <- new, second-layer HOCBF constraint
           u in U

The impact barrier and the sink-rate barrier overlap in intent. Keep both. The
sink-rate barrier is a simple kinematic cap. The impact barrier is the structural
constraint that accounts for dead rise and attitude, and it dominates near the
ground. Treat the sink-rate cap as the looser, always-on guard and the impact
barrier as the sharp terminal constraint.

---

## 5. Activation: making the barrier touchdown-only

The barrier must enforce the load limit at contact and be inactive earlier, since
contact acceleration at altitude does not matter.

### 5.1 Why the naive barrier is always-on

`h = n_limit - n_peak(tau, gamma0, ydot0)` depends only on attitude and velocity,
not altitude. During a normal descent `ydot > 0`, so `n_peak` can exceed `n_limit`
and `h` goes negative at altitude. A CBF enforces forward invariance of `{h >= 0}`
for all time, so the filter would fight the descent the whole way down. Scheduling
the class-K constant `alpha(h)` does not help, because the issue is that `h` itself
is negative, not that it is decreasing too fast.

A forward-integrated predictive or backup-set CBF would remove the need for a height
term, since a predicted touchdown state already includes the planned flare. That is
more rigorous but requires a short forward integration each step. It is heavier than
needed here and is not pursued. The two mechanisms below are the chosen approach.

### 5.2 Option A: height-relaxed barrier (primary)

Introduce height above water `z` (AGL, `z >= 0`, contact at `z = 0`) and relax the
barrier by a nonnegative budget that vanishes at the surface:

    h(x) = [ n_limit - n_peak(tau, gamma0, ydot0) ] + Phi(z),
    Phi(0) = 0,  Phi'(z) >= 0.

Properties:

- At `z = 0`, `Phi = 0`, so `h = n_limit - n_peak`. Maintaining `h >= 0` to the
  surface enforces the true load constraint exactly at contact.
- At altitude, `Phi(z)` is large, `h >= 0` holds trivially, and the barrier is
  inactive. This encodes "do not care about contact load up high."

The activation falls out of the CBF condition. With `zdot = -ydot`,

    hdot = -ndot_peak - Phi'(z) * ydot,

and the CBF condition `hdot >= -alpha(h)` rearranges to

    ndot_peak <= alpha(h) - Phi'(z) * ydot.

Read this as a budget. The allowable rate of increase of predicted load is
`alpha(h)` minus the rate at which the altitude budget shrinks, `Phi'(z) * ydot`.
At altitude, `h` is large, `alpha(h)` is large, and `Phi'` is small, so the bound
is loose and the barrier does nothing. As `z -> 0`, `Phi -> 0`, `h` collapses to the
true margin, `alpha(h)` shrinks, and `Phi'(z) * ydot` forces `ndot_peak` down. The
flare emerges from the shrinking budget, with no hard switch.

Choosing `Phi(z)`. Use a bounded, smooth function of `z` alone:

    Phi(z)  = N_b * (1 - exp(-z / z_s))
    Phi'(z) = (N_b / z_s) * exp(-z / z_s)

`N_b` is the largest counterfactual excess load the prediction may show at altitude
before the barrier reacts. `z_s` is the altitude scale over which the flare
authority operates. `Phi'(0) = N_b / z_s` is the maximum budget-shrink sensitivity
at the surface and must not exceed the rate at which the flare can actually reduce
`n_peak`. Keeping `Phi` a function of `z` alone keeps the barrier gradient simple,
adding only `dh/dz = Phi'(z)`.

Design rule. Pick `N_b` large enough that the barrier is inactive along the nominal
descent profile, and `z_s` so that `Phi'(0)` is matched by the flare authority.
Validate by confirming the barrier stays inactive in cruise and descent and engages
only inside the flare window.

### 5.3 Option C: height-scheduled slack (backstop)

Soften the impact row with a slack `delta >= 0` and a penalty in the QP cost:

    impact row relaxed by delta,   add  rho(z) * delta^2  to the cost.

Schedule `rho(z)` small at altitude (barrier cheap to relax, effectively off) and
large near the ground (barrier strongly enforced). This guarantees QP feasibility
and absorbs a momentary shortfall in flare authority near the surface as a logged,
bounded slack rather than an infeasible solve. It is a soft mechanism with no
invariance guarantee, so it is a backstop to Option A, not a primary activation.
Combined, Option A provides the anticipatory engagement and the formal terminal
constraint, and Option C guarantees the solve stays feasible when authority is
briefly exceeded.

### 5.4 Altitude assumption (current scope)

Assume perfect altitude information for now. `z` is the true vertical height of the
keel above the contact surface, used directly in `Phi(z)` and its gradient with no
biasing, smoothing, or sensor model. In smooth water `z` is the keel height above
mean water level. In rough water, reference `z` to the local wave surface, which is a
modeling choice independent of how `z` is measured. Under this assumption Option C is
a feasibility backstop for momentary flare-authority shortfall only, not for
estimator error.

Deferred. When realistic altitude estimation is reintroduced, the engagement timing
of Option A becomes sensitive to error in `z`. Overestimating `z` inflates `Phi(z)`,
relaxes the barrier too much, and is the dangerous direction, so the barrier should
then consume a lower-bound estimate, `Phi'` should be desensitized to noise through a
larger `z_s` and light bounded smoothing, and any rangefinder beam range should be
converted to vertical height using attitude. Until then this is out of scope.

### 5.5 Recommendation

Use Option A as the primary activation. The barrier engages from the shrinking
altitude budget, which is the natural formulation. Add Option C as a height-scheduled
slack backstop for feasibility, using the true height `z` per Section 5.4. Avoid pure
altitude gating, which is discontinuous. For compute, the
impact row may be skipped above a height where `Phi(z)` is within a small tolerance
of `N_b`, since the row is inactive there. That gate is an efficiency choice, not the
activation mechanism.

---

## 6. Implementation requirements

1. Precompute `Cl_max(kappa)` offline by solving eq 27 for the velocity ratio and
   evaluating eq 25, stored as a monotone interpolant. Verify `Cl_max(0) = 0.6123`.
2. Evaluate `n_peak` from `(tau, gamma0, ydot0)` and the hull coefficients for
   AHAB's beta, mass, and water density.
3. Form `h = (n_limit - n_peak) + Phi(z)` with
   `Phi(z) = N_b (1 - exp(-z / z_s))`. Carry `z` (true vertical height of the keel
   above the contact surface) in the barrier state.
4. Assemble the second-layer HOCBF constraint from the control-affine longitudinal
   dynamics, including the `dPhi/dz` term in the barrier gradient.
5. Add the impact constraint as a separate QP row, with the Option C slack and its
   height-scheduled penalty.
6. Expose tuning knobs: `n_limit`, `N_b`, `z_s`, the two class-K gains, and the
   slack penalty schedule.
7. Clamp `kappa` to the practical range and guard the singularities at
   `gamma0 -> 0`, `tau -> 0`, and `phi_A -> 0`.

Validation in SITL:

- Barrier inactive in cruise and descent, engages only in the flare window.
- Touchdown `n_peak <= n_limit` across a sweep of approaches (vary sink rate, trim,
  flight-path angle, weight).
- `kappa` stays within the practical range through the flare.
- Slack stays at zero on nominal approaches and activates only when flare authority
  is momentarily exceeded.

Tuning order: set `n_limit` from structure, set `N_b` from the worst nominal descent
excess, set `z_s` so the surface budget-shrink sensitivity is matched by flare
authority, tune the class-K gains against the flare timescale, then set the slack
penalty schedule as the backstop.

---

## 7. Caveats

- **Single planing surface (float _or_ hull), full weight — scope, verified against
  TN 1516.** The paper puts the *entire* seaplane weight on one keeled V-bottom ("a
  keeled float or hull") and says nothing about twin floats, load sharing, or
  lateral/asymmetric contact. This barrier inherits that single-surface, full-weight
  scope. On a twin-float airframe (the DHC-2 Beaver plant; see `paper_readiness.md`
  §6) it is used as a **representative, physically-grounded terminal constraint**, not
  a validated floatplane load prediction — a faithful twin-float model (per-float W/N
  split, inter-float interference, one-float-first asymmetry) is new, unvalidated
  modeling and is deliberately out of scope.
- The bound is conservative because chine immersion is ignored. The true peak is at
  or below the prediction.
- Set `n_limit` from the structural allowable with knockdowns. The prediction is CG
  load factor normal to the water, not local hull pressure.
- Use AHAB's actual hull dead rise `beta` and mass. Absolute load level is dominated
  by dead rise and contact sink rate, so placeholder values can be off by large
  factors.
- For rough water, define `tau` and `gamma0` relative to the local wave slope and
  reference `z` to that surface.
- Guard the singularities at `gamma0 -> 0`, `tau -> 0`, and `phi_A -> 0`.

---

## 8. Vehicle & float parameters (DHC-2 Beaver on Wipline 6000 floats)

The plant is the flight-validated DHC-2 Beaver (see `paper_readiness.md` §6), which flies on
**twin Wipline 6000 floats**. The impact barrier therefore uses the **W/2-per-float** conservative
parameterization: `n_surfaces = 2`, effective weight per surface `W/2`. Because
`n_peak ∝ (W/n)^(−1/3)`, the two-float split raises the predicted CG load factor by
2^(1/3) ≈ 1.26× over a single surface at full weight — the conservative (safe) direction for
symmetric contact (§7 covers the single-surface scope this rides on). `n_surfaces = 1` (default)
recovers the single-hull (flying-boat / AHAB) behavior unchanged.

| Parameter | Knob | Value | Basis |
|---|---|---|---|
| Planing surfaces | `n_surfaces` | **2** | twin float |
| Aircraft mass | `mass` | 2313 kg (5100 lb std gross) | Beaver spec |
| Effective weight / float | W/2 | ~11.3 kN | mass·g / n_surfaces |
| Dead-rise (forebody) | `beta_deg` | **22.5°** | representative; float 15–40°, ≥20° practical, NACA Model 57 family 20–25° |
| Water density | `rho_water` | 1000 fresh / 1025 sea | standard |
| Keel incidence | `tau_keel_deg` | 0 (confirm ~0–3° float rigging) | representative |
| Float length | — | 7.49 m | Wipline 6000 (context; not in peak-load formula) |
| Float beam/width | — | 0.99 m | Wipline 6000 (context; chine-immersion / local-load check) |
| Buoyancy / float (fresh) | — | 5664 lb | Wipline 6000 (sanity: ~2× reserve vs gross) |

**Documented** (Wipaire spec/service manuals): float length, beam, buoyancy, gross weight.
**Representative** (float makers do not publish these): dead-rise and keel incidence — acceptable
because the barrier is a *representative* constraint (§7) and dead-rise is the only parameter that
moves the peak load. Conservatism note: `n_peak ∝ W^(−1/3)`, so a *lower* assumed weight is the
conservative direction; the standard 5100 lb gross is used (a heavier STC gross would predict a
slightly lower load). Sources: Wipaire DHC-2 Beaver float data; Gudmundsson *GA Aircraft Design*
App. C3; USNA float-design notes.

---

## 9. Equation reference

| Quantity | Expression | Paper eq |
|---|---|---|
| Approach parameter | `kappa = sin(tau)/sin(gamma0) * cos(tau + gamma0)` | 20 |
| Hull coefficient | `alpha = f_beta^2 * phi_A * rho * pi / (6 sin(tau) cos(tau)^2)` | 12a |
| Dead-rise function | `f_beta = pi/(2 beta) - 1` | 45 |
| End-flow correction | `phi_A = 1 - tan(tau)/(2 tan(beta))` | 49 |
| Load-factor coefficient | `Cl = -(yddot/ydot0^2) (W/(alpha g))^(1/3)` | 14 |
| Velocity ratio at max accel | root of eq 27 in `u = ydot/ydot0` | 27 |
| Draft coeff cubed at max accel | `y^3 (alpha g / W) = 2u/(7u + 6 kappa)` | 23a |
| Peak load-factor coefficient | `Cl_max = [2u(u+kappa)^2/(3u+2kappa)] [(7u+6kappa)/(2u)]^(1/3)` | 25 |
| Anchor at kappa = 0 | `u = 7/9`, `Cl_max = 0.6123` | 39, 40 |
| Peak load factor | `n_peak = Cl_max(kappa) ydot0^2 (alpha/(W g^2))^(1/3)` | from 15 |
| Height-relaxed barrier | `h = (n_limit - n_peak) + Phi(z)`, `Phi(0)=0`, `Phi' >= 0` | this spec |
