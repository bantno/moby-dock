# Safety-Critical Autonomous Water Landing via Control Barrier Functions: Mathematical Formulation and Derivations

This document outlines the mathematical foundation and Control Barrier Function (CBF) formulations for an autonomous seaplane landing system. The system employs a high-order CBF-QP (Quadratic Program) safety filter over a nominal flight controller to guarantee hull-safe touchdown sink rates, airspeed stall margins, and hydrodynamic impact-force limits.

## Nomenclature

All quantities are **SI** (m, s, rad, kg, N) unless noted; angles in radians. An overdot is a
time derivative ($\dot x = \mathrm{d}x/\mathrm{d}t$); a hat denotes a nondimensional rate
($\hat q$). The **Code** column gives the corresponding identifier in `lon_augmented.hpp` /
`impact_barrier.hpp` where one exists.

> **Overloaded symbols** (disambiguated by context): $\alpha$ — angle of attack *vs.* extended
> class-$\mathcal{K}$ functions $\alpha_i(\cdot)$ *vs.* the hull coefficient $\alpha_\text{hull}$;
> $W$ — weight $mg$ (§3.5) *vs.* the QP cost-weighting matrix (§4); $u$ — control vector (§1) *vs.*
> the impact velocity-ratio $\dot y/\dot y_0$ (§3.5); $A$ — QP constraint matrix (§4) *vs.* the
> aspect-ratio argument of $\phi(A)$ (§3.5); $\bar c$ — mean aerodynamic chord *vs.* class-$\mathcal{K}$
> gains $c_i$.

### State and control variables

| Symbol | Code | Units | Meaning |
|---|---|---|---|
| $h$ | `LH` | m | Altitude above the water surface |
| $V$ | `LV` | m/s | True airspeed (magnitude of the velocity vector) |
| $\gamma$ | `LGAM` | rad | Flight-path angle (velocity-vector inclination, $+$ up) |
| $\theta$ | `LTH` | rad | Pitch attitude (body-axis inclination) |
| $q$ | `LQ` | rad/s | Body pitch rate, $q=\dot\theta$ |
| $T$ | `LT` | N | Thrust (integral-augmented to a state) |
| $\dot T$ | `LTDOT` | N/s | Thrust rate (augmented state) |
| $\delta_e$ | `LDE` | rad | Elevator deflection (virtual, through the mixing map) |
| $\ddot T$ | `LTDDOT` | N/s² | Thrust second derivative (control input) |
| $\alpha$ | `alpha` | rad | Angle of attack, $\alpha=\theta-\gamma$ |
| $x,\ u$ | — | — | Unaugmented state $[h,V,\gamma,\theta,q]$ / control $[\delta_e,T]$ |
| $X,\ U$ | — | — | Augmented state $[h,V,\gamma,\theta,q,T,\dot T]$ / control $[\delta_e,\ddot T]$ |

### Aircraft, geometry, and aerodynamics

| Symbol | Code | Units | Meaning |
|---|---|---|---|
| $m$ | `mass` | kg | Aircraft mass |
| $I_{yy}$ | `Iyy` | kg·m² | Pitch moment of inertia |
| $S$ | `Sref` | m² | Wing reference area |
| $\bar c$ | `cref` | m | Mean aerodynamic chord (moment reference length) |
| $z_{cp}$ | `zcp` | m | Thrust-line vertical offset from the c.g. (thrust pitch moment) |
| $L,\ D$ | `Lift`,`Drag` | N | Wind-axis lift and drag |
| $C_L,\ C_D$ | `CL`,`CD` | — | Wind-axis lift/drag coefficients |
| $C_{Fx},C_{Fz}$ | `CFx`,`CFz` | — | Body-axis force coefficients (VSPAero deck) |
| $C_m$ | `CMy` | — | Body-axis pitching-moment coefficient |
| $C_{m0}$ | `off_CMy` | — | Zero-$\alpha$ pitching-moment coefficient |
| $C_{m\alpha}$ | `dAlpha_CMy` | 1/rad | Static pitch-stability derivative $\partial C_m/\partial\alpha$ |
| $C_{mq}$ | `dQ_CMy` | — | Pitch-damping derivative $\partial C_m/\partial\hat q$ |
| $C_{mM}$ | `dMach_CMy` | — | Mach derivative $\partial C_m/\partial\mathrm{Ma}$ |
| $C_{m\delta_e}$ | `dDe_CMy` | 1/rad | Elevator control-effectiveness derivative |
| $\hat q$ | `qhat` | — | Nondimensional pitch rate, $\hat q=q\bar c/(2V)$ |
| $\mathrm{Ma}$ | `mach` | — | Mach number, $V/a_\text{sound}$ |
| $a_\text{sound}$ | `a_sound` | m/s | Speed of sound |
| $\bar q$ | `qbar` | Pa | Dynamic pressure, $\tfrac12\rho_a V^2$ |
| $C_{L,\max}$ | `CL_max` | — | Maximum (stall) lift coefficient — config placeholder |
| $C_{D,\text{maxlift}}$ | — | — | Drag coefficient at the extrapolated max-lift $\alpha$ |

### Environment and physical constants

| Symbol | Code | Units | Meaning |
|---|---|---|---|
| $g$ | `g` | m/s² | Gravitational acceleration |
| $\rho_a$ | `rho` | kg/m³ | Air density |
| $\rho_w$ | `rho_water` | kg/m³ | Water density |

### Descent-rate and airspeed barriers (§3.1–3.2)

| Symbol | Code | Units | Meaning |
|---|---|---|---|
| $b$ | — | m/s | Descent-rate barrier value |
| $b_V$ | — | m/s | Airspeed (stall) barrier, $V-V_\text{min}$ |
| $b_{V,\max}$ | — | m/s | Airspeed (over-speed) barrier, $V_\text{max}-V$ |
| $v_\text{safe}$ | `v_safe` | m/s | Hull-safe touchdown sink rate |
| $a_\text{brk}$ | — | m/s² | Available braking (upward) acceleration at the current state |
| $V_\text{min}$ | `Vmin` | m/s | Stall-margin minimum airspeed |
| $V_\text{max}$ | `Vmax_air` | m/s | Never-exceed / over-speed airspeed (placeholder) |

### Contact-force barrier — conceptual (§3.3)

| Symbol | Units | Meaning |
|---|---|---|
| $b_\text{comp},\ b_F$ | m/s | Composite attitude-coupled contact-force barrier |
| $F_\text{max}$ | N | Allowable peak structural impact force |
| $C_s$ | — | Slamming (impact-pressure) coefficient |
| $\beta_\text{eff}$ | rad | Effective deadrise angle, $\beta_\text{hull}+(\theta-\theta_\text{surf})$ |
| $\beta_\text{hull}$ | rad | Hull deadrise angle |
| $\theta_\text{surf}$ | rad | Water-surface slope angle |
| $A_\text{ref}$ | m² | Reference wetted area |

### Impact-load barrier — implemented, NACA TN 1516 (§3.5)

| Symbol | Code | Units | Meaning |
|---|---|---|---|
| $b_\text{imp}$ | — | — | Impact-load barrier value |
| $n_\text{peak}$ | `n_peak` | g | Peak c.g. load factor a contact at the current state would produce |
| $n_\text{limit}$ | `n_limit` | g | Allowable load-factor limit |
| $\tau$ | `tau` | rad | Trim angle, $\tau=\theta-\theta_\text{keel}$ |
| $\theta_\text{keel}$ | — | rad | Keel/hull reference-line angle |
| $\gamma_0$ | — | rad | Approach flight-path angle, $\gamma_0=-\gamma$ |
| $\dot y_0$ | `ydot0` | m/s | Sink rate at contact, $\dot y_0=-V\sin\gamma$ |
| $\kappa$ | `kappa` | — | Approach parameter (eq 20), $\tfrac{\sin\tau}{\sin\gamma_0}\cos(\tau+\gamma_0)$ |
| $C_{lf}$ | `Clf` | — | Load-factor coefficient (eqs 25/27) — **not** the aerodynamic $C_L$ |
| $u$ | — | — | Velocity ratio $\dot y/\dot y_0$ at peak acceleration (eq 27 root) |
| $\alpha_\text{hull}$ | `alpha_hull` | — | Hull impact coefficient |
| $f(\beta)$ | `f_beta` | — | Deadrise function $\tfrac{\pi}{2\beta}-1$ (eq 45) |
| $\phi(A)$ | `phi_A` | — | Aspect-ratio/immersion correction $1-\tfrac{\tan\tau}{2\tan\beta}$ (eq 49) |
| $W$ | — | N | Weight, $W=mg$ |
| $K_0$ | `K0` | g·s²/m² | Frozen hull coefficient, $(\alpha_\text{hull}/(Wg^2))^{1/3}$ |
| $z$ | — | m | Height above water (impact-barrier argument) |
| $\Phi(z)$ | `Phi` | — | Height-relaxation budget, $N_b(1-e^{-z/z_s})$ |
| $N_b$ | `Nb` | — | Relaxation saturation budget (aloft value of $\Phi$) |
| $z_s$ | `zs` | m | Relaxation height scale |
| $z_\text{gate}$ | `z_gate` | m | Model-valid gating height (row assembled only below it) |
| $C_{lf,0},\ \tfrac{\mathrm{d}C_{lf}}{\mathrm{d}\kappa}$ | — | — | Frozen local-affine $C_{lf}(\kappa)$ model, anchored at $\kappa_0$ |

### CBF / QP machinery (§3–4)

| Symbol | Meaning |
|---|---|
| $\alpha_i(\cdot)$ | Extended class-$\mathcal{K}$ functions (linear here, $\alpha_i(s)=c_i s$) |
| $c_i$ | Class-$\mathcal{K}$ gains (per barrier; e.g. descent $[c_1,c_2,c_3]$) |
| $\psi_i$ | HOCBF cascade functions, $\psi_0=b,\ \psi_i=\dot\psi_{i-1}+\alpha_i(\psi_{i-1})$ |
| $L_F^k b$ | $k$-th Lie derivative of $b$ along the drift $F$ (drift stack) |
| $L_G L_F^{r-1}b$ | Control Lie derivative (the QP authority/row) |
| $r$ | Relative degree of a barrier |
| $A,\ \mathbf{b}_{qp}$ | QP constraint matrix and right-hand-side vector |
| $W$ | QP cost-weighting matrix (in $\tfrac12\|U-U_\text{nom}\|_W^2$) |
| $U_\text{nom}$ | Nominal control from the nominal controller |
| $U^\star$ | QP-optimal (safety-filtered) control |
| $K_p,\ K_d$ | Nominal thrust-tracking PD gains (§2) |

---

## 1. System Dynamics and State Space

The aircraft is modeled in the longitudinal (pitch-plane) utilizing a dynamic extension (integral augmentation) on the thrust channel to align the relative degrees of the control inputs.

### 1.1 Original System (Unaugmented)
* **State Vector:** $x = [h, V, \gamma, \theta, q]^T$
* **Control Vector:** $u = [\delta_e, T]^T$

**Equations of Motion:**
$$
\begin{aligned}
\dot h &= V\sin\gamma \\
\dot V &= \frac{1}{m}\big(T\cos\alpha - D(\alpha,V)\big) - g\sin\gamma \\
\dot\gamma &= \frac{1}{mV}\big(L(\alpha,V) + T\sin\alpha\big) - \frac{g\cos\gamma}{V} \\
\dot\theta &= q \\
\dot q &= \frac{1}{I_{yy}} \left( \frac{1}{2}\rho_a V^2 S \bar{c} \right) \left( C_{m0} + C_{m\alpha}\alpha + C_{mq}\frac{\bar{c} q}{2V} + C_{m\delta_e}\delta_e \right)
\end{aligned}
$$
*Assumption:* Angle of attack $\alpha = \theta - \gamma$.

### 1.2 Augmented System (Aligned Relative Degree)
To solve the mixed-relative-degree problem between elevator ($\delta_e$) and thrust ($T$), thrust is augmented via two integrators, becoming a state variable.

* **Augmented State Vector:** $X = [h, V, \gamma, \theta, q, T, \dot{T}]^T$
* **Augmented Control Vector:** $U = [\delta_e, \ddot{T}]^T$

---

## 2. Nominal Control Strategy (Constant-Thrust Powered Approach)

The nominal controller's only job is to put the aircraft on a reasonable glide slope; the
CBF-QP owns all safety-critical shaping. In particular **the flare is not designed into the
nominal** — it emerges from the degree-3 descent-rate barrier near the surface (§3.1). We
therefore use a deliberately minimal, platform-agnostic law: a **constant-thrust powered
approach** with a **cascade elevator loop** that holds a constant (negative) flight-path angle.
This replaces the PX4 TECS of earlier drafts. Code: `include/autoland/lon_nominal.hpp`; wiring
and trim seeding: `src/lon_sim.cpp`.

### 2.1 Thrust channel (PD on the augmented thrust state)

The system is augmented, so the QP optimizes $\ddot{T}$, but the operator wants to command a
thrust *level* $T_\text{set}$. To avoid differentiating a setpoint, we let the augmented thrust
state track $T_\text{set}$ through a PD law whose output is the control $\ddot T$:
$$
\ddot{T}_\text{nom} = K_{p,T}\,(T_\text{set} - T) - K_{d,T}\,\dot{T}.
$$
Thrust is held *constant* (no airspeed feedback in the nominal). The approach is **powered**
on purpose: carrying thrust keeps airspeed energy available for the emergent flare. Airspeed
itself is protected by the airspeed CBF (§3.2), not the nominal.

### 2.2 Elevator channel (flight-path-angle cascade)

A two-loop cascade converts a flight-path-angle reference into an elevator command.

**Outer loop — PI on $\gamma$ to an attitude command** (with an anti-windup clamp):
$$
\theta_\text{cmd} = \theta_\text{trim}
  + K_{p,\gamma}\,(\gamma_\text{ref} - \gamma)
  + K_{i,\gamma}\!\int (\gamma_\text{ref} - \gamma)\,\mathrm{d}t,
\qquad
\big|\theta_\text{cmd} - \theta_\text{trim}\big| \le \theta_{\text{cmd},\max}.
$$
The integrator is conditionally frozen (back-calculated) whenever the clamp is active, so it
does not wind up against the attitude limit.

**Inner loop — PD on $\theta$ to elevator** (with pitch-rate damping):
$$
\delta_{e,\text{nom}} = K_{p,\theta}\,(\theta_\text{cmd} - \theta) - K_q\,q.
$$

The nominal control vector passed to the QP is
$U_\text{nom} = [\delta_{e,\text{nom}},\ \ddot{T}_\text{nom}]^T$.

### 2.3 Feedforward / trim seeding

The feedforward terms are seeded from a 6-DOF trim solve at the approach condition
$(V_\text{app}, \gamma_\text{app})$ before the run, then optionally overridden by config:
$\theta_\text{trim}$ from the trim attitude, $\gamma_\text{ref} = \gamma_\text{app}$, and
$T_\text{set}$ from the trim thrust. Because $\theta_\text{trim}$ and $T_\text{set}$ only feed a
*feedback* controller, the small mismatch between the body-axis trim and the longitudinal model
is harmless. (The run itself starts from a steady level-flight trim of the longitudinal model,
`lonTrim()`; the nominal then commands $\gamma_\text{ref}$ and the aircraft pushes over into the
approach — see design doc §8.4.)

### 2.4 Runtime configuration (`data/lon_scenario.yaml`)

| Parameter | Symbol | Value | Note |
|---|---|---|---|
| Approach airspeed | $V_\text{app}$ | 18 m/s | sets the trim seed |
| Approach path angle | $\gamma_\text{app}=\gamma_\text{ref}$ | $-3^\circ$ | shallow (clean high-L/D aero) |
| Thrust setpoint | $T_\text{set}$ | 2.0 N | powered approach (config override of trim seed) |
| Thrust PD | $K_{p,T},\,K_{d,T}$ | 4.0, 4.0 | |
| $\gamma\!\to\!\theta_\text{cmd}$ PI | $K_{p,\gamma},\,K_{i,\gamma}$ | 2.0, 0.5 | |
| Attitude-command clamp | $\theta_{\text{cmd},\max}$ | $25^\circ$ | anti-windup limit |
| $\theta\!\to\!\delta_e$ PD | $K_{p,\theta},\,K_q$ | 1.0, 0.25 | scaled down ~8× for the real mass |

These are tuning values for the present airframe/scenario, not first-principles constants; see
the design doc and `TODO.md` for outstanding calibration items.

---

## 3. High-Order Control Barrier Functions (HOCBFs)

The safety filter enforces forward invariance on a set of barriers. The formulation uses extended class-$\mathcal{K}$ functions, chosen here as linear functions: $\alpha_i(x) = c_i x$.

### 3.1 Descent-Rate Barrier ($b$)
Ensures sufficient altitude to arrest the current sink rate before touchdown.
$$
b(X) = V\sin\gamma + \sqrt{v_\text{safe}^2 + 2\, a_\text{brk}(V,\gamma)\, h} \ge 0
$$

**Speed/path-angle-dependent braking acceleration.** The braking authority is the maximum
upward acceleration available at the current airspeed by pulling to the (effective) max-lift
condition — lift and drag minus gravity:
$$
a_\text{brk}(V,\gamma) = \frac{\rho_a V^2 S}{2m}\big[C_{L,\max}\cos\gamma - C_{D,\text{maxlift}}\sin\gamma\big] - g .
$$
**Thrust is deliberately omitted.** Putting the thrust state $T$ (or $\theta$) into $b$ would make
$b$ depend on $T$/$\theta$ directly and drop the barrier's relative degree below 3, breaking the
augmented HOCBF alignment (that $\theta$-coupling is exactly the §3.3 mixed-degree problem). Lift
and drag depend only on $(V,\gamma)$, which are already in $b$, so **relative degree 3 is
preserved** and the machinery below is unchanged. $C_{L,\max}$ is a config input (placeholder
pending calibration); $C_{D,\text{maxlift}}$ is evaluated at the extrapolated max-lift $\alpha$.
The implementation differentiates this $a_\text{brk}(V,\gamma)$ **exactly** via autodiff, so the
control-affine *closed forms* written below (derived for a *constant* $a_\text{brk}$) are now
reference only — the live Lie derivatives carry the extra $\partial b/\partial V$, $\partial b/\partial\gamma$ terms.

**Derivation & Authority:**
* $\dot{b}$ yields vertical acceleration ($\ddot{h}$), exposing $T$ but no controls.
* $\ddot{b}$ exposes $\dot{T}$ and $q$.
* $\dddot{b}$ exposes the controls $\ddot{T}$ and $\dot{q}$ (which contains $\delta_e$).

**Relative Degree:** $3$ for both $\delta_e$ and $\ddot{T}$.
**Control Affine Form:**
$$
\dddot{b} = L_F^3 b + \underbrace{\left[ \left( \frac{T}{m}\cos\theta + \frac{1}{m}\left( \frac{\partial L}{\partial \alpha}\cos\gamma - \frac{\partial D}{\partial \alpha}\sin\gamma \right) \right) \left( \frac{\rho_a V^2 S \bar{c}}{2 I_{yy}} \right) C_{m\delta_e} \right]}_{\text{Elevator Authority}} \delta_e + \underbrace{\left[ \frac{\sin\theta}{m} \right]}_{\text{Thrust Authority}} \ddot{T}
$$

### 3.2 Airspeed Barrier ($b_V$)
Prevents aerodynamic stall during the flare maneuver.
$$
b_V(X) = V - V_\text{min} \ge 0
$$

**Upper (over-speed) barrier.** A symmetric maximum-airspeed barrier
$$
b_{V,\max}(X) = V_\text{max} - V \ge 0
$$
bounds high-energy water impact / structural limits. It has the same relative
degree (3) and control-affine structure as $b_V$, with the authority signs
flipped: since $b_{V,\max}$ depends on $V$ only, $L_f^k b_{V,\max} = -L_f^k b_V$
for $k \ge 1$ and $L_g L_f^2 b_{V,\max} = -L_g L_f^2 b_V$ (only the $k=0$ value
differs, $V_\text{max}-V$ vs $V-V_\text{min}$). It therefore reuses the same
degree-3 QP machinery unchanged. Implemented as `AirspeedUpperBarrier` in
`hocbf.hpp` and wired into the filter via the `airspeed_upper` flag + `Vmax_air`
in `LonCBFConfig` (soft by default, mirroring the lower barrier). `V_max` is a
placeholder pending a real never-exceed / structural speed.

**Derivation & Authority:**
* $\dot{b}_V = \dot{V}$, exposing $T$.
* $\ddot{b}_V$ exposes $\dot{T}$ and $q$.
* $\dddot{b}_V$ exposes $\ddot{T}$ and $\dot{q}$ (containing $\delta_e$).

**Relative Degree:** $3$ for both $\delta_e$ and $\ddot{T}$.
**Control Affine Form:**
$$
\dddot{b}_V = L_F^3 b_V + \underbrace{\left[ -\frac{1}{m}\left( T\sin\alpha + \frac{\partial D}{\partial \alpha} \right) \left( \frac{\rho_a V^2 S \bar{c}}{2 I_{yy}} \right) C_{m\delta_e} \right]}_{\text{Elevator Authority}} \delta_e + \underbrace{\left[ \frac{\cos\alpha}{m} \right]}_{\text{Thrust Authority}} \ddot{T}
$$

### 3.3 Contact-Force Barrier ($b_F$ / $b_{comp}$)
Couples the allowable sink rate to pitch attitude based on von Kármán/Wagner water-entry slamming theory to bound the peak structural load.

Allowable velocity becomes a function of attitude $\theta$:
$$
v_\text{safe}^2(\theta) = \frac{2\,F_\text{max}}{\rho_w\,C_s(\beta_\text{eff}(\theta))\,A_\text{ref}}
$$
Where effective deadrise $\beta_\text{eff}(\theta) = \beta_\text{hull} + (\theta - \theta_\text{surf})$.

**Composite Barrier:**
$$
b_{comp}(X) = V\sin\gamma + \sqrt{v_\text{safe}^2(\theta) + 2 a_\text{brk} h} \ge 0
$$

**Derivation & Authority:**
Because $v_\text{safe}^2$ is a function of $\theta$, taking $\dot{b}_{comp}$ introduces $\dot{\theta} = q$.
Consequently, taking the *second* derivative ($\ddot{b}_{comp}$) introduces $\dot{q}$, exposing the elevator control $\delta_e$ one derivative earlier than thrust.

**Relative Degree:** $2$ for $\delta_e$, $3$ for $\ddot{T}$.
*Architecture Decision:* Enforced as a mixed-degree formulation at Degree 2. The elevator ($\delta_e$) acts alone to enforce this specific constraint, accurately reflecting the high-bandwidth nature of attitude adjustments prior to hydrodynamic impact.

**QP Constraint Formulation (Degree 2):**
$$
\ddot{b}_{comp} = L_F^2 b_{comp} + L_{G_{\delta_e}} L_F b_{comp} \cdot \delta_e
$$
$$
L_{G_{\delta_e}} L_F b_{comp} \cdot \delta_e \ge - \left( L_F^2 b_{comp} + c_1 \dot{b}_{comp} + c_2(\dot{b}_{comp} + c_1 b_{comp}) \right)
$$

### 3.4 Actuator Limit Barriers
Because thrust is an augmented state, its physical limits must be enforced via state-constrained CBFs.

**Minimum Thrust ($b_1 = T \ge 0$):**
* Relative Degree: 2 (with respect to $\ddot{T}$)
* Constraint: $\ddot{T} \ge -(c_{1,1} + c_{1,2})\dot{T} - c_{1,1}c_{1,2}T$

**Maximum Thrust ($b_2 = T_\text{max} - T \ge 0$):**
* Relative Degree: 2 (with respect to $\ddot{T}$)
* Constraint: $\ddot{T} \le -(c_{2,1} + c_{2,2})\dot{T} + c_{2,1}c_{2,2}(T_\text{max} - T)$

### 3.5 Impact-Load Barrier (Implemented — NACA TN 1516)
The rigorous realization of the §3.3 contact-load idea, using the closed(ish)-form peak
load factor of Milwitzky's V-bottom step-landing theory (NACA TN 1516) instead of a
hand-shaped $v_\text{safe}(\theta)$. Code: `include/autoland/impact_barrier.hpp`,
`scripts/precompute_impact_clf.py`. Reference: `documentation/impact_load_barrier_spec.md`.

**Peak load factor at contact.** With trim $\tau=\theta-\theta_\text{keel}$, flight-path
$\gamma_0=-\gamma$, sink $\dot y_0 = -V\sin\gamma$, and the approach parameter
$$
\kappa = \frac{\sin\tau}{\sin\gamma_0}\cos(\tau+\gamma_0) \quad\text{(eq 20)},
$$
the peak CG load factor (g) is
$$
n_\text{peak} = C_{lf}(\kappa)\,\dot y_0^{\,2}\,\Big(\tfrac{\alpha_\text{hull}}{W g^2}\Big)^{1/3},
\qquad
\alpha_\text{hull} = \frac{f(\beta)^2\,\phi(A)\,\rho_w\,\pi}{6\sin\tau\cos^2\tau},
$$
with $f(\beta)=\tfrac{\pi}{2\beta}-1$ (eq 45), $\phi(A)=1-\tfrac{\tan\tau}{2\tan\beta}$ (eq 49),
and $W=mg$. The **load-factor coefficient** $C_{lf}(\kappa)$ (NOT the aerodynamic lift
coefficient $C_{L,\max}$) has no closed form — the velocity ratio $u=\dot y/\dot y_0$ at
maximum acceleration is the root of the transcendental eq 27,
$$
\kappa\!\left[\tfrac{1}{1+\kappa}-\tfrac{1}{u+\kappa}\right]
= \ln\frac{(9u+6\kappa)(u+\kappa)}{(7u+6\kappa)(1+\kappa)},
\qquad
C_{lf}=\frac{2u(u+\kappa)^2}{3u+2\kappa}\Big(\tfrac{7u+6\kappa}{2u}\Big)^{1/3}\ \text{(eq 25)},
$$
so it is precomputed offline over $\kappa\in[0.2,10]$ and interpolated (anchor:
$C_{lf}(0)=0.6123$ at $u=7/9$).

**Barrier (height-relaxed, Option A).** Height above water $z$ enters only through a
nonnegative budget that vanishes at the surface:
$$
b_\text{imp}(X) = \big(n_\text{limit}-n_\text{peak}(\tau,\gamma_0,\dot y_0)\big) + \Phi(z),
\qquad \Phi(z)=N_b\big(1-e^{-z/z_s}\big).
$$
At $z=0$ this is the true load constraint $n_\text{limit}-n_\text{peak}\ge0$; aloft $\Phi\to N_b$
makes it trivially satisfied (touchdown-only). $K_0=(\alpha_\text{hull}/Wg^2)^{1/3}$ and a
local-affine $C_{lf}(\kappa)\approx C_{lf,0}+\tfrac{dC_{lf}}{d\kappa}(\kappa-\kappa_0)$ are
**frozen** at the evaluation point so the templated barrier is smooth (no $\sqrt[3]{\cdot}$ in
the Taylor jet) while the attitude coupling stays live through $\kappa$.

**Relative degree & the affine-row obstruction.** $b_\text{imp}$ is relative degree **2 via
the elevator** ($\theta\!\to\!q\!\to\!\delta_e$, since $n_\text{peak}$ depends on $\theta$
through $\kappa$) and **3 via thrust** ($V\!\to\!T\!\to\!\dot T\!\to\!\ddot T$). A single
control-affine row cannot mix the two: the exact third derivative
$$
\dddot b = L_F^3 b + (L_{G}L_F^2 b)u + L_F(L_{G}L_F b)u + (L_{G}L_{G}L_F b)u^2 + (L_{G}L_F b)\dot u
$$
carries non-affine $u^2$ and $\dot u$ terms once the elevator has entered at degree 2.
It is therefore enforced as a **clean degree-2 HOCBF** (the §3.3 architecture decision): the
thrust column $L_{G_{\ddot T}}L_F b_\text{imp}=0$ drops out and the elevator/flare acts alone,
$$
\underbrace{L_{G_{\delta_e}}L_F b_\text{imp}}_{a}\,\delta_e \;\ge\; -\big(L_F^2 b_\text{imp}
+ (c_1+c_2)L_F b_\text{imp} + c_1 c_2\, b_\text{imp}\big).
$$
Thrust still bounds impact load through the descent/sink-rate barrier (§3.1). The row is
assembled only in the model-valid window (below $z_\text{gate}$, descending, positive trim,
$\kappa$ clamped to $[0.2,10]$) and softened with a height-scheduled slack (Option C).

---

## 4. Final QP Architecture

The safety filter operates continuously, solving the following Quadratic Program at each timestep:

$$
U^\star = \arg\min_{U \in \mathbb{R}^2}\ \frac{1}{2}\|U - U_\text{nom}\|_W^2
$$

**Subject to the matrix inequality $A U \le \textbf{b}_{qp}$:**

1.  **Descent-Rate (Degree 3):** Shared authority between $\delta_e$ and $\ddot{T}$.
2.  **Airspeed (Degree 3):** Shared authority between $\delta_e$ and $\ddot{T}$.
3.  **Impact-Load (Degree 2, §3.5):** Enforced exclusively by $\delta_e$; height-relaxed
    via $\Phi(z)$ and softened with a height-scheduled slack. Realizes the §3.3 contact-load idea.
4.  **Min Thrust (Degree 2):** Enforced exclusively by $\ddot{T}$. (Hard constraint.)
5.  **Max Thrust (Degree 2):** Enforced exclusively by $\ddot{T}$. (Hard contraint.)

*Post-Processing:* The optimized $U^\star = [\delta_e^\star, \ddot{T}^\star]^T$ is returned. $\delta_e^\star$ is passed directly to the elevator servos. $\ddot{T}^\star$ is integrated twice onboard the flight computer to generate the final throttle setpoint passed to the PX4 control interface.