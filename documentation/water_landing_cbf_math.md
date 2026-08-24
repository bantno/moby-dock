# Safety-Critical Autonomous Water Landing via Control Barrier Functions: Mathematical Formulation and Derivations

This document outlines the mathematical foundation and Control Barrier Function (CBF) formulations for an autonomous seaplane landing system. The system employs a high-order CBF-QP (Quadratic Program) safety filter over a nominal flight controller to guarantee hull-safe touchdown impact loads, angle-of-attack (stall) margin, a nose-up touchdown attitude, and a bounded touchdown energy/speed.

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

## 2. Nominal Control Strategy

The nominal controller is responsible for the baseline glide-slope capture (currently a cascade holding a constant flight-path angle $\gamma$; slated to move to an exponential-altitude flare $h(t)=h_0 e^{-t/\tau}$). Because the system is augmented, the QP optimizes $\ddot{T}$, but the nominal outputs a raw thrust setpoint $T_\text{set}$. To avoid differentiating it, we track $T_\text{set}$ with the augmented thrust state via a PD law:

$$
\ddot{T}_\text{nom} = K_p (T_\text{set} - T) - K_d \dot{T}
$$

The nominal elevator command $\delta_{e,\text{nom}}$ comes from the pitch cascade. The nominal control vector passed to the QP is $U_\text{nom} = [\delta_{e,\text{nom}}, \ddot{T}_\text{nom}]^T$. **The barriers below depend only on the instantaneous state, so they are independent of this nominal / glide-path choice.**

---

## 3. High-Order Control Barrier Functions (HOCBFs)

The safety filter enforces forward invariance on a set of barriers, using extended class-$\mathcal{K}$ functions chosen as linear: $\alpha_i(s)=c_i s$. In the current ("recovery") set the **impact-load barrier (§3.5) is the only hard (no-slack) safety row**; the stall, nose-up, and energy barriers are soft (penalized slack) with a strict priority set by their slack weights,
$$
w_\text{stall}\,(10^5) \;>\; w_\text{energy}\,(10^4) \;>\; w_\text{noseup}\,(10^3),
$$
so under a conflict on the shared elevator the QP sacrifices nose-up first, then energy, and protects the stall guard. The min/max-thrust barriers (§3.4) are also hard, as **actuator-effectiveness / HOCBF-validity guards**. There is deliberately **no airspeed-floor barrier** — stall is handled directly in angle of attack by §3.1.

### 3.1 Stall (Angle-of-Attack) Barrier ($b_\text{stall}$)
Prevents aerodynamic stall and produces the pilot low-altitude stall recovery. With $\alpha=\theta-\gamma$,
$$
b_\text{stall}(X) = \alpha_\text{lim} - (\theta-\gamma) \ge 0,
\qquad \alpha_\text{lim} = \alpha_\text{stall} - \alpha_\text{margin},
$$
where $\alpha_\text{stall}$ is the NACA 4414 wing-stall AoA ($11^\circ$) and $\alpha_\text{margin}$ a buffer. As $\alpha\to\alpha_\text{lim}$ the HOCBF drives the elevator **nose-down** (saturating at the boundary) — i.e. *pitch down + full elevator*, exactly the pilot recovery.

**Derivation & Authority:**
* $\dot b_\text{stall} = -\dot\theta + \dot\gamma = -q + \dot\gamma$ — no control ($\dot\gamma$ carries $T$ but no input).
* $\ddot b_\text{stall}$ exposes $\dot q$, which contains $\delta_e$.

**Relative Degree:** $2$, via the elevator ($\theta\!\to\!q\!\to\!\delta_e$); $\ddot T$ is absent from the degree-2 control row. **Control-affine form** (linear in the state, exact via autodiff):
$$
\ddot b_\text{stall} = L_F^2 b_\text{stall} + \big(L_{G_{\delta_e}}L_F b_\text{stall}\big)\,\delta_e,
\qquad
L_{G_{\delta_e}}L_F b_\text{stall} = -\left(\frac{\rho_a V^2 S \bar c}{2 I_{yy}}\right)C_{m\delta_e}.
$$

### 3.2 Nose-Up (Attitude) Barrier ($b_\text{nose}$)
Holds the pitch attitude nose-up in the final metres so the impact model (§3.5) stays valid — its gate requires $\tau=\theta-\theta_\text{keel}>0$ — and so the touchdown is nose-up.
$$
b_\text{nose}(X) = \theta - \theta_\text{min} \ge 0,
\qquad \theta_\text{min} \ge \theta_\text{keel}.
$$
It is **$\theta$-based, not AoA-based**: a floor on $\theta$ directly guarantees the impact gate (an AoA floor would not, since a steep $\gamma$ lets $\theta=\alpha+\gamma$ go negative), and gives a clean degree-2 with no $\gamma$-coupling. The row is assembled only below $h_\text{noseup}$ (the final $\sim$3 m). As the softest row it **yields to the stall guard**: when a steep $\gamma$ makes the demanded nose-up conflict with the stall limit, §3.1 wins and the nose drops — the recovery.

**Relative Degree:** $2$, via the elevator ($\theta\!\to\!q\!\to\!\delta_e$). **Control-affine form:**
$$
\dot b_\text{nose} = q,\qquad
\ddot b_\text{nose} = L_F^2 b_\text{nose} + \big(L_{G_{\delta_e}}L_F b_\text{nose}\big)\,\delta_e,
\qquad
L_{G_{\delta_e}}L_F b_\text{nose} = \left(\frac{\rho_a V^2 S \bar c}{2 I_{yy}}\right)C_{m\delta_e}.
$$

### 3.3 Total-Energy Barrier ($b_E$)
Bounds the touchdown kinetic energy / horizontal speed via a height-scheduled **total-specific-energy ceiling**. With specific energy $E=\tfrac12 V^2 + g h$ and a ceiling $E_\text{cap}(h)=\tfrac12 V_\text{td,max}^2 + g_\text{eff}\, h$,
$$
b_E(X) = E_\text{cap}(h) - E = \tfrac12\big(V_\text{td,max}^2 - V^2\big) + (g_\text{eff}-g)\,h \ge 0
\;\Longleftrightarrow\; V \le \sqrt{V_\text{td,max}^2 + 2(g_\text{eff}-g)\,h}.
$$
A *constant* energy cap is useless for bounding touchdown speed — potential energy converts to kinetic on descent, so a constant-$E$ glider arrives fast — hence the height schedule, which makes it a **descending airspeed cap**. It is a **loose never-exceed ceiling**: $V_\text{td,max}$ is the hull/structural touchdown limit ($\approx 1.4\,V_\text{stall}$), and $g_\text{eff}$ is sized so the ceiling clears the whole approach and *starts satisfied*,
$$
g_\text{eff} \ge g + \frac{V_0^2 - V_\text{td,max}^2}{2 h_0}.
$$
$b_E$ is a **polynomial in $(V,h)$** ($C^\infty$, no $\sqrt{\cdot}$ — the airspeed cap is $b_E\ge0$ solved for $V$, never formed); only $V,h$ enter (not $\theta$/$T$), so the degree-3 alignment is preserved.

**Derivation & Authority:**
* $\dot b_E = -V\dot V + (g_\text{eff}-g)\dot h$ — exposes $T$ (through $\dot V$) but no control.
* $\ddot b_E$ exposes $\dot T$ and $q$.
* $\dddot b_E$ exposes both controls $\ddot T$ and $\dot q$ (containing $\delta_e$).

**Relative Degree:** $3$ for **both** $\delta_e$ (through drag via $\alpha$) and $\ddot T$ (through the $T,\dot T$ chain) — the natural energy lever is thrust, with the elevator a weaker drag lever. Enforced with the shared degree-3 machinery ($L_F^3 b_E$ plus a two-column control row), class-$\mathcal K$ of size 3.

### 3.3b Total-Energy Floor (opt-in — the powered stall recovery)
Dual of §3.3, **off by default** (the recovery set deliberately carries no always-on airspeed
floor; the AoA row is the stall guard). When enabled it adds the *throttle* half of the pilot
stall recovery — the channel the degree-2 stall row cannot command ($\ddot T$'s column of that
row is identically zero):
$$
b_{E,\text{fl}}(X) = E - E_\text{floor}(h) = \tfrac12\big(V^2 - V_\text{floor}^2\big) + (g-g_\text{floor})\,h \ge 0 ,
$$
degree 3 for both controls like the ceiling but with **opposite authority signs**: when it binds
the QP throttles **up** and noses **down**. $g_\text{floor}=g$ degenerates to a pure airspeed
floor $V\ge V_\text{floor}$ (gentle, stays near-tangent); $g_\text{floor}<g$ relaxes it with
altitude. Two hard-won usage rules (`data/lon_stall_recovery_cbf_efloor.yaml`):
* **Price thrust by its actuator range** ($w_{\ddot T}\sim 1/\ddot T_\text{max}^2$, elevator at
  unit price). Under the default identity control cost the elevator's row coefficients
  ($\sim10^3$–$10^4\times$ thrust's) make the QP satisfy the floor by *diving* while altitude
  lasts, and thrust moves only at elevator saturation — late and slammed.
* **Keep the floor gently binding.** A deeply-violated floor (e.g. a true-energy
  $g_\text{floor}=0$ form sized to bind mid-departure) produces raw $\psi$ units $\sim100\times$
  the other rows'; its quadratic slack then out-muscles the nominally higher-priority stall row
  (slack weights compare **raw units**, not intent) and the $\ddot T$ slamming drives the thrust
  chain into genuine thrust-row/box infeasibility. The soft-row priority ladder is only
  meaningful while every row stays near-tangent — a known sharp edge of this formulation.

### 3.4 Actuator Limit / HOCBF-Validity Barriers
Because thrust is an augmented state, its physical limits must be enforced via state-constrained CBFs. These are kept **hard** and serve double duty as the **actuator-effectiveness / HOCBF-validity guards**: every elevator-driven barrier's control coefficient scales with dynamic pressure ($L_{G_{\delta_e}}L_F b \propto \rho_a V^2$), so keeping the augmented thrust chain bounded (and the aircraft flying) preserves the relative-degree structure the whole formulation rests on. This is why no separate airspeed-floor CBF is needed.

**Minimum Thrust ($b_1 = T \ge 0$):**
* Relative Degree: 2 (with respect to $\ddot{T}$)
* Constraint: $\ddot{T} \ge -(c_{1,1} + c_{1,2})\dot{T} - c_{1,1}c_{1,2}T$

**Maximum Thrust ($b_2 = T_\text{max} - T \ge 0$):**
* Relative Degree: 2 (with respect to $\ddot{T}$)
* Constraint: $\ddot{T} \le -(c_{2,1} + c_{2,2})\dot{T} + c_{2,1}c_{2,2}(T_\text{max} - T)$

### 3.5 Impact-Load Barrier (the only hard safety row — NACA TN 1516)
The rigorous hydrodynamic touchdown-load barrier, using the closed(ish)-form peak
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
makes it trivially satisfied (touchdown-only). $K_0(\tau)=(\alpha_\text{hull}(\tau)/Wg^2)^{1/3}$ and
$C_{lf}(\kappa)$ are **frozen** at the evaluation point as **local-affine models** (value *and*
first derivative: $K_0\approx K_{0,0}+\tfrac{dK_0}{d\tau}(\tau-\tau_0)$ with
$\tfrac{dK_0}{d\tau}=\tfrac{K_0}{3}\big(-\cot\tau+2\tan\tau-\tfrac{\sec^2\tau}{2\tan\beta\,\phi(A)}\big)$,
and $C_{lf}\approx C_{lf,0}+\tfrac{dC_{lf}}{d\kappa}(\kappa-\kappa_0)$), so the templated barrier is
smooth (no $\sqrt[3]{\cdot}$ in the Taylor jet) while the full first-order attitude coupling stays
live. Freezing the $K_0$ *value* alone is not enough: at landing trims
$\tfrac{d\ln K_0}{d\tau}\approx-\tfrac{1}{3}\cot\tau$ partially cancels the retained $\kappa$
channel and its omission inflates $\partial n_\text{peak}/\partial\theta$ (hence the hard row's
elevator coefficient) by $\sim2\times$.

**Relative degree & the affine-row obstruction.** $b_\text{imp}$ is relative degree **2 via
the elevator** ($\theta\!\to\!q\!\to\!\delta_e$, since $n_\text{peak}$ depends on $\theta$
through $\kappa$) and **3 via thrust** ($V\!\to\!T\!\to\!\dot T\!\to\!\ddot T$). A single
control-affine row cannot mix the two: the exact third derivative
$$
\dddot b = L_F^3 b + (L_{G}L_F^2 b)u + L_F(L_{G}L_F b)u + (L_{G}L_{G}L_F b)u^2 + (L_{G}L_F b)\dot u
$$
carries non-affine $u^2$ and $\dot u$ terms once the elevator has entered at degree 2.
It is therefore enforced as a **clean degree-2 HOCBF** (elevator-only): the
thrust column $L_{G_{\ddot T}}L_F b_\text{imp}=0$ drops out and the elevator/flare acts alone,
$$
\underbrace{L_{G_{\delta_e}}L_F b_\text{imp}}_{a}\,\delta_e \;\ge\; -\big(L_F^2 b_\text{imp}
+ (c_1+c_2)L_F b_\text{imp} + c_1 c_2\, b_\text{imp}\big).
$$
Thrust still bounds the approach energy/speed through the total-energy barrier (§3.3); the impact
row itself is elevator-only. The row is assembled only in the model-valid window (below
$z_\text{gate}$, descending, positive trim, $\kappa$ clamped to $[0.2,10]$) and is now enforced
**hard** — the only hard safety row — with the height gate $\Phi(z)$ keeping it inactive except
near contact.

---

## 4. Final QP Architecture

The safety filter operates continuously, solving the following Quadratic Program at each timestep:

$$
U^\star = \arg\min_{U \in \mathbb{R}^2}\ \frac{1}{2}\|U - U_\text{nom}\|_W^2 + \sum_i \tfrac12 w_i\,\delta_i^2
$$

where $\delta_i\ge0$ is the slack on soft row $i$ (quadratically penalized by $w_i$; hard rows carry no slack). If the hard set is infeasible the filter falls back to a best-effort minimum-violation solve.

**Subject to the matrix inequality $A U \le \textbf{b}_{qp}$:**

1.  **Impact-Load (Degree 2, §3.5):** Enforced exclusively by $\delta_e$; height-relaxed via $\Phi(z)$. The only **hard** safety row.
2.  **Stall / AoA (Degree 2, §3.1):** Enforced by $\delta_e$. Soft, highest weight ($10^5$) — produces the pitch-down recovery.
3.  **Total-Energy (Degree 3, §3.3):** Shared authority between $\delta_e$ and $\ddot{T}$. Soft ($10^4$).
4.  **Nose-Up (Degree 2, §3.2):** Enforced by $\delta_e$, gated to the final metres. Soft, lowest weight ($10^3$) — yields to stall.
5.  **Min Thrust (Degree 2, §3.4):** Enforced exclusively by $\ddot{T}$. Hard (actuator/validity guard).
6.  **Max Thrust (Degree 2, §3.4):** Enforced exclusively by $\ddot{T}$. Hard (actuator/validity guard).

*Post-Processing:* The optimized $U^\star = [\delta_e^\star, \ddot{T}^\star]^T$ is returned. $\delta_e^\star$ is passed directly to the elevator servos. $\ddot{T}^\star$ is integrated twice onboard the flight computer to generate the final throttle setpoint passed to the PX4 control interface.

---

## 5. Beaver plant: single-integrator power and the degree-2 mixed-actuator impact barrier

This section derives the propulsion architecture used with the flight-validated DHC-2 Beaver plant (LR-556; see `paper_readiness.md` §6). It is **implemented and tested** in `include/autoland/beaver_lon.hpp` + `test/test_beaver_lon.cpp`, parallel to the AHAB path.

### 5.1 The change: power as a single-integrator state

The Beaver has **no separate thrust force** — net thrust and slipstream enter the aerodynamics through the normalized propeller total-pressure rise $d_{pt}$, which appears in the force coefficients $C_X, C_Z$ (and the pitching moment $C_m$ via $C_{m,d_{pt}}$). We therefore make **engine power $P$ a state** with control the **power rate** $u_P = \dot P$ — a *single* integrator — replacing the AHAB thrust chain $[T,\dot T]$ / $\ddot T$ (two integrators):

$$
X = [\,h,\ V,\ \gamma,\ \theta,\ q,\ P\,]^T,\qquad U = [\,\delta_e,\ u_P\,]^T,\qquad \dot P = u_P,\qquad d_{pt}=0.08696+\tfrac{2\cdot191.18}{\rho}\,\frac{P}{V^3}.
$$

The elevator is routed through the pitch **moment only** ($C_{m,\delta_e}$); the Beaver's direct elevator-lift term $C_{Z,\delta_e}=-0.398$ is dropped from the drift/$g$, keeping $\delta_e$ cleanly relative degree 2 (the same modeling choice the AHAB path makes).

### 5.2 Relative degree of the impact barrier to each actuator

The impact barrier $b = \big(n_\text{lim}-n_\text{peak}(\tau,\gamma_0,\dot y_0)\big)+\Phi(z)$ depends on $(\theta,\gamma,V,h)$ but **not** on $q$ or $P$. Differentiating along $\dot X = f(X)+g(X)U$:

$$
\dot b = \frac{\partial b}{\partial h}\dot h+\frac{\partial b}{\partial V}\dot V+\frac{\partial b}{\partial \gamma}\dot\gamma+\frac{\partial b}{\partial\theta}\,q .
$$

No control appears yet ($\delta_e$ enters only $\dot q$; $u_P$ enters only $\dot P$), so relative degree $\ge 2$. At the second derivative:

- **Elevator (moment channel):** $\dfrac{\partial b}{\partial\theta}\,\dot q$ with $\dot q = M/I_y$, $M\supset C_{m,\delta_e}\,\delta_e$ — so $\delta_e$ appears in $\ddot b$.
- **Power (force channel):** $\dfrac{\partial b}{\partial V}\ddot V$ (and $\tfrac{\partial b}{\partial\gamma}\ddot\gamma$) with $\ddot V \supset \dfrac{\partial \dot V}{\partial P}\,\dot P = \dfrac{\partial \dot V}{\partial P}\,u_P$ — so $u_P$ appears in $\ddot b$.

Both controls first appear at $\ddot b$: **uniform relative degree 2**. The control row therefore has **two** nonzero columns,

$$
L_gL_f\,b = \big[\,\underbrace{c_{\delta_e}}_{\text{moment}}\ \ \underbrace{c_{u_P}}_{\text{force}}\,\big],\qquad c_{\delta_e}\ne 0,\ c_{u_P}\ne 0,
$$

a **degree-2 mixed-actuator** row — vs. the AHAB double integrator, where thrust ($\ddot T$) enters the impact barrier only at degree $\ge 3$ (force channel) / $4$ (moment) and drops out of the degree-2 row. The power's degree-2 authority comes from the *force channel* (reducing sink rate); the $C_{m,d_{pt}}$ pitch coupling is a bonus that enriches the drift $L_f^2 b$ but reaches the control only at degree 3.

*(Verified in `test_beaver_lon.cpp` with the exact flow-Taylor Lie jet: at degree 1 both columns are $0$; at degree 2 both are nonzero.)*

### 5.3 Complementary authority vs. speed

Scaling the two columns:

$$
c_{\delta_e} \propto C_{m,\delta_e}\,\bar q\,S\bar c/I_y \ \propto\ \rho V^2 \quad(\text{weak at low speed}),\qquad
c_{u_P} \propto \frac{\partial \dot V}{\partial P}\ \propto\ \bar q\,\frac{\partial C_X}{\partial d_{pt}}\frac{\partial d_{pt}}{\partial P}\ \propto\ \rho V^2\cdot\frac{1}{\rho V^3}=\frac{1}{V}\quad(\text{strong at low speed}).
$$

The two actuators are **anti-correlated in speed**: exactly at the low-$\bar q$ flare/touchdown — where the elevator-only impact barrier is weakest (the "controllability guard" gap) — the power authority is strongest. *(Verified: the $|c_{u_P}|/|c_{\delta_e}|$ ratio is larger at 25 m/s than at 40 m/s.)*

### 5.4 Why this helps, and the cost

**Benefit.** The hard impact row is enforced by two independent control directions instead of one: it survives elevator-authority collapse near the water, improves QP feasibility (fewer best-effort recoveries), and — because both actuators are uniform degree 2 — realizes the *thrust-jointly-enforces-impact* construction (previously "highest research value") with the **standard** HOCBF, no non-affine $u^2/\dot u$ terms. The energy barrier likewise drops to degree 2. The system stays **control-affine** in $[\delta_e, u_P]$: $d_{pt}$'s nonlinearity in $(P,V)$ lives in the drift and is handled exactly by the frozen-affine + Taylor/autodiff Lie machinery.

**Cost.** Power becomes $C^0$ (rate-limited, not $C^1$) — physically reasonable for engine spool-up, but $u_P$ must be **bounded to the real spool rate** or the barrier leans on power it cannot deliver instantly. Because the hard row now depends on power authority, the $d_{pt}(P,V)$ map and spool bound must reflect the real engine, not be optimistic.
