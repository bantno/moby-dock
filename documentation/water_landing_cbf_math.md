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
