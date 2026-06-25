# Safety-Critical Autonomous Water Landing via Control Barrier Functions: Mathematical Formulation and Derivations

This document outlines the mathematical foundation and Control Barrier Function (CBF) formulations for an autonomous seaplane landing system. The system employs a high-order CBF-QP (Quadratic Program) safety filter over a nominal flight controller to guarantee hull-safe touchdown sink rates, airspeed stall margins, and hydrodynamic impact-force limits.

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

## 2. Nominal Control Strategy (PX4 TECS Integration)

The nominal controller is responsible for the baseline glide-slope capture. We utilize the PX4 Total Energy Control System (TECS). 

Because the system is augmented, the QP optimizes $\ddot{T}$, but TECS outputs a raw thrust/throttle setpoint ($T_\text{TECS}$). To avoid the noise amplification of differentiating the TECS output, we treat $T_\text{TECS}$ as a tracking reference for our augmented thrust state using a Proportional-Derivative (PD) control law:

$$
\ddot{T}_\text{nom} = K_p (T_\text{TECS} - T) - K_d \dot{T}
$$

The nominal elevator command $\delta_{e,\text{nom}}$ is taken directly from the TECS pitch control output. The nominal control vector passed to the QP is $U_\text{nom} = [\delta_{e,\text{nom}}, \ddot{T}_\text{nom}]^T$.

---

## 3. High-Order Control Barrier Functions (HOCBFs)

The safety filter enforces forward invariance on a set of barriers. The formulation uses extended class-$\mathcal{K}$ functions, chosen here as linear functions: $\alpha_i(x) = c_i x$.

### 3.1 Descent-Rate Barrier ($b$)
Ensures sufficient altitude to arrest the current sink rate before touchdown.
$$
b(X) = V\sin\gamma + \sqrt{v_\text{safe}^2 + 2 a_\text{brk} h} \ge 0
$$

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

---

## 4. Final QP Architecture

The safety filter operates continuously, solving the following Quadratic Program at each timestep:

$$
U^\star = \arg\min_{U \in \mathbb{R}^2}\ \frac{1}{2}\|U - U_\text{nom}\|_W^2
$$

**Subject to the matrix inequality $A U \le \textbf{b}_{qp}$:**

1.  **Descent-Rate (Degree 3):** Shared authority between $\delta_e$ and $\ddot{T}$.
2.  **Airspeed (Degree 3):** Shared authority between $\delta_e$ and $\ddot{T}$.
3.  **Contact-Force (Degree 2):** Enforced exclusively by $\delta_e$.
4.  **Min Thrust (Degree 2):** Enforced exclusively by $\ddot{T}$. (Hard constraint.)
5.  **Max Thrust (Degree 2):** Enforced exclusively by $\ddot{T}$. (Hard contraint.)

*Post-Processing:* The optimized $U^\star = [\delta_e^\star, \ddot{T}^\star]^T$ is returned. $\delta_e^\star$ is passed directly to the elevator servos. $\ddot{T}^\star$ is integrated twice onboard the flight computer to generate the final throttle setpoint passed to the PX4 control interface.