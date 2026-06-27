> **ARCHIVED — historical reference.** This is the original brainstorming-era research
> pitch that started the project. The formalized math lives in
> [`../water_landing_cbf_math.md`](../water_landing_cbf_math.md) and the current design
> in [`../water_landing_cbf_design.md`](../water_landing_cbf_design.md). Kept for history;
> do not edit.

# Safety-Critical Autonomous Water Landing via Control Barrier Functions
*Working title — research pitch*

**In one line:** A CBF-based safety filter that guarantees a seaplane touches down
within hull-safe sink-rate and airspeed limits — and, optionally, within a
hydrodynamic impact-force limit — demonstrated in simulation and on a PX4 aircraft.

---

## Motivation

Autonomous fixed-wing landing research overwhelmingly targets runways. Water
landing is harder and far less studied: the surface is compliant and moving, and
touchdown becomes a safety-critical event governed by sink rate and attitude.
CBFs provide formal forward-invariance guarantees and act as a minimally-invasive
filter over *any* nominal controller — a natural fit that has not been applied to
the water-landing / slamming-load problem.

## Core idea

Work in the **longitudinal (pitch-plane) model** with elevator as the primary
control. A nominal glide-slope-capture + flare controller produces the commands; a
**CBF-QP safety filter** minimally adjusts them to keep the aircraft in a provably
safe set defined by two barriers:

- **Descent-rate CBF** — the soft-landing envelope: always retain enough authority
  to arrest to an acceptable touchdown sink rate. Under direct surface control this
  barrier is high relative degree, so it is enforced as an **HOCBF**.
- **Airspeed CBF** — stall avoidance. Couples with the descent barrier during flare
  (pitching up to slow the descent bleeds the airspeed you rely on).

Both are solved together in a single QP with actuator limits, so the filter passes
the nominal command through untouched until a constraint is about to be violated.

## Optional extension (the novel hook)

- **Contact-force CBF** — bound the *predicted peak slamming load* from water-entry
  theory (von Kármán / Wagner): the load scales as $\rho\, v_\text{rel}^2$ with a
  strong dependence on effective deadrise / trim. Treated as a pre-contact state
  constraint, it generalizes the descent-rate barrier so the **allowable sink rate
  becomes a function of touchdown attitude**. This ties the CBF safe set directly to
  impact physics and is the contribution most likely to make this a paper rather
  than an application note.

## Validation

- **Simulation** — demonstrate set-invariance and the guarantees against an
  unfiltered baseline across sea states; report touchdown sink rate, airspeed
  margin, and predicted impact load.
- **Hardware** — fixed-wing aircraft running **PX4**, with the QP filter inserted
  between guidance and PX4's control interface. (Whether the filter commands at the
  actuator level or through a rate/attitude setpoint sets the relative degree, and
  is itself a design point worth reporting.)

## Contributions

1. A CBF safety-filter formulation for autonomous water landing on longitudinal
   aircraft dynamics.
2. A contact-force barrier that grounds the safe set in hydrodynamic slamming
   theory, coupling sink rate and attitude.
3. HOCBF design demonstrated both in simulation and on a real PX4 platform.

## Scope / deferred to future work

Lateral and heading dynamics, the wave-direction (track-parallel-to-crests)
strategy, predictive wave-phase exploitation, and post-touchdown porpoising / the
on-water regime.

---

# Appendix: Model and Approach

## A.1 Longitudinal model

State $x = (h, V, \gamma, \theta, q)$ with angle of attack $\alpha = \theta - \gamma$;
control $u = \delta_e$ (elevator), with thrust $T$ as an optional second input.

$$
\begin{aligned}
\dot h &= V\sin\gamma \\
\dot V &= \tfrac{1}{m}\big(T\cos\alpha - D(\alpha,V)\big) - g\sin\gamma \\
\dot\gamma &= \tfrac{1}{mV}\big(L(\alpha,V) + T\sin\alpha\big) - \tfrac{g\cos\gamma}{V} \\
\dot\theta &= q \\
\dot q &= \tfrac{1}{I_{yy}}\, M(\alpha, q, \delta_e)
\end{aligned}
$$

with $L = \tfrac12\rho_a V^2 S\,C_L(\alpha)$, $D = \tfrac12\rho_a V^2 S\,C_D(\alpha)$,
$M = \tfrac12\rho_a V^2 S\bar c\,C_m(\alpha,q,\delta_e)$ and
$C_m = C_{m0} + C_{m\alpha}\alpha + C_{mq}\tfrac{\bar c q}{2V} + C_{m\delta_e}\delta_e$.
Air density $\rho_a$; water density $\rho_w$ enters in A.6.

## A.2 Nominal controller

The nominal layer only has to put the aircraft on a reasonable approach; the
CBF filter does the safety-critical shaping. So we keep it minimal and
platform-agnostic: **command a constant approach flight-path angle**
$\gamma_\text{ref}$ — a fixed glide slope (e.g. $-3^\circ$) — tracked by a simple
inner law (PID or LQR on the longitudinal states) that drives the elevator to hold
the corresponding pitch attitude, with throttle holding a nominal approach airspeed.
No flare is designed here.

The flare instead **emerges from the descent-rate barrier (A.3)**: as the
constant-glide-slope trajectory nears the envelope boundary close to the surface, the
CBF-QP bends the command upward to keep $b\ge0$, which is precisely a flare. In the
$(h,\dot h)$ phase portrait this is the straight constant-descent line being curved
into a gentle arrest near $h=0$. This is the cleaner division of labor — the nominal
is trivial, and the safety filter owns the part that matters.

The split also decouples simulation from hardware: in simulation the reference-angle
controller is a few lines of PID; on PX4 the same $\gamma_\text{ref}$ becomes an
attitude (or rate) setpoint through the existing interface. Either way the filter
(A.5) sees the same nominal $(\delta_e^\text{nom}, T^\text{nom})$ and is agnostic to
how they were produced.

## A.3 Descent-rate (soft-landing) barrier

The barrier is the constant-deceleration braking relation. Standing at height $h$,
descending at rate $v \equiv \dot h = V\sin\gamma < 0$, with a reliably-commandable
upward acceleration $a_\text{brk}>0$, apply maximum braking to touchdown. The
time-independent kinematic relation $s_f^2 = s_i^2 - 2a_\text{brk}\,\Delta x$
(i.e. $v^2 = u^2 + 2ax$ with $a=-a_\text{brk}$, $\Delta x = h$) gives the touchdown
descent speed from the current one $s_i = |v|$:

$$
s_f^2 = v^2 - 2 a_\text{brk} h .
$$

Requiring a hull-safe touchdown $s_f \le v_\text{safe}$ gives
$v^2 \le v_\text{safe}^2 + 2 a_\text{brk} h$, and since $v<0$,

$$
b(h,v) = v + \sqrt{v_\text{safe}^2 + 2\,a_\text{brk}\,h} \;\ge\; 0,
\qquad \mathcal C = \{b\ge0\}.
$$

**Reading the terms.**

- $v$ — the current sink rate, the quantity being lower-bounded.
- $v_\text{safe}^2$ — the touchdown-speed budget. It is what the radicand collapses to
  at $h=0$, making $b\ge0 \Leftrightarrow |v|\le v_\text{safe}$ there. Setting
  $v_\text{safe}=0$ demands a perfectly soft kiss and leaves the bare braking
  parabola $v+\sqrt{2a_\text{brk}h}$.
- $2a_\text{brk}h$ — the **velocity-squared margin the remaining altitude buys**. The
  factor 2 and the linearity in $h$ are clearest in energy form. Halving the squared
  inequality:

$$
\tfrac12 v^2 - a_\text{brk} h \;\le\; \tfrac12 v_\text{safe}^2 .
$$

$\tfrac12 v^2$ is the specific kinetic energy of the descent; $a_\text{brk}h$ is the
specific **work** the braking acceleration can do over the remaining drop (force ×
distance, hence linear in $h$). So: descent KE minus the energy braking can still
bleed off must not exceed the allowable touchdown KE. Un-halved, $2a_\text{brk}h$ is
twice that brakeable work — the $\Delta(v^2)=2a\,\Delta x$ budget the height
provides. High up it is large and a fast descent is permitted; as $h\to0$ it
vanishes and the envelope clamps to $v_\text{safe}$.

**Two notes.** (i) This curve is the boundary of the *maximal control-invariant set*
for the idealized double integrator with bounded deceleration — the largest set from
which a soft landing is recoverable. The real system cannot produce $a_\text{brk}$
instantly — building lift takes a pitch/attitude change — so its true recoverable set
is smaller; hence $a_\text{brk}$ is chosen conservatively and/or backed by the backup
set in A.5. (ii) Keeping $v_\text{safe}>0$ holds the radicand away from zero, so
$\partial b/\partial h \sim a_\text{brk}/\sqrt{\,\cdot\,}$ stays bounded at touchdown;
with $v_\text{safe}=0$ the gradient blows up at $h=0$ and the QP misbehaves exactly
when it matters. Since real braking authority degrades with speed and density, the
honest form is $a_\text{brk}=a_\text{brk}(V)$.

## A.4 Airspeed (stall) barrier

$$
b_V(V) = V - V_\text{min}, \qquad V_\text{min} = V_\text{stall} + \Delta V.
$$

Near flare this competes with the descent barrier (A.3) — pitching up to reduce sink
rate sheds the airspeed it depends on — so the descent coefficient is made
speed-dependent, $a_\text{brk} = a_\text{brk}(V)$, and the QP arbitrates the two.
Throttle, if used as a second input, is the natural airspeed authority and eases the
conflict directly.

## A.5 Safety-filter QP

Each barrier $b_i \in \{b,\, b_V\}$ (and optionally $b_F$) defines a safe set to stay
within. Given the nominal command $u^\text{nom}$ from A.2, the filter solves

$$
u^\star = \arg\min_u\ \tfrac12\|u - u^\text{nom}\|_W^2
\quad\text{s.t.}\quad
\dot b_i(x,u) + \alpha_i\big(b_i(x)\big) \ge 0\ \ \forall i,
\quad u \in \mathcal U,
$$

with $u = \delta_e$ (and throttle as an optional second input), each $\alpha_i$ a
class-$\mathcal K$ function, and $\mathcal U$ the actuator limits. The filter passes
the nominal command through whenever the state sits comfortably inside every safe set
and bends it only near a boundary; the condition above is the basic CBF form, with any
high-relative-degree handling left to a later stage.

Control invariance under the actuator limits is not automatic: size $a_\text{brk}$
conservatively, or define the safe set via a committed-flare backup controller, so the
descent set $\mathcal C$ stays invariant — a property to prove rather than tune.

## A.6 Contact-force barrier (extension)

Let $h_\text{rel} = h - \eta$ and $v_\text{rel} = \dot h - \dot\eta$ be clearance and
closing rate to the (possibly moving) surface $\eta$; the distance sensor reads
$h_\text{rel}$, and in still water $\eta = 0$, $v_\text{rel} = \dot h$. The normal impact
velocity at contact is $v_n = -v_\text{rel}\big|_{h_\text{rel}=0}$.

Water-entry slamming theory (von Kármán momentum / Wagner wetted-surface expansion)
gives a peak load

$$
F_\text{pred} = \tfrac12\,\rho_w\,C_s(\beta_\text{eff})\,A_\text{ref}\,v_n^2,
$$

where the slamming coefficient $C_s$ grows sharply as the effective deadrise
$\beta_\text{eff} \to 0$ (flat impacts are worst, V-bottoms cushion), and

$$
\beta_\text{eff} = \beta_\text{hull} + (\theta - \theta_\text{surf})
$$

is the hull deadrise adjusted by trim relative to the local surface slope
$\theta_\text{surf}$. The barrier

$$
b_F = F_\text{max} - F_\text{pred}(v_n, \beta_\text{eff}) \ge 0
$$

is enforced as a constraint on the **pre-contact state** (the impact itself is impulsive
and uncontrollable, so we bound the predicted load that *would* result if contact
occurred). Equivalently, it tightens A.3 into an attitude-dependent sink-rate limit,

$$
v_\text{safe}(\beta_\text{eff}) = \sqrt{\frac{2\,F_\text{max}}{\rho_w\,C_s(\beta_\text{eff})\,A_\text{ref}}},
$$

substituted into $b$. So the contact-force barrier is a strict generalization of the
descent-rate barrier in which the allowable touchdown sink rate depends on attitude —
and, through $\theta_\text{surf}$, on where in the wave the contact occurs.
