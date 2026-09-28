# Flipping the Reach-Energy CBF into a Ceiling: Design Notes for Gentle Touchdown

*Design discussion (no implementation yet), 2026-09-09. Context: the reach-energy FLOOR
h_E >= 0 is implemented and evaluated (see `energy_reach_cbf.md`). The eventual goal is
to flip the signs into an energy CEILING and pair floor + ceiling + the impact CBF to
enforce a gentle touchdown within a desired region — the way a seaplane pilot manages an
approach. These notes capture the barriers in the current implementation and how to
address them, before any derivation is done.*

## 1. The trap in the naive sign flip

The floor prices the remaining flight at the **current** alpha's L/D — that is the whole
point of the reformulation, and it is what gives the alpha-dot channel its RD-1 gradient.
Flip the sign with the same pricing and two things go wrong:

1. **Same conservation, same myopia.** h_ceil = (½V_td² + g·d·Ψ(α)) − E is conserved
   along any on-polar glide, exactly like the floor (`energy_reach_cbf.md` §4.1). A
   conserved ceiling never binds until a transient pushes it — the opposite of the pilot
   behavior wanted (start bleeding energy *early*, arrive on speed). A ceiling should
   *decay* when you fly clean with excess energy, so it progressively forces dissipation.
   That happens naturally if the ceiling is priced at **maximum dissipation capability**
   (worst achievable L/D in the current configuration) instead of at the current alpha:
   flying clean then makes the cap shrink faster than you are spending, the row binds
   gradually, and the funnel is automatic.
2. **Current-alpha pricing degenerates at low alpha for a ceiling.** Ψ(α̃) → ∞ as
   α̃ → 0, which for the floor was safely conservative ("you can't glide at zero lift").
   Flipped, it reads "diving at tiny alpha gives unbounded dissipation per ground meter"
   — the QP's cheapest lever for shedding energy becomes pushing alpha DOWN, i.e. a
   screaming dive at the water. Legal per the math (the steady-glide assumption behind
   g·d·(D/L) is being abused outside its validity), catastrophic per the airplane.

**Proposal:** price the ceiling's Ψ at the worst L/D achievable *within the trusted band
given the current drag configuration* — a function of sideslip/bank and flap, not of
instantaneous alpha. The ceiling loses its alpha gradient, but gains a clean actuator
allocation with no floor/ceiling tug-of-war on alpha: the floor owns pitch (fly
efficiently enough to reach), the ceiling owns the drag devices and the path (slip, bank,
S-turn, throttle-down) — exactly how a pilot splits the problem.

## 2. The three named enablers

### Bank-to-turn (biggest bang, least work)

χ̇ = g·tanμ/(V·cosγ) gives ~9°/s at 30° bank vs the skid's ~2°/s (turn radius ~220 m
instead of ~1.2 km), making both target capture *and* path-lengthening S-turns real
maneuvers instead of gestures. Two bonuses:

* Banking adds **modeled dissipation** — induced drag scales with n² = 1/cos²μ (a 45°
  bank doubles it) — so μ is a ceiling actuator too.
* Realization is *easier* than the rudder sideslip loop: the 6-DOF roll/cross-track
  machinery already exists, and the special wings-level lateral mode goes away.

One modeling obligation: bank tilts lift, so γ̇ gets a cosμ factor and the wings-level
point-mass assumption must go.

### The impact CBF (necessary — energy alone cannot produce "gentle")

At z = 0, E = ½V²; sink is V·sinγ, and γ is completely unconstrained by any energy
quantity — you can arrive inside a perfect energy corridor at 6 m/s sink. Gentle
touchdown is the flare, and the flare is the TN 1516 impact row's job (it already
produces exactly this behavior in the lon sim).

The barrier to using it here is **composition**: the impact row is RD-2 in the elevator
on the augmented lon model, while the reach layer is RD-1 in (α̇, β) on the point-mass —
different models, different actuation spaces — and the current layer overrides theta_cmd
wholesale, so the two filters would fight over the same pitch axis. Clean resolution:
derive everything on the 6-DOF model so all rows share one QP over (de, da, dr, dT),
following the lon filter's existing stacking pattern (impact hard, others soft, gates
like z_gate). Enabler already in place: `BeaverDynamics::xdotT` is templated for
autodiff, so exact HOCBF Lie derivatives on the *actual plant* are available. Wrinkle: a
few aero terms are not control-affine (the dpt² and de·β² terms) — freeze them locally
the way `AeroLocal` already does.

### The energy ceiling at touchdown

The repo already has the skeleton: `EnergyBarrier` in `hocbf.hpp` is literally
E <= ½V_td² + g_eff·h. What the flipped reach barrier adds is the d/(L/D) structure and
the heading dependence — and the heading channel's meaning flips delightfully: for the
ceiling, *increasing* d raises the allowed energy, so the barrier discovers S-turns and
360s on its own. That is the pilot's "extend the pattern to lose altitude," and it falls
out of the arc-distance term for free.

## 3. Barriers not on the original list

* **The target-blind nominal is the biggest gap.** *(Addressed 2026-09-26: the pattern
  guidance of `corridor_landing_roadmap.md` Phase 0 — Dubins path to a final-approach fix,
  L1 course tracking, glideslope cone — is exactly this "pilot".)* In the current runs the CBF was doing
  guidance, which is why even the wins missed by ~250 m. Floor + ceiling define a
  corridor; corridors *guard* a trajectory, they do not generate one. A nominal that
  tracks the aim point with a speed schedule (the "pilot") turns both barriers into rare
  interventions and dissolves most of the feasibility pressure. Do this before either
  barrier gets more sophisticated.
* **Dissipation levers don't exist in the model yet.** The ceiling's actuators —
  sideslip drag (the plant charges ~25% extra drag at β = 10° that the point-mass
  ignores), bank-induced drag, flaps, throttle-to-idle — are all unmodeled or absent.
  Related sign-flip hazard: unmodeled thrust was *conservative* for the floor but is
  *anti-conservative* for the ceiling (energy arriving that the barrier can't see), so
  the ceiling needs thrust in the drift or an idle-when-binding rule.
* **Corridor feasibility near the spot.** Floor and ceiling converge to the band
  [½V_min², ½V_td²] as d → 0; with bounded inputs and the measured myopia, the squeeze
  can go infeasible. The ck lesson carries over (both rows need small class-K gains so
  the funnel forms early), and capability-pricing the ceiling helps because it removes
  the alpha-axis conflict entirely.
* **Point target vs region.** "Within a desired region" falls out naturally if the floor
  references the *near* edge (must be able to reach it) and the ceiling the *far* edge
  (must be able to get down by it) — same d function, different target parameters per
  row. That also sidesteps the degenerate bearing-at-the-point geometry. Lateral
  containment comes back for free with bank-to-turn (the cross-track loop returns).
* **Envelope honesty.** The gentle-touchdown regime (slow, high-alpha, in ground effect)
  sits below the LR-556 30 m/s validity floor, on a plant with no stall aero; and the
  ceiling's dive escape hatch argues for placard rows (dynamic-pressure/V ceiling, keep
  the AoA band). State as a limitation now, consider a stall blend later (cf.
  `paper_readiness.md`).

## 4. Proposed sequencing

1. **Bank-to-turn** in the point-mass layer (μ replaces β; existing roll loops realize
   it).
2. **Target-tracking nominal** so the barriers stop doing guidance.
3. **Ceiling row** with capability pricing + modeled dissipation levers, region
   semantics (near/far edges), thrust handled.
4. **Impact row on the 6-DOF plant** (templated dynamics makes this mostly
   Lie-machinery work), gated near the surface; energy corridor hands off to it.
5. Only then the **full 6-DOF derivation** of the pair — steps 1–4 will have settled the
   design questions that derivation needs answered.

The pleasing property is the 1:1 map to airmanship: floor = "don't get low and slow away
from the spot," ceiling = "don't arrive high and hot" (slip / S-turn / flaps / idle),
impact row = the flare. If the math lands anywhere else, it is probably wrong.
