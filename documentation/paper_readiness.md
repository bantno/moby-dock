# Paper-Readiness Assessment

**Status:** living assessment — *what we have, what a reviewer will attack, and what remains to
reach a defensible conference submission.* Update it as gaps close. Written 2026-08-19.

| Doc | Job |
|---|---|
| **`paper_readiness.md`** (this) | Publishability: assets, gaps, candidate framings, roadmap. |
| `water_landing_cbf_design.md` | Current design: theory, approach, assumptions, decisions. |
| `water_landing_cbf_math.md` | Formal derivations. |
| `CHANGELOG.md` | Chronological "what changed, when, who, why." |
| `../TODO.md` | Granular backlog. |

---

## 1. What this is

A fast closed-loop **design sandbox** (not onboard code) for the autoland of the AHAB V-tail
flying-boat seaplane UAV. The research artifact is a **CBF-QP safety filter for water landing**,
built on a single nonlinear 6-DOF body-axis EOM (`Dynamics::xdot`) derived from a real VSPAERO
deck and real vehicle mass properties. Trim, linearization, and both sims all derive from that
one EOM.

## 2. What is solid (the paper assets)

| Asset | State | Paper value |
|---|---|---|
| Longitudinal CBF-QP filter (OSQP) | Implemented, tested | Core |
| Exact Lie derivatives — Taylor-jet + autodiff, **no finite differences** | Implemented, oracle-checked | Genuine method point |
| Barrier set: impact-load (hard), stall/AoA, nose-up, energy ceiling, thrust guards | Implemented | Core |
| **Emergent stall-recovery suite** — 6 scenarios + noisy seeds + recovery-floor heatmap, all A/B vs CBF-off, all PASS | Complete | **Strongest, cleanest result** |
| **Hydrodynamic impact-load HOCBF** (NACA TN 1516 slam load) | Implemented | **Most novel** |
| Plant realism: MIL-F-8785C gust, JONSWAP/Airy waves (STANAG 4194 / USACE CEM), NACA 4414 viscous-stall overlay | Implemented | Supporting |
| 6-DOF straight-in sim (cascaded PID, crosswind crab, waves) | Implemented — **no CBF, no flare/decrab** | Future-work base |
| Flight-log (PX4 ULog) validation of aero/thrust placeholders | Done | Credibility |

64/64 tests pass; real vehicle data.

## 3. Gaps a reviewer will attack

Stated bluntly (most are also in `TODO.md` and the design doc §8):

1. **Self-consistent plant = zero model mismatch.** The sim plant *is* the CBF model. Clean
   barrier-invariance partly reflects structure, not robustness — this validates the *math*, not
   robustness to model error. **This is the #1 objection.** (See §6 for the chosen fix.)
2. **The hard safety barrier never binds.** On nominal approaches `n_peak ≈ 0.17 g` vs the 3 g
   limit; even forced hard at `n_limit = 0.05 g` the trajectory is unchanged. Per TODO: *"it works
   is indistinguishable from it's off."* The headline safety barrier has never been shown to be the
   active constraint doing its job.
3. **Impact physics are placeholders** — `n_limit`, dead-rise `beta`, `Nb`, `zs`, `c_impact` are
   uncalibrated. The "18.8 g wave-referenced" figure is explicitly *not a calibrated load*.
4. **CBF is longitudinal-only.** The 6-DOF sim has no filter, no flare, no decrab. Crabbed /
   wing-down touchdowns — the worst *local* load case — are entirely outside the barrier set.
5. **Filter is wave-blind.** The measured flat-vs-wave gap (~6× over `n_limit`) is a strong
   *motivating* result but shows the filter as-is is unsafe in waves; the wave-aware barrier
   (Stage 2) is deliberately not attempted.
6. **Idealized estimation / actuators** — perfect full-state feedback, no servo lag; the impact
   row is the most noise-sensitive (elevator coeff ~3e4) and has not been exercised with sensor
   noise + lag.

## 4. Candidate paper framings

The choice drives the remaining work.

- **A. "Emergent stall recovery via a CBF safety filter for seaplane autoland"** — *closest to
  ready.* The stall suite is a complete A/B study with noise seeds, gust/shear, and a
  recovery-floor sweep. Novelty: an AoA-direct barrier producing pilot-like pitch-down recovery as
  *emergent* behavior with no explicit recovery logic. Remaining work is modest. Lower novelty
  ceiling.
- **B. "A hydrodynamic slam-load control barrier function for autonomous water landing"** — *most
  novel, most work.* Nobody has an HOCBF from TN 1516 slam physics. But it is the least validated
  (gaps #2, #3, #5).
- **C. Combined water-landing CBF safety filter** — broadest, needs the most gaps closed.

## 5. Roadmap

**Must-do for any venue** (closes objection #1, the fatal one):

- A **model-mismatch / robustness** demonstration — the plant is perturbed off the CBF's model and
  the guarantee is shown to degrade gracefully. This single result separates "validates the math"
  from "publishable control result." **See §6 — the chosen route is a validated third-party
  plant.**
- A crisp **contribution statement + baseline** (CBF-off is the A/B; the stall suite already has
  this).

**Framing A (recommended first submission):** run the robustness demonstration on the stall
scenarios; write up the recovery-floor heatmap as the headline figure (already generated).

**Framing B/C (impact-load CBF):** construct a scenario where the impact row **actually binds** and
holds `n_peak ≤ n_limit`; **calibrate or defensibly bound** the impact physics (or reframe them as
normalized/parametric); exercise the impact row with sensor noise + servo lag; Stage-2 wave-aware
barrier is a large lift, likely future work.

**Nice-to-have (any framing):** exponential-altitude nominal flare (nominal-only change, barriers
already glide-path agnostic); PX4 insertion-point discussion (the relative-degree design point is
of genuine interest to a GNC audience).

## 6. Strategic direction — validated third-party plant

**Decision (2026-08-19, final):** replace the self-authored VSPAERO-based plant with the
**flight-validated 6-DOF model of the de Havilland DHC-2 Beaver** (see model selection below).
Rationale:

- It **eliminates the aero model as a confound.** What the paper reports becomes *the CBF
  machinery* evaluated against an independent plant, not a controller matched to its own plant —
  the direct, decisive fix for gap #1.
- It creates **genuine, characterizable model mismatch**: the CBF keeps a reduced control-affine
  model (for the exact Lie derivatives) while the plant is the full validated nonlinear model. The
  robustness result the paper needs falls out of the architecture instead of being staged.
- It reuses the plant-agnostic machinery: `lie_taylor.hpp`, `lon_cbf_filter`, `impact_barrier`,
  `hocbf`, and the QP are independent of *which* aircraft supplies the coefficients.

**What changes:**

- **Plant:** the validated model (co-simulated, or its aero exported to a table), replacing
  `Dynamics::xdot` as ground truth.
- **CBF internal model:** a control-affine coefficient-table model of the *same* aircraft (so the
  mismatch is "reduced-order / quasi-steady vs full plant" — honest and quantifiable), or a
  deliberately simpler model to stress the filter.
- **Safety metric:** the existing barrier set (impact-load TN 1516, stall/AoA, energy) evaluated on
  the new plant's state.

**Resolved decisions (2026-08-19):**

- **Aero-only plant.** No float hydrodynamics in the plant. The TN 1516 impact model remains the
  touchdown safety evaluation; the plant only needs validated aerodynamics good enough to fly a
  credible approach to the contact state (V, γ, θ, sink) the impact model consumes.
- **JSBSim set aside.** JSBSim is a simulation *engine* — we already have a tested nonlinear 6-DOF
  engine wired to the CBF/QP/wave/impact machinery, so its engine + external-control API are
  redundant, and its GA aircraft aero is community-tuned, not rigorously flight-validated. The
  independence we want comes from the **data provenance**, not the engine. So: keep our engine,
  swap in a validated aero **dataset**.

**Selected model: de Havilland DHC-2 Beaver** (Tjee & Mulder, *Stability and Control Derivatives
of the DHC-2 "Beaver" Aircraft*, TU-Delft Report **LR-556**, 1988; nonlinear implementation in
Rauw's **FDC** Simulink toolbox). It is the canonical flight-test-identified nonlinear GA aero
model — a complete polynomial set of stability & control derivatives, freely available, decades of
citations — and the Beaver is the archetypal floatplane. Fits "aero-only, seaplane, thoroughly
validated, third-party."

**Caveats (both honest, both manageable):**

1. The validated derivatives are for the **landplane** configuration; floats add drag / side-area
   increments. Options: treat the Beaver as the seaplane airframe and add a documented float-drag
   increment, or accept it as a stated limitation.
2. The polynomial model is valid over the **attached-flow envelope only** — no validated
   post-stall aero. This couples the plant swap to **framing B (impact-load)**, whose landing
   approach stays attached; it does *not* supply validated post-stall aero for framing A (stall
   recovery), which would still rest on a modeled overlay.

**Validity band vs. the landing phase (verified 2026-08-19).** LR-556 Table 3 is flight-validated
only over **30–55 m/s IAS**. The Beaver's float approach (~33.5 m/s / 75 mph, flaps 35°, ~400 fpm)
sits *inside* that band, but the flare and touchdown occur near stall (Vso ≈ 20 m/s / 45 mph;
touchdown ~20–22 m/s), **below the 30 m/s floor** — so the exact contact state the impact barrier
consumes is on *extrapolated* aero, and there is no validated stall aero. This is unavoidable (every
aircraft lands near stall) and is a **stated limitation, not a defect**: design the scenario approach
at ~33–35 m/s to keep the bulk of the descent inside the validated envelope, accept that the final
seconds extrapolate below 30 m/s, and say so. The impact barrier is mechanically unaffected — it
takes (V, γ, θ, sink) from the sim regardless of aero-model provenance.

**Airframe vs. impact model — reconciled (a floatplane is fine).** TN 1516 (verified against the
paper text, 2026-08-19) is a **single planing-surface** theory — it puts the *entire* seaplane
weight on one V-bottom, "a keeled float **or** hull," and says nothing about two floats, load
sharing, or lateral/asymmetric contact. A faithful twin-float impact model would be new,
*unvalidated* modeling that the paper does not support — so we **do not** extend it. Resolution:
keep the single-surface TN 1516 model and present the impact barrier as a **representative,
physically-grounded terminal safety constraint** (single-surface scope stated), *not* a validated
floatplane load prediction. This decouples the two things we actually validate — the **aero**
(Beaver) and the **CBF machinery** — and needs no twin-float extension.

**The flying-boat alternative was investigated and rejected (2026-08-19).** A single-hull flying
boat would make the airframe and TN 1516 the same object, but there is **no complete, flight-
validated, publicly-available flying-boat derivative set**. The best-validated ones (AG600, US-2,
Be-200) are OEM/state-proprietary; the public NACA flying-boat corpus is **hull-drag / hull-
increment / hydrodynamic (tank) data, not a vehicle model** — verified by inventorying NTRS
20050019300 (hull-in-presence-of-wing increments) and 20050019278 (planing-tail hull hydro). It
supplies a static-longitudinal hull skeleton at best; **all control derivatives, all damping
derivatives, the wing/tail buildup, and power effects would be estimated** (DATCOM), and even the
raw numbers are trapped in poorly-OCR'd 1940s scans. That collapses the "validated plant" claim to
"hull increments validated, the rest estimated" — defeating the pivot. So: **Beaver, not a flying
boat.**

**Impact-model parameters** (dead-rise β) can use a **representative float/hull dead-rise** — real
numbers replacing the placeholders — consistent with the representative-constraint framing above.

**Integration style.** Export the Beaver aero (base + derivatives) into the existing
`AeroTable`/`Dynamics` control-affine path — the engine, wind/wave forcing, and impact barrier all
stay; only the coefficient source changes. The CBF's *internal* model is a reduced version of the
same Beaver aero, so the plant/model gap is characterizable.

**Effort:** multi-week — key the Beaver aero into our table format, verify trim/handling against
the published model, build the CBF's reduced Beaver model, re-tune. High value: converts the
project's biggest weakness into its cleanest strength.
