# Changelog

Chronological log of meaningful changes to this project. **Read the top entry at the start
of a session** to catch up on what your collaborator did since you last worked.

This doc's job — and how it differs from the others:

| Doc | Job |
|---|---|
| **`CHANGELOG.md`** (this file) | *When / who / what / why* — chronological, append-only. |
| `../TODO.md` | *Backlog* — granular open tasks/ideas to pick up. |
| `water_landing_cbf_design.md` | *Current state* — the living theory/approach/assumptions/decisions. |
| `water_landing_cbf_math.md` | *The math* — formal derivations. |
| `archive/` | Superseded historical docs. |

## Conventions

- **Newest entry on top.** Append above the existing entries; never rewrite past ones.
- One entry per working session (or per logically-grouped change). Keep it skimmable.
- Use the template below. Date is `YYYY-MM-DD`; author is a first name so collaborators can
  tell each other apart.

```
## YYYY-MM-DD — <author> — <short title>
**Branch/commit:** <branch> (<short-sha> if committed)
**What changed:** <the concrete change>
**Why:** <motivation / what problem it solves>
**Follow-ups / notes for collaborator:** <anything they should know or pick up>
**Files touched:** <key paths>
```

### Collaborator setup (one-time, per developer)

Each developer keeps their own local `CLAUDE.md` (gitignored, not shared via git). Add this
line to it so your Claude session reads this changelog at startup:

> At the start of each session, read `documentation/CHANGELOG.md` to see what changed since
> the last session before doing any work.

---

## 2026-06-29 — Brian — Altitude-sensor model + LPF, best-effort QP fallback, per-constraint slack weights, gain Monte-Carlo
**Branch/commit:** water-impact-cbf (sensor model + best-effort fallback + Lie test committed in `acdc27b`; per-constraint-slack refactor + gain study in the working tree)
**What changed:**
1. **Altitude-measurement sensor model + low-pass filter** (`lon_cbf_filter.hpp`, `lon_sim.cpp`).
   The CBF now acts on a noisy, low-pass-filtered altitude `h_filt(h + N(0,σ²))` while the
   plant still integrates the TRUE `h`. New YAML knobs: `h_meas_stddev` [m], `h_lpf_tau` [s]
   (first-order LPF `α=dt/(τ+dt)`, `f_c=1/2πτ`), `h_meas_seed` (RNG seed, exposed for the
   Monte-Carlo below). CSV gained `h_meas, h_filt`. Defaults (σ=0, τ=0) reproduce the
   perfect-altitude behavior exactly.
2. **Best-effort QP infeasibility fallback** (`lon_cbf_filter.cpp`). When the hard-constrained
   QP is infeasible the filter no longer softens-and-tracks-nominal (which deletes the
   constraint's effect); it solves a **minimum-violation** QP that penalizes each hard row's
   slack at a dominant weight (`kBestEffortHardPenalty=1e6`), spending the actuators to drive
   the hard-barrier violation toward zero. Last resort is the box-clamped nominal.
3. **QP objective refactor — per-constraint slack weights + identity control**
   (`lon_cbf_filter.{hpp,cpp}`, `lon_sim.cpp`, `lon_scenario.yaml`). Replaced the single global
   `slack_penalty` (and the impact Option-C **height-scheduled** slack) with per-constraint,
   YAML-settable weights `w_slack_{descent,airspeed,airspeed_upper,impact}`; dropped the control
   weights `w_de`/`w_Tddot` for an identity control cost `½‖u−u_nom‖²`. Hard rows still carry no
   slack. Kept the penalty **quadratic** `½wδ²` (not linear `wδ`): a pure-linear penalty leaves
   `P` singular on the slacks and OSQP then fails to converge on ~hundreds of feasible steps,
   firing best-effort spuriously (verified: identity+linear → 316 recoveries on the no-noise
   case; quadratic → 0, same landing).
4. **Lie-derivative verification test** (`test/test_lie_derivatives.cpp`, +CMake). 5 textbook
   systems (double/triple integrator, pendulum, `ẋ=x²`) with hand-computed `L_f^k b`,
   `L_gL_f^{r-1}b`, and HOCBF rows, cross-checking the Taylor-jet Lie engine and `hocbfRow`'s
   class-K (elementary-symmetric) coefficient placement to ~1e-10. Complements the existing
   finite-difference flow oracle.
5. **Gain robustness study** (`results/cbf_gain_{sweep,montecarlo}.{csv,md}`). A single-seed
   `c_descent` sweep picked `[1,2,128]`, but it's a razor-sharp single-realization optimum
   (overfit — `c3=128` ≈ a deadbeat pole at the timestep). A proper gain×seed **Monte-Carlo
   (40 common-random-number seeds)** showed `[1,2,128]` averages 0.55 m/s (not 0.024) across
   seeds and identified the robust optimum **`c_descent=[1,2,32]`** (two slow poles + one
   moderate-fast pole; mean 0.52, worst-case 1.17 m/s) — now set in `lon_scenario.yaml`.
**Why:** Study the CBF under a realistic noisy altimeter; make the infeasibility fallback
actually drive toward safety instead of deleting the constraint; give per-constraint control
over slack firmness; and put gain selection on a multi-seed footing instead of overfitting one
noise draw. The analysis figure `figures/lon_alt_noise_lpf.png` (ideal / filtered / no-CBF,
8 panels) tracks all of this.
**Follow-ups / notes for collaborator:**
- The per-constraint-slack refactor + gain study (`results/`) are **uncommitted**; commit when ready.
- Slack kept quadratic for OSQP robustness; if you want the L1/exact-penalty (linear) form it
  needs a slack rescale + looser OSQP tolerances.
- The impact barrier's height-scheduled slack (Option C) was **removed** for a constant
  `w_slack_impact`; re-add as a scheduled weight if you want cheap-aloft/firm-near-surface back.
- No gain keeps the descent barrier ≥ −0.5 in >~63% of seeds at σ=0.5 m — that's a
  noise/barrier-steepness limit (tune σ, `v_safe`, `τ`), not a gain-tuning one.
**Files touched:** `include/autoland/lon_cbf_filter.hpp`, `src/lon_cbf_filter.cpp`,
`src/lon_sim.cpp`, `data/lon_scenario.yaml`, `test/test_lie_derivatives.cpp`, `CMakeLists.txt`,
`results/cbf_gain_{sweep,montecarlo}.{csv,md}`, `figures/lon_alt_noise_lpf.png`.

---

## 2026-06-29 — Brian — Fix silent NaN-drop of the descent (sink-rate) barrier + un-mask the diagnostics
**Branch/commit:** water-impact-cbf
**What changed:** Three coupled fixes to a latent failure in the descent barrier
`b = V sinγ + √(v_safe² + 2·a_brk(V,γ)·h)` — the **only hard sink-rate guarantee**:
1. **Non-finite-proof radicand** (`hocbf.hpp`). When `v_safe² + 2·a_brk·h < 0` the raw
   `√` returned NaN, the filter's finiteness guard then **silently dropped the row** from
   the QP (coeffs→0, rhs→+∞), and the NaN was masked in the ψ-minima (`std::min` ignores
   NaN). `DescentBarrier::operator()` now floors the radicand with a smooth C∞ positive map
   `arg_pos = ½(arg + √(arg² + 4·eps_r²))` (new `eps_r = 0.01 m²/s²` field, a numerical
   regularizer — *not* a tuning knob). `b` stays finite and smooth, so the Lie jet and QP
   row stay well-posed and the barrier degrades to "flare as hard as possible" instead of
   vanishing.
2. **Dropped rows annunciated** (`lon_cbf_filter.{hpp,cpp}`). Removed the once-only
   `static bool warned`. The filter now exposes per-call `lastDroppedRows()`,
   `lastHardDropped()`, and `lastDescentInfeasible()` (the `a_brk ≤ 0` root-cause signal).
3. **NaN-guarded diagnostics** (`lon_sim.cpp`). ψ-minima go through a finiteness-guarded
   `gmin` with per-barrier non-finite-sample counters; the summary adds a **"CBF
   row-assembly faults"** line (so `recoveries=0` can no longer imply health — dropping the
   binding row makes the QP trivially feasible) and tags any ψ line with excluded samples.
   CSV gained `n_rows_dropped, hard_dropped, desc_infeasible` (appended at end).
**Why:** The trigger (`a_brk ≤ 0`, i.e. no available braking authority — low airspeed near
the surface) is **common-mode** with the soft `Vmin` barrier yielding and the elevator-only
impact barrier losing authority at low q̄, so all three sink protections could degrade
together with **no fault flag** and every health metric reading clean. Two reachable cases:
the **touchdown sample of every nominal run** (`h` integrates slightly < 0 with `v_safe²≈0`),
and any **`CL_max` recalibration downward** (an uncalibrated placeholder) that moves the
`a_brk=0` crossover into the flown envelope.
**Verified:** all **16 `[lon_cbf]` cases pass (143 assertions)**; the lone full-suite failure
(`test_cbf.cpp:218`, missing `data/AHAB_sweep.stab`) is **pre-existing** (confirmed identical
on HEAD with the patch stashed) and unrelated (body-axis path). Before→after:
- **Baseline** (CL_max=1.2): descent non-finite CSV cells **4→0**; touchdown sink
  **0.8023→0.8024 m/s** (unchanged to 4 s.f. — the `eps_r` bias is ~1e-8, no trajectory
  change); faults line reports `0 / 0 HARD / 0 a_brk≤0`.
- **Repro** (CL_max=0.30 → `a_brk≤0` over the descent): descent non-finite cells
  **2327 (100% of flight)→0**; the barrier now stays in the QP and **flares** — touchdown
  **1.63 m/s slam → 0.725 m/s** — and the summary loudly reports `a_brk≤0 on 5989 steps`
  with the *real* `descent ψ = −3.4 / −4.83` instead of a swallowed NaN. Previously this case
  ran with the barrier absent the entire flight while printing `recoveries: 0` and a positive
  `min ψ1`.
**Follow-ups / notes for collaborator:** This makes the failure **visible and graceful**, not
**safe** — when `a_brk ≤ 0` the envelope is genuinely infeasible and the barrier still can't
deliver a guarantee it physically lacks. The proper terminal fix is the **controllability /
handoff guard** (declare "cannot guarantee" + hold/abort rather than emit a control as if it
could) — new `TODO.md` item, related to the existing elevator-authority guard. CSV columns
were **appended** to preserve positions; name-keyed plot scripts are unaffected, but any
positional reader asserting an exact column count needs a +3 bump. A regression test asserting
the summary `min ψ` agrees with the per-sample CSV would lock this in.
**Files touched:** `include/autoland/{hocbf,lon_cbf_filter}.hpp`,
`src/{lon_cbf_filter,lon_sim}.cpp`, `test/test_lon_cbf.cpp`, `TODO.md`,
`documentation/CHANGELOG.md`.

## 2026-06-28 — Brian — Add hydrodynamic impact-load HOCBF barrier (NACA TN 1516)
**Branch/commit:** water-impact-cbf
**What changed:** New independent CBF-QP row that bounds the **peak CG load factor at
water touchdown**, derived from NACA TN 1516 (Milwitzky 1948). The barrier is
`b = (n_limit − n_peak(θ,γ,V)) + Φ(z)` with `n_peak = K0·ẏ₀²·Clf(κ)` (peak load factor,
g), the approach parameter `κ = sin τ/sin γ₀ · cos(τ+γ₀)`, and a height-relaxed term
`Φ(z) = Nb(1−e^{−z/zs})` that makes it touchdown-only. New files: an offline precompute
`scripts/precompute_impact_clf.py` (solves the transcendental eq 27 for the velocity
ratio, evaluates eq 25, verifies the `Clf(0)=0.6123` anchor) → generated table
`include/autoland/impact_clf_table.hpp`; the barrier itself in
`include/autoland/impact_barrier.hpp` (`clfLookup`, `ImpactLoadBarrier`,
`makeImpactLoadBarrier`). Wired into `LonCBFFilter` behind an `impact`/`impact_hard`
flag set with new config (`n_limit`, `beta`, `rho_water`, `Nb`, `zs`, `tau_keel`,
`z_gate`, `eps_g0`, `impact_slack_lo/hi`, `c_impact`) in `LonCBFConfig` +
`lon_scenario.yaml` + the `lon_sim` loader. Added an `exp()` overload to the Taylor-jet
Lie engine (`lie_taylor.hpp`) for `Φ`. The sim logs `b_impact,n_peak,kappa_imp` and the
degree-2 `psi1/psi2_imp` + `res_imp`, and prints the impact ψ-minima.
**Why:** Bound the structural *slam load* at hull contact — which depends on dead rise,
trim, and sink rate, not just sink rate (the descent barrier's kinematic cap). This is
the attitude-coupled contact-load constraint flagged as the §3.3 research item.
**Design choices:**
- **Degree-2 HOCBF, elevator-enforced.** The barrier is relative degree 2 via the
  elevator (θ→q→δe) and 3 via thrust (T→Tdot→Tddot). A single affine row can't mix the
  two (the elevator's degree-2 entry pushed to degree 3 yields non-affine u²/u̇ terms),
  so it's built as a clean degree-2 row (`barrierLie<2>` + 2 class-K gains): thrust's
  column in `L_gL_f b` is exactly 0 and drops out; the flare enforces it. Thrust still
  bounds impact load via the existing descent/sink-rate barrier (spec §4).
- **Frozen K0 + local-affine `Clf(κ)`** at the eval point (mirrors `makeDescentBarrier`),
  so the templated barrier is smooth with no cbrt/pow in the Taylor path while the
  attitude coupling stays live through `κ`.
- **Activation:** Option A height-relaxed `Φ(z)` (primary) + Option C height-scheduled
  slack (`impact_slack_lo→hi`) backstop. Plus a **model-validity gate** — the row is only
  assembled below `z_gate` AND while descending with positive trim (NACA TN 1516's
  validity domain); outside it the prediction is meaningless (κ<0) and would feed a
  huge-coefficient row to the QP.
- **Placeholders:** `n_limit=3 g`, `beta=22.5°`, `rho_water=1000` (fresh), `Nb=10`,
  `zs=2 m`. Soft + enabled by default and sized **non-binding** (additive safety, zero
  trajectory change), mirroring how the upper-airspeed barrier was introduced.
**Verified:** **34/35 tests pass** — incl. 6 new (clfLookup monotone/clamp, closed-form
value, relative-degree thrust-drops-out, drift flow-oracle for the degree-2 stack +
`exp()`, Φ height-relaxation, hard-enforcement under a violating nominal). The 1
remaining failure (`#19`, `test_cbf.cpp`) is **pre-existing** (old body-axis path, throws
unrelated to this work). Default sim run **unchanged**: touchdown 0.0137 m/s, 0
recoveries; `impact: min psi1=5.84, psi2=14.11 ≥ 0`; over the active window
`n_peak ≤ 0.17 g « n_limit=3`, `κ∈[0.001,1.149]`, slack ≡ 0 (row never binds).
**Follow-ups / notes for collaborator:** `n_limit`/`beta`/`Nb`/`zs`/`c_impact` need
calibration/tuning (new `TODO.md` items). The full mixed-degree (2/3) treatment and a
predictive/backup-set CBF remain open (§3.3). CSV gained 6 columns — name-keyed plot
scripts are unaffected. Source paper added at `documentation/19930082553.pdf`.
**Files touched:** `include/autoland/{lie_taylor,impact_barrier,impact_clf_table,
lon_cbf_filter}.hpp`, `src/{lon_cbf_filter,lon_sim}.cpp`, `scripts/precompute_impact_clf.py`,
`data/lon_scenario.yaml`, `test/test_lon_cbf.cpp`,
`documentation/{impact_load_barrier_spec,water_landing_cbf_design,CHANGELOG}.md`, `TODO.md`.

## 2026-06-27 — Jack — Add maximum-airspeed (over-speed) HOCBF barrier

**Branch/commit:** corbin-dev
**What changed:** Wired the previously-stubbed `AirspeedUpperBarrier` (`b = V_max − V`) into the
CBF-QP filter — a symmetric over-speed / high-energy-impact guard alongside the stall barrier.
New config: `airspeed_upper` / `airspeed_upper_hard` flags, `Vmax_air`, `c_airspeed_upper`
(`LonCBFConfig` + `lon_scenario.yaml` + loader). The sim now logs `b_airspeed_upper`, its ψ-minima
(console), and a `res_airup` residual column. Because the barrier is a function of `V` only, it
reuses the lower barrier's degree-3 machinery unchanged: its drift Lie stack and control row are
the **sign-flip** of the stall barrier's (`L_f^k b_up = −L_f^k b_lo` for `k≥1`).
**Why:** Bound high-energy water impact / structural limits — slam loads grow with impact speed
on a flying-boat hull. Closes the `TODO.md` "Upper airspeed barrier" item.
**Design choices:** `V_max` is a placeholder (`27` m/s ≈ 1.5·V_app) set generously so it is
protective but **non-binding** on the current approach; soft by default (mirrors the stall
barrier); enabled by default (additive safety, zero trajectory change).
**Verified:** all **29 tests** pass — incl. 2 new (oracle SECTION revalidating the upper
barrier's drift stack; a sign-flip identity test vs the lower barrier; a hard-enforcement test
under an over-speed nominal). Default sim run **unchanged** (touchdown 0.0137 m/s, 0 recoveries);
`b_airspeed_upper` stays in [8.67, 10.48] (V_max=27), upper ψ-minima ≈ 17/33 — confirmed
non-binding. Spot-check: lowering `V_max` / setting it hard makes the row actively cap `V`.
**Follow-ups / notes for collaborator:** `V_max` needs calibration (new `TODO.md` item, under
Modeling/data); both airspeed barriers remain soft (the "harden the airspeed barrier" item now
applies to the pair). CSV gained a column — name-keyed plot scripts are unaffected.
**Files touched:** `include/autoland/{hocbf,lon_cbf_filter}.hpp`, `src/{lon_cbf_filter,lon_sim}.cpp`,
`data/lon_scenario.yaml`, `test/test_lon_cbf.cpp`,
`documentation/{water_landing_cbf_math,water_landing_cbf_design,CHANGELOG}.md`, `TODO.md`.

## 2026-06-27 — Jack — Steady level-flight trim for the initial condition

**Branch/commit:** corbin-dev
**What changed:** Added `lonTrim()` (`src/lon_sim.cpp`) — a Newton solve on the lon EOM for a
steady **level-flight** equilibrium (`γ=0`, `V̇=γ̇=q̇=0`, mirroring the body-axis `trim()`). The
sim IC now uses it, so the run starts cruising level and the nominal pushes over into the −3°
approach (more realistic; no body-axis-seed startup transient).
**Why:** Asked for a realistic IC; also a step on the "true longitudinal trim" backlog item.
**Finding (corrects an earlier note):** the level IC did **not** change the residual descent
`ψ < 0` — so that dip is **not** the initial condition. It is a **terminal** effect in the last
~0.4 s as `h → 0` at touchdown; the actual barrier `b` stays ≥ 0 in flight (the realized
trajectory is safe). Design doc §6/§8.4 updated accordingly.
**Verified:** all 27 tests pass; default run lands at 0.014 m/s within budget, 0 recoveries,
`b_descent ≥ 0.10` airborne. (Also generated a class-K gain comparison figure showing the four
`b` barriers; stiff `[100,100,100]` gains under-flare and violate `b_descent` at touchdown,
gentle `[1,1,1]` stay safe.)
**Files touched:** `src/lon_sim.cpp`, `documentation/water_landing_cbf_design.md`, `TODO.md`,
`documentation/CHANGELOG.md`.

## 2026-06-27 — Jack — Force-based `a_brk(V,γ)` for the descent barrier; real deck by default

**Branch/commit:** corbin-dev
**What changed:**
- Replaced the constant `a_brk = 3.0` with a speed/path-angle-dependent
  `a_brk(V,γ) = (ρV²S/2m)[CL_max cosγ − CD_maxlift sinγ] − g` (lift + drag + gravity) in
  `DescentBarrier`, built by `makeDescentBarrier()` and differentiated exactly by autodiff.
  **Thrust deliberately omitted** — it would drop the barrier's relative degree below 3 and
  break the augmented HOCBF alignment (backlog, tied to §3.3).
- `C_L,max` is a new config field (`CL_max`, placeholder 1.2); `C_D,maxlift` is computed from the
  frozen aero at the extrapolated max-lift α. Added a positivity warning in the filter.
- Pointed `lon_autoland_sim`'s default aero deck at the **real `AHAB_combined.stab`** (was the
  `example.stab` placeholder — the source of earlier "violations").
- Tests: oracle test now auto-revalidates the new barrier's drift stack; replaced the descent
  closed-form authority check (invalid for state-dependent `a_brk`) with an `a_brk(V,γ)` formula
  value test. **All 27 pass.**
**Why:** Make the soft-landing envelope honestly track available lift (tighten at low speed,
loosen at high) instead of a hand-picked constant.
**Verified (real AHAB deck):** touchdown sink **0.011 m/s** (within the 0.1 budget), **0 QP
recoveries** — a big improvement over the constant-`a_brk` baseline (sink 2.18 m/s, 11
recoveries on the same deck). The residual descent `ψ < 0` at the start is the non-trim initial
condition, **not** control saturation (smaller class-K gains made it *worse*).
**Follow-ups / notes for collaborator:** `CL_max` needs calibration; thrust term + lag-aware
`a_brk` + a true lon trim are in `TODO.md`. The unit tests still use `example.stab` (deck-
agnostic math); removing `example.stab` is a backlog item.
**Files touched:** `include/autoland/{hocbf,lon_cbf_filter}.hpp`,
`src/{lon_cbf_filter,lon_sim}.cpp`, `apps/lon_autoland_sim.cpp`, `data/lon_scenario.yaml`,
`test/test_lon_cbf.cpp`, `documentation/{water_landing_cbf_math,water_landing_cbf_design}.md`,
`TODO.md`, `documentation/CHANGELOG.md`.

## 2026-06-27 — Jack — Add `TODO.md` backlog; retire stale README TODO

**Branch/commit:** corbin-dev
**What changed:** Added a root-level `TODO.md` (granular backlog, grouped by area) seeded from
the design doc's open items. Replaced the stale README "⚠️ TODO" section (placeholder-data
items resolved by the AHAB refactor) with a pointer to it. Wired `TODO.md` into the
collaboration convention (CLAUDE.md startup note + changelog doc-map).
**Why:** Give a low-friction, collaborator-visible place to capture tasks/ideas, distinct from
the changelog (done) and the design doc §9 (rationale for big threads).
**Files touched:** `TODO.md`, `README.md`, `CLAUDE.md` (gitignored), `documentation/CHANGELOG.md`.

## 2026-06-27 — Jack — Verify drift Lie stack (`L_f³b`) + make build MSVC-portable

**Branch/commit:** corbin-dev
**What changed:**
- Added a test-only **finite-difference flow oracle** that independently cross-checks the
  drift Lie stack `{b, L_f b, L_f²b, L_f³b}` (descent + airspeed barriers) against the
  bespoke Taylor-jet engine in `lie_taylor.hpp`. Engine vs oracle agree to ~2.5e-4 on
  `L_f³b`; a seeded 1% engine error is caught. Closes the prior "L_f³b not independently
  verified" gap (only the control authorities were checked before).
- Made the build work under **MSVC** (the only compiler on this Windows machine): guarded the
  GCC/Clang `-Wall -Wextra` flags behind `if(NOT MSVC)` and defined `_USE_MATH_DEFINES` for
  `M_PI`. Strictly additive — does not change the GCC/Clang build.
**Why:** Cheap insurance that the high-order Lie machinery is correct (the autodiff *primitive*
was never in doubt; the custom Picard/Taylor *iteration* was the untested part). Build fix was
needed to compile/run the suite locally on Windows.
**Verified:** all 26 tests pass (`ctest -C Release`); confirmed the new test fails on a 1%
perturbation, then reverted.
**Follow-ups / notes for collaborator:** to build on Windows/MSVC, configure with the Visual
Studio generator (`cmake -S . -B build -G "Visual Studio 17 2022" -A x64`) — Ninja hits an
OSQP `osqp.lib` multiple-rules error.
**Files touched:** `test/test_lon_cbf.cpp`, `CMakeLists.txt`,
`documentation/water_landing_cbf_design.md`, `documentation/CHANGELOG.md`.

## 2026-06-27 — Jack — Set up `corbin-dev` branch + collaboration docs

**Branch/commit:** corbin-dev (off `lon-cbf-water-landing` @ `2b610a3`)
**What changed:**
- Created `corbin-dev` from the latest `lon-cbf-water-landing` (the augmented-longitudinal
  HOCBF water-landing module with exact autodiff/Taylor-jet Lie derivatives).
- Gitignored personal Claude context (`CLAUDE.md`, `PROJECT_CONTEXT.md`).
- Added this changelog and a living design doc (`water_landing_cbf_design.md`) that supersedes
  the dated implementation notes.
- Moved historical docs (`water_landing_cbf.md`, `water_landing_cbf_implementation_notes.md`)
  into `documentation/archive/`.
**Why:** Establish a shared working branch and a lightweight way for two parallel
Claude-using collaborators to stay in sync, plus a single maintained source of truth for the
design instead of a point-in-time session snapshot.
**Follow-ups / notes for collaborator:**
- Add the one-line changelog pointer to your own local `CLAUDE.md` (see Collaborator setup).
- Open design questions captured in `water_landing_cbf_design.md`: force-based `a_brk(V)`,
  the §3.3 contact-force barrier, and an independent `L_f³b` drift check.
**Files touched:** `.gitignore`, `documentation/CHANGELOG.md`,
`documentation/water_landing_cbf_design.md`, `documentation/archive/`.
