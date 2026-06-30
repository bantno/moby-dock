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

## 2026-06-29 — Jack — Math-spec nomenclature + real nominal-controller section

**Branch/commit:** corbin-dev (committed this session)
**What changed:** Added a **Nomenclature** section to `water_landing_cbf_math.md` — grouped
symbol / units / code-identifier tables (state & control, aircraft & aero, environment, the
descent/airspeed/contact-force/impact-load barriers, and the CBF/QP machinery), plus an
overloaded-symbol callout (α, W, u, A, c̄). Rewrote **§2 (Nominal Control Strategy)** from the
stale PX4-TECS placeholder to the controller actually in the code: a constant-thrust powered
approach (PD on the augmented thrust state) + a γ→θ→δe cascade (outer PI, inner PD), with trim
feedforward seeding and a runtime-gain table tied to `lon_scenario.yaml`.
**Why:** Make the math spec presentable to an advisor and consistent with the implementation —
the TECS text never matched the code, and symbols were undefined for an external reader.
**Follow-ups / notes for collaborator:** §3.1's descent control-affine "closed form" is still
labelled reference-only (live code differentiates `a_brk(V,γ)` via autodiff); §3.3 vs §3.5 still
read as two contact-force barriers (3.5 is the implemented NACA TN 1516 one). Separate
experimental work this session — energy-reachability CBF, nominal hold-off flare, parameter-sweep
tool (`scripts/sweep.py`) — lives on branch `experiments-corbin` (pushed), not here. Key finding
there: touchdown airspeed is drag-limited (clean VSPAero deck can't decelerate via pitch); real
speed reduction needs landing-config drag.
**Files touched:** documentation/water_landing_cbf_math.md, documentation/CHANGELOG.md

## 2026-06-28 — Jack — Merge `water-impact-cbf` into `corbin-dev`

**Branch/commit:** corbin-dev (merge of `origin/water-impact-cbf` `8e2770a`)
**What changed:** Fetched and merged the collaborator's `water-impact-cbf` branch — Brian's
hydrodynamic impact-load HOCBF barrier (NACA TN 1516); see his entry below for the substance.
The branch forked from `eaa3076` (before the `AHAB_combined` consolidation), so this is a real
3-way merge. **All code/data merged cleanly**; the only conflicts were the two append-only docs
(`CHANGELOG.md`, `TODO.md`), resolved by **keeping both sides' entries** (Brian's impact-barrier
work folded in alongside Jack's deck-consolidation + over-speed-barrier work). The impact module
adds two new aero decks of its own (`AHAB_alpha_beta_sweep.stab`,
`AHAB_control_surface_effectiveness.stab`) — distinct filenames, no collision with the
`AHAB_combined` consolidation.
**Why:** Bring the deferred §3.3 contact-force / slam-load barrier (the highest research-value
thread) onto `corbin-dev`.
**Follow-ups / notes for collaborator:** verify build + tests on the merged tree (the impact
branch reported 34/35, with `#19` a pre-existing body-axis failure unrelated to this work).
**Files touched (merge):** conflict resolution in `documentation/CHANGELOG.md`, `TODO.md`; all
other impact-module files merged without conflict.

## 2026-06-27 — Jack — Merge `origin/main`; reconcile deck cleanup with collaborator

**Branch/commit:** corbin-dev (`106159e` merge; `73098b9` consolidation)
**What changed:** Fetched and merged `origin/main` into `corbin-dev`. Discovered the collaborator
(bantno) had independently pushed `eaa3076 "remove deprecated files"` (already merged to `main` via
PRs #1–#2) that deleted the same two decks (`AHAB_sweep.stab`, `AHAB 2.stab`) **plus four stale
figure PNGs** (`autoland_response`, `cbf_detail`, `cbf_on_vs_off`, `gain_compare`) — but **did not
repoint the code**, leaving `main`'s `autoland_sim.cpp`/`test_cbf.cpp` referencing the deleted
`AHAB_sweep.stab` (broken build/tests on `main`). The merge was clean: the two `.stab` deletions
resolved as both-deleted (no conflict), the 4 PNG deletions came in from `main`, and my repoint
commit (`73098b9`) supplies the fix `main` was missing. Net: `corbin-dev` now removes the
deprecated decks **and** builds/tests green.
**Why:** Keep the branch current with `main` and fold the two parallel deck-cleanup efforts into one
consistent state.
**Verified:** clean merge, no conflicts; rebuilt + **29/29 tests pass** (Release) on the merged tree.
Only `AHAB_combined.stab` + `example.stab` remain in `data/`; `figures/` is now just `README.md`.
**Follow-ups / notes for collaborator:** `figures/README.md` still documents the 3 deleted PNGs
(its "Regenerate" recipe — repointed to `AHAB_combined` — is still valid). Regenerate them on the
new deck or trim the README to the recipe (open question). Pushing `corbin-dev` up fixes the broken
deck references currently on `main`.
**Files touched (merge):** removed `figures/*.png` (from `main`); no source changes in the merge
commit itself.

## 2026-06-27 — Jack — Consolidate on the `AHAB_combined.stab` aero deck

**Branch/commit:** corbin-dev (`73098b9`)
**What changed:** Removed the two deprecated aero decks `data/AHAB_sweep.stab` (α −20…0, β 0…20)
and `data/AHAB 2.stab` (5-case fragment), and repointed everything that referenced them onto the
real full deck `data/AHAB_combined.stab` (α −20…20 × β 0…20, single Mach 0.059):
- `apps/autoland_sim.cpp` default deck `AHAB_sweep.stab` → `AHAB_combined.stab` (a superset in α,
  so the body-axis app now runs on-grid at positive α where it used to extrapolate).
- `test/test_cbf.cpp` "CBF runs on the real linearized longitudinal model" likewise.
- `figures/README.md` regen commands + caption.
**Kept `example.stab`** intentionally — the deck-agnostic unit tests (`test_aero_table`,
`test_linear_model`, `test_trim`, `test_lon_cbf`) still use it as a small synthetic fixture
(symmetric β ∈ [−10,10]); migrating them to the real deck is a separate, riskier change (β grid is
no longer symmetric; real airframe is near-neutral/slightly unstable so some sign-sense assertions
could legitimately flip). Left as a TODO, not done here.
**Why:** `AHAB_combined.stab` is the authoritative vehicle deck; the others were outdated/incorrect.
**Verified:** all **29 tests** pass (Release); body-axis `autoland_sim` runs to touchdown on the
combined deck (exit 0). No code besides the deck path changed.
**Follow-ups / notes for collaborator:** `main` had *already* removed these decks (collaborator's
`eaa3076`) but without repointing the code — see the merge entry above. `TODO.md` "Remove
example.stab" item revised to reflect the keep-decision.
**Files touched:** `apps/autoland_sim.cpp`, `test/test_cbf.cpp`, `figures/README.md`,
`data/` (removed `AHAB_sweep.stab`, `AHAB 2.stab`), `TODO.md`, `documentation/CHANGELOG.md`.

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
