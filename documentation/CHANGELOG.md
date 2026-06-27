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
