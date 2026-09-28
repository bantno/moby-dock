## What was built, in one paragraph

A linear Airy surface-wave field as **plant-side truth** with a **wave-blind CBF filter**
(your two decisions). One JONSWAP spectrum implementation covers every requested sea:
`gamma = 1` degenerates exactly to the two-parameter Bretschneider spectrum that
**NATO STANAG 4194** prescribes with its Annex D sea states (the mil-spec tier),
`gamma ~ 3.3` is the fetch-limited developing sea for lakes with (Hs, Tp) from the
**USACE CEM** fetch-growth laws, and `regular: true` collapses to one deterministic
component for tests/demos. Realization is the standard sum-of-sinusoids with seeded
random phases and within-bin frequency jitter (St. Denis & Pierson 1953; Fossen ch. 8;
ITTC 7.5-02-07-01.1), deep-water dispersion `k = w^2/g`.

## Files

| File | Status | What |
|---|---|---|
| `include/autoland/water_waves.hpp` | new | spectrum + realization + field evaluation (header-only, no sim state) |
| `include/autoland/impact_barrier.hpp` | modified | `impactK0` / `impactNPeakExact` factored out of the factory |
| `include/autoland/lon_sim.hpp` | modified | `LonScenario.waves`, `LonTouchdown` contact record |
| `src/lon_sim.cpp` | modified | YAML parse, altimeter re-reference, touchdown at `h = eta`, ground track, CSV, console |
| `test/test_water_waves.cpp` | new | 4 test cases, tag `[waves]` |
| `CMakeLists.txt` | modified | registers the test file |
| `data/lon_scenario.yaml` | modified | `waves:` schema block, **disabled** (flat-water default preserved) |
| `data/lon_landing_waves_lake.yaml` | new | JONSWAP lake demo (STANAG SS3 ocean alternative in comments) |
| `data/lon_landing_wave_regular.yaml` | new | locked regular-wave ctest scenario |
| `scripts/plot_wave_landing.py` | new | 4-panel wave-landing figure |
| `figures/lon_landing_waves.png` | new | lake-demo output |
| `README.md`, `TODO.md`, `documentation/CHANGELOG.md` | modified | run recipe; smooth-water TODO rewritten; changelog entry |

## The load-bearing code, with pointers

**Spectrum** — `water_waves.hpp:61`. `S_B = (5/16) Hs^2 wp^4 / w^5 exp(-1.25 (wp/w)^4)`;
JONSWAP multiplies by `gamma^b` with `sigma = 0.07/0.09` and the Hs-preserving
`A_gamma = 1 - 0.287 ln(gamma)` (DNV-RP-C205). `gamma <= 1` returns pure Bretschneider.

**Realization** — `makeWaveField`, `water_waves.hpp:139`. Band `[0.5, 5] * wp`
(tail outside holds < 0.3% of m0), `A_i = sqrt(2 S(w_i) dw)` at a frequency drawn
uniformly within each bin (non-periodic realization), seeded `mt19937` so a scenario
is one locked, reproducible sea.

**Hull-mean slope** — `slopeMean`, `water_waves.hpp:132` + `contact_len` at `:88`.
The point slope's spectrum is `k^2 S(w)` — dominated by sub-hull-length ripples that
produced spurious ±25 deg "faces" in the first demo run. The contact diagnostics use
the chord slope across the wetted length (0.4 m default) instead: components shorter
than the keel average out, exactly as the hull bridges them. **This is my judgment
call, not from a source — review it.**

**Sim wiring** — `src/lon_sim.cpp`:
- `:266` YAML `waves:` parse (direction as string `head|following`).
- `:330` field built once per run; `h_filt` seeded at initial *clearance*.
- `:377-384` per-step `eta`/`slopeMean`; **altimeter measures `h - eta`** (radar
  clearance). The filter flies the surface-relative altitude but keeps its flat-water
  barrier model — that is the whole "wave-blind" contract. The nominal controller
  still gets the true state (pre-existing convention, unchanged).
- `:548` touchdown test is now `X[LH] <= eta_now` (was `<= 0`).
- `:551-566` contact record: `sink_rel = sink + eta_x*xdot + eta_t` (closure of
  `h - eta`), and the TN 1516 truth pair `n_peak_flat` / `n_peak_wave` — the wave one
  tilts tau/gamma0 by `atan(slopeMean)` and closes at `sink_rel` (the spec's own
  rough-water referencing recommendation).
- `:579` earth-frame ground track `x_pos` advanced by trapezoid **outside** RK4 —
  legitimate because waves never force the airborne dynamics; x feeds only
  `eta(x, t)` and diagnostics. Error ~mm against ~5 m wavelengths.
- CSV appends 3 columns at the end: `x, eta, eta_slope` (hull-mean slope).

**Impact-barrier refactor** — `impact_barrier.hpp:74` (`impactK0`) and `:94`
(`impactNPeakExact`). Pure extraction: the factory now calls `impactK0` with the same
expressions in the same order (the 1e-3-tight frozen-gradient test still passes).
`impactNPeakExact` is the *unfrozen* evaluator the touchdown truth uses; it mirrors
the factory's kappa/tau clamps.

## Flat-water behavior is unchanged

- `waves.enabled` defaults false; `WaveField::eta` returns exactly `0.0` when disabled,
  so `h_meas = X[LH] - 0.0` and the touchdown test degrade to the old code path.
- All 53 pre-existing tests pass; the flagship flat touchdown still reads
  t=68.54 s, sink=0.2095, V=13.371, theta=2.997 — the documented baseline.
- One mid-session failure of the flagship test (`recoveries == 0`) was the known
  OSQP feasibility-recovery jitter (wall-time-adaptive rho), 8/8 green on rerun and
  green in the final full run. Not introduced by this change, but be aware it exists.

## Verification

```
cd build && ctest            # 57/57 (4 new [waves] cases)
```
- Spectrum: closed form at wp; gamma=1 == Bretschneider across the band;
  m0 = integral S dw == Hs^2/16 (1% for gamma=1; 6% for gamma=3.3 — the A_gamma
  normalization is approximate by construction).
- Realization: sum A_i^2/2 == band-limited m0 (8%, within-bin jitter scatter);
  k = w^2/g exact per component; same-seed reproducible, different-seed differs;
  time-variance of eta at a fixed point ~ sum A_i^2/2 (25%, single realization).
- Regular wave: exact eta/slope/etaDot closed forms, the advection identity
  `eta_t = -c_p eta_x`, crest-riding for both directions, disabled-field zeros.
- End-to-end (`lon_landing_wave_regular.yaml`): touchdown AT the surface
  (`|h_td - eta_td| < 2 cm`, `|eta_td| > 5 mm` so it isn't a zero-crossing),
  envelope holds, truth pair populated and differing.

Demos (CSV in git-ignored `results/`):
```
./build/lon_autoland_sim data/AHAB_combined.stab data/aircraft.yaml \
    data/lon_landing_waves_lake.yaml results/lon_waves_lake.csv
python3 scripts/plot_wave_landing.py results/lon_waves_lake.csv figures/lon_landing_waves.png
```
Lake result: filter meets its flat-water spec (sink 0.209 m/s, V 13.37 — matches the
flat baseline) while the truth is contact on a **7.7 deg rising face at 2.36 m/s
closure**: `n_peak` wave-referenced **18.8 g vs 0.32 g** flat-referenced (~59x, ~6x
over `n_limit = 3`). Regular-wave variant: 3.6 deg face, 1.25 m/s, 5.3 g (~16x).

## Caveats I want you to see (also in CHANGELOG/TODO)

1. **One seeded sea = one draw.** The 59x headline depends on where in the encounter
   phase the keel lands. A seed/phase sweep (suite-style) is the natural follow-up
   before quoting statistics.
2. **`n_peak_wave` leaves the TN 1516 validated range** on steep faces: tau relative
   to the surface goes negative and clamps at 1 deg. Read it as "far outside the
   flat-water design point", not a calibrated load.
3. Deep-water dispersion only (fine for lake chop; shallow water would need
   `w^2 = gk tanh(kd)`).
4. The altimeter measures the *point* surface under the aircraft; the touchdown
   diagnostics use the *hull-mean* slope. Both choices are argued in comments, both
   are reviewable judgment calls.
5. Stage 2 (wave-aware barrier: tau/gamma0/Phi referenced to the surface inside the
   filter) is deliberately not attempted — recorded in TODO.
