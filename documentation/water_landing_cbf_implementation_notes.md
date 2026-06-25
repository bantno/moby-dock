# Augmented Longitudinal CBF-QP: Implementation Notes & Session Summary

**Date:** 2026-06-25
**Scope:** Implementation of the augmented-longitudinal Control-Barrier-Function (CBF) safety
filter specified in `water_landing_cbf_math.md`, **excluding section 3.3** (contact-force
barrier, deferred), plus verification, closed-loop simulation, and analysis.

This document is the implementation companion to the math spec. It records what was built, how
the math was verified, the results, and — importantly — the assumptions/simplifications and the
analysis of *why* the closed-loop behaves as it does.

---

## 1. What was built

A new, self-contained **longitudinal** module (`lon_*`) living alongside the existing 11-state
body-axis code, which it leaves untouched. New files:

| File | Role |
|---|---|
| `include/autoland/lie_taylor.hpp` | Exact Lie-derivative engine (flow Taylor-jet + autodiff) |
| `include/autoland/lon_augmented.hpp` / `src/lon_augmented.cpp` | 7-state augmented EOM + wind-frame aero adapter |
| `include/autoland/hocbf.hpp` | Barriers, per-barrier Lie provider, generic HOCBF constraint builder, actuator rows |
| `include/autoland/lon_cbf_filter.hpp` / `src/lon_cbf_filter.cpp` | The §4 QP over OSQP, with feasibility recovery |
| `include/autoland/lon_nominal.hpp` | Nominal controller (constant thrust + cascade γ→θ→pitch) |
| `include/autoland/lon_sim.hpp` / `src/lon_sim.cpp` | Self-consistent RK4 sim + CSV logging + diagnostics |
| `apps/lon_autoland_sim.cpp` | Executable |
| `data/lon_scenario.yaml` | Scenario/config |
| `test/test_lon_cbf.cpp` | Unit tests |
| `scripts/plot_lon_results.py`, `scripts/plot_lon_compare.py` | Plotting |

New dependency: **autodiff** (autodiff.github.io, header-only C++17), fetched via CMake
FetchContent and exposed as an INTERFACE include target (its own CMake `find_package(Eigen3)`
is bypassed since we vendor Eigen as headers).

### Augmented system (math spec §1)
- State `X = [h, V, γ, θ, q, T, Ṫ]` (NXA=7); control `U = [δe, T̈]` (NUA=2).
- Thrust augmented with two integrators (`T`, `Ṫ` states; `T̈` control) to align the relative
  degrees of elevator and thrust.

### Nominal controller (replaces the spec's TECS, per directive)
- **Thrust:** constant setpoint with a PD law on the augmented thrust state,
  `T̈_nom = Kp_T (T_set − T) − Kd_T Ṫ`.
- **Elevator:** cascade holding a constant negative flight-path angle γ — outer PI
  `γ → θ_cmd` (limited), inner PD `θ → δe` with pitch-rate damping.

### Barriers implemented
- Descent-rate (§3.1), relative degree 3.
- Airspeed (§3.2), relative degree 3.
- Min/Max thrust actuator barriers (§3.4), relative degree 2 (closed form).
- **§3.3 contact-force barrier deferred** (not implemented).

---

## 2. The exact Lie-derivative engine (no finite differences)

The HOCBF rows need `L_f^k b` and the control row `L_g L_f^{r-1} b` over the augmented
dynamics. Per directive these are computed **exactly**, not by central differences.

**Approach:** the flow Taylor-jet. For `Ẋ = f(X)`, the Taylor coefficients of `t ↦ b(X(t))`
give `L_f^k b = k! · [t^k] b(X(t))`. We build the order-r jet by Picard iteration in
truncated-Taylor arithmetic (each pass evaluates `f` and fixes one more order). The barrier and
EOM are templated on the coefficient scalar.

- Drift Lie derivatives use scalar `double`.
- The control row `L_g L_f^{r-1} b = ∇(L_f^{r-1}b)·g` is the directional derivative of the
  field along a g-column; computed by promoting the coefficient scalar to `autodiff::dual` and
  seeding the jet base point along that column.

A **hand-rolled `Taylor<N,S>`** is used (not autodiff's `real`) because `real` requires an
arithmetic coefficient type and so cannot nest over `dual` for the control row; also
autodiff's `along()` gives fixed-direction derivatives `(v·∇)^k b`, which are **not** the Lie
derivatives `L_f^k b` (they differ from order 2 because `f` varies along the flow).

**Aero for autodiff:** `AeroLocal` freezes the VSPAero table slopes at the evaluation point
(the table model is piecewise-linear in α within a cell), so the templated EOM are smooth and
the derivatives are exact w.r.t. that local affine model.

---

## 3. Aero frame (verified)

The spec's EOM are in wind/path axes (need `L`, `D`, `∂L/∂α`, `∂D/∂α`), while VSPAero's primary
data is body-axis. We **rotate the body-axis coefficients (and their derivatives) into the wind
frame analytically** — standard 6DOF practice — rather than consuming VSPAero's pre-rotated
`CL/CD` columns.

Verified transform (OpenVSP Google Group, Rob McDonald; numerically checked against
`example.stab`), applied to the **file-frame** coefficients (before dynamics' `kFrameSign`
flip), β=0:

```
C_D = CFx·cosα + CFz·sinα + parasite_CD0
C_L = −CFx·sinα + CFz·cosα
```

At the α=−10°/β=−10° node this reproduces the file's `CL=−0.5834, CD=+0.0218` exactly. Note the
trap: rotating the *post-flip* (flight-dynamics) coefficients instead gives the wrong sign. The
pitch moment stays body-axis (`CMm == CMy`). The wind-frame rotation lives entirely in the new
module (`makeAeroLocal` reads a raw `AeroTable` lookup); **`dynamics.cpp` was not modified.**

---

## 4. Math verification (the "check all the math" task)

Verified the spec's control-authority closed forms against exact autodiff:

- **§3.1 descent thrust authority `sinθ/m`** and **§3.2 airspeed thrust authority `cosα/m`** —
  match exactly (to ~1e-10) with the full aero. ✓
- **§3.4 actuator HOCBFs** — exact. ✓
- **§3.1/§3.2 elevator authorities** — match the spec exactly **only in the no-rate-aero limit**
  the spec implicitly assumes (`L=L(α,V)`, `D=D(α,V)`).

**Finding:** the spec's elevator-authority closed form **omits the aerodynamic pitch-rate
dependence** of lift/drag (`∂C_L/∂q̂`, `∂C_D/∂q̂`). The autodiff differentiates the full
tabulated aero and so captures those terms; the two differ by a few percent. The autodiff
implementation is the more complete quantity. This is recorded as a dedicated test
(`full aero elevator authority diverges from the doc`).

---

## 5. Results

All 25 unit tests pass (23 prior + the new `[lon_cbf]` set). The closed-loop sim
(`lon_autoland_sim`, default `lon_scenario.yaml`, `example.stab`):

- **Touchdown sink = 0.006 m/s** (hull budget `v_safe = 0.6`), V = 14.2 m/s (≥ Vmin), θ = 3.8°,
  **zero QP feasibility recoveries**.
- **HOCBF nested-function minima** over the run (the real forward-invariance condition):
  `descent ψ1=1.16, ψ2=1.50`; `airspeed ψ1=1.13, ψ2=1.30` — all ≥ 0.
- **All four barriers stay strictly positive:** `b_descent` min +0.59, `b_airspeed` +0.74,
  `b_thrust_min (T)` +1.83, `b_thrust_max (Tmax−T)` +4.17.

### CBF on vs off (nominal-only)
Same nominal controller, identical trajectory until the last ~3 m, then:

| | touchdown sink | vs `v_safe`=0.6 | θ |
|---|---|---|---|
| CBF off (nominal only) | **0.879 m/s** | EXCEEDED | −1.3° (no flare, glides straight in) |
| CBF on | **0.006 m/s** | WITHIN | +3.8° (flared) |

The filter intervenes only in the final seconds, converting a hull-budget violation into a
safe touchdown — exactly the safety-filter contract.

---

## 6. Behavioral analysis (why it does what it does)

### The flare is produced by the high-order descent CBF, not the terminal bound
`b` has **relative degree 3**: `b, ḃ, b̈` are control-free; the elevator first appears in `b⃛`.
So `ψ₀=b, ψ₁, ψ₂` are all control-free and the QP can only act through the degree-3 row
`ψ₃ = ψ̇₂ + c₃ψ₂ ≥ 0`. Instrumentation confirms `res_desc = 0` (ψ₃ binding) through the flare,
with `ψ₁, ψ₂ ≈ 1.2–1.5` (slack). `ψ₃ ≥ 0` is *predictive*: through the class-K cascade it keeps
`ψ₂, ψ₁, b ≥ 0` in the future, so the QP must start pitching up (reducing sink) **early**, while
`b` is still far from 0. That early pitch-up *is* the flare. A relative-degree-1 view of `b`
cannot flare (`L_g b = 0`) — the HOCBF is what makes the constraint enforceable.

### Why touchdown sink ≈ 0 rather than v_safe
At the surface `b = v_safe − sink`, so `sink_td = v_safe − b_td`. With linear class-K the
cascade drives `b → 0` only **asymptotically** (rate ~ class-K gains). With the gentle default
`c=[2,2,2]`, `b` is still ≈ v_safe when h reaches 0, leaving `sink ≈ 0`. The CBF *allows* up to
v_safe; the gentle gains simply don't consume the allowance. Stiffening the gains recovers it:

| `c_descent` | touchdown sink |
|---|---|
| [2,2,2] | 0.006 m/s |
| [5,5,5] | 0.338 m/s |
| [12,12,12] | 0.512 m/s |
| [30,30,30] | 0.573 m/s → ≈ v_safe |

**Design implication:** a safety *filter* yields `sink ≤ v_safe` with slack; the exact
touchdown sink is a tuning consequence of the class-K gains, not a property of the barrier. A
*defined* touchdown sink belongs in the nominal trajectory, with the CBF as backstop.

### Should `h ≥ 0` be added as a CBF? — No.
`h` and the descent barrier act on the same channel (`ḣ = V sinγ`, relative degree 3). An
`h ≥ 0` CBF with any class-K `α` gives `sink ≤ α(h)`, and every class-K has `α(0)=0`, so it
forces `sink → 0` at the surface — i.e. it would **forbid touchdown** (permanent hover). The
descent-rate barrier is precisely the fix: it replaces the collapsing `α(h)` with the kinematic
envelope `√(v_safe²+2a_brk h)`, whose value at h=0 is `v_safe > 0`, allowing a controlled
touchdown. The descent barrier *is* the ground-proximity safety condition in the only form
compatible with landing. The worthwhile additions are instead **§3.3 (attitude-coupled impact
safety)** and, if a defined touchdown sink is wanted, a flare reference in the nominal.

---

## 7. Assumptions, simplifications, and limitations

Honest caveats on the validity of the quantitative results:

1. **Placeholder vehicle data.** `aircraft.yaml` is TODO-marked (mass 4 kg, made-up inertia and
   thrust); aero is the test `example.stab`. The results demonstrate the **method**, not a real
   vehicle — specific numbers are not physically calibrated.
2. **Self-consistent plant (directive).** The plant *is* the CBF model (zero model mismatch), so
   the clean barrier-invariance result is partly structural; it validates the math, not
   robustness to model error.
3. **`L_f³b` drift not independently verified.** The tests check the *control* authorities vs the
   spec; the spec never writes `L_f³b`, so a drift error would pass current tests. (Suggested
   follow-up: a test-only finite-difference cross-check as an oracle.)
4. **Frozen local-affine aero.** Autodiff exactness is relative to a C0 piecewise-linear table
   surrogate; 2nd/3rd Lie derivatives zero out true aero curvature, and `example.stab` is
   effectively single-Mach (mild extrapolation).
5. **Initial condition is not a true longitudinal-model trim.** `θ_trim/T_set/X0` are seeded
   from the existing 11-state body-axis trim, a different formulation; V drifts 18→16.8 before
   the flare. (Follow-up: solve a Newton trim on the longitudinal model itself.)
6. **No airspeed control in the nominal** (follows from the constant-thrust directive): airspeed
   is uncontrolled, so the approach decelerates toward the stall floor and leans on the airspeed
   CBF. With a high `Vmin` the energy-limited glider parks in a safe hover instead of landing;
   the default `Vmin≈0.75·V_app` gives margin to touch down.
7. **Soft airspeed barrier** (slack-penalized) — not a hard guarantee. Descent + thrust are hard.
8. **Discrete-time / ZOH** at dt=0.01 approximates the continuous-time CBF guarantee.
9. **No controllability guard** on `L_g L_f²b` (elevator authority → 0 at very low dynamic
   pressure would make the descent barrier unenforceable by elevator; did not occur here).
10. **Idealized scope:** longitudinal-only, perfect full-state feedback, no wind/gust/sensor
    noise, instantaneous elevator (no actuator lag).
11. **§3.3 deferred ⇒ "hull-safe" is incomplete.** `v_safe` is a flat constant and touchdown
    attitude is unconstrained; real slamming load is attitude-dependent.

**Bottom line:** the HOCBF-QP is internally correct and behaves as designed on an idealized,
self-consistent, placeholder-parameter model. Priority follow-ups: independent `L_f³b` check,
a true longitudinal trim, real vehicle data, and §3.3.

---

## 8. How to reproduce

```bash
cmake -S . -B build && cmake --build build -j
ctest --test-dir build --output-on-failure          # 25 tests
./build/lon_autoland_sim                             # default landing -> lon_autoland_log.csv
python3 scripts/plot_lon_results.py  lon_autoland_log.csv lon_trajectory.png 0.6 3.0 13.5 6.0
# CBF on vs off:
sed 's/^cbf_enabled: true/cbf_enabled: false/' data/lon_scenario.yaml > /tmp/off.yaml
./build/lon_autoland_sim data/example.stab data/aircraft.yaml /tmp/off.yaml off.csv
python3 scripts/plot_lon_compare.py lon_autoland_log.csv off.csv lon_compare.png 0.6 3.0 13.5
```
