# Autoland development environment — AHAB V-tail flying-boat UAV

A standalone, fast closed-loop simulation for **control-law design** of a V-tail
flying-boat seaplane UAV autoland. This is a design sandbox, **not** onboard code.
It is structured so a Control Barrier Function (CBF) QP safety filter drops in
later without restructuring.

## Build & test (no system installs)

Every dependency — Eigen, yaml-cpp, Catch2, OSQP — is fetched by CMake `FetchContent`,
so a network connection is needed on first configure but nothing must be installed
system-wide.

```bash
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Run the water-landing CBF sim + plot the results

The primary target is the **augmented-longitudinal CBF-QP water-landing sim**
(`lon_autoland_sim`). Its safety filter enforces the *recovery* barrier set:
hydrodynamic **impact-load** (the only hard safety row) plus soft **stall/AoA**,
**nose-up attitude**, and a **total-energy** ceiling, with the thrust actuator/validity
guards. Barriers and gains live under the `cbf:` block of `data/lon_scenario.yaml`; the
AoA guard reproduces the pilot low-altitude stall recovery (pitch down + full elevator).

### Option A — run-and-plot harness (recommended)

`scripts/harness.py` rebuilds the binary if the sources are stale, applies inline
scenario overrides, runs the sim, prints the touchdown summary + a barrier-minimum
table, and plots the trace — with **every plot annotation (`V_td_max`, `alpha_stall`,
`theta_min`, `Tmax`, …) read from the resolved scenario**, so the figure can never drift
out of sync with the run.

```bash
# single run (defaults: data/lon_scenario.yaml on data/AHAB_combined.stab)
python3 scripts/harness.py run

# override any scenario field (dotted key; value parsed as YAML)
python3 scripts/harness.py run --name slow_td --set cbf.V_td_max=11 --set stall.enabled=true
python3 scripts/harness.py run --set cbf.enabled=false -o figures/nominal_only.png

# overlay N labeled cases
python3 scripts/harness.py compare \
    --case "CBF on:cbf.enabled=true" \
    --case "CBF off:cbf.enabled=false"
python3 scripts/harness.py compare --name egains \
    --case "e2:cbf.c_energy=[2,2,2]" \
    --case "e5:cbf.c_energy=[5,5,5]"
```

Each run writes `runs/<name>.{csv,resolved.yaml,log,png}` (git-ignored). The barrier
table marks the **hard** rows (impact + thrust) `OK`/`VIOLATION` and reports the **soft**
rows' minima (stall / nose-up / energy — they may dip transiently). The harness also
surfaces NaN barrier values and the run-to-run jitter the QP shows in feasibility-recovery.

### Option B — run the binary directly, then plot

```bash
./build/lon_autoland_sim data/AHAB_combined.stab data/aircraft.yaml data/lon_scenario.yaml lon_log.csv

# standard 6-panel trace:   args = [csv] [out.png] [V_td_max] [g_eff] [Tmax]
python3 scripts/plot_lon_results.py lon_log.csv figures/lon_landing.png 14 16 50
# last N s before touchdown: args = ... [V_td_max] [g_eff] [Tmax] [window_s]
python3 scripts/plot_lon_zoom.py    lon_log.csv figures/lon_zoom.png    14 16 50 5
# impact-load barrier detail
python3 scripts/plot_lon_impact.py  lon_log.csv figures/lon_impact.png
```

`lon_autoland_sim [stab] [aircraft.yaml] [lon_scenario.yaml] [out.csv]` — all optional
(defaults to the bundled `data/` files).

> A separate legacy body-axis sim (`autoland_sim`, `src/sim.cpp`) exists as an older
> design sandbox and is **not** part of the CBF workflow.

## Architecture

```
aero_table   VSPAero .stab parser + (alpha,beta,Mach) trilinear interpolation
dynamics     ONE nonlinear 6-DOF EOM  Dynamics::xdot(x,u)  <-- single source of truth
trim         Newton solve on xdot for steady descent
linear_model central-difference linearization of xdot; split lon/lat sub-models
mixing       virtual [de,da,dr] -> physical control groups (config-driven)
controller   cascaded PID, frontside technique (+ flare, decrab)
cbf          CBF safety-filter interface; pass-through stub (OSQP later)
sim          load -> trim -> linearize -> closed loop to touchdown -> CSV
config       aircraft.yaml / scenario.yaml loaders
```

The trim solver and the linear model **both derive from `Dynamics::xdot`** — the
nonlinear EOM is the foundation. The closed loop currently runs against the
linear plant; switching `Sim`'s `plant` lambda to `Dynamics::xdot` runs the full
nonlinear plant with no other change.

## Key modeling facts (verified against `data/example.stab`)

- **Body axes throughout.** Forces/moments are assembled from body-axis
  coefficients `CFx,CFy,CFz,Cl=CMx,Cm=CMy,Cn=CMz`, never from wind-axis CL/CD.
  Full sign-convention block in `include/autoland/linear_model.hpp`.
- **Sweep grid:** 5×5 `(alpha, beta)` at a single Mach 0.059. Off-grid queries
  clamp and warn — never silent extrapolation.
- **Control groups** are read from the file header (not hardcoded):
  `ConGrp_1=Ailerons`, `ConGrp_2=Elevator`, `ConGrp_3=Rudder` — the V-tail
  ruddervator mixing is already baked into the OpenVSP groups, so the default
  mixing matrix is identity-style.
- **Rate derivatives** use VSPAero's nondimensional convention
  `phat=p·b/2V, qhat=q·c/2V, rhat=r·b/2V`; angle derivatives are per-radian.
  Both verified numerically against the perturbation rows.
- **Units:** model geometry is SI metres/seconds (Mach 0.059 = Vinf 20 / a forces
  m/s). The file's `Rho_=0.002377` is a stale imperial default and is ignored;
  air density comes from `aircraft.yaml` (default 1.225 kg/m³).
- **Stall:** the inviscid `.stab` cannot stall (lift is linear to ±20°). An optional
  NACA 4414 viscous-stall overlay (NeuralFoil-derived) splices a real post-stall lift
  drop / drag rise / pitch break onto the deck — OFF by default, enabled via
  `aircraft.yaml`'s `stall:` block. See `documentation/stall_model_spec.md`.

## TODO / backlog

Open tasks, data to calibrate, and ideas live in **[`TODO.md`](TODO.md)** (root). The
heavyweight design rationale is in `documentation/water_landing_cbf_design.md`; recent
changes are in `documentation/CHANGELOG.md`.

## CBF safety filter

The **longitudinal CBF-QP is implemented** (OSQP-backed):
`include/autoland/{hocbf,lon_cbf_filter,impact_barrier}.hpp`, `src/lon_cbf_filter.cpp`,
with exact Lie derivatives from `include/autoland/lie_taylor.hpp`. Design rationale is in
`documentation/water_landing_cbf_design.md` and the derivations (barrier definitions,
relative degrees, QP) in `documentation/water_landing_cbf_math.md`.

The older body-axis filter (`include/autoland/cbf.hpp`, `src/cbf.cpp`) remains a
pass-through stub for the legacy `autoland_sim` path only.
