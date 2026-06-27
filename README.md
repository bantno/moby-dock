# Autoland development environment — AHAB V-tail flying-boat UAV

A standalone, fast closed-loop simulation for **control-law design** of a V-tail
flying-boat seaplane UAV autoland. This is a design sandbox, **not** onboard code.
It is structured so a Control Barrier Function (CBF) QP safety filter drops in
later without restructuring.

## Build (no system installs)

Every dependency — Eigen, yaml-cpp, Catch2 — is fetched by CMake `FetchContent`,
so a network connection is needed on first configure but nothing must be
installed system-wide.

```bash
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/autoland_sim                 # uses bundled data/, writes autoland_log.csv
python3 scripts/plot_results.py autoland_log.csv
```

`autoland_sim [stab] [aircraft.yaml] [scenario.yaml] [out.csv]` — all optional.

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

## TODO / backlog

Open tasks, data to calibrate, and ideas live in **[`TODO.md`](TODO.md)** (root). The
heavyweight design rationale is in `documentation/water_landing_cbf_design.md`; recent
changes are in `documentation/CHANGELOG.md`.

## Next step: CBF QP

`include/autoland/cbf.hpp` defines the final filter signature and the QP it will
solve; `src/cbf.cpp` marks the OSQP `TODO`. Candidate barriers are listed there
and stubbed in `Sim::buildCandidateBarriers`: minimum airspeed above stall, bank
limit tightening near the surface, and sink rate bounded as a function of height.
