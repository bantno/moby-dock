# Figures

Closed-loop autoland figures, generated from `autoland_sim` CSV logs
(aero: `data/AHAB_sweep.stab`).

| File | What it shows |
|------|---------------|
| `autoland_response.png` | Standard 6-panel closed-loop response at the nominal approach (V_app=18 m/s): altitude vs glideslope, airspeed, pitch/alpha, sink rate, lateral capture, control surfaces + throttle. CBF is pass-through here (the airspeed barrier never binds). |
| `cbf_on_vs_off.png` | CBF safety filter **ON vs OFF**, approach pushed to the stall margin (V_app=15.5 m/s, barrier V_min=15). ON holds airspeed at 15.0; OFF sags to 14.57. Shows the filter raising throttle to defend the barrier while flying essentially the same glidepath (minimally invasive). |
| `cbf_detail.png` | CBF **internals** (V_app=15.5 m/s): descent rate (vs v_safe), pitch attitude, **nominal vs safe** elevator and throttle, and the two **barrier values** (airspeed `b_V` and descent `b`) with the b=0 safety boundary. The throttle panel shows the filter stepping above the nominal command to defend the airspeed barrier. |

## Regenerate

```bash
# nominal response
./build/autoland_sim                      # writes autoland_log.csv (AHAB_sweep default)
python3 scripts/plot_results.py autoland_log.csv --save figures/autoland_response.png

# CBF on/off study at V_app=15.5 (edit V_app + cbf.enabled in two scenario copies)
./build/autoland_sim data/AHAB_sweep.stab data/aircraft.yaml scenario_on.yaml  cbf_on.csv
./build/autoland_sim data/AHAB_sweep.stab data/aircraft.yaml scenario_off.yaml cbf_off.csv
python3 scripts/cbf_compare.py cbf_on.csv cbf_off.csv --vmin 15 --save figures/cbf_on_vs_off.png

# CBF internals (descent rate, pitch, nominal-vs-safe controls, barrier values)
./build/autoland_sim data/AHAB_sweep.stab data/aircraft.yaml scenario_on.yaml cbf_detail.csv
python3 scripts/plot_cbf_detail.py cbf_detail.csv --vsafe 0.6 --save figures/cbf_detail.png
```

CSV columns added for these plots: `de_nom_deg`, `dT_nom` (the nominal command
before filtering) and `b_airspeed`, `b_descent` (barrier values, b >= 0 = safe).
