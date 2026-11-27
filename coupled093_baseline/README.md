# Baseline runs re-evaluated with Photochem 0.9.3's coupled mode, at the saved CO2

Step 1 of the plan to compare the whole sweep against `photochem.extensions.coupled.CoupledPhotochemistryClimate`:
each of the 37 baseline runs from `../saved_inputs/` (all steady runs without a CO experiment; the haze-network version at the hot GJ 876
end) is solved once in the new coupled mode **at its saved CO2 pressure**. No CO2 search is done here, so this shows how far the sweep's
temperature is from the new mode's at the same CO2, not the CO2 needed for 288 K.

## Method

Same setup as `../photochem_v093_spotcheck/` (script: `coupled093_point.py <run_id> <output_root>`, run in a Photochem 0.9.3 environment):

- Settings, stellar flux and mechanism from the saved run; CO2 lower boundary at the saved pressure; each case starts from scratch
  (`initial_guess`, with CH4 and H2 supplied at the sweep's surface pressures, tropopause pressure from the saved run), then the default
  `FixedPointSolver` (temperature tolerance 1e-3, pressure tolerance 1e-3 dex, up to 30 iterations).
- Kzz as in the sweep, integrator `atol = 1e-18`, 50 climate layers, climate top and photochemical top at 1 dyn/cm2.
- The climate model sees all photochemical gases and their vertical profiles (the sweep's climate model had only H2O, CO2, N2, CH4, H2).

## Results (`map_table.csv`, `map_dT.png`)

T(0.9.3 coupled) minus T(0.9.1 sweep) at the same CO2, in K. 34 of 37 converged; 3 failed (below).

| Regime | Instellation range (W/m2) | Difference (K) | Notes |
|---|---|---|---|
| Sun, warm to hot | 575-1400 | +0.4 to +1.7 | all converged in 3-6 iterations; CO within about 15% |
| GJ 876, mid-range | 437.5-987.5 | -1.2 to +2.3 | CO within about 5% of the sweep up to about 800 W/m2 |
| GJ 876, hot, no haze network | 1056.25-1067.71 | +2.4 to +3.0 | CO about 1.6-1.7 times lower |
| GJ 876, hot, haze network | 1079.17 | +4.0 | CO 1.7 times lower |
| GJ 876, hot, haze network | 1102.08-1262.5 | **+12.4 to +14.8** | CO 2 to 8 times lower; see the fuller-climate folder |
| Sun, cold (CO2 capped at 10 bar or max-greenhouse) | 300-437.5 | +1.4 to +3.1 | surface pressure differs strongly (see below) |
| Sun, 506.25 | 506.25 | **+17.1** | the sweep's maximum-greenhouse case: 268.7 K in the sweep, 285.8 K here at the same 5.3 bar CO2 |
| GJ 876, 300 W/m2 | 300 | +13.1 | cold CO2-capped case |

Failures (initial climate solve did not converge, `ClimaException`): GJ 876 at **1331.25 and 1400 W/m2** (haze network, CO2 at the 1e-5 bar floor)
and at **368.75 W/m2** (cold, 10 bar CO2). Their `result.json` records the error; there is no atmosphere or state for them.

Other points from the table:

- At the cold end the total surface pressure differs a lot: Sun 368.75 W/m2, 8.58 bar in the sweep and 3.33 bar here; Sun 437.5 W/m2, 11.65 and 8.30 bar;
  GJ 876 300 W/m2, 2.89 and 1.69 bar. This looks like a different treatment of CO2 condensation; I did not check which is closer to the saturation limit.
- The haze-network GJ 876 results at 1102.08 (303.19 K), 1125 (301.80 K), 1193.75 (300.67 K) and 1262.5 (302.53 K) reproduce the first evaluations of the
  fuller-climate CO2 search in `../fuller_climate_haze/`.
- At the cold end the 0.9.1 runs were not at 288 K (CO2 capped, or maximum-greenhouse CO2), so differences there mix the method difference with a very different climate state.

## What is saved for each point (`<run_id>/`)

- `result.json`: status, iteration count, solver residual history, surface composition (CO2, CO, O2, H2O, CH4, H2, N2), minor gases above 1e-7, the 0.9.1 values for comparison, and settings used.
- `photochem_atmosphere.txt`: the final photochemical atmosphere in Photochem's format.
- `climate_profile.txt`: the final radiative-convective profile (pressure, temperature, altitude).
- `model_state.pkl`: pickled output of `CoupledPhotochemistryClimate.model_state()` for converged points, to warm-start a new coupled solve. It is only valid with Photochem 0.9.3 or a compatible version, and a pickle should only be loaded from a source you trust (these were written here).

```python
import pickle
from photochem.extensions.coupled import CoupledPhotochemistryClimate
# build the model exactly as in coupled093_point.py (same mechanism, settings, flux), then:
state = pickle.load(open("<run_id>/model_state.pkl", "rb"))
m.initialize_solver(state, solver_options=dict(max_iterations=30))
```

Checked on one point (Sun, 1400 W/m2): restarting from `model_state.pkl` converged after 1 iteration (about 5 minutes including model construction) to T_surf 288.266 K (saved 288.265 K) and the same surface CO (1.629e-6). `restart_test.py` is that check.

The run IDs match the folders in `../saved_inputs/`, which hold the settings and flux for each case.

## Files

`coupled093_point.py` (the run script), `restart_test.py`, `summarize_map.py` (builds the table and figure), `map_table.csv`, `map_dT.png`, and one folder per run.
