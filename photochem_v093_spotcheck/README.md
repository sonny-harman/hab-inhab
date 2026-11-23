# Spot-check of the 0.9.1 sweep against Photochem 0.9.3's coupled mode

Photochem 0.9.3 adds `photochem.extensions.coupled.CoupledPhotochemistryClimate`, a dedicated coupled
climate-photochemistry solver. This folder checks seven saved runs from the 0.9.1 sweep against it.

## What was compared

The sweep (Photochem 0.9.1) solved CO2 for a 288 K surface by alternating `AdiabatClimate` and `EvoAtmosphere`.
The new mode does not solve for CO2: it takes a fixed CO2 pressure and solves for temperature and surface pressure.
So each check **fixes the saved CO2 pressure and the saved boundary conditions** (the saved `*_photochem_settings.yaml`,
mechanism and stellar flux) and asks whether the new mode returns the sweep's temperature and surface composition.

- Each case starts from scratch in the new mode (not from the saved state), so it is an independent check.
- Same Kzz profile as the sweep, integrator `atol = 1e-18` (as in the sweep), tropopause pressure from the saved run,
  climate top at 1 dyn/cm2, 50 climate layers, solver defaults (`temperature_tol = 1e-3`, 30 iterations max).
- The new mode's climate model includes **all 92 photochemical gases**; the sweep's climate model used only
  H2O, CO2, N2, CH4 and H2.

## Results (`spotcheck_table.csv`)

| Case | pCO2 (bar) | T 0.9.1 (K) | T 0.9.3 coupled (K) | difference (K) | converged | iterations |
|---|---|---|---|---|---|---|
| GJ 876, 643.75 W/m2 | 0.98 | 288.11 | 288.26 | +0.15 | yes | 3 |
| Sun, 1262.5 W/m2 | 1.09e-3 | 287.70 | 288.41 | +0.71 | yes | 4 |
| Sun, 850 W/m2 | 0.306 | 287.45 | 289.18 | +1.73 | yes | 6 |
| Sun, 850 W/m2, CO vdep 1.2e-4, flux 0 | 0.372 | 288.11 | 288.87 | +0.75 | yes | 4 |
| Sun, 368.75 W/m2 (cold, CO2 capped) | 10 | 203.95 | 205.41 | +1.46 | yes | 11 |
| GJ 876, 1125 W/m2 (haze network) | 8.94e-4 | 288.11 | 301.69 | **+13.58** | yes | 7 |
| GJ 876, 1400 W/m2 (haze network, CO2 floor) | 1e-5 | 305.68 | -- | -- | **failed** | -- |

Surface composition in the converged cases (old -> new, mixing ratios):

- CH4: +33% at Sun 850 W/m2 (1.98e-4 -> 2.63e-4), +19% at GJ 876 643.75 W/m2 (3.64e-4 -> 4.34e-4), and within 4% at the other four converged cases.
- CO: within 12% in the warm cases (Sun 850: 1.52e-2 -> 1.71e-2) and a factor 2.7 higher in the cold Sun case (5.3e-3 -> 1.4e-2); at GJ 876 1125 W/m2 it is a factor 3 lower (1.85e-2 -> 6.2e-3), where the surface is 13 K warmer.
- H2O follows the surface temperature (it is 1.9x higher at GJ 876 1125 W/m2) and is within 25% in the other warm cases.

## Findings

1. **At five of the six runnable cases the two methods agree to within 0.15-1.7 K.** The new mode reproduces
   the sweep's surface temperature and its CH4, H2 and total pressure to that level, including a case with the fast CO sink.
2. **GJ 876, 1125 W/m2 (haze network) differs by 13.6 K.** This is not the climate solver: with the same 5 gases at the same
   surface composition, the 0.9.3 climate code returns 288.11 K, identical to the 0.9.1 result
   (`climate_only.py`). It is also not CH4, H2 or CO: the new-mode composition of those gases gives 288.8-289.9 K in the
   sweep-style climate (`decompose.py`). Giving the 0.9.3 climate model the **full 0.9.1 composition (all gases and their vertical profiles)**
   gives 300.65 K (`full_comp.py`). So what warms this case is the gases beyond the sweep's five (and/or their
   vertical structure). The 0.9.1 surface atmosphere here holds ppm levels of C3H4, C2H2, C2H6, C4H2 and C2H4, which are
   the obvious candidates, but I did not isolate which gases are responsible.
3. **The new mode failed at GJ 876, 1400 W/m2 (305 K, CO2 floor).** `initial_guess` raised
   `ClimaException: Climate did not converge`, with and without CH4/H2 in the initial guess. Giving the climate model the full 0.9.1 composition
   for this point also did not converge (T = 555 K, not converged). The sweep's own approximate-then-RCE solve converged at 305.7 K.
4. **Cold case (Sun, 368.75 W/m2):** the temperatures agree to 1.5 K, but the surface pressure does not: 3.35 bar in the new mode
   against 8.58 bar in the sweep. The new mode's surface CO2 partial pressure is 2.2 bar, against 7.5 bar in the sweep;
   both are below the nominal 10 bar boundary value, which suggests that condensation limits CO2 in both, by very different amounts.
   I have not checked the saturation vapour pressure explicitly, so which treatment is closer to the saturation limit is not established here.
   (N2 is 0.31 vs 0.12 and CO 1.4e-2 vs 5.3e-3 as a result.)

## What this implies for the earlier results

- **The haze-chemistry GJ 876 points (1079-1400 W/m2) are the ones most likely to be affected.** At 1125 W/m2 the more complete climate gives roughly
  12-13 K warmer at the saved CO2, so a consistent solution would need less CO2 and the range where 288 K can be held would shift. This is
  established for that one case only; the other haze points were not checked.
- **For the CO deposition experiments**, the CO flux/sink changes CH4 and (in the hot GJ 876 case) hydrocarbons, so the temperature
  response to CO deposition velocity may be incomplete where hydrocarbons matter. At the Sun cases checked it is not (differences <= 1.7 K).
- **Cold, CO2-capped cases** (below about 575 W/m2 for the Sun) should be treated as indicative only, as the sweep README already says;
  this check shows the surface pressure and CO2 partial pressure there depend strongly on the method.

## Caveats

- Seven runs out of 120 were checked, and not randomly: they were chosen to span warm/cold, low/high CO2, the haze network and the fast CO sink.
- Differences between 0.9.1 and 0.9.3 other than the coupled mode (mechanism data, `conv_longdy`, solver defaults) were not isolated.
- The "full composition" climate tests use the 0.9.1 atmosphere files; they test the climate model only, not the chemistry.

## Files

- `spot.py`: run one saved case in the new coupled mode (`python spot.py <run_id>` in a Photochem 0.9.3 environment).
- `climate_only.py`, `decompose.py`, `full_comp.py`: the climate-only tests behind finding 2; their outputs are in `decompose.txt` and `full_comp.txt`.
- `results/*.json`: per-case outputs, including the solver's residual history.
- `spotcheck_table.csv`: the table above with all compared quantities.
