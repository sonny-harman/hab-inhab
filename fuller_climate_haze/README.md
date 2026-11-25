# Hot GJ 876 haze points rerun with a fuller climate model (Photochem 0.9.3)

The seven GJ 876 haze-chemistry points (1079-1400 W/m2) of the main sweep were solved for CO2 with a climate model that
contained only H2O, CO2, N2, CH4 and H2. The spot-check in `../photochem_v093_spotcheck/` showed that at 1125 W/m2 the full
photochemical composition (including ppm-level hydrocarbons) warms the surface by about 13 K at the same CO2. Here the CO2 needed
for a 288 K surface is solved again with Photochem 0.9.3's `CoupledPhotochemistryClimate`, in which the radiative-convective model
sees **all photochemical gases and their vertical profiles**.

## Method

For each instellation (`solve_fuller.py <S>`, run in a Photochem 0.9.3 environment):

1. Take the settings, stellar flux, haze mechanism and starting CO2 from the sweep's haze run at that instellation.
2. Build the coupled state (`initial_guess`, with CH4 and H2 supplied) and converge the coupled solver
   (temperature and surface pressure) at the starting CO2. Solver defaults: `temperature_tol = 1e-3`, `log10pressure_tol = 1e-3`, 30 iterations.
3. Search for the CO2 giving T_surf = 288 K (secant in log CO2, tolerance 0.5 K, bounds 1e-5 to 10 bar, as in the sweep). After each CO2 change
   the coupled solver restarts from the previous converged state.
4. If a CO2 evaluation failed to converge, bisect back toward the last converged CO2.

Same Kzz profile and integrator `atol = 1e-18` as the sweep. The tropopause pressure used for the photochemical water cold trap is held at
the sweep run's value (it does not follow the changing CO2).

## Results (`comparison.csv`, `comparison.png`)

| S (W/m2) | CO2 sweep (bar) | CO2 fuller climate (bar) | T fuller (K) | CO sweep -> fuller | CH4 sweep -> fuller |
|---|---|---|---|---|---|
| 1079.17 | 2.5e-3 | 2.8e-4 | 288.08 | 2.8e-2 -> 3.4e-3 | 3.5e-5 -> 1.1e-5 |
| 1102.08 | 1.5e-3 | 1.7e-4 | 288.12 | 2.5e-2 -> 1.8e-3 | 3.0e-5 -> 8.2e-6 |
| 1125.00 | 8.9e-4 | 9.1e-5 | 288.18 | 1.8e-2 -> 7.3e-4 | 2.2e-5 -> 5.9e-6 |
| 1193.75 | 1.7e-4 | 1e-5 (floor) | 289.57 (288 K not reached) | 6.1e-3 -> 9.9e-5 | 1.0e-5 -> 4.2e-6 |
| 1262.50 | 2.2e-5 | 1e-5 (floor) | 299.83 (288 K not reached) | 2.4e-4 -> 8.7e-5 | 3.7e-6 -> 5.4e-6 |
| 1331.25 | 1e-5 (floor) | none found (see below) | 306.7 at 1e-4 bar; 326.4 at 5e-5 bar | -- | -- |
| 1400.00 | 1e-5 (floor) | -- | -- | -- | -- |

- **The CO2 needed for 288 K is about 9-10 times lower than in the sweep** at 1079-1125 W/m2.
- **The 288 K boundary moves to lower instellation.** The sweep held 288 K up to 1262.5 W/m2 (2.2e-5 bar); with the fuller climate 288 K
  is held at 1125 W/m2 but not at 1193.75 W/m2, where T_surf is 1.6 K above target even at the 1e-5 bar floor. So the boundary lies between those two points.
- **CO and CH4 are lower** (CO by factors of 8 to 25 at the three points that hold 288 K, CH4 by factors of 3 to 4) because the atmosphere now
  needs far less CO2, which feeds the CO source. This is an inference from the trend; I did not run a separate test of it.
- The temperatures at the saved sweep CO2 were 301.8 K (1125), 300.7 K (1193.75), 303.2 K (1102.08) and 302.5 K (1262.5), i.e. 12-15 K warmer than the sweep's 288 K.

## Points without a result

- **1331.25 W/m2: 288 K not reached at any CO2 that could be evaluated.** Fresh starts at three CO2 values all converged:

  | pCO2 (bar) | T_surf (K) | CH4 | H2O | other gases above 1e-7 |
  |---|---|---|---|---|
  | 2.0e-4 | 313.1 | 6.1e-6 | 3.2e-2 | none |
  | 1.0e-4 | 306.7 | 4.8e-6 | 2.2e-2 | none |
  | 5.0e-5 | 326.4 | 5.2e-5 | 6.2e-2 | C2H6 6.6e-6 |

  Temperature is **not monotonic in CO2**: the 5e-5 bar state (reproduced by a continuation run from 1e-4 bar: 326.3 K) is a methane-rich state, with CH4 about 10 times
  higher and C2H6 present, so I take it to be a hydrocarbon-greenhouse state, but I did not test that explanation. The 1e-5 bar floor could not be evaluated: the new mode's initial
  climate solve fails there, and earlier attempts to step down from 1e-4 bar either failed in the coupling iterations (`Climate did not converge during _fixed_point_map`), were
  lost to a container restart, or stalled with the integrator restarting from tiny time steps. Because of the non-monotonic behaviour I can no longer say whether 288 K is reachable
  at lower CO2; the evaluated points are all 307 K or warmer.
- **1400 W/m2: no solution.** The new mode's climate solve failed at every starting CO2 tried (1e-5, 1e-4, 1e-3, 1e-2, 0.1 and 1 bar). I did not establish whether
  this reflects a runaway greenhouse or a solver limitation. The sweep's own solve had converged at 305.7 K (1e-5 bar).

## Caveats

- Only the haze GJ 876 points were rerun. The other sweep points (including the Sun and the cooler GJ 876 points) were not.
- Haze particles act in the photochemistry only; their optical effects are still not in the climate model.
- The tropopause pressure is fixed at the sweep value; surface-pressure and temperature profiles come from the coupled solver.
- The two 288 K-unreachable converged points sit at the CO2 floor, so their composition is for that floor value.
- The CO2 boundary is not re-derived for the CO deposition experiments here; those results were produced with the five-gas climate.

## Files

- `GJ876_S<S>_FullClimate093/`: final photochemical atmosphere (`photochem_atmosphere.txt`), climate profile (`climate_profile.txt`) and
  `result.json` (full outer-search history, composition, convergence flags) for each completed point.
- `GJ876_S1331.25_FullClimate093_probe*/`: the three single-CO2 fresh evaluations at 1331.25 W/m2 (suffix = log10 of the starting CO2 in bar).
- `solve_fuller.py`: the solver (environment variables `START_LOG10`, `MAX_STEP_DEX`, `MAX_OUTER`, `TAG` were added for the 1331.25 probes); `summarize.py`: builds the table and figure.
- `comparison.csv`, `comparison.png`.
