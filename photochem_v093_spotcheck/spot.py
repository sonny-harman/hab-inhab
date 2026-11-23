"""Spot-check one saved run against Photochem 0.9.3's CoupledPhotochemistryClimate.
usage: python spot.py <run_id>   (run in the ph93 env)"""
import sys, os, json, time, numpy as np
os.environ["OMP_NUM_THREADS"] = "1"
import photochem
from photochem.extensions.coupled import CoupledPhotochemistryClimate
REPO = "/home/user/hab-inhab"; rid = sys.argv[1]
d = f"{REPO}/saved_inputs/{rid}"; meta = json.load(open(f"{d}/{rid}_run.json")); f = meta["saved_inputs"]["files"]
haze = meta.get("haze", {}).get("enabled", False)
mech = f"{REPO}/saved_inputs/shared/{'mechanism_haze.yaml' if haze else 'mechanism.yaml'}"
t0 = time.time()
m = CoupledPhotochemistryClimate(mech, f"{d}/{f['settings']}", f"{d}/{f['flux']}", climate_condensates=("H2O", "CO2"),
                                 climate_number_of_layers=50, climate_top_pressure=1.0)
m.pc.var.atol = 1.0e-18        # same integrator tolerance as the 0.9.1 sweep (default 1e-23 is very slow for cold CO2-rich cases)
def kzz(P_bar): return np.where(P_bar > 0.1, 1.0e5, np.minimum(1.0e5*(0.1/P_bar)**0.5, 1.0e7))   # same Kzz as the sweep
P = np.logspace(np.log10(2.0e7), 0.0, 80)                 # dyn/cm^2, decreasing
trop_p = meta["climate"]["P_trop_bar"] * 1e6          # dyn/cm^2: the tropopause pressure the sweep used for this run
try:
    guess = m.initial_guess(kzz(P/1e6), P, target_pressure=1.0, tropopause_pressure=trop_p)
except Exception as e:
    print(f"[{rid}] FAILED in initial_guess: {type(e).__name__}: {e}", flush=True)
    json.dump({"run_id": rid, "failed": f"{type(e).__name__}: {e}", "stage": "initial_guess", "T_surf_K_old": meta["T_surf_K"],
               "pCO2_bar": meta["pCO2_bar"], "star": meta["star"], "S_Wm2": meta["S_Wm2"]}, open(f"out/{rid}.json", "w"), indent=1)
    sys.exit(0)
print(f"[{rid}] initial guess built, T_surf={m.c.T_surf:.2f} K", flush=True)
m.initialize_solver(guess, solver_options=dict(max_iterations=30, verbose=True))
while True:
    give_up, conv = m.iterate()
    h = m.solver.solver_history[-1]
    print(f"[{rid}] it {m.solver.niterations}: dT/T={h['temperature_residual']:.2e} dlogP={h['log10pressure_residual']:.2e} chem_conv={h['chemistry_converged']} T_surf={m.c.T_surf:.2f} ({time.time()-t0:.0f}s)", flush=True)
    if give_up or conv: break
a = m.pc.mole_fraction_dict()
new = {k: float(a[k][0]) for k in ["CO2", "CO", "O2", "H2O", "CH4", "H2", "N2"]}
out = {"run_id": rid, "photochem_version": photochem.__version__, "fixed_point_converged": bool(conv), "give_up": bool(give_up),
       "chemistry_converged": m.solver.chemistry_converged, "iterations": m.solver.niterations, "wall_s": time.time() - t0,
       "T_surf_K_new": float(m.c.T_surf), "P_total_bar_new": float(m.pc.wrk.surface_pressure), "vmr_new": new,
       "T_surf_K_old": meta["T_surf_K"], "P_total_bar_old": meta["P_total_bar"], "vmr_old": meta["vmr"],
       "pCO2_bar": meta["pCO2_bar"], "star": meta["star"], "S_Wm2": meta["S_Wm2"], "haze": haze,
       "co_experiment": meta.get("co_experiment"), "residual_history": m.solver.solver_history}
json.dump(out, open(f"out/{rid}.json", "w"), indent=1, default=str)
print(f"[{rid}] DONE converged={conv} T_new={out['T_surf_K_new']:.2f} T_old={out['T_surf_K_old']:.2f} ({out['wall_s']:.0f}s)", flush=True)
