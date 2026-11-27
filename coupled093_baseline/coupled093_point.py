"""Run one saved 0.9.1 run through Photochem 0.9.3's CoupledPhotochemistryClimate at the saved CO2 pressure.
usage: python coupled093_point.py <run_id> <output_root>        (Photochem 0.9.3 environment)
Writes <output_root>/<run_id>/{result.json, photochem_atmosphere.txt, climate_profile.txt, model_state.pkl}."""
import sys, os, json, time, pickle, numpy as np
os.environ["OMP_NUM_THREADS"] = "1"
import photochem
from photochem.extensions.coupled import CoupledPhotochemistryClimate
REPO = "/home/user/hab-inhab"; rid = sys.argv[1]; root = sys.argv[2]
d = f"{REPO}/saved_inputs/{rid}"; meta = json.load(open(f"{d}/{rid}_run.json")); f = meta["saved_inputs"]["files"]
haze = meta.get("haze", {}).get("enabled", False)
mech = f"{REPO}/saved_inputs/shared/{'mechanism_haze.yaml' if haze else 'mechanism.yaml'}"
out = f"{root}/{rid}"; os.makedirs(out, exist_ok=True)
t0 = time.time(); log = lambda s: print(f"[{rid}] {s} ({time.time()-t0:.0f}s)", flush=True)
base = dict(run_id=rid, star=meta["star"], S_Wm2=meta["S_Wm2"], haze=haze, pCO2_bar=meta["pCO2_bar"], photochem_version=photochem.__version__,
            old=dict(T_surf_K=meta["T_surf_K"], P_total_bar=meta["P_total_bar"], vmr=meta["vmr"], photochem_steady=meta["photochem_steady"], target_status=meta["target_status"]),
            settings=dict(saved_inputs_run=rid, tropopause_pressure_dyn_cm2=meta["climate"]["P_trop_bar"] * 1e6, atol=1e-18, target_pressure_dyn_cm2=1.0,
                          climate_layers=50, climate_top_pressure_dyn_cm2=1.0, solver="FixedPointSolver defaults, max_iterations=30"))
def kzz(P_bar): return np.where(P_bar > 0.1, 1.0e5, np.minimum(1.0e5*(0.1/P_bar)**0.5, 1.0e7))
def finish(rec):
    rec["wall_s"] = time.time() - t0; json.dump(rec, open(f"{out}/result.json", "w"), indent=1, default=str)
try:
    m = CoupledPhotochemistryClimate(mech, f"{d}/{f['settings']}", f"{d}/{f['flux']}")
    m.pc.var.atol = 1.0e-18; m.c.verbose = False
    sp = {k: v for k, v in meta["climate"]["surface_partial_pressures_dyn_cm2"].items() if k in ("CO2", "N2", "H2O", "CH4", "H2")}
    P = np.logspace(np.log10(2.0e7), 0.0, 80)
    guess = m.initial_guess(kzz(P/1e6), P, target_pressure=1.0, tropopause_pressure=base["settings"]["tropopause_pressure_dyn_cm2"], surface_partial_pressures=sp)
except Exception as e:
    log(f"FAILED in initial_guess: {type(e).__name__}: {e}"); finish({**base, "status": "failed", "stage": "initial_guess", "error": f"{type(e).__name__}: {e}"}); sys.exit(0)
log(f"initial guess built, T_surf={m.c.T_surf:.2f} K")
try:
    m.initialize_solver(guess, solver_options=dict(max_iterations=30, verbose=False))
    while True:
        give_up, conv = m.iterate()
        if give_up or conv: break
except Exception as e:
    log(f"FAILED in coupling iterations: {type(e).__name__}: {e}")
    finish({**base, "status": "failed", "stage": "coupling", "error": f"{type(e).__name__}: {e}", "iterations_completed": m.solver.niterations,
            "solver_history": m.solver.solver_history}); sys.exit(0)
a = m.pc.mole_fraction_dict()
res = {**base, "status": "converged" if conv else "gave_up", "fixed_point_converged": bool(conv), "chemistry_converged": m.solver.chemistry_converged,
       "iterations": m.solver.niterations, "T_surf_K": float(m.c.T_surf), "P_total_bar": float(m.pc.wrk.surface_pressure),
       "vmr": {k: float(a[k][0]) for k in ("CO2", "CO", "O2", "H2O", "CH4", "H2", "N2")},
       "minor_gases_gt_1e7": {k: float(a[k][0]) for k in m.pc.dat.species_names[m.pc.dat.np:m.pc.dat.nq] if a[k][0] > 1e-7 and k not in ("N2", "CO2", "H2O", "H2", "CO", "CH4", "He")},
       "dT_K_vs_0.9.1": float(m.c.T_surf) - meta["T_surf_K"], "solver_history": m.solver.solver_history}
m.pc.out2atmosphere_txt(f"{out}/photochem_atmosphere.txt", number_of_decimals=10, overwrite=True)
c = m.c
np.savetxt(f"{out}/climate_profile.txt", np.c_[c.P / 1e6, c.T, c.z / 1e5], fmt="%.10e",
           header=f"columns: pressure_bar temperature_K altitude_km (top of atmosphere last)\nT_surf_K {c.T_surf:.10f}\nP_surf_bar {c.P_surf/1e6:.10e}")
if conv:
    try: pickle.dump(m.model_state(), open(f"{out}/model_state.pkl", "wb"))        # warm start for CoupledPhotochemistryClimate.initialize_solver
    except Exception as e: res["model_state_error"] = f"{type(e).__name__}: {e}"
finish(res)
log(f"DONE {res['status']}: T_surf={res['T_surf_K']:.2f} K (0.9.1: {meta['T_surf_K']:.2f}), {res['iterations']} iterations")
