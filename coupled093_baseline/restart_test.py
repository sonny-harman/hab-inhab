import sys, json, pickle, time, os
os.environ["OMP_NUM_THREADS"] = "1"
from photochem.extensions.coupled import CoupledPhotochemistryClimate
REPO = "/home/user/hab-inhab"; rid = sys.argv[1]
d = f"{REPO}/saved_inputs/{rid}"; meta = json.load(open(f"{d}/{rid}_run.json")); f = meta["saved_inputs"]["files"]
res = json.load(open(f"{REPO}/coupled093_baseline/{rid}/result.json"))
mech = f"{REPO}/saved_inputs/shared/{'mechanism_haze.yaml' if meta.get('haze', {}).get('enabled') else 'mechanism.yaml'}"
t0 = time.time()
m = CoupledPhotochemistryClimate(mech, f"{d}/{f['settings']}", f"{d}/{f['flux']}"); m.pc.var.atol = 1e-18; m.c.verbose = False
state = pickle.load(open(f"{REPO}/coupled093_baseline/{rid}/model_state.pkl", "rb"))
m.initialize_solver(state, solver_options=dict(max_iterations=30))
while True:
    give_up, conv = m.iterate()
    if give_up or conv: break
a = m.pc.mole_fraction_dict()
print(f"restart from model_state.pkl: converged={conv} after {m.solver.niterations} iteration(s), {time.time()-t0:.0f}s; T_surf {m.c.T_surf:.3f} K (saved {res['T_surf_K']:.3f}); CO {a['CO'][0]:.3e} (saved {res['vmr']['CO']:.3e})")
