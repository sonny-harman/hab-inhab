"""Solve CO2 for a 288 K surface with Photochem 0.9.3's coupled mode (climate sees all photochemical gases).
usage: python solve_fuller.py <S_Wm2>      (GJ 876, haze network; run in the ph93 env)"""
import sys, os, json, glob, time, numpy as np
os.environ["OMP_NUM_THREADS"] = "1"
import photochem
from photochem.extensions.coupled import CoupledPhotochemistryClimate
REPO = "/home/user/hab-inhab"; S = float(sys.argv[1]); T_TARGET = 288.0; TOL_K = 0.5
PMIN, PMAX = 1e-5, 10.0                                  # bar, same bounds as the sweep
# the baseline haze run at this S supplies settings, flux, mechanism and the starting CO2
base = [json.load(open(p)) for p in glob.glob(f"{REPO}/saved_inputs/GJ876_S{S:07.2f}_*_haze/*_run.json")
        if "COvdep" not in p and "hazeR" not in p]
assert len(base) == 1, base
meta = base[0]; rid = meta["run_id"]; d = f"{REPO}/saved_inputs/{rid}"; f = meta["saved_inputs"]["files"]
mech = f"{REPO}/saved_inputs/shared/mechanism_haze.yaml"
trop_p = meta["climate"]["P_trop_bar"] * 1e6
sp = meta["climate"]["surface_partial_pressures_dyn_cm2"]
def kzz(P_bar): return np.where(P_bar > 0.1, 1.0e5, np.minimum(1.0e5*(0.1/P_bar)**0.5, 1.0e7))
P = np.logspace(np.log10(2.0e7), 0.0, 80)
SOLVER = dict(max_iterations=30, verbose=False)
t0 = time.time(); log = lambda s: print(f"[S={S}] {s} ({time.time()-t0:.0f}s)", flush=True)

def new_model(): 
    m = CoupledPhotochemistryClimate(mech, f"{d}/{f['settings']}", f"{d}/{f['flux']}")
    m.pc.var.atol = 1.0e-18; m.c.verbose = False
    return m
def set_co2(m, pbar): m.pc.set_lower_bc("CO2", "press", press=pbar * 1e6)
def build(pbar):
    "initial climate-consistent state at this CO2 (None if the climate solve fails)"
    m = new_model(); set_co2(m, pbar)
    spp = {k: v for k, v in sp.items() if k in ("N2", "H2O", "CH4", "H2")}; spp["CO2"] = pbar * 1e6
    try:
        g = m.initial_guess(kzz(P/1e6), P, target_pressure=1.0, tropopause_pressure=trop_p, surface_partial_pressures=spp)
    except Exception as e:
        log(f"initial_guess failed at pCO2={pbar:.3g} bar: {type(e).__name__}"); return None, None
    return m, g
def converge(m, state, pbar):
    m.initialize_solver(state, solver_options=SOLVER); set_co2(m, pbar)       # restart from `state` with the (new) CO2 boundary
    while True:
        give_up, conv = m.iterate()
        if give_up or conv: break
    return conv

history = []; m = None
x = np.log10(meta["pCO2_bar"])
if os.environ.get("START_LOG10"): x = float(os.environ["START_LOG10"])      # optional: start the CO2 search from a higher CO2
# find a CO2 at which the new mode's climate converges to start from (hot points can fail at low CO2)
cands = [x] + [min(np.log10(PMAX), x + k) for k in (1, 2, 3)]
for xc in cands:
    m, g = build(10**xc)
    if m is not None: x = xc; break
if m is None:
    json.dump({"S": S, "failed": "initial_guess failed at all starting CO2 values", "tried_pCO2_bar": [10**c for c in cands]}, open(f"out/S{S:07.2f}.json", "w"), indent=1); sys.exit(0)
state = g
for it in range(8):
    pbar = 10 ** x
    try:
        conv = converge(m, state, pbar)
    except Exception as e:            # a climate solve failed during the coupling iterations: bisect back toward the last good CO2
        nfail = globals().get("nfail", 0) + 1; globals()["nfail"] = nfail
        log(f"coupling failed at pCO2={pbar:.3g} bar ({type(e).__name__}); bisecting toward the last converged CO2")
        if not history or nfail > 3: log("giving up: no converged state at the lowest CO2 reached"); break
        x = 0.5 * (x + np.log10(history[-1]["pCO2_bar"])); continue
    T = float(m.c.T_surf); state = m.model_state()
    a = m.pc.mole_fraction_dict()
    history.append(dict(pCO2_bar=pbar, T_surf_K=T, fixed_point_converged=bool(conv), chemistry_converged=m.solver.chemistry_converged,
                        coupling_iterations=m.solver.niterations, P_total_bar=float(m.pc.wrk.surface_pressure),
                        vmr={k: float(a[k][0]) for k in ("CO2", "CO", "O2", "H2O", "CH4", "H2", "N2")}, wall_s=time.time() - t0))
    log(f"outer {it}: pCO2={pbar:.4g} bar -> T_surf={T:.2f} K (converged={conv}, {m.solver.niterations} coupling its), CO={a['CO'][0]:.2e}")
    err = T - T_TARGET
    if abs(err) < TOL_K and conv: break
    if len(history) >= 2 and abs(np.log10(history[-1]["pCO2_bar"]) - np.log10(history[-2]["pCO2_bar"])) > 1e-6:
        slope = (history[-1]["T_surf_K"] - history[-2]["T_surf_K"]) / (np.log10(history[-1]["pCO2_bar"]) - np.log10(history[-2]["pCO2_bar"]))
        if not (2.0 < slope < 150.0): slope = 25.0
    else: slope = 25.0                                        # K per decade of CO2, first guess
    xn = float(np.clip(x - err / slope, np.log10(PMIN), np.log10(PMAX)))
    xn = x + float(np.clip(xn - x, -1.5, 1.5))
    if abs(xn - x) < 1e-3:
        log(f"at the CO2 bound ({10**x:.3g} bar) with T-288 = {err:+.2f} K: target not reachable"); break
    x = xn
if not history:
    json.dump({"S": S, "failed": "no converged coupled state", "tried_pCO2_bar": [10**c for c in cands]}, open(f"out/S{S:07.2f}.json", "w"), indent=1); sys.exit(0)
final = history[-1]; reachable = abs(final["T_surf_K"] - T_TARGET) < 1.0
out = dict(S_Wm2=S, star="GJ876", haze=True, photochem_version=photochem.__version__, baseline_run=rid,
           baseline_pCO2_bar=meta["pCO2_bar"], baseline_T_surf_K=meta["T_surf_K"], baseline_vmr=meta["vmr"],
           final=final, reached_288K=bool(reachable), history=history, wall_s=time.time() - t0,
           climate="Photochem 0.9.3 CoupledPhotochemistryClimate: all photochemical gases in the radiative-convective model; tropopause pressure fixed at the baseline run's value")
d_out = f"{REPO}/fuller_climate_haze/GJ876_S{S:07.2f}_FullClimate093"; os.makedirs(d_out, exist_ok=True)
m.pc.out2atmosphere_txt(f"{d_out}/photochem_atmosphere.txt", number_of_decimals=10, overwrite=True)
c = m.c
np.savetxt(f"{d_out}/climate_profile.txt", np.c_[c.P / 1e6, c.T, c.z / 1e5], fmt="%.10e",
           header=f"columns: pressure_bar temperature_K altitude_km (top of atmosphere last)\nT_surf_K {c.T_surf:.10f}\nP_surf_bar {c.P_surf/1e6:.10e}")
json.dump(out, open(f"{d_out}/result.json", "w"), indent=1)
json.dump(out, open(f"out/S{S:07.2f}.json", "w"), indent=1)
log(f"DONE T_surf={final['T_surf_K']:.2f} K at pCO2={final['pCO2_bar']:.4g} bar reached_288K={reachable}")
