"""For each finished spot-check: climate-only T at the NEW composition with (a) my 5 gases, (b) 5 gases + CO."""
import sys, json, glob, numpy as np
sys.path.insert(0, "/home/user/hab-inhab")
from photochem.extensions.coupled import CoupledPhotochemistryClimate
REPO = "/home/user/hab-inhab"; shared = f"{REPO}/saved_inputs/shared"
rows = []
for p in sorted(glob.glob("out/*.json")):
    r = json.load(open(p))
    if "failed" in r: continue
    rid = r["run_id"]; d = f"{REPO}/saved_inputs/{rid}"; meta = json.load(open(f"{d}/{rid}_run.json")); f = meta["saved_inputs"]["files"]
    mech = f"{shared}/{'mechanism_haze.yaml' if r['haze'] else 'mechanism.yaml'}"
    m = CoupledPhotochemistryClimate(mech, f"{d}/{f['settings']}", f"{d}/{f['flux']}"); c = m.c; c.verbose = False
    names = list(c.species_names); v = r["vmr_new"]; Pt = r["P_total_bar_new"] * 1e6
    def P_i(gases):
        out = np.full(len(names), 1e-15)
        for i, n in enumerate(names):
            if n == "H2O": out[i] = 270e6
            elif n in gases: out[i] = v[n] * Pt
        return out
    res = {}
    for label, gases in (("5gas", ("CO2", "N2", "CH4", "H2")), ("5gas+CO", ("CO2", "N2", "CH4", "H2", "CO"))):
        ok = c.RCE_robust(P_i(gases)); res[label] = (float(c.T_surf), bool(ok))
    print(f"{rid}\n   0.9.1 sweep T={r['T_surf_K_old']:.2f}  | new coupled mode T={r['T_surf_K_new']:.2f}  | climate-only at new composition: 5 gases {res['5gas'][0]:.2f}, +CO {res['5gas+CO'][0]:.2f}   (CO vmr new {v['CO']:.2e})", flush=True)
    rows.append(dict(run_id=rid, T_old=r["T_surf_K_old"], T_new=r["T_surf_K_new"], T_5gas_newcomp=res["5gas"][0], T_5gas_CO_newcomp=res["5gas+CO"][0], CO_new=v["CO"]))
json.dump(rows, open("decompose.json", "w"), indent=1)
