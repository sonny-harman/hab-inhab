"""Climate-only test at fixed surface composition: GJ876 S=1125 haze run (pCO2 8.944e-4 bar, CH4 2.24e-5, H2 3.88e-3)."""
import sys, json, numpy as np
sys.path.insert(0, "/home/user/hab-inhab")
import photochem
print("photochem", photochem.__version__)
from photochem.clima import AdiabatClimate
import sweeplib
rid = "GJ876_S1125.00_CO2-8.944e-04bar_haze"; d = f"/home/user/hab-inhab/saved_inputs/{rid}"
meta = json.load(open(f"{d}/{rid}_run.json")); f = meta["saved_inputs"]["files"]
shared = "/home/user/hab-inhab/saved_inputs/shared"
# (A) the sweep's own climate model/species (5 gases) in this photochem version
pt = sweeplib.Point(shared, f"{d}/{f['flux']}")
T, ok = pt.rce_T(meta["pCO2_bar"], meta["vmr"]["CH4"], meta["vmr"]["H2"])
print(f"(A) sweep climate, 5 gases, this version: T_surf={T:.2f} K ok={ok}  (0.9.1 result saved in run: {meta['T_surf_K']:.2f} K)")
# (B) the coupled mode's climate object (all photochemical gases), only the same 5 gases populated
from photochem.extensions.coupled import CoupledPhotochemistryClimate
m = CoupledPhotochemistryClimate(f"{shared}/mechanism_haze.yaml", f"{d}/{f['settings']}", f"{d}/{f['flux']}")
c = m.c; names = list(c.species_names); print("coupled-mode climate species:", len(names), names[:12], "...")
c.verbose = False
pp = pt.partial_pressures(meta["pCO2_bar"], meta["vmr"]["CH4"], meta["vmr"]["H2"])   # dyn/cm2, in the 5-species order
base = dict(zip(pt.names, pp))
def run(extra=None, label=""):
    P_i = np.full(len(names), 1e-15)
    for i, n in enumerate(names):
        if n in base: P_i[i] = base[n]
        if extra and n in extra: P_i[i] = extra[n]
    ok = c.RCE_robust(P_i)
    print(f"(B){label}: T_surf={c.T_surf:.2f} K ok={ok}", flush=True)
run(label=" coupled-mode climate, same 5 gases")
tot = sum(base.values())
run({"CO": 6.2e-3*tot}, " + CO at 6.2e-3")
print("RH:", dict(zip(names[:3], c.RH[:3])), "T_trop/solve_for_T_trop:", getattr(c, "solve_for_T_trop", None), getattr(c, "T_trop", None))
