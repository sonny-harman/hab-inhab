"""Climate-only: feed the coupled-mode climate model the FULL 0.9.1 photochemical composition (all gases, vertical profiles)."""
import sys, json, numpy as np
from photochem.extensions.coupled import CoupledPhotochemistryClimate
REPO = "/home/user/hab-inhab"; shared = f"{REPO}/saved_inputs/shared"
for rid in sys.argv[1:]:
    d = f"{REPO}/saved_inputs/{rid}"; meta = json.load(open(f"{d}/{rid}_run.json")); f = meta["saved_inputs"]["files"]
    haze = meta.get("haze", {}).get("enabled", False)
    m = CoupledPhotochemistryClimate(f"{shared}/{'mechanism_haze.yaml' if haze else 'mechanism.yaml'}", f"{d}/{f['settings']}", f"{d}/{f['flux']}")
    c = m.c; c.verbose = False; names = list(c.species_names)
    hdr = open(f"{d}/{f['atmosphere']}").readline().split(); arr = np.loadtxt(f"{d}/{f['atmosphere']}", skiprows=1)
    col = {h: arr[:, i] for i, h in enumerate(hdr)}
    press = col["press"] * 1e6                                   # dyn/cm2, surface first
    gas = np.maximum(np.array([col[n] for n in names]), 1e-30); gas /= gas.sum(axis=0)
    Psurf = meta["P_total_bar"] * 1e6
    fixed = {"H2O": 270e6, "CO2": meta["pCO2_bar"] * 1e6}
    dry = np.array([n not in ("H2O", "CO2") for n in names])
    P_i = np.array([fixed[n] if n in fixed else gas[i, 0] * Psurf for i, n in enumerate(names)])
    dry_mix = gas[dry] / gas[dry].sum(axis=0)
    cdm = {"pressure": press.copy()}
    for n, v in zip(np.asarray(names)[dry], dry_mix): cdm[str(n)] = v.copy()
    ok = c.RCE_robust(P_i, cdm)
    top = sorted(((gas[i, 0], n) for i, n in enumerate(names) if n not in ("N2", "H2O", "CO2", "He")), reverse=True)[:8]
    print(f"{rid}\n   0.9.1 sweep T={meta['T_surf_K']:.2f} | coupled-mode climate on the full 0.9.1 composition: T={c.T_surf:.2f} ok={ok}\n   top surface gases: " +
          ", ".join(f"{n} {v:.1e}" for v, n in top), flush=True)
