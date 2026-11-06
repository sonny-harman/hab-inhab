"""Haze formation rates for the saved haze runs -> haze_diagnostics.csv.
Restarts each saved run (seconds), reads Photochem's column-integrated production of the three haze
particles, and converts to carbon atoms using the particles' formation reactions
(HCaer1 from C2H + C4H2 = 6 C, HCaer2 from H2CN + HCN = 2 C, HCaer3 from C4H + HCCCN = 7 C)."""
import csv, glob, json, os
from sweeplib import *
CARBON = {"HCaer1": 6, "HCaer2": 2, "HCaer3": 7}
C_SOURCE = 3.0e10 + 1.0e9 + 1.7e8 + 2*1.3e8      # CH4 + CO + HCN + CH3CN surface fluxes, C atoms/cm2/s   # reduced-carbon input; haze carbon can exceed it because fixed-pCO2 CO2 also supplies carbon
def haze_rows(saved_dir="saved_inputs"):
    rows = []
    for p in sorted(glob.glob(f"{saved_dir}/*_haze*/*_run.json")):
        m = json.load(open(p)); rid = m["run_id"]
        pc, rep = reproduce_run(rid, saved_dir, max_wall=900, log=lambda *a: None)
        prod = {n: float(pc.production_and_loss(n, pc.wrk.usol).integrated_total_production) for n in CARBON}
        c = sum(CARBON[n] * v for n, v in prod.items())
        rows.append({"run_id": rid, "S_Wm2": m["S_Wm2"], "pCO2_bar": m["pCO2_bar"], "T_surf_K": m["T_surf_K"],
                     "particle_radius_cm": m["haze"]["particle_radius_cm"],
                     **{f"{n}_prod_cm-2_s-1": v for n, v in prod.items()},
                     "haze_carbon_atoms_cm-2_s-1": c, "reduced_carbon_input_atoms_cm-2_s-1": C_SOURCE,
                     "haze_carbon_over_reduced_carbon_input": c / C_SOURCE})
    return rows
if __name__ == "__main__":
    rows = haze_rows()
    with open("haze_diagnostics.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    for r in rows: print(r["run_id"], f'{r["haze_carbon_atoms_cm-2_s-1"]:.2e}', f'{r["haze_carbon_over_reduced_carbon_input"]:.3f}')
