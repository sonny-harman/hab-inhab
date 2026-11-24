import json, glob, csv, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = "/tmp/claude-0/-home-user-hab-inhab/e6bce99a-9bee-504b-9b23-ed0404f973d7/scratchpad/fuller/out"
REPO = "/home/user/hab-inhab/fuller_climate_haze"
res = {}
for p in sorted(glob.glob(f"{OUT}/S*.json")):
    r = json.load(open(p)); res[float(p.split("/S")[-1][:-5])] = r
rows = []
for S, r in sorted(res.items()):
    if "failed" in r:
        rows.append(dict(S=S, note=r["failed"])); continue
    f = r["final"]; b = r["baseline_vmr"]
    rows.append(dict(S=S, pCO2_sweep=r["baseline_pCO2_bar"], T_sweep=r["baseline_T_surf_K"], pCO2_fuller=f["pCO2_bar"], T_fuller=f["T_surf_K"],
                     reached_288K=r["reached_288K"], CO_sweep=b["CO"], CO_fuller=f["vmr"]["CO"], CH4_sweep=b["CH4"], CH4_fuller=f["vmr"]["CH4"],
                     H2O_sweep=b["H2O"], H2O_fuller=f["vmr"]["H2O"], P_total_fuller=f["P_total_bar"], evaluations=len(r["history"]),
                     first_T_at_sweep_pCO2=r["history"][0]["T_surf_K"] if abs(r["history"][0]["pCO2_bar"]-r["baseline_pCO2_bar"])/r["baseline_pCO2_bar"] < 1e-3 else "", note=""))
cols = ["S","pCO2_sweep","T_sweep","pCO2_fuller","T_fuller","reached_288K","first_T_at_sweep_pCO2","CO_sweep","CO_fuller","CH4_sweep","CH4_fuller","H2O_sweep","H2O_fuller","P_total_fuller","evaluations","note"]
with open(f"{REPO}/comparison.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(cols)
    for r in rows: w.writerow([f"{r.get(c):.6g}" if isinstance(r.get(c), float) else r.get(c, "") for c in cols])
ok = [r for r in rows if "pCO2_fuller" in r]
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
x = [r["S"] for r in ok]
ax[0].plot(x, [r["pCO2_sweep"] for r in ok], "s--", color="#7a7a7a", label="sweep (5-gas climate, Photochem 0.9.1)")
ax[0].plot([r["S"] for r in ok if r["reached_288K"]], [r["pCO2_fuller"] for r in ok if r["reached_288K"]], "o-", color="#2a78d6", label="fuller climate, holds 288 K")
ax[0].plot([r["S"] for r in ok if not r["reached_288K"]], [r["pCO2_fuller"] for r in ok if not r["reached_288K"]], "o", mfc="white", color="#2a78d6", label="fuller climate, 288 K not reached (CO2 floor)")
ax[0].axhline(1e-5, color="gray", lw=.8, ls=":"); ax[0].set_yscale("log"); ax[0].set_xlabel("Instellation (W/m$^2$)"); ax[0].set_ylabel("CO$_2$ partial pressure (bar)"); ax[0].legend(fontsize=8)
ax[1].plot(x, [r["T_sweep"] for r in ok], "s--", color="#7a7a7a", label="sweep")
ax[1].plot(x, [r["T_fuller"] for r in ok], "o-", color="#2a78d6", label="fuller climate"); ax[1].axhline(288, color="gray", lw=.8, ls=":")
ax[1].set_xlabel("Instellation (W/m$^2$)"); ax[1].set_ylabel("Surface temperature (K)"); ax[1].legend(fontsize=8)
fig.savefig(f"{REPO}/comparison.png", dpi=160)
for r in rows: print(r)
