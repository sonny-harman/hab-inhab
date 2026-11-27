import json, glob, csv
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = "/home/user/hab-inhab/coupled093_baseline"
rs = [json.load(open(p)) for p in glob.glob(f"{ROOT}/*/result.json")]
rs.sort(key=lambda r: (r["star"], -r["S_Wm2"]))
cols = ["run_id","star","S_Wm2","haze","pCO2_bar","status","iterations","T_old_K","T_new_K","dT_K","CO_old","CO_new","CH4_old","CH4_new","H2O_old","H2O_new","P_old_bar","P_new_bar","minor_gases","wall_s","error"]
with open(f"{ROOT}/map_table.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(cols)
    for r in rs:
        o = r["old"]; ok = r["status"] != "failed"
        w.writerow([r["run_id"], r["star"], r["S_Wm2"], r["haze"], f"{r['pCO2_bar']:.4g}", r["status"], r.get("iterations", ""), f"{o['T_surf_K']:.3f}",
                    f"{r['T_surf_K']:.3f}" if ok else "", f"{r['dT_K_vs_0.9.1']:.3f}" if ok else "",
                    f"{o['vmr']['CO']:.3e}", f"{r['vmr']['CO']:.3e}" if ok else "", f"{o['vmr']['CH4']:.3e}", f"{r['vmr']['CH4']:.3e}" if ok else "",
                    f"{o['vmr']['H2O']:.3e}", f"{r['vmr']['H2O']:.3e}" if ok else "", f"{o['P_total_bar']:.3f}", f"{r['P_total_bar']:.3f}" if ok else "",
                    ";".join(f"{k}={v:.1e}" for k, v in r.get("minor_gases_gt_1e7", {}).items()), f"{r['wall_s']:.0f}", r.get("error", "")[:80]])
fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.3), sharey=True, constrained_layout=True)
for ax, star, ttl in [(axs[0], "Sun", "Sun at 3 Ga"), (axs[1], "GJ876", "GJ 876")]:
    rr = [r for r in rs if r["star"] == star and r["status"] != "failed"]
    ax.axhline(0, color="gray", lw=.8, ls=":")
    for r in rr:
        capped = r["pCO2_bar"] >= 9.99 or r["old"]["target_status"] != "ok"
        ax.plot(r["S_Wm2"], r["dT_K_vs_0.9.1"], "o", color="#2a78d6" if not r["haze"] else "#eb6834", mfc="white" if abs(r["old"]["T_surf_K"] - 288) > 1.5 else ("#2a78d6" if not r["haze"] else "#eb6834"), ms=6)
    for r in [r for r in rs if r["star"] == star and r["status"] == "failed"]:
        ax.plot(r["S_Wm2"], 0, "x", color="k", ms=9, mew=1.5)
    ax.set_title(ttl); ax.set_xlabel("Instellation (W/m$^2$)"); ax.grid(alpha=.2)
axs[0].set_ylabel("T(0.9.3 coupled) $-$ T(0.9.1 sweep)  (K)")
from matplotlib.lines import Line2D
axs[1].legend(handles=[Line2D([], [], marker="o", ls="", color="#2a78d6", label="no-haze network, held 288 K"), Line2D([], [], marker="o", ls="", color="#eb6834", label="haze network (hot GJ 876)"),
              Line2D([], [], marker="o", ls="", color="gray", mfc="white", label="0.9.1 run did not hold 288 K"), Line2D([], [], marker="x", ls="", color="k", label="0.9.3 failed")], fontsize=8)
fig.savefig(f"{ROOT}/map_dT.png", dpi=160)
for r in rs:
    ok = r["status"] != "failed"
    print(f"{r['star']:5s} {r['S_Wm2']:8.2f} haze={str(r['haze'])[0]} pCO2={r['pCO2_bar']:9.3e} " + (f"T {r['old']['T_surf_K']:7.2f} -> {r['T_surf_K']:7.2f} dT {r['dT_K_vs_0.9.1']:+6.2f}  CO {r['old']['vmr']['CO']:.1e}->{r['vmr']['CO']:.1e} P {r['old']['P_total_bar']:.2f}->{r['P_total_bar']:.2f} {r['status']} {r['iterations']}its" if ok else f"FAILED {r['stage']}: {r['error'][:70]}"))
