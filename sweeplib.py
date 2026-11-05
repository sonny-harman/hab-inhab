"""Helpers for the instellation sweep notebook (written by the notebook, importable).

For a given star and instellation S, `run_point` finds the CO2 partial pressure that holds the
surface temperature at T_target in a coupled climate-photochemistry model of an anoxic planet.
"""
import os, json, time, shutil, sys
import numpy as np
import yaml
from scipy.optimize import brentq
from photochem import EvoAtmosphere
from photochem.clima import AdiabatClimate
from photochem.utils import (stars, zahnle_rx_and_thermo_files,
                             species_file_for_climate, settings_file_for_climate)

planet_mass, planet_radius, surface_albedo = 5.972e27, 6.371e8, 0.2
relative_humidity = 0.5
T_target = 288.0                      # K
ATOL = 1.0e-18                        # photochemical integrator absolute tolerance (mixing ratio)
TRACKED = ["CO", "O2", "H2O", "CH4", "H2", "CO2"]   # species whose surface values define "steady"
STAR_TAG = {"Sun": "Sun3Ga", "GJ876": "GJ876"}
SHARED_FILES = ["mechanism.yaml", "climate_species.yaml", "climate_settings.yaml"]
pN2_bar = 1.0
pCO2_min_bar, pCO2_max_bar = 1.0e-5, 10.0

def setup_files(workdir):
    os.makedirs(workdir, exist_ok=True)
    mech = f"{workdir}/mechanism.yaml"
    zahnle_rx_and_thermo_files(atoms_names=["H","He","N","O","C","S"], rxns_filename=mech,
                               thermo_filename=None, remove_reaction_particles=True)
    def bc(name, lower, upper=None):
        return {"name": name, "lower-boundary": lower, "upper-boundary": upper or {"type":"veff","veff":0.0}}
    def src(vdep, flux): return {"type":"vdep + dist flux","vdep":vdep,"flux":flux,"height":-1.0}
    bcs = [bc("H2O",{"type":"press","press":270.0e6}), bc("N2",{"type":"press","press":pN2_bar*1e6}),
           bc("CO2",{"type":"press","press":0.01e6}), bc("O2",{"type":"flux","flux":0.0}),
           bc("H2",src(0.0,3.0e10)), bc("CH4",{"type":"flux","flux":3.0e10}),
           bc("CO",src(1.0e-8,1.0e9)), bc("SO2",src(1.0,1.0e9)), bc("H2S",src(0.015,1.0e8)),
           bc("OCS",src(0.003,5.0e7)), bc("NO",src(0.01,3.0e8)),
           bc("NO2",{"type":"vdep","vdep":0.05}), bc("NO3",{"type":"vdep","vdep":0.05}),
           bc("HNO3",{"type":"vdep","vdep":4.0}), bc("NH3",src(1.0,1.5e10)),
           bc("H2SO4",src(1.0,7.0e8)), bc("HCN",src(0.13,1.7e8)), bc("CH3CN",src(0.13,1.3e8)),
           bc("N2O",src(0.0,0.0))]
    settings = {"atmosphere-grid":{"number-of-layers":100},
                "planet":{"planet-mass":planet_mass,"planet-radius":planet_radius,"surface-albedo":surface_albedo,
                          "solar-zenith-angle":60.0,"hydrogen-escape":{"type":"diffusion limited"},
                          "water":{"gas-rainout":True,"rainfall-rate":1,"tropopause-altitude":1.1e6}},
                "particles":[{"name":"H2Oaer","RH-condensation":relative_humidity}],
                "boundary-conditions":bcs}
    with open(f"{workdir}/photochem_settings.yaml","w") as f: yaml.safe_dump(settings,f,sort_keys=False)
    species_file_for_climate(f"{workdir}/climate_species.yaml",["H2O","CO2","N2","CH4","H2"],["H2O","CO2"])
    settings_file_for_climate(f"{workdir}/climate_settings.yaml",planet_mass,planet_radius,surface_albedo)
    return mech

def cache_stellar_spectrum(star, workdir):
    "Download an observed MUSCLES spectrum once and cache it (needs internet access)."
    path = f"{workdir}/{star}_muscles_raw.npz"
    if not os.path.exists(path):
        wv, F = stars.get_muscles_spectrum(star, False, 1000)
        np.savez(path, wv=wv, F=F)
    return path

def write_flux_file(star, S, path, workdir=None):
    "Spectrum of `star` scaled so total instellation is S (W/m^2)."
    if star == "Sun":
        stars.solar_spectrum(outputfile=path, age=3.0, stellar_flux=S, scale_before_age=False)
    else:
        d = np.load(cache_stellar_spectrum(star, workdir))
        F = stars.scale_spectrum_to_planet(d["wv"], d["F"], None, S)
        wv, F = stars.rebin_to_needed_resolution(d["wv"], F)
        stars.save_photochem_spectrum(wv, F, path, scale_to_planet=False)

def eddy_diffusion(P_bar):
    return np.where(P_bar > 0.1, 1.0e5, np.minimum(1.0e5*(0.1/P_bar)**0.5, 1.0e7))

class Point:
    "Climate model + helpers for one instellation."
    def __init__(self, workdir, flux_file):
        self.clim = AdiabatClimate(f"{workdir}/climate_species.yaml", f"{workdir}/climate_settings.yaml", flux_file)
        self.clim.verbose = False
        self.names = self.clim.species_names
        self.clim.RH = np.full(len(self.names), relative_humidity)
        self.clim.T_trop = 200.0
    def partial_pressures(self, pCO2_bar, f_CH4, f_H2):
        p_dry = (pN2_bar + pCO2_bar) / (1.0 - f_CH4 - f_H2)
        d = {"N2":pN2_bar, "CO2":pCO2_bar, "CH4":f_CH4*p_dry, "H2":f_H2*p_dry, "H2O":270.0}
        return np.array([d.get(n,1e-15) for n in self.names])*1e6
    def approx_T(self, pCO2_bar, f_CH4, f_H2):
        "Fast approximate climate (adiabat + isothermal stratosphere, solved tropopause). None if the solver fails."
        c = self.clim; c.solve_for_T_trop = True
        P_i = self.partial_pressures(pCO2_bar, f_CH4, f_H2)
        for guess in (288.0, 230.0, 340.0, 260.0):
            for T_trop in (200.0, 150.0):
                c.T_trop = T_trop
                try:
                    return c.surface_temperature(P_i, T_guess=guess)
                except Exception:
                    pass
        return None
    def rce_T(self, pCO2_bar, f_CH4, f_H2):
        "Full radiative-convective equilibrium started from the approximate solution. (T_surf, converged)."
        T0 = self.approx_T(pCO2_bar, f_CH4, f_H2)
        if T0 is None: return None, False
        c = self.clim; P_i = self.partial_pressures(pCO2_bar, f_CH4, f_H2)
        try:
            ok = c.RCE(P_i, c.T_surf, c.T, c.convecting_with_below)
        except Exception:
            ok = False
        if not ok:
            c.surface_temperature(P_i, T_guess=T0)   # restore the approximate state
            return T0, False
        return c.T_surf, True
    def solve_pCO2(self, f_CH4, f_H2, tol_K=0.25, max_refine=5):
        "pCO2 (bar) giving T_surf = T_target in the full RCE model. Returns (pCO2, status)."
        lo, hi = np.log10(pCO2_min_bar), np.log10(pCO2_max_bar)
        def fa(x):
            T = self.approx_T(10**x, f_CH4, f_H2)
            return np.nan if T is None else T - T_target
        xs = np.linspace(lo, hi, 25); ts = np.array([fa(x) for x in xs])
        good = np.isfinite(ts)
        if good[0] and ts[0] > 0:
            return 10**lo, "target not reachable: too warm even at minimum CO2"
        idx = [i for i in range(len(xs)-1) if good[i] and good[i+1] and ts[i] <= 0 <= ts[i+1]]
        if not idx:
            # 288 K is out of reach: use the CO2 that gives the warmest surface (maximum greenhouse)
            xbest = xs[int(np.nanargmax(np.where(good, ts, -np.inf)))]
            xf = np.linspace(max(lo, xbest-0.3), min(hi, xbest+0.3), 7)
            tf = np.array([fa(x) for x in xf]); xbest = xf[int(np.nanargmax(np.where(np.isfinite(tf), tf, -np.inf)))]
            return 10**xbest, "target not reachable: maximum-greenhouse CO2 used"
        i = idx[0]
        x = brentq(fa, xs[i], xs[i+1], xtol=0.002)
        # refine against the full RCE model with Newton-like steps
        for _ in range(max_refine):
            T, ok = self.rce_T(10**x, f_CH4, f_H2)
            if T is None: break
            err = T - T_target
            if abs(err) < tol_K: break
            d = 0.05
            Ta, Tb = self.approx_T(10**(x-d), f_CH4, f_H2), self.approx_T(10**(x+d), f_CH4, f_H2)
            slope = (Tb - Ta) / (2*d) if (Ta is not None and Tb is not None) else 40.0
            slope = max(slope, 5.0)
            x = float(np.clip(x - err/slope, lo, hi))
        return 10**x, "ok"

def photochem_steady(pc, tracked, tol=0.01, trace_floor=1e-12, trace_dex=0.05, min_time=1e12, decade=10.0,
                     check_every=25, max_wall=7200, runaway_vmr=5e-3,
                     expected=("N2", "CO2", "H2O", "H2", "CO", "CH4", "He")):
    """Integrate to a practical steady state.
    Steady = Photochem's own criterion, OR the surface mixing ratios of the `tracked` species
    changed by less than `tol` (1%) over the last decade in model time (and t > min_time).
    Trace species (mixing ratio < trace_floor, e.g. O2 at 1e-16) are instead compared in log space
    (change < trace_dex dex), since their relative noise is large and they do not matter physically.
    Also stops early with reason "runaway: <species>" if a gas other than the expected major/minor
    species exceeds `runaway_vmr` at the surface: in these runs hydrocarbons such as C4H4 can pile up
    without bound when there is little oxidant, because haze formation is not modelled, so no steady
    state exists.
    Returns (reason, model time, steps, wall time, species limiting convergence)."""
    pc.var.nsteps_before_giveup = 1000000       # the default (1e5) can be exhausted before steady state
    pc.initialize_robust_stepper(pc.wrk.usol)
    hist = []; n = 0; t0 = time.time(); reason = "max wall time"; worst = ""
    def surf():
        a = pc.mole_fraction_dict(); return np.array([a[s][0] for s in tracked])
    while True:
        give_up, conv = pc.robust_step(); n += 1
        if conv: reason = "photochem criterion"; break
        if give_up: reason = "give_up"; break
        if n % check_every == 0:
            tn = pc.wrk.tn; v = surf(); hist.append((tn, v))
            if n % (4*check_every) == 0:
                a = pc.mole_fraction_dict()
                bad = [(a[k][0], k) for k in pc.dat.species_names[:pc.dat.nq]
                       if k not in expected and not k.endswith("aer") and a[k][0] > runaway_vmr]
                if bad:
                    worst = max(bad)[1]; reason = f"runaway: {worst} accumulating"; break
            if tn > min_time:
                old = [h for h in hist if h[0] <= tn/decade]
                if old:
                    v0 = old[-1][1]
                    tiny = (np.abs(v) < trace_floor) & (np.abs(v0) < trace_floor)
                    rel = np.abs(v - v0) / np.maximum(np.maximum(np.abs(v), np.abs(v0)), 1e-300)
                    dex = np.abs(np.log10(np.maximum(np.abs(v), 1e-300)) - np.log10(np.maximum(np.abs(v0), 1e-300)))
                    score = np.where(tiny, dex/trace_dex, rel/tol)      # < 1 means converged
                    worst = tracked[int(np.argmax(score))]
                    if score.max() < 1.0: reason = "composition steady over last decade"; break
        if time.time() - t0 > max_wall: break
    pc.destroy_stepper()
    return reason, pc.wrk.tn, n, time.time() - t0, worst

def surface_summary(pc):
    a = pc.mole_fraction_dict()
    return {k: float(a[k][0]) for k in ["CO2","CO","O2","H2O","CH4","H2","N2"]}


# ---------------------------------------------------------------------------
# Saving and reproducing runs
# ---------------------------------------------------------------------------
def run_id(star, S, pCO2_bar):
    """Name that identifies a run, e.g. Sun3Ga_S1400.00_CO2-4.651e-05bar
    (star and spectrum, instellation in W/m^2, CO2 surface partial pressure in bar; e6.003e00 means 6.003 bar)."""
    return f"{STAR_TAG[star]}_S{S:07.2f}_CO2-{pCO2_bar:.3e}bar".replace("e+", "e")

def is_done(star, S, results_dir, saved_dir):
    "A point is done when its result file exists and its saved inputs are present."
    p = f"{results_dir}/{star}_S{S:.1f}.json"
    if not os.path.exists(p): return False
    rid = json.load(open(p)).get("run_id")
    return bool(rid) and os.path.exists(f"{saved_dir}/{rid}/{rid}_run.json")

def write_saved_inputs(saved_dir, rid, star, S, pCO2, pc, pt, workdir, pdir, flux_file, P_prof, T_prof, trop_p):
    """Write everything needed to restart this run and to plot its vertical structure.
    Must be called after the final photochemical solve and after the final climate solve."""
    d = f"{saved_dir}/{rid}"; os.makedirs(d, exist_ok=True)
    shared = f"{saved_dir}/shared"; os.makedirs(shared, exist_ok=True)
    for f in SHARED_FILES:
        if not os.path.exists(f"{shared}/{f}"): shutil.copy(f"{workdir}/{f}", f"{shared}/{f}")
    raw = f"{workdir}/{star}_muscles_raw.npz"
    if os.path.exists(raw) and not os.path.exists(f"{shared}/{os.path.basename(raw)}"):
        shutil.copy(raw, f"{shared}/{os.path.basename(raw)}")
    files = {k: f"{rid}_{v}" for k, v in dict(
        atmosphere="photochem_atmosphere.txt", settings="photochem_settings.yaml", flux="flux.txt",
        profile="PTKzz_profile.txt", climate="climate_profile.txt").items()}
    # 1. final photochemical state: altitude, pressure, density, T, Kzz, every species and particle (Photochem's own format)
    pc.out2atmosphere_txt(f"{d}/{files['atmosphere']}", number_of_decimals=10, overwrite=True)
    # 2. photochemical settings with this run's CO2 lower boundary (the sweep sets it at run time)
    st = yaml.safe_load(open(f"{pdir}/photochem_settings.yaml"))
    for bc in st["boundary-conditions"]:
        if bc["name"] == "CO2": bc["lower-boundary"] = {"type": "press", "press": float(pCO2 * 1e6)}
    with open(f"{d}/{files['settings']}", "w") as f: yaml.safe_dump(st, f, sort_keys=False)
    # 3. stellar flux used (spectrum scaled to this instellation)
    shutil.copy(flux_file, f"{d}/{files['flux']}")
    # 4. persistent pressure-temperature-eddy profile that the last photochemical solve used
    np.savetxt(f"{d}/{files['profile']}", np.c_[P_prof, T_prof, eddy_diffusion(P_prof / 1e6)], fmt="%.12e",
               header=("columns: pressure (dyn/cm^2, surface first), temperature (K), eddy diffusion Kzz (cm^2/s)\n"
                       f"tropopause_pressure_dyn_cm2 {trop_p:.12e}\ntarget_toa_pressure_dyn_cm2 {P_prof[-1]:.12e}"))
    # 5. final climate (radiative-convective) profile for the final composition
    c = pt.clim
    cols = ["pressure_bar", "temperature_K", "altitude_km", "convecting_with_below"] + [f"vmr_{n}" for n in c.species_names]
    arr = np.c_[c.P / 1e6, c.T, c.z / 1e5, np.asarray(c.convecting_with_below, float), np.asarray(c.f_i)]
    np.savetxt(f"{d}/{files['climate']}", arr, fmt="%.10e",
               header=("climate model layers, top of atmosphere last. columns: " + " ".join(cols) + "\n"
                       f"T_surf_K {c.T_surf:.10f}\nP_surf_bar {c.P_surf / 1e6:.10e}\nP_trop_bar {c.P_trop / 1e6:.10e}"))
    return files

def read_profile(path):
    "Read a *_PTKzz_profile.txt: returns P (dyn/cm^2), T (K), Kzz (cm^2/s), tropopause pressure, target TOA pressure."
    hdr = [l[1:].split() for l in open(path) if l.startswith("#")]
    meta = {h[0]: float(h[1]) for h in hdr if len(h) == 2 and h[0].endswith("dyn_cm2")}
    P, T, K = np.loadtxt(path).T
    return P, T, K, meta["tropopause_pressure_dyn_cm2"], meta["target_toa_pressure_dyn_cm2"]

def read_climate_profile(path):
    "Read a *_climate_profile.txt into a dict of arrays (plus T_surf_K, P_surf_bar, P_trop_bar)."
    lines = [l[1:].strip() for l in open(path) if l.startswith("#")]
    names = lines[0].split("columns:")[1].split()
    arr = np.loadtxt(path)
    out = {n: arr[:, i] for i, n in enumerate(names)}
    for l in lines[1:]:
        k, v = l.split(); out[k] = float(v)
    return out

def read_atmosphere(path):
    """Read a *_photochem_atmosphere.txt into a dict of arrays keyed by column name: 'alt' (km), 'press' (bar),
    'den' (molecules/cm^3), 'temp' (K), 'eddy' (cm^2/s), the mixing ratio of every particle and gas (e.g. 'CO2', 'CO'),
    and particle radii ('H2Oaer_r', cm). Layers run from the surface upward. No model is needed."""
    with open(path) as f:
        names = f.readline().split()
    arr = np.loadtxt(path, skiprows=1)
    return {n: arr[:, i] for i, n in enumerate(names)}

def reproduce_run(rid, saved_dir="saved_inputs", max_wall=1800, log=print):
    """Rebuild a saved run from its files and let the photochemistry settle, which is quick because the
    saved state is already converged. Returns (pc, report): `pc` is the model in its final state and
    `report` compares the reproduced surface composition with the saved one."""
    d, shared = f"{saved_dir}/{rid}", f"{saved_dir}/shared"
    meta = json.load(open(f"{d}/{rid}_run.json"))
    f = meta["saved_inputs"]["files"]
    pc = EvoAtmosphere(f"{shared}/mechanism.yaml", f"{d}/{f['settings']}", f"{d}/{f['flux']}", f"{d}/{f['atmosphere']}")
    pc.var.verbose = 0
    pc.var.atol = meta["photochem"]["atol"]
    P, T, K, trop_p, target = read_profile(f"{d}/{f['profile']}")
    pc.set_press_temp_edd_profile(P, T, K, trop_p=trop_p, hydro_pressure=True,
                                  maintain_toa_pressure=True, target_pressure=target)
    reason, tn, n, wall, worst = photochem_steady(pc, TRACKED, max_wall=max_wall)
    now, saved = surface_summary(pc), meta["vmr"]
    rel = {k: abs(now[k] - saved[k]) / max(abs(saved[k]), 1e-300) for k in saved}
    report = {"run_id": rid, "stop": reason, "wall_s": wall, "model_time_s": tn,
              "reproduced": now, "saved": saved, "max_rel_diff_tracked": max(rel[k] for k in TRACKED),
              "rel_diff": rel, "P_total_bar": float(pc.wrk.surface_pressure), "P_total_bar_saved": meta["P_total_bar"],
              "saved_run_was_steady": bool(meta["photochem_steady"])}
    note = "" if meta["photochem_steady"] else " (this saved run had no steady state, so it keeps running away and differences are expected)"
    log(f"{rid}: {reason} in {wall:.0f} s; largest relative difference in a tracked species "
        f"{report['max_rel_diff_tracked']:.1e}{note}")
    return pc, report

def load_climate(rid, saved_dir="saved_inputs"):
    """Rebuild the saved climate model state for a run. First restarts the full radiative-convective solve from the
    saved profile; if the solver cannot (it can fail for cold, CO2-rich cases), it repeats the call sequence the sweep
    itself used, which reproduces the saved state. Returns (AdiabatClimate, converged)."""
    d, shared = f"{saved_dir}/{rid}", f"{saved_dir}/shared"
    meta = json.load(open(f"{d}/{rid}_run.json")); f = meta["saved_inputs"]["files"]
    pt = Point(shared, f"{d}/{f['flux']}"); c = pt.clim
    prof = read_climate_profile(f"{d}/{f['climate']}")
    P_i = np.array([meta["climate"]["surface_partial_pressures_dyn_cm2"][n] for n in c.species_names])
    # The climate model needs its internal state set up by a first (approximate) solve, as in the sweep itself;
    # then the full radiative-convective solve restarts from the saved profile.
    try:
        c.solve_for_T_trop = True
        c.surface_temperature(P_i, T_guess=prof["T_surf_K"])
        ok = c.RCE(P_i, prof["T_surf_K"], prof["temperature_K"], prof["convecting_with_below"].astype(bool))
    except Exception:
        ok = False
    if not ok:      # the sweep's own sequence (several starting guesses for the approximate solve, then RCE)
        T, ok = pt.rce_T(meta["pCO2_bar"], meta["vmr"]["CH4"], meta["vmr"]["H2"])
        ok = bool(ok) and T is not None
    return c, ok

def make_readme(saved_dir="saved_inputs"):
    """Write <saved_dir>/README.md (naming key, file formats, how to restart, and a full description of every
    saved simulation) and <saved_dir>/run_index.csv, generated from the *_run.json files."""
    import csv, glob
    metas = [json.load(open(p)) for p in sorted(glob.glob(f"{saved_dir}/*/*_run.json"))]
    metas.sort(key=lambda m: (m["star"], -m["S_Wm2"]))
    if not metas: raise RuntimeError("no saved runs found")
    ex = metas[0]; pl = ex["planet"]; v = ex["versions"]
    sci = lambda x: f"{x:.2e}"
    def status(m):
        if m["photochem_stop"].startswith("runaway"):
            return "NO STEADY STATE", f"hydrocarbons accumulate without bound ({m['photochem_stop'].split(': ')[1]}); this chemistry has no haze sink. The saved state is the one at the moment the runaway was detected, so its composition is not meaningful and it is saved only for diagnosis."
        if m["target_status"] != "ok" or abs(m["T_surf_K"] - m["T_target_K"]) > 1.5:
            return "288 K NOT REACHABLE", f"the model could not reach {m['T_target_K']:.0f} K with CO2 up to the 10 bar cap, so it used the CO2 pressure giving the warmest surface; the surface settles at {m['T_surf_K']:.0f} K (part of the CO2 may condense)."
        return "holds 288 K", "the coupled model holds the surface at the target temperature in a steady state."
    lines = []
    A = lines.append
    A("# Saved inputs for the instellation sweep\n")
    A("Each simulation in `archean_instellation_sweep.ipynb` has its own folder here, holding everything needed to\n"
      "(1) restart the photochemical model from its converged state, (2) reproduce that state quickly, and\n"
      "(3) plot the vertical structure (temperature, pressure, eddy diffusion and the mixing ratio of every gas). "
      "This file is generated by `make_readme()` from the `*_run.json` files, so it always matches what is saved.\n")
    A("## The scenario (same for every run)\n")
    A(f"- **Planet:** Earth's mass ({pl['mass_g']:.3e} g) and radius ({pl['radius_cm']:.3e} cm), surface albedo {pl['surface_albedo']}, "
      "solar zenith angle 60 degrees.")
    A(f"- **Atmosphere:** anoxic. N2 fixed at {ex['N2_bar']:g} bar at the surface; CO2 fixed at the run's value (the lever the sweep adjusts); "
      f"water vapour limited to {ex['relative_humidity']*100:.0f}% relative humidity with condensation above the cold trap; hydrogen escape at the diffusion limit.")
    A("- **Surface sources and sinks** (molecules cm^-2 s^-1; `vdep` is a deposition velocity in cm/s): H2 flux 3e10 (volcanic), CH4 flux 3e10 (biological), "
      "CO flux 1e9 with vdep 1e-8, SO2 flux 1e9, H2S flux 1e8, OCS flux 5e7, NO flux 3e8, NH3 flux 1.5e10, H2SO4 flux 7e8, HCN flux 1.7e8, CH3CN flux 1.3e8, "
      "O2 flux 0 (no oxygenic photosynthesis). The full list is in each run's `*_photochem_settings.yaml`.")
    A(f"- **What a run does:** for a given star and instellation, it finds the CO2 partial pressure that holds the surface at {ex['T_target_K']:.0f} K in the "
      "coupled climate (`AdiabatClimate`, radiative-convective equilibrium) and photochemistry (`EvoAtmosphere`) models, iterating until the CO2 pressure "
      "changes by less than about 5%. Climate and photochemistry exchange the surface CH4, H2, CO2 and N2 abundances and the temperature profile.")
    A(f"- **Software:** Photochem {v['photochem']}, Python {v['python']}. The chemistry is a Zahnle-network subset with H, He, N, O, C and S "
      "(`shared/mechanism.yaml`); the climate model uses H2O, CO2, N2, CH4 and H2 (`shared/climate_species.yaml`).")
    A(f"- **Solver settings that matter for a restart:** integrator absolute tolerance {ex['photochem']['atol']:g} (Photochem's default is 1e-23, which is "
      "millions of times slower for cold, CO2-rich cases), maximum steps before giving up raised to 1,000,000, and a persistent pressure-based "
      "temperature and eddy-diffusion profile (see `*_PTKzz_profile.txt`). A run counts as steady when Photochem's own criterion is met or the surface "
      "mixing ratios of CO, O2, H2O, CH4, H2 and CO2 change by less than 1% over the last decade of model time (trace species below 1e-12 are compared in log space).\n")
    A("## File names\n")
    A("A run is identified by a name like `Sun3Ga_S1400.00_CO2-4.651e-05bar`:\n")
    A("| Part | Meaning |\n|---|---|")
    A("| `Sun3Ga` | the Sun as it was 3 billion years ago (`Sun` at age 3 Ga: Photochem's bundled ATLAS-1 solar spectrum, scaled in the past with the Claire et al. 2012 young-Sun routine, with a 5778 K blackbody added at long wavelengths) |")
    A("| `GJ876` | the M dwarf GJ 876: observed MUSCLES spectrum (MAST archive) with a blackbody added at long wavelengths |")
    A("| `S1400.00` | **S**, the total stellar instellation at the planet in W/m^2. The spectrum's shape is fixed and its scale is set so its integral equals this value |")
    A("| `CO2-4.651e-05bar` | the CO2 partial pressure at the surface, in bar, in scientific notation (`6.003e00` means 6.003 bar; `4.651e-05` means 0.00004651 bar), as the lower boundary condition for CO2. For runs where 288 K is not reachable this is the capped or maximum-greenhouse value; for runs with no steady state it is the value that was attempted |\n")
    A("Each run folder `saved_inputs/<run name>/` holds these files, each starting with the run name:\n")
    A("| Suffix | Contents |\n|---|---|")
    A("| `_run.json` | metadata: star, instellation, CO2, outcome flags, the iteration history, surface mixing ratios, parameters, versions |")
    A("| `_photochem_atmosphere.txt` | final photochemical state in Photochem's atmosphere format: `alt` (km), `press` (bar), `den` (molecules/cm^3), `temp` (K), `eddy` (cm^2/s), then the mixing ratio of every particle and gas, then particle radii (`*_r`, cm). Layers run from the surface upward. Mixing ratios are saved to 10 decimals |")
    A("| `_photochem_settings.yaml` | Photochem settings file with this run's CO2 lower boundary (planet, grid, boundary conditions, particles) |")
    A("| `_flux.txt` | stellar flux used (wavelength in nm, flux in mW/m^2/nm) |")
    A("| `_PTKzz_profile.txt` | the pressure-based temperature and eddy-diffusion profile the last photochemical solve ran on: pressure (dyn/cm^2, surface first), temperature (K), Kzz (cm^2/s); the header records the tropopause pressure and target top-of-atmosphere pressure |")
    A("| `_climate_profile.txt` | the final radiative-convective profile from the climate model for the final composition: layer pressure (bar), temperature (K), altitude (km), whether the layer convects with the one below, and the mixing ratio of each climate species; the header gives the surface temperature, surface pressure and tropopause pressure |\n")
    A("`saved_inputs/shared/` holds files common to all runs: `mechanism.yaml`, `climate_species.yaml`, `climate_settings.yaml`, and the downloaded GJ 876 spectrum (`GJ876_muscles_raw.npz`, so the flux can be rebuilt offline).\n")
    dT_hold, dT_unreach = [], []
    for m in metas:
        if m["photochem_stop"].startswith("runaway"): continue
        rid_ = m["run_id"]; T0 = read_profile(f"{saved_dir}/{rid_}/{m['saved_inputs']['files']['profile']}")[1][0]
        m["_T0"] = float(T0); m["_dT"] = abs(T0 - m["T_surf_K"])
        (dT_hold if status(m)[0] == "holds 288 K" else dT_unreach).append(m["_dT"])
    A("The photochemical profile (`_PTKzz_profile.txt`) and the final climate profile (`_climate_profile.txt`) can differ. The photochemistry ran on the profile "
      "computed *before* its final composition was known; the climate profile is then recomputed for the final composition, and its surface temperature is the one "
      f"reported as `T_surf`. For the {len(dT_hold)} runs that hold 288 K the two surface temperatures agree to within {max(dT_hold):.2f} K (median {np.median(dT_hold):.2f} K). "
      f"For the {len(dT_unreach)} steady runs where 288 K is not reachable they differ by up to {max(dT_unreach):.0f} K ({sum(d > 1 for d in dT_unreach)} of them by more than 1 K): the CO2 is pinned at the cap, "
      "so the loop stops changing CO2 and counts as converged while the climate recomputed for the final composition is still moving. Treat the compositions of those runs as indicative only.\n")
    A("## Restarting a run\n")
    example = next((m["run_id"] for m in metas if m["star"] == "Sun" and abs(m["S_Wm2"] - 850) < 1), metas[0]["run_id"])
    A("```python\nfrom sweeplib import reproduce_run, load_climate, read_profile\n\n"
      f"rid = '{example}'          # any folder name in saved_inputs/\n"
      "pc, report = reproduce_run(rid)                   # rebuilds the model from the saved files; converged state in seconds (3 to 81 s across the saved runs)\n"
      "print(report['max_rel_diff_tracked'])             # how closely it matches the saved surface composition\n"
      "atm = pc.mole_fraction_dict()                     # vertical profiles: 'pressure', 'temp', 'CO', 'CH4', ... for plotting\n"
      "clim, ok = load_climate(rid)                      # rebuilds the climate model state (restarts the RCE solve)\n```\n")
    A("Runs flagged `NO STEADY STATE` have nothing to converge to: restarting one simply keeps running away, which `reproduce_run` reports. "
      "Their files are kept for diagnosis only.\n")
    A("`sweeplib.py` is written by the sweep notebook (cell 2 in `archean_instellation_sweep.ipynb`) and is also committed next to it. "
      "The notebook `reproduce_saved_runs.ipynb` shows the whole workflow, including plotting vertical profiles straight from the text files without the model.\n")
    A("## Index of saved runs\n")
    A("| Run | Star | S (W/m2) | CO2 (bar) | T_surf (K) | P_total (bar) | State |\n|---|---|---|---|---|---|---|")
    for m in metas:
        st, _ = status(m)
        if st == "NO STEADY STATE":
            A(f"| `{m['run_id']}` | {STAR_TAG[m['star']]} | {m['S_Wm2']:g} | {m['pCO2_bar']:.3e} | -- | -- | {st} |")
        else:
            A(f"| `{m['run_id']}` | {STAR_TAG[m['star']]} | {m['S_Wm2']:g} | {m['pCO2_bar']:.3e} | {m['T_surf_K']:.1f} | {m['P_total_bar']:.3f} | {st} |")
    A("\n## Each simulation in detail\n")
    for m in metas:
        st, note = status(m); star_name = "the Sun at 3 Ga" if m["star"] == "Sun" else "the M dwarf GJ 876"
        A(f"### `{m['run_id']}`\n")
        A(f"- **Setup:** {star_name}, spectrum scaled to {m['S_Wm2']:g} W/m^2 at the planet; CO2 fixed at {m['pCO2_bar']:.3e} bar, N2 at {m['N2_bar']:g} bar.")
        A(f"- **State:** {st}: {note}")
        it = m["iterations"]
        A(f"- **Coupling:** {it} iteration{'s' if it != 1 else ''}; CO2 pressures tried: " + ", ".join(f"{h['pCO2_bar']:.3e}" for h in m["history"]) + " bar. "
          f"Photochemistry stopped by: {m['photochem_stop']}. Wall time {m['wall_s']/60:.0f} min.")
        if st != "NO STEADY STATE":
            vm = m["vmr"]
            A(f"- **Result:** surface temperature {m['T_surf_K']:.2f} K (the photochemistry ran on a profile whose surface temperature is {m['_T0']:.2f} K), "
              f"total surface pressure {m['P_total_bar']:.3f} bar, tropopause near {m['climate']['P_trop_bar']:.3f} bar. "
              f"Surface mixing ratios: CO2 {sci(vm['CO2'])}, CO {sci(vm['CO'])}, O2 {sci(vm['O2'])}, H2O {sci(vm['H2O'])}, CH4 {sci(vm['CH4'])}, H2 {sci(vm['H2'])}.")
        A(f"- **Files:** `{m['run_id']}/` (see the file-name table above).\n")
    open(f"{saved_dir}/README.md", "w").write("\n".join(lines) + "\n")
    with open(f"{saved_dir}/run_index.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["run_id", "star", "S_Wm2", "pCO2_bar", "T_surf_K", "P_total_bar", "state", "iterations", "photochem_stop"])
        for m in metas:
            st, _ = status(m); ns = st == "NO STEADY STATE"
            w.writerow([m["run_id"], m["star"], m["S_Wm2"], m["pCO2_bar"], "" if ns else round(m["T_surf_K"], 3),
                        "" if ns else round(m["P_total_bar"], 4), st, m["iterations"], m["photochem_stop"]])
    return len(metas)

def run_point(star, S, workdir, results_dir, max_iter=6, dlog_tol=0.02, max_wall_s=10800, log=print,
              saved_dir="saved_inputs"):
    """One sweep point: find the CO2 that gives T_surf = T_target for a coupled
    climate-photochemistry model at instellation S, starting from scratch.
    Result is checkpointed to results_dir/{star}_S{S}.json and returned."""
    os.makedirs(results_dir, exist_ok=True)
    path = f"{results_dir}/{star}_S{S:.1f}.json"
    if is_done(star, S, results_dir, saved_dir):
        return json.load(open(path))
    t_point = time.time()
    tag = f"{star} S={S:.1f}"
    pdir = f"{workdir}/{star}_S{S:.1f}"; os.makedirs(pdir, exist_ok=True)
    for f in ["mechanism.yaml", "photochem_settings.yaml", "climate_species.yaml", "climate_settings.yaml"]:
        if not os.path.exists(f"{pdir}/{f}"): os.symlink(os.path.abspath(f"{workdir}/{f}"), f"{pdir}/{f}")
    flux_file = f"{pdir}/flux.txt"; write_flux_file(star, S, flux_file, workdir)
    pt = Point(pdir, flux_file)
    tracked = TRACKED
    f_CH4, f_H2 = 1.0e-5, 1.0e-3          # first guess for the radiatively active trace gases
    pc = None; history = []
    for it in range(max_iter):
        pCO2, status = pt.solve_pCO2(f_CH4, f_H2)
        T_rce, rce_ok = pt.rce_T(pCO2, f_CH4, f_H2)
        c = pt.clim; P = np.append(c.P_surf, c.P); T = np.append(c.T_surf, c.T)
        P_used, T_used, trop_used = P.copy(), T.copy(), float(c.P_trop)   # profile this iteration's photochemistry runs on
        if pc is None:
            pc = EvoAtmosphere(f"{pdir}/mechanism.yaml", f"{pdir}/photochem_settings.yaml", flux_file)
            pc.var.verbose = 0
            # Absolute tolerance in mixing ratio. The default (1e-23) forces tiny time steps in cold, CO2-rich
            # atmospheres (10^6 times slower to reach steady state); 1e-18 is far below anything reported.
            pc.var.atol = ATOL
            pc.set_lower_bc("CO2", "press", press=pCO2*1e6)
            pc.initialize_atmosphere_p(P, T, eddy_diffusion(P/1e6), mix=None,
                particle_radius={"H2Oaer": np.full(P.size, 1.0e-3)}, persistent=True,
                tropopause_pressure=float(c.P_trop), target_pressure=P[-1])
        else:
            pc.set_lower_bc("CO2", "press", press=pCO2*1e6)
            pc.set_press_temp_edd_profile(P, T, eddy_diffusion(P/1e6), trop_p=float(c.P_trop),
                hydro_pressure=True, maintain_toa_pressure=True, target_pressure=P[-1])
        reason, tn, nsteps, wall, worst = photochem_steady(pc, tracked, max_wall=max_wall_s)
        steady = reason in ("photochem criterion", "composition steady over last decade")
        sm = surface_summary(pc); f_CH4, f_H2 = sm["CH4"], sm["H2"]
        history.append({"pCO2_bar": pCO2, "T_rce_K": T_rce, "photochem_stop": reason, "photochem_t_s": tn, "photochem_wall_s": wall, "limiting_species": worst})
        log(f"[{tag}] iter {it}: pCO2={pCO2:.4g} bar ({status}), T(RCE)={T_rce:.2f} K, photochem: {reason} "
            f"(t={tn:.1e} s, {wall:.0f} s wall, limiting: {worst}), CO={sm['CO']:.2e} CH4={sm['CH4']:.2e} H2={sm['H2']:.2e}")
        if reason.startswith("runaway"):
            pCO2_next = pCO2; break          # no steady state exists, so there is nothing to iterate on
        pCO2_next, status_next = pt.solve_pCO2(f_CH4, f_H2)
        if abs(np.log10(pCO2_next/pCO2)) < dlog_tol: break
    converged = (not reason.startswith("runaway")) and abs(np.log10(pCO2_next/pCO2)) < dlog_tol
    T_final, rce_ok = pt.rce_T(pCO2, f_CH4, f_H2)      # climate for the final composition, at the CO2 used
    sm = surface_summary(pc)
    rid = run_id(star, S, pCO2)
    files = write_saved_inputs(saved_dir, rid, star, S, pCO2, pc, pt, workdir, pdir, flux_file, P_used, T_used, trop_used)
    c = pt.clim
    result = {"run_id": rid, "star": star, "S_Wm2": S, "T_surf_K": float(T_final), "rce_converged": bool(rce_ok),
              "pCO2_bar": float(pCO2), "P_total_bar": float(pc.wrk.surface_pressure), "vmr": sm,
              "target_status": status, "coupling_converged": bool(converged), "iterations": it+1,
              "photochem_stop": reason, "photochem_steady": bool(steady), "wall_s": time.time()-t_point, "history": history,
              "T_target_K": T_target, "relative_humidity": relative_humidity, "N2_bar": pN2_bar,
              "planet": {"mass_g": planet_mass, "radius_cm": planet_radius, "surface_albedo": surface_albedo},
              "photochem": {"atol": ATOL, "nsteps_before_giveup": 1000000, "trop_p_dyn_cm2": trop_used,
                            "target_toa_pressure_dyn_cm2": float(P_used[-1]), "tracked_species": list(TRACKED)},
              "climate": {"species": list(c.species_names),
                          "surface_partial_pressures_dyn_cm2": dict(zip(c.species_names, map(float, pt.partial_pressures(pCO2, f_CH4, f_H2)))),
                          "T_surf_K": float(c.T_surf), "P_surf_bar": float(c.P_surf / 1e6), "P_trop_bar": float(c.P_trop / 1e6)},
              "versions": {"photochem": __import__("photochem").__version__, "python": sys.version.split()[0]},
              "saved_inputs": {"directory": rid, "files": files}}
    json.dump(result, open(path, "w"), indent=1)
    json.dump(result, open(f"{saved_dir}/{rid}/{rid}_run.json", "w"), indent=1)
    log(f"[{tag}] DONE: T={T_final:.2f} K pCO2={pCO2:.4g} bar Ptot={pc.wrk.surface_pressure:.3f} bar "
        f"converged={converged} ({time.time()-t_point:.0f} s)")
    return result

def worker(args):
    "Run one (star, S) point in a worker process; progress is appended to <workdir>/sweep.log."
    os.environ["OMP_NUM_THREADS"] = "1"
    star, S, workdir, results_dir = args[:4]
    saved_dir = args[4] if len(args) > 4 else "saved_inputs"
    def log(msg):
        print(msg, flush=True)
        with open(f"{workdir}/sweep.log", "a") as f: f.write(msg + "\n")
    return run_point(star, S, workdir, results_dir, log=log, saved_dir=saved_dir)
