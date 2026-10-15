"""
Stage 22 runner: equal habitat-scale wetting (50% of resolved pores wet in both soils).
Sand and Vertisol ONLY; no Loess; no measured data. Only the wetting rule differs from Stage 21.

    python stage22_run.py step0     wetting map + OFF-switch reproduction of Stage 21
    python stage22_run.py step1     S_eq, V_eq at OM 1.0 seed 0 + distances + figure
    python stage22_run.py step2     seeds + shape factorial (only if step1 passes)
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time

import numpy as np

import model3d as M
import packing3d as P
import stage21_run as S21
from model import load_yaml

OUT = os.path.join("results", "stage22")
OM = 1.0
SIGMA = 0.072  # N/m, water-air surface tension


def save(name, obj):
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, name), "w") as f:
        json.dump(obj, f, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))


def load(name):
    with open(os.path.join(OUT, name)) as f:
        return json.load(f)


def psi_kpa(d_um):
    return -4 * SIGMA / (d_um * 1e-6) / 1e3


def wetting_map(st, hp, xp):
    pore = st.phase == P.PORE
    lt = st.lt.astype(xp.float64)
    npore = float(pore.sum())
    Kn_hab = M.native_K(st, hp, xp) > 0
    out = {}
    for rule in ("matric", "equal50"):
        wet = (pore & (lt < hp["d_cut_um"])) if rule == "matric" else M.equal_wet_mask(st, xp)
        matrix = st.phase == P.MATRIX
        solid = st.phase == P.SOLID
        Sw_res = float(wet.sum()) / npore
        Sw_vox = float(wet.sum() + matrix.sum()) / float((~solid).sum())
        d_air = float(lt[wet].max()) if bool(wet.any()) else 0.0
        out[rule] = {
            "wet_fraction_resolved_pores": Sw_res,
            "wet_fraction_habitat": float((wet & Kn_hab).sum()) / float(Kn_hab.sum()),
            "habitat_volume_fraction": float(Kn_hab.sum()) / float(st.phase.size),
            "S_w_voxel": Sw_vox,
            "S_w_used_for_film": Sw_vox if rule == "matric" else Sw_res,
            "f_rel": (Sw_vox if rule == "matric" else Sw_res) ** (10.0 / 3.0),
            "air_entry_diameter_um": d_air,
            "air_entry_psi_kPa": psi_kpa(d_air) if d_air > 0 else None,
            "LT_median_um": float(xp.median(lt[pore])),
        }
    return out


# ------------------------------------------------------------- jobs
def jk(kind, seed, rule):
    return f"{kind}_s{seed}_{rule}"


def job(key):
    S21.check_biology()
    kind, s, rule = key.split("_")
    seed = int(s[1:])
    xp = P.get_xp()
    hp = M.load_inherited()
    bio = load_yaml("configs/biology.yaml")
    K_density = S21.load("step1.json")["K_density"]  # Stage 21 value, kept fixed
    st = S21.get_structure(kind, seed, S21.N_FULL, S21.VOX_FULL, xp)
    wetting = {"eq": "equal50", "off": "matric"}[rule]
    t0 = time.time()
    out, res = M.run_structure(st, hp, bio, xp, K_density, OM, wetting=wetting)
    res.update({"key": key, "seconds": time.time() - t0, "K_density": K_density, "wetting": wetting})
    os.makedirs(os.path.join(OUT, "runs"), exist_ok=True)
    np.savez_compressed(os.path.join(OUT, "runs", key + ".npz"), **{k: out[k] for k in M.STAT_NAMES})
    with open(os.path.join(OUT, "runs", key + ".json"), "w") as f:
        json.dump(res, f, indent=1)


def run_jobs(keys, ngpu=2):
    pending = [k for k in dict.fromkeys(keys) if not os.path.exists(os.path.join(OUT, "runs", k + ".json"))]
    running = {}
    while pending or running:
        for slot in range(ngpu):
            if slot not in running and pending:
                k = pending.pop(0)
                env = dict(os.environ, CUDA_VISIBLE_DEVICES=str(slot), PYTHONDONTWRITEBYTECODE="1")
                running[slot] = (k, subprocess.Popen([sys.executable, "stage22_run.py", "job", k], env=env,
                                                     stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True))
                print("start", k, "gpu", slot, flush=True)
        time.sleep(2)
        for slot, (k, pr) in list(running.items()):
            if pr.poll() is not None:
                if pr.returncode != 0:
                    raise RuntimeError(f"job {k} failed: {pr.stderr.read()[-2000:]}")
                print("done", k, flush=True)
                del running[slot]


def load_run(key):
    with open(os.path.join(OUT, "runs", key + ".json")) as f:
        r = json.load(f)
    z = np.load(os.path.join(OUT, "runs", key + ".npz"))
    r["R"], r["B"] = z["R"], z["B"]
    return r


def load21(kind, seed=0):
    return S21.load_run(S21.jkey(kind, seed, OM))


def dist(a, b):
    return S21.dist(a, b)


# ------------------------------------------------------------- steps
def step0():
    t0 = time.time()
    xp = P.get_xp()
    hp = M.load_inherited()
    res = {"maps": {}}
    st21 = S21.load("step1.json")
    for k in ("S", "V"):
        st = S21.get_structure(k, 0, S21.N_FULL, S21.VOX_FULL, xp)
        res["maps"][k] = wetting_map(st, hp, xp)
        res["maps"][k]["native_K_density_equal50"] = float(M.native_K(st, hp, xp, wetting="equal50").sum()) / st.phase.size
        res["maps"][k]["native_K_density_stage21"] = st21["native_K_density"][k]
    res["K_density_used"] = st21["K_density"]
    # OFF-switch reproduction: rerun Stage 21 S and V with the default rule on the new code
    run_jobs([jk("S", 0, "off"), jk("V", 0, "off")])
    rep = {}
    for k in ("S", "V"):
        a, b = load_run(jk(k, 0, "off")), load21(k)
        rep[k] = {"R_rel_err_max": float(np.max(np.abs(a["R"] - b["R"])) / np.max(np.abs(b["R"]))),
                  "B_rel_err_max": float(np.max(np.abs(a["B"] - b["B"])) / np.max(np.abs(b["B"])))}
    res["off_switch_reproduction"] = rep
    res["seconds"] = time.time() - t0
    save("step0.json", res)
    print(json.dumps(res, indent=1))


KEYS = ("alive_primary", "alive_secondary", "B_peak_over_B0", "R_peak_over_floor", "growth_share_at_peak",
        "peak_frac", "R_end_over_peak", "fuel_remaining_final", "diffusion_into_wet_habitat",
        "diffusion_into_drained_habitat", "exchange_into_wet_habitat", "exchange_into_drained_habitat",
        "mass_balance_rel", "seconds")


def plot1(Se, Ve, S0, V0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for a, norm in zip(ax, (False, True)):
        for r, lab, c, ls in ((Se, "S_eq", "C0", "-"), (Ve, "V_eq", "C3", "-"),
                              (S0, "S_21", "C0", "--"), (V0, "V_21", "C3", "--")):
            R = np.asarray(r["R"], float)
            t = np.linspace(0, 1, len(R))
            a.plot(t, R / R.max() if norm else R, ls, color=c, label=lab)
        a.set_xlabel("t / T"); a.set_ylabel("R / own peak" if norm else "R"); a.legend()
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "step1_curves.png"), dpi=110); plt.close(fig)


def step1():
    t0 = time.time()
    run_jobs([jk("S", 0, "eq"), jk("V", 0, "eq")])
    Se, Ve = load_run("S_s0_eq"), load_run("V_s0_eq")
    S0, V0 = load21("S"), load21("V")
    d = {"S_eq_vs_V_eq": dist(Se, Ve), "S_eq_vs_S_21": dist(Se, S0), "V_eq_vs_V_21": dist(Ve, V0),
         "S_eq_vs_V_21": dist(Se, V0), "V_eq_vs_S_21": dist(Ve, S0), "S_21_vs_V_21": dist(S0, V0)}
    both = bool(Se["alive_primary"] and Ve["alive_primary"])
    reading = ("collapsed: Stage 21 distinction carried by the wetting asymmetry" if d["S_eq_vs_V_eq"] < 0.3
               else ("material carries distinction under equal wetting" if both else "distinct but a soil is not alive"))
    res = {"runs": {k: {x: r[x] for x in KEYS} for k, r in (("S_eq", Se), ("V_eq", Ve), ("S_21", S0), ("V_21", V0))},
           "distances": d, "both_primary_alive": both, "reading": reading,
           "proceed_step2": bool(d["S_eq_vs_V_eq"] >= 0.3 and both), "seconds": time.time() - t0}
    save("step1.json", res)
    plot1(Se, Ve, S0, V0)
    print(json.dumps(res, indent=1))


def step2():
    t0 = time.time()
    keys = [jk(k, s, "eq") for k in ("S", "V") for s in (1, 2)] + [jk("Sp", 0, "eq"), jk("Vp", 0, "eq")]
    run_jobs(keys)
    R = {(k, s): load_run(jk(k, s, "eq")) for k in ("S", "V") for s in (1, 2)}
    R[("S", 0)], R[("V", 0)] = load_run("S_s0_eq"), load_run("V_s0_eq")
    Sp, Vp = load_run("Sp_s0_eq"), load_run("Vp_s0_eq")
    within = {k: [dist(R[(k, a)], R[(k, b)]) for a, b in ((0, 1), (0, 2), (1, 2))] for k in ("S", "V")}
    between = [dist(R[("S", a)], R[("V", b)]) for a in range(3) for b in range(3)]
    noise = max(max(v) for v in within.values())
    alive = all(r["alive_primary"] for r in R.values())
    a = {"within_soil": within, "between_all": between, "noise_max": noise, "between_min": min(between),
         "noise_lt_half_between": bool(noise < 0.5 * min(between)), "all_alive_primary": alive,
         "mass_balance_max": max(r["mass_balance_rel"] for r in R.values()),
         "per_seed": {f"{k}{s}": {x: R[(k, s)][x] for x in KEYS} for (k, s) in R}}
    Se, Ve = R[("S", 0)], R[("V", 0)]
    b = {"runs": {k: {x: r[x] for x in KEYS} for k, r in (("Sp_eq", Sp), ("Vp_eq", Vp))},
         "S_vs_V": dist(Se, Ve), "shape_sand_S_vs_Sp": dist(Se, Sp), "shape_vertisol_Vp_vs_V": dist(Vp, Ve),
         "material_sphere_S_vs_Vp": dist(Se, Vp), "material_plate_Sp_vs_V": dist(Sp, Ve), "Sp_vs_Vp": dist(Sp, Vp)}
    save("step2.json", {"a_robustness": a, "b_factorial": b, "seconds": time.time() - t0})
    print(json.dumps({"a": {k: v for k, v in a.items() if k != "per_seed"}, "b": b}, indent=1))


if __name__ == "__main__":
    globals()[sys.argv[1]](*sys.argv[2:])
