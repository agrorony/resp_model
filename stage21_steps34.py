"""Stage 21 steps 3 and 4 (imports stage21_run; python stage21_steps34.py step3|step4)."""
from __future__ import annotations

import itertools
import json
import os
import sys
import time

import numpy as np

import packing3d as P
from stage21_run import (RESULTS_DIR, THRESHOLD, dist, get_structure, jkey, load, run_jobs, save)


def _alive(r):
    return {"alive_primary": r["alive_primary"], "alive_secondary": r["alive_secondary"],
            "B_peak_over_B0": r["B_peak_over_B0"], "R_peak_over_floor": r["R_peak_over_floor"],
            "growth_share_at_peak": r["growth_share_at_peak"], "peak_frac": r["peak_frac"],
            "R_end_over_peak": r["R_end_over_peak"],
            "diffusion_into_wet_habitat": r["diffusion_into_wet_habitat"],
            "diffusion_into_drained_habitat": r["diffusion_into_drained_habitat"],
            "exchange_into_wet_habitat": r["exchange_into_wet_habitat"],
            "exchange_into_drained_habitat": r["exchange_into_drained_habitat"],
            "mass_balance_rel": r["mass_balance_rel"]}


def step3():
    s2 = load("step2.json")
    om = s2["selected_OM_density"]
    assert om is not None, "Step 2 selected no density"
    xp = P.get_xp()
    for k in ("S", "V", "Sp", "Vp", "Val"):
        get_structure(k, 0, xp=xp)
    keys = {"S": jkey("S", 0, om), "V": jkey("V", 0, om), "Sp": jkey("Sp", 0, om),
            "Vp": jkey("Vp", 0, om), "Val": jkey("Val", 0, om),
            "S_off": jkey("S", 0, om, dscale=1e-9), "V_off": jkey("V", 0, om, dscale=1e-9)}
    runs = dict(zip(keys, run_jobs(list(keys.values())).values()))
    d = lambda a, b: dist(runs[a], runs[b])
    res = {"OM_density": om, "runs": {k: _alive(r) for k, r in runs.items()}}
    res["a_factorial"] = {
        "S_vs_V (total)": d("S", "V"),
        "shape effect, sand material: S vs Sp": d("S", "Sp"),
        "shape effect, vertisol material: Vp vs V": d("Vp", "V"),
        "material effect, sphere shape: S vs Vp": d("S", "Vp"),
        "material effect, plate shape: Sp vs V": d("Sp", "V"),
        "Sp_vs_Vp (both off-diagonal)": d("Sp", "Vp"),
    }
    res["b_orientation"] = {"V_vs_Val": d("V", "Val"), "S_vs_Val": d("S", "Val")}
    res["c_transport_off"] = {
        "S_vs_V_off": d("S_off", "V_off"), "S_vs_S_off": d("S", "S_off"),
        "V_vs_V_off": d("V", "V_off"),
        "Rpeak_ratio_S": runs["S_off"]["R_peak_over_floor"] / runs["S"]["R_peak_over_floor"],
        "Rpeak_ratio_V": runs["V_off"]["R_peak_over_floor"] / runs["V"]["R_peak_over_floor"],
    }
    res["distinctness"] = d("S", "V")
    res["both_alive_primary"] = bool(runs["S"]["alive_primary"] and runs["V"]["alive_primary"])
    res["pass"] = bool(res["both_alive_primary"] and res["distinctness"] > THRESHOLD)
    save("step3.json", res)
    plot3(runs)
    return res


def plot3(runs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    for k, c in (("S", "#d62728"), ("Sp", "#ff9896"), ("V", "#4b0082"), ("Vp", "#9467bd"), ("Val", "#17becf")):
        r = runs[k]; ax[0].plot(r["R"] / (0.02 * r["B"][0]), c=c, label=k)
    ax[0].set_title("factorial / orientation"); ax[0].legend(); ax[0].set_ylabel("R/(m0 B0)")
    for k, c, ls in (("S", "#d62728", "-"), ("S_off", "#d62728", ":"), ("V", "#4b0082", "-"), ("V_off", "#4b0082", ":")):
        r = runs[k]; ax[1].plot(r["R"] / (0.02 * r["B"][0]), c=c, ls=ls, label=k)
    ax[1].set_title("transport off (D_scale 1e-9)"); ax[1].legend()
    fig.tight_layout(); fig.savefig(os.path.join(RESULTS_DIR, "step3_curves.png"), dpi=110); plt.close(fig)


def step4():
    s3 = load("step3.json")
    assert s3["pass"], "Step 3 did not pass"
    om = s3["OM_density"]
    xp = P.get_xp()
    seeds = (0, 1, 2)
    for k in ("S", "V"):
        for s in seeds:
            get_structure(k, s, xp=xp)
    keys = {(k, s): jkey(k, s, om) for k in ("S", "V") for s in seeds}
    runs = dict(zip(keys, run_jobs(list(keys.values())).values()))
    within = {k: [dist(runs[(k, a)], runs[(k, b)]) for a, b in itertools.combinations(seeds, 2)]
              for k in ("S", "V")}
    between = [dist(runs[("S", a)], runs[("V", b)]) for a in seeds for b in seeds]
    same_seed = [dist(runs[("S", s)], runs[("V", s)]) for s in seeds]
    noise = max(max(v) for v in within.values())
    res = {"OM_density": om, "within_soil_seed_distance": within,
           "between_soil_distances_all_pairs": between, "same_seed_pair_distance": same_seed,
           "seed_noise_max": noise, "between_min": min(between),
           "noise_lt_half_between": bool(noise < 0.5 * min(between)),
           "all_seeds_alive_primary": bool(all(r["alive_primary"] for r in runs.values())),
           "all_seed_pairs_distinct": bool(min(between) > THRESHOLD),
           "mass_balance_max": max(r["mass_balance_rel"] for r in runs.values()),
           "runs": {f"{k}_seed{s}": _alive(r) for (k, s), r in runs.items()}}
    res["pass"] = bool(res["noise_lt_half_between"] and res["all_seeds_alive_primary"]
                       and res["all_seed_pairs_distinct"] and res["mass_balance_max"] < 1e-9)
    save("step4.json", res)
    return res


if __name__ == "__main__":
    t0 = time.time()
    r = globals()[sys.argv[1]]()
    print(json.dumps(r, indent=1, default=str)[:5000])
    print(f"[{sys.argv[1]}] {time.time() - t0:.1f} s")
