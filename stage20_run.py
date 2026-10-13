"""
Stage 20 runner -- film transport (F) x wetting rule (W), 2x2 factorial.
TWO SOILS ONLY: Sand and Vertisol (docs/stage_results/STAGE20_RESULTS.md).

Everything except the two factors is Stage 19's, imported from stage19_run
(same mapping cleaning, densities, coating fractions, liveness, assertions):
  F0  constant film D_film x D_scale = 4.70 in drained macropores (Stage 19)
  F1  D_film(i,j) = D(d_ij) x D_scale x S_w^(10/3)  (Millington-Quirk)
  W0  matric_regions, d_cut 50 um (Stage 19)
  W1  equal saturation: 50% of non-solid cells wet, smallest per-cell d first
Arms: A0 = F0+W0 (Stage 19, reused; one point verified byte-identical),
      A1 = F1+W0 (candidate), A2 = F0+W1 (diagnostic only), A3 = F1+W1 (candidate).

    python stage20_run.py step0|step1|step2|step3
"""
from __future__ import annotations

import copy
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import dual_porosity
import stage19_run as s19
from model import load_yaml

SOILS = s19.SOILS
RESULTS_DIR = os.path.join("results", "stage20")
OM_DENSITIES = s19.OM_DENSITIES
K_DENSITY = 0.1095862986462395       # Stage 19's rule value, kept (R5)
W1_SW = 0.5
ARMS = {"A0": ("F0", "W0"), "A1": ("F1", "W0"), "A2": ("F0", "W1"), "A3": ("F1", "W1")}
CANDIDATES = ("A1", "A3")
EQUAL_CF = 0.5


def structure_cfg(soil, n, wetting, coating_fraction=None, seed_offset=0):
    cfg = s19.structure_cfg(soil, n, coating_fraction=coating_fraction, seed_offset=seed_offset)
    if wetting == "W1":
        cfg["saturation"] = {"mode": "equal_saturation", "S_w": W1_SW}
    return cfg


def mapping_for(om_density, film, D_scale=None):
    m = s19.mapping_for(om_density, K_DENSITY, D_scale=D_scale)
    if film == "F1":
        m["stage20_film_mode"] = "millington_quirk"
    return m


def spec(key, soil, n, arm, om, coating_fraction=None, D_scale=None, seed_offset=0, shuffle=False):
    film, wet = ARMS[arm]
    cfg = structure_cfg(soil, n, wet, coating_fraction, seed_offset)
    m = mapping_for(om, film, D_scale)
    if shuffle:
        cfg["pore"]["shuffle_seed"] = 9000 + seed_offset
        m["stage20_conserve_om_pools"] = True   # see dual_porosity / STAGE20_RESULTS
    return {"key": key, "soil": soil, "cfg": cfg, "m": m}


def secondary_alive(r):
    return bool(r["R_peak_over_floor"] >= 2.0 and r["growth_share_at_peak"] > 0.5)


def save(name, obj):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(os.path.join(RESULTS_DIR, name), "w") as f:
        json.dump(obj, f, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))


def load(name):
    with open(os.path.join(RESULTS_DIR, name)) as f:
        return json.load(f)


# ------------------------------------------------------------------ step 0

def step0():
    s19.check_biology()
    m_base = mapping_for(0.1, "F0")
    D_scale = float(m_base["stage7_D_scale"])
    k_cl = float(m_base.get("stage4_variantC_cluster_K_scale", 0.75))
    d_cl = float(m_base.get("stage4_variantC_cluster_D_scale", 0.75))
    rows = []
    for s in SOILS:
        for n in (40, 128):
            for wet in ("W0", "W1"):
                g = dual_porosity.build_dual_grids(structure_cfg(s, n, wet), m_base)
                macro, water = g["macro_mask"], g["water_mask"]
                dry = macro & ~water
                wetm = macro & water
                S_w = float(g["S_w"])
                f_rel = S_w ** (10.0 / 3.0)
                Dbase = g["D_full"] * D_scale
                _, _, scale = dual_porosity.habitat_properties(macro)
                D_wet = Dbase * (1.0 + d_cl * scale)
                F0 = float(m_base["D_film"]) * D_scale
                F1 = Dbase * f_rel
                K_t = np.where(macro, g["K"] * (1.0 + k_cl * scale), 0.0)
                med = lambda a: float(np.median(a)) if a.size else float("nan")
                mean = lambda a: float(np.mean(a)) if a.size else float("nan")
                rows.append({
                    "soil": s, "n": n, "wetting": wet, "S_w": S_w, "f_rel": f_rel,
                    "wet_frac_habitat": float(wetm.sum() / max(macro.sum(), 1)),
                    "wet_frac_matrix": float((g["porous_matrix_mask"] & water).sum()
                                             / max(g["porous_matrix_mask"].sum(), 1)),
                    "drained_D_F0": F0,
                    "drained_D_F1_median": med(F1[dry]), "drained_D_F1_mean": mean(F1[dry]),
                    "wet_macro_D_median": med(D_wet[wetm]), "wet_macro_D_mean": mean(D_wet[wetm]),
                    "ratio_drained_over_wet_F0": F0 / med(D_wet[wetm]) if wetm.any() else float("nan"),
                    "ratio_drained_over_wet_F1": med(F1[dry]) / med(D_wet[wetm]) if wetm.any() and dry.any() else float("nan"),
                    "coating_cells": int(g["coating_mask"].sum()), "trapped_cells": int(g["trapped_mask"].sum()),
                    "native_K_t_density": float(K_t.sum() / K_t.size),
                })
    save("stage20_step0_budget.json", {"rows": rows, "K_density_kept": K_DENSITY,
                                       "D_film_F0": float(m_base["D_film"]) * D_scale})
    for r in rows:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})


# ------------------------------------------------------------------ step 1

def step1():
    s19.check_biology()
    specs = [spec(f"{arm}|{x:.4f}", s, 40, arm, x) for arm in ("A1", "A2", "A3") for x in OM_DENSITIES for s in SOILS]
    res = s19.run_many(specs)
    rows, curves = [], {}
    for arm in ("A1", "A2", "A3"):
        for x in OM_DENSITIES:
            rs = [r for r in res if r["key"] == f"{arm}|{x:.4f}"]
            p = s19.pair(rs)
            curves[(arm, x)] = {r["soil"]: r["R_t"] for r in rs}
            for r in rs:
                rows.append(dict(s19.flat(r), arm=arm, OM_density=x, distance=p["distance"],
                                 alive_primary=r["alive"], alive_secondary=secondary_alive(r),
                                 both_alive_primary=p["both_alive"]))
    # A0 from Stage 19's map (verified byte-identical at one point)
    a0 = load_s19_map()
    for r in a0:
        rows.append(dict(r, arm="A0", alive_primary=r["alive"], alive_secondary=secondary_alive(r),
                         both_alive_primary=r["both_alive"]))
    # control (a) at every candidate density with distance > 0.3
    ctrl_specs = []
    for arm in CANDIDATES:
        for x in OM_DENSITIES:
            d = [r["distance"] for r in rows if r["arm"] == arm and r["OM_density"] == x][0]
            if d > 0.3:
                ctrl_specs += [spec(f"ctrlA|{arm}|{x:.4f}", s, 40, arm, x, coating_fraction=EQUAL_CF) for s in SOILS]
    ctrl = {}
    if ctrl_specs:
        cres = s19.run_many(ctrl_specs)
        for k in sorted({r["key"] for r in cres}):
            rs = [r for r in cres if r["key"] == k]
            p = s19.pair(rs)
            _, arm, x = k.split("|")
            ctrl[f"{arm}|{x}"] = {"distance": p["distance"], "both_alive": p["both_alive"],
                                  "soils": {r["soil"]: s19.flat(r) for r in rs}}
    # selection per candidate arm
    selected = {}
    for arm in CANDIDATES:
        ok = [x for x in OM_DENSITIES
              if all(r["alive_primary"] for r in rows if r["arm"] == arm and r["OM_density"] == x)]
        if ok:
            best = max(ok, key=lambda x: [r["distance"] for r in rows if r["arm"] == arm and r["OM_density"] == x][0])
            selected[arm] = best
    save("stage20_step1_map.json", {"rows": rows, "control_a": ctrl, "selected": selected})
    s19_csv(rows)
    plot_map(curves)
    for arm in ("A0", "A1", "A2", "A3"):
        for x in OM_DENSITIES:
            for s in SOILS:
                r = [q for q in rows if q["arm"] == arm and abs(q["OM_density"] - x) < 1e-9 and q["soil"] == s][0]
                print(f"{arm} OM={x:.4f} {s:8s} P={r['alive_primary']!s:5} S={r['alive_secondary']!s:5} "
                      f"B/B0={r['B_peak_over_B0']:.3f} R/fl={r['R_peak_over_floor']:.2f} gs={r['growth_share_at_peak']:.2f} "
                      f"pk@{r['peak_frac']:.2f} dry_hab={r['dry_habitat_frac']:.3f} "
                      f"del_wet={r.get('delivered_wet_total', float('nan')):.1f} del_film={r.get('delivered_dry_total', float('nan')):.1f} "
                      f"dist={r['distance']:.3f}")
    print("control (a):", {k: round(v["distance"], 3) for k, v in ctrl.items()})
    print("selected:", selected)


def load_s19_map():
    with open(os.path.join("results", "stage19", "stage19_step1_map.json")) as f:
        return json.load(f)["rows"]


def s19_csv(rows):
    import csv
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(RESULTS_DIR, "stage20_step1_map.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})


def plot_map(curves):
    arms = ("A1", "A2", "A3")
    fig, axes = plt.subplots(len(arms), len(OM_DENSITIES), figsize=(4 * len(OM_DENSITIES), 3.3 * len(arms)),
                             squeeze=False)
    for i, arm in enumerate(arms):
        for j, x in enumerate(OM_DENSITIES):
            ax = axes[i, j]
            for s, R in curves[(arm, x)].items():
                ax.plot(np.arange(len(R)) * 0.05, R / R.max(), color=s19.COLORS[s], label=f"{s} pk={R.max():.2f}")
            ax.set_title(f"{arm} {ARMS[arm]} OM={x:.3f} (R/peak)", fontsize=8)
            ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, "stage20_step1_map.png"), dpi=110)
    plt.close(fig)


# ------------------------------------------------------------------ steps 2-3

def step2():
    s19.check_biology()
    sel = load("stage20_step1_map.json")["selected"]
    assert sel, "no candidate arm has a density with both soils alive -- stop"
    out = {}
    for arm, x in sel.items():
        rs = s19.run_many([spec(f"{arm}", s, 128, arm, x) for s in SOILS])
        p = s19.pair(rs)
        out[arm] = {"OM_density": x, "distance": p["distance"], "both_alive": p["both_alive"], "pass": p["pass"],
                    "soils": {s: s19.flat(r) for s, r in p["by"].items()},
                    "curves": {s: r["R_t"] for s, r in p["by"].items()}}
        print(arm, x, "distance", round(p["distance"], 3), "both alive", p["both_alive"],
              {s: (r["alive"], round(r["B_peak_over_B0"], 3), round(r["R_peak_over_floor"], 2)) for s, r in p["by"].items()})
    save("stage20_step2_n128.json", out)


def step3():
    s19.check_biology()
    s2 = load("stage20_step2_n128.json")
    out = {}
    for arm, base in s2.items():
        x = base["OM_density"]
        specs = []
        for s in SOILS:
            specs.append(spec("a_equal_cf", s, 128, arm, x, coating_fraction=EQUAL_CF))
            specs.append(spec("b_transport_off", s, 128, arm, x, D_scale=1e-9))
        res = s19.run_many(specs)
        o = {"baseline": {"distance": base["distance"], "both_alive": base["both_alive"]}}
        for key in ("a_equal_cf", "b_transport_off"):
            p = s19.pair([r for r in res if r["key"] == key])
            o[key] = {"distance": p["distance"], "both_alive": p["both_alive"],
                      "soils": {s: s19.flat(r) for s, r in p["by"].items()},
                      "curves": {s: p["by"][s]["R_t"] for s in SOILS}}
        o["c_delivery"] = {s: {k: base["soils"][s][k] for k in ("dry_habitat_frac", "delivered_wet_total",
                                                                 "delivered_dry_total")} for s in SOILS}
        out[arm] = o
        print(arm, {k: (round(v["distance"], 3), v["both_alive"]) for k, v in o.items() if k != "c_delivery"},
              o["c_delivery"])
    save("stage20_step3_controls.json", out)


# ------------------------------------------------------------------ step 4

# Pre-registered Step 4 criteria (fixed before the run):
#   S3a  all 3 seed sets: both soils PRIMARY-alive and distance > 0.3
#   S3b  mean within-soil seed-to-seed distance < 0.5 x min between-soil distance
#   S4   shuffling the pore field (same PSD, arrangement destroyed) lowers the
#        distance by >= 0.1 (Stage 17's arrangement criterion)
#   MASS transport step (implicit diffusion + leak, the model's own solver
#        settings, on the configuration's real D field) conserves total
#        substrate to relative error < 1e-9
S4_MIN_DROP = 0.1
MASS_TOL = 1e-9


def mass_check(arm, om, n=128):
    """One transport step on the real grid of each soil: random S, the
    configuration's own D field (wet macropores with cluster bonus, film,
    matrix) and leak, with the model's iteration cap and tolerance."""
    film, wet = ARMS[arm]
    m = mapping_for(om, film)
    out = {}
    rng = np.random.default_rng(0)
    for s in SOILS:
        g = dual_porosity.build_dual_grids(structure_cfg(s, n, wet), m)
        D_scale = float(m["stage7_D_scale"])
        macro, water = g["macro_mask"], g["water_mask"]
        _, _, scale = dual_porosity.habitat_properties(macro)
        Dbase = g["D_full"] * D_scale
        f_rel = float(g["S_w"]) ** (10.0 / 3.0)
        D_film = Dbase * f_rel if film == "F1" else float(m["D_film"]) * D_scale
        D_macro = np.where(water, Dbase * (1.0 + float(m["stage4_variantC_cluster_D_scale"]) * scale), 0.0)
        D_macro = np.where(macro & ~water, D_film, D_macro)
        D_matrix = np.where(water & g["porous_matrix_mask"], float(m["stage4_variantC_matrix_diffusion"]) * D_scale, 0.0)
        D = D_macro + D_matrix
        S = rng.uniform(0.0, 5.0, D.shape)
        S_d = dual_porosity.implicit_diffusion_step(S, D, 1.0, 0.05, int(m["stage6_implicit_iterations"]),
                                                    float(m["stage6_implicit_relaxation"]),
                                                    tol=float(m["stage11_implicit_tol"]))
        S_new = S_d + 0.05 * dual_porosity.leak_exchange(S_d, g["porous_matrix_mask"], macro, float(m["k_leak"]))
        out[s] = float(abs(S_new.sum() - S.sum()) / S.sum())
    return out


def step4():
    import metrics
    s19.check_biology()
    s3 = load("stage20_step3_controls.json")
    s2 = load("stage20_step2_n128.json")
    passing = [a for a, v in s2.items() if v["pass"]]
    assert passing, "no configuration passed Step 3 -- Step 4 not run"
    out = {}
    for arm in passing:
        x = s2[arm]["OM_density"]
        res = {}
        for n in (128, 40):
            specs = [spec(f"seed{so}", s, n, arm, x, seed_offset=so) for so in (0, 10, 20) for s in SOILS]
            specs += [spec("shuffle", s, n, arm, x, shuffle=True) for s in SOILS]
            rs = s19.run_many(specs)
            sets = {so: s19.pair([r for r in rs if r["key"] == f"seed{so}"]) for so in (0, 10, 20)}
            shuf = s19.pair([r for r in rs if r["key"] == "shuffle"])
            within = [float(metrics.pairwise_distance(sets[a]["by"][s]["R_t"], sets[b]["by"][s]["R_t"]))
                      for s in SOILS for a, b in ((0, 10), (0, 20), (10, 20))]
            min_between = min(v["distance"] for v in sets.values())
            res[n] = {
                "per_seed_set": {so: {"distance": v["distance"], "both_alive": v["both_alive"],
                                      "soils": {s: {k: v["by"][s][k] for k in ("alive", "B_peak_over_B0",
                                               "R_peak_over_floor", "peak_frac", "R_end_over_peak")} for s in SOILS}}
                                 for so, v in sets.items()},
                "S3a": all(v["pass"] for v in sets.values()),
                "within_mean": float(np.mean(within)), "min_between": min_between,
                "S3b": float(np.mean(within)) < 0.5 * min_between,
                "shuffle": {"distance": shuf["distance"], "both_alive": shuf["both_alive"],
                            "soils": {s: {k: shuf["by"][s][k] for k in ("alive", "B_peak_over_B0",
                                     "R_peak_over_floor", "dry_habitat_frac", "peak_frac",
                                     "coating_fraction")} for s in SOILS}},
                "S4_drop": sets[0]["distance"] - shuf["distance"],
                "S4": (sets[0]["distance"] - shuf["distance"]) >= S4_MIN_DROP,
            }
        mc = mass_check(arm, x)
        out[arm] = {"OM_density": x, "grids": res, "mass_rel_error": mc,
                    "MASS": all(v < MASS_TOL for v in mc.values())}
        print(json.dumps(out[arm], indent=1, default=str))
    save("stage20_step4_robustness.json", out)


if __name__ == "__main__":
    {"step0": step0, "step1": step1, "step2": step2, "step3": step3, "step4": step4}[sys.argv[1]]()
