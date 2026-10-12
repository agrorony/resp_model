"""
Stage 19 runner -- Stage 8 revived (equal OM + per-soil coating fraction),
made alive. TWO SOILS ONLY: Sand and Vertisol (docs/stage_results/
STAGE19_RESULTS.md).

Rules (fixed in this file before any run):
  R1 Sand + Vertisol only.             R2 literature PSDs (psd_parametric).
  R3 biology frozen (SHA-256); ONE dissolution rate (k_dis_slow = k_dis);
     no enzyme; om_fine_bias = 0 (uniform OM within each pool).
  R4 om_mode coating_trapped, coating_thickness 2; coating_fraction is the
     only per-soil OM parameter (configs/soil_<s>_stage19.yaml: 0.8 / 0.3).
  R5 equal DENSITIES per cell, asserted: sum(OM)=OM_density*N,
     sum(B0)=0.0125*N, sum(K_t)=K_density*N (effective capacity), identical
     for both soils and at n=40 / n=128.
  R6 shared, inherited from Stage 8: D_scale = stage7_D_scale (12.2582),
     D_film, film K scale, matric wetting (d_cut 50 um), T = 1500, dt.
     Every per-soil Stage 9-16 key is stripped so nothing per-soil leaks in.

Steps (stop early on failure):  python stage19_run.py step0|step1|step2|step3|step4
"""
from __future__ import annotations

import copy
import csv
import hashlib
import json
import os
import sys
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import dual_porosity
import metrics
from model import load_yaml

SOILS = ("sand", "vertisol")
COLORS = {"sand": "#d62728", "vertisol": "#4b0082"}
RESULTS_DIR = os.path.join("results", "stage19")
BIOLOGY_SHA256 = "c76d9869dc5847f830ac3caadcce4c68986a9180875f3333504c2482a4b0bf5f"
B0_DENSITY = 0.0125                      # 20 / 1600, the v2 design density
OM_DENSITIES = tuple(float(x) for x in np.geomspace(0.037, 1.0, 5))   # Step 1 map
THRESHOLD = 0.3
EQUAL_CF = 0.5                           # Step 3 control (a)
T_SHARED = 1500


def check_biology():
    with open("configs/biology.yaml", "rb") as f:
        assert hashlib.sha256(f.read()).hexdigest() == BIOLOGY_SHA256, "R3: biology.yaml modified"


def structure_cfg(soil, n, coating_fraction=None, seed_offset=0, shuffle=False):
    cfg = copy.deepcopy(load_yaml(f"configs/soil_{soil}_stage19.yaml"))
    cfg["grid"]["n"] = int(n)
    cfg["T"] = T_SHARED
    cfg["pore"]["seed"] = int(cfg["pore"]["seed"]) + seed_offset
    if coating_fraction is not None:
        cfg["pore"]["coating_fraction"] = float(coating_fraction)
    return cfg


def base_mapping():
    """main's mapping with every per-soil / later-stage key removed, so the
    shared Stage 7/8 values are what `_variantc_*` resolve to."""
    m = copy.deepcopy(load_yaml("configs/mapping.yaml"))
    bio = load_yaml("configs/biology.yaml")
    for k in list(m):
        if k.startswith(("stage9_", "stage10_", "stage16_", "stage8_om_total", "stage7_om_fraction")):
            del m[k]
        elif k.startswith("stage11_") and k != "stage11_implicit_tol":
            del m[k]
    m["om_mode"] = dual_porosity.OM_MODE_COATING
    m["k_dis_slow"] = float(bio["k_dis"])          # R3: one dissolution rate
    m["om_fine_bias"] = 0.0                         # R3: no d^-bias weighting
    m["stage19_B0_density"] = B0_DENSITY
    return m


def K_density_rule(m):
    """R5's shared habitat density, fixed by rule (not tuned): the mean of the
    two soils' NATIVE effective capacity density sum(K_t)/N at n=128, their
    configured seeds. Used identically at n=40 and n=128."""
    vals = {}
    for s in SOILS:
        g = dual_porosity.build_dual_grids(structure_cfg(s, 128), dict(m, stage19_om_density=0.1))
        _, _, scale = dual_porosity.habitat_properties(g["macro_mask"])
        k_cl = float(m.get("stage4_variantC_cluster_K_scale", 0.75))
        K_t = np.where(g["macro_mask"], g["K"] * (1.0 + k_cl * scale), 0.0)
        vals[s] = float(K_t.sum()) / K_t.size
    return float(np.mean(list(vals.values()))), vals


def mapping_for(om_density, K_density, D_scale=None):
    m = base_mapping()
    m["stage19_om_density"] = float(om_density)
    m["stage19_K_density"] = float(K_density)
    if D_scale is not None:
        m["stage7_D_scale"] = float(D_scale)
    return m


# ------------------------------------------------------------- simulation

def run_job(spec):
    bio = load_yaml("configs/biology.yaml")
    out = dual_porosity.simulate_dual_porosity(bio, spec["cfg"], spec["m"], wicking_enabled=False)
    g = out["grids"]
    n = g["n"]
    N = n * n
    # R5 assertions, every run
    assert abs(g["om_total_placed"] / (spec["m"]["stage19_om_density"] * N) - 1) < 1e-9, "R5 OM"
    assert abs(g["B0"].sum() / (B0_DENSITY * N) - 1) < 1e-9, "R5 B0"
    assert abs(out["K_effective_sum"] / (spec["m"]["stage19_K_density"] * N) - 1) < 1e-9, "R5 K_t"
    R = out["R_t"]
    pk = int(np.argmax(R))
    floor = float(bio["m0"]) * float(g["B0"].sum())
    gsh = out["growth_term_t"][pk] / max(out["growth_term_t"][pk] + out["maintenance_term_t"][pk], 1e-300)
    res = {
        "key": spec["key"], "soil": spec["soil"], "n": n,
        "R_t": R, "delivered_wet_t": out["delivered_wet_t"], "delivered_dry_t": out["delivered_dry_t"],
        "wet_habitat_frac_t": out["wet_habitat_frac_t"],
        "B_peak_over_B0": float(out["B_total_t"].max() / g["B0"].sum()),
        "B_final_over_B0": float(out["B_total_t"][-1] / g["B0"].sum()),
        "R_peak_over_floor": float(R.max() / floor),
        "growth_share_at_peak": float(gsh),
        "fuel_remaining_final": float(out["fuel_remaining_t"][-1]),
        "peak_frac": pk / (len(R) - 1), "R_end_over_peak": float(R[-1] / R.max()),
        "coating_fraction": float(g["coating_fraction"]),
        "dry_habitat_frac": float(g["dry_macro_fraction_of_habitat"]),
        "delivered_wet_total": float(out["delivered_wet_t"].sum()),
        "delivered_dry_total": float(out["delivered_dry_t"].sum()),
        "diffusion_into_wet_habitat": float(out["diff_wet_t"].sum()) if "diff_wet_t" in out else float("nan"),
        "diffusion_into_drained_habitat": float(out["diff_dry_t"].sum()) if "diff_dry_t" in out else float("nan"),
        "leak_into_wet_habitat": float(out["leak_wet_t"].sum()) if "leak_wet_t" in out else float("nan"),
        "leak_into_drained_habitat": float(out["leak_dry_t"].sum()) if "leak_dry_t" in out else float("nan"),
        "shape": metrics.describe_shape(R)["shape"],
    }
    res["alive"] = bool(res["B_peak_over_B0"] > 1.0 and res["R_peak_over_floor"] >= 2.0
                        and res["growth_share_at_peak"] > 0.5)
    return res


def run_many(specs, workers=4):
    with Pool(workers) as p:
        return list(p.map(run_job, specs, chunksize=1))


def pair(results):
    by = {r["soil"]: r for r in results}
    d = float(metrics.pairwise_distance(by["sand"]["R_t"], by["vertisol"]["R_t"]))
    return {"distance": d, "both_alive": bool(by["sand"]["alive"] and by["vertisol"]["alive"]),
            "pass": bool(by["sand"]["alive"] and by["vertisol"]["alive"] and d > THRESHOLD), "by": by}


def flat(r):
    return {k: v for k, v in r.items() if not isinstance(v, np.ndarray)}


def save(name, obj):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(os.path.join(RESULTS_DIR, name), "w") as f:
        json.dump(obj, f, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))


def load(name):
    with open(os.path.join(RESULTS_DIR, name)) as f:
        return json.load(f)


# ------------------------------------------------------------- step 0

def step0():
    check_biology()
    bio = load_yaml("configs/biology.yaml")
    K_d, native = K_density_rule(base_mapping())
    f_star = (bio["m0"] + bio["m_s"]) / (bio["r_max"] + bio["m_s"])
    S_star = bio["Ks"] * f_star / (1 - f_star)
    m = mapping_for(0.1, K_d)
    D_scale = float(m["stage7_D_scale"])
    D_matrix = float(m["stage4_variantC_matrix_diffusion"]) * D_scale
    rows = []
    for s in SOILS:
        for n in (40, 128):
            g = dual_porosity.build_dual_grids(structure_cfg(s, n), m)
            allowed = ~g["solid_mask"]
            dist = dual_porosity.bfs_geodesic_distance(g["macro_mask"], allowed)
            dtr = dist[g["trapped_mask"]]
            dtr = dtr[np.isfinite(dtr)]
            cf = g["coating_fraction"]
            n_coat, n_trap = int(g["coating_mask"].sum()), int(g["trapped_mask"].sum())
            L = float(np.median(dtr)) if dtr.size else float("nan")
            rows.append({
                "soil": s, "n": n, "N": n * n, "coating_fraction": cf,
                "native_K_t_density": native[s] if n == 128 else None,
                "n_regions": int(g["region_labels"].max()),
                "habitat_cells": int(g["macro_mask"].sum()),
                "porous_matrix_cells": int(g["porous_matrix_mask"].sum()),
                "coating_cells": n_coat, "trapped_cells": n_trap,
                "dry_habitat_frac": float(g["dry_macro_fraction_of_habitat"]),
                "trapped_median_geodesic_cells": L,
                "trapped_median_geodesic_um": L * float(g["voxel_um"]),
                "trapped_unreachable": int((~np.isfinite(dist[g["trapped_mask"]])).sum()),
                "t_diffuse_trapped": L * L / D_matrix if np.isfinite(L) else float("nan"),
                # coating OM per coating cell / trapped OM per trapped cell, per unit OM_density
                "coating_OM_per_cell_per_density": cf * n * n / max(n_coat, 1),
                "trapped_OM_per_cell_per_density": (1 - cf) * n * n / max(n_trap, 1),
            })
    supply = [{"OM_density": x, "OM_per_habitat_capacity": x / K_d,
               **{f"supply_over_maint_{s}": bio["k_dis"] * r["coating_fraction"] * x / (bio["m0"] * B0_DENSITY)
                  for s in SOILS for r in rows if r["soil"] == s and r["n"] == 128}}
              for x in OM_DENSITIES]
    timescales = {"t_growth_max_rate": 1.0 / bio["r_max"], "t_dissolve": 1.0 / bio["k_dis"],
                  "t_leak": 1.0 / float(m["k_leak"]), "D_matrix": D_matrix,
                  "D_film_eff": float(m["D_film"]) * D_scale, "f_star": f_star, "S_star": S_star,
                  "window_time": T_SHARED * bio["dt"]}
    out = {"K_density": K_d, "native_K_t_density_n128": native, "B0_density": B0_DENSITY,
           "OM_densities": OM_DENSITIES, "grids": rows, "supply": supply, "timescales": timescales}
    save("stage19_step0_budget.json", out)
    print(json.dumps(out, indent=1, default=str))


# ------------------------------------------------------------- step 1

def step1():
    check_biology()
    K_d = load("stage19_step0_budget.json")["K_density"]
    specs = [{"key": f"om{x:.4f}", "soil": s, "cfg": structure_cfg(s, 40), "m": mapping_for(x, K_d)}
             for x in OM_DENSITIES for s in SOILS]
    res = run_many(specs)
    rows, curves = [], {}
    for x in OM_DENSITIES:
        rs = [r for r in res if r["key"] == f"om{x:.4f}"]
        p = pair(rs)
        curves[x] = {r["soil"]: r["R_t"] for r in rs}
        for r in rs:
            rows.append(dict(flat(r), OM_density=x, distance=p["distance"], both_alive=p["both_alive"]))
    cand = [x for x in OM_DENSITIES
            if all(r["both_alive"] for r in rows if r["OM_density"] == x)]
    sel = max(cand, key=lambda x: [r["distance"] for r in rows if r["OM_density"] == x][0]) if cand else None
    save("stage19_step1_map.json", {"rows": rows, "selected_OM_density": sel, "alive_densities": cand})
    _write_csv("stage19_step1_map.csv", rows)
    _plot_map(curves, "stage19_step1_map.png")
    for r in rows:
        print(f"OM={r['OM_density']:.4f} {r['soil']:8s} alive={r['alive']} B_pk/B0={r['B_peak_over_B0']:.3f} "
              f"R_pk/floor={r['R_peak_over_floor']:.2f} growth_share={r['growth_share_at_peak']:.2f} "
              f"fuel_end={r['fuel_remaining_final']:.3f} peak@{r['peak_frac']:.2f} dist={r['distance']:.3f}")
    print("alive densities:", cand, "selected:", sel)


# ------------------------------------------------------------- steps 2-4

def step2():
    check_biology()
    K_d = load("stage19_step0_budget.json")["K_density"]
    sel = load("stage19_step1_map.json")["selected_OM_density"]
    assert sel is not None, "Step 1 found no density with both soils alive -- stop"
    res = run_many([{"key": "n128", "soil": s, "cfg": structure_cfg(s, 128), "m": mapping_for(sel, K_d)}
                    for s in SOILS])
    p = pair(res)
    save("stage19_step2_n128.json", {"OM_density": sel, "distance": p["distance"],
                                      "both_alive": p["both_alive"], "pass": p["pass"],
                                      "soils": {s: flat(r) for s, r in p["by"].items()},
                                      "curves": {s: r["R_t"] for s, r in p["by"].items()},
                                      "delivery": {s: {"wet": r["delivered_wet_t"], "dry": r["delivered_dry_t"],
                                                       "wet_frac": r["wet_habitat_frac_t"]}
                                                   for s, r in p["by"].items()}})
    for s, r in p["by"].items():
        print(f"n=128 {s}: alive={r['alive']} B_pk/B0={r['B_peak_over_B0']:.3f} R_pk/floor={r['R_peak_over_floor']:.2f} "
              f"growth_share={r['growth_share_at_peak']:.2f} fuel_end={r['fuel_remaining_final']:.3f} peak@{r['peak_frac']:.2f}")
    print(f"distance={p['distance']:.3f} pass={p['pass']}")


def step3():
    check_biology()
    K_d = load("stage19_step0_budget.json")["K_density"]
    s2 = load("stage19_step2_n128.json")
    sel = s2["OM_density"]
    specs = []
    for s in SOILS:
        specs.append({"key": "a_equal_cf", "soil": s, "cfg": structure_cfg(s, 128, coating_fraction=EQUAL_CF),
                      "m": mapping_for(sel, K_d)})
        specs.append({"key": "b_transport_off", "soil": s, "cfg": structure_cfg(s, 128),
                      "m": mapping_for(sel, K_d, D_scale=1e-9)})
    res = run_many(specs)
    out = {"baseline": {"distance": s2["distance"], "both_alive": s2["both_alive"]}}
    for key in ("a_equal_cf", "b_transport_off"):
        p = pair([r for r in res if r["key"] == key])
        out[key] = {"distance": p["distance"], "both_alive": p["both_alive"],
                    "soils": {s: flat(r) for s, r in p["by"].items()},
                    "curves": {s: p["by"][s]["R_t"] for s in SOILS},
                    "max_dR_over_peak": {s: float(np.abs(np.asarray(p["by"][s]["R_t"]) - np.asarray(s2["curves"][s])).max()
                                              / max(np.asarray(s2["curves"][s]).max(), 1e-300)) for s in SOILS}}
    out["c_delivery"] = {s: {"wet_habitat_frac": float(np.mean(s2["delivery"][s]["wet_frac"])),
                             "delivered_wet_total": float(np.sum(s2["delivery"][s]["wet"])),
                             "delivered_dry_total": float(np.sum(s2["delivery"][s]["dry"]))} for s in SOILS}
    save("stage19_step3_controls.json", out)
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    panels = [("n=128 baseline (OM_density=%.3f; cf 0.8 / 0.3)" % sel, s2["curves"], s2["distance"]),
              ("(a) equal coating_fraction = 0.5", out["a_equal_cf"]["curves"], out["a_equal_cf"]["distance"]),
              ("(b) transport off (D_scale = 1e-9)", out["b_transport_off"]["curves"], out["b_transport_off"]["distance"])]
    for ax, (title, cs, d) in zip(axes, panels):
        for s in SOILS:
            R = np.asarray(cs[s])
            ax.plot(np.arange(len(R)) * 0.05, R, color=COLORS[s], label=f"{s} (R_peak={R.max():.2f})")
        ax.axhline(0.02 * B0_DENSITY * 128 * 128, color="gray", ls=":", lw=1, label="maintenance floor m0*B0")
        ax.set_title(f"{title}\ndistance = {d:.3f}", fontsize=9)
        ax.set_xlabel("time")
        ax.set_ylabel("R(t)")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, "stage19_step3_controls.png"), dpi=130)
    plt.close(fig)
    print(json.dumps({k: (v if k != "a_equal_cf" and k != "b_transport_off" else
                          {kk: vv for kk, vv in v.items() if kk not in ("soils", "curves")})
                      for k, v in out.items()}, indent=1, default=lambda o: o.tolist()))


def step4():
    raise SystemExit("Step 4 runs only if Step 3 passes (alive AND distinct at n=128).")


# ------------------------------------------------------------- io helpers

def _write_csv(name, rows):
    keys = list(rows[0])
    with open(os.path.join(RESULTS_DIR, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})


def _plot_map(curves, name):
    fig, axes = plt.subplots(2, len(curves), figsize=(4.2 * len(curves), 7.5), squeeze=False)
    for j, (x, cs) in enumerate(curves.items()):
        for s, R in cs.items():
            t = np.arange(len(R)) * 0.05
            axes[0, j].plot(t, R, color=COLORS[s], label=s)
            axes[1, j].plot(t, R / R.max(), color=COLORS[s], label=s)
        axes[0, j].set_title(f"OM_density={x:.3f}, n=40: R(t)", fontsize=9)
        axes[1, j].set_title("R(t) / own peak", fontsize=9)
        for ax in axes[:, j]:
            ax.set_xlabel("time")
            ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, name), dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    {"step0": step0, "step1": step1, "step2": step2, "step3": step3, "step4": step4}[sys.argv[1]]()
