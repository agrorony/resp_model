"""
Stage 17 runner -- reconciling the v3 branch's equal-totals / robustness
machinery with main's Stage 1-16 model, data and SS4 success battery.

See docs/stage_results/STAGE17_RESULTS.md for the verdict. Everything below
(rules, ladders, arms, selection rule, enzyme constants) is fixed IN THIS FILE
before any result is read; nothing is re-tuned after the fact.

Sources of truth (main): dual_porosity.py physics, configs/, Stage 16's SS4
criteria reused UNCHANGED by import from stage16_run.py (sand_pass_fail,
vertisol_pass_fail, recruitment_timing / recruitment_ok with its frac@2% <
0.90 and t_final > 0.10 bars, hotspot_selectivity, distinctness > 0.3), at
n=128 (primary) and n=40.

Ported from v3 (claude/brave-brahmagupta-edjq14, MODEL_SPEC SS13-14):
  R1 biology frozen (SHA-256 asserted)
  R2 equal totals across soils: ONE om_total, ONE D_scale, ONE T, and total
     habitat sum(K) renormalized to a shared density (`stage17_K_density`)
  R3 no OM hoarding: max cell OM <= 30 x mean porous-matrix OM
  S3 seed robustness (3 seed sets), S4 shuffle control (`pore.shuffle_seed`)
  P3 enzyme-mediated depolymerization (`mapping.enzyme`), S6 +-50% sensitivity

Tracks (user decision, see STAGE17_RESULTS.md SS1):
  L (primary)  main's operative Stage-16 setup -- literature PSDs
               (psd_parametric), saturation matric_regions, Sand + Vertisol
               under SS4; Loess reported descriptively.
  R (as asked) real measured PSDs (psd_measured), saturation global_theta
               (theta=0.6, configs/soil_<s>_measured.yaml), all three soils.

numpy-only physics, as always. Usage:  python stage17_run.py
"""
from __future__ import annotations

import copy
import csv
import hashlib
import json
import os
import time
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from model import load_yaml
import dual_porosity
import metrics
import stage16_run as s16

RESULTS_DIR = os.path.join("results", "stage17")
GRIDS = (128, 40)
SOILS_SS4 = ("sand", "vertisol")
BIOLOGY_SHA256 = "c76d9869dc5847f830ac3caadcce4c68986a9180875f3333504c2482a4b0bf5f"  # R1: main's configs/biology.yaml, unchanged

# ---------------------------------------------------------------- fixed rules
T_SHARED = 3000              # R2: one window for every soil. 3000 is the longest
                             # window main ever used (Sand, Stage 10) -- fixed a
                             # priori, never tuned against SS4.
OM_LADDER = (3000.0, 12247.0, 50000.0)   # Stage 16's own anchors (Sand working
D_LADDER = (0.7, 7.5, 80.0)              # dose, geometric midpoint, Vertisol's)
R3_MAX_OVER_MEAN = 30.0
OM_BIAS_MAIN = 2.5           # main's shared_fine_pore exponent (mapping.yaml)
OM_BIAS_R3 = 0.75            # largest exponent on a 0.25 grid with worst-case
                             # max/mean <= 30 over 3 soils x 2 grids x 3 seeds
                             # (21.1 at 0.75, 42.1 at 1.0) -- derived from R3
                             # alone, before any SS4 run.
# Enzyme constants: v3's a_E, d_E, f_E; K_E ported as v3's DIMENSIONLESS ratio
# kappa = K_E / E*, E* = a_E * (B0_total / n_cells) / d_E (the enzyme level the
# initial biomass density sustains). v3: K_E=0.05, E*=0.125 -> kappa=0.4.
ENZ = {"a_E": 1.0, "d_E": 0.1, "f_E": 0.1, "kappa": 0.4}

ARMS = {
    "F_off": {"bias": OM_BIAS_MAIN, "enzyme": False,
              "label": "main's shared_fine_pore (bias 2.5), enzyme off"},
    "F_on": {"bias": OM_BIAS_MAIN, "enzyme": True,
             "label": "main's shared_fine_pore (bias 2.5), enzyme on"},
    "U_off": {"bias": OM_BIAS_R3, "enzyme": False,
              "label": "R3-compliant shared rule (bias 0.75), enzyme off"},
    "U_on": {"bias": OM_BIAS_R3, "enzyme": True,
             "label": "R3-compliant shared rule (bias 0.75), enzyme on"},
}


def check_biology():
    with open("configs/biology.yaml", "rb") as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    assert digest == BIOLOGY_SHA256, "R1 violated: configs/biology.yaml was modified"
    return digest


# -------------------------------------------------------------- config build

def native_K_density(mapping_cfg):
    """R2's shared habitat density: the mean, over the three soils, of each
    soil's NATIVE sum(K)/n_cells at that grid (literature PSD, Track L). A
    rule, not a fit: it keeps the overall habitat scale where main had it."""
    out = {}
    m = s16.shared_mapping(mapping_cfg, {"stage16_om_total_shared": 10000.0},
                           {"stage16_D_scale_shared": 1.0})
    for n in GRIDS:
        vals = [dual_porosity.build_dual_grids(s16.structure_cfg(s, n, 10), m)["K"].sum() / (n * n)
                for s in ("sand", "vertisol", "loess")]
        out[n] = float(np.mean(vals))
    return out


def mapping_for(mapping_cfg, n, om=None, D=None, bias=OM_BIAS_MAIN, enzyme=False,
                K_density=None, per_soil=None, kappa=None, f_E=None):
    """Shared (R2) mapping unless `per_soil` = {"om": {...}, "D": {...}}
    (only used to reproduce/diagnose Stage 16's non-R2 rung-3 point)."""
    if per_soil is None:
        m = s16.shared_mapping(mapping_cfg, {"stage16_om_total_shared": float(om)},
                               {"stage16_D_scale_shared": float(D)})
    else:
        m = s16.shared_mapping(
            mapping_cfg, {f"stage16_om_total_{s}": v for s, v in per_soil["om"].items()},
            {f"stage16_D_scale_{s}": v for s, v in per_soil["D"].items()})
    m["om_fine_bias"] = float(bias)
    if K_density is not None:
        m["stage17_K_density"] = float(K_density)
    if enzyme:
        k = ENZ["kappa"] if kappa is None else kappa
        E_star = ENZ["a_E"] * (float(mapping_cfg["B0_total"]) / (n * n)) / ENZ["d_E"]
        m["enzyme"] = {"enabled": True, "a_E": ENZ["a_E"], "d_E": ENZ["d_E"],
                       "K_E": k * E_star, "f_E": ENZ["f_E"] if f_E is None else f_E}
    return m


def structure_for(soil, n, T, track="L", seed_offset=0, shuffle=False):
    if track == "L":
        cfg = s16.structure_cfg(soil, n, T)
    else:
        cfg = copy.deepcopy(load_yaml(f"configs/soil_{soil}_measured.yaml"))
        cfg["grid"]["n"] = n
        cfg["T"] = int(T)
    cfg["pore"]["seed"] = int(cfg["pore"]["seed"]) + int(seed_offset)
    if shuffle:
        cfg["pore"]["shuffle_seed"] = 12345 + int(seed_offset)
    return cfg


# -------------------------------------------------------------- job runner

def run_job(spec):
    """One soil, one grid, one configuration. Returns a picklable summary."""
    biology_cfg = load_yaml("configs/biology.yaml")
    t0 = time.time()
    out = dual_porosity.simulate_dual_porosity(
        biology_cfg, spec["structure"], spec["mapping"], wicking_enabled=False,
        track_region_biomass=True, debug_growth_distribution=spec.get("graded", False))
    g = out["grids"]
    soil = spec["soil"]
    res = {k: spec[k] for k in ("key", "soil", "n")}
    res.update({
        "R_t": out["R_t"], "runtime_s": time.time() - t0,
        "om_total_init": g["om_total_init"], "om_total_placed": g["om_total_placed"],
        "om_hoard": g["om_hoard"], "sumK": float(g["K"].sum()), "theta": g["theta"],
        "n_regions": int(g["region_labels"].max()),
        "fuel_end": float(out["fuel_remaining_t"][-1]),
        "R_peak": float(out["R_t"].max()),
        "B_peak_over_B0": float(out["B_total_t"].max() / g["B0"].sum()) if g["B0"].sum() > 0 else float("nan"),
        "timing": s16.recruitment_timing(out),
        "hotspot": s16.hotspot_selectivity(out),
        "shape": metrics.describe_shape(out["R_t"]),
    })
    if soil in s16.PASS_FAIL:
        checks, vals = s16.PASS_FAIL[soil](out)
        res["checks"] = {k: bool(v) for k, v in checks.items()}
        res["vals"] = vals
        res["shape_ok"] = bool(all(checks.values()))
    if spec.get("graded"):
        hm = s16.per_region_half_max(out)
        res["half_max"] = ({"p25": float(np.percentile(hm, 25)), "median": float(np.median(hm)),
                            "p75": float(np.percentile(hm, 75)), "max": float(hm.max())}
                           if hm.size else None)
    return res


def run_all(specs, workers=4):
    with Pool(workers) as p:
        return list(p.map(run_job, specs, chunksize=1))


# -------------------------------------------------------------- scoring

def score(res_by_soil):
    """SS4 for one configuration at one grid. (a) shape per soil, (b)
    Vertisol recruitment timing, (d) distinctness; (c) reported."""
    s, v = res_by_soil["sand"], res_by_soil["vertisol"]
    dist = float(metrics.pairwise_distance(s16.normalized_time_resample(s["R_t"]),
                                           s16.normalized_time_resample(v["R_t"])))
    b_ok = s16.recruitment_ok(v["timing"])
    return {"a_sand": s["shape_ok"], "a_vert": v["shape_ok"], "b_recruit": bool(b_ok),
            "d_dist": dist > s16.DISTINCTNESS_THRESHOLD, "dist": dist,
            "n_pass": int(s["shape_ok"]) + int(v["shape_ok"]) + int(b_ok)
            + int(dist > s16.DISTINCTNESS_THRESHOLD),
            "frac2": v["timing"]["frac_at_2pct"], "t_final": v["timing"]["t_final"],
            "hotspot_vert": v["hotspot"]["hotspot_pct"],
            "r3_ok": bool(max(s["om_hoard"], v["om_hoard"]) <= R3_MAX_OVER_MEAN),
            "om_hoard_max": float(max(s["om_hoard"], v["om_hoard"]))}


def _write_csv(path, rows):
    if not rows:
        return
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})


def _flat(res):
    r = {k: res[k] for k in ("key", "soil", "n", "om_total_init", "om_total_placed", "om_hoard",
                             "sumK", "theta", "n_regions", "fuel_end", "R_peak",
                             "B_peak_over_B0", "runtime_s")}
    r.update({f"timing_{k}": v for k, v in res["timing"].items()})
    r["hotspot_pct"] = res["hotspot"]["hotspot_pct"]
    r["shape_desc"] = res["shape"]["shape"]
    if "checks" in res:
        r["shape_ok"] = res["shape_ok"]
        r["checks_passed"] = f"{sum(res['checks'].values())}/{len(res['checks'])}"
        r["checks_failed"] = "; ".join(k for k, v in res["checks"].items() if not v)
        r.update({f"val_{k}": v for k, v in res["vals"].items()})
    return r


# -------------------------------------------------------------- S5 check

def mass_conservation_check():
    """S5: the implicit Jacobi transport (+tol early stop) conserves total
    substrate only up to its convergence error -- measured, not assumed."""
    rng = np.random.default_rng(0)
    rows = []
    for tol in (None, 1e-7):
        D = rng.uniform(0, 3, (64, 64))
        D[rng.random(D.shape) < 0.3] = 0.0
        u = rng.uniform(0, 5, D.shape)
        un = dual_porosity.implicit_diffusion_step(u, D, 1.0, 0.05, 90, 0.85, tol=tol)
        rows.append({"tol": str(tol), "mass_before": float(u.sum()), "mass_after": float(un.sum()),
                     "relative_error": float(abs(un.sum() - u.sum()) / u.sum())})
    return rows


# -------------------------------------------------------------- main

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    t_start = time.time()
    bio_sha = check_biology()
    biology_cfg = load_yaml("configs/biology.yaml")
    mapping_cfg = load_yaml("configs/mapping.yaml")
    summary = {"biology_sha256": bio_sha}

    Kd = native_K_density(mapping_cfg)
    summary["K_density"] = Kd
    print("R2 shared habitat density (mean native sum(K)/N):", Kd)
    summary["S5_mass_conservation"] = mass_conservation_check()
    print("S5:", summary["S5_mass_conservation"])

    # ---- Part A: audit Stage 16's best attempt (probe 9, rung 3) against R1-R3
    rung3 = {"om": {"sand": 3000.0, "vertisol": 50000.0}, "D": {"sand": 80.0, "vertisol": 0.7}}
    specs = []
    for n in GRIDS:
        m = mapping_for(mapping_cfg, n, per_soil=rung3)
        for soil in SOILS_SS4:
            specs.append({"key": f"A_rung3_n{n}", "soil": soil, "n": n, "mapping": m,
                          "structure": structure_for(soil, n, mapping_cfg[f"stage10_T_{soil}"])})
    # ---- Part D: enzyme at Stage 16's own working point (NOT R2) -- the most
    # direct test of "does enzyme gating fix (b) where Stage 16 failed it"
    for n in GRIDS:
        for tag, kap in (("v3raw", None), ("kappa", ENZ["kappa"])):
            m = mapping_for(mapping_cfg, n, per_soil=rung3, enzyme=True)
            if tag == "v3raw":
                m["enzyme"]["K_E"] = 0.05   # v3's literal constant
            for soil in SOILS_SS4:
                specs.append({"key": f"D_rung3_enz_{tag}_n{n}", "soil": soil, "n": n, "mapping": m,
                              "structure": structure_for(soil, n, mapping_cfg[f"stage10_T_{soil}"])})
    # ---- Part B: Track R (real measured PSDs, global_theta), R2 midpoint
    for n in GRIDS:
        m = mapping_for(mapping_cfg, n, OM_LADDER[1], D_LADDER[1], K_density=Kd[n])
        for soil in ("sand", "vertisol", "loess"):
            specs.append({"key": f"B_trackR_n{n}", "soil": soil, "n": n, "mapping": m,
                          "structure": structure_for(soil, n, T_SHARED, track="R")})
    # ---- Part C: the four R2 arms x 9-probe ladder (<=12 per arm, R7)
    for arm, a in ARMS.items():
        for om in OM_LADDER:
            for D in D_LADDER:
                for n in GRIDS:
                    m = mapping_for(mapping_cfg, n, om, D, bias=a["bias"], enzyme=a["enzyme"],
                                    K_density=Kd[n])
                    for soil in SOILS_SS4:
                        specs.append({"key": f"C_{arm}_om{int(om)}_D{D}_n{n}", "soil": soil, "n": n,
                                      "mapping": m, "structure": structure_for(soil, n, T_SHARED)})
    print(f"running {len(specs)} jobs (parts A-D)...")
    results = run_all(specs)
    _write_csv(os.path.join(RESULTS_DIR, "stage17_all_runs.csv"), [_flat(r) for r in results])
    by_key = {}
    for r in results:
        by_key.setdefault(r["key"], {})[r["soil"]] = r

    # score A, D, C
    scored = {}
    for key, rs in by_key.items():
        if "sand" in rs and "vertisol" in rs:
            scored[key] = score(rs)
    _write_csv(os.path.join(RESULTS_DIR, "stage17_ss4_scores.csv"),
               [dict(key=k, **v) for k, v in scored.items()])
    summary["scores"] = scored
    summary["trackR"] = {k: {s: _flat(r) for s, r in v.items()} for k, v in by_key.items()
                         if k.startswith("B_")}

    # ---- selection per arm (fixed in advance): most SS4 criteria at n=128,
    # then at n=40, then larger n=128 distinctness
    selected = {}
    for arm in ARMS:
        cands = []
        for om in OM_LADDER:
            for D in D_LADDER:
                base = f"C_{arm}_om{int(om)}_D{D}"
                s128, s40 = scored[f"{base}_n128"], scored[f"{base}_n40"]
                cands.append(((s128["n_pass"], s40["n_pass"], s128["dist"]), om, D, base))
        cands.sort(reverse=True)
        selected[arm] = {"om": cands[0][1], "D": cands[0][2], "base": cands[0][3]}
    summary["selected"] = selected
    print("selected:", selected)

    # ---- Part E: S3 seeds + S4 shuffle (+ Loess descriptive, + graded metric)
    # ---- Part F: S6 enzyme sensitivity on the enzyme arms' selected probes
    specs = []
    for arm, sel in selected.items():
        a = ARMS[arm]
        for n in GRIDS:
            m = mapping_for(mapping_cfg, n, sel["om"], sel["D"], bias=a["bias"],
                            enzyme=a["enzyme"], K_density=Kd[n])
            for soil in SOILS_SS4:
                for so in (10, 20):
                    specs.append({"key": f"E_{arm}_seed{so}_n{n}", "soil": soil, "n": n, "mapping": m,
                                  "structure": structure_for(soil, n, T_SHARED, seed_offset=so)})
                specs.append({"key": f"E_{arm}_shuffle_n{n}", "soil": soil, "n": n, "mapping": m,
                              "structure": structure_for(soil, n, T_SHARED, shuffle=True)})
            specs.append({"key": f"E_{arm}_loess_n{n}", "soil": "loess", "n": n, "mapping": m,
                          "structure": structure_for("loess", n, T_SHARED)})
            specs.append({"key": f"E_{arm}_graded_n{n}", "soil": "vertisol", "n": n, "mapping": m,
                          "structure": structure_for("vertisol", n, T_SHARED), "graded": True})
            if a["enzyme"]:
                for tag, kw in (("kappa0.5", {"kappa": 0.2}), ("kappa1.5", {"kappa": 0.6}),
                                ("fE0.5", {"f_E": 0.05}), ("fE1.5", {"f_E": 0.15})):
                    mm = mapping_for(mapping_cfg, n, sel["om"], sel["D"], bias=a["bias"],
                                     enzyme=True, K_density=Kd[n], **kw)
                    for soil in SOILS_SS4:
                        specs.append({"key": f"F_{arm}_{tag}_n{n}", "soil": soil, "n": n,
                                      "mapping": mm, "structure": structure_for(soil, n, T_SHARED)})
    print(f"running {len(specs)} jobs (parts E-F)...")
    results2 = run_all(specs)
    _write_csv(os.path.join(RESULTS_DIR, "stage17_robustness_runs.csv"), [_flat(r) for r in results2])
    by2 = {}
    for r in results2:
        by2.setdefault(r["key"], {})[r["soil"]] = r

    robust = {}
    for arm, sel in selected.items():
        for n in GRIDS:
            primary = by_key[f"{sel['base']}_n{n}"]
            sets = [primary] + [by2[f"E_{arm}_seed{so}_n{n}"] for so in (10, 20)]
            sc = [score(x) for x in sets]
            within = []
            for soil in SOILS_SS4:
                curves = [s16.normalized_time_resample(x[soil]["R_t"]) for x in sets]
                for i in range(3):
                    for j in range(i + 1, 3):
                        within.append(float(metrics.pairwise_distance(curves[i], curves[j])))
            shuf = score(by2[f"E_{arm}_shuffle_n{n}"])
            loess = by2[f"E_{arm}_loess_n{n}"]["loess"]
            d3 = {}
            for a_, b_ in (("sand", "vertisol"), ("sand", "loess"), ("vertisol", "loess")):
                ra = primary[a_]["R_t"] if a_ != "loess" else loess["R_t"]
                rb = primary[b_]["R_t"] if b_ != "loess" else loess["R_t"]
                d3[f"{a_}-{b_}"] = float(metrics.pairwise_distance(
                    s16.normalized_time_resample(ra), s16.normalized_time_resample(rb)))
            rob = {
                "per_seed_set": [{k: v for k, v in x.items()} for x in sc],
                "S3a_all_pass_same_criteria": [x["n_pass"] for x in sc],
                "S3b_within_mean": float(np.mean(within)),
                "S3b_min_between": float(min(x["dist"] for x in sc)),
                "S4_shuffled": shuf, "S4_dist_drop": float(sc[0]["dist"] - shuf["dist"]),
                "loess_shape": loess["shape"], "loess_R_peak": loess["R_peak"],
                "three_way_distances": d3,
                "graded_half_max_vertisol": by2[f"E_{arm}_graded_n{n}"]["vertisol"].get("half_max"),
            }
            if ARMS[arm]["enzyme"]:
                rob["S6"] = {tag: score(by2[f"F_{arm}_{tag}_n{n}"])
                             for tag in ("kappa0.5", "kappa1.5", "fE0.5", "fE1.5")}
            robust[f"{arm}_n{n}"] = rob
    summary["robustness"] = robust
    summary["runtime_min"] = (time.time() - t_start) / 60.0

    with open(os.path.join(RESULTS_DIR, "stage17_summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    plot_curves(by_key, selected)
    print(f"done in {summary['runtime_min']:.1f} min")


def plot_curves(by_key, selected):
    colors = {"sand": "#d62728", "vertisol": "#4b0082", "loess": "#2ca02c"}
    fig, axes = plt.subplots(2, 5, figsize=(24, 8.5))
    cols = [("A_rung3", "Stage 16 probe 9 (rung 3, NOT R2)")] + \
           [(selected[a]["base"], f"{a}: {ARMS[a]['label']}\nom={selected[a]['om']:.0f}, D={selected[a]['D']}")
            for a in ARMS]
    for j, (base, title) in enumerate(cols):
        for i, n in enumerate(GRIDS):
            ax = axes[i, j]
            for soil, r in by_key[f"{base}_n{n}"].items():
                R = r["R_t"]
                ax.plot(np.linspace(0, 1, len(R)), R / max(R.max(), 1e-12), color=colors[soil],
                        label=f"{soil} (R_peak={R.max():.2f})")
            ax.set_title(f"{title}\nn={n}", fontsize=8)
            ax.set_xlabel("fraction of window")
            ax.set_ylabel("R(t) / own peak")
            ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, "stage17_selected_curves.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
