"""
Stage 18 runner -- uniform-OM single-compartment architecture on real
measured PSDs (docs/stage_results/STAGE18_RESULTS.md, uniform_om_model.py).

Test procedure, in order, stopping early on failure (no parameter of D(d),
K(d) or the totals is tuned at any point):
  Step 1  distinctness, ONE seed (each soil's configured seed), Sand,
          Vertisol, Loess: pairwise distance (metrics.pairwise_distance,
          threshold 0.3) and non-triviality (main's is_nontrivial).
          Verdict grid = configured n=128; n=40 reported alongside.
  Step 2  only if step 1 passes: 3 seed sets, shuffle control (same PSD,
          arrangement destroyed), mass-conservation check.

Usage:  python stage18_run.py
"""
from __future__ import annotations

import copy
import csv
import itertools
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import metrics
import uniform_om_model as u
from model import load_yaml

SOILS = ("sand", "vertisol", "loess")
COLORS = {"sand": "#d62728", "vertisol": "#4b0082", "loess": "#2ca02c"}
THRESHOLD = 0.3
VERDICT_GRID = 128
GRIDS = (128, 40)
RESULTS_DIR = os.path.join("results", "stage18")


def soil_cfg(soil, n, seed_offset=0):
    cfg = copy.deepcopy(load_yaml(f"configs/soil_{soil}_measured.yaml"))
    cfg["grid"]["n"] = n
    cfg["pore"]["seed"] = int(cfg["pore"]["seed"]) + seed_offset
    return cfg


def check_K_mean_rule(m):
    native = [float(u.K_shape_of_d(u.diameter_field(soil_cfg(s, 128)), m).mean()) for s in SOILS]
    assert abs(np.mean(native) - m["K_mean"]) < 1e-12, "K_mean is not the stated rule's value"
    return dict(zip(SOILS, native))


def run_set(bio, m, n, seed_offset=0, shuffle=False):
    outs = {}
    for s in SOILS:
        cfg = soil_cfg(s, n, seed_offset)
        outs[s] = u.simulate(bio, cfg, m, shuffle_seed=(9000 + seed_offset) if shuffle else None)
    return outs


def summarize(outs):
    dist = {f"{a}-{b}": float(metrics.pairwise_distance(outs[a]["R_t"], outs[b]["R_t"]))
            for a, b in itertools.combinations(SOILS, 2)}
    per = {}
    for s, o in outs.items():
        g = o["grids"]
        per[s] = {
            "nontrivial": u.is_nontrivial(o),
            "R_peak": float(o["R_t"].max()), "peak_frac": float(np.argmax(o["R_t"]) / (len(o["R_t"]) - 1)),
            "R_end_over_peak": float(o["R_t"][-1] / o["R_t"].max()),
            "B_final_over_B0": float(o["B_final"].sum() / g["B0"].sum()),
            "B_peak_over_B0": float(o["B_t"].max() / g["B0"].sum()),
            "cum_co2": float(o["cum_co2"][-1]), "theta": float(g["theta"]),
            "sum_OM": float(g["OM"].sum()), "sum_B0": float(g["B0"].sum()), "sum_K": float(g["K"].sum()),
            "D_mean_wet": float(g["D"][g["water_mask"]].mean()),
            "d_median_um": float(np.median(g["d"])),
            "shape": metrics.describe_shape(o["R_t"])["shape"],
            "transport_mass_err_max": o["transport_mass_err_max"],
        }
    passed = all(v > THRESHOLD for v in dist.values()) and all(p["nontrivial"] for p in per.values())
    return {"distances": dist, "min_distance": min(dist.values()), "per_soil": per, "pass": passed}


def plot(outs_by_grid, path, title):
    fig, axes = plt.subplots(2, len(outs_by_grid), figsize=(6.5 * len(outs_by_grid), 8.5), squeeze=False)
    for j, (n, outs) in enumerate(outs_by_grid.items()):
        for s, o in outs.items():
            t = np.arange(len(o["R_t"])) * o["dt"]
            axes[0, j].plot(t, o["R_t"], color=COLORS[s], label=s)
            axes[1, j].plot(t, o["R_t"] / o["R_t"].max(), color=COLORS[s], label=s)
        axes[0, j].set_title(f"{title}, n={n}: R(t)")
        axes[1, j].set_title(f"n={n}: R(t) / own peak (what the distance metric sees)")
        for ax in axes[:, j]:
            ax.set_xlabel("time")
            ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    bio = load_yaml("configs/biology.yaml")
    m = load_yaml("configs/mapping_stage18.yaml")
    summary = {"K_mean_rule_native_means": check_K_mean_rule(m)}

    # ---------------- step 1
    step1, outs_by_grid = {}, {}
    for n in GRIDS:
        outs = run_set(bio, m, n)
        outs_by_grid[n] = outs
        step1[n] = summarize(outs)
        print(f"step 1, n={n}: distances {step1[n]['distances']}  pass={step1[n]['pass']}")
        for s, p in step1[n]["per_soil"].items():
            print(f"   {s}: nontrivial={p['nontrivial']} R_peak={p['R_peak']:.3f} peak@{p['peak_frac']:.2f} "
                  f"R_end/peak={p['R_end_over_peak']:.3f} B_final/B0={p['B_final_over_B0']:.3f} shape={p['shape']}")
    plot(outs_by_grid, os.path.join(RESULTS_DIR, "stage18_step1_curves.png"), "Step 1 (one seed)")
    with open(os.path.join(RESULTS_DIR, "stage18_step1_timeseries.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["grid_n", "soil", "step", "time", "R_t", "B_total"])
        for n, outs in outs_by_grid.items():
            for s, o in outs.items():
                for i in range(len(o["R_t"])):
                    w.writerow([n, s, i, f"{i * o['dt']:.4f}", f"{o['R_t'][i]:.6g}", f"{o['B_t'][i]:.6g}"])
    summary["step1"] = step1
    step1_pass = step1[VERDICT_GRID]["pass"]
    summary["step1_pass"] = step1_pass

    # ---------------- step 2 (only if step 1 passes)
    if step1_pass:
        step2 = {}
        for n in GRIDS:
            sets = [step1[n]] + [summarize(run_set(bio, m, n, seed_offset=so)) for so in (10, 20)]
            shuf = summarize(run_set(bio, m, n, shuffle=True))
            step2[n] = {"seed_sets": sets, "all_seed_sets_pass": all(x["pass"] for x in sets),
                        "shuffle": shuf,
                        "shuffle_drop": {k: step1[n]["distances"][k] - shuf["distances"][k]
                                         for k in shuf["distances"]}}
        summary["step2"] = step2
    else:
        summary["step2"] = "NOT RUN -- step 1 failed at the verdict grid (procedure: stop early)"

    with open(os.path.join(RESULTS_DIR, "stage18_summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=float)
    print("step1_pass:", step1_pass)
    return 0


if __name__ == "__main__":
    sys.exit(main())
