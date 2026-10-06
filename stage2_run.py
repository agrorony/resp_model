"""
Stage-2 runner: sweep theta (global_theta saturation mode) for each soil,
record respiration + connectivity diagnostics, and write the sweep figures
required by soil_respiration_prompt_stage2.md SS4-7.

Usage:
    python stage2_run.py                  # full theta sweep, all figures
    python stage2_run.py --theta 0.4       # also dump a single-theta water map
"""
from __future__ import annotations

import argparse
import copy
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from model import load_yaml, simulate
import metrics

SOILS = ["A", "B", "C"]
COLORS = {"A": "#1f77b4", "B": "#d62728", "C": "#2ca02c"}
THETA_SWEEP = [1.0, 0.8, 0.6, 0.4, 0.2]
MAP_THETA = 0.4  # theta chosen for the water-filled maps / single-theta curves


def run_one(soil, theta, biology_cfg, mapping_cfg):
    structure_cfg = copy.deepcopy(load_yaml(f"configs/soil_{soil}.yaml"))
    structure_cfg["saturation"] = {"mode": "global_theta", "theta": theta}
    out = simulate(biology_cfg, structure_cfg, mapping_cfg)
    conn = metrics.connectivity_diagnostics(out["grids"])
    return out, conn


def sweep(results_dir="results/stage2"):
    biology_cfg = load_yaml("configs/biology.yaml")
    mapping_cfg = load_yaml("configs/mapping.yaml")
    os.makedirs(results_dir, exist_ok=True)

    # sweep_data[soil][theta] = (out, conn)
    sweep_data = {s: {} for s in SOILS}
    for theta in THETA_SWEEP:
        for soil in SOILS:
            out, conn = run_one(soil, theta, biology_cfg, mapping_cfg)
            sweep_data[soil][theta] = (out, conn)

    write_sweep_csv(sweep_data, results_dir)
    plot_sweep_respiration(sweep_data, results_dir)
    plot_sweep_connectivity(sweep_data, results_dir)
    distinctness_by_theta = write_distinctness_csv(sweep_data, results_dir)

    # single-theta deliverables at MAP_THETA
    for soil in SOILS:
        out, conn = sweep_data[soil][MAP_THETA]
        plot_water_map(soil, out, results_dir)
    plot_curves_at_theta(sweep_data, MAP_THETA, results_dir)

    return sweep_data, distinctness_by_theta


def write_sweep_csv(sweep_data, results_dir):
    path = os.path.join(results_dir, "stage2_theta_sweep.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["soil", "theta", "R_peak", "R_end", "cum_co2_final",
                    "B_final_total", "nontrivial", "largest_cluster_frac",
                    "om_biomass_connectivity", "accessible_om_fraction"])
        for soil in SOILS:
            for theta in THETA_SWEEP:
                out, conn = sweep_data[soil][theta]
                w.writerow([
                    soil, theta,
                    float(out["R_t"].max()), float(out["R_t"][-1]),
                    float(out["cum_co2"][-1]), float(out["B_final"].sum()),
                    metrics.is_nontrivial(out),
                    conn["largest_cluster_frac"],
                    conn["om_biomass_connectivity"],
                    conn["accessible_om_fraction"],
                ])
    return path


def write_distinctness_csv(sweep_data, results_dir):
    path = os.path.join(results_dir, "stage2_distinctness_by_theta.csv")
    distinctness_by_theta = {}
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["theta", "AB", "AC", "BC", "min_distance", "success",
                    "nontrivial_A", "nontrivial_B", "nontrivial_C"])
        for theta in THETA_SWEEP:
            out_A, _ = sweep_data["A"][theta]
            out_B, _ = sweep_data["B"][theta]
            out_C, _ = sweep_data["C"][theta]
            success, detail = metrics.emergent_distinctness(out_A, out_B, out_C)
            distinctness_by_theta[theta] = (success, detail)
            d = detail["distances"]
            w.writerow([theta, d["AB"], d["AC"], d["BC"], detail["min_distance"],
                        success, detail["nontrivial"]["A"],
                        detail["nontrivial"]["B"], detail["nontrivial"]["C"]])
    return distinctness_by_theta


def plot_water_map(soil, out, results_dir):
    grids = out["grids"]
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    fields = [
        ("pore radius r", grids["r"], "viridis"),
        ("K(r)", grids["K"], "viridis"),
        ("D used", grids["D"], "viridis"),
        ("water-filled mask", grids["water_mask"].astype(float), "Blues"),
    ]
    for ax, (name, field, cmap) in zip(axes, fields):
        im = ax.imshow(field, cmap=cmap, origin="lower")
        ax.set_title(f"Soil {soil}: {name}")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle(f"Soil {soil} at theta={grids['theta']:.2f}")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, f"stage2_water_map_{soil}_theta{grids['theta']:.2f}.png"), dpi=150)
    plt.close(fig)


def plot_curves_at_theta(sweep_data, theta, results_dir):
    fig, ax = plt.subplots(figsize=(7, 5))
    for soil in SOILS:
        out, _ = sweep_data[soil][theta]
        dt = out["dt"]
        t = np.arange(len(out["R_t"])) * dt
        ax.plot(t, out["R_t"], label=f"Soil {soil}", color=COLORS[soil])
    ax.set_xlabel("time")
    ax.set_ylabel("R(t) - respiration rate")
    ax.set_title(f"Respiration curves at theta={theta:.2f}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, f"stage2_respiration_curves_theta{theta:.2f}.png"), dpi=150)
    plt.close(fig)


def plot_sweep_respiration(sweep_data, results_dir):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for soil in SOILS:
        peaks = [sweep_data[soil][th][0]["R_t"].max() for th in THETA_SWEEP]
        cums = [sweep_data[soil][th][0]["cum_co2"][-1] for th in THETA_SWEEP]
        axes[0].plot(THETA_SWEEP, peaks, "o-", label=f"Soil {soil}", color=COLORS[soil])
        axes[1].plot(THETA_SWEEP, cums, "o-", label=f"Soil {soil}", color=COLORS[soil])
    axes[0].set_xlabel("theta")
    axes[0].set_ylabel("R(t) peak")
    axes[0].set_title("Peak respiration vs theta")
    axes[0].invert_xaxis()
    axes[0].legend()
    axes[1].set_xlabel("theta")
    axes[1].set_ylabel("cumulative CO2 (final)")
    axes[1].set_title("Cumulative CO2 vs theta")
    axes[1].invert_xaxis()
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage2_respiration_vs_theta.png"), dpi=150)
    plt.close(fig)


def plot_sweep_connectivity(sweep_data, results_dir):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    labels = ["largest_cluster_frac", "om_biomass_connectivity", "accessible_om_fraction"]
    titles = ["Largest water cluster fraction", "OM-biomass connectivity", "Accessible-OM fraction"]
    for ax, key, title in zip(axes, labels, titles):
        for soil in SOILS:
            vals = [sweep_data[soil][th][1][key] for th in THETA_SWEEP]
            ax.plot(THETA_SWEEP, vals, "o-", label=f"Soil {soil}", color=COLORS[soil])
        ax.set_xlabel("theta")
        ax.set_ylabel(key)
        ax.set_title(title)
        ax.invert_xaxis()
        ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage2_connectivity_vs_theta.png"), dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", default="results/stage2")
    args = parser.parse_args()

    sweep_data, distinctness_by_theta = sweep(args.results_dir)

    print("theta   AB      AC      BC      min     success")
    for theta in THETA_SWEEP:
        success, detail = distinctness_by_theta[theta]
        d = detail["distances"]
        print(f"{theta:5.2f}  {d['AB']:.3f}   {d['AC']:.3f}   {d['BC']:.3f}   "
              f"{detail['min_distance']:.3f}   {success}")

    for soil in SOILS:
        for theta in THETA_SWEEP:
            out, conn = sweep_data[soil][theta]
            print(f"soil {soil} theta={theta:.2f}: R_peak={out['R_t'].max():.3f} "
                  f"cum={out['cum_co2'][-1]:.3f} B_final={out['B_final'].sum():.3f} "
                  f"largest_cluster={conn['largest_cluster_frac']:.3f} "
                  f"om_biomass_conn={conn['om_biomass_connectivity']:.3f} "
                  f"accessible_om={conn['accessible_om_fraction']:.3f}")


if __name__ == "__main__":
    main()
