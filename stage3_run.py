"""
Stage-3 runner (soil_respiration_prompt_stage3.md): theta sweep at BOTH
grid sizes (n=128, n=40) on the three soils' real empirical PSDs, the
outcome-classification + isolation-ratio diagnostics, and the isolation-
ratio-vs-PSD hypothesis test (SS9).

Usage:
    python stage3_run.py
"""
from __future__ import annotations

import copy
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from model import load_yaml, simulate
import metrics
import psd_data

SOILS = ["A", "B", "C"]
SOIL_PSD_DIR = {"A": "data/psd/rehovot", "B": "data/psd/nlm", "C": "data/psd/mishmar"}
SOIL_PSD_NAME = {"A": "rehovot", "B": "nlm", "C": "mishmar"}
COLORS = {"A": "#1f77b4", "B": "#d62728", "C": "#2ca02c"}
GRID_STYLE = {128: "-", 40: "--"}
GRID_MARKER = {128: "o", 40: "s"}
THETA_SWEEP = [1.0, 0.8, 0.6, 0.4, 0.2]
GRIDS = [128, 40]
MAP_THETA = 0.4


def run_one(soil, theta, grid_n, biology_cfg, mapping_cfg):
    structure_cfg = copy.deepcopy(load_yaml(f"configs/soil_{soil}.yaml"))
    structure_cfg["saturation"] = {"mode": "global_theta", "theta": theta}
    structure_cfg["grid"]["n"] = grid_n
    out = simulate(biology_cfg, structure_cfg, mapping_cfg)
    conn = metrics.connectivity_diagnostics(out["grids"])
    return out, conn


def sweep(results_dir="results/stage3"):
    biology_cfg = load_yaml("configs/biology.yaml")
    mapping_cfg = load_yaml("configs/mapping.yaml")
    os.makedirs(results_dir, exist_ok=True)

    psds = {s: psd_data.load_soil_psd(SOIL_PSD_DIR[s]) for s in SOILS}

    # sweep_data[grid_n][soil][theta] = (out, conn)
    sweep_data = {g: {s: {} for s in SOILS} for g in GRIDS}
    for grid_n in GRIDS:
        for theta in THETA_SWEEP:
            for soil in SOILS:
                out, conn = run_one(soil, theta, grid_n, biology_cfg, mapping_cfg)
                sweep_data[grid_n][soil][theta] = (out, conn)

    plot_loaded_psds(psds, results_dir)
    write_sweep_csv(sweep_data, results_dir)
    distinctness = write_distinctness_csv(sweep_data, results_dir)
    plot_respiration_vs_theta(sweep_data, results_dir)
    plot_connectivity_vs_theta(sweep_data, results_dir)
    plot_isolation_vs_theta(sweep_data, results_dir)

    for grid_n in GRIDS:
        for soil in SOILS:
            out, conn = sweep_data[grid_n][soil][MAP_THETA]
            plot_water_map(soil, grid_n, out, results_dir)
        plot_curves_at_theta(sweep_data, grid_n, MAP_THETA, results_dir)

    corr = correlation_analysis(sweep_data, psds, results_dir)

    return sweep_data, distinctness, corr


# ---------------------------------------------------------------- PSD load

def plot_loaded_psds(psds, results_dir):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for soil in SOILS:
        psd = psds[soil]
        label = f"Soil {soil} ({SOIL_PSD_NAME[soil]})"
        axes[0].plot(psd["bin_centers_um"], psd["differential_psd"], "-",
                     color=COLORS[soil], label=label)
        axes[1].plot(psd["bin_edges_um"], psd["cdf_edges"], "-",
                     color=COLORS[soil], label=label)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("pore diameter d (um, log scale)")
    axes[0].set_ylabel("differential PSD (volume fraction / um)")
    axes[0].set_title("Loaded empirical PSDs (density)")
    axes[0].axvspan(30, 150, color="gray", alpha=0.15, label="K(d) habitat window")
    axes[0].legend(fontsize=8)

    axes[1].set_xscale("log")
    axes[1].set_xlabel("pore diameter d (um, log scale)")
    axes[1].set_ylabel("cumulative volume fraction")
    axes[1].set_title("Empirical CDF (used for quantile mapping)")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage3_psd_loaded.png"), dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------- CSVs

def write_sweep_csv(sweep_data, results_dir):
    path = os.path.join(results_dir, "stage3_theta_sweep.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["grid_n", "soil", "theta", "R_peak", "R_end", "cum_co2_final",
                    "B_final_total", "outcome", "largest_cluster_frac",
                    "om_biomass_connectivity", "accessible_om_fraction",
                    "isolation_ratio", "n_clusters", "n_saturated_cells"])
        for grid_n in GRIDS:
            for soil in SOILS:
                for theta in THETA_SWEEP:
                    out, conn = sweep_data[grid_n][soil][theta]
                    w.writerow([
                        grid_n, soil, theta,
                        float(out["R_t"].max()), float(out["R_t"][-1]),
                        float(out["cum_co2"][-1]), float(out["B_final"].sum()),
                        metrics.classify_outcome(out),
                        conn["largest_cluster_frac"],
                        conn["om_biomass_connectivity"],
                        conn["accessible_om_fraction"],
                        conn["isolation_ratio"], conn["n_clusters"],
                        conn["n_saturated_cells"],
                    ])
    return path


def write_distinctness_csv(sweep_data, results_dir):
    path = os.path.join(results_dir, "stage3_distinctness_by_theta.csv")
    distinctness = {}
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["grid_n", "theta", "AB", "AC", "BC", "min_distance", "success",
                    "outcome_A", "outcome_B", "outcome_C"])
        for grid_n in GRIDS:
            for theta in THETA_SWEEP:
                out_A, _ = sweep_data[grid_n]["A"][theta]
                out_B, _ = sweep_data[grid_n]["B"][theta]
                out_C, _ = sweep_data[grid_n]["C"][theta]
                success, detail = metrics.emergent_distinctness(
                    out_A, out_B, out_C, require_nontrivial=False)
                distinctness[(grid_n, theta)] = (success, detail)
                d = detail["distances"]
                oc = detail["outcomes"]
                w.writerow([grid_n, theta, d["AB"], d["AC"], d["BC"],
                            detail["min_distance"], success,
                            oc["A"], oc["B"], oc["C"]])
    return distinctness


# ----------------------------------------------------------------- plots

def plot_water_map(soil, grid_n, out, results_dir):
    grids = out["grids"]
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    fields = [
        ("pore diameter d (um)", grids["r"], "viridis"),
        ("K(d)", grids["K"], "viridis"),
        ("D used", grids["D"], "viridis"),
        ("water-filled mask", grids["water_mask"].astype(float), "Blues"),
    ]
    for ax, (name, field, cmap) in zip(axes, fields):
        im = ax.imshow(field, cmap=cmap, origin="lower")
        ax.set_title(f"Soil {soil}: {name}")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle(f"Soil {soil} (n={grid_n}) at theta={grids['theta']:.2f}")
    fig.tight_layout()
    fname = f"stage3_water_map_{soil}_n{grid_n}_theta{grids['theta']:.2f}.png"
    fig.savefig(os.path.join(results_dir, fname), dpi=150)
    plt.close(fig)


def plot_curves_at_theta(sweep_data, grid_n, theta, results_dir):
    fig, ax = plt.subplots(figsize=(7, 5))
    for soil in SOILS:
        out, _ = sweep_data[grid_n][soil][theta]
        dt = out["dt"]
        t = np.arange(len(out["R_t"])) * dt
        ax.plot(t, out["R_t"], label=f"Soil {soil}", color=COLORS[soil])
    ax.set_xlabel("time")
    ax.set_ylabel("R(t) - respiration rate")
    ax.set_title(f"Respiration curves at theta={theta:.2f}, n={grid_n}")
    ax.legend()
    fig.tight_layout()
    fname = f"stage3_respiration_curves_theta{theta:.2f}_n{grid_n}.png"
    fig.savefig(os.path.join(results_dir, fname), dpi=150)
    plt.close(fig)


def plot_respiration_vs_theta(sweep_data, results_dir):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for grid_n in GRIDS:
        for soil in SOILS:
            peaks = [sweep_data[grid_n][soil][th][0]["R_t"].max() for th in THETA_SWEEP]
            cums = [sweep_data[grid_n][soil][th][0]["cum_co2"][-1] for th in THETA_SWEEP]
            label = f"Soil {soil} (n={grid_n})"
            axes[0].plot(THETA_SWEEP, peaks, GRID_MARKER[grid_n] + GRID_STYLE[grid_n],
                         label=label, color=COLORS[soil])
            axes[1].plot(THETA_SWEEP, cums, GRID_MARKER[grid_n] + GRID_STYLE[grid_n],
                         label=label, color=COLORS[soil])
    axes[0].set_xlabel("theta")
    axes[0].set_ylabel("R(t) peak")
    axes[0].set_title("Peak respiration vs theta (solid=n128, dashed=n40)")
    axes[0].invert_xaxis()
    axes[0].legend(fontsize=7)
    axes[1].set_xlabel("theta")
    axes[1].set_ylabel("cumulative CO2 (final)")
    axes[1].set_title("Cumulative CO2 vs theta (solid=n128, dashed=n40)")
    axes[1].invert_xaxis()
    axes[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage3_respiration_vs_theta.png"), dpi=150)
    plt.close(fig)


def plot_connectivity_vs_theta(sweep_data, results_dir):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    labels = ["largest_cluster_frac", "om_biomass_connectivity", "accessible_om_fraction"]
    titles = ["Largest water cluster fraction", "OM-biomass connectivity", "Accessible-OM fraction"]
    for ax, key, title in zip(axes, labels, titles):
        for grid_n in GRIDS:
            for soil in SOILS:
                vals = [sweep_data[grid_n][soil][th][1][key] for th in THETA_SWEEP]
                ax.plot(THETA_SWEEP, vals, GRID_MARKER[grid_n] + GRID_STYLE[grid_n],
                        label=f"Soil {soil} (n={grid_n})", color=COLORS[soil])
        ax.set_xlabel("theta")
        ax.set_ylabel(key)
        ax.set_title(title)
        ax.invert_xaxis()
        ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage3_connectivity_vs_theta.png"), dpi=150)
    plt.close(fig)


def plot_isolation_vs_theta(sweep_data, results_dir):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for grid_n in GRIDS:
        for soil in SOILS:
            iso = [sweep_data[grid_n][soil][th][1]["isolation_ratio"] for th in THETA_SWEEP]
            nclust = [sweep_data[grid_n][soil][th][1]["n_clusters"] for th in THETA_SWEEP]
            label = f"Soil {soil} (n={grid_n})"
            axes[0].plot(THETA_SWEEP, iso, GRID_MARKER[grid_n] + GRID_STYLE[grid_n],
                         label=label, color=COLORS[soil])
            axes[1].plot(THETA_SWEEP, nclust, GRID_MARKER[grid_n] + GRID_STYLE[grid_n],
                         label=label, color=COLORS[soil])
    axes[0].set_xlabel("theta")
    axes[0].set_ylabel("isolation_ratio = n_clusters / n_saturated_cells")
    axes[0].set_title("Isolation ratio vs theta")
    axes[0].invert_xaxis()
    axes[0].legend(fontsize=7)
    axes[1].set_xlabel("theta")
    axes[1].set_ylabel("n_isolated_clusters")
    axes[1].set_title("Raw isolated-cluster count vs theta")
    axes[1].invert_xaxis()
    axes[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage3_isolation_vs_theta.png"), dpi=150)
    plt.close(fig)


# ------------------------------------------------------- hypothesis test

def correlation_analysis(sweep_data, psds, results_dir):
    """soil_respiration_prompt_stage3.md SS9: does the isolation ratio
    predict respiration outcome better than a PSD summary stat does?"""
    psd_habitable = {s: psd_data.habitable_volume_fraction(psds[s]) for s in SOILS}
    psd_median = {s: psd_data.median_diameter_um(psds[s]) for s in SOILS}

    points = []
    for grid_n in GRIDS:
        for soil in SOILS:
            for theta in THETA_SWEEP:
                out, conn = sweep_data[grid_n][soil][theta]
                points.append({
                    "grid_n": grid_n, "soil": soil, "theta": theta,
                    "cum_co2": float(out["cum_co2"][-1]),
                    "R_peak": float(out["R_t"].max()),
                    "isolation_ratio": conn["isolation_ratio"],
                    "psd_habitable_frac": psd_habitable[soil],
                    "psd_median_diameter": psd_median[soil],
                })

    results = {}
    for outcome_key in ("cum_co2", "R_peak"):
        outcome_vals = np.array([p[outcome_key] for p in points])
        for pred_key in ("isolation_ratio", "psd_habitable_frac", "psd_median_diameter"):
            pred_vals = np.array([p[pred_key] for p in points])
            results[(outcome_key, pred_key)] = metrics.pearson_corr(outcome_vals, pred_vals)

    path = os.path.join(results_dir, "stage3_isolation_vs_psd_correlation.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["outcome", "predictor", "pearson_r", "n_points"])
        for (outcome_key, pred_key), r in results.items():
            w.writerow([outcome_key, pred_key, r, len(points)])

    plot_outcome_vs_predictor(points, "cum_co2", "isolation_ratio",
                               "cumulative CO2 (final)", "isolation_ratio",
                               "stage3_outcome_vs_isolation_ratio.png", results_dir)
    plot_outcome_vs_predictor(points, "cum_co2", "psd_habitable_frac",
                               "cumulative CO2 (final)", "PSD habitable fraction (30-150um)",
                               "stage3_outcome_vs_psd_stat.png", results_dir)

    return results, points


def plot_outcome_vs_predictor(points, outcome_key, pred_key, outcome_label,
                               pred_label, fname, results_dir):
    fig, ax = plt.subplots(figsize=(7, 5))
    for soil in SOILS:
        for grid_n in GRIDS:
            pts = [p for p in points if p["soil"] == soil and p["grid_n"] == grid_n]
            xs = [p[pred_key] for p in pts]
            ys = [p[outcome_key] for p in pts]
            ax.scatter(xs, ys, color=COLORS[soil], marker=GRID_MARKER[grid_n],
                       label=f"Soil {soil} (n={grid_n})")
    ax.set_xlabel(pred_label)
    ax.set_ylabel(outcome_label)
    r = metrics.pearson_corr(
        np.array([p[outcome_key] for p in points]),
        np.array([p[pred_key] for p in points]),
    )
    ax.set_title(f"{outcome_label} vs {pred_label} (pearson r={r:.3f})")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, fname), dpi=150)
    plt.close(fig)


def main():
    sweep_data, distinctness, (corr, points) = sweep()

    print("grid  theta   AB      AC      BC      min     success   outA outB outC")
    for grid_n in GRIDS:
        for theta in THETA_SWEEP:
            success, detail = distinctness[(grid_n, theta)]
            d = detail["distances"]
            oc = detail["outcomes"]
            print(f"{grid_n:4d}  {theta:5.2f}  {d['AB']:.3f}   {d['AC']:.3f}   {d['BC']:.3f}   "
                  f"{detail['min_distance']:.3f}   {success!s:5}   "
                  f"{oc['A']:8} {oc['B']:8} {oc['C']:8}")

    print("\ncorrelation of outcome with predictor (pearson r):")
    for (outcome_key, pred_key), r in corr.items():
        print(f"  {outcome_key:8} vs {pred_key:20} r={r:+.3f}")


if __name__ == "__main__":
    main()
