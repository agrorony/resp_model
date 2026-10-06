"""
Stage 6 (EXPLORATORY) runner
(prompts/soil_respiration_prompt_stage6_portback.md).

Ports two conclusions proved in the abstract PDE sandbox (`pde_sandbox/`,
its own RESULTS.md) into the real dual-porosity soil model, as left after
Stage 5:

  Change 1 -- implicit backward-Euler transport for the dissolved-substrate
              field (`dual_porosity.implicit_diffusion_step`), replacing the
              explicit conservative-flux update and its `max(D)*dt/dx^2<0.2`
              stability cap. Critically, this ALSO unifies the macropore and
              matrix D fields into one connected solve -- Stage 4/5 ran two
              separate explicit divergences (D_macro, D_matrix) that never
              overlapped, so matrix and macropore cells could never diffuse
              into each other directly at all; only the special-cased
              `k_leak` term crossed that boundary. Stage 6 fixes that too.
  Change 2/3 -- OM split into a fast COATING pool (near habitat, along a wet
              path, dissolves at k_dis) and a slow TRAPPED pool (interior,
              dissolves at k_dis_slow, a tail not a lock), per-soil coating
              fraction (the hypothesis-motivated input).
  Change 4 -- S_cap removed; wicking off by default (dual_porosity.py
              defaults now match this).

Diagnostics per soil: mean substrate reaching the habitat over time (the
Stage 3.5/4/5 killer metric -- this run's direct check of whether the port
actually helps), habitats-recruited-over-time, outcome class, R(t) shape,
pairwise distinctness. Target (qualitative): Sand sharp rise-then-fall,
Vertisol rising, Loess low/erratic. `MAX_ITERATIONS=12` not exercised.

Usage:
    python stage6_run.py
"""
from __future__ import annotations

import copy
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from model import load_yaml
import dual_porosity
import metrics
from stage35_run import (
    SOILS, LABEL, COLORS, EXPERIMENT, PERIOD_NAMES,
    period_means, shape_of_periods, normalized_periods,
)

RESULTS_DIR = os.path.join("results", "stage6")

GRID_N = 128
GRID_N_SMALL = 40

MATRIC_D_CUT = 50.0   # Stage 5's matric potential, kept unchanged (SS1)

MAX_ITERATIONS = 12   # not exercised -- a defined port + run, not a search


# ------------------------------------------------------------- run helpers

def _structure_cfg(soil, grid_n):
    cfg = copy.deepcopy(load_yaml(f"configs/soil_{soil}.yaml"))
    cfg["saturation"] = {"mode": "matric_regions", "matric_d_cut_um": MATRIC_D_CUT}
    cfg["grid"]["n"] = grid_n
    return cfg


def run_soil(soil, grid_n, biology_cfg, mapping_cfg):
    cfg = _structure_cfg(soil, grid_n)
    out = dual_porosity.simulate_dual_porosity(biology_cfg, cfg, mapping_cfg, wicking_enabled=False)
    conn = metrics.connectivity_diagnostics(out["grids"])
    return out, conn


# --------------------------------------------------------------- 1. headline

def headline_run(biology_cfg, mapping_cfg):
    return {soil: run_soil(soil, GRID_N, biology_cfg, mapping_cfg) for soil in SOILS}


def write_headline_csv(data, results_dir):
    rows = []
    for soil in SOILS:
        out, conn = data[soil]
        grids = out["grids"]
        pm = period_means(out["R_t"])
        rows.append({
            "soil": LABEL[soil],
            "emergent_theta": grids["theta"],
            "wet_habitat_fraction": conn["wet_habitat_fraction"],
            "habitat_wet_share": conn["habitat_wet_share"],
            "mean_S_habitat_final": metrics.mean_substrate_in_habitat(out["S_final"], grids["K"]),
            "mean_S_habitat_timeavg": float(out["S_habitat_mean_t"].mean()),
            "habitats_recruited_final": float(out["habitats_recruited_t"][-1]),
            "habitats_recruited_max": float(out["habitats_recruited_t"].max()),
            "habitats_recruited_mean": float(out["habitats_recruited_t"].mean()),
            "n_regions": int(grids["region_diam_um"].size),
            "outcome": metrics.classify_outcome(out),
            "curve_shape": metrics.describe_shape(out["R_t"])["shape"],
            "R_peak": float(out["R_t"].max()), "R_end": float(out["R_t"][-1]),
            "cum_co2_final": float(out["cum_co2"][-1]),
            "B_final_over_B0": float(out["B_final"].sum() / grids["B0"].sum()),
            "P1": pm[0], "P2": pm[1], "P3": pm[2], "P4": pm[3],
            "period_shape_model": shape_of_periods(pm),
            "period_shape_experiment": shape_of_periods(EXPERIMENT[soil]),
            "coating_fraction": grids["coating_fraction"],
            "om_total_init": grids["om_total_init"],
            "coating_cells": int(grids["coating_mask"].sum()),
            "trapped_cells": int(grids["trapped_mask"].sum()),
        })
    path = os.path.join(results_dir, "stage6_headline.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


def grid_check(biology_cfg, mapping_cfg, results_dir):
    rows = []
    for soil in SOILS:
        out, conn = run_soil(soil, GRID_N_SMALL, biology_cfg, mapping_cfg)
        grids = out["grids"]
        pm = period_means(out["R_t"])
        rows.append({
            "soil": LABEL[soil], "grid_n": GRID_N_SMALL,
            "emergent_theta": grids["theta"],
            "wet_habitat_fraction": conn["wet_habitat_fraction"],
            "mean_S_habitat_final": metrics.mean_substrate_in_habitat(out["S_final"], grids["K"]),
            "outcome": metrics.classify_outcome(out),
            "period_shape_model": shape_of_periods(pm),
            "period_shape_experiment": shape_of_periods(EXPERIMENT[soil]),
        })
    path = os.path.join(results_dir, "stage6_grid_check.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


# --------------------------------------------------------------- 2. figures

def plot_comparison_periods(data, results_dir):
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
    x = np.arange(4)
    for ax, soil in zip(axes, SOILS):
        out, _ = data[soil]
        model = normalized_periods(period_means(out["R_t"]))
        exp = normalized_periods(EXPERIMENT[soil])
        ax.plot(x, exp, "o-", color="black", lw=2, ms=7, label="experiment")
        ax.plot(x, model, "s--", color=COLORS[soil], lw=2, ms=7, label="model (Stage 6)")
        ax.axhline(1.0, color="gray", lw=0.6, ls=":")
        ax.set_xticks(x); ax.set_xticklabels(PERIOD_NAMES)
        ax.set_ylim(bottom=0)
        model_shape = shape_of_periods(period_means(out["R_t"]))
        exp_shape = shape_of_periods(EXPERIMENT[soil])
        match = "MATCH" if model_shape == exp_shape else "no match"
        ax.set_title(f"{LABEL[soil]}\nmodel: {model_shape} | experiment: {exp_shape} ({match})", fontsize=10)
        ax.set_xlabel("period (3/4/3/5 days)")
        ax.legend(fontsize=8)
    axes[0].set_ylabel("R, normalized to P1")
    fig.suptitle("Stage 6: implicit transport + coating/trapped OM vs measured respiration "
                 f"(matric_d_cut_um={MATRIC_D_CUT:.0f}, n={GRID_N}, wicking OFF)")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage6_comparison_periods.png"), dpi=150)
    plt.close(fig)


def plot_respiration_curves(data, results_dir):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for soil in SOILS:
        out, _ = data[soil]
        t = np.arange(len(out["R_t"])) * out["dt"]
        axes[0].plot(t, out["R_t"], color=COLORS[soil], label=LABEL[soil], lw=1.8)
        axes[1].plot(t, out["cum_co2"], color=COLORS[soil], label=LABEL[soil], lw=1.8)
    axes[0].set_xlabel("model time"); axes[0].set_ylabel("R(t)")
    axes[0].set_title("Respiration rate"); axes[0].legend()
    axes[1].set_xlabel("model time"); axes[1].set_ylabel("cumulative CO2")
    axes[1].set_title("Cumulative CO2"); axes[1].legend()
    fig.suptitle(f"Stage 6 -- constant moisture (matric_d_cut_um={MATRIC_D_CUT:.0f}), n={GRID_N}")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage6_respiration_curves.png"), dpi=150)
    plt.close(fig)


def plot_mean_s_habitat_vs_t(data, results_dir):
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    for soil in SOILS:
        out, _ = data[soil]
        t = np.arange(len(out["S_habitat_mean_t"])) * out["dt"]
        ax.plot(t, out["S_habitat_mean_t"], color=COLORS[soil], label=LABEL[soil], lw=1.8)
    ax.axhline(0.2, color="k", ls="--", lw=1, label="Ks=0.2 (biology half-saturation)")
    ax.set_xlabel("model time")
    ax.set_ylabel("mean substrate S in habitat cells")
    ax.set_title("Stage 6 killer metric: mean substrate reaching the habitat over time")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage6_mean_s_habitat_vs_t.png"), dpi=150)
    plt.close(fig)


def plot_habitats_recruited_vs_t(data, results_dir):
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    for soil in SOILS:
        out, _ = data[soil]
        t = np.arange(len(out["habitats_recruited_t"])) * out["dt"]
        ax.plot(t, out["habitats_recruited_t"], color=COLORS[soil], label=LABEL[soil], lw=1.8)
    ax.set_xlabel("model time")
    ax.set_ylabel("count of actively-growing macropore regions")
    ax.set_title("Stage 6: habitats recruited over time")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage6_habitats_recruited_vs_t.png"), dpi=150)
    plt.close(fig)


def plot_om_maps(soil, out, results_dir):
    grids = out["grids"]
    fig, axes = plt.subplots(1, 4, figsize=(19, 4.4))
    fields = [
        ("macropore(1)/matrix(0)/solid(-1)",
         np.where(grids["solid_mask"], -1.0, grids["macro_mask"].astype(float)), "PiYG"),
        ("coating mask (1=coating)", grids["coating_mask"].astype(float), "YlOrRd"),
        ("OM_coating (log10)", np.log10(grids["OM_coating"] + 1e-12), "magma"),
        ("OM_trapped (log10)", np.log10(grids["OM_trapped"] + 1e-12), "magma"),
    ]
    for ax, (name, field, cmap) in zip(axes, fields):
        im = ax.imshow(field, cmap=cmap, origin="lower")
        ax.set_title(name, fontsize=10)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle(f"Stage 6 -- {LABEL[soil]}: coating_fraction={grids['coating_fraction']:.2f}, "
                 f"coating cells={int(grids['coating_mask'].sum())}, "
                 f"trapped cells={int(grids['trapped_mask'].sum())}")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, f"stage6_om_map_{soil}.png"), dpi=150)
    plt.close(fig)


# --------------------------------------------------------- 3. distinctness

def write_distinctness_csv(data, results_dir):
    outs = [data[s][0] for s in SOILS]
    success, detail = metrics.emergent_distinctness(*outs, require_nontrivial=False)
    d = detail["distances"]
    row = {
        "Loess_vs_Sand": d["AB"], "Loess_vs_Vertisol": d["AC"], "Sand_vs_Vertisol": d["BC"],
        "min_distance": detail["min_distance"], "threshold": detail["threshold"],
        "success": success,
        "outcome_Loess": detail["outcomes"]["A"], "outcome_Sand": detail["outcomes"]["B"],
        "outcome_Vertisol": detail["outcomes"]["C"],
    }
    path = os.path.join(results_dir, "stage6_distinctness.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(row))
        w.writeheader()
        w.writerow(row)
    return row


# ------------------------------------------------------------------ main

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    biology_cfg = load_yaml("configs/biology.yaml")
    mapping_cfg = load_yaml("configs/mapping.yaml")
    # Stage 16 (soil_respiration_prompt_stage16_shared_om_rule.md SS1) moved
    # configs/mapping.yaml's own default to `om_mode: shared_fine_pore`.
    # This stage predates that change and is pinned to the Stage 6-15
    # coating/trapped path explicitly, so it stays independently
    # reproducible rather than silently inheriting the new rule.
    mapping_cfg["om_mode"] = dual_porosity.OM_MODE_COATING

    print("=" * 78)
    print("STAGE 6 (EXPLORATORY) -- port implicit transport + coating/trapped OM into Variant C")
    print("=" * 78)

    print(f"\n[1] headline: matric_d_cut_um={MATRIC_D_CUT:.0f}, n={GRID_N}, wicking OFF ...")
    data = headline_run(biology_cfg, mapping_cfg)
    rows = write_headline_csv(data, RESULTS_DIR)

    print(f"    {'soil':10} {'theta':>7} {'wet-hab':>8} {'mean S hab':>11} {'hab.recruit':>11} "
          f"{'outcome':>10} {'model shp':>10} {'exp shp':>10}")
    for r in rows:
        print(f"    {r['soil']:10} {r['emergent_theta']:7.3f} {r['wet_habitat_fraction']:8.3f} "
              f"{r['mean_S_habitat_final']:11.6f} {r['habitats_recruited_final']:11.0f} "
              f"{r['outcome']:>10} {r['period_shape_model']:>10} {r['period_shape_experiment']:>10}")

    for soil in SOILS:
        plot_om_maps(soil, data[soil][0], RESULTS_DIR)

    plot_comparison_periods(data, RESULTS_DIR)
    plot_respiration_curves(data, RESULTS_DIR)
    plot_mean_s_habitat_vs_t(data, RESULTS_DIR)
    plot_habitats_recruited_vs_t(data, RESULTS_DIR)

    dist_row = write_distinctness_csv(data, RESULTS_DIR)
    print(f"\n[2] distinctness (require_nontrivial=False): success={dist_row['success']}, "
          f"min_distance={dist_row['min_distance']:.4f}")

    print(f"\n[3] grid check (n={GRID_N_SMALL}) ...")
    grid_rows = grid_check(biology_cfg, mapping_cfg, RESULTS_DIR)
    for r in grid_rows:
        print(f"    {r['soil']:10} theta={r['emergent_theta']:.3f}  wet-hab={r['wet_habitat_fraction']:.3f}  "
              f"mean_S_hab={r['mean_S_habitat_final']:.6f}  outcome={r['outcome']:>10}  "
              f"shape={r['period_shape_model']:>10}")

    print(f"\n[4] all results written to {RESULTS_DIR}/")


if __name__ == "__main__":
    main()
