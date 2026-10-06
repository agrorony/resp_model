"""
Stage 5 (EXPLORATORY) runner
(prompts/soil_respiration_prompt_stage5_matric_inscribed.md).

Builds on the LATEST Variant C (dual_porosity.py, as evolved and run by
stage4_run.py: three phases macropore/porous-matrix/solid, per-soil OM
fractions, k_leak matrix->macropore leak, matrix internal diffusion, hydro
wicking, habitat-cluster K/D/S bonuses, S capacity) and changes ONLY two
things:

  Change 1 -- matric-potential wetting: replace the per-soil `global_theta`
              quantile with a SHARED air-entry diameter (`matric_d_cut_um`,
              default 50 um) -- matrix cells use their own PSD diameter
              (always wet, < 30 um), macropore cells use their REGION's
              inscribed diameter (Change 2), solid cells never wet. theta
              EMERGES per soil instead of being imposed.
  Change 2 -- geometric macropore sizing: retire the forced 80 um macropore
              cavity constant; label connected macropore regions and size
              each by its maximum inscribed circle (numpy-only iterative
              4-neighbor erosion distance -- dual_porosity.
              macropore_region_inscribed_diameter), fed into BOTH the
              wet/dry decision above and K(d)/D(d).

Why: Stage 4 showed Variant C's blocker was upstream of all its transport
machinery -- the water mask was still the per-soil global_theta quantile of
the FULL PSD, so Loess/Vertisol's habitat started bone dry (wet-habitat~0)
regardless of any leak/wicking/OM mechanics. This stage tests whether a
shared matric potential, combined with a geometrically realistic macropore
size, gets live wet habitat into those two soils.

A single matric condition is run (NOT a theta sweep -- theta emerges), plus
a 2-3-value matric_d_cut_um sensitivity (40/50/60), an n=40 grid check, and a
wicking on/off comparison (isolating wicking's own contribution from the
static matric map, since dynamic rewetting can blur Change 1's whole point).
`MAX_ITERATIONS=12` not exercised -- a defined run, not a search.

Usage:
    python stage5_run.py
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

RESULTS_DIR = os.path.join("results", "stage5")

GRID_N = 128
GRID_N_SMALL = 40                          # grid-robustness check

MATRIC_D_CUT_HEADLINE = 50.0                # mapping.yaml matric_d_cut_um
MATRIC_D_CUT_SENSITIVITY = [40.0, 50.0, 60.0]

MAX_ITERATIONS = 12   # not exercised -- a defined run, not a search


# ------------------------------------------------------------- run helpers

def _structure_cfg(soil, d_cut, grid_n):
    cfg = copy.deepcopy(load_yaml(f"configs/soil_{soil}.yaml"))
    cfg["saturation"] = {"mode": "matric_regions", "matric_d_cut_um": d_cut}
    cfg["grid"]["n"] = grid_n
    return cfg


def run_soil(soil, d_cut, grid_n, wicking_enabled, biology_cfg, mapping_cfg):
    """One (soil, matric_d_cut_um, grid, wicking on/off) run."""
    cfg = _structure_cfg(soil, d_cut, grid_n)
    out = dual_porosity.simulate_dual_porosity(
        biology_cfg, cfg, mapping_cfg, wicking_enabled=wicking_enabled)
    conn = metrics.connectivity_diagnostics(out["grids"])
    return out, conn


def region_stats(grids, d_cut):
    diam = grids["region_diam_um"]
    n_regions = int(diam.size)
    n_wet = int((diam < d_cut).sum())
    return {
        "n_regions": n_regions,
        "n_regions_wet": n_wet,
        "n_regions_dry": n_regions - n_wet,
        "mean_region_diam_um": float(diam.mean()) if n_regions else 0.0,
        "median_region_diam_um": float(np.median(diam)) if n_regions else 0.0,
        "max_region_diam_um": float(diam.max()) if n_regions else 0.0,
        "min_region_diam_um": float(diam.min()) if n_regions else 0.0,
    }


# --------------------------------------------------------------- 1. headline

def headline_sweep(biology_cfg, mapping_cfg):
    """Each soil, wicking ON and OFF, at the headline matric_d_cut_um, n=128."""
    data = {}
    for soil in SOILS:
        for wicking in (True, False):
            data[(soil, wicking)] = run_soil(
                soil, MATRIC_D_CUT_HEADLINE, GRID_N, wicking, biology_cfg, mapping_cfg)
    return data


def write_headline_csv(data, results_dir):
    rows = []
    for soil in SOILS:
        for wicking in (True, False):
            out, conn = data[(soil, wicking)]
            grids = out["grids"]
            pm = period_means(out["R_t"])
            row = {
                "soil": LABEL[soil], "wicking": wicking,
                "matric_d_cut_um": MATRIC_D_CUT_HEADLINE,
                "emergent_theta": grids["theta"],
                "wet_habitat_fraction": conn["wet_habitat_fraction"],
                "habitat_wet_share": conn["habitat_wet_share"],
                "isolation_ratio": conn["isolation_ratio"],
                "mean_S_habitat_final": metrics.mean_substrate_in_habitat(out["S_final"], grids["K"]),
                "mean_S_habitat_timeavg": float(out["S_habitat_mean_t"].mean()),
                "outcome": metrics.classify_outcome(out),
                "curve_shape": metrics.describe_shape(out["R_t"])["shape"],
                "P1": pm[0], "P2": pm[1], "P3": pm[2], "P4": pm[3],
                "period_shape_model": shape_of_periods(pm),
                "period_shape_experiment": shape_of_periods(EXPERIMENT[soil]),
                "matrix_fraction": grids.get("matrix_fraction"),
                "porous_matrix_fraction": grids.get("porous_matrix_fraction"),
                "solid_fraction": grids.get("solid_fraction"),
            }
            row.update(region_stats(grids, MATRIC_D_CUT_HEADLINE))
            rows.append(row)
    path = os.path.join(results_dir, "stage5_headline.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


def write_region_diameter_csv(data, results_dir):
    rows = []
    for soil in SOILS:
        out, _ = data[(soil, True)]
        grids = out["grids"]
        labels, sizes = dual_porosity._component_labels(grids["macro_mask"])
        diam = grids["region_diam_um"]
        for k in range(diam.size):
            rows.append({
                "soil": LABEL[soil], "region_id": k + 1,
                "size_cells": int(sizes[k]), "inscribed_diameter_um": float(diam[k]),
                "wet": bool(diam[k] < MATRIC_D_CUT_HEADLINE),
            })
    path = os.path.join(results_dir, "stage5_region_diameters.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


# ------------------------------------------------------ 2. matric sensitivity

def matric_sensitivity(biology_cfg, mapping_cfg, results_dir):
    rows = []
    for soil in SOILS:
        for d_cut in MATRIC_D_CUT_SENSITIVITY:
            out, conn = run_soil(soil, d_cut, GRID_N, True, biology_cfg, mapping_cfg)
            grids = out["grids"]
            pm = period_means(out["R_t"])
            row = {
                "soil": LABEL[soil], "matric_d_cut_um": d_cut,
                "emergent_theta": grids["theta"],
                "wet_habitat_fraction": conn["wet_habitat_fraction"],
                "mean_S_habitat_final": metrics.mean_substrate_in_habitat(out["S_final"], grids["K"]),
                "outcome": metrics.classify_outcome(out),
                "period_shape_model": shape_of_periods(pm),
                "period_shape_experiment": shape_of_periods(EXPERIMENT[soil]),
            }
            row.update(region_stats(grids, d_cut))
            rows.append(row)
    path = os.path.join(results_dir, "stage5_matric_sensitivity.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


# ------------------------------------------------------------ 3. grid check

def grid_check(biology_cfg, mapping_cfg, results_dir):
    rows = []
    for soil in SOILS:
        out, conn = run_soil(soil, MATRIC_D_CUT_HEADLINE, GRID_N_SMALL, True, biology_cfg, mapping_cfg)
        grids = out["grids"]
        pm = period_means(out["R_t"])
        row = {
            "soil": LABEL[soil], "grid_n": GRID_N_SMALL,
            "emergent_theta": grids["theta"],
            "wet_habitat_fraction": conn["wet_habitat_fraction"],
            "outcome": metrics.classify_outcome(out),
            "period_shape_model": shape_of_periods(pm),
            "period_shape_experiment": shape_of_periods(EXPERIMENT[soil]),
        }
        row.update(region_stats(grids, MATRIC_D_CUT_HEADLINE))
        rows.append(row)
    path = os.path.join(results_dir, "stage5_grid_check.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


# --------------------------------------------------------------- 4. figures

def plot_comparison_periods(data, results_dir):
    """The comparison figure: three soils' R(t) (period-binned) vs the
    measured P1-P4 pattern, wicking ON (solid) and OFF (dashed)."""
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8), sharey=False)
    x = np.arange(4)
    for ax, soil in zip(axes, SOILS):
        out_on, _ = data[(soil, True)]
        out_off, _ = data[(soil, False)]
        model_on = normalized_periods(period_means(out_on["R_t"]))
        model_off = normalized_periods(period_means(out_off["R_t"]))
        exp = normalized_periods(EXPERIMENT[soil])
        ax.plot(x, exp, "o-", color="black", lw=2, ms=7, label="experiment")
        ax.plot(x, model_on, "s--", color=COLORS[soil], lw=2, ms=7, label="model (wicking ON)")
        ax.plot(x, model_off, "^:", color=COLORS[soil], lw=1.6, ms=6, alpha=0.6,
                label="model (wicking OFF)")
        ax.axhline(1.0, color="gray", lw=0.6, ls=":")
        ax.set_xticks(x)
        ax.set_xticklabels(PERIOD_NAMES)
        ax.set_ylim(bottom=0)
        model_shape_on = shape_of_periods(period_means(out_on["R_t"]))
        exp_shape = shape_of_periods(EXPERIMENT[soil])
        match = "MATCH" if model_shape_on == exp_shape else "no match"
        ax.set_title(f"{LABEL[soil]}\nmodel(wick ON): {model_shape_on} | experiment: {exp_shape} ({match})",
                     fontsize=9.5)
        ax.set_xlabel("period (3/4/3/5 days)")
        ax.legend(fontsize=7.5)
    axes[0].set_ylabel("R, normalized to P1")
    fig.suptitle(f"Stage 5: matric-potential wetting + inscribed-circle macropores "
                 f"(matric_d_cut_um={MATRIC_D_CUT_HEADLINE:.0f}, n={GRID_N})")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage5_comparison_periods.png"), dpi=150)
    plt.close(fig)


def plot_region_map(soil, out, results_dir):
    grids = out["grids"]
    diam_field = grids["region_diam_field"]
    macro_mask = grids["macro_mask"]
    wet_dry = np.zeros(diam_field.shape)
    wet_dry[macro_mask & (diam_field < MATRIC_D_CUT_HEADLINE)] = 1.0   # wet macropore
    wet_dry[macro_mask & (diam_field >= MATRIC_D_CUT_HEADLINE)] = -1.0  # dry macropore
    # matrix stays 0 (background), solid gets its own marker
    wet_dry[grids["solid_mask"]] = -2.0

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    im0 = axes[0].imshow(np.log10(diam_field + 1e-6), cmap="viridis", origin="lower")
    axes[0].set_title("macropore region inscribed diameter (log10 um)\n(0 = matrix/solid)", fontsize=9.5)
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    im1 = axes[1].imshow(wet_dry, cmap="RdYlBu", origin="lower", vmin=-2, vmax=1)
    axes[1].set_title(f"wet(1,blue)/dry(-1,red) macropore, solid(-2,gray-ish)\n"
                       f"matric_d_cut_um={MATRIC_D_CUT_HEADLINE:.0f}", fontsize=9.5)
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    im2 = axes[2].imshow(grids["water_mask"].astype(float), cmap="Blues", origin="lower")
    axes[2].set_title(f"final water_mask (theta={grids['theta']:.3f})", fontsize=9.5)
    fig.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

    fig.suptitle(f"Stage 5 -- {LABEL[soil]}: {grids['region_diam_um'].size} macropore regions, "
                 f"{(grids['region_diam_um'] < MATRIC_D_CUT_HEADLINE).sum()} wet")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, f"stage5_region_map_{soil}.png"), dpi=150)
    plt.close(fig)


def plot_region_diameter_hist(data, results_dir):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), sharey=False)
    for ax, soil in zip(axes, SOILS):
        out, _ = data[(soil, True)]
        diam = out["grids"]["region_diam_um"]
        bins = np.logspace(np.log10(max(diam.min(), 1.0)), np.log10(diam.max() + 1), 20)
        ax.hist(diam, bins=bins, color=COLORS[soil], alpha=0.8)
        ax.axvline(MATRIC_D_CUT_HEADLINE, color="k", ls="--", lw=1.5,
                   label=f"matric_d_cut_um={MATRIC_D_CUT_HEADLINE:.0f}")
        ax.set_xscale("log")
        ax.set_xlabel("region inscribed diameter (um)")
        ax.set_title(f"{LABEL[soil]} ({diam.size} regions)")
        ax.legend(fontsize=8)
    axes[0].set_ylabel("region count")
    fig.suptitle("Stage 5: distribution of macropore-region inscribed diameters")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage5_region_diameter_hist.png"), dpi=150)
    plt.close(fig)


def plot_substrate_reaching_habitat(data, results_dir):
    fig, ax = plt.subplots(figsize=(8, 5))
    width = 0.35
    x = np.arange(len(SOILS))
    for i, wicking in enumerate((True, False)):
        vals = [metrics.mean_substrate_in_habitat(
                    data[(s, wicking)][0]["S_final"], data[(s, wicking)][0]["grids"]["K"])
                for s in SOILS]
        label = "wicking ON" if wicking else "wicking OFF"
        ax.bar(x + (i - 0.5) * width, vals, width, label=label,
               color=["#1f77b4", "#ff7f0e"][i])
    ax.set_xticks(x)
    ax.set_xticklabels([LABEL[s] for s in SOILS])
    ax.set_ylabel("mean substrate S in habitat cells (final)")
    ax.set_title("Stage 5: does carbon actually reach the habitat? (Stage-3.5/4 killer metric)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage5_substrate_reaching_habitat.png"), dpi=150)
    plt.close(fig)


def plot_matric_sensitivity(sens_rows, results_dir):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for soil in SOILS:
        soil_rows = [r for r in sens_rows if r["soil"] == LABEL[soil]]
        d_cuts = [r["matric_d_cut_um"] for r in soil_rows]
        wet_hab = [r["wet_habitat_fraction"] for r in soil_rows]
        theta = [r["emergent_theta"] for r in soil_rows]
        axes[0].plot(d_cuts, wet_hab, "o-", color=COLORS[soil], label=LABEL[soil])
        axes[1].plot(d_cuts, theta, "o-", color=COLORS[soil], label=LABEL[soil])
    axes[0].set_xlabel("matric_d_cut_um"); axes[0].set_ylabel("wet_habitat_fraction")
    axes[0].set_title("Wet-habitat fraction vs matric cutoff")
    axes[1].set_xlabel("matric_d_cut_um"); axes[1].set_ylabel("emergent theta")
    axes[1].set_title("Emergent theta vs matric cutoff")
    for ax in axes:
        ax.legend()
    fig.suptitle("Stage 5: matric_d_cut_um sensitivity (40/50/60 um), wicking ON, n=128")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage5_matric_sensitivity.png"), dpi=150)
    plt.close(fig)


# --------------------------------------------------------- 5. distinctness

def write_distinctness_csv(data, results_dir):
    rows = []
    for wicking in (True, False):
        outs = [data[(s, wicking)][0] for s in SOILS]
        success, detail = metrics.emergent_distinctness(*outs, require_nontrivial=False)
        d = detail["distances"]
        rows.append({
            "wicking": wicking,
            "Loess_vs_Sand": d["AB"], "Loess_vs_Vertisol": d["AC"], "Sand_vs_Vertisol": d["BC"],
            "min_distance": detail["min_distance"], "threshold": detail["threshold"],
            "success": success,
            "outcome_Loess": detail["outcomes"]["A"], "outcome_Sand": detail["outcomes"]["B"],
            "outcome_Vertisol": detail["outcomes"]["C"],
        })
    path = os.path.join(results_dir, "stage5_distinctness.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


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
    print("STAGE 5 (EXPLORATORY) -- matric-potential wetting + inscribed-circle macropores")
    print("=" * 78)

    print(f"\n[1] headline: matric_d_cut_um={MATRIC_D_CUT_HEADLINE:.0f}, n={GRID_N}, "
          f"wicking ON and OFF ...")
    data = headline_sweep(biology_cfg, mapping_cfg)
    rows = write_headline_csv(data, RESULTS_DIR)
    write_region_diameter_csv(data, RESULTS_DIR)

    print(f"    {'soil':10} {'wick':>5} {'theta':>7} {'wet-hab':>8} {'mean S hab':>11} "
          f"{'outcome':>10} {'model shp':>10} {'exp shp':>10} {'n_reg':>6} {'wet_reg':>8}")
    for r in rows:
        print(f"    {r['soil']:10} {str(r['wicking']):>5} {r['emergent_theta']:7.3f} "
              f"{r['wet_habitat_fraction']:8.3f} {r['mean_S_habitat_final']:11.5f} "
              f"{r['outcome']:>10} {r['period_shape_model']:>10} {r['period_shape_experiment']:>10} "
              f"{r['n_regions']:6d} {r['n_regions_wet']:8d}")

    for soil in SOILS:
        plot_region_map(soil, data[(soil, True)][0], RESULTS_DIR)

    plot_comparison_periods(data, RESULTS_DIR)
    plot_region_diameter_hist(data, RESULTS_DIR)
    plot_substrate_reaching_habitat(data, RESULTS_DIR)

    dist_rows = write_distinctness_csv(data, RESULTS_DIR)
    print(f"\n[2] distinctness (require_nontrivial=False):")
    for r in dist_rows:
        print(f"    wicking={r['wicking']}: success={r['success']}, min_distance={r['min_distance']:.4f}")

    print(f"\n[3] matric_d_cut_um sensitivity {MATRIC_D_CUT_SENSITIVITY} (wicking ON, n={GRID_N}) ...")
    sens_rows = matric_sensitivity(biology_cfg, mapping_cfg, RESULTS_DIR)
    plot_matric_sensitivity(sens_rows, RESULTS_DIR)
    for r in sens_rows:
        print(f"    {r['soil']:10} d_cut={r['matric_d_cut_um']:4.0f}  theta={r['emergent_theta']:.3f}  "
              f"wet-hab={r['wet_habitat_fraction']:.3f}  outcome={r['outcome']:>10}  "
              f"shape={r['period_shape_model']:>10}")

    print(f"\n[4] grid check (n={GRID_N_SMALL}, wicking ON) ...")
    grid_rows = grid_check(biology_cfg, mapping_cfg, RESULTS_DIR)
    for r in grid_rows:
        print(f"    {r['soil']:10} theta={r['emergent_theta']:.3f}  wet-hab={r['wet_habitat_fraction']:.3f}  "
              f"outcome={r['outcome']:>10}  shape={r['period_shape_model']:>10}")

    print(f"\n[5] all results written to {RESULTS_DIR}/")


if __name__ == "__main__":
    main()
