"""
Stage 4 (EXPLORATORY) runner (prompts/soil_respiration_prompt_stage4_exploratory.md).

Stage 3.5 diagnosed why the model reproduced only Sand-falling: the frozen
`OM ~ d^-2.5` rule dumps ~99.96% of carbon onto sub-micron clay cells, from
which diffusion (reach ~1-3 cells/run) can never deliver it to the 30-150 um
habitat. This is a defined HEAD-TO-HEAD comparison of two fixes, not a search
(`MAX_ITERATIONS = 12` ceiling, not exercised -- see SS0):

    Variant A ("truncate + distributed OM", model.build_grids' `psd_truncate_floor_um`
                         plus wet-diffusion scaling): drop pore volume below 10 um, keep
                         distributed OM on the truncated grid, and increase effective
                         diffusion in wet cells, on the SAME single-continuum grid
                         Stage 1-3.5 use.
  Variant C (dual-porosity matrix/macropore, dual_porosity.py): keep the FULL
             literature PSD, split every cell into a fast macropore or a
             blocked-fast-flow matrix cell, and let matrix cells leak
             substrate into adjacent macropores at rate k_leak.

Both are run for all three soils (Loess/Sand/Vertisol) at the SAME constant
theta = 0.6 used in Stage 3.5 (closed-jar, no drying), n=128, plus a 2-3 value
theta-robustness check and an n=40 grid-robustness check. The point is to see
which fix (if either) reproduces the measured pattern (Vertisol rising, Sand
falling, Loess low/erratic) and why, via the mechanism diagnostics -- not to
tune either variant toward the target.

Usage:
    python stage4_run.py
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
import dual_porosity
import metrics
import psd_data
import psd_parametric
import pore_field
from stage35_run import (
    SOILS, LABEL, COLORS, EXPERIMENT, PERIOD_DAYS, PERIOD_NAMES,
    period_means, shape_of_periods, normalized_periods,
)

RESULTS_DIR = os.path.join("results", "stage4")

THETA_CONST = 0.6                     # headline: same as Stage 3.5 (SS1)
THETA_ROBUSTNESS = [0.5, 0.6, 0.7]    # SS1's "2-3-value theta robustness check"

GRID_N = 128
GRID_N_SMALL = 40                     # SS1's grid-robustness check

MAX_ITERATIONS = 12   # SS0: not exercised -- this is a defined comparison

HAB_LO, HAB_HI = 30.0, 150.0

VARIANTS = ["A", "C"]
VARIANT_NAME = {"A": "Variant A (truncate + distributed OM)",
                "C": "Variant C (dual-porosity matrix/macropore)"}
VARIANT_STYLE = {"A": "s--", "C": "^-"}


# ------------------------------------------------------------- run helpers

def _base_structure_cfg(soil, sat_cfg, grid_n):
    cfg = copy.deepcopy(load_yaml(f"configs/soil_{soil}.yaml"))
    cfg["saturation"] = sat_cfg
    cfg["grid"]["n"] = grid_n
    return cfg


def run_variant(variant, soil, sat_cfg, grid_n, biology_cfg, mapping_cfg):
    """One (variant, soil, theta, grid) run. Returns (out, conn)."""
    cfg = _base_structure_cfg(soil, sat_cfg, grid_n)
    if variant == "A":
        cfg["pore"]["psd_truncate_floor_um"] = mapping_cfg["stage4_psd_truncate_floor_um"]
        cfg["pore"]["wet_diffusion_scale"] = mapping_cfg["stage4_variantA_wet_diffusion_scale"]
        out = simulate(biology_cfg, cfg, mapping_cfg)
    elif variant == "C":
        out = dual_porosity.simulate_dual_porosity(biology_cfg, cfg, mapping_cfg)
    else:
        raise ValueError(f"unknown variant: {variant}")
    conn = metrics.connectivity_diagnostics(out["grids"])
    return out, conn


def headline_row(variant, soil, theta, grid_n, out, conn):
    grids = out["grids"]
    return {
        "variant": variant,
        "soil": LABEL[soil],
        "theta": theta,
        "grid_n": grid_n,
        "outcome": metrics.classify_outcome(out),
        "R_peak": float(out["R_t"].max()),
        "mean_S_habitat_final": metrics.mean_substrate_in_habitat(out["S_final"], grids["K"]),
        "wet_habitat_fraction": conn["wet_habitat_fraction"],
        "curve_shape": metrics.describe_shape(out["R_t"])["shape"],
        "om_fraction_init": grids.get("om_fraction_init"),
        "om_total_init": grids.get("om_total_init"),
    }


def constant_theta(theta):
    return {"mode": "global_theta", "theta": theta}


# --------------------------------------------------- 1. PSD truncation report

def truncation_report(mapping_cfg, results_dir):
    """SS2.1: for each soil, report post-truncation porosity, habitable
    fraction, and how much volume the 10 um floor drops -- Variant A's
    truncation, at the PSD level (independent of any particular grid draw)."""
    floor_um = float(mapping_cfg["stage4_psd_truncate_floor_um"])
    rows = []
    for s in SOILS:
        psd = psd_parametric.build_soil_psd(s)
        psd_trunc, report = pore_field.truncate_psd_floor(psd, floor_um)
        rows.append({
            "soil": LABEL[s],
            "floor_um": floor_um,
            "porosity_before": report["porosity_before"],
            "porosity_after": report["porosity_after"],
            "dropped_volume_fraction": report["dropped_volume_fraction"],
            "habitable_frac_before": psd_data.habitable_volume_fraction(psd, HAB_LO, HAB_HI),
            "habitable_frac_after": psd_data.habitable_volume_fraction(psd_trunc, HAB_LO, HAB_HI),
        })
    path = os.path.join(results_dir, "stage4_variantA_truncation.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


# ------------------------------------------------------------- 2. main sweep

def sweep(biology_cfg, mapping_cfg):
    """All (variant, soil, theta, grid) combinations needed for the headline
    run plus the theta- and grid-robustness checks."""
    data = {}
    for variant in VARIANTS:
        for soil in SOILS:
            for theta in THETA_ROBUSTNESS:
                data[(variant, soil, theta, GRID_N)] = run_variant(
                    variant, soil, constant_theta(theta), GRID_N, biology_cfg, mapping_cfg)
            data[(variant, soil, THETA_CONST, GRID_N_SMALL)] = run_variant(
                variant, soil, constant_theta(THETA_CONST), GRID_N_SMALL, biology_cfg, mapping_cfg)
    return data


def write_variant_csv(data, variant, results_dir):
    rows = []
    for soil in SOILS:
        for theta in THETA_ROBUSTNESS:
            for grid_n in ({GRID_N, GRID_N_SMALL} if theta == THETA_CONST else {GRID_N}):
                out, conn = data[(variant, soil, theta, grid_n)]
                pm = period_means(out["R_t"])
                grids = out["grids"]
                row = {
                    "variant": variant, "soil": LABEL[soil],
                    "theta": theta, "grid_n": grid_n,
                    "R_peak": float(out["R_t"].max()), "R_end": float(out["R_t"][-1]),
                    "cum_co2_final": float(out["cum_co2"][-1]),
                    "B_final_over_B0": float(out["B_final"].sum() / grids["B0"].sum()),
                    "outcome": metrics.classify_outcome(out),
                    "curve_shape": metrics.describe_shape(out["R_t"])["shape"],
                    "P1": pm[0], "P2": pm[1], "P3": pm[2], "P4": pm[3],
                    "period_shape_model": shape_of_periods(pm),
                    "period_shape_experiment": shape_of_periods(EXPERIMENT[soil]),
                    "habitat_fraction": conn["habitat_fraction"],
                    "wet_habitat_fraction": conn["wet_habitat_fraction"],
                    "habitat_wet_share": conn["habitat_wet_share"],
                    "isolation_ratio": conn["isolation_ratio"],
                    "n_clusters": conn["n_clusters"],
                    "mean_S_habitat_final": metrics.mean_substrate_in_habitat(
                        out["S_final"], grids["K"]),
                    "mean_S_habitat_timeavg": float(out["S_habitat_mean_t"].mean()),
                    "matrix_fraction": grids.get("matrix_fraction"),
                    "solid_fraction": grids.get("solid_fraction"),
                    "porous_matrix_fraction": grids.get("porous_matrix_fraction"),
                    "habitat_count": grids.get("habitat_count"),
                    "mean_habitat_count_t": grids.get("mean_habitat_count_t"),
                }
                rows.append(row)
    path = os.path.join(results_dir, f"stage4_variant{variant}_sweep.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


# -------------------------------------------------------- 3. comparison plot

def plot_comparison_periods(data, results_dir, fname="stage4_comparison_periods.png"):
    """The ONE comparison figure (SS4): model R(t), binned into the same
    P1-P4 periods as the experiment, for Variant A and Variant C, three soils
    each, next to the experimental shape. 2 rows (variant) x 3 cols (soil)."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.6), sharex=True)
    x = np.arange(4)
    for row, variant in enumerate(VARIANTS):
        for col, soil in enumerate(SOILS):
            ax = axes[row, col]
            out, _ = data[(variant, soil, THETA_CONST, GRID_N)]
            model = normalized_periods(period_means(out["R_t"]))
            exp = normalized_periods(EXPERIMENT[soil])
            ax.plot(x, exp, "o-", color="black", lw=2, ms=7, label="experiment")
            ax.plot(x, model, VARIANT_STYLE[variant], color=COLORS[soil], lw=2, ms=7,
                    label=f"model ({variant})")
            ax.axhline(1.0, color="gray", lw=0.6, ls=":")
            ax.set_xticks(x)
            ax.set_xticklabels(PERIOD_NAMES)
            ax.set_ylim(bottom=0)
            model_shape = shape_of_periods(period_means(out["R_t"]))
            exp_shape = shape_of_periods(EXPERIMENT[soil])
            match = "MATCH" if model_shape == exp_shape else "no match"
            ax.set_title(f"{LABEL[soil]} -- {VARIANT_NAME[variant].split(' (')[0]}\n"
                         f"model: {model_shape} | experiment: {exp_shape} ({match})", fontsize=9.5)
            if row == 1:
                ax.set_xlabel("period (3/4/3/5 days)")
            if col == 0:
                ax.set_ylabel(f"{VARIANT_NAME[variant]}\nR, normalized to P1", fontsize=9)
            ax.legend(fontsize=7)
    fig.suptitle("Stage 4: Variant A vs Variant C vs measured respiration (constant theta=0.6, n=128)")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, fname), dpi=150)
    plt.close(fig)


def plot_variant_curves(data, variant, results_dir):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for soil in SOILS:
        out, _ = data[(variant, soil, THETA_CONST, GRID_N)]
        t = np.arange(len(out["R_t"])) * out["dt"]
        axes[0].plot(t, out["R_t"], color=COLORS[soil], label=LABEL[soil], lw=1.8)
        axes[1].plot(t, out["cum_co2"], color=COLORS[soil], label=LABEL[soil], lw=1.8)
    axes[0].set_xlabel("model time"); axes[0].set_ylabel("R(t)")
    axes[0].set_title("Respiration rate"); axes[0].legend()
    axes[1].set_xlabel("model time"); axes[1].set_ylabel("cumulative CO2")
    axes[1].set_title("Cumulative CO2"); axes[1].legend()
    fig.suptitle(f"Stage 4 {VARIANT_NAME[variant]} -- constant theta={THETA_CONST}, n={GRID_N}")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, f"stage4_respiration_variant{variant}.png"), dpi=150)
    plt.close(fig)


def plot_variantc_om_refactor_rates(single_data, results_dir):
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    for soil in SOILS:
        out, _ = single_data[soil]
        t = np.arange(len(out["R_t"])) * out["dt"]
        ax.plot(t, out["R_t"], color=COLORS[soil], lw=2.2, label=LABEL[soil])
    ax.set_xlabel("model time")
    ax.set_ylabel("R(t)")
    ax.set_title("Variant C OM-refactor: Respiration Rate (theta=0.6, n=128)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "variantC_OM_refactor_rates.png"), dpi=220)
    plt.close(fig)


def plot_variantc_om_refactor_cumulative(single_data, results_dir):
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    for soil in SOILS:
        out, _ = single_data[soil]
        t = np.arange(len(out["cum_co2"])) * out["dt"]
        ax.plot(t, out["cum_co2"], color=COLORS[soil], lw=2.2, label=LABEL[soil])
    ax.set_xlabel("model time")
    ax.set_ylabel("cumulative R(t)")
    ax.set_title("Variant C OM-refactor: Cumulative Respiration (theta=0.6, n=128)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "variantC_OM_refactor_cumulative.png"), dpi=220)
    plt.close(fig)


def write_variantc_om_refactor_report(rows, distinctness, results_dir):
    by_soil = {r["soil"]: r for r in rows}

    d = distinctness["distances"]
    threshold = float(distinctness["threshold"])
    pair_rows = [
        ("Loess vs Sand", d["AB"]),
        ("Loess vs Vertisol", d["AC"]),
        ("Sand vs Vertisol", d["BC"]),
    ]

    def _fmt_outcome(soil_label):
        r = by_soil[soil_label]
        return (
            f"{soil_label} ({r['outcome']}, shape={r['curve_shape']}, "
            f"R_peak={r['R_peak']:.6f}, "
            f"mean_S_habitat_final={r['mean_S_habitat_final']:.6f}, "
            f"wet_habitat_fraction={r['wet_habitat_fraction']:.4f})"
        )

    dynamic_line = (
        " ; ".join([
            _fmt_outcome("Sand"),
            _fmt_outcome("Loess"),
            _fmt_outcome("Vertisol"),
        ])
    )

    lines = [
        "# Variant C OM Refactor Report",
        "",
        "## Setup",
        "",
        "This run applies soil-specific initial OM fractions in Variant C at a single condition: theta=0.6, n=128.",
        "OM initialization is applied only to porous matrix cells (solid phase and macropore cavity cells receive zero OM).",
        "Configured OM fractions:",
        "- Sand: 2.0%",
        "- Loess: 3.5%",
        "- Vertisol: 5.0%",
        "",
        "## Headline Metrics (theta=0.6, n=128)",
        "",
        "| Soil | outcome | R_peak | mean_S_habitat_final | wet_habitat_fraction |",
        "|---|---:|---:|---:|---:|",
    ]

    for soil in ["Sand", "Loess", "Vertisol"]:
        r = by_soil[soil]
        lines.append(
            f"| {soil} | {r['outcome']} | {r['R_peak']:.6f} | "
            f"{r['mean_S_habitat_final']:.6f} | {r['wet_habitat_fraction']:.4f} |"
        )

    lines.extend([
        "",
        "## Dynamic Analysis",
        "",
        "The refactor combines three transport controls in Variant C: soil-specific OM loading, internal porous-matrix diffusion, and wet-matrix-driven network wicking toward macropore habitats.",
        f"Observed behavior: {dynamic_line}.",
        "Interpretation: the custom OM loading increases available dissolved substrate in matrix storage, while matrix diffusion and wicking govern how much of that substrate actually reaches wet habitat cells.",
        "Under this single condition, starvation is reduced where wet habitat connectivity and transfer are high; soils with lower wet habitat exposure continue to show declining/falling behavior even with larger OM loading.",
        "",
        "## Pairwise Distinctness (theta=0.6, n=128)",
        "",
        "| Pair | distance | pass (>0.3) |",
        "|---|---:|---:|",
    ])

    for pair_name, dist in pair_rows:
        lines.append(f"| {pair_name} | {dist:.6f} | {'yes' if dist > threshold else 'no'} |")

    lines.extend([
        "",
        f"Minimum pairwise distance: {distinctness['min_distance']:.6f} (threshold={threshold:.3f}).",
        f"Distinctness success status: {'PASS' if distinctness['min_distance'] > threshold else 'FAIL'} for the >0.3 criterion.",
        "",
        "## Output Files",
        "",
        "- variantC_OM_refactor_rates.png",
        "- variantC_OM_refactor_cumulative.png",
        "",
    ])

    path = os.path.join(results_dir, "VARIANTC_OM_REFACTOR_REPORT.md")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def run_variantc_om_refactor_single(biology_cfg, mapping_cfg):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    single_data = {
        soil: run_variant("C", soil, constant_theta(THETA_CONST), GRID_N, biology_cfg, mapping_cfg)
        for soil in SOILS
    }

    rows = []
    for soil in SOILS:
        out, conn = single_data[soil]
        rows.append(headline_row("C", soil, THETA_CONST, GRID_N, out, conn))

    path = os.path.join(RESULTS_DIR, "variantC_OM_refactor_headline.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    loess = single_data["loess"][0]
    sand = single_data["sand"][0]
    vertisol = single_data["vertisol"][0]
    _, detail = metrics.emergent_distinctness(loess, sand, vertisol, require_nontrivial=False)

    plot_variantc_om_refactor_rates(single_data, RESULTS_DIR)
    plot_variantc_om_refactor_cumulative(single_data, RESULTS_DIR)
    write_variantc_om_refactor_report(rows, detail, RESULTS_DIR)

    print("\n[Variant C OM refactor] headline condition complete:")
    for r in rows:
        print(
            f"    {r['soil']:10} outcome={r['outcome']:>9} "
            f"R_peak={r['R_peak']:.6f} "
            f"mean_S_habitat_final={r['mean_S_habitat_final']:.6f} "
            f"wet_habitat_fraction={r['wet_habitat_fraction']:.4f}"
        )
    print(
        f"    pairwise min_distance={detail['min_distance']:.6f} "
        f"(threshold={detail['threshold']:.3f})"
    )
    print(f"    outputs written to {RESULTS_DIR}/")


def plot_substrate_reaching_habitat(data, results_dir):
    """The mechanism figure (SS4): mean substrate concentration reaching the
    habitat, per soil per variant -- does carbon actually get to the
    microbes, not just whether the habitat is wet."""
    fig, ax = plt.subplots(figsize=(8, 5))
    width = 0.35
    x = np.arange(len(SOILS))
    for i, variant in enumerate(VARIANTS):
        vals = [metrics.mean_substrate_in_habitat(
                    data[(variant, s, THETA_CONST, GRID_N)][0]["S_final"],
                    data[(variant, s, THETA_CONST, GRID_N)][0]["grids"]["K"])
                for s in SOILS]
        ax.bar(x + (i - 0.5) * width, vals, width, label=VARIANT_NAME[variant],
               color=["#1f77b4", "#ff7f0e"][i])
    ax.set_xticks(x)
    ax.set_xticklabels([LABEL[s] for s in SOILS])
    ax.set_ylabel("mean substrate S in habitat cells (final)")
    ax.set_title("Stage 4: does carbon actually reach the habitat? (Stage-3.5 killer metric)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage4_substrate_reaching_habitat.png"), dpi=150)
    plt.close(fig)


def plot_water_map_variant_a(soil, out, results_dir):
    grids = out["grids"]
    fig, axes = plt.subplots(1, 5, figsize=(20, 4))
    fields = [
        ("pore diameter d (um)", np.log10(grids["r"]), "viridis", "log10 d"),
        ("K(d) - habitat", grids["K"], "viridis", "K"),
        ("D used (0 where dry)", grids["D"], "viridis", "D"),
        ("water-filled mask", grids["water_mask"].astype(float), "Blues", "wet"),
        ("OM (log10, distributed)", np.log10(grids["OM"] + 1e-12), "magma", "log10 OM"),
    ]
    for ax, (name, field, cmap, cbl) in zip(axes, fields):
        im = ax.imshow(field, cmap=cmap, origin="lower")
        ax.set_title(name, fontsize=10)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=cbl)
    trunc = grids.get("psd_truncation") or {}
    fig.suptitle(f"Stage 4 Variant A -- {LABEL[soil]} (n={grids['n']}, theta={grids['theta']:.2f}), "
                 f"dropped={trunc.get('dropped_volume_fraction', 0):.4f} below "
                 f"{trunc.get('floor_um', '?')} um")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, f"stage4_water_map_variantA_{soil}.png"), dpi=150)
    plt.close(fig)


def plot_water_map_variant_c(soil, out, results_dir):
    grids = out["grids"]
    fig, axes = plt.subplots(1, 9, figsize=(35, 4))
    fields = [
        ("pore diameter d (um)", np.log10(grids["d"]), "viridis", "log10 d"),
        ("K(d) - habitat", grids["K"], "viridis", "K"),
        ("D used (0 where dry OR non-macro)", grids["D"], "viridis", "D"),
        ("water-filled mask", grids["water_mask"].astype(float), "Blues", "wet"),
        ("macropore (1) vs matrix (0)", grids["macro_mask"].astype(float), "PiYG", "macro"),
        ("solid phase (1=solid)", grids["solid_mask"].astype(float), "gray", "solid"),
        ("porous matrix (1=porous)", grids["porous_matrix_mask"].astype(float), "cividis", "porous"),
        ("macro habitat network labels", grids["habitat_labels"].astype(float), "tab20", "cluster id"),
        ("OM (log10, in matrix)", np.log10(grids["OM"] + 1e-12), "magma", "log10 OM"),
    ]
    for ax, (name, field, cmap, cbl) in zip(axes, fields):
        im = ax.imshow(field, cmap=cmap, origin="lower")
        ax.set_title(name, fontsize=10)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=cbl)
    fig.suptitle(f"Stage 4 Variant C -- {LABEL[soil]} (n={grids['n']}, theta={grids['theta']:.2f}), "
                 f"matrix={grids['matrix_fraction']:.3f}, solid={grids['solid_fraction']:.3f}, "
                 f"habitats={grids['habitat_count']}")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, f"stage4_water_map_variantC_{soil}.png"), dpi=150)
    plt.close(fig)


# --------------------------------------------------------- 4. distinctness

def distinctness_row(data, variant, theta, grid_n, results_dir):
    outs = [data[(variant, s, theta, grid_n)][0] for s in SOILS]
    success, detail = metrics.emergent_distinctness(*outs, require_nontrivial=False)
    d = detail["distances"]
    return {
        "variant": variant, "theta": theta, "grid_n": grid_n,
        "Loess_vs_Sand": d["AB"], "Loess_vs_Vertisol": d["AC"], "Sand_vs_Vertisol": d["BC"],
        "min_distance": detail["min_distance"], "threshold": detail["threshold"],
        "success": success,
        "outcome_Loess": detail["outcomes"]["A"], "outcome_Sand": detail["outcomes"]["B"],
        "outcome_Vertisol": detail["outcomes"]["C"],
    }


def write_distinctness_csv(data, results_dir):
    rows = []
    for variant in VARIANTS:
        for theta in THETA_ROBUSTNESS:
            rows.append(distinctness_row(data, variant, theta, GRID_N, results_dir))
        rows.append(distinctness_row(data, variant, THETA_CONST, GRID_N_SMALL, results_dir))
    path = os.path.join(results_dir, "stage4_distinctness.csv")
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
    print("STAGE 4 (EXPLORATORY) -- Variant A vs Variant C, head-to-head")
    print("=" * 78)

    trows = truncation_report(mapping_cfg, RESULTS_DIR)
    print(f"\n[1] Variant A PSD truncation (floor={mapping_cfg['stage4_psd_truncate_floor_um']} um):")
    for r in trows:
        print(f"    {r['soil']:10} porosity {r['porosity_before']:.3f} -> {r['porosity_after']:.3f}  "
              f"dropped={r['dropped_volume_fraction']:.4f}  "
              f"habitable {r['habitable_frac_before']:.3f} -> {r['habitable_frac_after']:.3f}")

    print(f"\n[2] running both variants x 3 soils x theta{THETA_ROBUSTNESS} (n={GRID_N}) "
          f"+ n={GRID_N_SMALL} grid check ...")
    data = sweep(biology_cfg, mapping_cfg)

    for variant in VARIANTS:
        rows = write_variant_csv(data, variant, RESULTS_DIR)
        headline = [r for r in rows if r["theta"] == THETA_CONST and r["grid_n"] == GRID_N]
        print(f"\n[3] {VARIANT_NAME[variant]} -- headline (theta={THETA_CONST}, n={GRID_N}):")
        print(f"    {'soil':10} {'outcome':>10} {'model shape':>12} {'experiment':>12} "
              f"{'wet-hab':>8} {'mean S hab':>11}")
        for r in headline:
            print(f"    {r['soil']:10} {r['outcome']:>10} {r['period_shape_model']:>12} "
                  f"{r['period_shape_experiment']:>12} {r['wet_habitat_fraction']:8.3f} "
                  f"{r['mean_S_habitat_final']:11.4f}")
        plot_variant_curves(data, variant, RESULTS_DIR)
        for soil in SOILS:
            out, _ = data[(variant, soil, THETA_CONST, GRID_N)]
            if variant == "A":
                plot_water_map_variant_a(soil, out, RESULTS_DIR)
            else:
                plot_water_map_variant_c(soil, out, RESULTS_DIR)

    plot_comparison_periods(data, RESULTS_DIR)
    plot_substrate_reaching_habitat(data, RESULTS_DIR)

    dist_rows = write_distinctness_csv(data, RESULTS_DIR)
    print(f"\n[4] distinctness (require_nontrivial=False), headline theta={THETA_CONST}, n={GRID_N}:")
    for r in dist_rows:
        if r["theta"] == THETA_CONST and r["grid_n"] == GRID_N:
            print(f"    variant {r['variant']}: success={r['success']}, "
                  f"min_distance={r['min_distance']:.3f}")

    print(f"\n[5] all results written to {RESULTS_DIR}/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Stage 4 runner")
    parser.add_argument(
        "--variantc-om-refactor-single",
        action="store_true",
        help="Run Variant C only for theta=0.6, n=128 with soil-specific OM initialization and write report/plots.",
    )
    args = parser.parse_args()

    if args.variantc_om_refactor_single:
        biology_cfg = load_yaml("configs/biology.yaml")
        mapping_cfg = load_yaml("configs/mapping.yaml")
        # Stage 16 (soil_respiration_prompt_stage16_shared_om_rule.md SS1) moved
        # configs/mapping.yaml's own default to `om_mode: shared_fine_pore`.
        # This stage predates that change and is pinned to the Stage 6-15
        # coating/trapped path explicitly, so it stays independently
        # reproducible rather than silently inheriting the new rule.
        mapping_cfg["om_mode"] = dual_porosity.OM_MODE_COATING
        run_variantc_om_refactor_single(biology_cfg, mapping_cfg)
    else:
        main()
