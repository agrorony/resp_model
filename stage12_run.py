"""
Stage 12 runner (prompts/soil_respiration_prompt_stage12_significance.md).

An ANALYSIS stage, not a modeling stage. Builds on Stage 11: the two soils
(Sand, Vertisol) and their Stage-11 parameters (`om_total`, coating
fraction, `D_scale`, biology, `T`) are read from `configs/mapping.yaml`
EXACTLY as Stage 11 left them -- never re-tuned here. For each soil, an
ensemble of M random pore-field seeds is run (same structural class, same
everything except the seed that drives the pore-field realization), to
show statistically that each soil's characteristic respiration shape is
driven by its defined STRUCTURAL FEATURES (few concentrated coating-fed
regions for Sand; many distributed regions for Vertisol), not by a lucky
random layout. numpy-only (no scipy) -- Mann-Whitney U, effect size, and
bootstrap CIs are all implemented directly.

Usage:
    python stage12_run.py
"""
from __future__ import annotations

import copy
import csv
import datetime
import math
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from model import load_yaml
import dual_porosity
import metrics

SOILS = ["sand", "vertisol"]
LABEL = {"sand": "Sand", "vertisol": "Vertisol"}
COLORS = {"sand": "#d62728", "vertisol": "#4b0082"}

RESULTS_DIR = os.path.join("results", "stage12")
LOGBOOK_PATH = "LOGBOOK.md"

GRID_N = 128
MATRIC_D_CUT = 50.0
N_RESAMPLE = 500          # points on the common normalized time axis
N_BOOT = 5000             # bootstrap resamples for the within/between gap CI

M_SEEDS = 40              # ensemble size per soil (>= 30 required)
SEED_BASE = {"sand": 10000, "vertisol": 20000}   # disjoint seed ranges per soil


# ---------------------------------------------------------- numpy-only stats

def rankdata(a):
    """Average ranks (1-indexed), ties handled by averaging -- the standard
    rank transform Mann-Whitney needs. numpy-only (no scipy.stats.rankdata)."""
    a = np.asarray(a, dtype=float)
    sorter = np.argsort(a, kind="mergesort")
    inv = np.empty_like(sorter)
    inv[sorter] = np.arange(len(a))
    a_sorted = a[sorter]
    ranks_sorted = np.arange(1, len(a) + 1, dtype=float)
    unique_vals, first_idx, counts = np.unique(a_sorted, return_index=True, return_counts=True)
    for start, cnt in zip(first_idx, counts):
        if cnt > 1:
            ranks_sorted[start:start + cnt] = ranks_sorted[start:start + cnt].mean()
    return ranks_sorted[inv]


def _norm_cdf(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def mann_whitney_u(x, y):
    """Two-sided Mann-Whitney U test via the normal approximation with tie
    correction and continuity correction -- valid for the ensemble sizes
    used here (M>=30 per group). Returns (U_x, z, p_two_sided), where U_x is
    the U statistic for `x` (Sand) against `y` (Vertisol)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n1, n2 = len(x), len(y)
    combined = np.concatenate([x, y])
    ranks = rankdata(combined)
    r1 = ranks[:n1].sum()
    U1 = r1 - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    _, counts = np.unique(combined, return_counts=True)
    tie_term = float(np.sum(counts ** 3 - counts))
    N = n1 + n2
    sigma2 = (n1 * n2 / 12.0) * ((N + 1) - tie_term / (N * (N - 1)))
    sigma = math.sqrt(max(sigma2, 0.0))
    if sigma == 0:
        z = 0.0
    else:
        diff = U1 - mu
        cc = 0.5 if diff > 0 else (-0.5 if diff < 0 else 0.0)
        z = (diff - cc) / sigma
    p = 2.0 * (1.0 - _norm_cdf(abs(z)))
    p = min(max(p, 0.0), 1.0)
    return float(U1), float(z), float(p)


def cohens_d(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n1, n2 = len(x), len(y)
    s1, s2 = x.std(ddof=1), y.std(ddof=1)
    pooled = math.sqrt(((n1 - 1) * s1 ** 2 + (n2 - 1) * s2 ** 2) / (n1 + n2 - 2))
    if pooled == 0:
        return 0.0
    return float((x.mean() - y.mean()) / pooled)


def rank_biserial(U1, n1, n2):
    """Rank-biserial correlation from the Mann-Whitney U for group 1:
    r = 2*U1/(n1*n2) - 1, in [-1, 1]. |r| near 1 means near-total
    separation between the two groups' distributions."""
    return float(2.0 * U1 / (n1 * n2) - 1.0)


def bootstrap_gap_ci(within_dists, between_dists, n_boot=N_BOOT, seed=0):
    """Bootstrap CI on mean(between) - mean(within), resampling each
    distribution with replacement independently -- numpy-only (no scipy)."""
    rng = np.random.default_rng(seed)
    within_dists = np.asarray(within_dists, dtype=float)
    between_dists = np.asarray(between_dists, dtype=float)
    n_w, n_b = len(within_dists), len(between_dists)
    gaps = np.empty(n_boot)
    for i in range(n_boot):
        w_sample = within_dists[rng.integers(0, n_w, n_w)]
        b_sample = between_dists[rng.integers(0, n_b, n_b)]
        gaps[i] = b_sample.mean() - w_sample.mean()
    lo, hi = np.percentile(gaps, [2.5, 97.5])
    return float(gaps.mean()), float(lo), float(hi)


def normalized_time_resample(R_t, n_points=N_RESAMPLE):
    x_orig = np.linspace(0.0, 1.0, len(R_t))
    x_new = np.linspace(0.0, 1.0, n_points)
    return np.interp(x_new, x_orig, R_t)


# ------------------------------------------------------------- run helpers

def _T_for(soil, mapping_cfg):
    return int(mapping_cfg[f"stage10_T_{soil}"])


def _structure_cfg(soil, mapping_cfg, seed):
    cfg = copy.deepcopy(load_yaml(f"configs/soil_{soil}.yaml"))
    cfg["saturation"] = {"mode": "matric_regions", "matric_d_cut_um": MATRIC_D_CUT}
    cfg["grid"]["n"] = GRID_N
    cfg["T"] = _T_for(soil, mapping_cfg)
    cfg["pore"]["seed"] = int(seed)
    return cfg


def run_one(soil, seed, biology_cfg, mapping_cfg):
    cfg = _structure_cfg(soil, mapping_cfg, seed)
    out = dual_porosity.simulate_dual_porosity(biology_cfg, cfg, mapping_cfg, wicking_enabled=False)
    grids = out["grids"]
    R = out["R_t"]
    n3 = len(R) // 3
    resampled = normalized_time_resample(R)
    metrics_row = {
        "soil": soil, "seed": int(seed),
        "R_peak": float(R.max()),
        "peak_frac": float(np.argmax(R) / (len(R) - 1)),
        "R_end_over_peak": float(R[-1] / R.max()),
        "cum_co2_final": float(out["cum_co2"][-1]),
        "late_over_early": float(R[-n3:].mean() / R[:n3].mean()) if R[:n3].mean() > 0 else float("inf"),
        "n_regions": int(grids["region_diam_um"].size),
        "coating_cells": int(grids["coating_mask"].sum()),
        "dry_habitat_fraction_of_habitat": float(grids["dry_macro_fraction_of_habitat"]),
        "matrix_fraction": float(grids["matrix_fraction"]),
        "om_total": float(grids["om_total_init"]),
    }
    return metrics_row, resampled


def run_ensemble(soil, seeds, biology_cfg, mapping_cfg):
    rows = []
    curves = np.zeros((len(seeds), N_RESAMPLE))
    for i, seed in enumerate(seeds):
        row, resampled = run_one(soil, seed, biology_cfg, mapping_cfg)
        rows.append(row)
        curves[i] = resampled
    return rows, curves


# --------------------------------------------------------------- 1. figures

def plot_mean_bands(curves_by_soil, results_dir):
    x = np.linspace(0.0, 1.0, N_RESAMPLE)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for soil in SOILS:
        curves = curves_by_soil[soil]
        mean = curves.mean(axis=0)
        se = curves.std(axis=0, ddof=1) / math.sqrt(curves.shape[0])
        lo_p, hi_p = np.percentile(curves, [2.5, 97.5], axis=0)
        ax.plot(x, mean, color=COLORS[soil], lw=2.2, label=f"{LABEL[soil]} (mean, M={curves.shape[0]})")
        ax.fill_between(x, mean - 1.96 * se, mean + 1.96 * se, color=COLORS[soil], alpha=0.35,
                         label=f"{LABEL[soil]} mean +/- 1.96*SE")
        ax.plot(x, lo_p, color=COLORS[soil], lw=0.7, ls=":")
        ax.plot(x, hi_p, color=COLORS[soil], lw=0.7, ls=":", label=f"{LABEL[soil]} 2.5/97.5 percentile")
    ax.set_xlabel("fraction of own observation window")
    ax.set_ylabel("R(t)")
    ax.set_title(f"Stage 12: ensemble mean R(t) +/- 95% band (M={M_SEEDS} seeds/soil)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage12_mean_bands.png"), dpi=150)
    plt.close(fig)


def plot_metric_distributions(rows_by_soil, results_dir):
    metric_names = ["peak_frac", "R_end_over_peak", "late_over_early"]
    metric_labels = ["peak_frac (time-to-peak / window)", "R_end / R_peak", "late-third / early-third mean"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, mname, mlabel in zip(axes, metric_names, metric_labels):
        data = [
            np.array([r[mname] for r in rows_by_soil[s] if np.isfinite(r[mname])])
            for s in SOILS
        ]
        bp = ax.boxplot(data, tick_labels=[LABEL[s] for s in SOILS], patch_artist=True, showmeans=True)
        for patch, soil in zip(bp["boxes"], SOILS):
            patch.set_facecolor(COLORS[soil])
            patch.set_alpha(0.4)
        for i, d in enumerate(data, start=1):
            jitter = (np.random.default_rng(0).random(len(d)) - 0.5) * 0.15
            ax.scatter(np.full(len(d), i) + jitter, d, s=10, color="black", alpha=0.5, zorder=3)
        ax.set_title(mlabel, fontsize=10)
    fig.suptitle(f"Stage 12: per-seed shape metric distributions (M={M_SEEDS} seeds/soil)")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage12_metric_distributions.png"), dpi=150)
    plt.close(fig)


def plot_within_between(within_by_soil, between, results_dir):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    data = [within_by_soil["sand"], within_by_soil["vertisol"], between]
    labels = ["within Sand", "within Vertisol", "between Sand-Vertisol"]
    colors = [COLORS["sand"], COLORS["vertisol"], "#555555"]
    bp = ax.boxplot(data, tick_labels=labels, patch_artist=True, showmeans=True)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.4)
    ax.set_ylabel("pairwise normalized-curve distance")
    ax.set_title("Stage 12: within-soil vs between-soil curve distances")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage12_within_vs_between.png"), dpi=150)
    plt.close(fig)


# ----------------------------------------- optional causal feature sweep

# SS5 "optional causal check": vary ONE structural feature (Sand's coating
# fraction -- the lever that decides how much of Sand's OM sits on the
# fast, near-habitat coating pool vs. the slow trapped interior) over a
# few values, a few seeds each, and show the shape metrics move
# monotonically -- direct evidence the FEATURE itself is the driver, not
# just "this soil always bursts." Bounded: 4 values x 5 seeds = 20 runs.
SWEEP_COATING_VALUES = [0.5, 0.7, 0.85, 0.97]
SWEEP_SEEDS_PER_VALUE = 5
SWEEP_SEED_BASE = 30000


def run_coating_sweep(biology_cfg, mapping_cfg):
    rows = []
    for coat in SWEEP_COATING_VALUES:
        m = dict(mapping_cfg)
        m["stage11_coating_fraction_sand"] = coat  # Sand-only override, highest-priority tier
        for i in range(SWEEP_SEEDS_PER_VALUE):
            seed = SWEEP_SEED_BASE + int(coat * 1000) + i
            row, _ = run_one("sand", seed, biology_cfg, m)
            row["coating_fraction_swept"] = coat
            rows.append(row)
    return rows


def write_sweep_csv(rows, results_dir):
    path = os.path.join(results_dir, "stage12_coating_sweep.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


def plot_sweep(rows, results_dir):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, metric, mlabel in zip(
        axes, ["peak_frac", "R_end_over_peak"],
        ["peak_frac (time-to-peak / window)", "R_end / R_peak (crash depth)"],
    ):
        means, sds, xs = [], [], []
        for coat in SWEEP_COATING_VALUES:
            vals = np.array([r[metric] for r in rows if r["coating_fraction_swept"] == coat])
            xs.append(coat); means.append(vals.mean()); sds.append(vals.std(ddof=1))
        ax.errorbar(xs, means, yerr=sds, marker="o", color=COLORS["sand"], capsize=4, lw=1.8)
        ax.set_xlabel("Sand coating fraction (swept)")
        ax.set_ylabel(mlabel)
        ax.set_title(mlabel, fontsize=10)
    fig.suptitle(f"Stage 12 (optional): Sand's shape metrics vs. its own coating fraction "
                 f"({SWEEP_SEEDS_PER_VALUE} seeds/value)")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage12_coating_sweep.png"), dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------- 2. CSVs

def write_seed_metrics_csv(rows_by_soil, results_dir):
    all_rows = rows_by_soil["sand"] + rows_by_soil["vertisol"]
    path = os.path.join(results_dir, "stage12_seed_metrics.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(all_rows[0]))
        w.writeheader()
        w.writerows(all_rows)
    return all_rows


def write_metric_summary_csv(rows_by_soil, results_dir):
    metric_names = ["R_peak", "peak_frac", "R_end_over_peak", "cum_co2_final", "late_over_early"]
    rows = []
    for soil in SOILS:
        for m in metric_names:
            vals = np.array([r[m] for r in rows_by_soil[soil] if np.isfinite(r[m])])
            rows.append({"soil": LABEL[soil], "metric": m, "mean": float(vals.mean()),
                         "sd": float(vals.std(ddof=1)), "n": len(vals)})
    path = os.path.join(results_dir, "stage12_metric_summary.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


def write_significance_csv(rows_by_soil, results_dir):
    metric_names = ["peak_frac", "R_end_over_peak", "late_over_early", "R_peak", "cum_co2_final"]
    out_rows = []
    for m in metric_names:
        x = np.array([r[m] for r in rows_by_soil["sand"] if np.isfinite(r[m])])
        y = np.array([r[m] for r in rows_by_soil["vertisol"] if np.isfinite(r[m])])
        U1, z, p = mann_whitney_u(x, y)
        d = cohens_d(x, y)
        rb = rank_biserial(U1, len(x), len(y))
        out_rows.append({
            "metric": m, "n_sand": len(x), "n_vertisol": len(y),
            "mean_sand": float(x.mean()), "mean_vertisol": float(y.mean()),
            "U_sand": U1, "z": z, "p_two_sided": p,
            "cohens_d": d, "rank_biserial_r": rb,
        })
    path = os.path.join(results_dir, "stage12_significance.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0]))
        w.writeheader()
        w.writerows(out_rows)
    return out_rows


def write_within_between_csv(within_by_soil, between, gap_mean, gap_lo, gap_hi, results_dir):
    row = {
        "within_sand_mean": float(np.mean(within_by_soil["sand"])),
        "within_sand_n": len(within_by_soil["sand"]),
        "within_vertisol_mean": float(np.mean(within_by_soil["vertisol"])),
        "within_vertisol_n": len(within_by_soil["vertisol"]),
        "between_mean": float(np.mean(between)),
        "between_n": len(between),
        "gap_mean_boot": gap_mean, "gap_ci_lo_2.5pct": gap_lo, "gap_ci_hi_97.5pct": gap_hi,
        "gap_excludes_zero": bool(gap_lo > 0),
    }
    path = os.path.join(results_dir, "stage12_within_vs_between.csv")
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
    print("STAGE 12 -- statistical significance via seed ensembles")
    print("=" * 78)
    print(f"\nFrozen Stage-11 config (never re-tuned): "
          f"Sand om_total={mapping_cfg['stage10_om_total_sand']}, "
          f"coat={mapping_cfg['stage10_coating_fraction_sand']}, "
          f"D_scale={mapping_cfg['stage10_D_scale_sand']}, T={mapping_cfg['stage10_T_sand']} | "
          f"Vertisol om_total={mapping_cfg['stage11_om_total_vertisol']}, "
          f"coat={mapping_cfg['stage11_coating_fraction_vertisol']}, "
          f"D_scale={mapping_cfg['stage11_D_scale_vertisol']}, T={mapping_cfg['stage10_T_vertisol']}")

    t_start = time.time()

    rows_by_soil = {}
    curves_by_soil = {}
    for soil in SOILS:
        seeds = [SEED_BASE[soil] + i for i in range(M_SEEDS)]
        print(f"\n[1] running {M_SEEDS}-seed ensemble for {LABEL[soil]} (n={GRID_N}, T={_T_for(soil, mapping_cfg)}) ...")
        t0 = time.time()
        rows, curves = run_ensemble(soil, seeds, biology_cfg, mapping_cfg)
        elapsed = time.time() - t0
        print(f"    done in {elapsed:.1f}s ({elapsed / M_SEEDS:.2f}s/seed)")
        rows_by_soil[soil] = rows
        curves_by_soil[soil] = curves

    total_run_elapsed = time.time() - t_start
    print(f"\n[2] total ensemble runtime: {total_run_elapsed:.1f}s ({total_run_elapsed/60:.1f} min)")

    write_seed_metrics_csv(rows_by_soil, RESULTS_DIR)
    metric_summary = write_metric_summary_csv(rows_by_soil, RESULTS_DIR)
    print("\n[3] per-soil metric mean +/- SD:")
    for r in metric_summary:
        print(f"    {r['soil']:10} {r['metric']:18} mean={r['mean']:10.4f}  sd={r['sd']:10.4f}  n={r['n']}")

    plot_mean_bands(curves_by_soil, RESULTS_DIR)
    plot_metric_distributions(rows_by_soil, RESULTS_DIR)

    print("\n[4] significance (Mann-Whitney U, Sand vs Vertisol):")
    sig_rows = write_significance_csv(rows_by_soil, RESULTS_DIR)
    for r in sig_rows:
        print(f"    {r['metric']:18} U={r['U_sand']:8.1f}  z={r['z']:7.2f}  p={r['p_two_sided']:.2e}  "
              f"Cohen's d={r['cohens_d']:7.2f}  rank-biserial r={r['rank_biserial_r']:7.3f}")

    print("\n[5] within-soil vs between-soil curve distances ...")
    within_by_soil = {}
    for soil in SOILS:
        curves = curves_by_soil[soil]
        M = curves.shape[0]
        dists = []
        for i in range(M):
            for j in range(i + 1, M):
                dists.append(metrics.pairwise_distance(curves[i], curves[j]))
        within_by_soil[soil] = np.array(dists)

    between = []
    for i in range(M_SEEDS):
        for j in range(M_SEEDS):
            between.append(metrics.pairwise_distance(curves_by_soil["sand"][i], curves_by_soil["vertisol"][j]))
    between = np.array(between)

    gap_mean, gap_lo, gap_hi = bootstrap_gap_ci(
        np.concatenate([within_by_soil["sand"], within_by_soil["vertisol"]]), between)
    print(f"    within Sand:     mean={within_by_soil['sand'].mean():.4f}  n={len(within_by_soil['sand'])}")
    print(f"    within Vertisol: mean={within_by_soil['vertisol'].mean():.4f}  n={len(within_by_soil['vertisol'])}")
    print(f"    between:         mean={between.mean():.4f}  n={len(between)}")
    print(f"    bootstrap gap (between - within): mean={gap_mean:.4f}  95% CI=[{gap_lo:.4f}, {gap_hi:.4f}]  "
          f"excludes 0: {gap_lo > 0}")

    plot_within_between(within_by_soil, between, RESULTS_DIR)
    write_within_between_csv(within_by_soil, between, gap_mean, gap_lo, gap_hi, RESULTS_DIR)

    with open(LOGBOOK_PATH, "a", encoding="utf-8") as f:
        f.write(f"\n**Stage 12 headline: M={M_SEEDS} seeds/soil, total ensemble runtime "
                f"{total_run_elapsed:.1f}s ({total_run_elapsed/60:.1f} min). "
                f"Mann-Whitney (Sand vs Vertisol) peak_frac p={sig_rows[0]['p_two_sided']:.2e}, "
                f"R_end/R_peak p={sig_rows[1]['p_two_sided']:.2e}, late/early p={sig_rows[2]['p_two_sided']:.2e}. "
                f"Between-soil mean distance={between.mean():.3f} vs within-soil "
                f"({within_by_soil['sand'].mean():.3f}/{within_by_soil['vertisol'].mean():.3f}); "
                f"bootstrap gap 95% CI=[{gap_lo:.3f}, {gap_hi:.3f}] (excludes 0: {gap_lo>0}). "
                f"{datetime.datetime.now().isoformat(timespec='seconds')}. Full detail in "
                "docs/stage_results/STAGE12_RESULTS.md.**\n")

    print(f"\n[6] optional causal check: Sand coating-fraction sweep "
          f"({SWEEP_COATING_VALUES} x {SWEEP_SEEDS_PER_VALUE} seeds) ...")
    sweep_rows = run_coating_sweep(biology_cfg, mapping_cfg)
    write_sweep_csv(sweep_rows, RESULTS_DIR)
    plot_sweep(sweep_rows, RESULTS_DIR)
    for coat in SWEEP_COATING_VALUES:
        vals = [r for r in sweep_rows if r["coating_fraction_swept"] == coat]
        pf = np.array([r["peak_frac"] for r in vals])
        rp = np.array([r["R_end_over_peak"] for r in vals])
        print(f"    coating={coat:.2f}: peak_frac={pf.mean():.3f}+/-{pf.std(ddof=1):.3f}  "
              f"R_end/R_peak={rp.mean():.3f}+/-{rp.std(ddof=1):.3f}")

    print(f"\n[7] all results written to {RESULTS_DIR}/")
    print(f"\nTOTAL SCRIPT RUNTIME: {time.time() - t_start:.1f}s")


if __name__ == "__main__":
    main()
