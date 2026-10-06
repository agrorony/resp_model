"""
Stage 14 runner (prompts/soil_respiration_prompt_stage14_habitat_diagnostic.md).

A DIAGNOSTIC stage, not a retuning stage. Stage 11 reported, as an
explicitly unresolved limitation: `habitats_recruited(t)` reaches its final
count within the first ~2% of the observation window and stays completely
flat afterward, for Vertisol, even though `R(t)` itself keeps genuinely
rising. The suspected cause (never confirmed): the implicit backward-Euler
solver's per-step equilibration crosses the "recruited" activation
threshold (`growth > 1e-9`) almost everywhere almost instantly, independent
of transport speed. This stage confirms or refutes that directly, by
instrumenting `simulate_dual_porosity`'s new (purely observational)
`debug_growth_distribution` flag.

Uses Stage 11's frozen LITERATURE-PSD Vertisol config, not Stage 13's real
data: Stage 13 found real Vertisol data collapses to zero growth
everywhere (OM_total_placed=0.0, R_peak pinned at the trivial maintenance
floor) -- there is no growth distribution to diagnose there. Stage 11's
config is the most recent one where Vertisol actually grows.

Usage:
    python stage14_run.py
"""
from __future__ import annotations

import copy
import csv
import datetime
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from model import load_yaml
import dual_porosity

RESULTS_DIR = os.path.join("results", "stage14")
LOGBOOK_PATH = "LOGBOOK.md"

GRID_N = 128
MATRIC_D_CUT = 50.0
GROWTH_THRESHOLD = 1e-9
WINDOW_2PCT_FRAC = 0.02

# Stage 11's exact recorded reference numbers for Vertisol (n=128,
# results/stage11/stage11_headline.csv) -- SS3's required identity check
# compares the debug-instrumented run against these AND against a
# same-session non-debug run, so it catches both drift since Stage 11 and
# any accidental side effect of the new flag itself.
VERTISOL_REFERENCE = {
    "R_peak": 3.7561399861379976,
    "R_end": 3.7493280228711905,
    "B_peak_over_B0": 2.3372611216671735,
    "habitats_recruited_final": 55.0,
    "n_regions": 70,
}


def _structure_cfg(mapping_cfg):
    cfg = copy.deepcopy(load_yaml("configs/soil_vertisol.yaml"))
    cfg["saturation"] = {"mode": "matric_regions", "matric_d_cut_um": MATRIC_D_CUT}
    cfg["grid"]["n"] = GRID_N
    cfg["T"] = int(mapping_cfg["stage10_T_vertisol"])
    return cfg


def run_vertisol(biology_cfg, mapping_cfg, debug=False):
    cfg = _structure_cfg(mapping_cfg)
    return dual_porosity.simulate_dual_porosity(
        biology_cfg, cfg, mapping_cfg, wicking_enabled=False, debug_growth_distribution=debug)


# ------------------------------------------------------- 1/3. identity check

def identity_check(out_baseline, out_debug, results_dir):
    """SS3: confirm the debug instrumentation is purely observational --
    R(t), B(t), and all existing headline numbers must be identical
    before/after, not just noise-level (no stochastic element is touched by
    the flag, so exact equality is the correct, stronger check here)."""
    checks = {
        "R_t exactly equal": bool(np.array_equal(out_baseline["R_t"], out_debug["R_t"])),
        "B_total_t exactly equal": bool(np.array_equal(out_baseline["B_total_t"], out_debug["B_total_t"])),
        "habitats_recruited_t exactly equal": bool(np.array_equal(
            out_baseline["habitats_recruited_t"], out_debug["habitats_recruited_t"])),
        "fuel_remaining_t exactly equal": bool(np.array_equal(
            out_baseline["fuel_remaining_t"], out_debug["fuel_remaining_t"])),
        "R_peak matches Stage 11 reference": bool(
            abs(float(out_debug["R_t"].max()) - VERTISOL_REFERENCE["R_peak"]) / VERTISOL_REFERENCE["R_peak"] < 0.01),
        "B_peak/B0 matches Stage 11 reference": bool(
            abs(float(out_debug["B_total_t"].max() / out_debug["grids"]["B0"].sum())
                - VERTISOL_REFERENCE["B_peak_over_B0"]) / VERTISOL_REFERENCE["B_peak_over_B0"] < 0.01),
    }
    row = {**checks, "all_pass": all(checks.values())}
    path = os.path.join(results_dir, "stage14_identity_check.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(row))
        w.writeheader()
        w.writerow(row)
    return row


# ------------------------------------------------------- 2. mechanism confirmation

def region_recruitment_analysis(out_debug):
    """Per-region: first timestep any cell crosses GROWTH_THRESHOLD (boolean
    activation, matching count_active_habitats exactly), the region's own
    mean-growth trajectory, its eventual peak mean growth, and the first
    timestep that mean reaches half its own eventual peak (a graded,
    per-region 'time to half-max' -- SS2's candidate metric)."""
    gd = out_debug["growth_distribution_t"]          # (T, n_macro_cells)
    macro_mask = out_debug["macro_mask_debug"]
    region_labels = out_debug["region_labels_debug"]
    region_ids = region_labels[macro_mask]            # (n_macro_cells,), parallel to gd's columns
    n_regions = int(region_ids.max())
    T = gd.shape[0]

    first_active_t = np.full(n_regions, -1, dtype=int)
    first_half_max_t = np.full(n_regions, -1, dtype=int)
    eventual_peak_mean = np.zeros(n_regions)
    region_mean_series = np.zeros((T, n_regions))

    for k in range(1, n_regions + 1):
        cols = np.where(region_ids == k)[0]
        region_growth = gd[:, cols]                   # (T, cells_in_region_k)
        region_max_cell = region_growth.max(axis=1)
        region_mean = region_growth.mean(axis=1)
        region_mean_series[:, k - 1] = region_mean

        active = np.where(region_max_cell > GROWTH_THRESHOLD)[0]
        if active.size:
            first_active_t[k - 1] = int(active[0])

        peak = float(region_mean.max())
        eventual_peak_mean[k - 1] = peak
        if peak > 0:
            half = np.where(region_mean >= 0.5 * peak)[0]
            if half.size:
                first_half_max_t[k - 1] = int(half[0])

    return {
        "n_regions": n_regions, "T": T,
        "first_active_t": first_active_t,
        "first_half_max_t": first_half_max_t,
        "eventual_peak_mean": eventual_peak_mean,
        "region_mean_series": region_mean_series,
    }


def write_mechanism_csv(analysis, out_debug, results_dir):
    n_regions, T = analysis["n_regions"], analysis["T"]
    window_2pct_t = int(WINDOW_2PCT_FRAC * T)
    first_active = analysis["first_active_t"]
    activated_ever = first_active >= 0

    rows = []
    for k in range(n_regions):
        rows.append({
            "region": k + 1,
            "first_active_t": int(first_active[k]),
            "first_active_frac_of_window": float(first_active[k] / (T - 1)) if first_active[k] >= 0 else float("nan"),
            "first_half_max_t": int(analysis["first_half_max_t"][k]),
            "first_half_max_frac_of_window": (
                float(analysis["first_half_max_t"][k] / (T - 1)) if analysis["first_half_max_t"][k] >= 0 else float("nan")
            ),
            "eventual_peak_mean_growth": float(analysis["eventual_peak_mean"][k]),
        })
    path = os.path.join(results_dir, "stage14_region_recruitment.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    # cross-check against the existing boolean diagnostic: at every t, the
    # number of regions active by first_active_t<=t must equal
    # habitats_recruited_t[t] exactly (ties the new instrumentation to the
    # already-trusted metric, not a parallel, unverified computation).
    habitats_recruited_t = out_debug["habitats_recruited_t"]
    reconstructed = np.array([
        int(((first_active >= 0) & (first_active <= t)).sum()) for t in range(T)
    ])
    cross_check_exact = bool(np.array_equal(reconstructed, habitats_recruited_t.astype(int)))

    summary = {
        "n_regions": n_regions,
        "n_regions_ever_active": int(activated_ever.sum()),
        "n_regions_active_by_t0": int((first_active == 0).sum()),
        "n_regions_active_by_t2": int(((first_active >= 0) & (first_active <= 2)).sum()),
        "n_regions_active_by_2pct_window": int(((first_active >= 0) & (first_active <= window_2pct_t)).sum()),
        "window_2pct_t": window_2pct_t,
        "reconstructed_matches_habitats_recruited_t_exactly": cross_check_exact,
        "median_first_active_t": float(np.median(first_active[activated_ever])) if activated_ever.any() else float("nan"),
        "median_first_half_max_t": float(np.median(
            analysis["first_half_max_t"][analysis["first_half_max_t"] >= 0])) if (analysis["first_half_max_t"] >= 0).any() else float("nan"),
        "median_first_half_max_frac_of_window": float(np.median(
            analysis["first_half_max_t"][analysis["first_half_max_t"] >= 0] / (T - 1))) if (analysis["first_half_max_t"] >= 0).any() else float("nan"),
    }
    path2 = os.path.join(results_dir, "stage14_mechanism_summary.csv")
    with open(path2, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(summary))
        w.writeheader()
        w.writerow(summary)
    return rows, summary


# ------------------------------------------------------- 3. candidate graded metric

def graded_metric_comparison(out_debug, analysis, results_dir):
    """SS2: test a graded, growth-WEIGHTED 'active mass' metric
    (sum of raw growth over ALL macropore cells, every step) against the
    boolean recruited-count -- and against R(t)'s own known-genuine gradual
    rise (Stage 11: driven by per-region growth RATE building up, not new
    regions joining)."""
    gd = out_debug["growth_distribution_t"]
    T = gd.shape[0]
    active_mass_t = gd.sum(axis=1)                       # graded metric candidate 1
    higher_threshold = 0.1 * float(gd.max())              # candidate 2: 10% of the run's peak single-cell growth
    higher_thresh_count_t = (gd > higher_threshold).astype(int)
    region_labels_active = out_debug["region_labels_debug"][out_debug["macro_mask_debug"]]
    n_regions = analysis["n_regions"]
    higher_thresh_recruited_t = np.zeros(T, dtype=int)
    for t in range(T):
        active_cells = gd[t] > higher_threshold
        if active_cells.any():
            higher_thresh_recruited_t[t] = len(np.unique(region_labels_active[active_cells]))

    R_t = out_debug["R_t"]
    habitats_recruited_t = out_debug["habitats_recruited_t"]

    def norm(x):
        x = np.asarray(x, dtype=float)
        m = x.max()
        return x / m if m > 0 else x

    rows = []
    for t in range(T):
        rows.append({
            "t": t,
            "R_t": float(R_t[t]),
            "R_t_normalized": float(norm(R_t)[t]),
            "habitats_recruited_t": float(habitats_recruited_t[t]),
            "habitats_recruited_t_normalized": float(norm(habitats_recruited_t)[t]),
            "active_mass_t": float(active_mass_t[t]),
            "active_mass_t_normalized": float(norm(active_mass_t)[t]),
            "higher_threshold_recruited_t": int(higher_thresh_recruited_t[t]),
            "higher_threshold_recruited_t_normalized": float(norm(higher_thresh_recruited_t)[t]),
        })
    path = os.path.join(results_dir, "stage14_graded_metric_comparison.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    # correlation of each candidate (normalized) against normalized R(t) --
    # a simple, numpy-only agreement check (project invariant: no scipy).
    def pearson(a, b):
        a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
        a, b = a - a.mean(), b - b.mean()
        denom = np.sqrt((a ** 2).sum() * (b ** 2).sum())
        return float((a * b).sum() / denom) if denom > 0 else 0.0

    R_norm = norm(R_t)
    corr_summary = {
        "corr_habitats_recruited_vs_R": pearson(norm(habitats_recruited_t), R_norm),
        "corr_active_mass_vs_R": pearson(norm(active_mass_t), R_norm),
        "corr_higher_threshold_recruited_vs_R": pearson(norm(higher_thresh_recruited_t), R_norm),
        "higher_threshold_value": higher_threshold,
    }
    path2 = os.path.join(results_dir, "stage14_graded_metric_correlations.csv")
    with open(path2, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(corr_summary))
        w.writeheader()
        w.writerow(corr_summary)

    return active_mass_t, higher_thresh_recruited_t, corr_summary


# --------------------------------------------------------------- figures

def plot_growth_distribution_heatmap(out_debug, analysis, results_dir):
    gd = out_debug["growth_distribution_t"]
    region_ids = out_debug["region_labels_debug"][out_debug["macro_mask_debug"]]
    n_regions = analysis["n_regions"]
    T = gd.shape[0]
    region_mean_series = analysis["region_mean_series"]   # (T, n_regions)

    order = np.argsort(analysis["first_active_t"])
    mat = np.log10(np.clip(region_mean_series[:, order].T, 1e-15, None))

    fig, ax = plt.subplots(figsize=(10, 6.5))
    im = ax.imshow(mat, aspect="auto", origin="lower", cmap="viridis",
                    extent=[0, T, 0, n_regions])
    ax.set_xlabel("timestep"); ax.set_ylabel(f"region (sorted by first-active time, n={n_regions})")
    ax.set_title("Stage 14: per-region mean growth over time, log10 (Vertisol)")
    fig.colorbar(im, ax=ax, label="log10(mean growth in region)")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage14_growth_distribution_heatmap.png"), dpi=150)
    plt.close(fig)

    # zoomed early-window view: does activation happen near-simultaneously?
    zoom_t = max(int(0.05 * T), 10)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    im = ax.imshow(mat[:, :zoom_t], aspect="auto", origin="lower", cmap="viridis",
                    extent=[0, zoom_t, 0, n_regions])
    ax.set_xlabel("timestep (zoomed to first "
                  f"{zoom_t} of {T} steps, {100*zoom_t/T:.1f}% of window)")
    ax.set_ylabel(f"region (sorted by first-active time, n={n_regions})")
    ax.set_title("Stage 14: growth distribution, zoomed to confirm/refute near-instant activation")
    fig.colorbar(im, ax=ax, label="log10(mean growth in region)")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage14_growth_distribution_heatmap_zoomed.png"), dpi=150)
    plt.close(fig)


def plot_mechanism_summary(analysis, results_dir):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    first_active = analysis["first_active_t"]
    first_half = analysis["first_half_max_t"]
    T = analysis["T"]

    axes[0].hist(first_active[first_active >= 0] / (T - 1), bins=30, color="#4b0082")
    axes[0].set_xlabel("fraction of window at first activation (growth>1e-9)")
    axes[0].set_ylabel("number of regions")
    axes[0].set_title("Boolean activation: WHEN each region first crosses threshold")

    axes[1].hist(first_half[first_half >= 0] / (T - 1), bins=30, color="#d62728")
    axes[1].set_xlabel("fraction of window at time-to-half-max mean growth")
    axes[1].set_ylabel("number of regions")
    axes[1].set_title("Graded metric: WHEN each region reaches half its own eventual growth")

    fig.suptitle("Stage 14: boolean activation (near-instant) vs. graded half-max timing (spread out)")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage14_activation_vs_halfmax_histograms.png"), dpi=150)
    plt.close(fig)


def plot_graded_metric_vs_R(out_debug, active_mass_t, higher_thresh_recruited_t, results_dir):
    R_t = out_debug["R_t"]
    habitats_recruited_t = out_debug["habitats_recruited_t"]
    t = np.arange(len(R_t)) * out_debug["dt"]

    def norm(x):
        x = np.asarray(x, dtype=float)
        m = x.max()
        return x / m if m > 0 else x

    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    ax.plot(t, norm(R_t), color="black", lw=2, label="R(t), normalized")
    ax.plot(t, norm(habitats_recruited_t), color="#4b0082", lw=1.5, ls="--",
            label="habitats_recruited(t) / max (boolean count, saturates instantly)")
    ax.plot(t, norm(active_mass_t), color="#d62728", lw=1.5,
            label="growth-weighted active mass(t) / max (graded)")
    ax.plot(t, norm(higher_thresh_recruited_t), color="#2ca02c", lw=1.5, ls=":",
            label="recruited count at 10x-higher threshold / max")
    ax.set_xlabel("model time"); ax.set_ylabel("normalized")
    ax.set_title("Stage 14: candidate recruitment metrics vs. R(t)'s genuine gradual rise (Vertisol)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage14_graded_metric_vs_R.png"), dpi=150)
    plt.close(fig)


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
    print("STAGE 14 -- investigating the habitats_recruited(t) saturation artifact")
    print("=" * 78)

    print("\n[1] baseline run (no instrumentation) + debug-instrumented run, Vertisol, n=128 ...")
    out_baseline = run_vertisol(biology_cfg, mapping_cfg, debug=False)
    out_debug = run_vertisol(biology_cfg, mapping_cfg, debug=True)

    id_row = identity_check(out_baseline, out_debug, RESULTS_DIR)
    print(f"    identity check all_pass={id_row['all_pass']}: {id_row}")

    print("\n[2] per-region recruitment analysis (boolean activation vs. graded half-max) ...")
    analysis = region_recruitment_analysis(out_debug)
    _, mech_summary = write_mechanism_csv(analysis, out_debug, RESULTS_DIR)
    print(f"    n_regions={mech_summary['n_regions']}, ever active={mech_summary['n_regions_ever_active']}")
    print(f"    active by t=0: {mech_summary['n_regions_active_by_t0']}, "
          f"by t<=2: {mech_summary['n_regions_active_by_t2']}, "
          f"by 2% of window (t<={mech_summary['window_2pct_t']}): {mech_summary['n_regions_active_by_2pct_window']}")
    print(f"    reconstructed-vs-habitats_recruited_t exact match: "
          f"{mech_summary['reconstructed_matches_habitats_recruited_t_exactly']}")
    print(f"    median first-active t: {mech_summary['median_first_active_t']:.1f}  "
          f"median first-half-max t: {mech_summary['median_first_half_max_t']:.1f} "
          f"({mech_summary['median_first_half_max_frac_of_window']:.3f} of window)")

    print("\n[3] candidate graded metrics vs. R(t) ...")
    active_mass_t, higher_thresh_recruited_t, corr = graded_metric_comparison(out_debug, analysis, RESULTS_DIR)
    print(f"    correlation with normalized R(t): habitats_recruited={corr['corr_habitats_recruited_vs_R']:.4f}  "
          f"active_mass={corr['corr_active_mass_vs_R']:.4f}  "
          f"higher_threshold_recruited={corr['corr_higher_threshold_recruited_vs_R']:.4f}")

    plot_growth_distribution_heatmap(out_debug, analysis, RESULTS_DIR)
    plot_mechanism_summary(analysis, RESULTS_DIR)
    plot_graded_metric_vs_R(out_debug, active_mass_t, higher_thresh_recruited_t, RESULTS_DIR)

    confirmed = mech_summary["n_regions_active_by_2pct_window"] >= 0.9 * mech_summary["n_regions_ever_active"]

    with open(LOGBOOK_PATH, "a", encoding="utf-8") as f:
        f.write(f"\n**Stage 14 headline: equilibration hypothesis "
                f"{'CONFIRMED' if confirmed else 'NOT confirmed as stated'} -- "
                f"{mech_summary['n_regions_active_by_2pct_window']}/{mech_summary['n_regions_ever_active']} "
                f"active regions crossed growth>1e-9 within the first 2% of the window "
                f"(median first-active t={mech_summary['median_first_active_t']:.1f}), while the graded "
                f"time-to-half-max metric spreads across the window (median t={mech_summary['median_first_half_max_t']:.1f}, "
                f"{mech_summary['median_first_half_max_frac_of_window']:.1%} of window) and correlates with R(t) at "
                f"r={corr['corr_active_mass_vs_R']:.3f} (growth-weighted active mass) vs. "
                f"r={corr['corr_habitats_recruited_vs_R']:.3f} (boolean count). Identity check all_pass="
                f"{id_row['all_pass']} (debug instrumentation confirmed purely observational). "
                f"{datetime.datetime.now().isoformat(timespec='seconds')}. Full detail in "
                "docs/stage_results/STAGE14_RESULTS.md.**\n")

    print(f"\n[4] OVERALL: equilibration hypothesis {'CONFIRMED' if confirmed else 'NOT confirmed as stated'}, "
          f"identity check all_pass={id_row['all_pass']}")
    print(f"\n[5] all results written to {RESULTS_DIR}/")


if __name__ == "__main__":
    main()
