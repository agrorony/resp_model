"""
Stage-3.5 runner (soil_respiration_prompt_stage3_5.md).

Stage 3 fed the model the three soils' real, image-derived PSDs and could only
tell the fine soil apart from the two coarse ones. Stage 3.5 keeps the whole
pipeline and swaps the PSD source for literature-informed lognormal mixtures
for the three real soil types (Loess / Sand / Vertisol, configs/psd_literature
.yaml), fixes the D(d) fusion (SS3), and -- because the experiment ran at
CONSTANT moisture in closed humid jars, not as a drying run -- evaluates the
temporal shape of R(t) at a single fixed theta rather than sweeping theta.

The validation target is the measured respiration pattern, qualitatively only
(ranking and shape, not absolute mg CO2):

    Vertisol   rising     (~104 -> 142 -> 158 -> 173)
    Sand       falling    (~116 ->  94 ->  90 ->  74)
    Loess      low/erratic(~ 68 ->  25 ->  50 ->  11)     mg CO2 kg-1 day-1

What this script produces (all into results/stage35/):
  1. the three loaded literature PSDs + a habitable-fraction sanity check
  2. the D(d) fusion fix (SS3) and its separation / stability checks
  3. HEADLINE: the constant-theta time-evolution run at n=128 (SS4)
  4. the experiment overlay -- model vs measured temporal shape (SS5)
  5. the theta-sensitivity check at 2-3 nearby fixed theta (SS4, secondary)
  6. DIAGNOSTIC: the same three soils at a shared MATRIC POTENTIAL instead of
     a shared theta -- the headline run fails for two of three soils and this
     is what explains why (see STAGE3_5_RESULTS.md)
  7. isolation-ratio / distinctness diagnostics at the fixed theta (SS6)

Usage:
    python stage35_run.py
"""
from __future__ import annotations

import copy
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from model import load_yaml, simulate, build_grids
from pore_field import D_of_d
import metrics
import psd_data
import psd_parametric

SOILS = ["loess", "sand", "vertisol"]
LABEL = {"loess": "Loess", "sand": "Sand", "vertisol": "Vertisol"}
COLORS = {"loess": "#b8860b", "sand": "#d62728", "vertisol": "#4b0082"}

RESULTS_DIR = os.path.join("results", "stage35")

# SS4: a single constant moisture, the SAME fixed value for all soils, held for
# the whole simulation -- moist but aerated. This is the headline condition.
THETA_CONST = 0.6
THETA_SENSITIVITY = [0.5, 0.6, 0.7]     # secondary robustness check only
# a wider, purely diagnostic sweep used to show WHERE each soil's habitat wets
THETA_DIAGNOSTIC = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

GRID_N = 128
GRID_N_SMALL = 40        # optional second grid size (SS6)

# SS6: this is a defined comparison, not a search -- no iteration budget is
# spent nudging anything. Kept for the invariant's sake.
MAX_ITERATIONS = 12

# The measured experiment (SS0). Period means, mg CO2 kg-1 day-1.
EXPERIMENT = {
    "vertisol": [104.0, 142.0, 158.0, 173.0],
    "sand":     [116.0,  94.0,  90.0,  74.0],
    "loess":    [ 68.0,  25.0,  50.0,  11.0],
}
PERIOD_DAYS = [3, 4, 3, 5]              # P1..P4, 15 days total
PERIOD_NAMES = ["P1", "P2", "P3", "P4"]

# K(d) habitat window (mapping.yaml K_d_low/K_d_high), for the PSD plots
HAB_LO, HAB_HI = 30.0, 150.0


# --------------------------------------------------------------- runs

def run_one(soil, sat_cfg, grid_n, biology_cfg, mapping_cfg):
    """One soil, one saturation setting, one grid size."""
    structure_cfg = copy.deepcopy(load_yaml(f"configs/soil_{soil}.yaml"))
    structure_cfg["saturation"] = sat_cfg
    structure_cfg["grid"]["n"] = grid_n
    out = simulate(biology_cfg, structure_cfg, mapping_cfg)
    conn = metrics.connectivity_diagnostics(out["grids"])
    return out, conn


def constant_theta(theta):
    return {"mode": "global_theta", "theta": theta}


# ------------------------------------------------- 1. the literature PSDs

def psd_summary(mapping_cfg, results_dir):
    """SS8.1: plot the three literature PSDs and sanity-check the habitable
    (30-150 um) volume fraction of each."""
    psds = {s: psd_parametric.build_soil_psd(s) for s in SOILS}

    rows = []
    for s in SOILS:
        p = psds[s]
        rows.append({
            "soil": LABEL[s],
            "total_porosity": p["total_porosity"],
            "median_diameter_um": psd_data.median_diameter_um(p),
            "vol_frac_below_30um": psd_parametric.volume_fraction_between(p, 1e-4, HAB_LO),
            "habitable_vol_frac_30_150um": psd_parametric.volume_fraction_between(p, HAB_LO, HAB_HI),
            "vol_frac_above_150um": psd_parametric.volume_fraction_between(p, HAB_HI, 5000.0),
        })

    path = os.path.join(results_dir, "stage35_psd_summary.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for s in SOILS:
        p = psds[s]
        # plot the density against LOG diameter (dV/dlog10 d), which is how a
        # PSD's modes actually read as bumps -- dV/dd on a log axis would hide
        # the coarse modes entirely.
        dlog = np.diff(np.log10(p["bin_edges_um"]))
        axes[0].plot(p["bin_centers_um"], p["volume_count"] / dlog, "-",
                     color=COLORS[s], label=LABEL[s], lw=1.8)
        axes[1].plot(p["bin_edges_um"], p["cdf_edges"], "-",
                     color=COLORS[s], label=LABEL[s], lw=1.8)

    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel("pore diameter d (um, log scale)")
        ax.axvspan(HAB_LO, HAB_HI, color="gray", alpha=0.15)
    axes[0].set_ylabel("volume density  dV / dlog10(d)")
    axes[0].set_title("Literature-informed PSDs (grey = 30-150 um habitat window)")
    axes[0].legend(fontsize=9)
    axes[1].set_ylabel("cumulative volume fraction")
    axes[1].set_title("Empirical CDF (used for the quantile mapping)")
    axes[1].legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage35_psd_loaded.png"), dpi=150)
    plt.close(fig)

    return psds, rows


# --------------------------------------------- 2. the D(d) fusion fix (SS3)

def d_of_d_check(biology_cfg, mapping_cfg, results_dir):
    """SS3: Stage 3 fused the two coarse soils because D(d) saturated at
    D_d_ref=90 um -- every large pore looked alike to transport. Show that at
    the new D_d_ref, D keeps rising appreciably across ~30-300 um and that two
    coarse-but-different PSDs no longer collapse to the same mean D. Also
    re-check the numerical-stability bound."""
    d_ref_old, d_ref_new = 90.0, float(mapping_cfg["D_d_ref"])
    grid = np.logspace(0, 3.2, 300)

    def D_at(d, ref):
        m = dict(mapping_cfg)
        m["D_d_ref"] = ref
        return D_of_d(d, m)

    def mean_D_over_psd(psd, ref):
        d = np.interp(np.linspace(0, 1, 5000), psd["cdf_edges"], psd["bin_edges_um"])
        return float(D_at(d, ref).mean())

    # the two coarse-but-different soils Stage 3 fused, i.e. the real measured
    # PSDs of the two coarse samples
    coarse = {"nlm (coarse)": psd_data.load_soil_psd("data/psd/nlm"),
              "mishmar (coarse)": psd_data.load_soil_psd("data/psd/mishmar")}

    rows = []
    for ref in (d_ref_old, d_ref_new):
        spread = float(D_at(np.array([300.0]), ref)[0] / D_at(np.array([30.0]), ref)[0])
        a = mean_D_over_psd(coarse["nlm (coarse)"], ref)
        b = mean_D_over_psd(coarse["mishmar (coarse)"], ref)
        rows.append({
            "D_d_ref_um": ref,
            "D_at_30um": float(D_at(np.array([30.0]), ref)[0]),
            "D_at_150um": float(D_at(np.array([150.0]), ref)[0]),
            "D_at_300um": float(D_at(np.array([300.0]), ref)[0]),
            "D300_over_D30": spread,
            "mean_D_coarse_nlm": a,
            "mean_D_coarse_mishmar": b,
            "coarse_pair_separation": max(a, b) / min(a, b),
        })

    path = os.path.join(results_dir, "stage35_D_of_d_fix.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    # numerical stability, on the actual fields the run will use
    dt, worst = biology_cfg["dt"], 0.0
    for s in SOILS:
        cfg = load_yaml(f"configs/soil_{s}.yaml")
        g = build_grids(cfg, mapping_cfg)
        worst = max(worst, float(g["D_full"].max()) * dt / cfg["grid"]["dx"] ** 2)
    assert worst < 0.2, f"diffusion stability violated after the D(d) fix: {worst:.3f}"

    fig, ax = plt.subplots(figsize=(7.5, 5))
    for ref, style in ((d_ref_old, "--"), (d_ref_new, "-")):
        ax.plot(grid, D_at(grid, ref), style, lw=1.8,
                label=f"D_d_ref = {ref:.0f} um" + (" (Stage 3)" if ref == d_ref_old else " (Stage 3.5)"))
    ax.axvspan(HAB_LO, HAB_HI, color="gray", alpha=0.15, label="habitat window")
    ax.axvspan(30, 300, color="tab:blue", alpha=0.05)
    ax.set_xscale("log")
    ax.set_xlabel("pore diameter d (um, log scale)")
    ax.set_ylabel("D(d)")
    ax.set_title("D(d) fusion fix: D must keep discriminating across ~30-300 um\n"
                 f"D(300)/D(30): {rows[0]['D300_over_D30']:.1f}x -> {rows[1]['D300_over_D30']:.1f}x   |   "
                 f"coarse-pair separation: {rows[0]['coarse_pair_separation']:.2f}x -> "
                 f"{rows[1]['coarse_pair_separation']:.2f}x", fontsize=10)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage35_D_of_d_fix.png"), dpi=150)
    plt.close(fig)

    return rows, worst


# ------------------------------------- 3-4. the constant-theta headline run

def period_means(R_t):
    """Average the model's R(t) over four windows whose lengths are in the
    same 3/4/3/5-day proportion as the experiment's P1-P4, so the model and
    the experiment can be compared period-by-period."""
    R = np.asarray(R_t, dtype=float)
    total = sum(PERIOD_DAYS)
    edges = np.cumsum([0] + PERIOD_DAYS) / total * len(R)
    idx = np.round(edges).astype(int)
    return np.array([R[idx[k]:idx[k + 1]].mean() for k in range(4)])


def shape_of_periods(vals):
    """rising / falling / erratic, from the four period means: sign of the
    end-to-start change, unless the series reverses direction (non-monotone by
    more than a small tolerance), in which case it is erratic."""
    v = np.asarray(vals, dtype=float)
    scale = max(abs(v).max(), 1e-30)
    diffs = np.diff(v)
    ups = (diffs > 0.05 * scale).sum()
    downs = (diffs < -0.05 * scale).sum()
    if ups and downs:
        return "erratic"
    if v[-1] > v[0] * 1.05:
        return "rising"
    if v[-1] < v[0] * 0.95:
        return "falling"
    return "flat"


def normalized_periods(vals):
    """Each soil's four period means divided by its own P1 -- strips the
    absolute level so only the temporal SHAPE is compared (SS5: qualitative
    overlay, the point is shape/ranking, not absolute mg CO2)."""
    v = np.asarray(vals, dtype=float)
    return v / v[0] if v[0] > 0 else v


def headline_run(biology_cfg, mapping_cfg, results_dir):
    """SS4: each soil at ONE constant theta, held for the whole simulation."""
    data = {}
    for s in SOILS:
        data[s] = run_one(s, constant_theta(THETA_CONST), GRID_N, biology_cfg, mapping_cfg)
    return data


def write_headline_csv(data, results_dir, tag="stage35_constant_theta"):
    path = os.path.join(results_dir, f"{tag}.csv")
    rows = []
    for s in SOILS:
        out, conn = data[s]
        pm = period_means(out["R_t"])
        rows.append({
            "soil": LABEL[s], "theta": conn["theta"],
            "R_peak": float(out["R_t"].max()), "R_end": float(out["R_t"][-1]),
            "cum_co2_final": float(out["cum_co2"][-1]),
            "B_final_over_B0": float(out["B_final"].sum() / out["grids"]["B0"].sum()),
            "outcome": metrics.classify_outcome(out),
            "curve_shape": metrics.describe_shape(out["R_t"])["shape"],
            "P1": pm[0], "P2": pm[1], "P3": pm[2], "P4": pm[3],
            "period_shape_model": shape_of_periods(pm),
            "period_shape_experiment": shape_of_periods(EXPERIMENT[s]),
            "habitat_fraction": conn["habitat_fraction"],
            "wet_habitat_fraction": conn["wet_habitat_fraction"],
            "habitat_wet_share": conn["habitat_wet_share"],
            "isolation_ratio": conn["isolation_ratio"],
            "n_clusters": conn["n_clusters"],
            "largest_cluster_frac": conn["largest_cluster_frac"],
            "om_biomass_connectivity": conn["om_biomass_connectivity"],
            "accessible_om_fraction": conn["accessible_om_fraction"],
        })
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


def plot_curves(data, results_dir, title, fname):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for s in SOILS:
        out, _ = data[s]
        t = np.arange(len(out["R_t"])) * out["dt"]
        axes[0].plot(t, out["R_t"], color=COLORS[s], label=LABEL[s], lw=1.8)
        axes[1].plot(t, out["cum_co2"], color=COLORS[s], label=LABEL[s], lw=1.8)
    axes[0].set_xlabel("model time")
    axes[0].set_ylabel("R(t) - respiration rate")
    axes[0].set_title("Respiration rate")
    axes[0].legend()
    axes[1].set_xlabel("model time")
    axes[1].set_ylabel("cumulative CO2")
    axes[1].set_title("Cumulative CO2")
    axes[1].legend()
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, fname), dpi=150)
    plt.close(fig)


def plot_water_map(soil, out, results_dir, tag):
    grids = out["grids"]
    fig, axes = plt.subplots(1, 5, figsize=(20, 4))
    fields = [
        ("pore diameter d (um)", np.log10(grids["r"]), "viridis", "log10 d"),
        ("K(d) - habitat", grids["K"], "viridis", "K"),
        ("D used (0 where dry)", grids["D"], "viridis", "D"),
        ("water-filled mask", grids["water_mask"].astype(float), "Blues", "wet"),
        ("OM (log10)", np.log10(grids["OM"] + 1e-12), "magma", "log10 OM"),
    ]
    for ax, (name, field, cmap, cbl) in zip(axes, fields):
        im = ax.imshow(field, cmap=cmap, origin="lower")
        ax.set_title(name, fontsize=10)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=cbl)
    wet_hab = metrics.wet_habitat_fraction(grids["water_mask"], grids["K"])
    fig.suptitle(f"{LABEL[soil]} at theta={grids['theta']:.2f} (n={grids['n']}) -- "
                 f"habitat={metrics.habitat_fraction(grids['K']):.3f}, "
                 f"WET habitat={wet_hab:.3f}")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, f"{tag}_water_map_{soil}.png"), dpi=150)
    plt.close(fig)


def plot_experiment_overlay(data, results_dir, fname="stage35_experiment_overlay.png",
                            subtitle=""):
    """SS5: model vs measured TEMPORAL SHAPE, per soil. Both are normalized to
    their own P1, so what is compared is the shape (rising / falling /
    erratic), not the absolute mg CO2 -- which the model has no units for."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharex=True)
    x = np.arange(4)
    for ax, s in zip(axes, SOILS):
        out, _ = data[s]
        model = normalized_periods(period_means(out["R_t"]))
        exp = normalized_periods(EXPERIMENT[s])
        ax.plot(x, exp, "o-", color="black", lw=2, ms=7, label="experiment (measured)")
        ax.plot(x, model, "s--", color=COLORS[s], lw=2, ms=7, label="model")
        ax.axhline(1.0, color="gray", lw=0.6, ls=":")
        ax.set_xticks(x)
        ax.set_xticklabels(PERIOD_NAMES)
        ax.set_ylim(bottom=0)
        ax.set_title(f"{LABEL[s]}\nexperiment: {shape_of_periods(EXPERIMENT[s])}  |  "
                     f"model: {shape_of_periods(period_means(out['R_t']))}", fontsize=11)
        ax.set_xlabel("period (3 / 4 / 3 / 5 days)")
        ax.legend(fontsize=8)
    axes[0].set_ylabel("respiration, normalized to each soil's own P1")
    fig.suptitle("Stage 3.5: model vs measured respiration -- temporal shape only" +
                 (f"\n{subtitle}" if subtitle else ""))
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, fname), dpi=150)
    plt.close(fig)


# ---------------------------------------- 5. theta sensitivity (secondary)

def theta_sensitivity(biology_cfg, mapping_cfg, results_dir):
    """SS4 secondary check: is the ranking/shape knife-edge on the exact theta?
    Runs the nearby fixed values, plus a wider diagnostic sweep that exposes
    WHERE each soil's habitat actually becomes wet."""
    sweep = {}
    for theta in THETA_DIAGNOSTIC:
        for s in SOILS:
            sweep[(s, theta)] = run_one(s, constant_theta(theta), GRID_N,
                                        biology_cfg, mapping_cfg)

    path = os.path.join(results_dir, "stage35_theta_sensitivity.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["theta", "soil", "d_cut_um", "habitat_fraction",
                    "wet_habitat_fraction", "habitat_wet_share", "R_peak",
                    "cum_co2_final", "outcome", "period_shape_model",
                    "isolation_ratio", "n_clusters"])
        for theta in THETA_DIAGNOSTIC:
            for s in SOILS:
                out, conn = sweep[(s, theta)]
                d_cut = float(np.quantile(out["grids"]["r"], theta))
                w.writerow([theta, LABEL[s], d_cut, conn["habitat_fraction"],
                            conn["wet_habitat_fraction"], conn["habitat_wet_share"],
                            float(out["R_t"].max()), float(out["cum_co2"][-1]),
                            metrics.classify_outcome(out),
                            shape_of_periods(period_means(out["R_t"])),
                            conn["isolation_ratio"], conn["n_clusters"]])

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    for s in SOILS:
        wet_hab = [sweep[(s, th)][1]["wet_habitat_fraction"] for th in THETA_DIAGNOSTIC]
        cums = [sweep[(s, th)][0]["cum_co2"][-1] for th in THETA_DIAGNOSTIC]
        d_cut = [float(np.quantile(sweep[(s, th)][0]["grids"]["r"], th)) for th in THETA_DIAGNOSTIC]
        axes[0].plot(THETA_DIAGNOSTIC, d_cut, "o-", color=COLORS[s], label=LABEL[s])
        axes[1].plot(THETA_DIAGNOSTIC, wet_hab, "o-", color=COLORS[s], label=LABEL[s])
        axes[2].plot(THETA_DIAGNOSTIC, cums, "o-", color=COLORS[s], label=LABEL[s])

    axes[0].axhspan(HAB_LO, HAB_HI, color="gray", alpha=0.15, label="habitat window")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("wetting cutoff diameter d_cut (um)")
    axes[0].set_title("Which pores are wet at each theta\n(small pores fill first)", fontsize=10)
    axes[1].set_ylabel("wet habitat fraction (wet AND K>0)")
    axes[1].set_title("Fraction of cells that are BOTH wet and habitable", fontsize=10)
    axes[2].set_ylabel("cumulative CO2 (final)")
    axes[2].set_title("Total respiration vs theta", fontsize=10)
    for ax in axes:
        ax.axvline(THETA_CONST, color="k", ls=":", lw=1)
        ax.set_xlabel("theta (constant, shared by all soils)")
        ax.legend(fontsize=8)
    fig.suptitle("Stage 3.5 theta-sensitivity: the shared-theta rule wets the Vertisol's and "
                 "Loess's HABITAT last (dotted line = the headline theta)")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage35_theta_sensitivity.png"), dpi=150)
    plt.close(fig)

    # the nearby-theta panel the prompt actually asks for (SS4)
    fig, axes = plt.subplots(1, len(THETA_SENSITIVITY), figsize=(5 * len(THETA_SENSITIVITY), 4.4),
                             sharey=True)
    for ax, theta in zip(axes, THETA_SENSITIVITY):
        for s in SOILS:
            out, _ = sweep[(s, theta)]
            t = np.arange(len(out["R_t"])) * out["dt"]
            ax.plot(t, out["R_t"], color=COLORS[s], label=LABEL[s], lw=1.7)
        ax.set_title(f"theta = {theta:.1f}")
        ax.set_xlabel("model time")
        ax.legend(fontsize=8)
    axes[0].set_ylabel("R(t)")
    fig.suptitle("Stage 3.5 robustness: R(t) at nearby fixed theta (is the ranking knife-edge?)")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage35_theta_robustness.png"), dpi=150)
    plt.close(fig)

    return sweep


# ------------------------- 6. shared-matric-potential diagnostic (see SS9)

def matric_diagnostic(biology_cfg, mapping_cfg, results_dir):
    """NOT the headline run. The headline (SS4) imposes the SAME theta on all
    three soils. But a closed humid jar imposes the same MATRIC POTENTIAL, and
    three soils at one matric potential sit at three different theta -- that is
    exactly what a water-retention curve is. Under a shared theta the Vertisol's
    habitat macropores are the LAST pores to wet (68% of its volume is
    sub-micron), so they are bone dry at theta=0.6 and its biomass is cut off
    from substrate. This run holds the air-entry diameter (mapping.yaml
    matric_d_cut) fixed instead, lets each soil's theta emerge from its own PSD,
    and is reported alongside -- never instead of -- the headline."""
    d_cut = float(mapping_cfg["matric_d_cut"])
    sat = {"mode": "matric_d", "d_cut": d_cut}
    data = {s: run_one(s, sat, GRID_N, biology_cfg, mapping_cfg) for s in SOILS}
    return data, d_cut


# --------------------------------------------- 7. distinctness (SS6)

def distinctness(data, label, results_dir, rows_accum):
    out_l, out_s, out_v = (data[s][0] for s in SOILS)
    success, detail = metrics.emergent_distinctness(
        out_l, out_s, out_v, require_nontrivial=False)
    d = detail["distances"]
    rows_accum.append({
        "run": label,
        "Loess_vs_Sand": d["AB"], "Loess_vs_Vertisol": d["AC"],
        "Sand_vs_Vertisol": d["BC"],
        "min_distance": detail["min_distance"],
        "threshold": detail["threshold"], "success": success,
        "outcome_Loess": detail["outcomes"]["A"],
        "outcome_Sand": detail["outcomes"]["B"],
        "outcome_Vertisol": detail["outcomes"]["C"],
    })
    return success, detail


# ------------------------------------------------------------------ main

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    biology_cfg = load_yaml("configs/biology.yaml")
    mapping_cfg = load_yaml("configs/mapping.yaml")

    print("=" * 78)
    print("STAGE 3.5 -- literature-informed PSDs, constant moisture")
    print("=" * 78)

    # 1. PSDs
    psds, psd_rows = psd_summary(mapping_cfg, RESULTS_DIR)
    print("\n[1] literature PSDs (SS2) -- habitable = volume fraction in 30-150 um:")
    print(f"    {'soil':10} {'porosity':>9} {'median um':>10} {'<30um':>7} {'habitable':>10} {'>150um':>8}")
    for r in psd_rows:
        print(f"    {r['soil']:10} {r['total_porosity']:9.2f} {r['median_diameter_um']:10.2f} "
              f"{r['vol_frac_below_30um']:7.3f} {r['habitable_vol_frac_30_150um']:10.3f} "
              f"{r['vol_frac_above_150um']:8.3f}")

    # 2. D(d) fusion fix
    d_rows, worst_stability = d_of_d_check(biology_cfg, mapping_cfg, RESULTS_DIR)
    print(f"\n[2] D(d) fusion fix (SS3) -- D_d_ref {d_rows[0]['D_d_ref_um']:.0f} -> "
          f"{d_rows[1]['D_d_ref_um']:.0f} um:")
    for r in d_rows:
        print(f"    D_d_ref={r['D_d_ref_um']:5.0f}: D(30)={r['D_at_30um']:.3f} "
              f"D(150)={r['D_at_150um']:.3f} D(300)={r['D_at_300um']:.3f}  "
              f"D(300)/D(30)={r['D300_over_D30']:5.2f}x  "
              f"coarse-pair separation={r['coarse_pair_separation']:.3f}x")
    print(f"    numerical stability: max(D)*dt/dx^2 = {worst_stability:.4f} < 0.2  OK")

    dist_rows = []

    # 3-4. headline constant-theta run + experiment overlay
    print(f"\n[3] HEADLINE: constant theta = {THETA_CONST} for all soils, n={GRID_N} (SS4)")
    data = headline_run(biology_cfg, mapping_cfg, RESULTS_DIR)
    rows = write_headline_csv(data, RESULTS_DIR)
    plot_curves(data, RESULTS_DIR,
                f"Stage 3.5 headline: constant theta = {THETA_CONST} (n={GRID_N})",
                "stage35_respiration_constant_theta.png")
    for s in SOILS:
        plot_water_map(s, data[s][0], RESULTS_DIR, "stage35")

    print(f"    {'soil':10} {'habitat':>8} {'WET hab':>8} {'iso ratio':>10} {'outcome':>10} "
          f"{'model shape':>12} {'experiment':>12}")
    for r in rows:
        print(f"    {r['soil']:10} {r['habitat_fraction']:8.3f} {r['wet_habitat_fraction']:8.3f} "
              f"{r['isolation_ratio']:10.4f} {r['outcome']:>10} "
              f"{r['period_shape_model']:>12} {r['period_shape_experiment']:>12}")

    matched = [r["soil"] for r in rows
               if r["period_shape_model"] == r["period_shape_experiment"]]
    plot_experiment_overlay(
        data, RESULTS_DIR,
        subtitle=f"HEADLINE: shared constant theta = {THETA_CONST}. "
                 f"Model reproduces the measured shape for: "
                 f"{', '.join(matched) if matched else 'NONE'}")
    print(f"\n[4] experiment overlay (SS5): model reproduces the measured shape for "
          f"{matched if matched else 'NONE of the three'}")

    ok, detail = distinctness(data, f"headline_constant_theta_{THETA_CONST}",
                              RESULTS_DIR, dist_rows)
    print(f"    distinctness (require_nontrivial=False): success={ok}, "
          f"min_distance={detail['min_distance']:.3f}, distances={ {k: round(v,3) for k,v in detail['distances'].items()} }")

    # optional smaller grid (SS6)
    data40 = {s: run_one(s, constant_theta(THETA_CONST), GRID_N_SMALL,
                         biology_cfg, mapping_cfg) for s in SOILS}
    write_headline_csv(data40, RESULTS_DIR, tag="stage35_constant_theta_n40")
    plot_curves(data40, RESULTS_DIR,
                f"Stage 3.5 grid check: constant theta = {THETA_CONST} (n={GRID_N_SMALL})",
                "stage35_respiration_constant_theta_n40.png")
    distinctness(data40, f"headline_constant_theta_{THETA_CONST}_n40", RESULTS_DIR, dist_rows)

    # 5. theta sensitivity
    print(f"\n[5] theta sensitivity (SS4, secondary): {THETA_SENSITIVITY} "
          f"(+ diagnostic sweep {THETA_DIAGNOSTIC})")
    sweep = theta_sensitivity(biology_cfg, mapping_cfg, RESULTS_DIR)
    print(f"    {'theta':>6} " + " ".join(f"{LABEL[s]:>22}" for s in SOILS))
    print(f"    {'':>6} " + " ".join(f"{'wet-hab / shape':>22}" for _ in SOILS))
    for th in THETA_DIAGNOSTIC:
        cells = []
        for s in SOILS:
            out, conn = sweep[(s, th)]
            cells.append(f"{conn['wet_habitat_fraction']:.3f} / "
                         f"{shape_of_periods(period_means(out['R_t'])):>8}")
        print(f"    {th:6.2f} " + " ".join(f"{c:>22}" for c in cells))

    # 6. shared-matric-potential diagnostic
    mdata, d_cut = matric_diagnostic(biology_cfg, mapping_cfg, RESULTS_DIR)
    mrows = write_headline_csv(mdata, RESULTS_DIR, tag="stage35_matric_diagnostic")
    plot_curves(mdata, RESULTS_DIR,
                f"Stage 3.5 DIAGNOSTIC: shared matric potential (air-entry d_cut = {d_cut:.0f} um),\n"
                f"theta emerges per soil -- NOT the headline run",
                "stage35_respiration_matric_diagnostic.png")
    for s in SOILS:
        plot_water_map(s, mdata[s][0], RESULTS_DIR, "stage35_matric")
    mmatched = [r["soil"] for r in mrows
                if r["period_shape_model"] == r["period_shape_experiment"]]
    plot_experiment_overlay(
        mdata, RESULTS_DIR, fname="stage35_experiment_overlay_matric.png",
        subtitle=f"DIAGNOSTIC (not the headline): shared MATRIC POTENTIAL, air-entry "
                 f"d_cut = {d_cut:.0f} um, theta emerges per soil. "
                 f"Reproduces: {', '.join(mmatched) if mmatched else 'NONE'}")
    ok_m, detail_m = distinctness(mdata, f"diagnostic_matric_d_cut_{d_cut:.0f}um",
                                  RESULTS_DIR, dist_rows)

    print(f"\n[6] DIAGNOSTIC: shared matric potential (air-entry d_cut = {d_cut:.0f} um) -- "
          f"theta EMERGES per soil")
    print(f"    {'soil':10} {'theta':>7} {'WET hab':>8} {'outcome':>10} {'model shape':>12} {'experiment':>12}")
    for r in mrows:
        print(f"    {r['soil']:10} {r['theta']:7.3f} {r['wet_habitat_fraction']:8.3f} "
              f"{r['outcome']:>10} {r['period_shape_model']:>12} {r['period_shape_experiment']:>12}")
    print(f"    reproduces the measured shape for {mmatched if mmatched else 'NONE of the three'}; "
          f"distinctness success={ok_m}, min_distance={detail_m['min_distance']:.3f}")

    # 7. distinctness table
    path = os.path.join(RESULTS_DIR, "stage35_distinctness.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(dist_rows[0]))
        w.writeheader()
        w.writerows(dist_rows)

    print(f"\n[7] all results written to {RESULTS_DIR}/")


if __name__ == "__main__":
    main()
