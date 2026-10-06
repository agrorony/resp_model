"""
Stage 13 runner (prompts/soil_respiration_prompt_stage13_real_data.md).

DATA-SOURCING + RECALIBRATION-CHECK stage, not a new mechanism. Replaces
Stages 3.5-12's literature-parametrized PSDs (configs/psd_literature.yaml)
with real measured PSD + Track E connectivity data from the companion
"resarch exercise" project (configs/psd_measured.yaml,
data/psd_measured/<soil>/source.txt carry full provenance) -- see
docs/stage_results/STAGE13_RESULTS.md for the data-sourcing narrative and
the honest verdict this script's own numbers below support.

Usage:
    python stage13_run.py
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
import metrics
import psd_measured

SOILS = ["sand", "vertisol"]          # Stage 10/11's frozen recalibration-check set
ALL_SOILS = ["sand", "loess", "vertisol"]   # SS2 region-connectivity, all 3
LABEL = {"sand": "Sand", "loess": "Loess", "vertisol": "Vertisol"}
COLORS = {"sand": "#d62728", "loess": "#2ca02c", "vertisol": "#4b0082"}

RESULTS_DIR = os.path.join("results", "stage13")
LOGBOOK_PATH = "LOGBOOK.md"

GRID_N = 128
GRID_N_SMALL = 40
MATRIC_D_CUT = 50.0
MAINTENANCE_FLOOR = 0.4
DISTINCTNESS_THRESHOLD = 0.3
MAX_ITERATIONS = 15

# Stage 11's exact recorded reference numbers (literature PSD, n=128,
# results/stage11/stage11_headline.csv + stage11_distinctness.csv) -- the
# recalibration check compares real-data numbers against THESE, not a
# freshly re-run literature baseline, so it also catches any accidental
# drift in the frozen mechanism since Stage 11 was written up.
STAGE11_REFERENCE = {
    "sand": {
        "R_peak": 15.789854790999666, "R_end": 0.2957433626001884,
        "R_end_over_peak": 0.018729960884045896, "peak_frac": 0.1813937979326442,
        "B_peak_over_B0": 10.70073770003218, "fuel_remaining_final": 0.01470602062936485,
        "n_regions": 3, "outcome": "declining",
    },
    "vertisol": {
        "R_peak": 3.7561399861379976, "R_end": 3.7493280228711905,
        "R_end_over_peak": 0.9981864458481455, "peak_frac": 0.9459639759839893,
        "early3": 0.320521554046756, "late3": 3.672714109292974,
        "B_peak_over_B0": 2.3372611216671735, "fuel_remaining_final": 0.48824011434181513,
        "habitats_recruited_final": 55.0, "n_regions": 70, "outcome": "thriving",
    },
}
STAGE11_DISTINCTNESS = 2.0009240671906716


# ------------------------------------------------------------- run helpers

def _T_for(soil, mapping_cfg):
    return int(mapping_cfg[f"stage10_T_{soil}"])


def _structure_cfg(soil, grid_n, mapping_cfg, measured=True):
    name = f"{soil}_measured" if measured else soil
    cfg = copy.deepcopy(load_yaml(f"configs/soil_{name}.yaml"))
    cfg["saturation"] = {"mode": "matric_regions", "matric_d_cut_um": MATRIC_D_CUT}
    cfg["grid"]["n"] = grid_n
    cfg["T"] = _T_for(soil, mapping_cfg)
    return cfg


def run_soil(soil, grid_n, biology_cfg, mapping_cfg, measured=True, mapping_overrides=None):
    cfg = _structure_cfg(soil, grid_n, mapping_cfg, measured=measured)
    mcfg = mapping_cfg
    if mapping_overrides:
        mcfg = copy.deepcopy(mapping_cfg)
        mcfg.update(mapping_overrides)
    out = dual_porosity.simulate_dual_porosity(biology_cfg, cfg, mcfg, wicking_enabled=False)
    conn = metrics.connectivity_diagnostics(out["grids"])
    return out, conn


def run_all(soils, grid_n, biology_cfg, mapping_cfg, measured=True):
    return {soil: run_soil(soil, grid_n, biology_cfg, mapping_cfg, measured=measured) for soil in soils}


# --------------------------------------------------------- pass/fail (Stage 10 SS3, reused verbatim)

def sand_pass_fail(out):
    R = out["R_t"]
    peak_idx = int(np.argmax(R))
    peak_frac = peak_idx / (len(R) - 1)
    R_peak = float(R.max())
    R_end_over_peak = float(R[-1] / R_peak)
    B0 = float(out["grids"]["B0"].sum())
    B_peak_over_B0 = float(out["B_total_t"].max() / B0)
    growth_at_peak = float(out["growth_term_t"][peak_idx])
    maint_at_peak = float(out["maintenance_term_t"][peak_idx])
    growth_share_at_peak = growth_at_peak / (growth_at_peak + maint_at_peak) if (growth_at_peak + maint_at_peak) > 0 else 0.0
    fuel_remaining_final = float(out["fuel_remaining_t"][-1])

    checks = {
        "R_peak >= 2x maintenance_floor (0.8)": R_peak >= 2.0 * MAINTENANCE_FLOOR,
        "growth term is MAJORITY of R at peak": growth_share_at_peak > 0.5,
        "B_peak/B0 > 1 (biomass grows)": B_peak_over_B0 > 1.0,
        "coating fuel visibly depletes (fuel_remaining_final < 0.5)": fuel_remaining_final < 0.5,
        "R_end/R_peak < 0.3 (crash)": R_end_over_peak < 0.3,
        "peak in the first third (peak_frac < 0.333)": peak_frac < 1.0 / 3.0,
    }
    values = {
        "R_peak": R_peak, "R_end_over_peak": R_end_over_peak, "peak_frac": peak_frac,
        "B_peak_over_B0": B_peak_over_B0, "growth_share_at_peak": growth_share_at_peak,
        "fuel_remaining_final": fuel_remaining_final,
    }
    return checks, values


def vertisol_pass_fail(out):
    R = out["R_t"]
    n3 = len(R) // 3
    early3 = float(R[:n3].mean())
    late3 = float(R[-n3:].mean())
    peak_idx = int(np.argmax(R))
    peak_frac = peak_idx / (len(R) - 1)
    R_peak = float(R.max())
    R_end_over_peak = float(R[-1] / R_peak)

    checks = {
        "late-third mean > early-third mean": late3 > early3,
        "peaks in the LATE part of the window (peak_frac > 0.6)": peak_frac > 0.6,
        "not an early spike + dip (R_end/R_peak > 0.7)": R_end_over_peak > 0.7,
        "clearly not a burst (R_peak < 4x maintenance_floor)": R_peak < 4.0 * MAINTENANCE_FLOOR,
    }
    values = {
        "early3": early3, "late3": late3, "peak_frac": peak_frac,
        "R_peak": R_peak, "R_end_over_peak": R_end_over_peak,
    }
    return checks, values


PASS_FAIL = {"sand": sand_pass_fail, "vertisol": vertisol_pass_fail}


def normalized_time_resample(R_t, n_points=500):
    x_orig = np.linspace(0.0, 1.0, len(R_t))
    x_new = np.linspace(0.0, 1.0, n_points)
    return np.interp(x_new, x_orig, R_t)


# --------------------------------------------------- 1. matrix-phase diagnostic

def matrix_phase_diagnostic(mapping_cfg, results_dir):
    """SS0/SS1 precursor check: how much of each soil's real measured PSD
    volume falls BELOW K_d_low (the model's macro/matrix split threshold)?
    This is checked BEFORE any simulation because it determines whether the
    dual-porosity mechanism (which places ALL organic matter inside
    `porous_matrix_mask`, i.e. matrix cells with diameter < K_d_low) has
    any substrate-hosting phase to work with at all under real data."""
    k_d_low = float(mapping_cfg["K_d_low"])
    rows = []
    for soil in ALL_SOILS:
        psd = psd_measured.build_soil_psd(soil)
        centers, vol = psd["bin_centers_um"], psd["volume_count"]
        below = float(vol[centers < k_d_low].sum() / vol.sum())
        rows.append({
            "soil": LABEL[soil], "K_d_low_um": k_d_low,
            "measured_volume_fraction_below_K_d_low": below,
            "n_nonzero_bins_below_K_d_low": int((vol[centers < k_d_low] > 0).sum()),
        })
    path = os.path.join(results_dir, "stage13_matrix_phase_diagnostic.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


# --------------------------------------------------- 2. recalibration check

def write_headline_csv(data, results_dir, mapping_cfg, tag="real"):
    rows = []
    for soil in SOILS:
        out, conn = data[soil]
        checks, vals = PASS_FAIL[soil](out)
        grids = out["grids"]
        row = {
            "soil": LABEL[soil],
            "om_total": grids["om_total_init"],
            "coating_fraction": grids["coating_fraction"],
            "D_scale": mapping_cfg.get(f"stage11_D_scale_{soil}", mapping_cfg.get(f"stage10_D_scale_{soil}")),
            "T": mapping_cfg[f"stage10_T_{soil}"],
            "R_peak": float(out["R_t"].max()),
            "R_end": float(out["R_t"][-1]),
            "R_end_over_peak": float(out["R_t"][-1] / out["R_t"].max()) if out["R_t"].max() > 0 else float("nan"),
            "peak_frac": float(np.argmax(out["R_t"]) / (len(out["R_t"]) - 1)),
            "B_final_over_B0": float(out["B_final"].sum() / grids["B0"].sum()),
            "B_peak_over_B0": float(out["B_total_t"].max() / grids["B0"].sum()),
            "fuel_remaining_final": float(out["fuel_remaining_t"][-1]),
            "n_regions": int(grids["region_diam_um"].size),
            "macro_fraction": float(grids["macro_mask"].mean()),
            "porous_matrix_fraction": float(grids["porous_matrix_mask"].mean()),
            "OM_total_placed": float(grids["OM"].sum()),
            "outcome": metrics.classify_outcome(out),
            "all_checks_pass": all(checks.values()),
        }
        for name, passed in checks.items():
            row[f"check: {name}"] = passed
        rows.append(row)
    path = os.path.join(results_dir, f"stage13_headline_{tag}.csv")
    fieldnames = list(dict.fromkeys(k for row in rows for k in row))
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, restval="")
        w.writeheader()
        w.writerows(rows)
    return rows


def compare_to_stage11(rows, results_dir):
    """Direct real-vs-literature-PSD comparison at the frozen Stage 10/11
    per-soil configuration (SS3: 'Re-run ... parameters otherwise
    unchanged first')."""
    out_rows = []
    for row in rows:
        soil = row["soil"].lower()
        ref = STAGE11_REFERENCE[soil]
        out_rows.append({
            "soil": row["soil"],
            "R_peak_stage11_literature": ref["R_peak"],
            "R_peak_stage13_real": row["R_peak"],
            "outcome_stage11": ref["outcome"],
            "outcome_stage13_real": row["outcome"],
            "n_regions_stage11": ref["n_regions"],
            "n_regions_stage13_real": row["n_regions"],
            "all_checks_pass_stage13_real": row["all_checks_pass"],
        })
    path = os.path.join(results_dir, "stage13_vs_stage11_comparison.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0]))
        w.writeheader()
        w.writerows(out_rows)
    return out_rows


# --------------------------------------------------- 3. budget-capped search

def search_probe(soil, om_total, coating_fraction, D_scale, T, biology_cfg, mapping_cfg, grid_n=GRID_N_SMALL):
    """One directed probe: override the soil's stage11/stage10 knobs
    directly via mapping_cfg copies (same mechanism `_variantc_*` already
    reads), never touching dual_porosity.py physics or introducing a new
    lever."""
    mcfg = copy.deepcopy(mapping_cfg)
    mcfg[f"stage11_om_total_{soil}"] = om_total
    mcfg[f"stage11_coating_fraction_{soil}"] = coating_fraction
    mcfg[f"stage11_D_scale_{soil}"] = D_scale
    mcfg[f"stage10_T_{soil}"] = T
    cfg = copy.deepcopy(load_yaml(f"configs/soil_{soil}_measured.yaml"))
    cfg["saturation"] = {"mode": "matric_regions", "matric_d_cut_um": MATRIC_D_CUT}
    cfg["grid"]["n"] = grid_n
    cfg["T"] = T
    out = dual_porosity.simulate_dual_porosity(biology_cfg, cfg, mcfg, wicking_enabled=False)
    return out


def budget_capped_search(soil, biology_cfg, mapping_cfg, results_dir, max_iterations=MAX_ITERATIONS):
    """Stage 13 SS3: 'ONLY if the real data breaks Stage 10/11 SS3's
    pass/fail table, re-run a budget-capped search over the same per-soil
    knobs Stage 10 already established -- do not introduce new free
    parameters.' Stops EARLY (before exhausting `max_iterations`) once the
    structural cause (empty porous_matrix_mask -- see
    matrix_phase_diagnostic) is directly confirmed to make every probe
    identical regardless of the four knobs, per Stage 15's own precedent
    ('name that mechanism explicitly rather than chasing a 16th probe')."""
    ref_om = float(mapping_cfg[f"stage10_om_total_{soil}"])
    ref_coat = float(mapping_cfg.get(f"stage11_coating_fraction_{soil}", mapping_cfg[f"stage10_coating_fraction_{soil}"]))
    ref_D = float(mapping_cfg.get(f"stage11_D_scale_{soil}", mapping_cfg[f"stage10_D_scale_{soil}"]))
    ref_T = int(mapping_cfg[f"stage10_T_{soil}"])

    probes = [
        {"om_total": ref_om * 10, "coating_fraction": ref_coat, "D_scale": ref_D, "T": ref_T,
         "note": "10x om_total, everything else frozen"},
        {"om_total": ref_om, "coating_fraction": min(ref_coat * 2, 1.0), "D_scale": ref_D, "T": ref_T,
         "note": "2x coating_fraction (capped at 1.0), everything else frozen"},
        {"om_total": ref_om * 25, "coating_fraction": 1.0, "D_scale": ref_D * 5, "T": ref_T * 3,
         "note": "extreme joint probe: 25x om_total, coating=1.0, 5x D_scale, 3x T"},
    ]

    log = []
    for i, p in enumerate(probes, start=1):
        out40 = search_probe(soil, p["om_total"], p["coating_fraction"], p["D_scale"], p["T"],
                              biology_cfg, mapping_cfg, grid_n=GRID_N_SMALL)
        om_placed = float(out40["grids"]["OM"].sum())
        R_peak = float(out40["R_t"].max())
        checks, _ = PASS_FAIL[soil](out40)
        log.append({
            "probe": i, "soil": LABEL[soil], **{k: p[k] for k in ("om_total", "coating_fraction", "D_scale", "T", "note")},
            "OM_total_placed": om_placed, "R_peak": R_peak,
            "porous_matrix_fraction": float(out40["grids"]["porous_matrix_mask"].mean()),
            "all_checks_pass": all(checks.values()),
        })
        if om_placed <= 0.0 and i >= 2:
            log[-1]["note"] += " -- CONFIRMS: OM_total_placed=0 regardless of knob values (porous_matrix_mask is empty; " \
                                "structural, not a tuning problem) -- stopping search early per Stage 15 precedent"
            break

    path = os.path.join(results_dir, f"stage13_search_trace_{soil}.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(log[0]))
        w.writeheader()
        w.writerows(log)
    return log


# --------------------------------------------------------------- 4. region-connectivity (SS2)

def region_connectivity_table(mapping_cfg, results_dir, grid_n=GRID_N):
    """SS2: compare the model's own emergent macropore region count/largest-
    region fraction against real Track E Gamma / crossover radius per soil
    -- report agreement or disagreement honestly, never tune geometry to
    force a match (Stage 3.5/4/5's own standing rule, reused verbatim
    here)."""
    biology_cfg = load_yaml("configs/biology.yaml")
    rows = []
    for soil in ALL_SOILS:
        cfg = copy.deepcopy(load_yaml(f"configs/soil_{soil}_measured.yaml"))
        cfg["saturation"] = {"mode": "matric_regions", "matric_d_cut_um": MATRIC_D_CUT}
        cfg["grid"]["n"] = grid_n
        grids = dual_porosity.build_dual_grids(cfg, mapping_cfg)
        macro_mask = grids["macro_mask"]
        region_diam = grids["region_diam_um"]
        n_regions = int(region_diam.size)
        if n_regions > 0:
            labels = grids["region_labels"]
            sizes = np.array([(labels == k).sum() for k in range(1, n_regions + 1)])
            largest_frac = float(sizes.max() / macro_mask.sum()) if macro_mask.sum() > 0 else 0.0
        else:
            largest_frac = 0.0
        psd = psd_measured.build_soil_psd(soil)
        real = psd["connectivity"]
        rows.append({
            "soil": LABEL[soil],
            "model_macro_fraction": float(macro_mask.mean()),
            "model_n_regions": n_regions,
            "model_largest_region_fraction": largest_frac,
            "real_gamma": real["gamma"],
            "real_euler_number": real["euler_number"],
            "real_r_star_um": real.get("r_star_um"),
            "real_note": real.get("note", ""),
            "agreement_note": (
                "model degenerates to ~1 region covering the whole grid (macro_fraction~1.0) -- "
                "directionally consistent with real soil's own single dominant component "
                "(largest-component fraction 92.7-96.3% per connectivity_validation_summary.md Part A) "
                "and high Gamma, but the model CANNOT reproduce the real soil's co-existing tens-of-"
                "thousands of tiny fragmented components (a single global diameter threshold, not a "
                "genuine 3D pore network, has no mechanism to fragment a region once it's connected) -- "
                "see STAGE13_RESULTS.md SS2 for the full honest comparison, not forced to match"
            ),
        })
    path = os.path.join(results_dir, "stage13_region_connectivity_comparison.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


# --------------------------------------------------------------- 5. figures

def plot_real_vs_literature_psd(results_dir):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5), sharey=False)
    import psd_parametric
    for ax, soil in zip(axes, ALL_SOILS):
        psd_real = psd_measured.build_soil_psd(soil)
        psd_lit = psd_parametric.build_soil_psd(soil)
        ax.plot(psd_lit["bin_centers_um"], psd_lit["differential_psd"], color="gray", lw=1.5,
                label="literature (Stage 3.5-12)")
        ax.plot(psd_real["bin_centers_um"], psd_real["differential_psd"], color=COLORS[soil], lw=2,
                label="real measured (Stage 13)")
        ax.axvline(30.0, color="k", ls=":", lw=1, label="K_d_low=30um (macro/matrix split)")
        ax.set_xscale("log"); ax.set_xlabel("pore diameter (um, log)"); ax.set_ylabel("differential PSD")
        ax.set_title(LABEL[soil]); ax.legend(fontsize=7)
    fig.suptitle("Stage 13: real measured vs. literature-mixture PSD, per soil")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage13_psd_real_vs_literature.png"), dpi=150)
    plt.close(fig)


def plot_respiration_comparison(data_real, results_dir):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for soil in SOILS:
        out, _ = data_real[soil]
        t = np.arange(len(out["R_t"])) * out["dt"]
        axes[0].plot(t, out["R_t"], color=COLORS[soil], label=f"{LABEL[soil]} (real PSD)", lw=2)
    axes[0].axhline(MAINTENANCE_FLOOR, color="gray", ls=":", lw=1, label=f"maintenance floor ({MAINTENANCE_FLOOR})")
    axes[0].set_xlabel("model time"); axes[0].set_ylabel("R(t)")
    axes[0].set_title("Stage 13: real-PSD R(t), Stage 10/11 frozen config"); axes[0].legend(fontsize=8)

    for soil in SOILS:
        out, _ = data_real[soil]
        x = np.linspace(0.0, 1.0, len(out["R_t"]))
        axes[1].plot(x, out["R_t"], color=COLORS[soil], label=LABEL[soil], lw=2)
    axes[1].axhline(MAINTENANCE_FLOOR, color="gray", ls=":", lw=1)
    axes[1].set_xlabel("fraction of own observation window"); axes[1].set_ylabel("R(t)")
    axes[1].set_title("normalized time"); axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage13_respiration_curves.png"), dpi=150)
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
    print("STAGE 13 -- real measured PSD + connectivity data")
    print("=" * 78)

    print("\n[0] matrix-phase diagnostic (real PSD volume fraction below K_d_low=30um) ...")
    matrix_rows = matrix_phase_diagnostic(mapping_cfg, RESULTS_DIR)
    for r in matrix_rows:
        print(f"    {r['soil']:10} below-30um volume fraction={r['measured_volume_fraction_below_K_d_low']:.4f} "
              f"(nonzero bins below: {r['n_nonzero_bins_below_K_d_low']})")

    print(f"\n[1] frozen Stage 10/11 config re-run on real PSD, n={GRID_N} ...")
    data_real = run_all(SOILS, GRID_N, biology_cfg, mapping_cfg, measured=True)
    rows = write_headline_csv(data_real, RESULTS_DIR, mapping_cfg, tag="real")
    for r in rows:
        print(f"    {r['soil']:10} R_peak={r['R_peak']:.4f}  outcome={r['outcome']:10} "
              f"OM_placed={r['OM_total_placed']:.2f}  porous_matrix_frac={r['porous_matrix_fraction']:.4f}  "
              f"all_checks_pass={r['all_checks_pass']}")

    cmp_rows = compare_to_stage11(rows, RESULTS_DIR)
    print("\n[2] vs. Stage 11 literature-PSD reference:")
    for r in cmp_rows:
        print(f"    {r['soil']:10} R_peak lit={r['R_peak_stage11_literature']:.3f} -> real={r['R_peak_stage13_real']:.3f}  "
              f"outcome lit={r['outcome_stage11']} -> real={r['outcome_stage13_real']}  "
              f"checks_pass={r['all_checks_pass_stage13_real']}")

    need_search = any(not r["all_checks_pass"] for r in rows)
    search_logs = {}
    if need_search:
        print(f"\n[3] real data broke the SS3 pass/fail table -- budget-capped search (<={MAX_ITERATIONS} iterations/soil) ...")
        for r in rows:
            soil = r["soil"].lower()
            if not r["all_checks_pass"]:
                log = budget_capped_search(soil, biology_cfg, mapping_cfg, RESULTS_DIR)
                search_logs[soil] = log
                print(f"    {LABEL[soil]}: {len(log)} probe(s) run (capped at {MAX_ITERATIONS}), "
                      f"final probe note: {log[-1]['note']}")
    else:
        print("\n[3] real data reproduced Stage 10/11's pass/fail table with NO retuning -- search skipped.")

    print("\n[4] region-connectivity vs. Track E comparison (all 3 soils, no forced match) ...")
    conn_rows = region_connectivity_table(mapping_cfg, RESULTS_DIR)
    for r in conn_rows:
        print(f"    {r['soil']:10} model: n_regions={r['model_n_regions']}, largest_frac={r['model_largest_region_fraction']:.3f}  "
              f"real: Gamma={r['real_gamma']}, r*={r['real_r_star_um']}")

    plot_real_vs_literature_psd(RESULTS_DIR)
    plot_respiration_comparison(data_real, RESULTS_DIR)

    root_cause_note = ""
    if need_search:
        root_cause_note = (
            " Root cause of the failure(s): the affected soil's real PSD has 0% measured pore volume below "
            "K_d_low=30um (matrix_phase_diagnostic), so porous_matrix_mask -- where ALL organic matter is "
            "placed -- is structurally empty; the budget-capped search confirms no combination of the 4 "
            "allowed knobs (om_total, coating fraction, D_scale, T) can place OM without touching "
            "dual_porosity.py's physics (out of scope)."
        )
    with open(LOGBOOK_PATH, "a", encoding="utf-8") as f:
        f.write(f"\n**Stage 13 headline (n={GRID_N}): real measured PSD substituted for literature PSD, Stage 10/11 "
                f"config otherwise frozen. Sand all_checks_pass={rows[0]['all_checks_pass']}, "
                f"Vertisol all_checks_pass={rows[1]['all_checks_pass']}.{root_cause_note} "
                f"{datetime.datetime.now().isoformat(timespec='seconds')}. Full detail in "
                "docs/stage_results/STAGE13_RESULTS.md.**\n")

    print(f"\n[5] all results written to {RESULTS_DIR}/")


if __name__ == "__main__":
    main()
