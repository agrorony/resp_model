"""
Stage 11 runner (prompts/soil_respiration_prompt_stage11_vertisol.md).

Continues from Stage 10. Stage 10 gave Sand a strong early burst-then-crash
(kept exactly as-is here -- Sand is FROZEN, no Sand parameter touched) but
its Vertisol only ever reached biomass PARITY with its initial total
(`B_peak/B0=1.000`, never actual growth) even though its R(t) shape already
read as a plausible slow rise. This stage: (1) fixes a real per-timestep
performance bug in `simulate_dual_porosity` (`habitat_properties` and the
K/habitable fields derived from it were recomputed every timestep despite
being loop-invariant), verified to change nothing; (2) raises Vertisol's
total OM and slows its transport magnitude further (Vertisol-only levers)
so its biomass clearly grows, not just clears a trivial floor comparison.
n=128 ONLY (no n=40 anywhere, per this stage's own restriction).

Usage:
    python stage11_run.py
"""
from __future__ import annotations

import copy
import csv
import datetime
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

RESULTS_DIR = os.path.join("results", "stage11")
LOGBOOK_PATH = "LOGBOOK.md"

GRID_N = 128
MATRIC_D_CUT = 50.0
MAINTENANCE_FLOOR = 0.4   # m0 * B0_total = 0.02 * 20
DISTINCTNESS_THRESHOLD = 0.3

# Stage 10's own recorded reference numbers for Sand (n=128, T=3000, the
# pre-Stage-11, fixed-90-iteration, per-timestep-recomputed baseline) --
# the identity check compares against these directly, not just "before vs
# after in this same run", so the check also catches any accidental drift
# introduced since Stage 10 was written up.
SAND_REFERENCE_R_PEAK = 15.794899991768297
SAND_REFERENCE_R_END = 0.29810401858904995


def _T_for(soil, mapping_cfg):
    return int(mapping_cfg[f"stage10_T_{soil}"])


def _structure_cfg(soil, mapping_cfg):
    cfg = copy.deepcopy(load_yaml(f"configs/soil_{soil}.yaml"))
    cfg["saturation"] = {"mode": "matric_regions", "matric_d_cut_um": MATRIC_D_CUT}
    cfg["grid"]["n"] = GRID_N
    cfg["T"] = _T_for(soil, mapping_cfg)
    return cfg


def run_soil(soil, biology_cfg, mapping_cfg):
    cfg = _structure_cfg(soil, mapping_cfg)
    t0 = time.time()
    out = dual_porosity.simulate_dual_porosity(biology_cfg, cfg, mapping_cfg, wicking_enabled=False)
    elapsed = time.time() - t0
    conn = metrics.connectivity_diagnostics(out["grids"])
    return out, conn, elapsed


def run_all(biology_cfg, mapping_cfg):
    results = {}
    for soil in SOILS:
        out, conn, elapsed = run_soil(soil, biology_cfg, mapping_cfg)
        results[soil] = (out, conn, elapsed)
    return results


# --------------------------------------------------- 1. performance + identity

def sand_identity_check(data):
    out, _, elapsed = data["sand"]
    R_peak = float(out["R_t"].max())
    R_end = float(out["R_t"][-1])
    peak_diff = abs(R_peak - SAND_REFERENCE_R_PEAK) / SAND_REFERENCE_R_PEAK
    end_diff = abs(R_end - SAND_REFERENCE_R_END) / SAND_REFERENCE_R_END
    row = {
        "R_peak_stage10_reference": SAND_REFERENCE_R_PEAK,
        "R_peak_stage11": R_peak,
        "R_peak_relative_diff": peak_diff,
        "R_end_stage10_reference": SAND_REFERENCE_R_END,
        "R_end_stage11": R_end,
        "R_end_relative_diff": end_diff,
        "elapsed_seconds_stage11": elapsed,
        "identical_within_numerical_noise": bool(peak_diff < 0.01 and end_diff < 0.01),
    }
    return row


# --------------------------------------------------------------- 2. figures

def plot_respiration_curves(data, results_dir):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for soil in SOILS:
        out, _, _ = data[soil]
        t = np.arange(len(out["R_t"])) * out["dt"]
        axes[0].plot(t, out["R_t"], color=COLORS[soil], label=LABEL[soil], lw=2)
    axes[0].axhline(MAINTENANCE_FLOOR, color="gray", ls=":", lw=1, label=f"maintenance floor ({MAINTENANCE_FLOOR})")
    axes[0].set_xlabel("model time (absolute)"); axes[0].set_ylabel("R(t)")
    axes[0].set_title("Respiration rate, own absolute time axis"); axes[0].legend(fontsize=8)

    for soil in SOILS:
        out, _, _ = data[soil]
        x = np.linspace(0.0, 1.0, len(out["R_t"]))
        axes[1].plot(x, out["R_t"], color=COLORS[soil], label=LABEL[soil], lw=2)
    axes[1].axhline(MAINTENANCE_FLOOR, color="gray", ls=":", lw=1)
    axes[1].set_xlabel("fraction of own observation window")
    axes[1].set_ylabel("R(t)")
    axes[1].set_title("Respiration rate, normalized time (0=start, 1=end)")
    axes[1].legend(fontsize=8)
    fig.suptitle("Stage 11: Sand frozen (burst-then-crash) vs Vertisol re-tuned (slow rise), n=128")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage11_respiration_curves.png"), dpi=150)
    plt.close(fig)


def plot_vertisol_diagnostics(data, results_dir):
    out = data["vertisol"][0]
    t = np.arange(len(out["R_t"])) * out["dt"]
    B0 = out["grids"]["B0"].sum()
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    axes[0, 0].plot(t, out["R_t"], color=COLORS["vertisol"], lw=1.8)
    axes[0, 0].axhline(MAINTENANCE_FLOOR, color="gray", ls=":", lw=1, label="maintenance floor")
    axes[0, 0].set_title("R(t)"); axes[0, 0].set_xlabel("model time"); axes[0, 0].legend(fontsize=8)

    axes[0, 1].plot(t, out["B_total_t"] / B0, color=COLORS["vertisol"], lw=1.8)
    axes[0, 1].axhline(1.0, color="gray", ls=":", lw=1, label="B0 (initial)")
    axes[0, 1].set_title("biomass(t) / B0"); axes[0, 1].set_xlabel("model time"); axes[0, 1].legend(fontsize=8)

    axes[1, 0].plot(t, out["habitats_recruited_t"], color=COLORS["vertisol"], lw=1.8)
    axes[1, 0].set_title("habitats_recruited(t)"); axes[1, 0].set_xlabel("model time")
    axes[1, 0].set_ylim(bottom=0)

    axes[1, 1].plot(t, out["S_habitat_mean_t"], color=COLORS["vertisol"], lw=1.8)
    axes[1, 1].set_title("mean_S_habitat(t)"); axes[1, 1].set_xlabel("model time")

    fig.suptitle("Stage 11: Vertisol diagnostics")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage11_vertisol_diagnostics.png"), dpi=150)
    plt.close(fig)


# --------------------------------------------------------- 3. headline CSV

def write_headline_csv(data, results_dir, mapping_cfg):
    rows = []
    for soil in SOILS:
        out, conn, elapsed = data[soil]
        grids = out["grids"]
        R = out["R_t"]
        n3 = len(R) // 3
        B0 = grids["B0"].sum()
        rec = out["habitats_recruited_t"]
        rows.append({
            "soil": LABEL[soil],
            "om_total": grids["om_total_init"],
            "coating_fraction": grids["coating_fraction"],
            "D_scale": dual_porosity._variantc_D_scale(mapping_cfg, soil),
            "T": mapping_cfg[f"stage10_T_{soil}"],
            "elapsed_seconds": elapsed,
            "R_peak": float(R.max()),
            "R_end": float(R[-1]),
            "R_end_over_peak": float(R[-1] / R.max()),
            "peak_frac": float(np.argmax(R) / (len(R) - 1)),
            "early_third_mean": float(R[:n3].mean()),
            "late_third_mean": float(R[-n3:].mean()),
            "R_peak_over_maintenance_floor": float(R.max() / MAINTENANCE_FLOOR),
            "B_peak_over_B0": float(out["B_total_t"].max() / B0),
            "B_final_over_B0": float(out["B_final"].sum() / B0),
            "fuel_remaining_final": float(out["fuel_remaining_t"][-1]),
            "habitats_recruited_final": float(rec[-1]),
            "habitats_recruited_at_2pct_window": float(rec[max(int(0.02 * len(rec)), 0)]),
            "n_regions": int(grids["region_diam_um"].size),
            "outcome": metrics.classify_outcome(out),
        })
    path = os.path.join(results_dir, "stage11_headline.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


def write_identity_csv(row, results_dir):
    path = os.path.join(results_dir, "stage11_sand_identity_check.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(row))
        w.writeheader()
        w.writerow(row)
    return row


def write_distinctness_csv(data, results_dir):
    r_sand = data["sand"][0]["R_t"]
    r_vert = data["vertisol"][0]["R_t"]
    x_sand = np.linspace(0.0, 1.0, len(r_sand))
    x_vert = np.linspace(0.0, 1.0, len(r_vert))
    n_points = 500
    x_new = np.linspace(0.0, 1.0, n_points)
    r_sand_rs = np.interp(x_new, x_sand, r_sand)
    r_vert_rs = np.interp(x_new, x_vert, r_vert)
    d = metrics.pairwise_distance(r_sand_rs, r_vert_rs)
    row = {
        "Sand_vs_Vertisol_distance": d,
        "threshold": DISTINCTNESS_THRESHOLD,
        "distinct": d > DISTINCTNESS_THRESHOLD,
        "note": "computed on time-normalized (fraction-of-own-window) resampled curves",
    }
    path = os.path.join(results_dir, "stage11_distinctness.csv")
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
    print("STAGE 11 -- freeze Sand, fix Vertisol, optimize the loop (n=128 only)")
    print("=" * 78)

    print(f"\n[1] two-soil run at n={GRID_N} (Sand frozen from Stage 10; "
          f"Vertisol: om={mapping_cfg['stage11_om_total_vertisol']}, "
          f"coat={mapping_cfg['stage11_coating_fraction_vertisol']}, "
          f"D_scale={mapping_cfg['stage11_D_scale_vertisol']}) ...")
    data = run_all(biology_cfg, mapping_cfg)

    identity_row = sand_identity_check(data)
    write_identity_csv(identity_row, RESULTS_DIR)
    print(f"\n[2] Sand identity check (performance fix must not alter results):")
    print(f"    R_peak: Stage10={identity_row['R_peak_stage10_reference']:.6f}  "
          f"Stage11={identity_row['R_peak_stage11']:.6f}  "
          f"rel.diff={identity_row['R_peak_relative_diff']:.2e}")
    print(f"    R_end:  Stage10={identity_row['R_end_stage10_reference']:.6f}  "
          f"Stage11={identity_row['R_end_stage11']:.6f}  "
          f"rel.diff={identity_row['R_end_relative_diff']:.2e}")
    print(f"    identical within numerical noise: {identity_row['identical_within_numerical_noise']}")
    print(f"    Sand elapsed (with fix): {identity_row['elapsed_seconds_stage11']:.2f}s")

    rows = write_headline_csv(data, RESULTS_DIR, mapping_cfg)
    print(f"\n[3] headline:")
    for r in rows:
        print(f"    {r['soil']:10} R_peak={r['R_peak']:8.3f}  R_peak/floor={r['R_peak_over_maintenance_floor']:6.2f}x  "
              f"B_peak/B0={r['B_peak_over_B0']:6.3f}  R_end/pk={r['R_end_over_peak']:.3f}  "
              f"peak_frac={r['peak_frac']:.3f}  elapsed={r['elapsed_seconds']:.2f}s  outcome={r['outcome']}")

    plot_respiration_curves(data, RESULTS_DIR)
    plot_vertisol_diagnostics(data, RESULTS_DIR)

    dist_row = write_distinctness_csv(data, RESULTS_DIR)
    print(f"\n[4] distinctness (time-normalized): distance={dist_row['Sand_vs_Vertisol_distance']:.4f}, "
          f"distinct={dist_row['distinct']}")

    vert_row = rows[SOILS.index("vertisol")]
    vertisol_success = (
        vert_row["late_third_mean"] > vert_row["early_third_mean"]
        and vert_row["R_peak_over_maintenance_floor"] > 1.0
        and vert_row["B_peak_over_B0"] > 1.0
        and vert_row["peak_frac"] > 0.6
        and vert_row["R_end_over_peak"] > 0.7
    )
    overall = identity_row["identical_within_numerical_noise"] and vertisol_success and dist_row["distinct"]

    with open(LOGBOOK_PATH, "a", encoding="utf-8") as f:
        f.write(f"\n**Stage 11 headline (n={GRID_N}): Sand identical to Stage 10 within numerical noise="
                f"{identity_row['identical_within_numerical_noise']} (elapsed {identity_row['elapsed_seconds_stage11']:.1f}s), "
                f"Vertisol slow-rise success={vertisol_success} "
                f"(R_peak={vert_row['R_peak']:.2f} = {vert_row['R_peak_over_maintenance_floor']:.1f}x floor, "
                f"B_peak/B0={vert_row['B_peak_over_B0']:.2f}), distinct={dist_row['distinct']} "
                f"(distance={dist_row['Sand_vs_Vertisol_distance']:.4f}). Overall success={overall}. "
                f"{datetime.datetime.now().isoformat(timespec='seconds')}. Full detail in "
                "docs/stage_results/STAGE11_RESULTS.md.**\n")

    print(f"\n[5] OVERALL: Sand identical={identity_row['identical_within_numerical_noise']}, "
          f"Vertisol slow-rise={vertisol_success}, distinct={dist_row['distinct']} -> success={overall}")
    print(f"\n[6] all results written to {RESULTS_DIR}/")


if __name__ == "__main__":
    main()
