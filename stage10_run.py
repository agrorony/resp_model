"""
Stage 10 (FINAL) runner (prompts/soil_respiration_prompt_stage10_final.md).

The last stage. Two soils only -- Vertisol and Sand (Loess dropped
entirely). Stage 9 got the qualitative direction right but never produced
a genuine growth burst for Sand (`R_peak` sat AT the trivial maintenance
floor `m0*B0=0.4`, i.e. pure decay of the seeded biomass, not growth) or a
genuine slow rise for Vertisol (its "rise" was really an early spike that
then dipped, because the shared, fast `D_scale` recruits a soil's whole
habitat almost at once). Stage 10 explicitly sets aside the
shared/ranked-scalar purity of Stages 7-9 and tunes PER-SOIL physical
parameters directly (never arbitrary curve-drawing -- the mechanism stays
the same finite-fuel + Monod-limitation + implicit-transport physics kept
from Stage 6 onward):

  Sand     -- concentrated, fast, finite: a large total, mostly on the fast
              coating, a raised per-soil D_scale so the burst arrives early
              despite Sand's huge single macropore region, and a LONGER
              observation window (the same candidate crashes cleanly by
              t=1500 at n=40, but needs materially longer at n=128 -- a
              grid-scale timescale effect, checked directly below).
  Vertisol -- distributed, slow, sustained: a large total, mostly in the
              slow trapped pool, and its own per-soil D_scale SLOWED well
              below the shared Stage 7-9 value so its ~70 dispersed
              macropore regions recruit gradually across the whole window.

Usage:
    python stage10_run.py
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

SOILS = ["sand", "vertisol"]
LABEL = {"sand": "Sand", "vertisol": "Vertisol"}
COLORS = {"sand": "#d62728", "vertisol": "#4b0082"}
EXPERIMENT_SHAPE = {"sand": "falling", "vertisol": "rising"}

RESULTS_DIR = os.path.join("results", "stage10")
LOGBOOK_PATH = "LOGBOOK.md"

GRID_N = 128
GRID_N_SMALL = 40
MATRIC_D_CUT = 50.0

MAINTENANCE_FLOOR = 0.4   # m0 * B0_total = 0.02 * 20 (configs/biology.yaml, mapping.yaml)
DISTINCTNESS_THRESHOLD = 0.3


# ------------------------------------------------------------- run helpers

def _T_for(soil, mapping_cfg):
    return int(mapping_cfg[f"stage10_T_{soil}"])


def _structure_cfg(soil, grid_n, mapping_cfg):
    cfg = copy.deepcopy(load_yaml(f"configs/soil_{soil}.yaml"))
    cfg["saturation"] = {"mode": "matric_regions", "matric_d_cut_um": MATRIC_D_CUT}
    cfg["grid"]["n"] = grid_n
    cfg["T"] = _T_for(soil, mapping_cfg)
    return cfg


def run_soil(soil, grid_n, biology_cfg, mapping_cfg):
    cfg = _structure_cfg(soil, grid_n, mapping_cfg)
    out = dual_porosity.simulate_dual_porosity(biology_cfg, cfg, mapping_cfg, wicking_enabled=False)
    conn = metrics.connectivity_diagnostics(out["grids"])
    return out, conn


def run_all(grid_n, biology_cfg, mapping_cfg):
    return {soil: run_soil(soil, grid_n, biology_cfg, mapping_cfg) for soil in SOILS}


def normalized_time_resample(R_t, n_points=500):
    """Resample R(t) onto a common [0, 1] fraction-of-window time axis
    (linear interpolation) -- needed because Sand and Vertisol use
    DIFFERENT observation windows (T_sand != T_vertisol), so their raw
    R_t arrays cannot be compared element-wise. Used only for the
    distinctness check; all other diagnostics use each soil's own,
    honestly-reported absolute time axis."""
    x_orig = np.linspace(0.0, 1.0, len(R_t))
    x_new = np.linspace(0.0, 1.0, n_points)
    return np.interp(x_new, x_orig, R_t)


# --------------------------------------------------------- pass/fail (SS3)

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
    growth_share_at_peak = growth_at_peak / (growth_at_peak + maint_at_peak)
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


# --------------------------------------------------------------- 1. headline

def write_headline_csv(data, results_dir, mapping_cfg):
    rows = []
    sand_checks, sand_vals = sand_pass_fail(data["sand"][0])
    vert_checks, vert_vals = vertisol_pass_fail(data["vertisol"][0])
    for soil, checks, vals in (("sand", sand_checks, sand_vals), ("vertisol", vert_checks, vert_vals)):
        out, grids_conn = data[soil]
        grids = out["grids"]
        row = {
            "soil": LABEL[soil],
            "om_total": grids["om_total_init"],
            "coating_fraction": grids["coating_fraction"],
            "D_scale": mapping_cfg[f"stage10_D_scale_{soil}"],
            "T": mapping_cfg[f"stage10_T_{soil}"],
            "R_peak": float(out["R_t"].max()),
            "R_end": float(out["R_t"][-1]),
            "R_end_over_peak": float(out["R_t"][-1] / out["R_t"].max()),
            "peak_frac": float(np.argmax(out["R_t"]) / (len(out["R_t"]) - 1)),
            "B_final_over_B0": float(out["B_final"].sum() / grids["B0"].sum()),
            "B_peak_over_B0": float(out["B_total_t"].max() / grids["B0"].sum()),
            "fuel_remaining_final": float(out["fuel_remaining_t"][-1]),
            "habitats_recruited_final": float(out["habitats_recruited_t"][-1]),
            "n_regions": int(grids["region_diam_um"].size),
            "outcome": metrics.classify_outcome(out),
            "all_checks_pass": all(checks.values()),
        }
        for name, passed in checks.items():
            row[f"check: {name}"] = passed
        rows.append(row)
    path = os.path.join(results_dir, "stage10_headline.csv")
    # Sand and Vertisol have DIFFERENT §3 check names, so each row's key set
    # differs -- fieldnames must be the union (in first-seen order), not
    # just rows[0]'s keys, or DictWriter raises on the second soil's row.
    fieldnames = list(dict.fromkeys(k for row in rows for k in row))
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, restval="")
        w.writeheader()
        w.writerows(rows)
    return rows, (sand_checks, sand_vals), (vert_checks, vert_vals)


def grid_check(biology_cfg, mapping_cfg, results_dir):
    rows = []
    data40 = run_all(GRID_N_SMALL, biology_cfg, mapping_cfg)
    sand_checks, sand_vals = sand_pass_fail(data40["sand"][0])
    vert_checks, vert_vals = vertisol_pass_fail(data40["vertisol"][0])
    for soil, checks, vals in (("sand", sand_checks, sand_vals), ("vertisol", vert_checks, vert_vals)):
        out, _ = data40[soil]
        row = {
            "soil": LABEL[soil], "grid_n": GRID_N_SMALL,
            "R_peak": float(out["R_t"].max()),
            "R_end_over_peak": float(out["R_t"][-1] / out["R_t"].max()),
            "peak_frac": float(np.argmax(out["R_t"]) / (len(out["R_t"]) - 1)),
            "fuel_remaining_final": float(out["fuel_remaining_t"][-1]),
            "outcome": metrics.classify_outcome(out),
            "all_checks_pass": all(checks.values()),
        }
        rows.append(row)
    path = os.path.join(results_dir, "stage10_grid_check.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows, data40


# --------------------------------------------------------------- 2. figures

def plot_respiration_curves(data, results_dir, mapping_cfg):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for soil in SOILS:
        out, _ = data[soil]
        t = np.arange(len(out["R_t"])) * out["dt"]
        axes[0].plot(t, out["R_t"], color=COLORS[soil], label=f"{LABEL[soil]} (T={mapping_cfg[f'stage10_T_{soil}']})", lw=2)
    axes[0].axhline(MAINTENANCE_FLOOR, color="gray", ls=":", lw=1, label=f"maintenance floor ({MAINTENANCE_FLOOR})")
    axes[0].set_xlabel("model time (absolute)"); axes[0].set_ylabel("R(t)")
    axes[0].set_title("Respiration rate, own absolute time axis"); axes[0].legend(fontsize=8)

    for soil in SOILS:
        out, _ = data[soil]
        x = np.linspace(0.0, 1.0, len(out["R_t"]))
        axes[1].plot(x, out["R_t"], color=COLORS[soil], label=LABEL[soil], lw=2)
    axes[1].axhline(MAINTENANCE_FLOOR, color="gray", ls=":", lw=1)
    axes[1].set_xlabel("fraction of own observation window")
    axes[1].set_ylabel("R(t)")
    axes[1].set_title("Respiration rate, normalized time (0=start, 1=end)")
    axes[1].legend(fontsize=8)
    fig.suptitle("Stage 10 (FINAL): Sand real burst-then-crash vs Vertisol slow rise, n=128")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage10_respiration_curves.png"), dpi=150)
    plt.close(fig)


def plot_growth_vs_maintenance(data, results_dir):
    """Sand's growth-vs-maintenance split over time -- the direct evidence
    that the peak is growth-driven, not just decay of the seeded biomass."""
    out = data["sand"][0]
    t = np.arange(len(out["R_t"])) * out["dt"]
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    ax.stackplot(t, out["growth_term_t"], out["maintenance_term_t"],
                 labels=["growth-associated ((1-Y)*growth)", "maintenance (m0*B)"],
                 colors=["#d62728", "#f4a582"])
    ax.plot(t, out["R_t"], color="black", lw=1.2, ls="--", label="R(t) total")
    ax.axhline(MAINTENANCE_FLOOR, color="gray", ls=":", lw=1, label=f"maintenance floor ({MAINTENANCE_FLOOR})")
    ax.set_xlabel("model time"); ax.set_ylabel("respiration")
    ax.set_title("Stage 10: Sand's growth-vs-maintenance split over time")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage10_sand_growth_vs_maintenance.png"), dpi=150)
    plt.close(fig)


def plot_mean_s_habitat_and_fuel(data, results_dir):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for soil in SOILS:
        out, _ = data[soil]
        t = np.arange(len(out["S_habitat_mean_t"])) * out["dt"]
        axes[0].plot(t, out["S_habitat_mean_t"], color=COLORS[soil], label=LABEL[soil], lw=1.8)
        axes[1].plot(t, out["fuel_remaining_t"], color=COLORS[soil], label=LABEL[soil], lw=1.8)
    axes[0].set_xlabel("model time"); axes[0].set_ylabel("mean substrate S in habitat")
    axes[0].set_title("mean_S_habitat(t)"); axes[0].legend()
    axes[1].set_xlabel("model time"); axes[1].set_ylabel("fraction of initial OM undissolved")
    axes[1].set_title("fuel_remaining(t)"); axes[1].legend()
    fig.suptitle("Stage 10: fed-habitat substrate and fuel depletion, own absolute time axes")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage10_mean_s_habitat_and_fuel.png"), dpi=150)
    plt.close(fig)


def plot_biomass(data, results_dir):
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    for soil in SOILS:
        out, _ = data[soil]
        t = np.arange(len(out["B_total_t"])) * out["dt"]
        B0 = out["grids"]["B0"].sum()
        ax.plot(t, out["B_total_t"] / B0, color=COLORS[soil], label=LABEL[soil], lw=1.8)
    ax.axhline(1.0, color="gray", ls=":", lw=1, label="B0 (initial)")
    ax.set_xlabel("model time"); ax.set_ylabel("total biomass / B0")
    ax.set_title("Stage 10: biomass(t), own absolute time axes")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage10_biomass_vs_t.png"), dpi=150)
    plt.close(fig)


# --------------------------------------------------------- 3. distinctness

def write_distinctness_csv(data, results_dir):
    r_sand = normalized_time_resample(data["sand"][0]["R_t"])
    r_vert = normalized_time_resample(data["vertisol"][0]["R_t"])
    d = metrics.pairwise_distance(r_sand, r_vert)
    row = {
        "Sand_vs_Vertisol_distance": d,
        "threshold": DISTINCTNESS_THRESHOLD,
        "distinct": d > DISTINCTNESS_THRESHOLD,
        "note": "computed on time-normalized (fraction-of-own-window) resampled curves, "
                "since Sand and Vertisol use different observation windows T",
    }
    path = os.path.join(results_dir, "stage10_distinctness.csv")
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
    print("STAGE 10 (FINAL) -- Sand real burst-then-crash, Vertisol slow rise")
    print("=" * 78)

    print(f"\n[1] headline: n={GRID_N}, wicking OFF, "
          f"Sand(om={mapping_cfg['stage10_om_total_sand']}, coat={mapping_cfg['stage10_coating_fraction_sand']}, "
          f"D_scale={mapping_cfg['stage10_D_scale_sand']}, T={mapping_cfg['stage10_T_sand']}), "
          f"Vertisol(om={mapping_cfg['stage10_om_total_vertisol']}, coat={mapping_cfg['stage10_coating_fraction_vertisol']}, "
          f"D_scale={mapping_cfg['stage10_D_scale_vertisol']}, T={mapping_cfg['stage10_T_vertisol']}) ...")
    data = run_all(GRID_N, biology_cfg, mapping_cfg)
    rows, (sand_checks, sand_vals), (vert_checks, vert_vals) = write_headline_csv(data, RESULTS_DIR, mapping_cfg)

    print("\n    SAND pass/fail:")
    for name, passed in sand_checks.items():
        print(f"      [{'PASS' if passed else 'FAIL'}] {name}")
    print(f"      values: {sand_vals}")

    print("\n    VERTISOL pass/fail:")
    for name, passed in vert_checks.items():
        print(f"      [{'PASS' if passed else 'FAIL'}] {name}")
    print(f"      values: {vert_vals}")

    plot_respiration_curves(data, RESULTS_DIR, mapping_cfg)
    plot_growth_vs_maintenance(data, RESULTS_DIR)
    plot_mean_s_habitat_and_fuel(data, RESULTS_DIR)
    plot_biomass(data, RESULTS_DIR)

    dist_row = write_distinctness_csv(data, RESULTS_DIR)
    print(f"\n[2] distinctness (time-normalized): distance={dist_row['Sand_vs_Vertisol_distance']:.4f}, "
          f"distinct={dist_row['distinct']}")

    print(f"\n[3] grid check (n={GRID_N_SMALL}) ...")
    grid_rows, data40 = grid_check(biology_cfg, mapping_cfg, RESULTS_DIR)
    for r in grid_rows:
        print(f"    {r['soil']:10} R_peak={r['R_peak']:.3f}  R_end/pk={r['R_end_over_peak']:.3f}  "
              f"peak_frac={r['peak_frac']:.3f}  all_checks_pass={r['all_checks_pass']}")

    sand_all_pass = all(sand_checks.values())
    vert_all_pass = all(vert_checks.values())
    overall = sand_all_pass and vert_all_pass and dist_row["distinct"]

    with open(LOGBOOK_PATH, "a", encoding="utf-8") as f:
        f.write(f"\n**Stage 10 headline (n={GRID_N}): Sand all checks pass={sand_all_pass}, "
                f"Vertisol all checks pass={vert_all_pass}, distinct={dist_row['distinct']} "
                f"(distance={dist_row['Sand_vs_Vertisol_distance']:.4f}). Overall success={overall}. "
                f"{datetime.datetime.now().isoformat(timespec='seconds')}. Full detail in "
                "docs/stage_results/STAGE10_RESULTS.md.**\n")

    print(f"\n[4] OVERALL: Sand pass={sand_all_pass}, Vertisol pass={vert_all_pass}, "
          f"distinct={dist_row['distinct']} -> success={overall}")
    print(f"\n[5] all results written to {RESULTS_DIR}/")


if __name__ == "__main__":
    main()
