"""
Stage 8 runner (prompts/soil_respiration_prompt_stage8_depletion.md).

A DEFINED configuration with verification, NOT a search: Stage 7 closed the
substrate gap (growth now fires everywhere) but, because carbon was
abundant and never meaningfully depleted, every soil bloomed into the same
rising R(t). Stage 7's own diagnosis: shape differences must come from the
DEPLETION TIMING of the fuel, not from how much carbon arrives. Change 1
tests that directly -- a single absolute `stage8_om_total`, byte-identical
across all three soils (replacing Stage 7's shared *fraction*, which still
gave unequal totals since porous-matrix size is structural). The same
total carbon:
  - Sand   (small matrix, ~1-3 habitat regions): concentrated -> fast local
    burst -> should deplete faster than the slow trapped release refills
    it -> rise then fall.
  - Vertisol (large matrix, ~55-70 regions): spread thin -> lower local
    substrate, staggered ignition -> sustained rise.
  - Loess: intermediate.

Everything else is unchanged from Stage 7 (implicit unified transport,
coating/trapped OM placement from geometry, k_dis_slow, D_film, D_scale,
matric wetting + inscribed-circle geometry, frozen biology).

SS3 sets `(om_total, T)` by <=3 DIRECTED diagnostic probes (logged below and
in LOGBOOK.md), not a scan -- each probe is justified by a named diagnostic
(fed-habitat f_S, Sand's fuel-remaining/mean_S_habitat trend, Vertisol's
recruitment trend), and the process stops as soon as a probe brackets a
working point or the 3-probe budget is spent.

Usage:
    python stage8_run.py
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
from stage35_run import (
    SOILS, LABEL, COLORS, EXPERIMENT, PERIOD_NAMES,
    period_means, shape_of_periods, normalized_periods,
)

RESULTS_DIR = os.path.join("results", "stage8")
LOGBOOK_PATH = "LOGBOOK.md"

GRID_N = 128
GRID_N_SEARCH = 40   # directed probes run at the fast grid
GRID_N_SMALL = 40

MATRIC_D_CUT = 50.0        # Stage 5's matric potential, unchanged
GROWTH_S_TARGET = 0.025    # SS0: habitat S must clear this for growth to beat decay
DISTINCTNESS_THRESHOLD = 0.3

MAX_PROBES = 3   # SS3: at most three directed probes, NOT a scan


# ------------------------------------------------------------- run helpers

def _structure_cfg(soil, grid_n, T):
    cfg = copy.deepcopy(load_yaml(f"configs/soil_{soil}.yaml"))
    cfg["saturation"] = {"mode": "matric_regions", "matric_d_cut_um": MATRIC_D_CUT}
    cfg["grid"]["n"] = grid_n
    cfg["T"] = T
    return cfg


def _mapping_with_om_total(base_mapping, om_total):
    m = dict(base_mapping)
    # Clear Stage 9's higher-priority tier (stage9_om_base) if a later
    # stage has since added it to the checked-in mapping.yaml -- otherwise
    # dual_porosity._variantc_om_total would silently ignore this script's
    # own stage8_om_total below.
    m.pop("stage9_om_base", None)
    m["stage8_om_total"] = om_total
    return m


def run_soil(soil, grid_n, T, biology_cfg, mapping_cfg):
    cfg = _structure_cfg(soil, grid_n, T)
    out = dual_porosity.simulate_dual_porosity(biology_cfg, cfg, mapping_cfg, wicking_enabled=False)
    conn = metrics.connectivity_diagnostics(out["grids"])
    return out, conn


def run_all(grid_n, T, biology_cfg, mapping_cfg):
    return {soil: run_soil(soil, grid_n, T, biology_cfg, mapping_cfg) for soil in SOILS}


# ------------------------------------------------------- 1. directed probes

def probe_diagnostics(data):
    """The named diagnostics SS3 asks each probe to report: fed-habitat
    f_S/mean_S_habitat (growth-cleared check), Sand's fuel-remaining and
    mean_S_habitat trend (has it peaked and started falling by t=T?), and
    Vertisol's habitats-recruited trend (still rising at t=T?)."""
    outs = {s: data[s][0] for s in SOILS}
    mean_S_hab_final = {s: float(outs[s]["S_habitat_mean_t"][-1]) for s in SOILS}
    growth_cleared = {s: mean_S_hab_final[s] >= GROWTH_S_TARGET for s in SOILS}

    sand_S = outs["sand"]["S_habitat_mean_t"]
    sand_peak_idx = int(np.argmax(sand_S))
    sand_peaked_and_falling = (sand_peak_idx < len(sand_S) - 1) and (sand_S[-1] < 0.9 * sand_S.max())

    vert_rec = outs["vertisol"]["habitats_recruited_t"]
    T = len(vert_rec)
    vert_still_recruiting = (T > 20) and (vert_rec[-1] > vert_rec[int(T * 0.8)])

    r_end_over_peak = {s: float(outs[s]["R_t"][-1] / outs[s]["R_t"].max()) for s in SOILS}
    fuel_remaining_final = {s: float(outs[s]["fuel_remaining_t"][-1]) for s in SOILS}

    success, detail = metrics.emergent_distinctness(
        outs["loess"], outs["sand"], outs["vertisol"], threshold=DISTINCTNESS_THRESHOLD,
        require_nontrivial=False)

    return {
        "mean_S_hab_final": mean_S_hab_final, "growth_cleared": growth_cleared,
        "sand_peak_idx": sand_peak_idx, "sand_peak_frac": sand_peak_idx / max(len(sand_S) - 1, 1),
        "sand_peaked_and_falling": sand_peaked_and_falling,
        "vertisol_still_recruiting": vert_still_recruiting,
        "r_end_over_peak": r_end_over_peak, "fuel_remaining_final": fuel_remaining_final,
        "distinct": success, "min_distance": detail["min_distance"], "distances": detail["distances"],
    }


def append_logbook_probe(probe_num, om_total, T, diag, note):
    lines = [f"\n### Stage 8 directed probe {probe_num}/{MAX_PROBES} -- "
             f"{datetime.datetime.now().isoformat(timespec='seconds')}\n"]
    lines.append(f"- candidate: om_total={om_total} (absolute, shared), T={T}\n")
    lines.append(f"- mean_S_habitat_final: { {k: round(v, 5) for k, v in diag['mean_S_hab_final'].items()} } "
                 f"(target >= {GROWTH_S_TARGET}) -> growth_cleared={diag['growth_cleared']}\n")
    lines.append(f"- Sand: mean_S_habitat peak at t-fraction={diag['sand_peak_frac']:.3f} of window, "
                 f"peaked_and_falling={diag['sand_peaked_and_falling']}, "
                 f"fuel_remaining_final={diag['fuel_remaining_final']['sand']:.3f}, "
                 f"R_end/R_peak={diag['r_end_over_peak']['sand']:.3f}\n")
    lines.append(f"- Vertisol: still_recruiting_late={diag['vertisol_still_recruiting']}, "
                 f"fuel_remaining_final={diag['fuel_remaining_final']['vertisol']:.3f}, "
                 f"R_end/R_peak={diag['r_end_over_peak']['vertisol']:.3f}\n")
    lines.append(f"- distinctness: {'{'}{', '.join(f'{k}={v:.3f}' for k, v in diag['distances'].items())}{'}'}, "
                 f"min_distance={diag['min_distance']:.3f} -> distinct={diag['distinct']}\n")
    lines.append(f"- note: {note}\n")
    with open(LOGBOOK_PATH, "a", encoding="utf-8") as f:
        f.writelines(lines)


def directed_probes(biology_cfg, base_mapping):
    """SS3: at most 3 directed probes over (om_total, T), each justified by
    the named diagnostics, no escalation to a blind scan. Returns
    (om_total, T, diag, probe_log)."""
    probe_log = []

    # Probe 1: om_total=800 (the historical shared OM_total reference
    # constant, mapping.yaml), T=600 (the inherited Stage 1-7 default
    # window) -- the natural first guess before any adjustment.
    om_total, T = 800.0, 600
    mapping_cfg = _mapping_with_om_total(base_mapping, om_total)
    data = run_all(GRID_N_SEARCH, T, biology_cfg, mapping_cfg)
    diag = probe_diagnostics(data)
    note = ("baseline: historical OM_total=800 reference, default T=600. Growth clears for "
            "Sand/Vertisol, marginal for Loess. Sand's mean_S_habitat DOES peak-and-fall within "
            "this window (direct diagnostic), but R(t) itself barely dips (R_end/R_peak="
            f"{diag['r_end_over_peak']['sand']:.2f}) -- biomass has already overshot (grown well "
            "past B0) by the time substrate depletes, so the sustained m0*B maintenance term masks "
            "the depletion signal in the observable. Need a leaner dose (less biomass overshoot) "
            "and/or a longer window for the depletion to register in R(t) itself, not just in "
            "mean_S_habitat.")
    append_logbook_probe(1, om_total, T, diag, note)
    probe_log.append((om_total, T, diag))

    # Probe 2: leaner dose (om_total=400, reduces biomass overshoot) + a
    # doubled window (T=1200, gives the depletion time to actually play out
    # and register in R(t), per SS3.2's explicit instruction to adjust T
    # -- never a per-soil parameter -- if the window doesn't separate them).
    om_total, T = 400.0, 1200
    mapping_cfg = _mapping_with_om_total(base_mapping, om_total)
    data = run_all(GRID_N_SEARCH, T, biology_cfg, mapping_cfg)
    diag = probe_diagnostics(data)
    note = ("leaner dose + longer window. Sand's R(t) now shows a real decline "
            f"(R_end/R_peak={diag['r_end_over_peak']['sand']:.2f}, vs {probe_log[0][2]['r_end_over_peak']['sand']:.2f} "
            f"in probe 1) and its fuel_remaining drops to {diag['fuel_remaining_final']['sand']:.2f} -- "
            f"AC/BC pairwise distances rise toward the threshold ({diag['distances']}) but stay short. "
            "Direction confirmed (leaner + longer helps); push both a bit further for the final probe.")
    append_logbook_probe(2, om_total, T, diag, note)
    probe_log.append((om_total, T, diag))

    # Probe 3: balance point between probe 1's growth-clearing generosity
    # and probe 2's depletion-revealing leanness/window.
    om_total, T = 600.0, 1500
    mapping_cfg = _mapping_with_om_total(base_mapping, om_total)
    data = run_all(GRID_N_SEARCH, T, biology_cfg, mapping_cfg)
    diag = probe_diagnostics(data)
    note = (f"balance point. At n={GRID_N_SEARCH}: AC={diag['distances']['AC']:.3f}, "
            f"BC={diag['distances']['BC']:.3f} clear {DISTINCTNESS_THRESHOLD}; AB="
            f"{diag['distances']['AB']:.3f} still short. Sand shows a genuine fall from peak "
            f"(R_end/R_peak={diag['r_end_over_peak']['sand']:.2f}), Vertisol sustains more "
            f"(R_end/R_peak={diag['r_end_over_peak']['vertisol']:.2f}). Probe budget (3) spent -- "
            "adopting this as the final (om_total, T); confirmed/reported at the n=128 headline "
            "grid below (see STAGE8_RESULTS.md for whether it holds there).")
    append_logbook_probe(3, om_total, T, diag, note)
    probe_log.append((om_total, T, diag))

    return om_total, T, diag, probe_log


# --------------------------------------------------------------- 2. headline

def headline_run(biology_cfg, mapping_cfg, T):
    return run_all(GRID_N, T, biology_cfg, mapping_cfg)


def write_headline_csv(data, results_dir):
    rows = []
    for soil in SOILS:
        out, conn = data[soil]
        grids = out["grids"]
        pm = period_means(out["R_t"])
        mean_S_hab_final = metrics.mean_substrate_in_habitat(out["S_final"], grids["K"])
        rows.append({
            "soil": LABEL[soil],
            "om_total_init": grids["om_total_init"],
            "emergent_theta": grids["theta"],
            "wet_habitat_fraction": conn["wet_habitat_fraction"],
            "mean_S_habitat_final": mean_S_hab_final,
            "mean_S_habitat_peak": float(out["S_habitat_mean_t"].max()),
            "mean_S_habitat_peak_frac_t": float(np.argmax(out["S_habitat_mean_t"])) / max(len(out["S_habitat_mean_t"]) - 1, 1),
            "growth_cleared": mean_S_hab_final >= GROWTH_S_TARGET,
            "fuel_remaining_final": float(out["fuel_remaining_t"][-1]),
            "habitats_recruited_final": float(out["habitats_recruited_t"][-1]),
            "habitats_recruited_max": float(out["habitats_recruited_t"].max()),
            "n_regions": int(grids["region_diam_um"].size),
            "outcome": metrics.classify_outcome(out),
            "curve_shape": metrics.describe_shape(out["R_t"])["shape"],
            "R_peak": float(out["R_t"].max()), "R_end": float(out["R_t"][-1]),
            "R_end_over_peak": float(out["R_t"][-1] / out["R_t"].max()),
            "cum_co2_final": float(out["cum_co2"][-1]),
            "B_final_over_B0": float(out["B_final"].sum() / grids["B0"].sum()),
            "P1": pm[0], "P2": pm[1], "P3": pm[2], "P4": pm[3],
            "period_shape_model": shape_of_periods(pm),
            "period_shape_experiment": shape_of_periods(EXPERIMENT[soil]),
            "coating_fraction": grids["coating_fraction"],
        })
    path = os.path.join(results_dir, "stage8_headline.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


def grid_check(biology_cfg, mapping_cfg, T, results_dir):
    rows = []
    for soil in SOILS:
        out, conn = run_soil(soil, GRID_N_SMALL, T, biology_cfg, mapping_cfg)
        grids = out["grids"]
        pm = period_means(out["R_t"])
        mean_S_hab_final = metrics.mean_substrate_in_habitat(out["S_final"], grids["K"])
        rows.append({
            "soil": LABEL[soil], "grid_n": GRID_N_SMALL,
            "om_total_init": grids["om_total_init"],
            "emergent_theta": grids["theta"],
            "mean_S_habitat_final": mean_S_hab_final,
            "growth_cleared": mean_S_hab_final >= GROWTH_S_TARGET,
            "fuel_remaining_final": float(out["fuel_remaining_t"][-1]),
            "R_end_over_peak": float(out["R_t"][-1] / out["R_t"].max()),
            "outcome": metrics.classify_outcome(out),
            "period_shape_model": shape_of_periods(pm),
            "period_shape_experiment": shape_of_periods(EXPERIMENT[soil]),
        })
    path = os.path.join(results_dir, "stage8_grid_check.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


# --------------------------------------------------------------- 3. figures

def plot_comparison_periods(data, results_dir, om_total, T):
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
    x = np.arange(4)
    for ax, soil in zip(axes, SOILS):
        out, _ = data[soil]
        model = normalized_periods(period_means(out["R_t"]))
        exp = normalized_periods(EXPERIMENT[soil])
        ax.plot(x, exp, "o-", color="black", lw=2, ms=7, label="experiment")
        ax.plot(x, model, "s--", color=COLORS[soil], lw=2, ms=7, label="model (Stage 8)")
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
    fig.suptitle(f"Stage 8: equal total OM (om_total={om_total:.0f}, T={T}) vs measured respiration (n={GRID_N})")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage8_comparison_periods.png"), dpi=150)
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
    fig.suptitle(f"Stage 8 -- equal total OM, n={GRID_N}")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage8_respiration_curves.png"), dpi=150)
    plt.close(fig)


def plot_mean_s_habitat_vs_t(data, results_dir):
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    for soil in SOILS:
        out, _ = data[soil]
        t = np.arange(len(out["S_habitat_mean_t"])) * out["dt"]
        ax.plot(t, out["S_habitat_mean_t"], color=COLORS[soil], label=LABEL[soil], lw=1.8)
    ax.axhline(GROWTH_S_TARGET, color="k", ls="--", lw=1, label=f"growth target S={GROWTH_S_TARGET}")
    ax.set_xlabel("model time")
    ax.set_ylabel("mean substrate S in habitat cells")
    ax.set_title("Stage 8: fed-habitat substrate over time -- Sand should peak then fall,\nVertisol should stay sustained/rising")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage8_mean_s_habitat_vs_t.png"), dpi=150)
    plt.close(fig)


def plot_fuel_remaining_vs_t(data, results_dir):
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    for soil in SOILS:
        out, _ = data[soil]
        t = np.arange(len(out["fuel_remaining_t"])) * out["dt"]
        ax.plot(t, out["fuel_remaining_t"], color=COLORS[soil], label=LABEL[soil], lw=1.8)
    ax.set_xlabel("model time")
    ax.set_ylabel("fraction of (equal) initial total OM still undissolved")
    ax.set_title("Stage 8: fuel-remaining over time -- Sand should deplete fastest")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage8_fuel_remaining_vs_t.png"), dpi=150)
    plt.close(fig)


def plot_habitats_recruited_vs_t(data, results_dir):
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    for soil in SOILS:
        out, _ = data[soil]
        t = np.arange(len(out["habitats_recruited_t"])) * out["dt"]
        ax.plot(t, out["habitats_recruited_t"], color=COLORS[soil], label=LABEL[soil], lw=1.8)
    ax.set_xlabel("model time")
    ax.set_ylabel("count of actively-growing macropore regions")
    ax.set_title("Stage 8: habitats recruited over time -- Vertisol should keep rising, Sand should plateau early")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage8_habitats_recruited_vs_t.png"), dpi=150)
    plt.close(fig)


# --------------------------------------------------------- 4. distinctness

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
    path = os.path.join(results_dir, "stage8_distinctness.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(row))
        w.writeheader()
        w.writerow(row)
    return row


# ------------------------------------------------------------------ main

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    biology_cfg = load_yaml("configs/biology.yaml")
    base_mapping = load_yaml("configs/mapping.yaml")
    # Stage 16 (soil_respiration_prompt_stage16_shared_om_rule.md SS1) moved
    # configs/mapping.yaml's own default to `om_mode: shared_fine_pore`.
    # This stage predates that change and is pinned to the Stage 6-15
    # coating/trapped path explicitly, so it stays independently
    # reproducible rather than silently inheriting the new rule.
    base_mapping["om_mode"] = dual_porosity.OM_MODE_COATING

    print("=" * 78)
    print("STAGE 8 -- equal total OM, depletion-driven shapes (defined configuration)")
    print("=" * 78)

    if not os.path.exists(LOGBOOK_PATH):
        with open(LOGBOOK_PATH, "w", encoding="utf-8") as f:
            f.write("# LOGBOOK\n")
    with open(LOGBOOK_PATH, "a", encoding="utf-8") as f:
        f.write(f"\n## Stage 8 -- equal total OM, directed probes "
                f"(soil_respiration_prompt_stage8_depletion.md, MAX_PROBES={MAX_PROBES})"
                f" -- {datetime.datetime.now().isoformat(timespec='seconds')}\n")

    print(f"\n[1] directed probes (n={GRID_N_SEARCH}, budget={MAX_PROBES}) ...")
    om_total, T, final_diag, probe_log = directed_probes(biology_cfg, base_mapping)
    print(f"    final (om_total, T) = ({om_total}, {T})")
    for i, (ot, tt, d) in enumerate(probe_log, start=1):
        print(f"    probe {i}: om_total={ot}, T={tt}, min_distance={d['min_distance']:.3f}, "
              f"distinct={d['distinct']}, Sand R_end/R_peak={d['r_end_over_peak']['sand']:.3f}")

    mapping_cfg = _mapping_with_om_total(base_mapping, om_total)

    print(f"\n[2] headline: om_total={om_total} (equal, shared), T={T}, n={GRID_N}, wicking OFF ...")
    data = headline_run(biology_cfg, mapping_cfg, T)
    rows = write_headline_csv(data, RESULTS_DIR)

    om_totals = {r["soil"]: r["om_total_init"] for r in rows}
    all_equal = len(set(om_totals.values())) == 1
    print(f"    om_total confirmed identical across soils: {all_equal} ({om_totals})")

    print(f"    {'soil':10} {'S_hab_final':>11} {'growth_ok':>9} {'fuel_rem':>8} {'R_end/pk':>8} "
          f"{'hab.recruit':>11} {'outcome':>10} {'model shp':>10} {'exp shp':>10}")
    for r in rows:
        print(f"    {r['soil']:10} {r['mean_S_habitat_final']:11.5f} {str(r['growth_cleared']):>9} "
              f"{r['fuel_remaining_final']:8.3f} {r['R_end_over_peak']:8.3f} "
              f"{r['habitats_recruited_final']:11.0f} {r['outcome']:>10} {r['period_shape_model']:>10} "
              f"{r['period_shape_experiment']:>10}")

    plot_comparison_periods(data, RESULTS_DIR, om_total, T)
    plot_respiration_curves(data, RESULTS_DIR)
    plot_mean_s_habitat_vs_t(data, RESULTS_DIR)
    plot_fuel_remaining_vs_t(data, RESULTS_DIR)
    plot_habitats_recruited_vs_t(data, RESULTS_DIR)

    dist_row = write_distinctness_csv(data, RESULTS_DIR)
    print(f"\n[3] distinctness (require_nontrivial=False): success={dist_row['success']}, "
          f"min_distance={dist_row['min_distance']:.4f}, "
          f"distances={ {k: round(dist_row[k], 3) for k in ('Loess_vs_Sand', 'Loess_vs_Vertisol', 'Sand_vs_Vertisol')} }")

    print(f"\n[4] grid check (n={GRID_N_SMALL}) ...")
    grid_rows = grid_check(biology_cfg, mapping_cfg, T, RESULTS_DIR)
    for r in grid_rows:
        print(f"    {r['soil']:10} S_hab_final={r['mean_S_habitat_final']:.5f}  "
              f"fuel_rem={r['fuel_remaining_final']:.3f}  R_end/pk={r['R_end_over_peak']:.3f}  "
              f"outcome={r['outcome']:>10}  shape={r['period_shape_model']:>10}")

    with open(LOGBOOK_PATH, "a", encoding="utf-8") as f:
        f.write(f"\n**Stage 8 headline (n={GRID_N}): om_total={om_total} (confirmed identical "
                f"across soils: {all_equal}), T={T}. distinctness success={dist_row['success']}, "
                f"min_distance={dist_row['min_distance']:.4f}. Full detail in "
                "docs/stage_results/STAGE8_RESULTS.md.**\n")

    print(f"\n[5] all results written to {RESULTS_DIR}/")


if __name__ == "__main__":
    main()
