"""
Stage 7 runner (prompts/soil_respiration_prompt_stage7_calibration.md).

Stage 6 proved the ported implicit-transport + coating/trapped mechanism
works *structurally* (dozens of macropore regions recruit for Loess/
Vertisol) but left fed-habitat substrate ~25x too low for growth to ever
beat maintenance decay (`r_max=1, m0=0.02, m_s=0.1, Ks=0.2` => growth needs
`f_S > ~0.11`, i.e. habitat `S >~ 0.025`). This is a quantitative
CALIBRATION of three SHARED scalars, not a new mechanism:

  Change 1 -- a single shared `stage7_om_fraction` replaces the per-soil OM
              fractions (`dual_porosity._variantc_om_fraction`).
  Change 2 -- `stage7_D_scale`, a shared post-scale multiplier on the
              transport magnitude (raises sqrt(D*T) reach without touching
              `D_d_ref`, Stage 3.5's soil-discriminating SHAPE of D(d)).
  Change 3 -- `D_film`/`stage7_film_K_scale`: drained (dry) macropore cells
              get a small shared film diffusivity instead of D=0, hosting
              activity at reduced capacity.

SS5's calibration loop (`MAX_ITERATIONS=15`) searches these SHARED scalars
only (never per-soil) at the fast grid (n=40) and scores each candidate on
two axes: whether fed-habitat substrate clears the ~0.025 threshold, and
whether the three curves stay pairwise distinct (>0.3) instead of
collapsing into the degenerate "everything blooms and rises together"
outcome. The winning candidate is then run at the headline grid (n=128)
with the full diagnostic suite + an n=40 check.

Usage:
    python stage7_run.py
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

RESULTS_DIR = os.path.join("results", "stage7")
LOGBOOK_PATH = "LOGBOOK.md"

GRID_N = 128
GRID_N_SMALL = 40
GRID_N_SEARCH = 40   # calibration loop runs at the fast grid

MATRIC_D_CUT = 50.0          # Stage 5's matric potential, unchanged
GROWTH_S_TARGET = 0.025      # SS0: habitat S must clear this for growth to beat decay
DISTINCTNESS_THRESHOLD = 0.3

MAX_ITERATIONS = 15   # SS0: budget for the shared calibration loop

# Fixed (not searched) Change-3 constant -- Change 3's own D_film magnitude
# IS searched; the capacity-reduction factor for a film cell is held fixed
# so the search stays a manageable 3 axes (om_fraction, D_scale, D_film).
STAGE7_FILM_K_SCALE = 0.45


# ------------------------------------------------------------- run helpers

def _structure_cfg(soil, grid_n):
    cfg = copy.deepcopy(load_yaml(f"configs/soil_{soil}.yaml"))
    cfg["saturation"] = {"mode": "matric_regions", "matric_d_cut_um": MATRIC_D_CUT}
    cfg["grid"]["n"] = grid_n
    return cfg


def _mapping_with_candidate(base_mapping, candidate):
    m = dict(base_mapping)
    # Clear any higher-priority om_total tier a LATER stage may have added to
    # the checked-in mapping.yaml since this script was written (Stage 8's
    # stage8_om_total, Stage 9's stage9_om_base) -- otherwise
    # dual_porosity._variantc_om_total would silently ignore this script's
    # own stage7_om_fraction candidate below.
    m.pop("stage8_om_total", None)
    m.pop("stage9_om_base", None)
    m["stage7_om_fraction"] = candidate["om_fraction"]
    m["stage7_D_scale"] = candidate["D_scale"]
    m["D_film"] = candidate["D_film"]
    m["stage7_film_K_scale"] = STAGE7_FILM_K_SCALE
    return m


def run_soil(soil, grid_n, biology_cfg, mapping_cfg):
    cfg = _structure_cfg(soil, grid_n)
    out = dual_porosity.simulate_dual_porosity(biology_cfg, cfg, mapping_cfg, wicking_enabled=False)
    conn = metrics.connectivity_diagnostics(out["grids"])
    return out, conn


def run_all(grid_n, biology_cfg, mapping_cfg):
    return {soil: run_soil(soil, grid_n, biology_cfg, mapping_cfg) for soil in SOILS}


# ------------------------------------------------------- 1. calibration loop

def score_candidate(data):
    """Two axes (SS5): growth cleared (min fed-habitat S over soils, vs
    GROWTH_S_TARGET) and shape distinctness (metrics.emergent_distinctness,
    require_nontrivial=False -- a structurally-caused die-off is still a
    legitimate, distinct outcome per the existing convention). Also flags
    the degenerate "everything blooms and rises together" outcome: all
    three soils `thriving` with B_final > 2x B0 AND all three period-shapes
    "rising" -- distinct in magnitude only, not in the qualitative pattern
    the prompt asks for."""
    outs = {s: data[s][0] for s in SOILS}
    mean_S_hab = {s: metrics.mean_substrate_in_habitat(outs[s]["S_final"], outs[s]["grids"]["K"])
                  for s in SOILS}
    min_S_hab = min(mean_S_hab.values())
    growth_cleared = min_S_hab >= GROWTH_S_TARGET

    success, detail = metrics.emergent_distinctness(
        outs["loess"], outs["sand"], outs["vertisol"], threshold=DISTINCTNESS_THRESHOLD,
        require_nontrivial=False)

    shapes = {s: shape_of_periods(period_means(outs[s]["R_t"])) for s in SOILS}
    b_ratios = {s: float(outs[s]["B_final"].sum() / outs[s]["grids"]["B0"].sum()) for s in SOILS}
    degenerate_bloom = (all(shapes[s] == "rising" for s in SOILS)
                        and all(b_ratios[s] > 2.0 for s in SOILS))

    return {
        "mean_S_hab": mean_S_hab, "min_S_hab": min_S_hab, "growth_cleared": growth_cleared,
        "distinct": success, "min_distance": detail["min_distance"],
        "distances": detail["distances"], "shapes": shapes, "b_ratios": b_ratios,
        "degenerate_bloom": degenerate_bloom,
    }


def append_logbook(iteration, candidate, score, note):
    lines = [f"\n### Stage 7 calibration -- iteration {iteration} -- "
             f"{datetime.datetime.now().isoformat(timespec='seconds')}\n"]
    lines.append(f"- candidate: om_fraction={candidate['om_fraction']:.3f}, "
                 f"D_scale={candidate['D_scale']:.3f}, D_film={candidate['D_film']:.3f}\n")
    lines.append(f"- min fed-habitat S over soils = {score['min_S_hab']:.4f} "
                 f"(target >= {GROWTH_S_TARGET}) -> growth_cleared={score['growth_cleared']}\n")
    lines.append(f"- pairwise distances: "
                 f"{ {k: round(v, 3) for k, v in score['distances'].items()} }, "
                 f"min={score['min_distance']:.3f} -> distinct={score['distinct']}\n")
    lines.append(f"- shapes: {score['shapes']}, B_final/B0: "
                 f"{ {k: round(v, 2) for k, v in score['b_ratios'].items()} }, "
                 f"degenerate_bloom={score['degenerate_bloom']}\n")
    lines.append(f"- note: {note}\n")
    with open(LOGBOOK_PATH, "a", encoding="utf-8") as f:
        f.writelines(lines)


def calibration_loop(biology_cfg, base_mapping):
    """SS5: small budget-capped search over the shared scalars only. Starts
    from a deliberately LOW dose (so the first iterations demonstrate the
    starvation end of the tradeoff) and geometrically raises om_fraction/
    D_scale/D_film together until fed-habitat substrate clears
    GROWTH_S_TARGET for the weakest soil; then, if not yet distinct, tries a
    few D_film-only nudges (the one Change-3 lever that acts differentially
    per soil, since soils differ sharply in how much of their habitat is
    dry: Sand ~100%, Loess/Vertisol ~91-92%) before exhausting the budget.
    Returns (winning_candidate, winning_score, history)."""
    candidate = {"om_fraction": 0.30, "D_scale": 1.5, "D_film": 0.05}
    growth_step = 1.35
    history = []
    best_growth_cleared = None
    best_growth_cleared_score = None

    for iteration in range(1, MAX_ITERATIONS + 1):
        mapping_cfg = _mapping_with_candidate(base_mapping, candidate)
        data = run_all(GRID_N_SEARCH, biology_cfg, mapping_cfg)
        score = score_candidate(data)
        history.append((dict(candidate), score))

        if score["growth_cleared"]:
            if best_growth_cleared is None or score["min_distance"] > best_growth_cleared_score:
                best_growth_cleared, best_growth_cleared_score = dict(candidate), score["min_distance"]

        if score["growth_cleared"] and score["distinct"]:
            note = "growth cleared AND distinct -- full success, stopping."
            append_logbook(iteration, candidate, score, note)
            return candidate, score, history

        if not score["growth_cleared"]:
            note = (f"fed-habitat substrate still below target "
                    f"({score['min_S_hab']:.4f} < {GROWTH_S_TARGET}) -- raising "
                    f"om_fraction/D_scale/D_film together by {growth_step:.2f}x.")
            candidate = {
                "om_fraction": candidate["om_fraction"] * growth_step,
                "D_scale": candidate["D_scale"] * growth_step,
                "D_film": candidate["D_film"] * growth_step,
            }
        else:
            # Growth cleared but curves are not distinct (or degenerately
            # blooming together). Try nudging D_film alone -- the only axis
            # that touches soils asymmetrically (their dry-habitat SHARE
            # differs) -- down on odd iterations, up on even, a small local
            # probe rather than another geometric escalation (SS5: "do not
            # tune per-soil to force shapes; only the shared scalars move").
            direction = -1.0 if iteration % 2 else 1.0
            note = (f"growth cleared (min_S_hab={score['min_S_hab']:.4f}) but not distinct "
                    f"(min_distance={score['min_distance']:.3f}, degenerate_bloom="
                    f"{score['degenerate_bloom']}) -- probing D_film "
                    f"{'down' if direction < 0 else 'up'}.")
            candidate = dict(candidate)
            candidate["D_film"] = max(candidate["D_film"] * (1.0 + direction * 0.25), 1e-4)

        append_logbook(iteration, candidate, score, note)

    # Budget exhausted. SS5's primary target is the growth threshold (SS0);
    # among growth-cleared candidates, keep the one with the best (even if
    # still failing) distinctness rather than the last candidate probed.
    if best_growth_cleared is not None:
        mapping_cfg = _mapping_with_candidate(base_mapping, best_growth_cleared)
        data = run_all(GRID_N_SEARCH, biology_cfg, mapping_cfg)
        final_score = score_candidate(data)
        with open(LOGBOOK_PATH, "a", encoding="utf-8") as f:
            f.write(f"\n**Stage 7 calibration: budget of {MAX_ITERATIONS} iterations exhausted "
                    "without clearing BOTH axes. Falling back to the best growth-cleared "
                    f"candidate found: {best_growth_cleared} "
                    f"(min_distance={final_score['min_distance']:.3f}).**\n")
        return best_growth_cleared, final_score, history

    with open(LOGBOOK_PATH, "a", encoding="utf-8") as f:
        f.write(f"\n**Stage 7 calibration: budget of {MAX_ITERATIONS} iterations exhausted; "
                "fed-habitat substrate never cleared the growth target for all three soils. "
                f"Reporting the final candidate probed: {candidate}.**\n")
    return candidate, history[-1][1], history


# --------------------------------------------------------------- 2. headline

def headline_run(biology_cfg, mapping_cfg):
    return run_all(GRID_N, biology_cfg, mapping_cfg)


def write_headline_csv(data, results_dir):
    rows = []
    for soil in SOILS:
        out, conn = data[soil]
        grids = out["grids"]
        pm = period_means(out["R_t"])
        rows.append({
            "soil": LABEL[soil],
            "om_total_init": grids["om_total_init"],
            "emergent_theta": grids["theta"],
            "wet_habitat_fraction": conn["wet_habitat_fraction"],
            "dry_macro_fraction_of_habitat": grids["dry_macro_fraction_of_habitat"],
            "mean_S_habitat_final": metrics.mean_substrate_in_habitat(out["S_final"], grids["K"]),
            "mean_S_habitat_timeavg": float(out["S_habitat_mean_t"].mean()),
            "fS_habitat_final": metrics.mean_substrate_in_habitat(out["S_final"], grids["K"])
                                / (metrics.mean_substrate_in_habitat(out["S_final"], grids["K"]) + 0.20),
            "growth_cleared": metrics.mean_substrate_in_habitat(out["S_final"], grids["K"]) >= GROWTH_S_TARGET,
            "habitats_recruited_final": float(out["habitats_recruited_t"][-1]),
            "habitats_recruited_max": float(out["habitats_recruited_t"].max()),
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
        })
    path = os.path.join(results_dir, "stage7_headline.csv")
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
            "growth_cleared": metrics.mean_substrate_in_habitat(out["S_final"], grids["K"]) >= GROWTH_S_TARGET,
            "outcome": metrics.classify_outcome(out),
            "period_shape_model": shape_of_periods(pm),
            "period_shape_experiment": shape_of_periods(EXPERIMENT[soil]),
        })
    path = os.path.join(results_dir, "stage7_grid_check.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


# --------------------------------------------------------------- 3. figures

def plot_comparison_periods(data, results_dir, mapping_cfg):
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
    x = np.arange(4)
    for ax, soil in zip(axes, SOILS):
        out, _ = data[soil]
        model = normalized_periods(period_means(out["R_t"]))
        exp = normalized_periods(EXPERIMENT[soil])
        ax.plot(x, exp, "o-", color="black", lw=2, ms=7, label="experiment")
        ax.plot(x, model, "s--", color=COLORS[soil], lw=2, ms=7, label="model (Stage 7)")
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
    fig.suptitle("Stage 7: calibrated shared substrate budget/reach vs measured respiration "
                 f"(om_fraction={mapping_cfg['stage7_om_fraction']:.2f}, "
                 f"D_scale={mapping_cfg['stage7_D_scale']:.2f}, D_film={mapping_cfg['D_film']:.2f}, n={GRID_N})")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage7_comparison_periods.png"), dpi=150)
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
    fig.suptitle(f"Stage 7 -- constant moisture, n={GRID_N}")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage7_respiration_curves.png"), dpi=150)
    plt.close(fig)


def plot_mean_s_habitat_vs_t(data, results_dir):
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    for soil in SOILS:
        out, _ = data[soil]
        t = np.arange(len(out["S_habitat_mean_t"])) * out["dt"]
        ax.plot(t, out["S_habitat_mean_t"], color=COLORS[soil], label=LABEL[soil], lw=1.8)
    ax.axhline(GROWTH_S_TARGET, color="k", ls="--", lw=1, label=f"growth target S={GROWTH_S_TARGET}")
    ax.axhline(0.2, color="gray", ls=":", lw=1, label="Ks=0.2")
    ax.set_xlabel("model time")
    ax.set_ylabel("mean substrate S in habitat cells")
    ax.set_title("Stage 7: fed-habitat substrate over time (the gap this stage closes)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage7_mean_s_habitat_vs_t.png"), dpi=150)
    plt.close(fig)


def plot_habitats_recruited_vs_t(data, results_dir):
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    for soil in SOILS:
        out, _ = data[soil]
        t = np.arange(len(out["habitats_recruited_t"])) * out["dt"]
        ax.plot(t, out["habitats_recruited_t"], color=COLORS[soil], label=LABEL[soil], lw=1.8)
    ax.set_xlabel("model time")
    ax.set_ylabel("count of actively-growing macropore regions")
    ax.set_title("Stage 7: habitats recruited over time")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage7_habitats_recruited_vs_t.png"), dpi=150)
    plt.close(fig)


def plot_calibration_history(history, results_dir):
    its = np.arange(1, len(history) + 1)
    min_s = [s["min_S_hab"] for _, s in history]
    min_d = [s["min_distance"] for _, s in history]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    axes[0].plot(its, min_s, "o-", color="tab:blue")
    axes[0].axhline(GROWTH_S_TARGET, color="k", ls="--", label=f"target={GROWTH_S_TARGET}")
    axes[0].set_xlabel("calibration iteration"); axes[0].set_ylabel("min fed-habitat S over soils")
    axes[0].set_title("Growth-cleared axis"); axes[0].legend()
    axes[1].plot(its, min_d, "o-", color="tab:red")
    axes[1].axhline(DISTINCTNESS_THRESHOLD, color="k", ls="--", label=f"threshold={DISTINCTNESS_THRESHOLD}")
    axes[1].set_xlabel("calibration iteration"); axes[1].set_ylabel("min pairwise distance")
    axes[1].set_title("Distinctness axis"); axes[1].legend()
    fig.suptitle("Stage 7 calibration loop: the two scoring axes across iterations")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage7_calibration_history.png"), dpi=150)
    plt.close(fig)


def write_calibration_history_csv(history, results_dir):
    path = os.path.join(results_dir, "stage7_calibration_history.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["iteration", "om_fraction", "D_scale", "D_film", "min_S_hab",
                    "growth_cleared", "min_distance", "distinct", "degenerate_bloom",
                    "shapes"])
        for i, (cand, score) in enumerate(history, start=1):
            w.writerow([i, cand["om_fraction"], cand["D_scale"], cand["D_film"],
                        score["min_S_hab"], score["growth_cleared"], score["min_distance"],
                        score["distinct"], score["degenerate_bloom"], score["shapes"]])


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
    path = os.path.join(results_dir, "stage7_distinctness.csv")
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
    print("STAGE 7 -- calibration to close the substrate gap (uniform OM %)")
    print("=" * 78)

    if not os.path.exists(LOGBOOK_PATH):
        with open(LOGBOOK_PATH, "w", encoding="utf-8") as f:
            f.write("# LOGBOOK\n")
    with open(LOGBOOK_PATH, "a", encoding="utf-8") as f:
        f.write(f"\n## Stage 7 -- shared calibration search "
                f"(soil_respiration_prompt_stage7_calibration.md, MAX_ITERATIONS={MAX_ITERATIONS})"
                f" -- {datetime.datetime.now().isoformat(timespec='seconds')}\n")

    print(f"\n[1] calibration loop (n={GRID_N_SEARCH}, budget={MAX_ITERATIONS}) ...")
    winning_candidate, winning_score, history = calibration_loop(biology_cfg, base_mapping)
    print(f"    winning candidate: {winning_candidate}")
    print(f"    min fed-habitat S = {winning_score['min_S_hab']:.4f} "
          f"(growth_cleared={winning_score['growth_cleared']}), "
          f"min_distance={winning_score['min_distance']:.3f} (distinct={winning_score['distinct']}), "
          f"degenerate_bloom={winning_score['degenerate_bloom']}")
    write_calibration_history_csv(history, RESULTS_DIR)
    plot_calibration_history(history, RESULTS_DIR)

    # The search's winning candidate becomes the calibration used below --
    # configs/mapping.yaml's checked-in stage7_* values are this same
    # winning candidate, frozen there for every OTHER module to read (never
    # auto-written by this script, matching every prior stage's convention
    # -- LOGBOOK.md's header: "mapping.yaml ... NEVER edited by the search
    # loop").
    mapping_cfg = _mapping_with_candidate(base_mapping, winning_candidate)

    print(f"\n[2] headline: matric_d_cut_um={MATRIC_D_CUT:.0f}, n={GRID_N}, wicking OFF, "
          f"om_fraction={mapping_cfg['stage7_om_fraction']}, D_scale={mapping_cfg['stage7_D_scale']}, "
          f"D_film={mapping_cfg['D_film']} ...")
    data = headline_run(biology_cfg, mapping_cfg)
    rows = write_headline_csv(data, RESULTS_DIR)

    print(f"    {'soil':10} {'om_total':>9} {'S_hab_final':>11} {'growth_ok':>9} {'hab.recruit':>11} "
          f"{'outcome':>10} {'model shp':>10} {'exp shp':>10}")
    for r in rows:
        print(f"    {r['soil']:10} {r['om_total_init']:9.1f} {r['mean_S_habitat_final']:11.5f} "
              f"{str(r['growth_cleared']):>9} {r['habitats_recruited_final']:11.0f} "
              f"{r['outcome']:>10} {r['period_shape_model']:>10} {r['period_shape_experiment']:>10}")

    plot_comparison_periods(data, RESULTS_DIR, mapping_cfg)
    plot_respiration_curves(data, RESULTS_DIR)
    plot_mean_s_habitat_vs_t(data, RESULTS_DIR)
    plot_habitats_recruited_vs_t(data, RESULTS_DIR)

    dist_row = write_distinctness_csv(data, RESULTS_DIR)
    print(f"\n[3] distinctness (require_nontrivial=False): success={dist_row['success']}, "
          f"min_distance={dist_row['min_distance']:.4f}")

    print(f"\n[4] grid check (n={GRID_N_SMALL}) ...")
    grid_rows = grid_check(biology_cfg, mapping_cfg, RESULTS_DIR)
    for r in grid_rows:
        print(f"    {r['soil']:10} S_hab_final={r['mean_S_habitat_final']:.5f}  "
              f"growth_ok={r['growth_cleared']}  outcome={r['outcome']:>10}  "
              f"shape={r['period_shape_model']:>10}")

    # D_film "live habitat gained" report (Change 3)
    print(f"\n[5] D_film live-habitat-opened-up per soil:")
    for soil in SOILS:
        out, _ = data[soil]
        g = out["grids"]
        print(f"    {LABEL[soil]:10} dry macro cells given a film pathway: "
              f"{g['dry_macro_cells']} / {int(g['macro_mask'].sum())} habitat cells "
              f"({g['dry_macro_fraction_of_habitat']*100:.1f}% of habitat)")

    print(f"\n[6] all results written to {RESULTS_DIR}/")


if __name__ == "__main__":
    main()
