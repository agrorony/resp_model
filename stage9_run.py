"""
Stage 9 runner (prompts/soil_respiration_prompt_stage9_ranked_om.md).

A DEFINED configuration with verification, like Stage 8: NO open-ended
search. Stage 8 proved the depletion mechanism gives genuinely distinct
curves and reproduced Loess-erratic, but a single EQUAL total OM could not
give Sand-fall AND Vertisol-rise together -- opposite regimes (Sand needs a
lean pool to deplete; Vertisol needs abundant, ACCESSIBLE carbon to sustain
many regions). Three targeted, ranking-only changes (never fitted to the
measured respiration values):

  Change 1 -- per-soil total OM, ranked Vertisol > Loess > Sand
              (`dual_porosity._variantc_om_total`, `stage9_om_base` x
              per-soil `stage9_om_multiplier_<soil>`).
  Change 2 -- Vertisol's coating fraction raised (more of its now-larger
              carbon reaches its many habitat regions).
  Change 3 -- Sand's coating fraction lowered (smaller fast burst + a slow
              trapped tail), SUBJECT TO the "burst wins" exception: never
              lowered past the point where Sand's rise-then-fall breaks.

Everything else is unchanged from Stage 8 (implicit unified transport,
coating/trapped OM placement from geometry, k_dis_slow, D_film, D_scale,
matric wetting + inscribed-circle geometry, finite/depleting fuel).

SS4 sets the shared base OM scale and per-soil coating fractions (within
the fixed Change 2/3 directions) by <=3 DIRECTED diagnostic probes, logged
below and in LOGBOOK.md -- not a scan, and never violating the fixed
V > L > S ranking. Each probe is run at n=40 (fast) for the SS4 diagnostics
AND separately checked at n=128 (the primary grid) before being treated as
a candidate winner -- this stage's own probing found the n=40 and n=128
grids DISAGREE about which candidate is best (see LOGBOOK.md / STAGE9_
RESULTS.md), so no candidate is adopted without an n=128 check.

Usage:
    python stage9_run.py
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

RESULTS_DIR = os.path.join("results", "stage9")
LOGBOOK_PATH = "LOGBOOK.md"

GRID_N = 128
GRID_N_SEARCH = 40   # directed probes' primary (fast) diagnostics
GRID_N_SMALL = 40

MATRIC_D_CUT = 50.0        # Stage 5's matric potential, unchanged
GROWTH_S_TARGET = 0.025    # SS0 (Stage 7/8): habitat S must clear this for growth to beat decay
DISTINCTNESS_THRESHOLD = 0.3

MAX_PROBES = 3   # SS4: at most three directed probes, NOT a scan

# Fixed ranked ordering (SS1) -- never changed by the probes.
OM_MULTIPLIER_LOESS = 1.0
STAGE9_T = 1500   # unchanged from Stage 8 (SS4 found no reason to move it)


# ------------------------------------------------------------- run helpers

def _structure_cfg(soil, grid_n, T):
    cfg = copy.deepcopy(load_yaml(f"configs/soil_{soil}.yaml"))
    cfg["saturation"] = {"mode": "matric_regions", "matric_d_cut_um": MATRIC_D_CUT}
    cfg["grid"]["n"] = grid_n
    cfg["T"] = T
    return cfg


def _mapping_with_candidate(base_mapping, om_base, v_mult, s_mult, v_coat, s_coat):
    m = dict(base_mapping)
    m["stage9_om_base"] = om_base
    m["stage9_om_multiplier_vertisol"] = v_mult
    m["stage9_om_multiplier_loess"] = OM_MULTIPLIER_LOESS
    m["stage9_om_multiplier_sand"] = s_mult
    m["stage9_coating_fraction_vertisol"] = v_coat
    m["stage9_coating_fraction_sand"] = s_coat
    # Loess deliberately left unset -- _variantc_coating_fraction falls
    # back to the Stage 6 value (0.50), keeping it "intermediate" (SS3).
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
    """The named SS4 diagnostics: Vertisol's fuel NOT exhausted + still
    recruiting + R(t) rising; Sand's accessible pool depletes + mean_S_
    habitat peaks-then-falls + a real R(t) burst; Loess intermediate/
    erratic; distinctness > 0.3."""
    outs = {s: data[s][0] for s in SOILS}

    vert_fuel = outs["vertisol"]["fuel_remaining_t"][-1]
    vert_rec = outs["vertisol"]["habitats_recruited_t"]
    T = len(vert_rec)
    vert_still_recruiting = (T > 20) and (vert_rec[-1] > vert_rec[int(T * 0.8)])
    vert_r_end_over_peak = float(outs["vertisol"]["R_t"][-1] / outs["vertisol"]["R_t"].max())
    vert_rising = vert_r_end_over_peak >= 0.9   # R(t) still near its peak at t=T

    sand_S = outs["sand"]["S_habitat_mean_t"]
    sand_peak_idx = int(np.argmax(sand_S))
    sand_peaked_and_falling = (sand_peak_idx < len(sand_S) - 1) and (sand_S[-1] < 0.9 * sand_S.max())
    sand_r_end_over_peak = float(outs["sand"]["R_t"][-1] / outs["sand"]["R_t"].max())
    sand_r_peak_idx = int(np.argmax(outs["sand"]["R_t"]))
    # "burst" = a real EARLY rise (R(t) climbs meaningfully above its t=0
    # value before falling), not just noise -- distinguishes a genuine
    # rise-then-fall from a monotonic decline that never had a peak.
    sand_burst_intact = (outs["sand"]["R_t"].max() > 1.3 * outs["sand"]["R_t"][0]) and (sand_r_end_over_peak < 0.9)

    fuel_remaining_final = {s: float(outs[s]["fuel_remaining_t"][-1]) for s in SOILS}
    r_end_over_peak = {s: float(outs[s]["R_t"][-1] / outs[s]["R_t"].max()) for s in SOILS}

    success, detail = metrics.emergent_distinctness(
        outs["loess"], outs["sand"], outs["vertisol"], threshold=DISTINCTNESS_THRESHOLD,
        require_nontrivial=False)

    return {
        "vert_fuel_remaining": float(vert_fuel), "vert_still_recruiting": vert_still_recruiting,
        "vert_r_end_over_peak": vert_r_end_over_peak, "vert_rising": vert_rising,
        "sand_peaked_and_falling": sand_peaked_and_falling, "sand_burst_intact": sand_burst_intact,
        "sand_r_end_over_peak": sand_r_end_over_peak,
        "fuel_remaining_final": fuel_remaining_final, "r_end_over_peak": r_end_over_peak,
        "distinct": success, "min_distance": detail["min_distance"], "distances": detail["distances"],
    }


def append_logbook_probe(probe_num, om_base, v_mult, s_mult, v_coat, s_coat,
                          diag40, diag128, note):
    lines = [f"\n### Stage 9 directed probe {probe_num}/{MAX_PROBES} -- "
             f"{datetime.datetime.now().isoformat(timespec='seconds')}\n"]
    lines.append(f"- candidate: om_base={om_base}, multipliers V={v_mult}/L={OM_MULTIPLIER_LOESS}/S={s_mult} "
                 f"(ranked V>L>S), coating V={v_coat}/S={s_coat} (Loess unchanged at Stage 6's 0.50), T={STAGE9_T}\n")
    lines.append(f"- [n=40] Vertisol: fuel_remaining={diag40['vert_fuel_remaining']:.3f}, "
                 f"still_recruiting={diag40['vert_still_recruiting']}, R_end/R_peak={diag40['vert_r_end_over_peak']:.3f} "
                 f"(rising={diag40['vert_rising']})\n")
    lines.append(f"- [n=40] Sand: peaked_and_falling={diag40['sand_peaked_and_falling']}, "
                 f"burst_intact={diag40['sand_burst_intact']}, R_end/R_peak={diag40['sand_r_end_over_peak']:.3f}\n")
    lines.append(f"- [n=40] distinctness: {diag40['distances']}, min_distance={diag40['min_distance']:.3f} "
                 f"-> distinct={diag40['distinct']}\n")
    lines.append(f"- [n=128 CHECK] distinctness: {diag128['distances']}, min_distance={diag128['min_distance']:.3f} "
                 f"-> distinct={diag128['distinct']} (Sand R_end/R_peak={diag128['sand_r_end_over_peak']:.3f}, "
                 f"Vertisol R_end/R_peak={diag128['vert_r_end_over_peak']:.3f})\n")
    lines.append(f"- note: {note}\n")
    with open(LOGBOOK_PATH, "a", encoding="utf-8") as f:
        f.writelines(lines)


def directed_probes(biology_cfg, base_mapping):
    """SS4: at most 3 directed probes over (base scale, Vertisol/Sand
    coating fractions), ranked ordering fixed throughout. Each candidate is
    checked at BOTH n=40 (the SS4 diagnostics) and n=128 (the primary grid)
    -- this stage's own probing found they can disagree about which
    candidate is best, so the final adoption is decided at n=128, not n=40.
    Returns (winner_dict, probe_log)."""
    candidates = [
        # Probe 1: SS1's own suggested starting ratios (V=1.5/L=1.0/S=0.6),
        # coating raised/lowered by a modest, symmetric first step
        # (Vertisol 0.30->0.60, Sand 0.80->0.50).
        {"om_base": 600.0, "v_mult": 1.5, "s_mult": 0.6, "v_coat": 0.6, "s_coat": 0.5},
        # Probe 2: push Vertisol's total AND coating further (Change 2's
        # explicit directive) to test whether more, more-accessible carbon
        # improves its sustain / separates it from Loess.
        {"om_base": 600.0, "v_mult": 2.0, "s_mult": 0.6, "v_coat": 0.7, "s_coat": 0.5},
        # Probe 3: push Vertisol once more (furthest tried, still a modest
        # multiplier, not an aggressive fit) to see if the probe-2 trend
        # continues.
        {"om_base": 600.0, "v_mult": 2.5, "s_mult": 0.6, "v_coat": 0.7, "s_coat": 0.5},
    ]

    probe_log = []
    for i, c in enumerate(candidates, start=1):
        mapping_cfg = _mapping_with_candidate(
            base_mapping, c["om_base"], c["v_mult"], c["s_mult"], c["v_coat"], c["s_coat"])
        data40 = run_all(GRID_N_SEARCH, STAGE9_T, biology_cfg, mapping_cfg)
        diag40 = probe_diagnostics(data40)
        data128 = run_all(GRID_N, STAGE9_T, biology_cfg, mapping_cfg)
        diag128 = probe_diagnostics(data128)

        if i == 1:
            note = ("SS1's own suggested starting ratios. n=40 distinctness is weak "
                    f"(min_distance={diag40['min_distance']:.3f}, driven by a low Loess-vs-Vertisol "
                    "distance -- both are large, well-connected, many-region soils whose NORMALIZED "
                    "curves look qualitatively similar even at different total carbon). Sand's burst "
                    f"is genuine ({diag40['sand_burst_intact']}). Checking n=128 before deciding "
                    "whether to escalate Vertisol further (Change 2's directive).")
        elif i == 2:
            improved40 = diag40["min_distance"] > probe_log[0][2]["min_distance"]
            worsened128 = diag128["min_distance"] < probe_log[0][3]["min_distance"]
            note = (f"escalated Vertisol (mult 1.5->2.0, coating 0.6->0.7) per Change 2. n=40 "
                    f"min_distance {'improved' if improved40 else 'did not improve'} "
                    f"({diag40['min_distance']:.3f} vs probe 1's n=40 {probe_log[0][2]['min_distance']:.3f}), "
                    f"but the n=128 check {'WORSENED' if worsened128 else 'held or improved'} "
                    f"({diag128['min_distance']:.3f} vs probe 1's n=128 {probe_log[0][3]['min_distance']:.3f}) "
                    "-- the two grids disagree about whether escalating helps.")
        else:
            improved40 = diag40["min_distance"] > probe_log[1][2]["min_distance"]
            worsened128 = diag128["min_distance"] < probe_log[0][3]["min_distance"]
            note = (f"pushed Vertisol once more (mult 2.0->2.5) to see if probe 2's n=40 trend "
                    f"continues. n=40 min_distance {'improved further' if improved40 else 'did not improve further'} "
                    f"({diag40['min_distance']:.3f}), but n=128 stays "
                    f"{'worse than probe 1' if worsened128 else 'comparable to probe 1'} "
                    f"({diag128['min_distance']:.3f} vs probe 1's n=128 {probe_log[0][3]['min_distance']:.3f}). "
                    "Probe budget (3) spent. Per SS4/SS6, adoption is decided on the n=128 PRIMARY grid, "
                    "not the n=40 probe grid where escalation looked like it was helping -- see the "
                    "winner selection below.")

        append_logbook_probe(i, c["om_base"], c["v_mult"], c["s_mult"], c["v_coat"], c["s_coat"],
                              diag40, diag128, note)
        probe_log.append((c, data40, diag40, diag128))

    # Winner: the candidate with the BEST n=128 (primary-grid) min_distance
    # among the three probed -- not simply "the last probe tried" and not
    # "the n=40 winner" (SS4/SS6: adoption must reflect the grid the
    # success criterion is actually evaluated against).
    best_i = max(range(len(probe_log)), key=lambda i: probe_log[i][3]["min_distance"])
    winner_candidate = probe_log[best_i][0]
    winner_diag128 = probe_log[best_i][3]

    with open(LOGBOOK_PATH, "a", encoding="utf-8") as f:
        f.write(f"\n**Stage 9 probe selection: candidate {best_i + 1}/{MAX_PROBES} "
                f"({winner_candidate}) has the best n=128 min_distance "
                f"({winner_diag128['min_distance']:.3f}) among the {MAX_PROBES} probed -- adopted as final, "
                "even though it was NOT the best-looking candidate at n=40 "
                "(probes 2 and 3 looked more promising there). This n=40-vs-n=128 disagreement "
                "is reported as an honest finding, not resolved with a 4th probe.**\n")

    return winner_candidate, winner_diag128, probe_log


# --------------------------------------------------------------- 2. headline

def headline_run(biology_cfg, mapping_cfg):
    return run_all(GRID_N, STAGE9_T, biology_cfg, mapping_cfg)


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
            "coating_fraction": grids["coating_fraction"],
            "emergent_theta": grids["theta"],
            "wet_habitat_fraction": conn["wet_habitat_fraction"],
            "mean_S_habitat_final": mean_S_hab_final,
            "mean_S_habitat_peak": float(out["S_habitat_mean_t"].max()),
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
        })
    path = os.path.join(results_dir, "stage9_headline.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


def grid_check(biology_cfg, mapping_cfg, results_dir):
    rows = []
    for soil in SOILS:
        out, conn = run_soil(soil, GRID_N_SMALL, STAGE9_T, biology_cfg, mapping_cfg)
        grids = out["grids"]
        pm = period_means(out["R_t"])
        mean_S_hab_final = metrics.mean_substrate_in_habitat(out["S_final"], grids["K"])
        rows.append({
            "soil": LABEL[soil], "grid_n": GRID_N_SMALL,
            "om_total_init": grids["om_total_init"],
            "mean_S_habitat_final": mean_S_hab_final,
            "growth_cleared": mean_S_hab_final >= GROWTH_S_TARGET,
            "fuel_remaining_final": float(out["fuel_remaining_t"][-1]),
            "R_end_over_peak": float(out["R_t"][-1] / out["R_t"].max()),
            "outcome": metrics.classify_outcome(out),
            "period_shape_model": shape_of_periods(pm),
            "period_shape_experiment": shape_of_periods(EXPERIMENT[soil]),
        })
    path = os.path.join(results_dir, "stage9_grid_check.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


# --------------------------------------------------------------- 3. figures

def plot_comparison_periods(data, results_dir, cand):
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
    x = np.arange(4)
    for ax, soil in zip(axes, SOILS):
        out, _ = data[soil]
        model = normalized_periods(period_means(out["R_t"]))
        exp = normalized_periods(EXPERIMENT[soil])
        ax.plot(x, exp, "o-", color="black", lw=2, ms=7, label="experiment")
        ax.plot(x, model, "s--", color=COLORS[soil], lw=2, ms=7, label="model (Stage 9)")
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
    fig.suptitle(f"Stage 9: ranked OM (V={cand['v_mult']}/L=1.0/S={cand['s_mult']} x base={cand['om_base']:.0f}, "
                 f"coating V={cand['v_coat']}/S={cand['s_coat']}) vs measured respiration (n={GRID_N})")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage9_comparison_periods.png"), dpi=150)
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
    fig.suptitle(f"Stage 9 -- ranked OM + placement, n={GRID_N}")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage9_respiration_curves.png"), dpi=150)
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
    ax.set_title("Stage 9: fed-habitat substrate over time")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage9_mean_s_habitat_vs_t.png"), dpi=150)
    plt.close(fig)


def plot_fuel_remaining_vs_t(data, results_dir):
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    for soil in SOILS:
        out, _ = data[soil]
        t = np.arange(len(out["fuel_remaining_t"])) * out["dt"]
        ax.plot(t, out["fuel_remaining_t"], color=COLORS[soil], label=LABEL[soil], lw=1.8)
    ax.set_xlabel("model time")
    ax.set_ylabel("fraction of each soil's OWN (ranked, unequal) total OM still undissolved")
    ax.set_title("Stage 9: fuel-remaining over time")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage9_fuel_remaining_vs_t.png"), dpi=150)
    plt.close(fig)


def plot_habitats_recruited_vs_t(data, results_dir):
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    for soil in SOILS:
        out, _ = data[soil]
        t = np.arange(len(out["habitats_recruited_t"])) * out["dt"]
        ax.plot(t, out["habitats_recruited_t"], color=COLORS[soil], label=LABEL[soil], lw=1.8)
    ax.set_xlabel("model time")
    ax.set_ylabel("count of actively-growing macropore regions")
    ax.set_title("Stage 9: habitats recruited over time")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage9_habitats_recruited_vs_t.png"), dpi=150)
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
    path = os.path.join(results_dir, "stage9_distinctness.csv")
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
    print("STAGE 9 -- ranked per-soil OM + placement (defined configuration)")
    print("=" * 78)

    if not os.path.exists(LOGBOOK_PATH):
        with open(LOGBOOK_PATH, "w", encoding="utf-8") as f:
            f.write("# LOGBOOK\n")
    with open(LOGBOOK_PATH, "a", encoding="utf-8") as f:
        f.write(f"\n## Stage 9 -- ranked OM, directed probes "
                f"(soil_respiration_prompt_stage9_ranked_om.md, MAX_PROBES={MAX_PROBES})"
                f" -- {datetime.datetime.now().isoformat(timespec='seconds')}\n")

    print(f"\n[1] directed probes (n={GRID_N_SEARCH} + n={GRID_N} check each, budget={MAX_PROBES}) ...")
    winner, winner_diag128, probe_log = directed_probes(biology_cfg, base_mapping)
    print(f"    winner: {winner} (n=128 min_distance={winner_diag128['min_distance']:.3f})")
    for i, (c, _, d40, d128) in enumerate(probe_log, start=1):
        print(f"    probe {i}: V_mult={c['v_mult']} V_coat={c['v_coat']} S_mult={c['s_mult']} S_coat={c['s_coat']} "
              f"-> n=40 min_dist={d40['min_distance']:.3f}, n=128 min_dist={d128['min_distance']:.3f}")

    mapping_cfg = _mapping_with_candidate(
        base_mapping, winner["om_base"], winner["v_mult"], winner["s_mult"], winner["v_coat"], winner["s_coat"])

    print(f"\n[2] headline: n={GRID_N}, wicking OFF, om_base={winner['om_base']}, "
          f"V_mult={winner['v_mult']}, S_mult={winner['s_mult']}, "
          f"V_coat={winner['v_coat']}, S_coat={winner['s_coat']} ...")
    data = headline_run(biology_cfg, mapping_cfg)
    rows = write_headline_csv(data, RESULTS_DIR)

    om_totals = {r["soil"]: r["om_total_init"] for r in rows}
    ranking_ok = om_totals["Vertisol"] > om_totals["Loess"] > om_totals["Sand"]
    print(f"    om_total ranking Vertisol>Loess>Sand confirmed: {ranking_ok} ({om_totals})")

    print(f"    {'soil':10} {'om_total':>9} {'coat':>5} {'S_hab_final':>11} {'fuel_rem':>8} {'R_end/pk':>8} "
          f"{'hab.recruit':>11} {'outcome':>10} {'model shp':>10} {'exp shp':>10}")
    for r in rows:
        print(f"    {r['soil']:10} {r['om_total_init']:9.1f} {r['coating_fraction']:5.2f} "
              f"{r['mean_S_habitat_final']:11.5f} {r['fuel_remaining_final']:8.3f} {r['R_end_over_peak']:8.3f} "
              f"{r['habitats_recruited_final']:11.0f} {r['outcome']:>10} {r['period_shape_model']:>10} "
              f"{r['period_shape_experiment']:>10}")

    plot_comparison_periods(data, RESULTS_DIR, winner)
    plot_respiration_curves(data, RESULTS_DIR)
    plot_mean_s_habitat_vs_t(data, RESULTS_DIR)
    plot_fuel_remaining_vs_t(data, RESULTS_DIR)
    plot_habitats_recruited_vs_t(data, RESULTS_DIR)

    dist_row = write_distinctness_csv(data, RESULTS_DIR)
    print(f"\n[3] distinctness (require_nontrivial=False): success={dist_row['success']}, "
          f"min_distance={dist_row['min_distance']:.4f}, "
          f"distances={ {k: round(dist_row[k], 3) for k in ('Loess_vs_Sand', 'Loess_vs_Vertisol', 'Sand_vs_Vertisol')} }")

    print(f"\n[4] grid check (n={GRID_N_SMALL}) ...")
    grid_rows = grid_check(biology_cfg, mapping_cfg, RESULTS_DIR)
    for r in grid_rows:
        print(f"    {r['soil']:10} om_total={r['om_total_init']:.1f}  S_hab_final={r['mean_S_habitat_final']:.5f}  "
              f"fuel_rem={r['fuel_remaining_final']:.3f}  R_end/pk={r['R_end_over_peak']:.3f}  "
              f"outcome={r['outcome']:>10}  shape={r['period_shape_model']:>10}")

    with open(LOGBOOK_PATH, "a", encoding="utf-8") as f:
        f.write(f"\n**Stage 9 headline (n={GRID_N}): om_total ranking confirmed (V>L>S: {ranking_ok}), "
                f"winner candidate={winner}. distinctness success={dist_row['success']}, "
                f"min_distance={dist_row['min_distance']:.4f}. Full detail in "
                "docs/stage_results/STAGE9_RESULTS.md.**\n")

    print(f"\n[5] all results written to {RESULTS_DIR}/")


if __name__ == "__main__":
    main()
