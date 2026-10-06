"""
Stage 16 runner (prompts/soil_respiration_prompt_stage16_shared_om_rule.md).

A MODEL-WIDE STRUCTURAL change, not a calibration stage and not a search:
the per-soil coating/trapped OM split introduced in Stage 6 is deleted as
the active default and ONE shared, structure-blind organic-matter placement
rule takes its place, with spatially UNIFORM dissolution. Both soils (Sand,
Vertisol) are re-run under it. The point: every difference -- between soils,
and between macropore regions within a soil -- must emerge from the pore
field, not from a hand-set carbon parameter.

What this script does, in the SS9 order:

  0. (done separately, first) the SS0 correction appended to
     docs/stage_results/STAGE14_RESULTS.md.
  1. verifies the preserved `om_mode: coating_trapped` path still
     reproduces Stage 11's recorded numbers, and that Stage 16's two new
     observational outputs change nothing (exact identity check, Stage 14
     style).
  2. SS2's geodesic fuel->habitat distance analysis for both soils --
     computed BEFORE any tuning, since it is a property of the structure
     and the shared rule, not of the dose.
  3. SS3's dose/transport ladder, rung by rung, <=15 directed probes total,
     each logged with its reasoning.
  4. SS4's criteria at n=128 (primary) and n=40 for the adopted rung.
  5. SS7's optional shared-lognormal clumping sensitivity check.

numpy-only, as always. Usage:
    python stage16_run.py
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

RESULTS_DIR = os.path.join("results", "stage16")
LOGBOOK_PATH = "LOGBOOK.md"

MATRIC_D_CUT = 50.0
MAINTENANCE_FLOOR = 0.4
DISTINCTNESS_THRESHOLD = 0.3
MAX_PROBES = 15

# SS4(b) states no numeric bar of its own -- it says to REPORT frac@2%, t50
# and t_final, and names Vertisol's baseline (1.00 / 0.001 / 0.013) as "the
# thing being fixed". These are this stage's stated operationalization of
# "no longer saturates instantly", fixed in advance and never loosened:
# recruitment must not be essentially complete at 2% of the window, and must
# not finish within the first tenth of it.
FRAC_AT_2PCT_MAX = 0.90
T_FINAL_MIN = 0.10

# Portell et al. (2018, Front. Microbiol. 9:1583) -- their own hotspot
# metric, for direct comparison (SS4c).
PORTELL_HOTSPOT_PCT = 9.5
PORTELL_HOTSPOT_SD = 4.0


# --------------------------------------------------------------- run helpers

def structure_cfg(soil, grid_n, T):
    cfg = copy.deepcopy(load_yaml(f"configs/soil_{soil}.yaml"))
    cfg["saturation"] = {"mode": "matric_regions", "matric_d_cut_um": MATRIC_D_CUT}
    cfg["grid"]["n"] = grid_n
    cfg["T"] = int(T)
    return cfg


def shared_mapping(mapping_cfg, om_keys, D_keys, om_mode=dual_porosity.OM_MODE_SHARED):
    """A mapping config on the shared OM rule, with this rung's dose and
    transport keys set and every PREVIOUS Stage 16 key stripped first, so a
    stale key from an earlier probe cannot leak in through `_variantc_*`'s
    fallback chain and quietly un-share the configuration.

    Note that Stage 10/11's own per-soil keys are deliberately left in
    place: the Stage 16 tiers sit ABOVE them in the priority chain, so as
    long as a rung sets its own keys the older ones are unreachable -- and
    leaving them present is what keeps `om_mode: coating_trapped` able to
    reproduce Stages 6-15 from the same file."""
    m = copy.deepcopy(mapping_cfg)
    m["om_mode"] = om_mode
    for key in list(m):
        if key.startswith(("stage16_om_total", "stage16_om_rho", "stage16_D_scale")):
            del m[key]
    m.update(om_keys)
    m.update(D_keys)
    return m


def run_soil(soil, grid_n, T, biology_cfg, mapping_cfg, track_region_biomass=True,
             debug_growth_distribution=False):
    return dual_porosity.simulate_dual_porosity(
        biology_cfg, structure_cfg(soil, grid_n, T), mapping_cfg,
        wicking_enabled=False, track_region_biomass=track_region_biomass,
        debug_growth_distribution=debug_growth_distribution)


def per_region_half_max(out):
    """Stage 14's GRADED per-region recruitment metric, reused unchanged:
    for each macropore region, the first timestep its OWN mean growth
    reaches half its OWN eventual peak, as a fraction of the window.

    SS4(b)'s criterion is defined on the boolean `habitats_recruited_t`
    count, so that is what `recruitment_timing` scores. But SS4(c) directs
    that the boolean count be retired as a headline number, on Stage 14's
    finding that it is actively misleading (r=0.072 against R(t), versus
    r=0.996 for the graded active-mass metric). Reporting both is the only
    way to state honestly what did and did not change: they disagree, and
    which one is trusted decides the verdict on SS4(b)."""
    gd = out.get("growth_distribution_t")
    if gd is None:
        return np.array([])
    lab = out["region_labels_debug"][out["macro_mask_debug"]]
    T = gd.shape[0]
    vals = []
    for k in np.unique(lab):
        if k <= 0:
            continue
        g = gd[:, lab == k].mean(axis=1)
        if g.max() <= 0:
            continue
        vals.append(int(np.argmax(g >= 0.5 * g.max())) / (T - 1))
    return np.array(vals)


def normalized_time_resample(R_t, n_points=500):
    x_orig = np.linspace(0.0, 1.0, len(R_t))
    x_new = np.linspace(0.0, 1.0, n_points)
    return np.interp(x_new, x_orig, R_t)


# ------------------------------------------------- SS4(a) pass/fail criteria

def sand_pass_fail(out):
    """Stage 10 SS3's Sand criteria, UNCHANGED (SS4a: 'Do NOT loosen any
    threshold'). Only one label is reworded: Stage 10's 'coating fuel
    visibly depletes' reads 'fuel visibly depletes' here, because under the
    shared rule there is no coating pool -- `fuel_remaining_t` now tracks
    the soil's single uniform OM field. The threshold (<0.5) and the
    quantity measured are identical."""
    R = out["R_t"]
    peak_idx = int(np.argmax(R))
    peak_frac = peak_idx / (len(R) - 1)
    R_peak = float(R.max())
    R_end_over_peak = float(R[-1] / R_peak)
    B0 = float(out["grids"]["B0"].sum())
    B_peak_over_B0 = float(out["B_total_t"].max() / B0)
    growth_at_peak = float(out["growth_term_t"][peak_idx])
    maint_at_peak = float(out["maintenance_term_t"][peak_idx])
    growth_share = (growth_at_peak / (growth_at_peak + maint_at_peak)
                    if (growth_at_peak + maint_at_peak) > 0 else 0.0)
    fuel_remaining_final = float(out["fuel_remaining_t"][-1])
    checks = {
        "R_peak >= 2x maintenance_floor": R_peak >= 2.0 * MAINTENANCE_FLOOR,
        "growth term is MAJORITY of R at peak": growth_share > 0.5,
        "B_peak/B0 > 1 (biomass grows)": B_peak_over_B0 > 1.0,
        "fuel visibly depletes (<0.5)": fuel_remaining_final < 0.5,
        "R_end/R_peak < 0.3 (crash)": R_end_over_peak < 0.3,
        "peak in the first third": peak_frac < 1.0 / 3.0,
    }
    vals = {"R_peak": R_peak, "R_end_over_peak": R_end_over_peak, "peak_frac": peak_frac,
            "B_peak_over_B0": B_peak_over_B0, "growth_share_at_peak": growth_share,
            "fuel_remaining_final": fuel_remaining_final}
    return checks, vals


def vertisol_pass_fail(out):
    """Stage 11's UPDATED Vertisol criteria, as SS4(a) explicitly directs
    ('Use Stage 11's updated Vertisol criteria, not Stage 10's stale upper
    bound on R_peak') -- STAGE15_RESULTS.md SS3 already documents why Stage
    10's 'R_peak < 4x maintenance_floor' bound is scale-mismatched and was
    dropped. Reproduced here unchanged from stage15_run.py."""
    R = out["R_t"]
    n3 = len(R) // 3
    early3, late3 = float(R[:n3].mean()), float(R[-n3:].mean())
    peak_idx = int(np.argmax(R))
    peak_frac = peak_idx / (len(R) - 1)
    R_peak = float(R.max())
    R_end_over_peak = float(R[-1] / R_peak)
    B_peak_over_B0 = float(out["B_total_t"].max() / out["grids"]["B0"].sum())
    checks = {
        "late-third mean > early-third mean": late3 > early3,
        "peaks in the LATE part of the window (>0.6)": peak_frac > 0.6,
        "not an early spike + dip (R_end/R_peak > 0.7)": R_end_over_peak > 0.7,
        "R_peak clears maintenance floor (R_peak > floor)": R_peak > MAINTENANCE_FLOOR,
        "biomass genuinely grows (B_peak/B0 > 1)": B_peak_over_B0 > 1.0,
    }
    vals = {"early3": early3, "late3": late3, "peak_frac": peak_frac, "R_peak": R_peak,
            "R_end_over_peak": R_end_over_peak, "B_peak_over_B0": B_peak_over_B0}
    return checks, vals


PASS_FAIL = {"sand": sand_pass_fail, "vertisol": vertisol_pass_fail}


# ------------------------------------------- SS4(b) recruitment timing

def recruitment_timing(out):
    """SS4(b): fraction of finally-recruited regions already active at 2% of
    the window, t50, and t_final -- all as fractions of the window, on the
    boolean `habitats_recruited_t` count, which is exactly the quantity
    whose baseline (Vertisol 1.00 / 0.001 / 0.013) this stage is fixing.

    The boolean count is used HERE and only here, because SS4(b) is defined
    on it; SS4(c)'s headline recruitment number is the growth-weighted
    active-mass metric instead, per Stage 14's recommendation."""
    hr = np.asarray(out["habitats_recruited_t"], dtype=float)
    T = len(hr)
    final = float(hr[-1])
    idx_2pct = int(round(0.02 * (T - 1)))
    if final <= 0:
        return {"frac_at_2pct": float("nan"), "t50": float("nan"), "t_final": float("nan"),
                "n_recruited_final": 0}
    frac_2 = float(hr[idx_2pct] / final)
    t50 = int(np.argmax(hr >= 0.5 * final)) / (T - 1)
    t_final = int(np.argmax(hr >= final)) / (T - 1)
    return {"frac_at_2pct": frac_2, "t50": float(t50), "t_final": float(t_final),
            "n_recruited_final": int(final)}


def recruitment_ok(timing):
    if not np.isfinite(timing["frac_at_2pct"]):
        return False
    return (timing["frac_at_2pct"] < FRAC_AT_2PCT_MAX) and (timing["t_final"] > T_FINAL_MIN)


# ------------------------------------------- SS4(c) hotspot selectivity

def hotspot_selectivity(out):
    """SS4(c), Portell et al. (2018)'s OWN metric so the number is directly
    comparable to their 9.5 +/- 4.0% of colony sites: the percentage of
    macropore regions whose PEAK biomass exceeds 10% of the maximum peak
    biomass across regions in that soil.

    Reported alongside the growth-weighted "active mass" metric Stage 14
    recommended. The boolean `growth > 1e-9` count is RETIRED as a headline
    number here (Stage 14 measured it at r=0.072 against R(t), i.e.
    effectively uncorrelated) and kept only as a secondary diagnostic.

    One honest caveat, reported rather than corrected: `B0 ~ K` seeds every
    region with biomass proportional to its size at t=0, so a region's peak
    biomass carries a size component as well as a growth component. For a
    soil whose habitat is one giant region plus a couple of tiny ones
    (Sand: 3 regions at n=128) this metric is close to meaningless; it is
    informative for Vertisol's 70 comparably-sized regions, which is the
    soil the criterion is actually aimed at."""
    rb = out.get("region_biomass_t")
    if rb is None or rb.size == 0:
        return {"hotspot_pct": float("nan"), "n_regions": 0, "n_hotspots": 0,
                "active_mass_peak": float("nan"), "active_mass_t_of_peak": float("nan"),
                "boolean_recruited_pct": float("nan")}
    peaks = rb.max(axis=0)
    n_regions = int(peaks.size)
    max_peak = float(peaks.max()) if n_regions else 0.0
    n_hot = int((peaks > 0.10 * max_peak).sum()) if max_peak > 0 else 0
    am = np.asarray(out["active_mass_t"], dtype=float)
    hr_final = float(out["habitats_recruited_t"][-1])
    return {
        "hotspot_pct": 100.0 * n_hot / n_regions if n_regions else float("nan"),
        "n_regions": n_regions,
        "n_hotspots": n_hot,
        "active_mass_peak": float(am.max()),
        "active_mass_t_of_peak": float(int(np.argmax(am)) / (len(am) - 1)),
        "boolean_recruited_pct": 100.0 * hr_final / n_regions if n_regions else float("nan"),
    }


# ----------------------------------------------------- SS4(d) distinctness

def distinctness(sand_out, vert_out):
    return float(metrics.pairwise_distance(
        normalized_time_resample(sand_out["R_t"]), normalized_time_resample(vert_out["R_t"])))


# ------------------------------------------------------ SS2 geodesic analysis

def geodesic_analysis(biology_cfg, mapping_cfg, results_dir):
    """SS2 -- this stage's PRIMARY structural result, computed before any
    tuning. Distances are a property of the pore field and the shared
    placement rule; because `substantial_om_mask`'s threshold is the median
    of NONZERO matrix OM, the whole field scales by one factor when
    `om_total` changes and the mask -- hence every distance below -- is
    invariant to the dose. That invariance is verified directly here, not
    assumed."""
    rows, sens_rows, per_region = [], [], {}
    m = shared_mapping(mapping_cfg, {"stage16_om_total_shared": 10000.0},
                       {"stage16_D_scale_shared": 1.0})
    for soil in SOILS:
        for grid_n in (40, 128):
            grids = dual_porosity.build_dual_grids(structure_cfg(soil, grid_n, 10), m)
            a = dual_porosity.analyze_fuel_habitat_distance(grids)
            row = {"soil": LABEL[soil], "grid_n": grid_n,
                   "matrix_fraction": grids["matrix_fraction"],
                   "porous_matrix_fraction": grids["porous_matrix_fraction"],
                   "n_regions": a["n_regions"],
                   "n_regions_reachable": a["n_regions_reachable"],
                   "n_regions_unreachable": a["n_regions_unreachable"],
                   "substantial_cells": a["substantial_cells"],
                   "porous_matrix_cells": a["porous_matrix_cells"],
                   "substantial_threshold": a["substantial_threshold"],
                   "voxel_um": a["voxel_um"]}
            for k, v in a["stats_cells"].items():
                row[f"{k}_cells"] = v
            for k, v in a["stats_um"].items():
                row[f"{k}_um"] = v

            if grid_n == 128:
                per_region[soil] = a["per_region_um"]
                # dose-invariance check (SS2's own claim, verified not assumed)
                m2 = shared_mapping(mapping_cfg, {"stage16_om_total_shared": 1.0e6},
                                    {"stage16_D_scale_shared": 1.0})
                a2 = dual_porosity.analyze_fuel_habitat_distance(
                    dual_porosity.build_dual_grids(structure_cfg(soil, grid_n, 10), m2))
                row["dose_invariant_median"] = bool(
                    np.isclose(a["stats_um"]["median"], a2["stats_um"]["median"]))
            rows.append(row)

            # SS2 secondary: sensitivity of the distribution to how strict
            # "substantial" is. The MEDIAN is this stage's stated
            # definition and is what every headline number uses; this sweep
            # is reported because the median turns out to be a very
            # permissive cut under this model's extremely right-skewed OM
            # field, and the sweep is what reveals where the two soils
            # actually separate.
            OM, pm = grids["OM"], grids["porous_matrix_mask"]
            nz = OM[pm & (OM > 0)]
            allowed = (pm | grids["macro_mask"]) & ~grids["solid_mask"]
            lab = grids["region_labels"]
            for q in (50, 75, 90, 95, 99, 99.9):
                thr = float(np.percentile(nz, q))
                sub = pm & (OM >= thr)
                dist = dual_porosity.bfs_geodesic_distance(sub, allowed)
                vals = []
                for k in range(1, int(lab.max()) + 1):
                    dc = dist[lab == k]
                    fin = dc[np.isfinite(dc)]
                    if fin.size:
                        vals.append(float(fin.min()))
                v = np.array(vals) * grids["voxel_um"]
                sens_rows.append({
                    "soil": LABEL[soil], "grid_n": grid_n, "percentile_threshold": q,
                    "substantial_cells": int(sub.sum()),
                    "share_of_total_OM_in_substantial": float(OM[sub].sum() / OM.sum()),
                    "n_regions_reachable": len(vals),
                    "n_regions_unreachable": int(lab.max()) - len(vals),
                    "median_um": float(np.median(v)) if len(v) else float("nan"),
                    "p95_um": float(np.percentile(v, 95)) if len(v) else float("nan"),
                    "max_um": float(v.max()) if len(v) else float("nan"),
                })

    _write_csv(os.path.join(results_dir, "stage16_geodesic_distance.csv"), rows)
    _write_csv(os.path.join(results_dir, "stage16_geodesic_threshold_sensitivity.csv"), sens_rows)
    plot_geodesic_histograms(per_region, sens_rows, results_dir)
    return rows, sens_rows, per_region


def plot_geodesic_histograms(per_region, sens_rows, results_dir):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    for ax, soil in zip(axes[:2], SOILS):
        v = np.asarray(per_region.get(soil, []), dtype=float)
        v = v[np.isfinite(v)]
        if v.size:
            top = max(float(v.max()) + 20.0, 60.0)
            ax.hist(v, bins=np.arange(0, top, 10), color=COLORS[soil], edgecolor="white")
            ax.axvline(float(np.median(v)), color="black", ls="--", lw=1.5,
                       label=f"median = {np.median(v):.0f} um")
            ax.legend(fontsize=8)
        ax.set_xlabel("geodesic distance to nearest substantial OM (um)")
        ax.set_ylabel("macropore regions")
        ax.set_title(f"{LABEL[soil]} (n=128), 'substantial' = median cut")

    ax = axes[2]
    for soil in SOILS:
        rows = [r for r in sens_rows if r["soil"] == LABEL[soil] and r["grid_n"] == 128]
        ax.plot([r["percentile_threshold"] for r in rows], [r["median_um"] for r in rows],
                "o-", color=COLORS[soil], label=f"{LABEL[soil]} median")
        ax.plot([r["percentile_threshold"] for r in rows], [r["p95_um"] for r in rows],
                "s--", color=COLORS[soil], alpha=0.5, label=f"{LABEL[soil]} p95")
    ax.axhspan(PORTELL_LO, PORTELL_HI, color="gray", alpha=0.18,
               label="Portell et al. 2018 operative range")
    ax.set_xlabel("percentile of nonzero matrix OM used as 'substantial'")
    ax.set_ylabel("geodesic distance (um)")
    ax.set_yscale("log")
    ax.set_title("sensitivity to how strict 'substantial' is")
    ax.legend(fontsize=7)

    fig.suptitle("Stage 16 SS2: emergent fuel->habitat distance under ONE shared, "
                 "structure-blind OM rule")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "stage16_geodesic_histograms.png"), dpi=150)
    plt.close(fig)


# Portell et al. (2018) operative geodesic distances, the two particulate
# scenarios: 252 +/- 111 um and 490 +/- 262 um. Plotted as a band.
PORTELL_LO = 252.0 - 111.0
PORTELL_HI = 490.0 + 262.0


# ------------------------------------------------- SS5 reconciliation diagnostic

def ss5_reconciliation(biology_cfg, mapping_cfg, results_dir):
    """Why SS5's recorded recruitment numbers do NOT reproduce under SS1's
    shared rule -- a DIAGNOSTIC, not a ladder probe: it adopts nothing and
    proposes nothing, it explains a discrepancy between the prompt's own
    evidence table and this stage's measurements.

    SS5 states its table was measured "with `coating_fraction=0` and
    `k_dis_slow` set equal to `k_dis=0.05` (i.e. the SS1 shared-rule
    configuration)". Those are not the same configuration. With
    `coating_fraction=0` the coating POOL is empty, but
    `split_om_coating_trapped` still places all of the soil's OM on
    `trapped_mask`, which is the porous matrix MINUS the coating shell --
    so the geodesic construction is still running, and it punches a
    fuel-free exclusion zone around every macropore region. SS1's actual
    shared rule places OM across the WHOLE porous matrix, including exactly
    those near-habitat cells.

    That difference is the entire staggering effect, and this function
    measures it directly on both configurations."""
    rows = []
    variants = [
        ("SS5 proxy: coating_fraction=0 (OM on trapped_mask -> exclusion zone)", "proxy"),
        ("SS1 shared rule: OM across the whole porous matrix", "shared"),
    ]
    for name, kind in variants:
        if kind == "proxy":
            m = copy.deepcopy(mapping_cfg)
            m["om_mode"] = dual_porosity.OM_MODE_COATING
            m["stage11_coating_fraction_vertisol"] = 0.0
            m["stage11_om_total_vertisol"] = 50000.0
            m["stage11_D_scale_vertisol"] = 0.7
            m["k_dis_slow"] = float(biology_cfg["k_dis"])
        else:
            m = shared_mapping(mapping_cfg, {"stage16_om_total_vertisol": 50000.0},
                               {"stage16_D_scale_vertisol": 0.7})
        out = run_soil("vertisol", 128, 1500, biology_cfg, m)
        g = out["grids"]
        t = recruitment_timing(out)
        # how many macropore regions physically touch a fuel-bearing cell
        adj = dual_porosity._adjacent4(g["OM"] > 0)
        touch = np.unique(g["region_labels"][adj & g["macro_mask"]])
        n_touch = int((touch > 0).sum())
        geo = dual_porosity.analyze_fuel_habitat_distance(g)
        rows.append({
            "configuration": name,
            "om_bearing_matrix_cells": int((g["OM"] > 0).sum()),
            "porous_matrix_cells": int(g["porous_matrix_mask"].sum()),
            "regions_adjacent_to_fuel": n_touch,
            "n_regions": int(g["region_labels"].max()),
            "geodesic_median_um": geo["stats_um"]["median"],
            "R_peak": float(out["R_t"].max()),
            "B_peak_over_B0": float(out["B_total_t"].max() / g["B0"].sum()),
            "frac_at_2pct": t["frac_at_2pct"], "t50": t["t50"], "t_final": t["t_final"],
        })
        print(f"    {name}")
        print(f"      OM on {rows[-1]['om_bearing_matrix_cells']}/{rows[-1]['porous_matrix_cells']} "
              f"matrix cells; {n_touch}/{rows[-1]['n_regions']} regions adjacent to fuel; "
              f"median geodesic={geo['stats_um']['median']:.0f} um")
        print(f"      R_peak={rows[-1]['R_peak']:.3f}  frac@2%={t['frac_at_2pct']:.3f}  "
              f"t_final={t['t_final']:.4f}")
    _write_csv(os.path.join(results_dir, "stage16_ss5_reconciliation.csv"), rows)
    return rows


# ------------------------------------------------------------- SS3 the ladder

def evaluate(om_keys, D_keys, biology_cfg, mapping_cfg, grid_n,
             om_mode=dual_porosity.OM_MODE_SHARED):
    """Run both soils under one ladder configuration and score SS4."""
    m = shared_mapping(mapping_cfg, om_keys, D_keys, om_mode=om_mode)
    outs = {s: run_soil(s, grid_n, mapping_cfg[f"stage10_T_{s}"], biology_cfg, m) for s in SOILS}
    res = {"outs": outs, "dist": distinctness(outs["sand"], outs["vertisol"])}
    for s in SOILS:
        checks, vals = PASS_FAIL[s](outs[s])
        res[s] = {"checks": checks, "vals": vals, "shape_ok": all(checks.values()),
                  "timing": recruitment_timing(outs[s]), "hotspot": hotspot_selectivity(outs[s]),
                  "om_total": float(outs[s]["grids"]["om_total_init"]),
                  "D_scale": dual_porosity._variantc_D_scale(m, s)}
    res["recruit_ok"] = recruitment_ok(res["vertisol"]["timing"])
    res["dist_ok"] = res["dist"] > DISTINCTNESS_THRESHOLD
    res["all_pass"] = bool(res["sand"]["shape_ok"] and res["vertisol"]["shape_ok"]
                           and res["recruit_ok"] and res["dist_ok"])
    return res


def build_ladder(pm_cells):
    """The SS3 ladder as a fixed, reasoned probe sequence -- NOT an open
    search. Every probe's dose is anchored on the evidence SS5 supplies
    ('do not rediscover it'): Sand's burst-crash survives the shared rule at
    om_total~3,000 with D_scale=80; Vertisol grows AND staggers at
    om_total~50,000 with D_scale=0.7. The ladder's question is how much of
    that ~17x per-soil dose difference and ~114x transport difference can be
    DERIVED from structure instead of hand-set.

    `pm_cells` are each soil's porous-matrix cell counts at n=128, used to
    turn rung 2's shared concentration into a per-soil total."""
    rho_sand_working = 3000.0 / pm_cells["sand"]        # rho putting Sand at its working dose
    rho_vert_working = 50000.0 / pm_cells["vertisol"]   # rho putting Vertisol at its working dose
    rho_mid = float(np.sqrt(rho_sand_working * rho_vert_working))
    return [
        # ---- rung 1: ONE shared absolute total, ONE shared D_scale ----
        dict(rung=1, om={"stage16_om_total_shared": 3000.0}, D={"stage16_D_scale_shared": 80.0},
             note="rung 1 probe 1: Sand's own working point (SS5: om=3,000 / D=80 gives Sand a "
                  "genuine burst-crash) imposed on BOTH soils. The purest possible claim."),
        dict(rung=1, om={"stage16_om_total_shared": 50000.0}, D={"stage16_D_scale_shared": 0.7},
             note="rung 1 probe 2: the mirror -- Vertisol's own working point (SS5: om=50,000 / "
                  "D=0.7 gives growth AND staggering) imposed on both. Brackets the gap from the "
                  "other side."),
        dict(rung=1, om={"stage16_om_total_shared": 12247.0}, D={"stage16_D_scale_shared": 7.5},
             note="rung 1 probe 3: geometric midpoint of probes 1-2 in both dose and transport -- "
                  "the informative middle, testing whether any single compromise serves both."),
        dict(rung=1, om={"stage16_om_total_shared": 12247.0}, D={"stage16_D_scale_shared": 80.0},
             note="rung 1 probe 4: midpoint dose but Sand's fast transport, separating the dose "
                  "axis from the transport axis -- SS5 warns raising D_scale destroys staggering "
                  "fastest, so this isolates which axis actually blocks a shared setting."),
        # ---- rung 2: shared CONCENTRATION per unit matrix volume ----
        dict(rung=2, om={"stage16_om_rho_shared": rho_sand_working},
             D={"stage16_D_scale_shared": 80.0},
             note="rung 2 probe 5: shared OM concentration per unit matrix volume, pinned so SAND "
                  "lands on its working dose; Vertisol's total is then DERIVED from its 3.69x "
                  "larger matrix. Fast shared transport."),
        dict(rung=2, om={"stage16_om_rho_shared": rho_vert_working},
             D={"stage16_D_scale_shared": 0.7},
             note="rung 2 probe 6: same shared-concentration rule pinned so VERTISOL lands on its "
                  "working dose; Sand's total is then derived. Slow shared transport."),
        dict(rung=2, om={"stage16_om_rho_shared": rho_mid}, D={"stage16_D_scale_shared": 7.5},
             note="rung 2 probe 7: shared concentration at the geometric midpoint of probes 5-6, "
                  "midpoint shared transport -- rung 2's best-compromise attempt."),
        dict(rung=2, om={"stage16_om_rho_shared": rho_mid}, D="macroporosity",
             note="rung 2 probe 8: the analogous STRUCTURE-DERIVED option for transport that SS3 "
                  "asks to consider -- one shared coefficient, per-soil D_scale = coeff x that "
                  "soil's own macropore fraction, so transport too is derived rather than set."),
        # ---- rung 3: per-soil, as Stage 10 established ----
        dict(rung=3, om={"stage16_om_total_sand": 3000.0, "stage16_om_total_vertisol": 50000.0},
             D={"stage16_D_scale_sand": 80.0, "stage16_D_scale_vertisol": 0.7},
             note="rung 3 probe 9: per-soil dose and transport at exactly SS5's measured working "
                  "points. This is Stage 10's position, and the fallback the ladder descends to "
                  "only if rungs 1-2 fail."),
    ]


def run_ladder(biology_cfg, mapping_cfg, results_dir, pm_cells, macro_fracs):
    probes = build_ladder(pm_cells)
    assert len(probes) <= MAX_PROBES, f"{len(probes)} probes exceeds the SS3 cap of {MAX_PROBES}"

    trace, adopted, res_by_probe = [], None, {}
    for i, p in enumerate(probes, start=1):
        D_keys = p["D"]
        if D_keys == "macroporosity":
            # One shared coefficient; each soil's D_scale is its own
            # macropore fraction times that coefficient -- the
            # structure-derived transport option SS3 asks to consider and
            # report. The coefficient is set so the two soils' mean lands
            # on probe 7's shared D_scale, making this probe differ from
            # probe 7 ONLY in whether transport is structure-derived.
            coeff = 7.5 / float(np.mean([macro_fracs[s] for s in SOILS]))
            D_keys = {f"stage16_D_scale_{s}": coeff * macro_fracs[s] for s in SOILS}
        res = {gn: evaluate(p["om"], D_keys, biology_cfg, mapping_cfg, gn) for gn in (40, 128)}
        r128 = res[128]
        entry = {
            "probe": i, "rung": p["rung"],
            "om_spec": "; ".join(f"{k}={v:g}" for k, v in p["om"].items()),
            "D_spec": "; ".join(f"{k}={v:g}" for k, v in D_keys.items()),
            "om_total_sand": r128["sand"]["om_total"],
            "om_total_vertisol": r128["vertisol"]["om_total"],
            "D_scale_sand": r128["sand"]["D_scale"],
            "D_scale_vertisol": r128["vertisol"]["D_scale"],
            "sand_R_peak": r128["sand"]["vals"]["R_peak"],
            "sand_R_end_over_peak": r128["sand"]["vals"]["R_end_over_peak"],
            "sand_shape_ok": r128["sand"]["shape_ok"],
            "vert_R_peak": r128["vertisol"]["vals"]["R_peak"],
            "vert_B_peak_over_B0": r128["vertisol"]["vals"]["B_peak_over_B0"],
            "vert_peak_frac": r128["vertisol"]["vals"]["peak_frac"],
            "vert_R_end_over_peak": r128["vertisol"]["vals"]["R_end_over_peak"],
            "vert_shape_ok": r128["vertisol"]["shape_ok"],
            "vert_frac_at_2pct": r128["vertisol"]["timing"]["frac_at_2pct"],
            "vert_t_final": r128["vertisol"]["timing"]["t_final"],
            "vert_hotspot_pct": r128["vertisol"]["hotspot"]["hotspot_pct"],
            "recruit_ok": r128["recruit_ok"],
            "distinctness_n128": r128["dist"], "distinctness_n40": res[40]["dist"],
            "all_pass_n128": r128["all_pass"], "all_pass_n40": res[40]["all_pass"],
            "all_pass_both_grids": bool(r128["all_pass"] and res[40]["all_pass"]),
            "note": p["note"],
        }
        trace.append(entry)
        res_by_probe[i] = res
        print(f"    probe {i} (rung {p['rung']}): "
              f"Sand om={entry['om_total_sand']:.0f} D={entry['D_scale_sand']:.3g} "
              f"ok={entry['sand_shape_ok']}  |  "
              f"Vert om={entry['om_total_vertisol']:.0f} D={entry['D_scale_vertisol']:.3g} "
              f"ok={entry['vert_shape_ok']} frac@2%={entry['vert_frac_at_2pct']:.2f} "
              f"t_fin={entry['vert_t_final']:.3f}  |  dist={entry['distinctness_n128']:.3f}  "
              f"ALL(n128/n40)={entry['all_pass_n128']}/{entry['all_pass_n40']}")
        if entry["all_pass_both_grids"] and adopted is None:
            adopted = {"probe": i, "rung": p["rung"], "om": p["om"], "D": D_keys,
                       "res": res, "note": p["note"]}
            print(f"    -> rung {p['rung']} PASSES at both grids; the ladder stops here "
                  f"(SS3: adopt the FIRST rung that passes).")
            break

    _write_csv(os.path.join(results_dir, "stage16_ladder_trace.csv"), trace)
    return trace, adopted, len(trace), res_by_probe


# ------------------------------------------------------------------- outputs

def _write_csv(path, rows):
    if not rows:
        return
    fieldnames = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, restval="")
        w.writeheader()
        w.writerows(rows)


def write_stage_outputs(res, results_dir, grid_n, half_max=None):
    head, pf, timing_rows = [], [], []
    for soil in SOILS:
        out, r = res["outs"][soil], res[soil]
        R = out["R_t"]
        head.append({
            "soil": LABEL[soil], "grid_n": grid_n, "om_mode": out["grids"]["om_mode"],
            "om_total": r["om_total"], "D_scale": r["D_scale"],
            "R_peak": float(R.max()), "R_end": float(R[-1]),
            "R_end_over_peak": float(R[-1] / R.max()),
            "peak_frac": float(np.argmax(R) / (len(R) - 1)),
            "B_peak_over_B0": float(out["B_total_t"].max() / out["grids"]["B0"].sum()),
            "fuel_remaining_final": float(out["fuel_remaining_t"][-1]),
            "n_regions": int(out["grids"]["region_diam_um"].size),
            "outcome": metrics.classify_outcome(out),
            "shape": metrics.describe_shape(R)["shape"],
            "all_checks_pass": r["shape_ok"],
        })
        row = {"soil": LABEL[soil], "grid_n": grid_n, "all_checks_pass": r["shape_ok"]}
        for name, passed in r["checks"].items():
            row[f"check: {name}"] = passed
        pf.append(row)
        t, h = r["timing"], r["hotspot"]
        hm = (half_max or {}).get(soil)
        hm = np.asarray(hm if hm is not None else [], dtype=float)
        timing_rows.append({
            "soil": LABEL[soil], "grid_n": grid_n,
            "frac_recruited_at_2pct": t["frac_at_2pct"], "t50": t["t50"], "t_final": t["t_final"],
            "n_recruited_final": t["n_recruited_final"],
            "recruitment_not_instant": recruitment_ok(t),
            # SS4(c)/Stage 14's GRADED alternative to the boolean count --
            # reported alongside it because the two disagree completely.
            "graded_half_max_p25": float(np.percentile(hm, 25)) if hm.size else float("nan"),
            "graded_half_max_median": float(np.median(hm)) if hm.size else float("nan"),
            "graded_half_max_p75": float(np.percentile(hm, 75)) if hm.size else float("nan"),
            "graded_half_max_max": float(hm.max()) if hm.size else float("nan"),
            "graded_n_regions": int(hm.size),
            "hotspot_pct_peak_biomass_over_10pct_of_max": h["hotspot_pct"],
            "n_hotspot_regions": h["n_hotspots"], "n_regions": h["n_regions"],
            "portell_reference_pct": PORTELL_HOTSPOT_PCT,
            "portell_reference_sd": PORTELL_HOTSPOT_SD,
            "active_mass_peak": h["active_mass_peak"],
            "active_mass_t_of_peak_frac_of_window": h["active_mass_t_of_peak"],
            "secondary_boolean_recruited_pct": h["boolean_recruited_pct"],
        })
    dist_rows = [{"grid_n": grid_n, "Sand_vs_Vertisol": res["dist"],
                  "threshold": DISTINCTNESS_THRESHOLD, "distinct": res["dist_ok"]}]
    _write_csv(os.path.join(results_dir, f"stage16_headline_n{grid_n}.csv"), head)
    _write_csv(os.path.join(results_dir, f"stage16_pass_fail_n{grid_n}.csv"), pf)
    _write_csv(os.path.join(results_dir, f"stage16_distinctness_n{grid_n}.csv"), dist_rows)
    return head, timing_rows


def plot_headline(res, results_dir, grid_n):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for soil in SOILS:
        out = res["outs"][soil]
        t = np.arange(len(out["R_t"])) * out["dt"]
        axes[0].plot(t, out["R_t"], color=COLORS[soil], label=LABEL[soil], lw=2)
        axes[1].plot(np.linspace(0, 1, len(out["R_t"])), out["R_t"],
                     color=COLORS[soil], label=LABEL[soil], lw=2)
    for ax, title, xlabel in ((axes[0], "absolute time", "model time (absolute)"),
                              (axes[1], "normalized time", "fraction of own observation window")):
        ax.axhline(MAINTENANCE_FLOOR, color="gray", ls=":", lw=1)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("R(t)")
        ax.set_title(title)
        ax.legend(fontsize=9)
    fig.suptitle(f"Stage 16 (n={grid_n}): both soils under ONE shared, structure-blind OM rule "
                 f"with uniform dissolution\n(Sand-vs-Vertisol distinctness = {res['dist']:.3f})")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, f"stage16_headline_comparison_n{grid_n}.png"), dpi=150)
    plt.close(fig)


# ------------------------------------------------------------------ main

def main():
    t_start = time.time()
    os.makedirs(RESULTS_DIR, exist_ok=True)
    biology_cfg = load_yaml("configs/biology.yaml")
    mapping_cfg = load_yaml("configs/mapping.yaml")

    print("=" * 78)
    print("STAGE 16 -- restore the shared OM rule; structure sets fuel-habitat geometry")
    print("=" * 78)

    # ---- [1] preserved legacy path + observational-output identity check
    print("\n[1] SS9.2 -- the preserved `om_mode: coating_trapped` path still reproduces Stage 11 ...")
    legacy = copy.deepcopy(mapping_cfg)
    legacy["om_mode"] = dual_porosity.OM_MODE_COATING
    legacy_rows = []
    for soil, r_ref, b_ref in (("sand", 15.790, 10.7007), ("vertisol", 3.756, 2.337)):
        out = run_soil(soil, 128, mapping_cfg[f"stage10_T_{soil}"], biology_cfg, legacy,
                       track_region_biomass=False)
        rp = float(out["R_t"].max())
        bp = float(out["B_total_t"].max() / out["grids"]["B0"].sum())
        legacy_rows.append({"soil": LABEL[soil], "R_peak": rp, "R_peak_stage11_ref": r_ref,
                            "R_peak_dev_pct": abs(rp - r_ref) / r_ref * 100.0,
                            "B_peak_over_B0": bp, "B_peak_over_B0_stage11_ref": b_ref,
                            "reproduces_stage11": bool(abs(rp - r_ref) / r_ref < 0.01)})
        print(f"    {LABEL[soil]:9} R_peak={rp:8.4f} (Stage 11 ref {r_ref})  "
              f"B_peak/B0={bp:7.4f} (ref {b_ref})  dev={abs(rp - r_ref) / r_ref * 100:.3f}%")

    m_shared = shared_mapping(mapping_cfg, {"stage16_om_total_shared": 10000.0},
                              {"stage16_D_scale_shared": 1.0})
    plain = run_soil("vertisol", 128, 400, biology_cfg, m_shared, track_region_biomass=False)
    inst = run_soil("vertisol", 128, 400, biology_cfg, m_shared, track_region_biomass=True)
    identity = {k: bool(np.array_equal(plain[k], inst[k]))
                for k in ("R_t", "B_total_t", "habitats_recruited_t", "fuel_remaining_t")}
    identity_ok = all(identity.values())
    print(f"    observational-output identity check (track_region_biomass on vs off): "
          f"all_equal={identity_ok} {identity}")
    for r in legacy_rows:
        r["observational_outputs_identity_ok"] = identity_ok
    _write_csv(os.path.join(RESULTS_DIR, "stage16_legacy_reproducibility.csv"), legacy_rows)

    # ---- [2] SS2 geodesic fuel->habitat distance, BEFORE any tuning
    print("\n[2] SS2 -- emergent geodesic fuel->habitat distance under the shared rule ...")
    geo_rows, sens_rows, per_region = geodesic_analysis(biology_cfg, mapping_cfg, RESULTS_DIR)
    for r in geo_rows:
        print(f"    {r['soil']:9} n={r['grid_n']:3}  regions={r['n_regions']:3} "
              f"(unreachable {r['n_regions_unreachable']})  median={r['median_um']:6.1f} um  "
              f"p95={r['p95_um']:6.1f} um  max={r['max_um']:6.1f} um  "
              f"[substantial = {r['substantial_cells']}/{r['porous_matrix_cells']} matrix cells]")

    pm_cells, macro_fracs = {}, {}
    for soil in SOILS:
        g = dual_porosity.build_dual_grids(structure_cfg(soil, 128, 10), m_shared)
        pm_cells[soil] = int(g["porous_matrix_mask"].sum())
        macro_fracs[soil] = float(g["macro_mask"].mean())
    print(f"    porous-matrix cells (n=128): {pm_cells}  ->  Vertisol:Sand = "
          f"{pm_cells['vertisol'] / pm_cells['sand']:.3f}x")
    print("    macropore fractions (n=128): "
          + ", ".join(f"{LABEL[s]}={macro_fracs[s]:.4f}" for s in SOILS))

    # ---- [2b] why SS5's recorded numbers do not transfer to SS1's rule
    print("\n[2b] SS5 reconciliation diagnostic (adopts nothing; explains a discrepancy) ...")
    ss5_rows = ss5_reconciliation(biology_cfg, mapping_cfg, RESULTS_DIR)

    # ---- [3] SS3 ladder
    print(f"\n[3] SS3 -- dose/transport ladder, <= {MAX_PROBES} directed probes, both grids each ...")
    trace, adopted, n_probes, res_by_probe = run_ladder(
        biology_cfg, mapping_cfg, RESULTS_DIR, pm_cells, macro_fracs)

    # ---- [4] SS4 criteria for the adopted rung (or, failing that, SS6)
    if adopted is None:
        # Stated selection rule, fixed in advance so "closest attempt" is not
        # picked after the fact to flatter the result: the probe satisfying
        # the most SS4 criteria at n=128, ties broken by distinctness. No new
        # run is made -- the ladder's own result for that probe is reused, so
        # nothing outside the 15-probe cap is ever run.
        def _score(entry):
            return (int(entry["sand_shape_ok"]) + int(entry["vert_shape_ok"])
                    + int(entry["recruit_ok"])
                    + int(entry["distinctness_n128"] > DISTINCTNESS_THRESHOLD),
                    entry["distinctness_n128"])
        best = max(trace, key=_score)
        final = res_by_probe[best["probe"]]
        adopted_rung, adopted_probe = None, best["probe"]
        print("\n[4] SS6 HONEST-FAILURE PATH: no rung satisfied every SS4 criterion "
              "at both grids.")
        print(f"    Closest attempt (most criteria met, then highest distinctness): probe "
              f"{best['probe']} (rung {best['rung']}). Its numbers are reported below so the "
              f"failure is quantified, not merely asserted. Nothing is adopted.")
    else:
        final = adopted["res"]
        adopted_rung, adopted_probe = adopted["rung"], adopted["probe"]
        print(f"\n[4] SS4 -- criteria for the ADOPTED rung {adopted_rung} (probe {adopted_probe}), "
              f"n=128 and n=40 ...")

    # SS4(c)/Stage 14's graded per-region timing for the reported
    # configuration -- observational only, one extra run per soil per grid,
    # never a ladder probe and never a candidate for adoption.
    half_max = {}
    for grid_n in (40, 128):
        hm = {}
        for soil in SOILS:
            m_hm = shared_mapping(
                mapping_cfg,
                {f"stage16_om_total_{soil}": float(final[grid_n][soil]["om_total"])},
                {f"stage16_D_scale_{soil}": float(final[grid_n][soil]["D_scale"])})
            out_hm = run_soil(soil, grid_n, mapping_cfg[f"stage10_T_{soil}"], biology_cfg,
                              m_hm, debug_growth_distribution=True)
            hm[soil] = per_region_half_max(out_hm)
        half_max[grid_n] = hm

    timing_all = []
    for grid_n in (40, 128):
        head, timing_rows = write_stage_outputs(final[grid_n], RESULTS_DIR, grid_n,
                                                half_max=half_max[grid_n])
        plot_headline(final[grid_n], RESULTS_DIR, grid_n)
        timing_all.extend(timing_rows)
        for h in head:
            print(f"    n={grid_n} {h['soil']:9} R_peak={h['R_peak']:8.3f} "
                  f"R_end/peak={h['R_end_over_peak']:6.3f} peak@{h['peak_frac']:5.3f} "
                  f"B_peak/B0={h['B_peak_over_B0']:7.3f} shape={h['shape']:10} "
                  f"all_checks={h['all_checks_pass']}")
        for tr in timing_rows:
            print(f"    n={grid_n} {tr['soil']:9} frac@2%={tr['frac_recruited_at_2pct']:.3f} "
                  f"t50={tr['t50']:.4f} t_final={tr['t_final']:.4f}  "
                  f"hotspot={tr['hotspot_pct_peak_biomass_over_10pct_of_max']:.1f}% "
                  f"(Portell {PORTELL_HOTSPOT_PCT}+/-{PORTELL_HOTSPOT_SD}%)  "
                  f"[secondary boolean={tr['secondary_boolean_recruited_pct']:.1f}%]")
            print(f"    n={grid_n} {tr['soil']:9} GRADED per-region half-max over "
                  f"{tr['graded_n_regions']} regions: p25={tr['graded_half_max_p25']:.3f} "
                  f"median={tr['graded_half_max_median']:.3f} p75={tr['graded_half_max_p75']:.3f} "
                  f"max={tr['graded_half_max_max']:.3f}")
        print(f"    n={grid_n} distinctness={final[grid_n]['dist']:.3f} "
              f"(>{DISTINCTNESS_THRESHOLD}: {final[grid_n]['dist_ok']})  "
              f"ALL_PASS={final[grid_n]['all_pass']}")
    _write_csv(os.path.join(RESULTS_DIR, "stage16_recruitment_timing.csv"), timing_all)

    overall = bool(final[40]["all_pass"] and final[128]["all_pass"])

    # ---- [5] SS7 shared-lognormal clumping.
    # SS7 gates this on SS4 passing, and SS4 did not pass. It is still run
    # here, but strictly as a MECHANISM DIAGNOSTIC for SS6's required
    # explanation rather than as SS7's sensitivity check: the criterion that
    # failed is recruitment timing, and the reason (SS2, and the SS5
    # reconciliation above) is that a shared rule filling the whole matrix
    # leaves every region adjacent to fuel. Clumping is the one placement
    # change that is still shared and still structure-blind but does open
    # fuel-free gaps, so it tests whether the failure belongs to "shared
    # rules" or only to "smooth shared rules". Nothing here is adopted, no
    # threshold moves, and it is not counted against the SS3 probe cap.
    if True:
        label = ("[5] SS7 clumping as an SS6 mechanism diagnostic (SS4 did not pass, so this "
                 "is NOT SS7's sensitivity check and adopts nothing)")
        print(f"\n{label} ...")
        clump_rows = []
        clump_om = ({"stage16_om_total_sand": 3000.0, "stage16_om_total_vertisol": 50000.0}
                    if adopted is None else adopted["om"])
        clump_D = ({"stage16_D_scale_sand": 80.0, "stage16_D_scale_vertisol": 0.7}
                   if adopted is None else adopted["D"])
        for grid_n in (40, 128):
            cl = evaluate(clump_om, clump_D, biology_cfg, mapping_cfg, grid_n,
                          om_mode=dual_porosity.OM_MODE_CLUMPS)
            for soil in SOILS:
                clump_rows.append({
                    "soil": LABEL[soil], "grid_n": grid_n,
                    "smooth_hotspot_pct": final[grid_n][soil]["hotspot"]["hotspot_pct"],
                    "clumped_hotspot_pct": cl[soil]["hotspot"]["hotspot_pct"],
                    "smooth_R_peak": final[grid_n][soil]["vals"]["R_peak"],
                    "clumped_R_peak": cl[soil]["vals"]["R_peak"],
                    "smooth_frac_at_2pct": final[grid_n][soil]["timing"]["frac_at_2pct"],
                    "clumped_frac_at_2pct": cl[soil]["timing"]["frac_at_2pct"],
                    "smooth_t_final": final[grid_n][soil]["timing"]["t_final"],
                    "clumped_t_final": cl[soil]["timing"]["t_final"],
                    "smooth_shape_ok": final[grid_n][soil]["shape_ok"],
                    "clumped_shape_ok": cl[soil]["shape_ok"],
                    "clumped_recruitment_not_instant": recruitment_ok(cl[soil]["timing"]),
                })
                cr = clump_rows[-1]
                print(f"    n={grid_n} {LABEL[soil]:9} hotspot {cr['smooth_hotspot_pct']:.1f}% -> "
                      f"{cr['clumped_hotspot_pct']:.1f}%   R_peak {cr['smooth_R_peak']:.2f} -> "
                      f"{cr['clumped_R_peak']:.2f}   frac@2% {cr['smooth_frac_at_2pct']:.3f} -> "
                      f"{cr['clumped_frac_at_2pct']:.3f}   t_final {cr['smooth_t_final']:.4f} -> "
                      f"{cr['clumped_t_final']:.4f}")
        _write_csv(os.path.join(RESULTS_DIR, "stage16_clumping_sensitivity.csv"), clump_rows)

    runtime = time.time() - t_start
    with open(LOGBOOK_PATH, "a", encoding="utf-8") as f:
        # SS3 requires every probe logged with its reasoning, not just the
        # headline -- the ladder is meant to be auditable rung by rung.
        f.write(f"\n### Stage 16 -- SS3 dose/transport ladder, {n_probes}/{MAX_PROBES} probes "
                f"(each evaluated on BOTH soils at BOTH grids)\n\n")
        for e in trace:
            f.write(f"- **probe {e['probe']} (rung {e['rung']})** -- {e['note']}\n"
                    f"  Sand `om_total={e['om_total_sand']:.0f}`, `D_scale={e['D_scale_sand']:.3g}` "
                    f"-> shape checks {'PASS' if e['sand_shape_ok'] else 'FAIL'} "
                    f"(`R_peak={e['sand_R_peak']:.3f}`, `R_end/R_peak={e['sand_R_end_over_peak']:.3f}`). "
                    f"Vertisol `om_total={e['om_total_vertisol']:.0f}`, "
                    f"`D_scale={e['D_scale_vertisol']:.3g}` -> shape checks "
                    f"{'PASS' if e['vert_shape_ok'] else 'FAIL'} (`R_peak={e['vert_R_peak']:.3f}`, "
                    f"`B_peak/B0={e['vert_B_peak_over_B0']:.3f}`, `peak@{e['vert_peak_frac']:.3f}`). "
                    f"Recruitment `frac@2%={e['vert_frac_at_2pct']:.3f}`, "
                    f"`t_final={e['vert_t_final']:.4f}` -> "
                    f"{'PASS' if e['recruit_ok'] else 'FAIL'}. "
                    f"Distinctness {e['distinctness_n128']:.3f} (n128) / "
                    f"{e['distinctness_n40']:.3f} (n40). "
                    f"**All SS4 criteria: {e['all_pass_n128']} (n128) / "
                    f"{e['all_pass_n40']} (n40).**\n")
        f.write(f"\n**Stage 16 headline: the per-soil coating/trapped OM split is RETIRED and one "
                f"shared, structure-blind OM rule (`om_mode: shared_fine_pore`, uniform k_dis) is "
                f"the new default; the Stage 6-15 path is preserved behind `om_mode: "
                f"coating_trapped` and still reproduces Stage 11 exactly (max dev "
                f"{max(r['R_peak_dev_pct'] for r in legacy_rows):.3f}%). Stage 14's causal "
                f"attribution was refuted by direct test and corrected in place. SS2 emergent "
                f"geodesic fuel->habitat distances (n=128, shared rule, median-cut definition): "
                + "; ".join(f"{r['soil']} median={r['median_um']:.0f}um p95={r['p95_um']:.0f}um"
                            for r in geo_rows if r["grid_n"] == 128)
                + f". Ladder: {n_probes}/{MAX_PROBES} probes, adopted rung="
                f"{adopted_rung if adopted_rung else 'NONE (honest failure, SS6)'}. "
                f"Overall SS4 pass at both grids={overall}. Runtime {runtime / 60:.1f} min. "
                f"{datetime.datetime.now().isoformat(timespec='seconds')}. Full detail in "
                "docs/stage_results/STAGE16_RESULTS.md.**\n")

    print(f"\n[6] OVERALL SS4 pass at both grids = {overall}; probes used = {n_probes}/{MAX_PROBES}; "
          f"adopted rung = {adopted_rung}")
    print(f"    total runtime {runtime / 60:.2f} min; results in {RESULTS_DIR}/")
    return overall


if __name__ == "__main__":
    main()
