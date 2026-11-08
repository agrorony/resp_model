"""Post-processing for the fixed Stage 23 respiration design."""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

import chem23
import metrics
import stage23_run as S23
import timemap
from model import load_yaml


SOIL_CROPS = {
    "S": [f"{scan}_c{index}" for scan in ("SA", "SB") for index in range(4)],
    "V": [f"{scan}_c{index}" for scan in ("VA", "VB") for index in range(4)],
    "L": [f"{scan}_c{index}" for scan in ("LA", "LB", "LC") for index in range(4)],
}
SPECIMEN_CROPS = {
    scan: [f"{scan}_c{index}" for index in range(4)]
    for scan in S23.SCAN_IDS
}
SOIL_NAMES = {"S": "Sand", "V": "Vertisol", "L": "Loess"}
PERIOD_NAMES = ("P1", "P2", "P3", "P4")
_ARRAY_CACHE: dict[str, dict[str, np.ndarray]] = {}
_RUNS: dict[str, dict] | None = None
_RUN_INDEX: dict[tuple[str, str, str], str] | None = None


def _read_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as file:
        return json.load(file)


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite Stage 23 analysis output: {path}")
    path.write_text(
        json.dumps(value, indent=2, default=S23._json_default) + "\n",
        encoding="utf-8",
    )


def _load_runs() -> tuple[dict[str, dict], dict[tuple[str, str, str], str]]:
    global _RUNS, _RUN_INDEX
    if _RUNS is None or _RUN_INDEX is None:
        step2 = _read_json(S23.OUT / "step2.json")
        _RUNS = step2["runs"]
        _RUN_INDEX = {}
        for key, summary in _RUNS.items():
            index = (summary["arm"], summary["crop_id"], summary["chemistry"])
            if index in _RUN_INDEX:
                raise ValueError(f"Duplicate Stage 23 run index: {index}")
            _RUN_INDEX[index] = key
    return _RUNS, _RUN_INDEX


def _get_run(arm: str, crop_id: str, chemistry: str) -> dict:
    runs, index = _load_runs()
    key = index.get((arm, crop_id, chemistry))
    if key is None:
        raise KeyError(f"Missing Stage 23 run: {arm}, {crop_id}, {chemistry}")
    return runs[key]


def _arrays(summary: dict) -> dict[str, np.ndarray]:
    key = summary["key"]
    if key not in _ARRAY_CACHE:
        with np.load(S23.RUNS / f"{key}.npz", allow_pickle=False) as data:
            _ARRAY_CACHE[key] = {name: data[name] for name in data.files}
    return _ARRAY_CACHE[key]


def _normalized_curve(summary: dict) -> np.ndarray:
    arrays = _arrays(summary)
    m0 = float(load_yaml("configs/biology.yaml")["m0"])
    denominator = m0 * float(arrays["B"][0])
    if denominator <= 0:
        raise ValueError(f"Nonpositive m0*B0 for {summary['key']}")
    return np.asarray(arrays["R"], dtype=float) / denominator


def _mean_curve(curves: list[np.ndarray]) -> np.ndarray:
    if not curves:
        raise ValueError("Cannot average an empty curve collection")
    return np.mean(np.stack(curves), axis=0)


def _distance(curve_a: np.ndarray, curve_b: np.ndarray, steps: int) -> float:
    return float(metrics.pairwise_distance(curve_a[:steps], curve_b[:steps]))


def _within_cross(
    soil_a: str,
    soil_b: str,
    curves: dict[str, np.ndarray],
    steps: int,
) -> dict:
    within_pairs = []
    for soil in (soil_a, soil_b):
        ids = SOIL_CROPS[soil]
        for first, second in itertools.combinations(ids, 2):
            within_pairs.append({
                "soil": soil,
                "a": first,
                "b": second,
                "distance": _distance(curves[first], curves[second], steps),
            })
    cross_pairs = []
    for first in SOIL_CROPS[soil_a]:
        for second in SOIL_CROPS[soil_b]:
            cross_pairs.append({
                "a": first,
                "b": second,
                "distance": _distance(curves[first], curves[second], steps),
            })
    return {
        "within_max": max(within_pairs, key=lambda pair: pair["distance"]),
        "cross_min": min(cross_pairs, key=lambda pair: pair["distance"]),
        "within_pairs": within_pairs,
        "cross_pairs": cross_pairs,
    }


def _distinctness_for_arm(arm: str, chemistry_by_soil: dict[str, str]) -> dict:
    summaries = {
        crop: _get_run(arm, crop, chemistry_by_soil[soil])
        for soil, crops in SOIL_CROPS.items()
        for crop in crops
    }
    curves = {crop: _normalized_curve(run) for crop, run in summaries.items()}
    mean_curves = {
        soil: _mean_curve([curves[crop] for crop in crops])
        for soil, crops in SOIL_CROPS.items()
    }
    alive = {
        soil: {
            "alive": sum(bool(summaries[crop]["alive_primary"]) for crop in crops),
            "total": len(crops),
            "fraction": float(np.mean([summaries[crop]["alive_primary"] for crop in crops])),
            "crops": {crop: bool(summaries[crop]["alive_primary"]) for crop in crops},
        }
        for soil, crops in SOIL_CROPS.items()
    }
    pair_results = {}
    for soil_a, soil_b in (("S", "V"), ("S", "L"), ("V", "L")):
        pair = f"{soil_a}-{soil_b}"
        values = {}
        for steps in (1500, 3000):
            separation = _within_cross(soil_a, soil_b, curves, steps)
            mean_distance = _distance(mean_curves[soil_a], mean_curves[soil_b], steps)
            values[str(steps)] = {
                "soil_mean_distance": mean_distance,
                "max_within": separation["within_max"],
                "min_cross": separation["cross_min"],
                "within_pair_count": len(separation["within_pairs"]),
                "cross_pair_count": len(separation["cross_pairs"]),
            }
        a_alive = alive[soil_a]["fraction"]
        b_alive = alive[soil_b]["fraction"]
        criteria = {
            "at_least_75_percent_alive_each": a_alive >= 0.75 and b_alive >= 0.75,
            "mean_distance_gt_0_3": values["1500"]["soil_mean_distance"] > 0.3,
            "max_within_lt_half_min_cross": (
                values["1500"]["max_within"]["distance"]
                < 0.5 * values["1500"]["min_cross"]["distance"]
            ),
        }
        pair_results[pair] = {
            "soil_a_alive_fraction": a_alive,
            "soil_b_alive_fraction": b_alive,
            "by_steps": values,
            "criteria": criteria,
            "distinct": all(criteria.values()),
        }
    return {
        "arm": arm,
        "alive_primary": alive,
        "soil_mean_curves": mean_curves,
        "pairs": pair_results,
    }


def step3() -> dict:
    if (S23.OUT / "step3.json").exists():
        raise FileExistsError("Stage 23 Step 3 output already exists")
    arm_a = _distinctness_for_arm("A", {soil: "chem_eq" for soil in SOIL_CROPS})
    stage22 = _read_json(Path("results/stage22_scan_structure/step3.json"))
    sv = arm_a["pairs"]["S-V"]
    stage22_comparison = {}
    for steps in ("1500", "3000"):
        current = sv["by_steps"][steps]
        prior = {
            "soil_mean_distance": stage22["soil_mean_distance"][steps],
            "max_within": stage22["within_soil"][steps]["max_all"],
            "min_cross": stage22["cross_soil"][steps]["min"],
        }
        observed = {
            "soil_mean_distance": current["soil_mean_distance"],
            "max_within": current["max_within"]["distance"],
            "min_cross": current["min_cross"]["distance"],
        }
        stage22_comparison[steps] = {
            "stage22": prior,
            "stage23_arm_a": observed,
            "absolute_differences": {
                name: abs(observed[name] - prior[name]) for name in prior
            },
            "equal_to_1e-8": all(
                abs(observed[name] - prior[name]) < 1e-8 for name in prior
            ),
            "equal_at_4dp": all(
                f"{observed[name]:.4f}" == f"{prior[name]:.4f}" for name in prior
            ),
            "K_density_note": (
                "Stage 23 uses the fixed 0.16445; Stage 22 used 0.16445337235927582."
            ),
        }
    result = {
        "stage": 23,
        "step": 3,
        "distance": "metrics.pairwise_distance on R/(m0*B[0])",
        "horizons_steps": [1500, 3000],
        "arm_a": {
            "alive_primary": arm_a["alive_primary"],
            "pairs": arm_a["pairs"],
            "soil_mean_curves": {
                soil: curve.tolist() for soil, curve in arm_a["soil_mean_curves"].items()
            },
        },
        "stage22_sand_vertisol_comparison": stage22_comparison,
    }
    _write_json(S23.OUT / "step3.json", result)
    print(json.dumps({
        "distinctness": {
            pair: {
                "1500": data["by_steps"]["1500"]["soil_mean_distance"],
                "3000": data["by_steps"]["3000"]["soil_mean_distance"],
                "distinct": data["distinct"],
            }
            for pair, data in arm_a["pairs"].items()
        },
        "stage22_SV_comparison": stage22_comparison,
    }, indent=2))
    return result


def _physiology(summary: dict, r_lit: float) -> dict:
    arrays = _arrays(summary)
    run_n = int(summary["n"]) ** 3
    initial_om = float(summary["OM_density"]) * run_n
    output = timemap.physiological_outputs(
        arrays["R"],
        arrays["B"],
        arrays["uptake_dt"],
        initial_om,
        r_lit,
        dt=S23.DT,
    )
    cumulative_volume = output["cumulative_respiration_model_units"] / run_n
    density = float(chem23.CHEMISTRY[summary["soil_code"]]["core_BD_g_cm3"])
    return {
        **output,
        "cumulative_R_per_volume": cumulative_volume,
        "cumulative_R_per_mass": cumulative_volume / density,
        "core_BD_g_cm3": density,
        "initial_OM_total": initial_om,
    }


def _cross_scalar(summary: dict, metric: str) -> float:
    if metric == "R_peak_over_floor":
        return float(summary["R_peak_over_floor"])
    if metric == "peak_time_fraction":
        return float(summary["peak_frac"])
    if metric == "R_end_over_peak":
        return float(summary["R_end_over_peak"])
    physiology = _physiology(summary, S23.R_LIT_DAYS["primary"])
    if metric == "cumulative_R_per_volume":
        return float(physiology["cumulative_R_per_volume"])
    if metric == "cumulative_R_per_mass":
        return float(physiology["cumulative_R_per_mass"])
    if metric == "qCO2_days10_15":
        return float(physiology["qCO2_per_day"])
    raise ValueError(f"Unknown cross scalar: {metric}")


def _mean_sd(values: list[float] | np.ndarray) -> dict:
    data = np.asarray(values, dtype=float)
    if not data.size or not np.all(np.isfinite(data)):
        raise ValueError("Expected finite values for mean and SD")
    return {
        "mean": float(data.mean()),
        "sd": float(data.std(ddof=1)) if data.size > 1 else 0.0,
        "n": int(data.size),
    }


def _mean_optional(values: list[float | None]) -> float | None:
    defined = [float(value) for value in values if value is not None]
    return float(np.mean(defined)) if defined else None


def _two_way_decomposition(values: dict[tuple[str, str], list[float]]) -> dict:
    rows = ("S", "V", "L")
    columns = ("chem_S", "chem_V", "chem_L")
    n_by_row = {row: len(values[(row, columns[0])]) for row in rows}
    if any(len(values[(row, column)]) != n_by_row[row] for row in rows for column in columns):
        raise ValueError("Cross is not balanced across chemistry within each structure row")
    observations = [
        (row, column, value)
        for row in rows
        for column in columns
        for value in values[(row, column)]
    ]
    grand = float(np.mean([value for _, _, value in observations]))
    row_means = {
        row: float(np.mean([value for r, _, value in observations if r == row]))
        for row in rows
    }
    column_means = {
        column: float(np.mean([value for _, c, value in observations if c == column]))
        for column in columns
    }
    cell_means = {
        (row, column): float(np.mean(values[(row, column)]))
        for row in rows
        for column in columns
    }
    structure_ss = sum(
        3 * n_by_row[row] * (row_means[row] - grand) ** 2 for row in rows
    )
    chemistry_ss = sum(
        sum(n_by_row.values()) * (column_means[column] - grand) ** 2
        for column in columns
    )
    interaction_ss = sum(
        n_by_row[row]
        * (
            cell_means[(row, column)]
            - row_means[row]
            - column_means[column]
            + grand
        ) ** 2
        for row in rows
        for column in columns
    )
    residual_ss = sum(
        (value - cell_means[(row, column)]) ** 2
        for row, column, value in observations
    )
    total_ss = sum((value - grand) ** 2 for _, _, value in observations)
    parts = {
        "structure": float(structure_ss),
        "chemistry": float(chemistry_ss),
        "interaction": float(interaction_ss),
        "residual": float(residual_ss),
    }
    if abs(sum(parts.values()) - total_ss) > max(total_ss, 1.0) * 1e-9:
        raise ArithmeticError("Two-way sums of squares do not partition the total")
    return {
        "sum_of_squares": parts,
        "total_sum_of_squares": float(total_ss),
        "share_of_total": {
            key: value / total_ss if total_ss > 0 else 0.0
            for key, value in parts.items()
        },
        "grand_mean": grand,
        "row_means": row_means,
        "column_means": column_means,
    }


def _cross_cell_curves() -> dict[tuple[str, str], list[tuple[str, np.ndarray]]]:
    result = {}
    for soil, crops in SOIL_CROPS.items():
        for chemistry in ("chem_S", "chem_V", "chem_L"):
            result[(soil, chemistry)] = [
                (crop, _normalized_curve(_get_run("B" if chemistry == f"chem_{soil}" else "C", crop, chemistry)))
                for crop in crops
            ]
    return result


def _write_cross_figure(cell_curves: dict, target: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 3, figsize=(13, 10), sharex=True, sharey=True)
    step_axis = np.arange(S23.T_STEPS)
    for row_index, soil in enumerate(("S", "V", "L")):
        for col_index, chemistry in enumerate(("chem_S", "chem_V", "chem_L")):
            ax = axes[row_index, col_index]
            crop_runs = []
            for crop in SOIL_CROPS[soil]:
                arm = "B" if chemistry == f"chem_{soil}" else "C"
                summary = _get_run(arm, crop, chemistry)
                values = _arrays(summary)["R"] / (int(summary["n"]) ** 3)
                ax.plot(step_axis, values, color="#8c9aa8", alpha=0.25, linewidth=0.65)
                crop_runs.append(values)
            mean = np.mean(np.stack(crop_runs), axis=0)
            ax.plot(step_axis, mean, color="#b0232a", linewidth=1.6)
            ax.set_title(f"{SOIL_NAMES[soil]} × {chemistry}", fontsize=9)
            ax.grid(alpha=0.2)
    for ax in axes[-1, :]:
        ax.set_xlabel("Model step")
    for ax in axes[:, 0]:
        ax.set_ylabel("R per voxel")
    fig.tight_layout()
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite Stage 23 figure: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=140)
    plt.close(fig)


def step4() -> dict:
    target = S23.OUT / "step4.json"
    if target.exists():
        raise FileExistsError(f"Stage 23 Step 4 output already exists: {target}")
    cell_curves = _cross_cell_curves()
    metrics_to_analyze = (
        "R_peak_over_floor",
        "peak_time_fraction",
        "cumulative_R_per_volume",
        "cumulative_R_per_mass",
        "R_end_over_peak",
        "qCO2_days10_15",
    )
    values_by_metric = {metric: {} for metric in metrics_to_analyze}
    cells = {}
    for soil in ("S", "V", "L"):
        for chemistry in ("chem_S", "chem_V", "chem_L"):
            runs = [
                _get_run("B" if chemistry == f"chem_{soil}" else "C", crop, chemistry)
                for crop in SOIL_CROPS[soil]
            ]
            key = f"{soil}×{chemistry}"
            cells[key] = {
                "structure_soil": soil,
                "chemistry": chemistry,
                "n_crops": len(runs),
                "metrics": {},
            }
            for metric in metrics_to_analyze:
                values = [_cross_scalar(summary, metric) for summary in runs]
                values_by_metric[metric][(soil, chemistry)] = values
                cells[key]["metrics"][metric] = _mean_sd(values)

    decomposition = {
        metric: _two_way_decomposition(values)
        for metric, values in values_by_metric.items()
    }
    curve_distances = {}
    for steps in (1500, 3000):
        structure_pairs = []
        chemistry_pairs = []
        for soil in ("S", "V", "L"):
            for first, second in itertools.combinations(("chem_S", "chem_V", "chem_L"), 2):
                curve_a = _mean_curve([curve for _, curve in cell_curves[(soil, first)]])
                curve_b = _mean_curve([curve for _, curve in cell_curves[(soil, second)]])
                structure_pairs.append(_distance(curve_a, curve_b, steps))
        for chemistry in ("chem_S", "chem_V", "chem_L"):
            for first, second in itertools.combinations(("S", "V", "L"), 2):
                curve_a = _mean_curve([curve for _, curve in cell_curves[(first, chemistry)]])
                curve_b = _mean_curve([curve for _, curve in cell_curves[(second, chemistry)]])
                chemistry_pairs.append(_distance(curve_a, curve_b, steps))
        curve_distances[str(steps)] = {
            "same_structure_different_chemistry": _mean_sd(structure_pairs),
            "same_chemistry_different_structure": _mean_sd(chemistry_pairs),
            "structure_pair_distances": structure_pairs,
            "chemistry_pair_distances": chemistry_pairs,
        }
    dominant = {}
    for metric, result in decomposition.items():
        shares = result["share_of_total"]
        dominant[metric] = (
            "structure" if shares["structure"] > shares["chemistry"] else "chemistry"
        )
    chemistry_dominant = [
        metric for metric in dominant if dominant[metric] == "chemistry"
    ]
    structure_dominant = [
        metric for metric in dominant if dominant[metric] == "structure"
    ]
    factors = []
    if chemistry_dominant:
        factors.append("chemistry for " + ", ".join(chemistry_dominant))
    if structure_dominant:
        factors.append("structure for " + ", ".join(structure_dominant))
    sentence = "The larger main-effect share is " + " and ".join(factors) + "."
    arm_b_distinctness = _distinctness_for_arm(
        "B", {"S": "chem_S", "V": "chem_V", "L": "chem_L"}
    )
    figure = S23.OUT / "step4_cross.png"
    _write_cross_figure(cell_curves, figure)
    result = {
        "stage": 23,
        "step": 4,
        "arms": "B and C",
        "cross_cell_means_and_sd": cells,
        "two_way_decomposition": decomposition,
        "curve_level_distances": curve_distances,
        "arm_b_distinctness": {
            "alive_primary": arm_b_distinctness["alive_primary"],
            "pairs": arm_b_distinctness["pairs"],
            "soil_mean_curves": {
                soil: curve.tolist()
                for soil, curve in arm_b_distinctness["soil_mean_curves"].items()
            },
        },
        "dominant_main_effect_by_output": dominant,
        "plain_sentence": sentence,
        "figure": str(figure),
    }
    _write_json(target, result)
    print(sentence)
    print(json.dumps({
        metric: {
            factor: result["share_of_total"][factor]
            for factor in ("structure", "chemistry", "interaction", "residual")
        }
        for metric, result in decomposition.items()
    }, indent=2))
    return result


def _delivery_split(summary: dict) -> dict:
    diffusion = (
        float(summary["diffusion_into_wet_habitat"])
        + float(summary["diffusion_into_drained_habitat"])
    )
    exchange = (
        float(summary["exchange_into_wet_habitat"])
        + float(summary["exchange_into_drained_habitat"])
    )
    total = diffusion + exchange
    result = {
        "diffusion_total": diffusion,
        "exchange_total": exchange,
        "total": total,
        "fractions_defined": abs(total) > 1e-300,
    }
    if result["fractions_defined"]:
        result["diffusion_fraction"] = diffusion / total
        result["exchange_fraction"] = exchange / total
    else:
        result["diffusion_fraction"] = None
        result["exchange_fraction"] = None
    return result


def _ordering_check(arm: str, sand_arm_d: str, chemistry: str) -> dict:
    if sand_arm_d == "stage21_pore_fill":
        sand_runs = [_get_run("D", crop, chemistry) for crop in SOIL_CROPS["S"]]
    else:
        sand_runs = [_get_run(arm, crop, chemistry) for crop in SOIL_CROPS["S"]]
    sand_curve = _mean_curve([_normalized_curve(run) for run in sand_runs])
    vertisol_chem = "chem_eq" if arm == "A" else "chem_V"
    loess_chem = "chem_eq" if arm == "A" else "chem_L"
    vertisol_curve = _mean_curve([
        _normalized_curve(_get_run(arm, crop, vertisol_chem))
        for crop in SOIL_CROPS["V"]
    ])
    loess_curve = _mean_curve([
        _normalized_curve(_get_run(arm, crop, loess_chem))
        for crop in SOIL_CROPS["L"]
    ])
    values = {}
    for steps in (1500, 3000):
        sv = _distance(sand_curve, vertisol_curve, steps)
        sl = _distance(sand_curve, loess_curve, steps)
        values[str(steps)] = {
            "S_V_distance": sv,
            "S_L_distance": sl,
            "closer_pair": "S-V" if sv < sl else "S-L",
        }
    return {
        "arm": arm,
        "sand_rule": sand_arm_d,
        "chemistry": chemistry,
        "by_steps": values,
    }


def step5() -> dict:
    target = S23.OUT / "step5.json"
    if target.exists():
        raise FileExistsError(f"Stage 23 Step 5 output already exists: {target}")
    arm_d = {}
    for chemistry in ("chem_eq", "chem_S"):
        base_arm = "A" if chemistry == "chem_eq" else "B"
        pairs = []
        for crop in SOIL_CROPS["S"]:
            base = _get_run(base_arm, crop, chemistry)
            alternative = _get_run("D", crop, chemistry)
            base_curve = _normalized_curve(base)
            alt_curve = _normalized_curve(alternative)
            pairs.append({
                "crop_id": crop,
                "curve_distance": {
                    str(steps): _distance(base_curve, alt_curve, steps)
                    for steps in (1500, 3000)
                },
                "base": {
                    "R_peak_over_floor": base["R_peak_over_floor"],
                    "peak_frac": base["peak_frac"],
                    **_delivery_split(base),
                },
                "pore_fill": {
                    "R_peak_over_floor": alternative["R_peak_over_floor"],
                    "peak_frac": alternative["peak_frac"],
                    **_delivery_split(alternative),
                },
            })
        arm_d[chemistry] = {
            "paired_crop_results": pairs,
            "mean_curve_distance": {
                str(steps): float(np.mean([pair["curve_distance"][str(steps)] for pair in pairs]))
                for steps in (1500, 3000)
            },
            "mean_base_peak_over_floor": float(np.mean([p["base"]["R_peak_over_floor"] for p in pairs])),
            "mean_pore_fill_peak_over_floor": float(np.mean([p["pore_fill"]["R_peak_over_floor"] for p in pairs])),
            "mean_base_peak_frac": float(np.mean([p["base"]["peak_frac"] for p in pairs])),
            "mean_pore_fill_peak_frac": float(np.mean([p["pore_fill"]["peak_frac"] for p in pairs])),
            "mean_base_delivery": {
                key: _mean_optional([p["base"][key] for p in pairs])
                for key in (
                    "diffusion_total",
                    "exchange_total",
                    "diffusion_fraction",
                    "exchange_fraction",
                )
            },
            "mean_pore_fill_delivery": {
                key: _mean_optional([p["pore_fill"][key] for p in pairs])
                for key in (
                    "diffusion_total",
                    "exchange_total",
                    "diffusion_fraction",
                    "exchange_fraction",
                )
            },
        }
    ordering = [
        _ordering_check("A", "stage22_solid_lt", "chem_eq"),
        _ordering_check("A", "stage21_pore_fill", "chem_eq"),
        _ordering_check("B", "stage22_solid_lt", "chem_S"),
        _ordering_check("B", "stage21_pore_fill", "chem_S"),
    ]
    for arm in ("A", "B"):
        base = next(row for row in ordering if row["arm"] == arm and row["sand_rule"] == "stage22_solid_lt")
        alternative = next(row for row in ordering if row["arm"] == arm and row["sand_rule"] == "stage21_pore_fill")
        for steps in ("1500", "3000"):
            alternative["by_steps"][steps]["ordering_changed"] = (
                alternative["by_steps"][steps]["closer_pair"]
                != base["by_steps"][steps]["closer_pair"]
            )
    domains = {}
    for scan in ("LA", "LB", "LC"):
        domain = _get_run("E", f"{scan}_domain256", "chem_eq")
        domain_curve = _normalized_curve(domain)
        crop_curves = [
            _normalized_curve(_get_run("A", crop, "chem_eq"))
            for crop in SPECIMEN_CROPS[scan]
        ]
        mean_100 = _mean_curve(crop_curves)
        domains[scan] = {
            "domain_run": domain["key"],
            "n_100_crops": len(crop_curves),
            "domain_vs_100_crop_mean_distance": {
                str(steps): _distance(domain_curve, mean_100, steps)
                for steps in (1500, 3000)
            },
            "domain_R_peak_over_floor": domain["R_peak_over_floor"],
            "domain_peak_frac": domain["peak_frac"],
        }
    result = {
        "stage": 23,
        "step": 5,
        "arm_d_sand_rule_check": arm_d,
        "sand_ordering_with_arm_d_substitution": ordering,
        "loess_domain_size_check": domains,
        "alternative_fines_rule_adopted": False,
    }
    _write_json(target, result)
    return result


def classify_shape_stage23(means: np.ndarray) -> str:
    values = np.asarray(means, dtype=float)
    if values.shape != (4,) or not np.all(np.isfinite(values)):
        raise ValueError("Four finite period means are required")
    maximum = float(values.max())
    maxima = set(np.flatnonzero(values == maximum).tolist())
    if maxima.intersection((1, 2)) and values[3] < 0.9 * maximum:
        return "rise-then-fall"
    if (
        maxima.intersection((2, 3))
        and values[3] >= 0.9 * maximum
        and maximum > 1.1 * values[0]
    ):
        return "rise"
    if 0 in maxima and values[3] < 0.9 * values[0]:
        return "fall"
    return "flat"


def _physiology_rows(arm: str, soil: str, r_lit: float, crops: list[str]) -> list[dict]:
    chemistry = "chem_eq" if arm == "A" else f"chem_{soil}"
    rows = []
    for crop in crops:
        summary = _get_run(arm, crop, chemistry)
        arrays = _arrays(summary)
        R = np.asarray(arrays["R"], dtype=float)
        run_n = int(summary["n"]) ** 3
        if R[0] <= 0:
            raise ValueError(f"Nonpositive R(t=0) for {summary['key']}")
        normalized = R / R[0]
        period_absolute = timemap.period_means(R / run_n, r_lit, dt=S23.DT)
        period_normalized = timemap.period_means(normalized, r_lit, dt=S23.DT)
        physiology = _physiology(summary, r_lit)
        density = float(chem23.CHEMISTRY[soil]["core_BD_g_cm3"])
        rows.append({
            "crop_id": crop,
            "chemistry": chemistry,
            "period_R_per_voxel": period_absolute.tolist(),
            "period_R_normalized_to_R0": period_normalized.tolist(),
            "shape_stage23": classify_shape_stage23(period_absolute),
            "shape_stage22": timemap.classify_period_shape(period_absolute),
            "cumulative_R_per_volume": physiology["cumulative_R_per_volume"],
            "cumulative_R_per_mass": physiology["cumulative_R_per_mass"],
            "P4_rate_per_mass_per_day": float(period_absolute[3] * r_lit / density),
            "qCO2_days10_15": physiology["qCO2_per_day"],
            "CUE": physiology["CUE"],
            "alive_primary": bool(summary["alive_primary"]),
        })
    return rows


def _summarize_group(rows: list[dict]) -> dict:
    absolute = np.asarray([row["period_R_per_voxel"] for row in rows], dtype=float)
    normalized = np.asarray(
        [row["period_R_normalized_to_R0"] for row in rows], dtype=float
    )
    mean_absolute = absolute.mean(axis=0)
    mean_normalized = normalized.mean(axis=0)
    stage23_shapes = {}
    stage22_shapes = {}
    for row in rows:
        stage23_shapes[row["shape_stage23"]] = stage23_shapes.get(row["shape_stage23"], 0) + 1
        stage22_shapes[row["shape_stage22"]] = stage22_shapes.get(row["shape_stage22"], 0) + 1
    return {
        "n_crops": len(rows),
        "period_R_per_voxel_mean_sd": [
            _mean_sd(absolute[:, index]) for index in range(4)
        ],
        "period_R_normalized_to_R0_mean_sd": [
            _mean_sd(normalized[:, index]) for index in range(4)
        ],
        "shape_stage23_mean_curve": classify_shape_stage23(mean_absolute),
        "shape_stage22_mean_curve": timemap.classify_period_shape(mean_absolute),
        "crop_shape_stage23_counts": stage23_shapes,
        "crop_shape_stage22_counts": stage22_shapes,
        "cumulative_R_per_volume": _mean_sd(
            [row["cumulative_R_per_volume"] for row in rows]
        ),
        "cumulative_R_per_mass": _mean_sd(
            [row["cumulative_R_per_mass"] for row in rows]
        ),
        "P4_rate_per_mass_per_day": _mean_sd(
            [row["P4_rate_per_mass_per_day"] for row in rows]
        ),
        "qCO2_days10_15": _mean_sd([row["qCO2_days10_15"] for row in rows]),
        "CUE": _mean_sd([row["CUE"] for row in rows]),
        "alive_primary_fraction": float(np.mean([row["alive_primary"] for row in rows])),
        "crop_results": rows,
    }


def _time_groups(arm: str, r_lit: float) -> dict:
    groups = {}
    for soil in ("S", "V", "L"):
        rows = _physiology_rows(arm, soil, r_lit, SOIL_CROPS[soil])
        groups[soil] = _summarize_group(rows)
    for specimen in S23.SCAN_IDS:
        soil = S23.SOIL_FOR_SCAN[specimen]
        rows = _physiology_rows(arm, soil, r_lit, SPECIMEN_CROPS[specimen])
        groups[specimen] = _summarize_group(rows)
    return groups


def _write_period_figure(time_results: dict, target: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 2, figsize=(12, 12), sharex=True, sharey=True)
    colors = ("#4e79a7", "#f28e2b", "#59a14f", "#e15759")
    for row_index, mapping in enumerate(("low", "primary", "high")):
        for col_index, arm in enumerate(("A", "B")):
            ax = axes[row_index, col_index]
            groups = time_results[arm][mapping]
            x = np.arange(3)
            width = 0.18
            for period_index, period in enumerate(PERIOD_NAMES):
                values = [
                    groups[soil]["period_R_normalized_to_R0_mean_sd"][period_index]["mean"]
                    for soil in ("S", "V", "L")
                ]
                errors = [
                    groups[soil]["period_R_normalized_to_R0_mean_sd"][period_index]["sd"]
                    for soil in ("S", "V", "L")
                ]
                positions = x + (period_index - 1.5) * width
                ax.bar(
                    positions,
                    values,
                    width,
                    yerr=errors,
                    color=colors[period_index],
                    alpha=0.8,
                    label=period if row_index == 0 and col_index == 0 else None,
                    capsize=2,
                )
            ax.set_xticks(x, [SOIL_NAMES[soil] for soil in ("S", "V", "L")])
            ax.set_title(f"{mapping} r_lit, Arm {arm}")
            ax.grid(axis="y", alpha=0.2)
    for ax in axes[-1, :]:
        ax.set_xlabel("Structure soil")
    for ax in axes[:, 0]:
        ax.set_ylabel("Period mean R / R(t=0)")
    axes[0, 0].legend(title="Period", fontsize=8)
    fig.tight_layout()
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite Stage 23 figure: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=140)
    plt.close(fig)


def step6() -> dict:
    target = S23.OUT / "step6.json"
    if target.exists():
        raise FileExistsError(f"Stage 23 Step 6 output already exists: {target}")
    result = {"stage": 23, "step": 6, "arms": {}}
    for arm in ("A", "B"):
        arm_result = {}
        for mapping, r_lit in S23.R_LIT_DAYS.items():
            arm_result[mapping] = _time_groups(arm, r_lit)
        result["arms"][arm] = arm_result
    rank_tables = {}
    for mapping in S23.R_LIT_DAYS:
        groups = result["arms"]["B"][mapping]
        rank_tables[mapping] = {}
        for metric in (
            "cumulative_R_per_mass",
            "P4_rate_per_mass_per_day",
            "qCO2_days10_15",
        ):
            ranked = sorted(
                ("S", "V", "L"),
                key=lambda soil: groups[soil][metric]["mean"],
                reverse=True,
            )
            rank_tables[mapping][metric] = [
                {
                    "rank": index + 1,
                    "soil": SOIL_NAMES[soil],
                    "mean": groups[soil][metric]["mean"],
                    "sd": groups[soil][metric]["sd"],
                }
                for index, soil in enumerate(ranked)
            ]
    result["arm_b_rank_tables"] = rank_tables
    result["figure"] = str(S23.OUT / "step6_periods.png")
    _write_period_figure(result["arms"], Path(result["figure"]))
    _write_json(target, result)
    return result


def main(argv: list[str] | None = None) -> None:
    args = list(sys.argv[1:] if argv is None else argv)
    steps = {"step3": step3, "step4": step4, "step5": step5, "step6": step6}
    if len(args) != 1 or args[0] not in steps:
        raise SystemExit("Usage: python stage23_analysis.py step3|step4|step5|step6")
    steps[args[0]]()


if __name__ == "__main__":
    main()
