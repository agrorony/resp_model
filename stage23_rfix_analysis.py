"""Reproduce the Stage 23 gates and analyze corrected respiration traces."""
from __future__ import annotations

import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import stage23_analysis as A
import model3d as M
import packing3d as P
import stage24_run as S24
import stage23_run as S23
import stage22_scan_structure_run as S22
import stage12_run as S12
import timemap
from model import load_yaml


OUT = Path("results/stage23_rfix")
GATE_DIR = OUT / "gates"
GATE_RUNS = (
    ("A", "SA_c0"),
    ("A", "VA_c0"),
    ("A", "LC_c0"),
    ("B", "SA_c0"),
    ("B", "VA_c0"),
    ("B", "LC_c0"),
)
STATE_FIELDS = (
    "B",
    "S_total",
    "uptake_dt",
    "fuel",
    "maint_R",
    "diff_wet",
    "diff_dry",
    "ex_wet",
    "ex_dry",
    "clip_added",
    "source_dt",
    "iters",
    "S_final",
    "B_final",
)
_CORRECTED_CACHE: dict[str, dict[str, np.ndarray]] = {}
_CORRECTED_SUMMARY_CACHE: dict[str, dict] = {}
_BIOLOGY: dict | None = None
_SOURCE_OUT = S23.OUT
_SOURCE_GET_RUN = None
_STAGE22_ARRAY_CACHE: dict[str, dict[str, np.ndarray]] = {}


def _read_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as file:
        return json.load(file)


def _write_json(path: Path, value: dict) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite growth-respiration gate output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def _max_relative_difference(actual: np.ndarray, reference: np.ndarray) -> float:
    actual = np.asarray(actual, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if actual.shape != reference.shape:
        raise ValueError(f"Array shapes differ: {actual.shape} versus {reference.shape}")
    scale = max(float(np.max(np.abs(reference))), 1e-300)
    denominator = np.maximum(np.abs(reference), scale * 1e-12)
    return float(np.max(np.abs(actual - reference) / denominator))


def _save_output(path: Path, output: dict) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite fresh gate trace: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        **{
            name: P.to_np(value)
            for name, value in output.items()
            if name in M.STAT_NAMES
        },
    )


def run_reproduction_gates() -> dict:
    GATE_DIR.mkdir(parents=True, exist_ok=False)
    started = datetime.now(timezone.utc).isoformat()
    xp = P.get_xp()
    biology = load_yaml("configs/biology.yaml")
    stage23 = _read_json(S23.OUT / "step2.json")
    checks = {}

    for arm, crop_id in GATE_RUNS:
        summary = next(
            row
            for row in stage23["runs"].values()
            if row["arm"] == arm
            and row["crop_id"] == crop_id
        )
        chemistry = summary["chemistry"]
        key = summary["key"]
        structure_id = summary["structure_id"]
        structure = S23._load_structure(structure_id, xp)
        hp = M.load_inherited()
        hp["B0_density"] = float(summary["B0_density"])
        stored_path = S23.RUNS / f"{key}.npz"
        with np.load(stored_path, allow_pickle=False) as reference_file:
            reference = {name: reference_file[name] for name in reference_file.files}
        required = {"R", "B", "growth_R", "maint_R", "S_total", "uptake_dt"}
        missing = sorted(required - reference.keys())
        if missing:
            raise KeyError(f"Stored Stage 23 trace {key} lacks required arrays: {missing}")

        outputs = {}
        run_checks = {}
        for mode, legacy in (("legacy", True), ("corrected", False)):
            output, run_summary = M.run_structure(
                structure,
                hp,
                biology,
                xp,
                float(summary["K_density"]),
                float(summary["OM_density"]),
                T=int(summary["T"]),
                legacy_growth_resp=legacy,
            )
            outputs[mode] = output
            _save_output(GATE_DIR / f"{key}.{mode}.npz", output)
            run_checks[mode] = {
                "mass_balance_rel": float(run_summary["mass_balance_rel"]),
                "closure_growth": float(output["closure_growth"]),
                "closure_growth_rel": float(output["closure_growth_rel"]),
            }

        legacy = outputs["legacy"]
        corrected = outputs["corrected"]
        legacy_differences = {
            name: _max_relative_difference(P.to_np(legacy[name]), value)
            for name, value in reference.items()
            if name in M.STAT_NAMES
        }
        max_legacy_difference = max(legacy_differences.values())
        state_differences = {
            name: _max_relative_difference(
                P.to_np(corrected[name]), P.to_np(legacy[name])
            )
            for name in STATE_FIELDS
        }
        expected_corrected_R = (
            P.to_np(legacy["growth_R"]) / float(biology["Y"])
            + P.to_np(legacy["maint_R"])
        )
        corrected_formula_difference = _max_relative_difference(
            P.to_np(corrected["R"]), expected_corrected_R
        )
        checks[key] = {
            "arm": arm,
            "crop_id": crop_id,
            "chemistry": chemistry,
            "stored_trace": str(stored_path),
            "fresh_trace_mtime_utc": datetime.fromtimestamp(
                (GATE_DIR / f"{key}.corrected.npz").stat().st_mtime,
                tz=timezone.utc,
            ).isoformat(),
            "legacy_max_relative_difference_by_array": legacy_differences,
            "legacy_max_relative_difference": max_legacy_difference,
            "corrected_vs_legacy_state_max_relative_difference_by_array": state_differences,
            "corrected_formula_max_relative_difference": corrected_formula_difference,
            "legacy": run_checks["legacy"],
            "corrected": run_checks["corrected"],
            "pass": (
                max_legacy_difference < 1e-8
                and max(state_differences.values()) < 1e-12
                and corrected_formula_difference < 1e-10
                and run_checks["corrected"]["closure_growth_rel"] < 1e-10
                and run_checks["legacy"]["mass_balance_rel"] <= 1e-10
                and run_checks["corrected"]["mass_balance_rel"] <= 1e-10
            ),
        }
        if not checks[key]["pass"]:
            _write_json(
                GATE_DIR / "gate.json",
                {
                    "started_utc": started,
                    "completed_utc": datetime.now(timezone.utc).isoformat(),
                    "passed": False,
                    "runs": checks,
                },
            )
            raise RuntimeError(f"Stage 23 respiration reproduction gate failed: {key}")
        print(
            f"PASS {key}: legacy={max_legacy_difference:.3g}, "
            f"state={max(state_differences.values()):.3g}, "
            f"formula={corrected_formula_difference:.3g}, "
            f"closure={run_checks['corrected']['closure_growth_rel']:.3g}, "
            f"mass={run_checks['corrected']['mass_balance_rel']:.3g}",
            flush=True,
        )

    result = {
        "stage": 23,
        "step": "growth_respiration_reproduction_gate",
        "started_utc": started,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "fresh_gpu_runs": True,
        "run_count": len(checks),
        "passed": len(checks) == 6 and all(row["pass"] for row in checks.values()),
        "runs": checks,
    }
    _write_json(GATE_DIR / "gate.json", result)
    return result


def _source_arrays(summary: dict) -> dict[str, np.ndarray]:
    path = S23.RUNS / f"{summary['key']}.npz"
    with np.load(path, allow_pickle=False) as data:
        arrays = {name: data[name] for name in data.files}
    required = {"R", "growth_R", "maint_R"}
    missing = sorted(required - arrays.keys())
    if missing:
        raise KeyError(f"Stored Stage 23 trace {summary['key']} lacks required arrays: {missing}")
    return arrays


def _corrected_arrays(summary: dict) -> dict[str, np.ndarray]:
    key = summary["key"]
    if key not in _CORRECTED_CACHE:
        if _BIOLOGY is None:
            raise RuntimeError("Stage 23 biology configuration is not initialized")
        arrays = _source_arrays(summary)
        legacy_growth = np.asarray(arrays["growth_R"], dtype=float)
        growth_R = legacy_growth / float(_BIOLOGY["Y"])
        arrays["R"] = growth_R + np.asarray(arrays["maint_R"], dtype=float)
        arrays["growth_R"] = growth_R
        _CORRECTED_CACHE[key] = arrays
    return _CORRECTED_CACHE[key]


def _corrected_summary(summary: dict) -> dict:
    return _summary_from_arrays(summary, _corrected_arrays(summary))


def _summary_from_arrays(summary: dict, arrays: dict[str, np.ndarray]) -> dict:
    key = summary["key"]
    if key in _CORRECTED_SUMMARY_CACHE:
        return _CORRECTED_SUMMARY_CACHE[key]
    if _BIOLOGY is None:
        raise RuntimeError("Growth-respiration biology configuration is not initialized")
    corrected = dict(summary)
    respiration = np.asarray(arrays["R"], dtype=float)
    peak = int(np.argmax(respiration))
    floor = float(_BIOLOGY["m0"]) * float(summary["B0_density"]) * int(summary["n"]) ** 3
    growth_R = float(arrays["growth_R"][peak])
    maintenance_R = float(arrays["maint_R"][peak])
    corrected.update({
        "R_peak": float(respiration[peak]),
        "R_peak_over_floor": float(respiration[peak] / floor),
        "growth_share_at_peak": growth_R / (growth_R + maintenance_R),
        "peak_frac": peak / (len(respiration) - 1),
        "R_end_over_peak": float(respiration[-1] / respiration[peak]),
        "peak_step": peak,
    })
    _CORRECTED_SUMMARY_CACHE[key] = corrected
    return corrected


def _persist_corrected_stage23_traces(stage23: dict) -> None:
    target_dir = OUT / "runs"
    target_dir.mkdir(parents=True, exist_ok=False)
    for summary in stage23["runs"].values():
        key = summary["key"]
        arrays = _corrected_arrays(summary)
        source_path = S23.RUNS / f"{key}.npz"
        source_mtime = datetime.fromtimestamp(source_path.stat().st_mtime, timezone.utc)
        corrected_path = target_dir / f"{key}.npz"
        if corrected_path.exists():
            raise FileExistsError(f"Refusing to overwrite corrected Stage 23 trace: {corrected_path}")
        np.savez_compressed(
            corrected_path,
            R=arrays["R"],
            growth_R=arrays["growth_R"],
        )
        _write_json(
            target_dir / f"{key}.json",
            {
                **_corrected_summary(summary),
                "respiration_formula": "R = growth_R_legacy / Y + maint_R",
                "source_trace": str(source_path),
                "source_trace_mtime_utc": source_mtime.isoformat(),
            },
        )


def _install_stage23_analysis_overlay(stage23_runs: dict) -> None:
    global _SOURCE_GET_RUN
    _, run_index = A._load_runs()
    A._RUNS = stage23_runs["runs"]
    A._RUN_INDEX = run_index
    A._ARRAY_CACHE = {}
    A._arrays = _corrected_arrays
    _SOURCE_GET_RUN = A._get_run

    def corrected_get_run(arm: str, crop_id: str, chemistry: str) -> dict:
        return _corrected_summary(_SOURCE_GET_RUN(arm, crop_id, chemistry))

    A._get_run = corrected_get_run
    S23.OUT = OUT


def _mean_sd(values) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    return float(values.mean()), float(values.std(ddof=1)) if values.size > 1 else 0.0


def _soil_summary(arm: str, soil: str, corrected: bool) -> dict:
    if _BIOLOGY is None:
        raise RuntimeError("Stage 23 biology configuration is not initialized")
    crops = A.SOIL_CROPS[soil]
    r_lit = S23.R_LIT_DAYS["primary"]
    rows = []
    for crop in crops:
        if _SOURCE_GET_RUN is None:
            raise RuntimeError("Stage 23 source run lookup is not initialized")
        source_summary = _SOURCE_GET_RUN(
            arm, crop, "chem_eq" if arm == "A" else f"chem_{soil}"
        )
        summary = _corrected_summary(source_summary) if corrected else source_summary
        arrays = _corrected_arrays(source_summary) if corrected else _source_arrays(source_summary)
        respiration = np.asarray(arrays["R"], dtype=float)
        growth_R = np.asarray(arrays["growth_R"], dtype=float)
        maintenance_R = np.asarray(arrays["maint_R"], dtype=float)
        peak = int(np.argmax(respiration))
        run_n = int(summary["n"]) ** 3
        physiology = timemap.physiological_outputs(
            respiration,
            arrays["B"],
            arrays["uptake_dt"],
            float(summary["OM_density"]) * run_n,
            r_lit,
            dt=S23.DT,
        )
        volume = physiology["cumulative_respiration_model_units"] / run_n
        density = float(A.chem23.CHEMISTRY[soil]["core_BD_g_cm3"])
        rows.append({
            "R_peak_over_floor": float(summary["R_peak_over_floor"]),
            "peak_frac": float(summary["peak_frac"]),
            "growth_share_at_peak": float(
                growth_R[peak] / (growth_R[peak] + maintenance_R[peak])
            ),
            "CUE": float(physiology["CUE"]),
            "qCO2": float(physiology["qCO2_per_day"]),
            "cumulative_R_per_volume": float(volume),
            "cumulative_R_per_mass": float(volume / density),
            "R_end_over_peak": float(summary["R_end_over_peak"]),
        })
    return {
        metric: _mean_sd([row[metric] for row in rows])
        for metric in rows[0]
    }


def _rank_codes(rows: list[dict]) -> list[str]:
    soil_codes = {"Sand": "S", "Vertisol": "V", "Loess": "L"}
    return [soil_codes[row["soil"]] for row in rows]


def _write_stage23_report(old: dict, new: dict) -> None:
    old3, old4, old6 = old["step3"], old["step4"], old["step6"]
    new3, new4, new6 = new["step3"], new["step4"], new["step6"]
    lines = [
        "# Stage 23 growth-respiration correction",
        "",
        "Legacy values are from the committed Stage 23 outputs; corrected values use `growth_R_legacy / Y + maint_R`.",
        "Only the six registered Arm A/B × SA_c0/VA_c0/LC_c0 runs were freshly rerun as gates. The remaining 125 traces were post-hoc corrected; no Stage 23 state was rerun.",
        "",
        "## Fresh reproduction and consistency gates",
        "",
        "See [`gates/gate.json`](../../results/stage23_rfix/gates/gate.json) and the paired fresh traces in `results/stage23_rfix/gates/`.",
        "",
        "## Arm A distinctness and peak order",
        "",
        "| Pair | Distance 0–1500, legacy → corrected | Distance 0–3000, legacy → corrected | Criteria, legacy → corrected | Verdict, legacy → corrected |",
        "|---|---:|---:|---|---|",
    ]
    for pair in ("S-V", "S-L", "V-L"):
        before = old3["arm_a"]["pairs"][pair]
        after = new3["arm_a"]["pairs"][pair]
        lines.append(
            f"| {pair} | {before['by_steps']['1500']['soil_mean_distance']:.6g} → {after['by_steps']['1500']['soil_mean_distance']:.6g} "
            f"| {before['by_steps']['3000']['soil_mean_distance']:.6g} → {after['by_steps']['3000']['soil_mean_distance']:.6g} "
            f"| `{before['criteria']}` → `{after['criteria']}` | {before['distinct']} → {after['distinct']} |"
        )
    height_before, height_after = _arm_peak_orders(old, new, "height")
    time_before, time_after = _arm_peak_orders(old, new, "time")
    lines.extend([
        "",
        f"- Peak-height order: **{' > '.join(height_before)} → {' > '.join(height_after)}**.",
        f"- Peak-time order (earliest to latest): **{' < '.join(time_before)} → {' < '.join(time_after)}**.",
        "",
        "Corrected Stage 22 scan-structure comparison for the shared Sand–Vertisol Arm A pair:",
    ])
    for horizon, comparison in new3["stage22_sand_vertisol_comparison"].items():
        lines.append(
            f"- 0–{horizon}: Stage 22 {comparison['stage22']['soil_mean_distance']:.6g}, "
            f"Stage 23 {comparison['stage23_arm_a']['soil_mean_distance']:.6g}; "
            f"equal at 4 decimals: {comparison['equal_at_4dp']}."
        )
    lines.extend([
        "",
        "## Arm B distinctness and rankings",
        "",
        "| Pair | Distance 0–1500, legacy → corrected | Distance 0–3000, legacy → corrected | Verdict, legacy → corrected |",
        "|---|---:|---:|---|",
    ])
    for pair in ("S-V", "S-L", "V-L"):
        before = old4["arm_b_distinctness"]["pairs"][pair]
        after = new4["arm_b_distinctness"]["pairs"][pair]
        lines.append(
            f"| {pair} | {before['by_steps']['1500']['soil_mean_distance']:.6g} → {after['by_steps']['1500']['soil_mean_distance']:.6g} "
            f"| {before['by_steps']['3000']['soil_mean_distance']:.6g} → {after['by_steps']['3000']['soil_mean_distance']:.6g} "
            f"| {before['distinct']} → {after['distinct']} |"
        )
    lines.extend([
        "",
        "| Mapping | Metric | Legacy rank | Corrected rank |",
        "|---|---|---|---|",
    ])
    for mapping in ("low", "primary", "high"):
        for metric in ("cumulative_R_per_mass", "P4_rate_per_mass_per_day", "qCO2_days10_15"):
            before = _rank_codes(old6["arm_b_rank_tables"][mapping][metric])
            after = _rank_codes(new6["arm_b_rank_tables"][mapping][metric])
            lines.append(f"| {mapping} | {metric} | {' > '.join(before)} | {' > '.join(after)} |")

    lines.extend([
        "",
        "## Shape classes",
        "",
        "| Arm | Mapping | Soil | Legacy | Corrected |",
        "|---|---|---|---|---|",
    ])
    for arm in ("A", "B"):
        for mapping in ("low", "primary", "high"):
            for soil in ("S", "V", "L"):
                before = old6["arms"][arm][mapping][soil]["shape_stage23_mean_curve"]
                after = new6["arms"][arm][mapping][soil]["shape_stage23_mean_curve"]
                lines.append(f"| {arm} | {mapping} | {soil} | {before} | {after} |")

    lines.extend([
        "",
        "## Two-way decomposition shares",
        "",
        "| Output | Component | Legacy | Corrected |",
        "|---|---|---:|---:|",
    ])
    for metric, decomposition in old4["two_way_decomposition"].items():
        updated = new4["two_way_decomposition"][metric]
        for component, share in decomposition["share_of_total"].items():
            lines.append(
                f"| {metric} | {component} | {share:.6g} | {updated['share_of_total'][component]:.6g} |"
            )

    lines.extend([
        "",
        "## Soil summary, Arm A/B primary mapping",
        "",
        "| Arm | Soil | Metric | Legacy mean ± SD | Corrected mean ± SD |",
        "|---|---|---|---:|---:|",
    ])
    soil_summary = {}
    for arm in ("A", "B"):
        for soil in ("S", "V", "L"):
            before = _soil_summary(arm, soil, False)
            after = _soil_summary(arm, soil, True)
            soil_summary[(arm, soil)] = {"legacy": before, "corrected": after}
            for metric in (
                "R_peak_over_floor", "peak_frac", "growth_share_at_peak",
                "R_end_over_peak", "CUE",
                "qCO2", "cumulative_R_per_volume", "cumulative_R_per_mass",
            ):
                old_mean, old_sd = before[metric]
                new_mean, new_sd = after[metric]
                lines.append(
                    f"| {arm} | {soil} | {metric} | {old_mean:.6g} ± {old_sd:.3g} | {new_mean:.6g} ± {new_sd:.3g} |"
                )

    lines.extend([
        "",
        "## Complete Arm B/C structure × chemistry outputs",
        "",
        "| Structure × chemistry | Metric | Legacy mean ± SD | Corrected mean ± SD |",
        "|---|---|---:|---:|",
    ])
    for cell, before_cell in old4["cross_cell_means_and_sd"].items():
        after_cell = new4["cross_cell_means_and_sd"][cell]
        for metric, before in before_cell["metrics"].items():
            after = after_cell["metrics"][metric]
            lines.append(
                f"| {cell} | {metric} | {before['mean']:.6g} ± {before['sd']:.3g} "
                f"| {after['mean']:.6g} ± {after['sd']:.3g} |"
            )
    lines.extend([
        "",
        "## Curve-level structure/chemistry distances",
        "",
        "| Horizon | Comparison | Legacy mean ± SD | Corrected mean ± SD |",
        "|---|---|---:|---:|",
    ])
    for horizon in ("1500", "3000"):
        for metric in (
            "same_structure_different_chemistry",
            "same_chemistry_different_structure",
        ):
            before = old4["curve_level_distances"][horizon][metric]
            after = new4["curve_level_distances"][horizon][metric]
            lines.append(
                f"| 0–{horizon} | {metric} | {before['mean']:.6g} ± {before['sd']:.3g} "
                f"| {after['mean']:.6g} ± {after['sd']:.3g} |"
            )
    lines.extend([
        "",
        "## Arm D and Loess domain side checks",
        "",
        "| Check | Legacy | Corrected |",
        "|---|---:|---:|",
    ])
    for chemistry, row in old["step5"]["arm_d_sand_rule_check"].items():
        updated = new["step5"]["arm_d_sand_rule_check"][chemistry]
        lines.append(
            f"| Arm D {chemistry} mean peak/floor, solid-fines → pore-fill | "
            f"{row['mean_base_peak_over_floor']:.6g} → {row['mean_pore_fill_peak_over_floor']:.6g} | "
            f"{updated['mean_base_peak_over_floor']:.6g} → {updated['mean_pore_fill_peak_over_floor']:.6g} |"
        )
        lines.append(
            f"| Arm D {chemistry} mean peak fraction, solid-fines → pore-fill | "
            f"{row['mean_base_peak_frac']:.6g} → {row['mean_pore_fill_peak_frac']:.6g} | "
            f"{updated['mean_base_peak_frac']:.6g} → {updated['mean_pore_fill_peak_frac']:.6g} |"
        )
    for specimen, row in old["step5"]["loess_domain_size_check"].items():
        updated = new["step5"]["loess_domain_size_check"][specimen]
        lines.append(
            f"| Loess {specimen} domain peak/floor | {row['domain_R_peak_over_floor']:.6g} "
            f"| {updated['domain_R_peak_over_floor']:.6g} |"
        )
        lines.append(
            f"| Loess {specimen} domain peak fraction | {row['domain_peak_frac']:.6g} "
            f"| {updated['domain_peak_frac']:.6g} |"
        )

    lines.extend([
        "",
        "## P1–P4 means by arm, mapping, and soil",
        "",
        "The first table is absolute R per voxel; the second is normalized to R(0).",
        "",
        "| Arm | Mapping | Soil | Legacy P1 / P2 / P3 / P4 | Corrected P1 / P2 / P3 / P4 |",
        "|---|---|---|---|---|",
    ])
    for arm in ("A", "B"):
        for mapping in ("low", "primary", "high"):
            for soil in ("S", "V", "L"):
                before = old6["arms"][arm][mapping][soil]["period_R_per_voxel_mean_sd"]
                after = new6["arms"][arm][mapping][soil]["period_R_per_voxel_mean_sd"]
                old_values = " / ".join(f"{row['mean']:.6g}" for row in before)
                new_values = " / ".join(f"{row['mean']:.6g}" for row in after)
                lines.append(
                    f"| {arm} | {mapping} | {soil} | {old_values} | {new_values} |"
                )
    lines.extend([
        "",
        "| Arm | Mapping | Soil | Legacy normalized P1 / P2 / P3 / P4 | Corrected normalized P1 / P2 / P3 / P4 |",
        "|---|---|---|---|---|",
    ])
    for arm in ("A", "B"):
        for mapping in ("low", "primary", "high"):
            for soil in ("S", "V", "L"):
                before = old6["arms"][arm][mapping][soil]["period_R_normalized_to_R0_mean_sd"]
                after = new6["arms"][arm][mapping][soil]["period_R_normalized_to_R0_mean_sd"]
                old_values = " / ".join(f"{row['mean']:.6g}" for row in before)
                new_values = " / ".join(f"{row['mean']:.6g}" for row in after)
                lines.append(
                    f"| {arm} | {mapping} | {soil} | {old_values} | {new_values} |"
                )

    if any(
        values["corrected"]["CUE"][0] >= float(_BIOLOGY["Y"])
        for values in soil_summary.values()
    ):
        raise ArithmeticError("Corrected Stage 23 CUE must remain below yield")

    flips = _list_flips(old, new)
    lines.extend([
        "",
        "## What flipped",
        "",
        *([f"- {flip}" for flip in flips] if flips else ["- Nothing flipped among the registered distinctness verdicts, ranks, peak orders, or shape classes."]),
        "",
        "## Numbers to update in the appendix",
        "",
        "Replace the Stage 23 Step 2 and Step 4 peak/floor, growth-share, peak-fraction, end/peak, cumulative-R, qCO2, and decomposition numbers with the complete old-to-new tables above. Replace the Step 5 Arm D and domain peak/floor and peak-time numbers, the Step 6 P1–P4 means, shape classes, CUE, qCO2 and cumulative-R/rank tables, and the peak-order/curve-distance figures. Biomass, habitat, delivery, and substrate-uptake results are unchanged.",
        "",
    ])
    path = Path("docs/stage_results/STAGE23_RFIX.md")
    if path.exists():
        first_line = path.read_text(encoding="utf-8").splitlines()[:1]
        if first_line != ["# Stage 23 growth-respiration correction"]:
            raise FileExistsError(f"Refusing to overwrite an unrelated report: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def _arm_peak_orders(old: dict, new: dict, metric: str) -> tuple[list[str], list[str]]:
    def order(corrected: bool) -> list[str]:
        values = {}
        for soil, crops in A.SOIL_CROPS.items():
            summaries = []
            for crop in crops:
                if _SOURCE_GET_RUN is None:
                    raise RuntimeError("Stage 23 source run lookup is not initialized")
                source = _SOURCE_GET_RUN("A", crop, "chem_eq")
                summary = _corrected_summary(source) if corrected else source
                arrays = _corrected_arrays(source) if corrected else _source_arrays(source)
                if metric == "height":
                    value = float(summary["R_peak_over_floor"])
                else:
                    value = float(summary["peak_frac"])
                summaries.append(value)
            values[soil] = float(np.mean(summaries))
        return sorted(values, key=values.get, reverse=(metric == "height"))

    return order(False), order(True)


def _list_flips(old: dict, new: dict) -> list[str]:
    flips = []
    height_before, height_after = _arm_peak_orders(old, new, "height")
    time_before, time_after = _arm_peak_orders(old, new, "time")
    if height_before != height_after:
        flips.append(
            f"Arm A peak-height order: {' > '.join(height_before)} → {' > '.join(height_after)}."
        )
    if time_before != time_after:
        flips.append(
            f"Arm A peak-time order: {' < '.join(time_before)} → {' < '.join(time_after)}."
        )
    for arm, before_section, after_section in (
        ("A", old["step3"]["arm_a"], new["step3"]["arm_a"]),
        ("B", old["step4"]["arm_b_distinctness"], new["step4"]["arm_b_distinctness"]),
    ):
        for pair in ("S-V", "S-L", "V-L"):
            before = before_section["pairs"][pair]
            after = after_section["pairs"][pair]
            if before["distinct"] != after["distinct"]:
                dist = before["by_steps"]["1500"]["soil_mean_distance"]
                flips.append(
                    f"Arm {arm} {pair} distinctness: {before['distinct']} → {after['distinct']}; "
                    f"legacy 0–1500 distance={dist:.6g} (criterion threshold 0.3; noise margin requires review)."
                )
    for mapping in ("low", "primary", "high"):
        for metric in ("cumulative_R_per_mass", "P4_rate_per_mass_per_day", "qCO2_days10_15"):
            before = _rank_codes(old["step6"]["arm_b_rank_tables"][mapping][metric])
            after = _rank_codes(new["step6"]["arm_b_rank_tables"][mapping][metric])
            if before != after:
                old_rows = {
                    row["soil"][0]: row
                    for row in old["step6"]["arm_b_rank_tables"][mapping][metric]
                }
                first, second = next(
                    (a, b) for a, b in zip(before, before[1:]) if after.index(a) > after.index(b)
                )
                gap = old_rows[first]["mean"] - old_rows[second]["mean"]
                combined_sd = float(np.hypot(old_rows[first]["sd"], old_rows[second]["sd"]))
                flips.append(
                    f"Arm B {mapping} {metric} rank: {' > '.join(before)} → {' > '.join(after)}; "
                    f"legacy {first}−{second} mean gap={gap:.6g}, combined crop SD={combined_sd:.6g} "
                    f"({'within' if abs(gap) <= combined_sd else 'outside'} one combined SD)."
                )
        for arm in ("A", "B"):
            for soil in ("S", "V", "L"):
                before = old["step6"]["arms"][arm][mapping][soil]["shape_stage23_mean_curve"]
                after = new["step6"]["arms"][arm][mapping][soil]["shape_stage23_mean_curve"]
                if before != after:
                    old_periods = old["step6"]["arms"][arm][mapping][soil]["period_R_per_voxel_mean_sd"]
                    new_periods = new["step6"]["arms"][arm][mapping][soil]["period_R_per_voxel_mean_sd"]
                    p1_p2_gap = old_periods[1]["mean"] - old_periods[0]["mean"]
                    p1_p2_sd = float(np.hypot(old_periods[1]["sd"], old_periods[0]["sd"]))
                    flips.append(
                        f"Arm {arm} {mapping} {soil} shape: {before} → {after} "
                        f"(legacy P2−P1={p1_p2_gap:.6g}, combined crop SD={p1_p2_sd:.6g}; "
                        f"old P1/P2 ordering {'within' if abs(p1_p2_gap) <= p1_p2_sd else 'outside'} one SD; "
                        f"corrected P1/P2 means={new_periods[0]['mean']:.6g}/{new_periods[1]['mean']:.6g})."
                    )
    for arm in ("A", "B"):
        old_order = next(
            row for row in old["step5"]["sand_ordering_with_arm_d_substitution"]
            if row["arm"] == arm and row["sand_rule"] == "stage21_pore_fill"
        )
        old_base = next(
            row for row in old["step5"]["sand_ordering_with_arm_d_substitution"]
            if row["arm"] == arm and row["sand_rule"] == "stage22_solid_lt"
        )
        new_order = next(
            row for row in new["step5"]["sand_ordering_with_arm_d_substitution"]
            if row["arm"] == arm and row["sand_rule"] == "stage21_pore_fill"
        )
        new_base = next(
            row for row in new["step5"]["sand_ordering_with_arm_d_substitution"]
            if row["arm"] == arm and row["sand_rule"] == "stage22_solid_lt"
        )
        for horizon in ("1500", "3000"):
            before = old_order["by_steps"][horizon]["ordering_changed"]
            after = new_order["by_steps"][horizon]["ordering_changed"]
            if before != after:
                before_gap = abs(
                    old_order["by_steps"][horizon]["S_V_distance"]
                    - old_order["by_steps"][horizon]["S_L_distance"]
                )
                after_gap = abs(
                    new_order["by_steps"][horizon]["S_V_distance"]
                    - new_order["by_steps"][horizon]["S_L_distance"]
                )
                base_old_pair = old_base["by_steps"][horizon]["closer_pair"]
                base_new_pair = new_base["by_steps"][horizon]["closer_pair"]
                old_pair = old_order["by_steps"][horizon]["closer_pair"]
                new_pair = new_order["by_steps"][horizon]["closer_pair"]
                flips.append(
                    f"Arm {arm} Stage 21 pore-fill ordering at {horizon} steps: "
                    f"{before} → {after}; base closest pair {base_old_pair} → {base_new_pair}, "
                    f"pore-fill closest pair {old_pair} → {new_pair}; legacy S–V/S–L distance gap "
                    f"{before_gap:.6g}, corrected gap {after_gap:.6g}. No registered noise interval "
                    "was calculated for this single mean-curve ordering."
                )
    return flips


def run_stage23_analysis() -> dict:
    gate = _read_json(GATE_DIR / "gate.json")
    if not gate.get("passed") or not gate.get("fresh_gpu_runs"):
        raise RuntimeError("Corrected analysis requires a passing fresh Stage 23 GPU gate")
    global _BIOLOGY
    _BIOLOGY = load_yaml("configs/biology.yaml")
    stage23 = _read_json(_SOURCE_OUT / "step2.json")
    OUT.mkdir(parents=True, exist_ok=True)
    _persist_corrected_stage23_traces(stage23)
    _install_stage23_analysis_overlay(stage23)
    old = {
        f"step{step}": _read_json(_SOURCE_OUT / f"step{step}.json")
        for step in (3, 4, 5, 6)
    }
    new = {
        "step3": A.step3(),
        "step4": A.step4(),
        "step5": A.step5(),
        "step6": A.step6(),
    }
    _write_stage23_report(old, new)
    return new


def refresh_stage23_report() -> str:
    global _BIOLOGY
    _BIOLOGY = load_yaml("configs/biology.yaml")
    stage23 = _read_json(_SOURCE_OUT / "step2.json")
    _install_stage23_analysis_overlay(stage23)
    old = {
        f"step{step}": _read_json(_SOURCE_OUT / f"step{step}.json")
        for step in (3, 4, 5, 6)
    }
    new = {
        f"step{step}": _read_json(OUT / f"step{step}.json")
        for step in (3, 4, 5, 6)
    }
    _write_stage23_report(old, new)
    return "docs/stage_results/STAGE23_RFIX.md"


def refresh_stage23_stage22_comparison() -> dict:
    target_path = OUT / "step3.json"
    result = _read_json(target_path)
    stage22 = _read_json(
        Path("results/stage22_scan_structure_rfix/step3.json")
    )
    pair = result["arm_a"]["pairs"]["S-V"]
    comparison = {}
    for horizon in ("1500", "3000"):
        current = pair["by_steps"][horizon]
        prior = {
            "soil_mean_distance": stage22["soil_mean_distance"][horizon],
            "max_within": stage22["within_soil"][horizon]["max_all"],
            "min_cross": stage22["cross_soil"][horizon]["min"],
        }
        observed = {
            "soil_mean_distance": current["soil_mean_distance"],
            "max_within": current["max_within"]["distance"],
            "min_cross": current["min_cross"]["distance"],
        }
        comparison[horizon] = {
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
    result["stage22_sand_vertisol_comparison"] = comparison
    target_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return comparison


def run_stage24_posthoc_analysis() -> dict:
    global _BIOLOGY
    gate = _read_json(GATE_DIR / "gate.json")
    if not gate.get("passed"):
        raise RuntimeError("Stage 24 correction requires the passing Stage 23 gate")
    _BIOLOGY = load_yaml("configs/biology.yaml")
    stage23 = _read_json(_SOURCE_OUT / "step2.json")
    _install_stage23_analysis_overlay(stage23)
    source_out = Path("results/stage24")
    target_out = Path("results/stage24_rfix")
    target_out.mkdir(parents=True, exist_ok=False)
    run_manifest = _read_json(source_out / "runs.json")
    geometry = _read_json(source_out / "geometry.json")
    _write_json(target_out / "runs.json", run_manifest)
    _write_json(target_out / "geometry.json", geometry)
    _write_json(
        target_out / "posthoc_manifest.json",
        {
            "source": str(source_out),
            "corrected_formula": "R = growth_R_legacy / Y + maint_R",
            "stage24_variant_runs_posthoc_corrected": len(run_manifest["runs"]),
            "fresh_stage24_runs": 0,
            "stage23_corrected_analysis": [
                str(OUT / f"step{step}.json") for step in (3, 4, 5, 6)
            ],
        },
    )

    variant_runs_dir = target_out / "runs"
    variant_runs_dir.mkdir()
    for key in run_manifest["runs"]:
        with np.load(S24.RUNS / f"{key}.npz", allow_pickle=False) as data:
            arrays = {name: data[name] for name in data.files}
        required = {"growth_R", "maint_R"}
        missing = sorted(required - arrays.keys())
        if missing:
            raise KeyError(f"Stored Stage 24 trace {key} lacks required arrays: {missing}")
        legacy_growth = np.asarray(arrays["growth_R"], dtype=float)
        corrected_growth = legacy_growth / float(_BIOLOGY["Y"])
        corrected_R = corrected_growth + np.asarray(arrays["maint_R"], dtype=float)
        npz_path = variant_runs_dir / f"{key}.npz"
        if npz_path.exists():
            raise FileExistsError(f"Refusing to overwrite Stage 24 corrected trace: {npz_path}")
        np.savez_compressed(npz_path, R=corrected_R, growth_R=corrected_growth)

    source_read_json = S24._read_json
    source_get_run = S24._get_run
    source_run_arrays = S24._run_arrays

    def corrected_read_json(path: Path) -> dict:
        if path == Path("results/stage23/step3.json"):
            return _read_json(OUT / "step3.json")
        if path == Path("results/stage23/step6.json"):
            return _read_json(OUT / "step6.json")
        return source_read_json(path)

    def corrected_run_arrays(summary: dict) -> dict[str, np.ndarray]:
        if summary.get("variant") is None:
            return _corrected_arrays(summary)
        arrays = source_run_arrays(summary)
        legacy_growth = np.asarray(arrays["growth_R"], dtype=float)
        corrected_growth = legacy_growth / float(_BIOLOGY["Y"])
        return {
            **arrays,
            "R": corrected_growth + np.asarray(arrays["maint_R"], dtype=float),
            "growth_R": corrected_growth,
        }

    def corrected_get_run(variant: str, arm: str, crop_id: str) -> dict:
        summary = source_get_run(variant, arm, crop_id)
        return _summary_from_arrays(summary, corrected_run_arrays(summary))

    S24._read_json = corrected_read_json
    S24._run_arrays = corrected_run_arrays
    S24._get_run = corrected_get_run
    S24.OUT = target_out
    result = S24.step4()
    return result


def _stage22_yield(key: str) -> float:
    config_path = S22.JOB_CONFIGS / f"{key}.json"
    config = _read_json(config_path)
    biology = load_yaml("configs/biology.yaml")
    biology.update(config.get("biology_overrides", {}))
    return float(biology["Y"])


def _stage22_arrays(key: str) -> dict[str, np.ndarray]:
    if key not in _STAGE22_ARRAY_CACHE:
        source = S22.RUNS / f"{key}.npz"
        with np.load(source, allow_pickle=False) as data:
            arrays = {name: data[name] for name in data.files}
        required = {"growth_R", "maint_R"}
        missing = sorted(required - arrays.keys())
        if missing:
            raise KeyError(f"Stored Stage 22 trace {key} lacks required arrays: {missing}")
        corrected_growth = (
            np.asarray(arrays["growth_R"], dtype=float) / _stage22_yield(key)
        )
        arrays["R"] = corrected_growth + np.asarray(arrays["maint_R"], dtype=float)
        arrays["growth_R"] = corrected_growth
        _STAGE22_ARRAY_CACHE[key] = arrays
    return _STAGE22_ARRAY_CACHE[key]


def _stage22_run(key: str) -> dict:
    summary = _read_json(S22.RUNS / f"{key}.json")
    arrays = _stage22_arrays(key)
    summary["R_peak"] = float(np.max(arrays["R"]))
    summary["peak_step"] = int(np.argmax(arrays["R"]))
    summary["peak_frac"] = summary["peak_step"] / (len(arrays["R"]) - 1)
    summary["R_end_over_peak"] = float(arrays["R"][-1] / summary["R_peak"])
    return {**summary, **arrays}


def run_stage22_sensitivity_posthoc() -> dict:
    old = _read_json(Path("results/stage22_scan_structure/step5.json"))
    representatives = old["representative_crops"]
    representative_runs = {}
    for soil, record in representatives.items():
        key = S22._run_key(record["structure_id"])
        representative_runs[soil] = _stage22_run(key)
    base_orders = S22._peak_order(
        representative_runs["sand"], representative_runs["vertisol"]
    )
    result = {
        "representative_crops": representatives,
        "perturbations": {},
        "n_runs": int(old["n_runs"]),
        "active_step2": "step2.json",
        "method": "Post-hoc correction of stored traces; no reruns.",
    }
    for parameter, factors in old["perturbations"].items():
        result["perturbations"][parameter] = {}
        for factor, before in factors.items():
            sand_key = before["run_keys"]["sand"]
            vertisol_key = before["run_keys"]["vertisol"]
            sand_run = _stage22_run(sand_key)
            vertisol_run = _stage22_run(vertisol_key)
            order = S22._peak_order(sand_run, vertisol_run)
            result["perturbations"][parameter][factor] = {
                "base_value": before["base_value"],
                "perturbed_value": before["perturbed_value"],
                "alive_primary_both": before["alive_primary_both"],
                "alive_primary": before["alive_primary"],
                "distance_01500": S22._distance(
                    sand_run["R"], vertisol_run["R"], 1500
                ),
                "peak_amplitude_order": {
                    "base": base_orders[0],
                    "perturbed": order[0],
                    "unchanged": base_orders[0] == order[0],
                },
                "peak_time_order": {
                    "base": base_orders[1],
                    "perturbed": order[1],
                    "unchanged": base_orders[1] == order[1],
                },
                "run_keys": before["run_keys"],
                "yield_values": {
                    "sand": _stage22_yield(sand_key),
                    "vertisol": _stage22_yield(vertisol_key),
                },
            }
    target = Path("results/stage22_scan_structure_rfix/step5_posthoc.json")
    _write_json(target, result)
    return result


def run_figure_a1_rfix() -> dict:
    target = Path("results/appendix_figure_rfix")
    target.mkdir(parents=True, exist_ok=False)
    biology = load_yaml("configs/biology.yaml")
    mapping = load_yaml("configs/mapping.yaml")
    mapping["om_mode"] = S12.dual_porosity.OM_MODE_COATING
    rows_by_soil = {}
    curves_by_soil = {}
    started = time.monotonic()
    for soil in S12.SOILS:
        rows = []
        curves = []
        seeds = [
            S12.SEED_BASE[soil] + index
            for index in range(S12.M_SEEDS)
        ]
        for index, seed in enumerate(seeds):
            row, curve = S12.run_one(soil, seed, biology, mapping)
            rows.append(row)
            curves.append(curve)
            np.save(target / f"curve_{soil}_{index:02}.npy", curve)
            if (index + 1) % 10 == 0:
                print(
                    f"{soil}: {index + 1}/{S12.M_SEEDS} seeds complete",
                    flush=True,
                )
        rows_by_soil[soil] = rows
        curves_by_soil[soil] = np.stack(curves)
        np.save(target / f"curves_{soil}.npy", curves_by_soil[soil])

    S12.write_seed_metrics_csv(rows_by_soil, str(target))
    S12.write_metric_summary_csv(rows_by_soil, str(target))
    significance = S12.write_significance_csv(rows_by_soil, str(target))
    source_script = Path("results/appendix_figure/plot_fig_A1.py").read_text(
        encoding="utf-8"
    )
    target_prefix = target.as_posix() + "/"
    output_script = source_script.replace(
        'D="results/appendix_figure/"', f'D="{target_prefix}"'
    ).replace(
        'open("results/stage12/stage12_significance.csv")',
        f'open("{target_prefix}stage12_significance.csv")',
    )
    if output_script == source_script:
        raise ValueError("Figure A1 source paths did not match the expected inputs")
    exec(compile(output_script, "plot_fig_A1.py", "exec"), {"__name__": "__main__"})

    changes = {}
    for soil in S12.SOILS:
        legacy = np.load(
            Path("results/appendix_figure") / f"curves_{soil}.npy",
            allow_pickle=False,
        )
        corrected = curves_by_soil[soil]
        if legacy.shape != corrected.shape:
            raise ValueError(
                f"Figure A1 {soil} ensemble shape changed: {legacy.shape} versus {corrected.shape}"
            )
        legacy_mean = legacy.mean(axis=0)
        corrected_mean = corrected.mean(axis=0)
        max_abs_change = float(np.max(np.abs(corrected_mean - legacy_mean)))
        changes[soil] = {
            "legacy_peak": float(legacy_mean.max()),
            "corrected_peak": float(corrected_mean.max()),
            "legacy_peak_fraction": float(np.argmax(legacy_mean) / (legacy_mean.size - 1)),
            "corrected_peak_fraction": float(
                np.argmax(corrected_mean) / (corrected_mean.size - 1)
            ),
            "max_absolute_mean_curve_change": max_abs_change,
            "max_change_fraction_of_legacy_peak": max_abs_change / float(legacy_mean.max()),
            "visible_change": max_abs_change > 0.05 * float(legacy_mean.max()),
        }
    peak_fraction_stats = next(row for row in significance if row["metric"] == "peak_frac")
    result = {
        "runs": 2 * S12.M_SEEDS,
        "fresh_runs": True,
        "elapsed_seconds": time.monotonic() - started,
        "source_script": "results/appendix_figure/plot_fig_A1.py",
        "source_script_changed_only_for_paths": True,
        "lab_files_read": False,
        "normalization": (
            "The Figure A1 curves are resampled in time but not divided by their initial R; "
            "the factor does not cancel. R(0) remains the maintenance floor."
        ),
        "panel_changes": changes,
        "corrected_peak_frac_significance": peak_fraction_stats,
    }
    _write_json(target / "rfix_manifest.json", result)
    return result


def calculate_starvation_loss() -> dict:
    biology = load_yaml("configs/biology.yaml")
    Y = float(biology["Y"])
    m_s = float(biology["m_s"])
    dt = float(S23.DT)
    step2 = _read_json(S23.OUT / "step2.json")
    per_soil = {}
    for soil, crop_prefix in (("S", "SA"), ("V", "VA"), ("L", "LC")):
        rows = []
        for summary in step2["runs"].values():
            if summary["arm"] != "A" or summary["soil_code"] != soil:
                continue
            arrays = _source_arrays(summary)
            biomass = np.asarray(arrays["B"], dtype=float)
            uptake_dt = np.asarray(arrays["uptake_dt"], dtype=float)
            maint_R = np.asarray(arrays["maint_R"], dtype=float)
            # The stored B trace ends before the final biomass update.
            starvation_partial = (
                Y * float(uptake_dt[:-1].sum())
                - dt * float(maint_R[:-1].sum())
                - float(biomass[-1] - biomass[0])
            )
            terminal_bound = m_s * float(biomass[-1]) * dt
            corrected_R = (
                np.asarray(arrays["growth_R"], dtype=float) / Y
                + maint_R
            )
            uptake_total = float(uptake_dt.sum())
            respiration_total = float(corrected_R.sum() * dt)
            rows.append({
                "key": summary["key"],
                "crop_id": summary["crop_id"],
                "starvation_partial": starvation_partial,
                "terminal_upper_bound": terminal_bound,
                "uptake_total": uptake_total,
                "corrected_respiration_total": respiration_total,
            })
        if not rows:
            raise ValueError(f"No stored Stage 23 Arm A traces for soil {soil}")
        total_partial = sum(row["starvation_partial"] for row in rows)
        total_bound = sum(row["terminal_upper_bound"] for row in rows)
        total_uptake = sum(row["uptake_total"] for row in rows)
        total_respiration = sum(row["corrected_respiration_total"] for row in rows)

        fresh_key = next(
            key for key in _read_json(GATE_DIR / "gate.json")["runs"]
            if key.startswith(f"A-{crop_prefix}_c0-")
        )
        fresh_summary = next(
            row for row in step2["runs"].values() if row["key"] == fresh_key
        )
        with np.load(GATE_DIR / f"{fresh_key}.corrected.npz", allow_pickle=False) as fresh:
            fresh_partial = (
                Y * float(fresh["uptake_dt"][:-1].sum())
                - dt * float(fresh["maint_R"][:-1].sum())
                - float(fresh["B"][-1] - fresh["B"][0])
            )
            fresh_terminal_bound = m_s * float(fresh["B"][-1]) * dt
            fresh_uptake = float(fresh["uptake_dt"].sum())
            fresh_respiration = float(fresh["R"].sum() * dt)
            fresh_result = {
                "key": fresh_key,
            }
            if "B_final" in fresh.files:
                fresh_result.update({
                    "fraction_of_cumulative_uptake_range": [
                        fresh_partial / fresh_uptake,
                        (fresh_partial + fresh_terminal_bound) / fresh_uptake,
                    ],
                    "fraction_of_corrected_cumulative_R_range": [
                        fresh_partial / fresh_respiration,
                        (fresh_partial + fresh_terminal_bound) / fresh_respiration,
                    ],
                })
                fresh_growth = Y * float(fresh["uptake_dt"].sum())
                fresh_maintenance = dt * float(fresh["maint_R"].sum())
                initial_biomass = (
                    float(fresh_summary["B0_density"]) * int(fresh_summary["n"]) ** 3
                )
                final_biomass = float(np.asarray(fresh["B_final"]).sum())
                fresh_starvation = (
                    fresh_growth - fresh_maintenance - (final_biomass - initial_biomass)
                )
                fresh_terminal = fresh_starvation - fresh_partial
                fresh_result.update({
                    "starvation_loss": fresh_starvation,
                    "fraction_of_cumulative_uptake": fresh_starvation / fresh_uptake,
                    "fraction_of_corrected_cumulative_R": fresh_starvation / fresh_respiration,
                    "last_step_loss": fresh_terminal,
                    "last_step_upper_bound": fresh_terminal_bound,
                })
                if fresh_terminal < -1e-10 or fresh_terminal > fresh_terminal_bound + 1e-8:
                    raise ArithmeticError(
                        f"Fresh {soil} c0 starvation residual is outside its final-step bound"
                    )
            else:
                retained = _read_json(OUT / "starvation_loss.json")["soils"][soil]["fresh_c0_run"]
                fresh_starvation = float(retained["starvation_loss"])
                fresh_terminal = fresh_starvation - fresh_partial
                if fresh_terminal < -1e-10 or fresh_terminal > fresh_terminal_bound + 1e-8:
                    raise ArithmeticError(
                        f"Retained {soil} c0 starvation residual is outside its final-step bound"
                    )
                fresh_result.update({
                    "starvation_loss": fresh_starvation,
                    "fraction_of_cumulative_uptake": fresh_starvation / fresh_uptake,
                    "fraction_of_corrected_cumulative_R": fresh_starvation / fresh_respiration,
                    "last_step_loss": fresh_terminal,
                    "last_step_upper_bound": fresh_terminal_bound,
                })

        per_soil[soil] = {
            "stored_run_count": len(rows),
            "stored_trace_starvation_lower_bound": total_partial,
            "stored_trace_final_step_upper_bound": total_bound,
            "fraction_of_cumulative_uptake_range": [
                total_partial / total_uptake,
                (total_partial + total_bound) / total_uptake,
            ],
            "fraction_of_corrected_cumulative_R_range": [
                total_partial / total_respiration,
                (total_partial + total_bound) / total_respiration,
            ],
            "fresh_c0_run": fresh_result,
        }
        if total_partial < 0:
            raise ArithmeticError(f"Negative stored-trace starvation loss for soil {soil}")

    result = {
        "m_s": m_s,
        "method": (
            "Stored-trace loss is reconstructed from biomass change, Y*uptake, and maintenance "
            "for steps 0..T-2; the final-step loss is bounded by m_s*B[T-1]*dt. "
            "The six fresh gate runs include B_final, enabling exact full-window Arm A c0 estimates."
        ),
        "biology_changed": False,
        "decision_owner": "Rony",
        "soils": per_soil,
        "interpretation": (
            "Starved biomass is removed from dB but its carbon is neither emitted as CO2 nor returned to S."
        ),
    }
    return result


def run_starvation_loss_report() -> dict:
    result = calculate_starvation_loss()
    _write_json(OUT / "starvation_loss.json", result)
    return result


def verify_starvation_loss_report() -> dict:
    calculated = calculate_starvation_loss()
    saved = _read_json(OUT / "starvation_loss.json")
    if calculated != saved:
        raise ArithmeticError("Reconstructed starvation-loss results differ from the saved report")
    return {"matches_saved_report": True, "soils": len(calculated["soils"])}


def _write_stage22_report(old: dict, new: dict, sensitivity: dict) -> None:
    before, after = old["step3"], new["step3"]
    lines = [
        "# Stage 22 scan-structure respiration correction",
        "",
        "All 58 stored Stage 22 scan-structure traces were post-hoc corrected from `growth_R` and `maint_R`; no Stage 22 or earlier simulation was rerun. The Stage 22 normalized-curve analysis was rerun from corrected traces.",
        "",
        "| Metric | Legacy | Corrected |",
        "|---|---:|---:|",
        f"| Sand–Vertisol mean distance, 0–1500 | {before['soil_mean_distance']['1500']:.6g} | {after['soil_mean_distance']['1500']:.6g} |",
        f"| Sand–Vertisol mean distance, 0–3000 | {before['soil_mean_distance']['3000']:.6g} | {after['soil_mean_distance']['3000']:.6g} |",
        f"| Maximum within-soil distance, 0–1500 | {before['within_soil']['1500']['max_all']:.6g} | {after['within_soil']['1500']['max_all']:.6g} |",
        f"| Minimum cross-soil distance, 0–1500 | {before['cross_soil']['1500']['min']:.6g} | {after['cross_soil']['1500']['min']:.6g} |",
        f"| Registered distinctness verdict | {before['pass']} | {after['pass']} |",
        "",
        "The Step 4 comparison against Stage 21 synthetic controls remains mixed-era: the Stage 22 scan curves are corrected, while Stage 21 control traces are retained unchanged, as the prompt permits correction of Stage 22 only. Its distances are not a same-formula comparison.",
        "",
        "## Step 4 curve-distance comparisons",
        "",
        "| Horizon | Comparison | Legacy | Corrected |",
        "|---|---|---:|---:|",
    ]
    for horizon in ("1500", "3000"):
        for metric, value in old["step4"]["distances"][horizon].items():
            lines.append(
                f"| 0–{horizon} | {metric} | {value:.6g} "
                f"| {new['step4']['distances'][horizon][metric]:.6g} |"
            )
    old_factorial = old["step4"]["factorial_2x2"]
    new_factorial = new["step4"]["factorial_2x2"]
    lines.extend([
        "",
        f"Mean scan-source effect at 15 µm: {old_factorial['mean_structure_source_effect_at_15um']:.6g} → {new_factorial['mean_structure_source_effect_at_15um']:.6g}.",
        f"Mean resolution effect in synthetic controls: {old_factorial['mean_resolution_effect_in_synthetic']:.6g} → {new_factorial['mean_resolution_effect_in_synthetic']:.6g}.",
        f"The larger shift remains **{old_factorial['larger_mean_curve_shift']} → {new_factorial['larger_mean_curve_shift']}**.",
        "",
        "## Time-mapped summaries",
        "",
        "| Mapping | Soil | Shape | CUE | qCO2/day | Period means R/R0 (P1–P4) |",
        "|---|---|---|---:|---:|---|",
    ])
    for mapping in ("low", "primary", "high"):
        for soil in ("sand", "vertisol"):
            old_row = old["step6"]["soil_results"][mapping][soil]
            new_row = new["step6"]["soil_results"][mapping][soil]
            old_periods = [row["mean"] for row in old_row["period_means_normalized_to_R0"]]
            new_periods = [row["mean"] for row in new_row["period_means_normalized_to_R0"]]
            lines.append(
                f"| {mapping} | {soil} | {old_row['shape']} → {new_row['shape']} "
                f"| {old_row['CUE']['mean']:.6g} → {new_row['CUE']['mean']:.6g} "
                f"| {old_row['qCO2_per_day']['mean']:.6g} → {new_row['qCO2_per_day']['mean']:.6g} "
                f"| {' / '.join(f'{v:.5g}' for v in old_periods)} → {' / '.join(f'{v:.5g}' for v in new_periods)} |"
            )
    if sensitivity:
        old_step5 = _read_json(Path("results/stage22_scan_structure/step5.json"))
        lines.extend([
            "",
            "## Sensitivity checks from stored Stage 22 traces",
            "",
            "The 32 existing sensitivity runs were post-hoc corrected using each job's configured Y; no sensitivity run was repeated.",
            "",
            "| Parameter | Factor | Legacy distance | Corrected distance | Peak amplitude order unchanged, legacy → corrected | Peak time order unchanged, legacy → corrected |",
            "|---|---:|---:|---:|---|---|",
        ])
        for parameter, factors in old_step5["perturbations"].items():
            for factor, before in factors.items():
                after = sensitivity["perturbations"][parameter][factor]
                lines.append(
                    f"| {parameter} | {factor} | {before['distance_01500']:.6g} | "
                    f"{after['distance_01500']:.6g} | {before['peak_amplitude_order']['unchanged']} → "
                    f"{after['peak_amplitude_order']['unchanged']} | "
                    f"{before['peak_time_order']['unchanged']} → {after['peak_time_order']['unchanged']} |"
                )
    path = Path("docs/stage_results/STAGE22_RFIX.md")
    if path.exists():
        first_line = path.read_text(encoding="utf-8").splitlines()[:1]
        if first_line != ["# Stage 22 scan-structure respiration correction"]:
            raise FileExistsError(f"Refusing to overwrite an unrelated report: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_stage22_posthoc_analysis() -> dict:
    source_out = Path("results/stage22_scan_structure")
    target_out = Path("results/stage22_scan_structure_rfix")
    target_out.mkdir(parents=True, exist_ok=False)
    _write_json(target_out / "step1.json", _read_json(source_out / "step1.json"))
    _write_json(target_out / "step2.json", _read_json(source_out / "step2.json"))
    corrected_runs_dir = target_out / "runs"
    corrected_runs_dir.mkdir()
    source_files = sorted(S22.RUNS.glob("*.npz"))
    for source_path in source_files:
        key = source_path.stem
        arrays = _stage22_arrays(key)
        corrected_path = corrected_runs_dir / source_path.name
        if corrected_path.exists():
            raise FileExistsError(f"Refusing to overwrite Stage 22 corrected trace: {corrected_path}")
        np.savez_compressed(
            corrected_path,
            R=arrays["R"],
            growth_R=arrays["growth_R"],
        )
    _write_json(
        target_out / "posthoc_manifest.json",
        {
            "source": str(source_out),
            "corrected_formula": "R = growth_R_legacy / Y + maint_R",
            "stage22_traces_posthoc_corrected": len(source_files),
            "fresh_stage22_runs": 0,
            "yield_overrides_read_from_job_configs": True,
        },
    )

    source_run_key = S22._run_key
    source_load_run = S22._load_run
    source_load_run_key = S22._load_run_key

    def corrected_run(key: str) -> dict:
        summary = _read_json(S22.RUNS / f"{key}.json")
        arrays = _stage22_arrays(key)
        return {**summary, **arrays}

    def corrected_load_run(structure_id: str) -> dict:
        return corrected_run(source_run_key(structure_id))

    def corrected_load_run_key(key: str) -> dict:
        return corrected_run(key)

    S22._load_run = corrected_load_run
    S22._load_run_key = corrected_load_run_key
    S22.OUT = target_out
    S22.step3()
    S22.step4()
    S22.step6()

    old = {
        f"step{step}": _read_json(source_out / f"step{step}.json")
        for step in (3, 4, 6)
    }
    new = {
        f"step{step}": _read_json(target_out / f"step{step}.json")
        for step in (3, 4, 6)
    }
    _write_stage22_report(old, new, {})
    return new


def refresh_stage22_report() -> str:
    old = {
        f"step{step}": _read_json(Path("results/stage22_scan_structure") / f"step{step}.json")
        for step in (3, 4, 6)
    }
    new = {
        f"step{step}": _read_json(Path("results/stage22_scan_structure_rfix") / f"step{step}.json")
        for step in (3, 4, 6)
    }
    sensitivity = _read_json(
        Path("results/stage22_scan_structure_rfix/step5_posthoc.json")
    )
    _write_stage22_report(old, new, sensitivity)
    return "docs/stage_results/STAGE22_RFIX.md"


def write_stage24_report() -> str:
    old = _read_json(Path("results/stage24/analysis.json"))["variants_x_outputs"]
    new = _read_json(Path("results/stage24_rfix/analysis.json"))["variants_x_outputs"]

    def display(value) -> str:
        return " / ".join(value) if isinstance(value, list) else str(value)

    lines = [
        "# Stage 24 post-hoc growth-respiration correction",
        "",
        "Stage 24 was present on the checked-out `agents/stage24-matrix-sensitivity` base. Its 184 variant traces were corrected post-hoc from their stored `growth_R` and `maint_R`; no Stage 24 simulation was rerun and the Stage 24 branch was not merged.",
        "The base variant was validated against the corrected Stage 23 analysis before comparing variants.",
        "",
        "| Variant | Output | Legacy value | Corrected value | Same/different from base, legacy | Same/different from base, corrected |",
        "|---|---|---|---|---|---|",
    ]
    output_keys = (
        "sand_vertisol_distinct",
        "peak_height_order",
        "peak_time_order",
        "arm_b_shape_primary",
        "arm_b_cumulative_rank_primary",
    )
    labels = {
        "sand_vertisol_distinct": "Arm A Sand–Vertisol distinctness",
        "peak_height_order": "Arm A peak-height order",
        "peak_time_order": "Arm A peak-time order",
        "arm_b_shape_primary": "Arm B primary shapes (S/V/L)",
        "arm_b_cumulative_rank_primary": "Arm B primary cumulative-R/mass rank",
    }
    for variant in old:
        for key in output_keys:
            old_same = old[variant]["same_as_base"][key]
            new_same = new[variant]["same_as_base"][key]
            lines.append(
                f"| {variant} | {labels[key]} | "
                f"{display(old[variant]['outputs'][key])} | "
                f"{display(new[variant]['outputs'][key])} | "
                f"{'same' if old_same else 'different'} | "
                f"{'same' if new_same else 'different'} |"
            )

    flips = []
    for variant in old:
        if (
            old[variant]["outputs"]["sand_vertisol_distinct"]
            != new[variant]["outputs"]["sand_vertisol_distinct"]
        ):
            flips.append(
                f"{variant}: Sand–Vertisol distinctness "
                f"{old[variant]['outputs']['sand_vertisol_distinct']} → "
                f"{new[variant]['outputs']['sand_vertisol_distinct']}."
            )
        if (
            old[variant]["outputs"]["arm_b_shape_primary"]
            != new[variant]["outputs"]["arm_b_shape_primary"]
        ):
            flips.append(
                f"{variant}: Arm B primary shapes "
                f"{display(old[variant]['outputs']['arm_b_shape_primary'])} → "
                f"{display(new[variant]['outputs']['arm_b_shape_primary'])}."
            )
        if (
            old[variant]["outputs"]["arm_b_cumulative_rank_primary"]
            != new[variant]["outputs"]["arm_b_cumulative_rank_primary"]
        ):
            flips.append(
                f"{variant}: Arm B primary cumulative-R/mass rank "
                f"{display(old[variant]['outputs']['arm_b_cumulative_rank_primary'])} → "
                f"{display(new[variant]['outputs']['arm_b_cumulative_rank_primary'])}."
            )
        if old[variant]["same_as_base"] != new[variant]["same_as_base"]:
            changed_comparisons = [
                key for key in old[variant]["same_as_base"]
                if old[variant]["same_as_base"][key] != new[variant]["same_as_base"][key]
            ]
            flips.append(
                f"{variant}: base comparison changed for "
                + ", ".join(labels[key] for key in changed_comparisons)
                + "."
            )
    lines.extend(["", "## What flipped", ""])
    lines.extend(f"- {flip}" for flip in flips)
    if not flips:
        lines.append("- No variant verdicts or base comparisons changed.")

    old_s16 = old["S16"]["outputs"]["sand_vertisol_distinct"]
    new_s16 = new["S16"]["outputs"]["sand_vertisol_distinct"]
    old_raw = _read_json(Path("results/stage24/analysis.json"))["variants"]["S16"]
    new_raw = _read_json(Path("results/stage24_rfix/analysis.json"))["variants"]["S16"]
    old_pair = old_raw["arm_a"]["pairwise_distances_and_distinctness"]["S-V"]["by_steps"]["1500"]
    new_pair = new_raw["arm_a"]["pairwise_distances_and_distinctness"]["S-V"]["by_steps"]["1500"]
    lines.extend([
        "",
        "## Notable criterion change",
        "",
        f"S16 Arm A Sand–Vertisol distinctness changed {old_s16} → {new_s16}. The 0–1500 soil-mean distance moved from {old_pair['soil_mean_distance']:.6g} to {new_pair['soil_mean_distance']:.6g} (criterion `> 0.3`); maximum within distance was {old_pair['max_within']:.6g} → {new_pair['max_within']:.6g}, and minimum cross distance was {old_pair['min_cross']:.6g} → {new_pair['min_cross']:.6g}.",
        "",
        "The full corrected variants table and figures are in [`results/stage24_rfix/`](../../results/stage24_rfix/).",
        "",
    ])
    target = Path("docs/stage_results/STAGE24_RFIX.md")
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite Stage 24 rfix report: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(lines), encoding="utf-8")
    return str(target)


if __name__ == "__main__":
    command = sys.argv[1] if len(sys.argv) > 1 else "gate"
    if command == "gate":
        result = run_reproduction_gates()
        print(json.dumps(result, indent=2))
    elif command == "analyze":
        result = run_stage23_analysis()
        print(json.dumps({"generated_steps": list(result), "report": "docs/stage_results/STAGE23_RFIX.md"}, indent=2))
    elif command == "report23":
        print(json.dumps({"report": refresh_stage23_report()}, indent=2))
    elif command == "refresh23":
        print(json.dumps(refresh_stage23_stage22_comparison(), indent=2))
    elif command == "stage24":
        result = run_stage24_posthoc_analysis()
        print(json.dumps(result["variants_x_outputs"], indent=2))
    elif command == "stage22":
        result = run_stage22_posthoc_analysis()
        print(json.dumps({
            "step3": result["step3"]["pass"],
            "report": "docs/stage_results/STAGE22_RFIX.md",
        }, indent=2))
    elif command == "stage22_sens":
        result = run_stage22_sensitivity_posthoc()
        print(json.dumps({
            "n_runs": result["n_runs"],
            "perturbation_groups": len(result["perturbations"]),
            "results": result["perturbations"],
        }, indent=2))
    elif command == "report22":
        print(json.dumps({"report": refresh_stage22_report()}, indent=2))
    elif command == "report24":
        print(json.dumps({"report": write_stage24_report()}, indent=2))
    elif command == "figure_a1":
        result = run_figure_a1_rfix()
        print(json.dumps(result, indent=2))
    elif command == "starvation":
        result = run_starvation_loss_report()
        print(json.dumps(result, indent=2))
    elif command == "starvation_check":
        print(json.dumps(verify_starvation_loss_report(), indent=2))
    else:
        raise SystemExit(
            "Usage: python stage23_rfix_analysis.py gate|analyze|stage22|"
            "stage22_sens|report22|stage24|report24|figure_a1|starvation|"
            "report23|refresh23"
        )
