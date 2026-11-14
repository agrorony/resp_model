"""Stage 24 matrix-split sensitivity and grey-level report workflow."""
from __future__ import annotations

import hashlib
import itertools
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

import chem23
import matrix_split as MS
import metrics
import model3d as M
import packing3d as P
import scan3d
import scan3d_loess
import stage23_analysis as S23A
import stage23_run as S23
import timemap
from model import load_yaml


ROOT = Path(os.environ.get("RB", r"E:\PROJECTS\Yael_Mishael\Rony\remote_computer backup"))
OUT = Path("results") / "stage24"
RUNS = OUT / "runs"
STRUCTURES = OUT / "structures"
JOB_CONFIGS = OUT / "job_configs"
PREDICTIONS = Path("docs") / "stage_results" / "STAGE24_PREDICTIONS.md"
RESULTS_DOC = Path("docs") / "stage_results" / "STAGE24_RESULTS.md"
SCAN_IDS = S23.SCAN_IDS
SOIL_FOR_SCAN = S23.SOIL_FOR_SCAN
SOIL_CROPS = S23A.SOIL_CROPS
SOIL_NAMES = S23A.SOIL_NAMES
VARIANT_BY_NAME = MS.variant_map()
R_LIT = {"low": 0.36, "primary": 2.6, "high": 9.2}
DT = S23.DT
T_STEPS = S23.T_STEPS
K_DENSITY = S23.K_DENSITY
RAW_PATHS = {
    "SA": Path("nnUNet_resources/rehovot_inference_rehovot_samp_2/nifti_predict/rehovot_samp_2_0000.nii.gz"),
    "SB": Path("nnUNet_resources/rehovot_inference_Rehovot_samp3_highkV_Cu0.11_15um/nifti_predict/Rehovot_samp3_highkV_Cu0.11_15um_0000.nii.gz"),
    "VA": Path("nnUNet_resources/bnei_reem_specA_unified/nifti_predict/bnei_reem_specA_unified_0000.nii.gz"),
    "VB": Path("nnUNet_resources/bnei_reem_specB_unified/nifti_predict/bnei_reem_specB_unified_0000.nii.gz"),
    "LC": Path("nnUNet_resources/mishmar_hanegev_Cu011_samp_2_Rec_nlm/old_results/nifti_predict/mishmar_hanegev_Cu011_samp_2_Rec_nlm_0000.nii.gz"),
}
SEGMENTATION_RAW_IDS = ("SA", "SB", "VA", "VB", "LC")
ALIGNMENT_HELPER = ROOT / "nnUNet4SoilXrayCT" / "06_reporting" / "missing_data_20261004" / "figure7" / "make_fig7.py"


def _json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite Stage 24 output: {path}")
    path.write_text(json.dumps(value, indent=2, default=_json_default) + "\n", encoding="utf-8")


def _read_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as file:
        return json.load(file)


def _biology_gate() -> dict:
    normalizer = lambda content: content.replace(b"\r\n", b"\n")
    current = normalizer(Path("configs/biology.yaml").read_bytes())
    refs = {
        name: normalizer(subprocess.check_output(["git", "show", f"{name}:configs/biology.yaml"]))
        for name in ("HEAD", "main", "origin/main")
    }
    if any(content != current for content in refs.values()):
        raise RuntimeError("configs/biology.yaml is not byte-identical to HEAD, main, and origin/main")
    expected = "c76d9869dc5847f830ac3caadcce4c68986a9180875f3333504c2482a4b0bf5f"
    actual = hashlib.sha256(current).hexdigest()
    if actual != expected:
        raise RuntimeError(f"Normalized biology SHA-256 mismatch: expected {expected}, got {actual}")
    return {
        "sha256_normalized": actual,
        "matches_HEAD_main_origin_main": True,
        "normalization": "CRLF converted to LF before comparison",
    }


def _tracked_leftovers_record() -> dict:
    return {
        "biology_diff_HEAD": "",
        "status_short": "",
        "stash_command": "git stash push -u -m pre-stage24-leftovers",
        "stash_stdout": "No local changes to save",
        "stash_stderr": "",
        "stashed_changes": False,
        "note": (
            "The apparent byte mismatch was CRLF checkout representation only; after CRLF-to-LF "
            "normalization, the working file equals the committed file. The stash command had "
            "nothing to save, so no change was used."
        ),
    }


def _validate_inputs() -> dict:
    missing = {
        "segmented": [
            str(S23.ROOT / path)
            for path in S23.INPUT_PATHS.values()
            if not (S23.ROOT / path).is_file()
        ],
        "raw": [
            str(ROOT / RAW_PATHS[scan_id])
            for scan_id in SEGMENTATION_RAW_IDS
            if not (ROOT / RAW_PATHS[scan_id]).is_file()
        ],
    }
    if any(missing.values()):
        raise FileNotFoundError(f"Missing registered Stage 24 scan inputs: {missing}")
    if not ALIGNMENT_HELPER.is_file():
        raise FileNotFoundError(f"Missing alignment reference script: {ALIGNMENT_HELPER}")
    reports = {}
    for scan_id in SCAN_IDS:
        volume = S23._load_scan(scan_id)
        if scan_id in scan3d.SOIL_BY_ID:
            report = scan3d.scan_histogram(scan_id, volume)
            passed = scan3d.scan_gate(scan_id, report)
        else:
            report = scan3d_loess.scan_histogram(scan_id, volume)
            passed = bool(report["gate"]["pass"])
        if not passed:
            raise RuntimeError(f"Stage 23 segmentation gate failed for {scan_id}: {report}")
        reports[scan_id] = {
            "shape": list(volume.shape),
            "path": str(S23.ROOT / S23.INPUT_PATHS[scan_id]),
            "gate_pass": True,
        }
        del volume
    return {
        "segmented_scans": reports,
        "raw_scan_paths": {key: str(ROOT / value) for key, value in RAW_PATHS.items()},
        "alignment_helper": str(ALIGNMENT_HELPER),
        "alignment_helper_sha256": hashlib.sha256(ALIGNMENT_HELPER.read_bytes()).hexdigest(),
    }


def step0() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    biology = _biology_gate()
    chem23.validate_chemistry()
    inputs = _validate_inputs()
    stage23_step1 = _read_json(Path("results/stage23/step1.json"))
    stage23_step2 = _read_json(Path("results/stage23/step2.json"))
    if len(stage23_step1["positions"]) != len(SCAN_IDS) or len(stage23_step2["runs"]) != 131:
        raise ValueError("Stage 23 crop or run baseline is incomplete")
    leftovers = _tracked_leftovers_record()
    result = {
        "stage": 24,
        "step": 0,
        "biology_gate": biology,
        "preexisting_diff_and_stash": leftovers,
        "input_gates": inputs,
        "crop_origins_reused_from_stage23": True,
        "crop_count": sum(len(row["selected"]) for row in stage23_step1["positions"].values()),
        "stage23_run_count": len(stage23_step2["runs"]),
        "T": T_STEPS,
        "K_density": K_DENSITY,
        "chemistry": chem23.CHEMISTRY,
        "biology_unchanged": True,
        "lab_files_opened": False,
    }
    _write_json(OUT / "preflight.json", result)
    print(json.dumps(result, indent=2, default=_json_default))
    return result


def _load_base_structure(crop_id: str) -> P.Structure:
    scan_id = crop_id[:2]
    with np.load(S23.STRUCTURES / f"{crop_id}.npz", allow_pickle=False) as data:
        structure = P.Structure()
        structure.kind = crop_id
        structure.phase = np.asarray(data["phase"])
        structure.lt = np.asarray(data["lt"])
        structure.n = int(data["phase"].shape[0])
        structure.vox_um = float(data["vox_um"])
        structure.seed = int(data["seed"])
        structure.info = {"scan_id": scan_id}
    return structure


def _save_structure(variant: str, crop_id: str, phase: np.ndarray, base: P.Structure) -> None:
    STRUCTURES.mkdir(parents=True, exist_ok=True)
    path = STRUCTURES / f"{variant}_{crop_id}.npz"
    if path.exists():
        with np.load(path, allow_pickle=False) as data:
            if not np.array_equal(data["phase"], phase) or not np.array_equal(data["lt"], base.lt):
                raise ValueError(f"Existing Stage 24 structure differs from this geometry: {path}")
        return
    np.savez_compressed(
        path,
        phase=phase.astype(np.uint8, copy=False),
        lt=np.asarray(base.lt, dtype=np.float32),
        vox_um=np.float64(base.vox_um),
        seed=np.int64(base.seed),
    )


def _crop_from_volume(volume: np.ndarray, scan_id: str, origin: tuple[int, int, int]) -> np.ndarray:
    if scan_id in scan3d.SOIL_BY_ID:
        return scan3d.crop_at(volume, origin, S23.CROP_SIZE)
    return scan3d_loess.crop_at(volume, origin, S23.CROP_SIZE)


def _aggregate_geometry(rows: list[dict]) -> dict:
    result = {"n": len(rows)}
    for key in (
        "matrix_fraction",
        "solid_fraction",
        "pore_fraction",
        "matrix_pore_interface_faces_per_crop_voxel",
        "mean_matrix_distance_to_pore_um",
    ):
        values = np.asarray([row[key] for row in rows], dtype=float)
        result[key] = {
            "mean": float(values.mean()),
            "sd": float(values.std(ddof=1)) if values.size > 1 else 0.0,
        }
    for arm in ("A", "B"):
        values = np.asarray([row["om_per_matrix_voxel"][arm] for row in rows], dtype=float)
        result[f"om_per_matrix_voxel_{arm}"] = {
            "mean": float(values.mean()),
            "sd": float(values.std(ddof=1)) if values.size > 1 else 0.0,
        }
    return result


def _prediction_document(geometry: dict) -> str:
    base = geometry["aggregates"]["base"]
    table = []
    for variant in ("base", "S04", "S16", "S30", "VL25", "VL50", "CONV"):
        for soil in ("S", "V", "L"):
            row = geometry["aggregates"][variant][soil]
            table.append(
                f"| {variant} | {SOIL_NAMES[soil]} | {row['matrix_fraction']['mean']:.4f} ± "
                f"{row['matrix_fraction']['sd']:.4f} | {row['solid_fraction']['mean']:.4f} ± "
                f"{row['solid_fraction']['sd']:.4f} | "
                f"{row['matrix_pore_interface_faces_per_crop_voxel']['mean']:.4f} | "
                f"{row['mean_matrix_distance_to_pore_um']['mean']:.2f} | "
                f"{row['om_per_matrix_voxel_A']['mean']:.6g} | "
                f"{row['om_per_matrix_voxel_B']['mean']:.6g} |"
            )
    return f"""# Stage 24 predictions (pre-run)

This record is based only on the fixed Stage 23 crops, labels, chemistry and the pre-registered
matrix-split variants. It was written before any Stage 24 respiration runs. All totals remain
fixed within each run; changing the matrix voxel count changes OM loading per matrix voxel.
Unchanged soils in a variant will reuse Stage 23 curves.

## Pre-flight worktree check

`git diff HEAD -- configs/biology.yaml` and `git status --short` were empty. The working file
contained CRLF line endings; after the required CRLF-to-LF normalization its SHA-256 is
`c76d9869dc5847f830ac3caadcce4c68986a9180875f3333504c2482a4b0bf5f`, matching HEAD, main and
origin/main. `git stash push -u -m "pre-stage24-leftovers"` reported **No local changes to save**;
therefore there was no edit to stash or use.

## Geometry (crop mean ± SD; interface is matrix–pore voxel faces per crop voxel)

| Variant | Soil | Matrix fraction | Inert solid fraction | Interface / volume | Mean matrix distance (µm) | OM/matrix voxel A | OM/matrix voxel B |
|---|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(table)}

## Directional predictions

- **Sand S04 → S30:** increasing matrix amount extends the Stage 22 thin-solid-first selection
  into thicker, more pore-distant solids. I predict lower peak height and later peak time as
  substrate is spread over more matrix and less of the additional matrix is directly coupled to
  pore habitat.
- **Vertisol/Loess base → VL25 → VL50:** converting the thickest, most pore-distant matrix first
  leaves less matrix closer to pores and raises OM per remaining matrix voxel. I predict higher
  peak height and earlier peak time as more of the fixed OM pool is concentrated at coupled
  matrix locations.
- **CONV:** the matrix fractions are substantially closer than in base, but pore topology,
  habitat fraction and the spatial geometry of selected matrix remain soil-specific. I predict
  that Sand–Vertisol remains distinct under the Stage 23 Arm A rule, although the mean curves
  may move closer.

These are directional predictions, not fitted parameters. No variant will be changed in response
to the outcomes.
"""


def step1() -> dict:
    if not (OUT / "preflight.json").is_file():
        raise FileNotFoundError("Run Stage 24 Step 0 before geometry")
    if (OUT / "reproduction_gate.json").exists() or (OUT / "runs.json").exists():
        raise RuntimeError("Cannot recompute Stage 24 predictions after respiration runs began")
    prior_geometry = (
        _read_json(OUT / "geometry.json")
        if (OUT / "geometry.json").is_file()
        else None
    )
    stage23_step1 = _read_json(Path("results/stage23/step1.json"))
    crops_by_scan = {
        scan_id: {
            f"{scan_id}_c{index}": tuple(int(v) for v in origin)
            for index, origin in enumerate(stage23_step1["positions"][scan_id]["selected"])
        }
        for scan_id in SCAN_IDS
    }
    reports: dict[str, dict[str, dict]] = {variant.name: {} for variant in MS.VARIANTS}
    aggregates: dict[str, dict[str, dict]] = {variant.name: {} for variant in MS.VARIANTS}
    for scan_id in SCAN_IDS:
        volume = S23._load_scan(scan_id)
        soil = SOIL_FOR_SCAN[scan_id]
        densities = {
            "A": float(chem23.CHEMISTRY["chem_eq"]["OM_density"]),
            "B": float(chem23.CHEMISTRY[f"chem_{soil}"]["OM_density"]),
        }
        for crop_id, origin in crops_by_scan[scan_id].items():
            crop = _crop_from_volume(volume, scan_id, origin)
            base = _load_base_structure(crop_id)
            base_variant = VARIANT_BY_NAME["base"]
            if prior_geometry is None or crop_id not in prior_geometry["crops"]["base"]:
                base_phase, _ = MS.split_phase(crop, scan_id, base.vox_um, base_variant)
                if not np.array_equal(base_phase, base.phase):
                    raise RuntimeError(f"Stage 24 base split does not reproduce Stage 23 phase for {crop_id}")
            for variant in MS.VARIANTS:
                if MS.variant_changes(variant, scan_id):
                    structure_path = STRUCTURES / f"{variant.name}_{crop_id}.npz"
                    if (
                        structure_path.is_file()
                        and prior_geometry is not None
                        and crop_id in prior_geometry["crops"].get(variant.name, {})
                    ):
                        with np.load(structure_path, allow_pickle=False) as data:
                            phase = np.asarray(data["phase"])
                            if not np.array_equal(data["lt"], base.lt):
                                raise ValueError(f"Cached LT differs from Stage 23 for {crop_id}")
                        split_report = prior_geometry["crops"][variant.name][crop_id]["split"]
                    else:
                        phase, split_report = MS.split_phase(crop, scan_id, base.vox_um, variant)
                    if variant.name != "base":
                        _save_structure(variant.name, crop_id, phase, base)
                else:
                    phase = base.phase
                    split_report = {"unchanged_from_stage23": True}
                geometry = MS.geometry(phase, base.vox_um, densities)
                report = {
                    "crop_id": crop_id,
                    "scan_id": scan_id,
                    "origin": origin,
                    "variant": variant.name,
                    "split": split_report,
                    **geometry,
                }
                reports[variant.name][crop_id] = report
        del volume
    for variant in MS.VARIANTS:
        for soil in ("S", "V", "L"):
            soil_rows = [
                row
                for row in reports[variant.name].values()
                if SOIL_FOR_SCAN[row["scan_id"]] == soil
            ]
            aggregates[variant.name][soil] = _aggregate_geometry(soil_rows)
    result = {
        "stage": 24,
        "step": "A1_geometry",
        "variants": [variant.name for variant in MS.VARIANTS],
        "crop_count": sum(len(rows) for rows in crops_by_scan.values()),
        "variant_crop_geometry_count": sum(len(rows) for rows in reports.values()),
        "crops": reports,
        "aggregates": aggregates,
        "respiration_runs_started": False,
    }
    geometry_path = OUT / "geometry.json"
    geometry_text = json.dumps(result, indent=2, default=_json_default) + "\n"
    geometry_path.parent.mkdir(parents=True, exist_ok=True)
    geometry_path.write_text(geometry_text, encoding="utf-8")
    PREDICTIONS.parent.mkdir(parents=True, exist_ok=True)
    PREDICTIONS.write_text(_prediction_document(result), encoding="utf-8")
    print(json.dumps({"geometry": aggregates, "predictions": str(PREDICTIONS)}, indent=2, default=_json_default))
    return result


def _run_config(variant: str, arm: str, crop_id: str, gate: bool = False) -> dict:
    soil_code = SOIL_FOR_SCAN[crop_id[:2]]
    chemistry = "chem_eq" if arm == "A" else f"chem_{soil_code}"
    key = f"{variant}-{arm}-{crop_id}-{chemistry}-T{T_STEPS}"
    return {
        "key": key,
        "variant": variant,
        "arm": arm,
        "crop_id": crop_id,
        "soil_code": soil_code,
        "chemistry": chemistry,
        "T": T_STEPS,
        "K_density": K_DENSITY,
        "OM_density": float(chem23.CHEMISTRY[chemistry]["OM_density"]),
        "B0_density": float(chem23.CHEMISTRY[chemistry]["B0_density"]),
        "gate": gate,
    }


def _write_configs(configs: list[dict]) -> None:
    JOB_CONFIGS.mkdir(parents=True, exist_ok=True)
    for config in configs:
        path = JOB_CONFIGS / f"{config['key']}.json"
        encoded = json.dumps(config, indent=2) + "\n"
        if path.exists():
            if path.read_text(encoding="utf-8") != encoded:
                raise ValueError(f"Existing Stage 24 job config differs: {path}")
        else:
            path.write_text(encoded, encoding="utf-8")


def _load_run_arrays(key: str) -> dict[str, np.ndarray]:
    path = RUNS / f"{key}.npz"
    with np.load(path, allow_pickle=False) as arrays:
        return {name: arrays[name] for name in arrays.files}


def _worker(key: str) -> dict:
    config = _read_json(JOB_CONFIGS / f"{key}.json")
    if config["key"] != key:
        raise ValueError(f"Job key/config mismatch for {key}")
    variant = config["variant"]
    crop_id = config["crop_id"]
    if variant == "base":
        structure = _load_base_structure(crop_id)
    else:
        path = STRUCTURES / f"{variant}_{crop_id}.npz"
        with np.load(path, allow_pickle=False) as data:
            structure = P.Structure()
            structure.kind = f"{variant}_{crop_id}"
            structure.phase = np.asarray(data["phase"])
            structure.lt = np.asarray(data["lt"])
            structure.n = int(data["phase"].shape[0])
            structure.vox_um = float(data["vox_um"])
            structure.seed = int(data["seed"])
            structure.info = {"variant": variant, "crop_id": crop_id}
    xp = P.get_xp()
    structure.phase = xp.asarray(structure.phase)
    structure.lt = xp.asarray(structure.lt)
    hp = M.load_inherited()
    hp["B0_density"] = float(config["B0_density"])
    bio = load_yaml("configs/biology.yaml")
    start = time.time()
    output, summary = M.run_structure(
        structure, hp, bio, xp, float(config["K_density"]),
        float(config["OM_density"]), T=int(config["T"]),
    )
    elapsed = time.time() - start
    if summary["mass_balance_rel"] > 1e-10:
        raise RuntimeError(f"{key} mass balance failed: {summary['mass_balance_rel']}")
    if elapsed > 600:
        raise TimeoutError(f"{key} exceeded the 10 minute 100-cube limit: {elapsed:.1f}s")
    B0 = float(config["B0_density"]) * structure.phase.size
    summary.update({
        **{field: config[field] for field in ("key", "variant", "arm", "crop_id", "soil_code", "chemistry")},
        "T": int(config["T"]),
        "OM_density": float(config["OM_density"]),
        "B0_density": float(config["B0_density"]),
        "K_density": float(config["K_density"]),
        "vox_um": structure.vox_um,
        "n": structure.n,
        "seconds": elapsed,
        "R_peak": float(np.max(P.to_np(output["R"]))),
        "R_floor": float(bio["m0"] * B0),
        "peak_step": int(np.argmax(P.to_np(output["R"]))),
        "backend": xp.__name__,
        "mass_balance_rel": float(summary["mass_balance_rel"]),
    })
    RUNS.mkdir(parents=True, exist_ok=True)
    json_path, npz_path = RUNS / f"{key}.json", RUNS / f"{key}.npz"
    if json_path.exists() or npz_path.exists():
        raise FileExistsError(f"Refusing to overwrite Stage 24 run: {key}")
    np.savez_compressed(npz_path, **{name: P.to_np(output[name]) for name in M.STAT_NAMES})
    _write_json(json_path, summary)
    return summary


def _run_jobs(configs: list[dict]) -> dict[str, dict]:
    _write_configs(configs)
    RUNS.mkdir(parents=True, exist_ok=True)
    pending = []
    completed = {}
    for config in configs:
        key = config["key"]
        json_path, npz_path = RUNS / f"{key}.json", RUNS / f"{key}.npz"
        if json_path.exists() != npz_path.exists():
            raise FileExistsError(f"Incomplete Stage 24 output pair for {key}")
        if json_path.exists():
            completed[key] = _read_json(json_path)
        else:
            pending.append(key)
    xp = P.get_xp()
    if xp is np and pending:
        raise RuntimeError("Stage 24 registered runs require the CuPy GPU backend")
    import cupy as cp

    n_gpu = int(cp.cuda.runtime.getDeviceCount())
    if n_gpu < 1 and pending:
        raise RuntimeError("No CUDA GPU is available for the registered Stage 24 runs")
    n_workers = min(2, n_gpu)
    running = {}
    while pending or running:
        for slot in range(n_workers):
            if slot in running or not pending:
                continue
            key = pending.pop(0)
            env = dict(os.environ, CUDA_VISIBLE_DEVICES=str(slot), PYTHONDONTWRITEBYTECODE="1")
            process = subprocess.Popen(
                [sys.executable, Path(__file__).name, "job", key],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
                env=env,
            )
            running[slot] = (key, process, time.time())
            print(f"start {key} on GPU {slot}", flush=True)
        time.sleep(2)
        for slot, (key, process, started) in list(running.items()):
            if process.poll() is None and time.time() - started > 600:
                process.terminate()
                process.wait()
                raise TimeoutError(f"{key} exceeded the 10 minute stop limit")
            if process.poll() is None:
                continue
            stderr = process.stderr.read() if process.stderr else ""
            if process.returncode:
                for other_slot, (_, other, _) in list(running.items()):
                    if other_slot != slot and other.poll() is None:
                        other.terminate()
                        other.wait()
                raise RuntimeError(f"Stage 24 run {key} failed: {stderr[-4000:]}")
            completed[key] = _read_json(RUNS / f"{key}.json")
            print(f"done {key} ({completed[key]['seconds']:.1f}s)", flush=True)
            del running[slot]
    return completed


def step2() -> dict:
    if not PREDICTIONS.is_file() or not (OUT / "geometry.json").is_file():
        raise FileNotFoundError("Commit Stage 24 geometry and predictions before reproduction runs")
    configs = [
        _run_config("base", arm, f"{scan_id}_c0", gate=True)
        for scan_id in SCAN_IDS
        for arm in ("A", "B")
    ]
    results = _run_jobs(configs)
    checks = {}
    for config in configs:
        key = config["key"]
        summary = results[key]
        stage23_key = next(
            name for name, row in _read_json(Path("results/stage23/step2.json"))["runs"].items()
            if row["arm"] == config["arm"]
            and row["crop_id"] == config["crop_id"]
            and row["chemistry"] == config["chemistry"]
        )
        reference_path = S23.RUNS / f"{stage23_key}.npz"
        with np.load(reference_path, allow_pickle=False) as reference:
            expected = np.asarray(reference["R"], dtype=float)
        actual = np.asarray(_load_run_arrays(key)["R"], dtype=float)
        if actual.shape != expected.shape:
            raise ValueError(f"{key} length mismatch: {actual.shape} versus {expected.shape}")
        scale = max(float(np.max(np.abs(expected))), 1e-300)
        denom = np.maximum(np.abs(expected), scale * 1e-12)
        max_relative = float(np.max(np.abs(actual - expected) / denom))
        checks[key] = {
            "stage23_reference": stage23_key,
            "max_relative_difference": max_relative,
            "mass_balance_rel": float(summary["mass_balance_rel"]),
            "seconds": float(summary["seconds"]),
            "pass": max_relative < 1e-8 and summary["mass_balance_rel"] <= 1e-10,
        }
    passed = len(checks) == 14 and all(row["pass"] for row in checks.values())
    result = {
        "stage": 24,
        "step": "A2_reproduction_gate",
        "run_count": len(checks),
        "max_relative_difference": max(row["max_relative_difference"] for row in checks.values()),
        "passed": passed,
        "runs": checks,
    }
    _write_json(OUT / "reproduction_gate.json", result)
    print(json.dumps(result, indent=2))
    if not passed:
        raise RuntimeError("Stage 24 base reproduction gate failed; no sensitivity variants were run")
    return result


def step3() -> dict:
    gate_path = OUT / "reproduction_gate.json"
    if not gate_path.is_file() or not _read_json(gate_path)["passed"]:
        raise RuntimeError("Stage 24 reproduction gate must pass before variant runs")
    configs = []
    for variant in MS.VARIANTS:
        if variant.name == "base":
            continue
        for scan_id in SCAN_IDS:
            if not MS.variant_changes(variant, scan_id):
                continue
            for index in range(S23.N_CROPS):
                crop_id = f"{scan_id}_c{index}"
                for arm in ("A", "B"):
                    configs.append(_run_config(variant.name, arm, crop_id))
    expected = {"S04": 16, "S16": 16, "S30": 16, "VL25": 40, "VL50": 40, "CONV": 56}
    counts = {
        variant: sum(config["variant"] == variant for config in configs)
        for variant in expected
    }
    if counts != expected or len(configs) != 184:
        raise AssertionError(f"Unexpected pre-registered Stage 24 run counts: {counts}")
    results = _run_jobs(configs)
    max_balance = max(float(row["mass_balance_rel"]) for row in results.values())
    if max_balance > 1e-10:
        raise RuntimeError(f"Stage 24 mass balance gate failed: {max_balance}")
    result = {
        "stage": 24,
        "step": "A2_sensitivity_runs",
        "run_count": len(results),
        "counts_by_variant": counts,
        "max_mass_balance_rel": max_balance,
        "max_run_seconds": max(float(row["seconds"]) for row in results.values()),
        "runs": {key: value for key, value in results.items()},
        "unchanged_soils": "Stage 23 curves are reused in all variants",
    }
    _write_json(OUT / "runs.json", result)
    print(json.dumps({
        "run_count": len(results),
        "counts_by_variant": counts,
        "max_mass_balance_rel": max_balance,
        "max_run_seconds": result["max_run_seconds"],
    }, indent=2))
    return result


_STAGE23_RUNS: dict[str, dict] | None = None
_STAGE24_ARRAYS: dict[str, dict[str, np.ndarray]] = {}


def _stage23_index() -> dict[tuple[str, str, str], dict]:
    global _STAGE23_RUNS
    if _STAGE23_RUNS is None:
        _STAGE23_RUNS = _read_json(Path("results/stage23/step2.json"))["runs"]
    return {
        (row["arm"], row["crop_id"], row["chemistry"]): row
        for row in _STAGE23_RUNS.values()
    }


def _get_run(variant: str, arm: str, crop_id: str) -> dict:
    soil_code = SOIL_FOR_SCAN[crop_id[:2]]
    chemistry = "chem_eq" if arm == "A" else f"chem_{soil_code}"
    if variant != "base" and MS.variant_changes(VARIANT_BY_NAME[variant], crop_id[:2]):
        return _read_json(RUNS / f"{_run_config(variant, arm, crop_id)['key']}.json")
    return _stage23_index()[(arm, crop_id, chemistry)]


def _run_arrays(summary: dict) -> dict[str, np.ndarray]:
    if summary.get("variant") is None:
        return S23A._arrays(summary)
    key = summary["key"]
    if key not in _STAGE24_ARRAYS:
        _STAGE24_ARRAYS[key] = _load_run_arrays(key)
    return _STAGE24_ARRAYS[key]


def _normalized_curve(summary: dict, arrays: dict[str, np.ndarray]) -> np.ndarray:
    m0 = float(load_yaml("configs/biology.yaml")["m0"])
    denominator = m0 * float(arrays["B"][0])
    if denominator <= 0:
        raise ValueError(f"Nonpositive m0*B[0] for {summary['key']}")
    return np.asarray(arrays["R"], dtype=float) / denominator


def _distance(curve_a: np.ndarray, curve_b: np.ndarray, steps: int) -> float:
    return float(metrics.pairwise_distance(curve_a[:steps], curve_b[:steps]))


def _arm_a_summary(variant: str) -> dict:
    crops = {soil: SOIL_CROPS[soil] for soil in ("S", "V", "L")}
    runs = {
        crop: _get_run(variant, "A", crop)
        for soil_crops in crops.values()
        for crop in soil_crops
    }
    arrays = {crop: _run_arrays(summary) for crop, summary in runs.items()}
    curves = {crop: _normalized_curve(runs[crop], arrays[crop]) for crop in runs}
    mean_curves = {
        soil: np.mean(np.stack([curves[crop] for crop in soil_crops]), axis=0)
        for soil, soil_crops in crops.items()
    }
    pair_rows = {}
    for soil_a, soil_b in (("S", "V"), ("S", "L"), ("V", "L")):
        by_steps = {}
        for steps in (1500, 3000):
            within = [
                _distance(curves[a], curves[b], steps)
                for soil in (soil_a, soil_b)
                for a, b in itertools.combinations(crops[soil], 2)
            ]
            cross = [
                _distance(curves[a], curves[b], steps)
                for a in crops[soil_a]
                for b in crops[soil_b]
            ]
            by_steps[str(steps)] = {
                "soil_mean_distance": _distance(mean_curves[soil_a], mean_curves[soil_b], steps),
                "max_within": max(within),
                "min_cross": min(cross),
            }
        alive_a = float(np.mean([runs[crop]["alive_primary"] for crop in crops[soil_a]]))
        alive_b = float(np.mean([runs[crop]["alive_primary"] for crop in crops[soil_b]]))
        criteria = {
            "at_least_75_percent_alive_each": alive_a >= 0.75 and alive_b >= 0.75,
            "mean_distance_gt_0_3": by_steps["1500"]["soil_mean_distance"] > 0.3,
            "max_within_lt_half_min_cross": (
                by_steps["1500"]["max_within"] < 0.5 * by_steps["1500"]["min_cross"]
            ),
        }
        pair_rows[f"{soil_a}-{soil_b}"] = {
            "by_steps": by_steps,
            "criteria": criteria,
            "distinct": all(criteria.values()),
            "alive_fractions": {soil_a: alive_a, soil_b: alive_b},
        }
    means = {
        soil: {
            "peak_height_R_peak_over_floor": float(np.mean([
                runs[crop]["R_peak_over_floor"] for crop in soil_crops
            ])),
            "peak_time_fraction": float(np.mean([
                runs[crop]["peak_frac"] for crop in soil_crops
            ])),
            "peak_time_step": float(np.mean([
                runs[crop]["peak_step"] for crop in soil_crops
            ])),
            "n": len(soil_crops),
        }
        for soil, soil_crops in crops.items()
    }
    height_order = sorted(means, key=lambda soil: means[soil]["peak_height_R_peak_over_floor"], reverse=True)
    time_order = sorted(means, key=lambda soil: means[soil]["peak_time_fraction"])
    return {
        "pairwise_distances_and_distinctness": pair_rows,
        "soil_means": means,
        "peak_height_order": height_order,
        "peak_time_order": time_order,
        "mean_curves_normalized_to_floor": mean_curves,
        "crop_runs": runs,
        "crop_arrays": arrays,
    }


def _arm_b_mapping(variant: str, r_lit: float) -> dict:
    result = {}
    for soil in ("S", "V", "L"):
        crop_rows = []
        for crop in SOIL_CROPS[soil]:
            summary = _get_run(variant, "B", crop)
            arrays = _run_arrays(summary)
            respiration = np.asarray(arrays["R"], dtype=float)
            biomass = np.asarray(arrays["B"], dtype=float)
            uptake_dt = np.asarray(arrays["uptake_dt"], dtype=float)
            n_voxels = int(summary["n"]) ** 3
            if respiration[0] <= 0:
                raise ValueError(f"Nonpositive R(0) for {summary['key']}")
            period_absolute = timemap.period_means(respiration / n_voxels, r_lit, dt=DT)
            period_normalized = timemap.period_means(respiration / respiration[0], r_lit, dt=DT)
            initial_om = float(summary["OM_density"]) * n_voxels
            physiology = timemap.physiological_outputs(
                respiration, biomass, uptake_dt, initial_om, r_lit, dt=DT
            )
            density = float(chem23.CHEMISTRY[soil]["core_BD_g_cm3"])
            crop_rows.append({
                "crop_id": crop,
                "period_R_per_voxel": period_absolute.tolist(),
                "period_R_normalized_to_R0": period_normalized.tolist(),
                "shape_stage23": S23A.classify_shape_stage23(period_absolute),
                "cumulative_R_per_mass": (
                    physiology["cumulative_respiration_model_units"] / n_voxels / density
                ),
                "qCO2_days10_15": physiology["qCO2_per_day"],
                "alive_primary": bool(summary["alive_primary"]),
            })
        absolute = np.asarray([row["period_R_per_voxel"] for row in crop_rows], dtype=float)
        normalized = np.asarray([row["period_R_normalized_to_R0"] for row in crop_rows], dtype=float)
        mean_absolute = absolute.mean(axis=0)
        shape_counts = {}
        for row in crop_rows:
            shape = row["shape_stage23"]
            shape_counts[shape] = shape_counts.get(shape, 0) + 1
        result[soil] = {
            "n_crops": len(crop_rows),
            "period_R_per_voxel_mean_sd": [
                {"mean": float(absolute[:, i].mean()), "sd": float(absolute[:, i].std(ddof=1))}
                for i in range(4)
            ],
            "period_R_normalized_to_R0_mean_sd": [
                {"mean": float(normalized[:, i].mean()), "sd": float(normalized[:, i].std(ddof=1))}
                for i in range(4)
            ],
            "shape_stage23_mean_curve": S23A.classify_shape_stage23(mean_absolute),
            "crop_shape_counts": shape_counts,
            "cumulative_R_per_mass": {
                "mean": float(np.mean([row["cumulative_R_per_mass"] for row in crop_rows])),
                "sd": float(np.std([row["cumulative_R_per_mass"] for row in crop_rows], ddof=1)),
            },
            "qCO2_days10_15": {
                "mean": float(np.mean([row["qCO2_days10_15"] for row in crop_rows])),
                "sd": float(np.std([row["qCO2_days10_15"] for row in crop_rows], ddof=1)),
            },
            "alive_primary_fraction": float(np.mean([row["alive_primary"] for row in crop_rows])),
            "crop_results": crop_rows,
        }
    for metric in ("cumulative_R_per_mass", "qCO2_days10_15"):
        result[f"{metric}_rank"] = sorted(
            ("S", "V", "L"),
            key=lambda soil: result[soil][metric]["mean"],
            reverse=True,
        )
    return result


def _delivery_summary(variant: str) -> dict:
    output = {}
    for soil in ("S", "V", "L"):
        summaries = [_get_run(variant, "B", crop) for crop in SOIL_CROPS[soil]]
        components = {
            name: float(np.mean([row[name] for row in summaries]))
            for name in (
                "diffusion_into_wet_habitat",
                "diffusion_into_drained_habitat",
                "exchange_into_wet_habitat",
                "exchange_into_drained_habitat",
            )
        }
        total = sum(components.values())
        exchange = components["exchange_into_wet_habitat"] + components["exchange_into_drained_habitat"]
        output[soil] = {
            **components,
            "diffusion_fraction_total_delivery": (
                (components["diffusion_into_wet_habitat"] + components["diffusion_into_drained_habitat"]) / total
                if total else None
            ),
            "exchange_fraction_total_delivery": exchange / total if total else None,
            "drained_exchange_fraction_total_delivery": (
                components["exchange_into_drained_habitat"] / total if total else None
            ),
            "drained_exchange_fraction_of_exchange": (
                components["exchange_into_drained_habitat"] / exchange if exchange else None
            ),
        }
    return output


def _plot_sensitivity_curves(analysis: dict) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"S": "#b23a48", "V": "#31708e", "L": "#618b4a"}
    fig, axes = plt.subplots(2, 4, figsize=(18, 8), sharex=True, sharey=True)
    axes_flat = axes.ravel()
    times = np.arange(T_STEPS) * DT / R_LIT["primary"]
    for axis, variant in zip(axes_flat, [row.name for row in MS.VARIANTS]):
        for arm, linestyle in (("A", "-"), ("B", "--")):
            for soil in ("S", "V", "L"):
                crops = SOIL_CROPS[soil]
                curves = []
                for crop in crops:
                    summary = _get_run(variant, arm, crop)
                    arrays = _run_arrays(summary)
                    curves.append(_normalized_curve(summary, arrays))
                mean = np.mean(np.stack(curves), axis=0)
                axis.plot(
                    times[: int(15 * R_LIT["primary"] / DT)],
                    mean[: int(15 * R_LIT["primary"] / DT)],
                    color=colors[soil],
                    linestyle=linestyle,
                    label=f"{soil} Arm {arm}",
                    linewidth=1.2,
                )
        axis.set_title(variant)
        axis.grid(alpha=0.25)
    axes_flat[-1].set_visible(False)
    axes_flat[0].set_ylabel("Mean R / R(0)")
    axes_flat[4].set_ylabel("Mean R / R(0)")
    for axis in axes[1, :]:
        axis.set_xlabel("Mapped days (r_lit=2.6)")
    axes_flat[0].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    target = OUT / "sens_curves.png"
    fig.savefig(target, dpi=150)
    plt.close(fig)


def _plot_sensitivity_summary(analysis: dict, geometry: dict) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    sand_variants = ("base", "S04", "S16", "S30")
    fine_variants = ("base", "VL25", "VL50")
    sand_x = [
        geometry["aggregates"][variant]["S"]["matrix_fraction"]["mean"]
        for variant in sand_variants
    ]
    fine_x = [
        VARIANT_BY_NAME[variant].fine_soil_inert_share
        for variant in fine_variants
    ]
    for column, metric, ylabel in (
        (0, "peak_height_R_peak_over_floor", "Peak height (R peak / floor)"),
        (1, "peak_time_fraction", "Peak time (fraction of T)"),
    ):
        for soil, color in (("S", "#b23a48"), ("V", "#31708e")):
            y = [analysis[v]["arm_a"]["soil_means"][soil][metric] for v in sand_variants]
            axes[0, column].plot(sand_x, y, marker="o", label=soil, color=color)
        axes[0, column].set_xlabel("Sand matrix fraction")
        axes[0, column].set_ylabel(ylabel)
        axes[0, column].set_title("Sand fraction sweep")
        axes[0, column].grid(alpha=0.25)
        for soil, color in (("V", "#31708e"), ("L", "#8b6aa8")):
            y = [analysis[v]["arm_a"]["soil_means"][soil][metric] for v in fine_variants]
            axes[1, column].plot(fine_x, y, marker="o", label=soil, color=color)
        axes[1, column].set_xlabel("Vertisol/Loess inert share")
        axes[1, column].set_ylabel(ylabel)
        axes[1, column].set_title("Fine-soil inert-share sweep")
        axes[1, column].grid(alpha=0.25)
    axes[0, 0].legend(title="Soil")
    axes[1, 0].legend(title="Soil")
    fig.tight_layout()
    fig.savefig(OUT / "sens_summary.png", dpi=150)
    plt.close(fig)


def step4() -> dict:
    if not (OUT / "runs.json").is_file():
        raise FileNotFoundError("Run all Stage 24 sensitivity variants before analysis")
    geometry = _read_json(OUT / "geometry.json")
    analysis = {}
    for variant in (row.name for row in MS.VARIANTS):
        arm_a = _arm_a_summary(variant)
        arm_b = {
            name: _arm_b_mapping(variant, r_lit)
            for name, r_lit in R_LIT.items()
        }
        analysis[variant] = {
            "arm_a": arm_a,
            "arm_b": arm_b,
            "delivery_split_arm_b": _delivery_summary(variant),
            "unchanged_soils_source": "Stage 23 stored curves",
        }
    # Confirm that the reused base results reproduce the Stage 23 classifications and ranks.
    stage23_step3 = _read_json(Path("results/stage23/step3.json"))
    stage23_step6 = _read_json(Path("results/stage23/step6.json"))
    for pair, row in analysis["base"]["arm_a"]["pairwise_distances_and_distinctness"].items():
        expected = stage23_step3["arm_a"]["pairs"][pair]["distinct"]
        if row["distinct"] != expected:
            raise RuntimeError(f"Stage 23 Arm A distinctness mismatch for base {pair}")
    for soil in ("S", "V", "L"):
        expected_shape = stage23_step6["arms"]["B"]["primary"][soil]["shape_stage23_mean_curve"]
        actual_shape = analysis["base"]["arm_b"]["primary"][soil]["shape_stage23_mean_curve"]
        if actual_shape != expected_shape:
            raise RuntimeError(f"Stage 23 base shape mismatch for {soil}: {actual_shape} != {expected_shape}")
    for metric in ("cumulative_R_per_mass", "qCO2_days10_15"):
        expected_rank = [
            row["soil"][0] for row in stage23_step6["arm_b_rank_tables"]["primary"][metric]
        ]
        if analysis["base"]["arm_b"]["primary"][f"{metric}_rank"] != expected_rank:
            raise RuntimeError(f"Stage 23 base {metric} rank did not reproduce")
    base = analysis["base"]
    base_summary = {
        "sand_vertisol_distinct": base["arm_a"]["pairwise_distances_and_distinctness"]["S-V"]["distinct"],
        "peak_height_order": base["arm_a"]["peak_height_order"],
        "peak_time_order": base["arm_a"]["peak_time_order"],
        "arm_b_shape_primary": [base["arm_b"]["primary"][soil]["shape_stage23_mean_curve"] for soil in ("S", "V", "L")],
        "arm_b_cumulative_rank_primary": base["arm_b"]["primary"]["cumulative_R_per_mass_rank"],
    }
    variants_summary = {}
    for variant, row in analysis.items():
        current = {
            "sand_vertisol_distinct": row["arm_a"]["pairwise_distances_and_distinctness"]["S-V"]["distinct"],
            "peak_height_order": row["arm_a"]["peak_height_order"],
            "peak_time_order": row["arm_a"]["peak_time_order"],
            "arm_b_shape_primary": [row["arm_b"]["primary"][soil]["shape_stage23_mean_curve"] for soil in ("S", "V", "L")],
            "arm_b_cumulative_rank_primary": row["arm_b"]["primary"]["cumulative_R_per_mass_rank"],
        }
        variants_summary[variant] = {
            "outputs": current,
            "same_as_base": {key: current[key] == base_summary[key] for key in current},
        }
    result = {
        "stage": 24,
        "step": "A2_analysis",
        "variants": analysis,
        "variants_x_outputs": variants_summary,
        "base_outputs": base_summary,
        "base_stage23_analysis_reproduction": True,
        "figures": ["sens_curves.png", "sens_summary.png"],
    }
    # Large internal arrays are kept out of the report artifact; the curves remain in NPZ files.
    public_variants = {}
    for variant, row in analysis.items():
        public_arm_a = {
            key: value
            for key, value in row["arm_a"].items()
            if key not in ("mean_curves_normalized_to_floor", "crop_runs", "crop_arrays")
        }
        public_variants[variant] = {
            "arm_a": public_arm_a,
            "arm_b": row["arm_b"],
            "delivery_split_arm_b": row["delivery_split_arm_b"],
            "unchanged_soils_source": row["unchanged_soils_source"],
        }
    serializable = {
        **{key: value for key, value in result.items() if key != "variants"},
        "variants": public_variants,
    }
    _write_json(OUT / "analysis.json", serializable)
    _plot_sensitivity_curves(analysis)
    _plot_sensitivity_summary(analysis, geometry)
    print(json.dumps({"variants_x_outputs": variants_summary, "base_reproduction": True}, indent=2))
    return serializable


def _correlation_pore_low_grey(pore: np.ndarray, grey: np.ndarray) -> float:
    if pore.shape != grey.shape:
        raise ValueError(f"Alignment shapes differ: {pore.shape} and {grey.shape}")
    n = 0
    sx = sy = sxx = syy = sxy = 0.0
    for z0 in range(0, grey.shape[2], 8):
        z1 = min(z0 + 8, grey.shape[2])
        x = pore[:, :, z0:z1].astype(np.float64, copy=False)
        y = -grey[:, :, z0:z1].astype(np.float64, copy=False)
        n += x.size
        sx += float(x.sum())
        sy += float(y.sum())
        sxx += float(np.square(x).sum())
        syy += float(np.square(y).sum())
        sxy += float((x * y).sum())
    cov = sxy - sx * sy / n
    vx = sxx - sx * sx / n
    vy = syy - sy * sy / n
    if vx <= 0 or vy <= 0:
        return 0.0
    return float(cov / np.sqrt(vx * vy))


def _slice_corr(pore: np.ndarray, grey: np.ndarray, axis: int, index: int) -> float:
    pore_slice = np.take(pore, index, axis=axis).astype(np.float64)
    grey_slice = -np.take(grey, index, axis=axis).astype(np.float64)
    if np.std(pore_slice) == 0 or np.std(grey_slice) == 0:
        return 0.0
    return float(np.corrcoef(pore_slice.ravel(), grey_slice.ravel())[0, 1])


def _aligned_scan(scan_id: str) -> tuple[np.ndarray, np.ndarray, dict]:
    import nibabel as nib

    raw_image = nib.load(str(ROOT / RAW_PATHS[scan_id]))
    raw = np.asanyarray(raw_image.dataobj)
    if scan_id in ("SA", "SB"):
        segmented = np.load(S23.ROOT / S23.INPUT_PATHS[scan_id], mmap_mode="r")
        best = None
        slice_index = min(325, raw.shape[2] - 1)
        for permutation in itertools.permutations(range(3)):
            candidate = np.transpose(segmented, permutation)
            if candidate.shape != raw.shape:
                continue
            score = abs(_slice_corr(candidate > 0, raw, 2, slice_index))
            if best is None or score > best[0]:
                best = (score, permutation, candidate)
        if best is None:
            raise ValueError(f"No NPY axis permutation aligns {scan_id}: {segmented.shape} vs {raw.shape}")
        aligned = best[2]
        pore = np.asarray(aligned > 0)
        alignment = {
            "method": "make_fig7.py permutation search, maximize absolute pore/low-grey slice correlation",
            "axis_permutation": list(best[1]),
            "reference_slice_axis2": slice_index,
            "slice_correlation_pore_low_grey": float(best[0]),
            "raw_shape": list(raw.shape),
            "segmentation_shape": list(segmented.shape),
        }
    else:
        segmented_image = nib.load(str(S23.ROOT / S23.INPUT_PATHS[scan_id]))
        segmented = np.asanyarray(segmented_image.dataobj)
        aligned = np.flip(segmented, axis=1)
        if aligned.shape != raw.shape:
            raise ValueError(f"Y-flipped segmentation shape mismatch for {scan_id}: {aligned.shape} vs {raw.shape}")
        pore = aligned == 5
        alignment = {
            "method": "make_fig7.py y-flip of nnU-Net segmentation",
            "axis_flip": 1,
            "raw_shape": list(raw.shape),
            "segmentation_shape": list(segmented.shape),
        }
    if raw.ndim != 3 or pore.shape != raw.shape:
        raise ValueError(f"Expected voxel-aligned 3-D arrays for {scan_id}")
    alignment["correlation_pore_low_grey"] = _correlation_pore_low_grey(pore, raw)
    alignment["gate_pass"] = alignment["correlation_pore_low_grey"] >= 0.5
    return raw, np.asarray(aligned), alignment


def _quantile_summary(values: np.ndarray) -> dict:
    if values.size == 0:
        raise ValueError("Cannot summarize an empty grey-level set")
    sample_limit = 500_000
    flattened = values.ravel()
    if flattened.size > sample_limit:
        sample = flattened[np.random.default_rng(24).choice(flattened.size, sample_limit, replace=False)]
    else:
        sample = flattened
    quantiles = np.quantile(sample, [0.01, 0.25, 0.5, 0.75, 0.95, 0.99])
    return {
        "n": int(values.size),
        "summary_sample_n": int(sample.size),
        "mean": float(np.mean(sample, dtype=np.float64)),
        "std": float(np.std(sample, dtype=np.float64)),
        "quantiles_sampled": bool(values.size > sample_limit),
        "q01": float(quantiles[0]),
        "q25": float(quantiles[1]),
        "median": float(quantiles[2]),
        "q75": float(quantiles[3]),
        "q95": float(quantiles[4]),
        "q99": float(quantiles[5]),
    }


def _gmm_test(values: np.ndarray) -> dict:
    from sklearn.mixture import GaussianMixture

    rng = np.random.default_rng(24)
    sample_limit = 500_000
    flattened = values.ravel()
    if flattened.size > sample_limit:
        sample = flattened[rng.choice(flattened.size, sample_limit, replace=False)]
    else:
        sample = flattened
    sample = sample.astype(np.float64, copy=False)
    if sample.size < 100:
        raise ValueError(f"Too few eroded voxels for GMM: {sample.size}")
    center = float(sample.mean())
    scale = float(sample.std())
    if scale <= 0:
        raise ValueError("Eroded grey-level sample has zero variance")
    standardized = ((sample - center) / scale).reshape(-1, 1)
    one = GaussianMixture(n_components=1, covariance_type="full", n_init=3, random_state=24, reg_covar=1e-6)
    two = GaussianMixture(n_components=2, covariance_type="full", n_init=5, random_state=24, reg_covar=1e-6)
    one.fit(standardized)
    two.fit(standardized)
    means = two.means_.ravel() * scale + center
    variances = two.covariances_.reshape(-1) * scale**2
    order = np.argsort(means)
    means, variances, weights = means[order], variances[order], two.weights_[order]
    ashman_d = float(np.sqrt(2.0) * abs(means[1] - means[0]) / np.sqrt(variances.sum()))
    delta_bic = float(one.bic(standardized) - two.bic(standardized))
    minor_weight = float(weights.min())
    qualifies = delta_bic > 10.0 and ashman_d > 2.0 and minor_weight >= 0.05
    return {
        "sample_n": int(sample.size),
        "sample_limit": sample_limit,
        "BIC_1_component": float(one.bic(standardized)),
        "BIC_2_component": float(two.bic(standardized)),
        "delta_BIC_1_minus_2": delta_bic,
        "Ashman_D": ashman_d,
        "means_grey": means.tolist(),
        "variances_grey": variances.tolist(),
        "weights_dark_to_bright": weights.tolist(),
        "minor_component_weight": minor_weight,
        "two_populations": bool(qualifies),
        "rule": "delta_BIC > 10, Ashman_D > 2, minor weight >= 0.05",
    }, two, (center, scale)


def _multiotsu_thresholds(grey: np.ndarray) -> tuple[list[float], np.ndarray, np.ndarray]:
    from skimage.filters import threshold_multiotsu

    values = np.asarray(grey)
    low = float(np.min(values))
    high = float(np.max(values))
    if low == high:
        raise ValueError("Raw scan has constant grey values")
    if np.issubdtype(values.dtype, np.integer) and high - low <= 65535:
        counts = np.zeros(int(high - low + 1), dtype=np.int64)
        for z0 in range(0, values.shape[2], 8):
            chunk = values[:, :, z0 : z0 + 8]
            shifted = chunk.astype(np.int64, copy=False) - int(low)
            counts += np.bincount(shifted.ravel(), minlength=counts.size)
        centers = np.arange(int(low), int(high) + 1, dtype=float)
    else:
        counts, edges = np.histogram(values, bins=4096, range=(low, high))
        centers = (edges[:-1] + edges[1:]) / 2.0
    thresholds = threshold_multiotsu(hist=(counts, centers), classes=3)
    return [float(value) for value in thresholds], counts, centers


def _histogram_panel(ax, whole_values: np.ndarray, raw_values: np.ndarray, eroded_values: np.ndarray, label1: np.ndarray | None, label1_eroded: np.ndarray | None, thresholds: list[float], gmm: dict | None, model, scale_info) -> None:
    high = float(np.quantile(raw_values, 0.999))
    low = float(np.quantile(raw_values, 0.001))
    bins = np.linspace(low, high, 240)
    ax.hist(whole_values, bins=bins, density=True, histtype="step", color="#c7c7c7", linewidth=0.8, label="whole volume")
    ax.hist(raw_values, bins=bins, density=True, histtype="step", color="#6e7d8b", linewidth=1.0, label="phase raw")
    ax.hist(eroded_values, bins=bins, density=True, histtype="step", color="#c47b36", linewidth=1.2, label="phase eroded")
    if label1 is not None and label1.size:
        ax.hist(label1, bins=bins, density=True, histtype="step", color="#7b5ca6", linewidth=0.9, linestyle="--", label="label 1 raw")
    if label1_eroded is not None and label1_eroded.size:
        ax.hist(label1_eroded, bins=bins, density=True, histtype="step", color="#7b5ca6", linewidth=1.1, linestyle=":", label="label 1 eroded")
    for index, threshold in enumerate(thresholds):
        ax.axvline(threshold, color="#168a8a", linestyle=":", linewidth=1, label="Multi-Otsu" if index == 0 else None)
    if gmm is not None and model is not None:
        center, scale = scale_info
        x = np.linspace(low, high, 500)
        standardized = (x - center) / scale
        weights = model.weights_
        means = model.means_.ravel()
        variances = model.covariances_.reshape(-1)
        density = np.zeros_like(x)
        for weight, mean, variance in zip(weights, means, variances):
            component = weight * np.exp(-0.5 * (standardized - mean) ** 2 / variance) / np.sqrt(2 * np.pi * variance) / scale
            density += component
            ax.plot(x, component, color="#26734d", linewidth=0.8, alpha=0.8)
        ax.plot(x, density, color="#26734d", linewidth=1.2, label="GMM-2 total")
    ax.set_xlim(low, high)
    ax.set_title(ax.get_title(), fontsize=9)
    ax.grid(alpha=0.2)


def step5() -> dict:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.ndimage import binary_erosion

    if not (OUT / "preflight.json").is_file():
        raise FileNotFoundError("Run Stage 24 Step 0 before grey-level analysis")
    reports = {}
    figure, axes = plt.subplots(1, len(SEGMENTATION_RAW_IDS), figsize=(20, 4.5))
    for axis, scan_id in zip(axes, SEGMENTATION_RAW_IDS):
        raw, aligned, alignment = _aligned_scan(scan_id)
        if not alignment["gate_pass"]:
            reports[scan_id] = {
                "alignment": alignment,
                "analysis_stopped": True,
                "reason": "pore-mask correlation with low grey values was below 0.5",
            }
            axis.set_title(f"{scan_id}: alignment gate failed")
            axis.axis("off")
            continue
        border = 10
        if any(size <= 2 * border for size in raw.shape):
            raise ValueError(f"{scan_id} is too small for the registered 10-voxel border")
        core = tuple(slice(border, -border) for _ in range(3))
        grey = raw[core]
        labels = aligned[core]
        if scan_id in ("SA", "SB"):
            phase_mask = labels == 0
            label1_mask = None
        else:
            phase_mask = np.isin(labels, (0, 1, 2, 3))
            label1_mask = labels == 1
        eroded = binary_erosion(phase_mask, iterations=2, border_value=0)
        raw_phase_values = np.asarray(grey[phase_mask])
        eroded_values = np.asarray(grey[eroded])
        label1_values = np.asarray(grey[label1_mask]) if label1_mask is not None else None
        label1_eroded = (
            binary_erosion(label1_mask, iterations=2, border_value=0)
            if label1_mask is not None
            else None
        )
        label1_eroded_values = (
            np.asarray(grey[label1_eroded])
            if label1_eroded is not None and label1_eroded.any()
            else None
        )
        gmm_report, gmm_model, scale_info = _gmm_test(eroded_values)
        thresholds, otsu_counts, otsu_centers = _multiotsu_thresholds(raw)
        middle = (grey > thresholds[0]) & (grey <= thresholds[1])
        middle_fraction = float(np.count_nonzero(middle & phase_mask) / np.count_nonzero(phase_mask))
        middle_eroded_fraction = float(np.count_nonzero(middle & eroded) / np.count_nonzero(eroded))
        report = {
            "alignment": alignment,
            "analysis_stopped": False,
            "core_shape_after_10_voxel_border": list(grey.shape),
            "histograms": {
                "whole_volume": _quantile_summary(raw),
                "phase_raw": _quantile_summary(raw_phase_values),
                "phase_eroded_2vox": _quantile_summary(eroded_values),
                "label1_raw": _quantile_summary(label1_values) if label1_values is not None and label1_values.size else None,
                "label1_eroded_2vox": (
                    _quantile_summary(label1_eroded_values)
                    if label1_eroded_values is not None
                    else None
                ),
            },
            "GMM": gmm_report,
            "multiotsu_3class": {
                "thresholds": thresholds,
                "histogram_bins_used": int(otsu_counts.size),
                "nonpore_middle_class_fraction_before_erosion": middle_fraction,
                "nonpore_middle_class_fraction_after_2vox_erosion": middle_eroded_fraction,
            },
        }
        if scan_id in ("SA", "SB"):
            sample = raw_phase_values.ravel()
            if sample.size > 500_000:
                sample = sample[np.random.default_rng(24).choice(sample.size, 500_000, replace=False)]
            q50, q90, q99, q999 = np.quantile(sample, [0.5, 0.9, 0.99, 0.999])
            report["sand_bright_tail"] = {
                "solid_q50": float(q50),
                "solid_q90": float(q90),
                "solid_q99": float(q99),
                "solid_q99_9": float(q999),
                "top_0_1_percent_range": [float(q999), float(np.max(raw_phase_values))],
                "interpretation": "A high-grey tail is descriptive only; it is not specific evidence of Fe-oxide coatings.",
            }
        if gmm_report["two_populations"]:
            means = gmm_report["means_grey"]
            threshold = 0.5 * (means[0] + means[1])
            darker = raw_phase_values < threshold
            report["grey_split_preview"] = {
                "dark_component_is_matrix": True,
                "decision_threshold_grey": threshold,
                "matrix_fraction_of_eroded_phase": float(darker.mean()),
                "matrix_fraction_of_core_crop": float(np.count_nonzero(phase_mask & (grey < threshold)) / phase_mask.size),
                "preview_png": str(OUT / f"greysplit_{scan_id}.png"),
            }
            mid = raw.shape[2] // 2
            core_z = mid - border
            image = raw[border:-border, border:-border, mid].T
            mask_slice = phase_mask[:, :, core_z].T & (grey[:, :, core_z].T < threshold)
            fig, ax = plt.subplots(figsize=(6, 6))
            ax.imshow(image, cmap="gray", interpolation="nearest")
            overlay = np.zeros((*mask_slice.shape, 4), dtype=float)
            overlay[mask_slice] = (1.0, 0.64, 0.0, 0.55)
            ax.imshow(overlay, interpolation="nearest")
            ax.set_title(f"{scan_id}: darker GMM population = candidate matrix")
            ax.set_axis_off()
            fig.tight_layout()
            fig.savefig(OUT / f"greysplit_{scan_id}.png", dpi=150)
            plt.close(fig)
        reports[scan_id] = report
        axis.set_title(scan_id)
        axis.set_xlabel("Grey value")
        axis.set_ylabel("Density")
        raw_plot = raw_phase_values
        eroded_plot = eroded_values
        _histogram_panel(axis, raw.ravel(), raw_plot, eroded_plot, label1_values, label1_eroded_values, thresholds, gmm_report, gmm_model, scale_info)
        axis.legend(fontsize=6)
        del raw, aligned, grey, labels
    figure.tight_layout()
    figure.savefig(OUT / "grey_histograms.png", dpi=150)
    plt.close(figure)
    result = {
        "stage": 24,
        "step": "B_grey_levels",
        "scans": reports,
        "alignment_helper": str(ALIGNMENT_HELPER),
        "figure": str(OUT / "grey_histograms.png"),
        "two_population_rule": "delta BIC > 10, Ashman D > 2, minor weight >= 0.05",
    }
    _write_json(OUT / "grey_levels.json", result)
    print(json.dumps({
        scan: {
            "correlation": report["alignment"]["correlation_pore_low_grey"],
            "two_populations": report.get("GMM", {}).get("two_populations"),
            "grey_split_preview": report.get("grey_split_preview"),
        }
        for scan, report in reports.items()
    }, indent=2))
    return result


def _percent(value: float | None) -> str:
    return "n/a" if value is None else f"{100 * value:.1f}%"


def step6() -> str:
    required = (
        OUT / "preflight.json",
        OUT / "geometry.json",
        OUT / "reproduction_gate.json",
        OUT / "runs.json",
        OUT / "analysis.json",
        OUT / "grey_levels.json",
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Stage 24 report inputs are missing: {missing}")
    preflight = _read_json(required[0])
    geometry = _read_json(required[1])
    reproduction = _read_json(required[2])
    run_summary = _read_json(required[3])
    analysis = _read_json(required[4])
    grey = _read_json(required[5])
    variants = [row.name for row in MS.VARIANTS]

    def curve_source(variant: str, soil: str) -> str:
        if variant != "base" and MS.variant_changes(VARIANT_BY_NAME[variant], soil):
            return "Stage 24"
        return "Stage 23 reused"

    geometry_rows = []
    for variant in variants:
        for soil in ("S", "V", "L"):
            row = geometry["aggregates"][variant][soil]
            geometry_rows.append(
                f"| {variant} | {SOIL_NAMES[soil]} | {row['matrix_fraction']['mean']:.4f} ± {row['matrix_fraction']['sd']:.4f} | "
                f"{row['solid_fraction']['mean']:.4f} ± {row['solid_fraction']['sd']:.4f} | "
                f"{row['pore_fraction']['mean']:.4f} ± {row['pore_fraction']['sd']:.4f} | "
                f"{row['matrix_pore_interface_faces_per_crop_voxel']['mean']:.4f} | "
                f"{row['mean_matrix_distance_to_pore_um']['mean']:.2f} µm | "
                f"{row['om_per_matrix_voxel_A']['mean']:.6g} | {row['om_per_matrix_voxel_B']['mean']:.6g} |"
            )
    summary_rows = []
    for variant in variants:
        row = analysis["variants_x_outputs"][variant]
        outputs, same = row["outputs"], row["same_as_base"]
        summary_rows.append(
            f"| {variant} | {outputs['sand_vertisol_distinct']} ({'same' if same['sand_vertisol_distinct'] else 'different'}) | "
            f"{' > '.join(outputs['peak_height_order'])} ({'same' if same['peak_height_order'] else 'different'}) | "
            f"{' < '.join(outputs['peak_time_order'])} ({'same' if same['peak_time_order'] else 'different'}) | "
            f"{'/'.join(outputs['arm_b_shape_primary'])} ({'same' if same['arm_b_shape_primary'] else 'different'}) | "
            f"{' > '.join(outputs['arm_b_cumulative_rank_primary'])} ({'same' if same['arm_b_cumulative_rank_primary'] else 'different'}) | "
            f"Stage 24 affected soils; Stage 23 reused otherwise |"
        )
    pair_rows = []
    for variant in variants:
        for pair, row in analysis["variants"][variant]["arm_a"]["pairwise_distances_and_distinctness"].items():
            one, full = row["by_steps"]["1500"], row["by_steps"]["3000"]
            soil_a, soil_b = pair.split("-")
            sources = f"{soil_a}: {curve_source(variant, soil_a)}; {soil_b}: {curve_source(variant, soil_b)}"
            pair_rows.append(
                f"| {variant} | {pair} | {one['soil_mean_distance']:.4f} | {one['max_within']:.4f} | "
                f"{one['min_cross']:.4f} | {row['distinct']} | {full['soil_mean_distance']:.4f} | {sources} |"
            )
    period_rows = []
    rank_rows = []
    for variant in variants:
        for mapping in ("low", "primary", "high"):
            arm_b = analysis["variants"][variant]["arm_b"][mapping]
            for soil in ("S", "V", "L"):
                row = arm_b[soil]
                periods = ", ".join(
                    f"{mean['mean']:.5g} ± {mean['sd']:.3g}"
                    for mean in row["period_R_per_voxel_mean_sd"]
                )
                period_rows.append(
                    f"| {variant} | {mapping} | {SOIL_NAMES[soil]} | {periods} | {row['shape_stage23_mean_curve']} | "
                    f"{curve_source(variant, soil)} |"
                )
            qco2_values = ", ".join(
                f"{soil} {arm_b[soil]['qCO2_days10_15']['mean']:.4g} ± "
                f"{arm_b[soil]['qCO2_days10_15']['sd']:.3g}"
                for soil in ("S", "V", "L")
            )
            rank_rows.append(
                f"| {variant} | {mapping} | {' > '.join(arm_b['cumulative_R_per_mass_rank'])} | "
                f"{qco2_values} | {' > '.join(arm_b['qCO2_days10_15_rank'])} | "
                f"Stage 24 affected soils; Stage 23 reused otherwise |"
            )
    delivery_rows = []
    for variant in variants:
        for soil in ("S", "V", "L"):
            row = analysis["variants"][variant]["delivery_split_arm_b"][soil]
            total = sum(
                row[name]
                for name in (
                    "diffusion_into_wet_habitat",
                    "diffusion_into_drained_habitat",
                    "exchange_into_wet_habitat",
                    "exchange_into_drained_habitat",
                )
            )
            delivery_rows.append(
                f"| {variant} | {SOIL_NAMES[soil]} | "
                f"{_percent(row['diffusion_into_wet_habitat'] / total if total else None)} | "
                f"{_percent(row['diffusion_into_drained_habitat'] / total if total else None)} | "
                f"{_percent(row['exchange_into_wet_habitat'] / total if total else None)} | "
                f"{_percent(row['exchange_into_drained_habitat'] / total if total else None)} | "
                f"{_percent(row['diffusion_fraction_total_delivery'])} | "
                f"{_percent(row['exchange_fraction_total_delivery'])} | "
                f"{_percent(row['drained_exchange_fraction_of_exchange'])} | {curve_source(variant, soil)} |"
            )
    grey_rows = []
    for scan_id in SEGMENTATION_RAW_IDS:
        row = grey["scans"][scan_id]
        align = row["alignment"]
        if row.get("analysis_stopped"):
            grey_rows.append(
                f"| {scan_id} | {align['correlation_pore_low_grey']:.3f} | stopped (alignment < 0.5) | — | — | — | — | — |"
            )
            continue
        gmm = row["GMM"]
        otsu = row["multiotsu_3class"]
        preview = row.get("grey_split_preview", {})
        hist = row["histograms"]
        grey_rows.append(
            f"| {scan_id} | {align['correlation_pore_low_grey']:.3f} | "
            f"ΔBIC {gmm['delta_BIC_1_minus_2']:.2f}; D {gmm['Ashman_D']:.2f}; "
            f"minor {gmm['minor_component_weight']:.3f}; two populations: {gmm['two_populations']} | "
            f"{otsu['thresholds'][0]:.4g}, {otsu['thresholds'][1]:.4g} | "
            f"{otsu['nonpore_middle_class_fraction_before_erosion']:.3f} / "
            f"{otsu['nonpore_middle_class_fraction_after_2vox_erosion']:.3f} | "
            f"{hist['phase_raw']['n']} / {hist['phase_eroded_2vox']['n']} | "
            f"{preview.get('matrix_fraction_of_core_crop', 'no qualifying split')} | "
            f"{'saved' if preview else 'not generated'} |"
        )

    base = analysis["variants"]["base"]
    sv_all = [analysis["variants"][variant]["arm_a"]["pairwise_distances_and_distinctness"]["S-V"]["distinct"] for variant in variants]
    distinct_answer = "passes all variants" if all(sv_all) else "flips in " + ", ".join(
        variant for variant, distinct in zip(variants, sv_all) if distinct != sv_all[0]
    )
    height_flips = [variant for variant in variants if analysis["variants_x_outputs"][variant]["same_as_base"]["peak_height_order"] is False]
    time_flips = [variant for variant in variants if analysis["variants_x_outputs"][variant]["same_as_base"]["peak_time_order"] is False]
    shape_flips = [variant for variant in variants if analysis["variants_x_outputs"][variant]["same_as_base"]["arm_b_shape_primary"] is False]
    rank_flips = [variant for variant in variants if analysis["variants_x_outputs"][variant]["same_as_base"]["arm_b_cumulative_rank_primary"] is False]
    sand_drained = {
        variant: analysis["variants"][variant]["delivery_split_arm_b"]["S"]["drained_exchange_fraction_total_delivery"]
        for variant in variants
    }
    sand_change = (
        sand_drained["S30"] - sand_drained["base"]
        if "S30" in sand_drained else 0.0
    )
    verdict = (
        f"Arm A Sand–Vertisol distinctness {distinct_answer}. Peak-height order differs from base in "
        f"{', '.join(height_flips) if height_flips else 'no variant'}; peak-time order differs in "
        f"{', '.join(time_flips) if time_flips else 'no variant'}. Arm B primary shape classes differ in "
        f"{', '.join(shape_flips) if shape_flips else 'no variant'}; cumulative-R/mass rank differs in "
        f"{', '.join(rank_flips) if rank_flips else 'no variant'}. Sand drained-exchange delivery share "
        f"changes from {_percent(sand_drained['base'])} in base to {_percent(sand_drained['S30'])} in S30 "
        f"({100*sand_change:+.1f} percentage points)."
    )
    elapsed_by_variant = {}
    for run in run_summary["runs"].values():
        name = run["variant"]
        elapsed_by_variant.setdefault(name, []).append(float(run["seconds"]))
    runtime_rows = [
        f"| {name} | {len(values)} | {np.mean(values):.1f} | {max(values):.1f} | {sum(values):.1f} |"
        for name, values in elapsed_by_variant.items()
    ]
    runtime_rows.insert(
        0,
        f"| base reproduction | {reproduction['run_count']} | "
        f"{np.mean([row['seconds'] for row in reproduction['runs'].values()]):.1f} | "
        f"{max(row['seconds'] for row in reproduction['runs'].values()):.1f} | "
        f"{sum(row['seconds'] for row in reproduction['runs'].values()):.1f} |",
    )
    alignment_helper_hash = preflight["input_gates"]["alignment_helper_sha256"]
    report = f"""# Stage 24 — matrix/inert-solid sensitivity and raw grey levels

## Pre-flight

- Biology gate: passed; normalized SHA-256 `{preflight['biology_gate']['sha256_normalized']}`, matching HEAD, main and origin/main.
- The requested initial command outputs were empty:
  - `git --no-pager diff HEAD -- configs/biology.yaml`: *(no output)*
  - `git status --short`: *(no output)*
- The working copy used CRLF line endings; after CRLF-to-LF normalization its bytes matched the committed config. `git stash push -u -m "pre-stage24-leftovers"` returned **No local changes to save**. No change was stashed or used; the apparent difference was only line-ending representation.
- The committed file history was `ce82d87 first commit`; HEAD and main both normalize to the expected hash. The branch was created from fetched `origin/main`.
- All seven Stage 23 segmentations, five raw scans and fixed Stage 23 run artifacts were present. The alignment helper copy used was `{ALIGNMENT_HELPER}` (SHA-256 `{alignment_helper_hash}`); it is byte-identical to the handoff copy of the referenced `make_fig7.py`.
- Stage 23 crops/origins, chemistry values, `K_density={K_DENSITY}`, `T={T_STEPS}`, the Stage 23 shape classifier, `r_lit=0.36/2.6/9.2`, and biology were unchanged. No lab files were opened.

## Reproduction gate

The 14 base runs (A and B, c0 for each of the seven specimens) were rerun and compared pointwise with Stage 23 R(t). The gate is maximum relative difference `<1e-8`; mass balance must remain `≤1e-10`.

| Runs | Max relative difference | Max mass-balance error | Gate |
|---:|---:|---:|---|
| {reproduction['run_count']} | {reproduction['max_relative_difference']:.3g} | {max(row['mass_balance_rel'] for row in reproduction['runs'].values()):.3g} | {'PASS' if reproduction['passed'] else 'FAIL'} |

## Part A — geometry first, then respiration

Every listed geometry is a mean ± SD across the fixed Stage 23 crop set for that soil. “Solid” is inert SOLID; the Sand target is a fraction of crop volume. Vertisol/Loess inert shares are fractions of the original matrix voxels converted to solid. Unchanged soils in sensitivity variants reuse Stage 23 curves.

| Variant | Soil | Matrix fraction | Inert solid fraction | Pore fraction | Matrix–pore faces / N | Mean matrix-to-pore distance | OM/matrix voxel A | OM/matrix voxel B |
|---|---|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(geometry_rows)}

### Predictions versus outcomes

| Prediction | Observed outcome |
|---|---|
| Sand peak height decreases and time becomes later from S04 to S30. | Sand Arm A height {analysis['variants']['S04']['arm_a']['soil_means']['S']['peak_height_R_peak_over_floor']:.4g} → {analysis['variants']['S30']['arm_a']['soil_means']['S']['peak_height_R_peak_over_floor']:.4g}; peak-time fraction {analysis['variants']['S04']['arm_a']['soil_means']['S']['peak_time_fraction']:.4g} → {analysis['variants']['S30']['arm_a']['soil_means']['S']['peak_time_fraction']:.4g}. |
| Vertisol/Loess peak height increases and time becomes earlier from base through VL25 to VL50. | Vertisol height {[analysis['variants'][name]['arm_a']['soil_means']['V']['peak_height_R_peak_over_floor'] for name in ('base','VL25','VL50')]} and peak-time fraction {[analysis['variants'][name]['arm_a']['soil_means']['V']['peak_time_fraction'] for name in ('base','VL25','VL50')]}; Loess height {[analysis['variants'][name]['arm_a']['soil_means']['L']['peak_height_R_peak_over_floor'] for name in ('base','VL25','VL50')]} and peak-time fraction {[analysis['variants'][name]['arm_a']['soil_means']['L']['peak_time_fraction'] for name in ('base','VL25','VL50')]}. |
| Sand–Vertisol remain distinct in CONV. | CONV Sand–Vertisol distinct: {analysis['variants']['CONV']['arm_a']['pairwise_distances_and_distinctness']['S-V']['distinct']}. |

### Arm A: distances and distinctness

The Stage 23 rule applies at 0–1500 steps: both soils at least 75% alive, soil-mean distance >0.3, and max within-soil distance < half the min cross-soil distance.

| Variant | Pair | Mean distance 0–1500 | Max within | Min cross | Distinct? | Mean distance 0–3000 | Curve source by soil |
|---|---|---:|---:|---:|---|---:|---|
{chr(10).join(pair_rows)}

| Variant | Peak-height order (high to low) | Peak-time order (early to late) |
|---|---|---|
{chr(10).join(f"| {variant} | {' > '.join(analysis['variants'][variant]['arm_a']['peak_height_order'])} | {' < '.join(analysis['variants'][variant]['arm_a']['peak_time_order'])} |" for variant in variants)}

### Arm B: periods, shape, cumulative rank, qCO₂

P1–P4 values are mean ± SD of per-voxel R, in model units. Shape classes use the unchanged Stage 23 classifier.

| Variant | r_lit | Soil | P1 | P2 | P3 | P4 | Shape | Curve source |
|---|---|---|---:|---:|---:|---:|---|---|
{chr(10).join(period_rows)}

| Variant | r_lit | Cumulative R/mass rank (high to low) | qCO₂ days 10–15 (S/V/L mean ± SD) | qCO₂ rank (high to low) | Curve source |
|---|---|---|---|---|---|
{chr(10).join(rank_rows)}

### Variants × outputs

| Variant | Sand–Vertisol distinct? | Peak-height order | Peak-time order | Arm B shapes at 2.6 (S/V/L) | Arm B cumulative rank at 2.6 | Curve source |
|---|---|---|---|---|---|---|
{chr(10).join(summary_rows)}

`same` / `different` compares the listed output with base. The primary mapping is `r_lit=2.6`.

### Delivery split (Arm B)

The first two delivery columns partition total net delivery. Wet-exchange is its share of total delivery; drained-exchange columns show Sand’s drained-pore exchange reliance explicitly.

| Variant | Soil | Diffusion wet / total | Diffusion drained / total | Exchange wet / total | Exchange drained / total | Total diffusion | Total exchange | Drained exchange / exchange | Curve source |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
{chr(10).join(delivery_rows)}

Sand’s drained-exchange share is {_percent(sand_drained['base'])} in base versus {_percent(sand_drained['S30'])} in S30.

## Part B — raw grey-level split evidence

Grey histograms cover the whole volume and the core after excluding a 10-voxel border; phase masks and label 1 are shown raw and after 2-voxel erosion. Multi-Otsu thresholds are fit to the whole-volume histogram. The fixed two-population rule is applied to the eroded phase sample. Preview fractions are candidate darker-population volume fractions, not adopted parameters.

| Scan | Corr(pore, low grey) | Two-population test (ΔBIC, Ashman D, minor wt, pass) | Multi-Otsu thresholds | Middle non-pore fraction before / after erosion | Phase voxels raw / eroded | Candidate matrix fraction of core | Overlay |
|---|---:|---|---|---|---:|---:|---|
{chr(10).join(grey_rows)}

For Sand, the solid histograms show the bright-tail quantiles in `grey_levels.json`; a bright tail can be consistent with dense/coated grains but is not specific evidence for Fe-oxide coatings. No scan failing the 0.5 alignment threshold was analyzed further.

## Verdict

{verdict}

CONV: Sand–Vertisol {'remain' if analysis['variants']['CONV']['arm_a']['pairwise_distances_and_distinctness']['S-V']['distinct'] else 'do not remain'} distinct even though their matrix fractions are similar.

The grey-level verdict is determined by the fixed GMM rule: {', '.join(f"{scan}={'two populations' if row.get('GMM', {}).get('two_populations') else 'no qualifying two-population split' if not row.get('analysis_stopped') else 'alignment stopped'}" for scan, row in grey['scans'].items())}. A split is only plausible where all three pre-registered GMM criteria pass; it is not used in the model.

## Limitations

- Matrix fractions are sensitivity inputs, not measured quantities; no variant is adopted.
- Crop replicates are subsamples of the seven registered scans, not independent field replicates.
- The dark/bright interpretation assumes darker material is porous matrix; grey-level populations may also reflect composition, density, reconstruction, noise or partial volume.
- Multi-Otsu is an intensity partition, not a mineralogical or porosity assay; erosion cannot eliminate every partial-volume effect.
- Stage 23 chemistry means and the fixed time mapping were not varied. No laboratory respiration data were opened.
- Sand’s possible bright tail is not diagnostic of Fe-oxide coatings.
- Stage 24 changes matrix amount while preserving total OM, intentionally changing OM per matrix voxel.

## Runtime and delivery log

| Run group | Runs | Mean seconds | Maximum seconds | Sum seconds |
|---|---:|---:|---:|---:|
{chr(10).join(runtime_rows)}

Figures: [sens_curves.png](../../results/stage24/sens_curves.png), [sens_summary.png](../../results/stage24/sens_summary.png), [grey_histograms.png](../../results/stage24/grey_histograms.png). Preview overlays are saved for qualifying scans only.

Branch: `agents/stage24-matrix-sensitivity`. Commits and PR: recorded in `LOGBOOK.md` after publication.
"""
    RESULTS_DOC.parent.mkdir(parents=True, exist_ok=True)
    if RESULTS_DOC.exists():
        raise FileExistsError(f"Refusing to overwrite Stage 24 report: {RESULTS_DOC}")
    RESULTS_DOC.write_text(report, encoding="utf-8")
    print(f"Wrote {RESULTS_DOC}")
    return report


def main(argv: list[str] | None = None) -> None:
    args = list(sys.argv[1:] if argv is None else argv)
    if args and args[0] == "job":
        if len(args) != 2:
            raise SystemExit("Usage: python stage24_run.py job <run-key>")
        _worker(args[1])
        return
    steps = {
        "step0": step0,
        "step1": step1,
        "step2": step2,
        "step3": step3,
        "step4": step4,
        "step5": step5,
        "step6": step6,
    }
    if len(args) != 1 or args[0] not in steps:
        raise SystemExit("Usage: python stage24_run.py step0|step1|step2|step3|step4|step5|step6")
    steps[args[0]]()


if __name__ == "__main__":
    main()
