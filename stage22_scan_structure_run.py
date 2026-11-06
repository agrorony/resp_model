"""Stage 22 scan-structure workflow.

Run geometry before respiration:
    python stage22_scan_structure_run.py step1
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

import model3d as M
import packing3d as P
import scan3d
import stage21_run as S21
import timemap
from model import load_yaml


ROOT = Path(os.environ.get("RB", r"E:\PROJECTS\Yael_Mishael\Rony\remote_computer backup"))
OUT = Path("results") / "stage22_scan_structure"
STRUCTURES = OUT / "structures"
RUNS = OUT / "runs"
JOB_CONFIGS = OUT / "job_configs"
SCAN_IDS = ("VA", "VB", "SA", "SB")
VOX_UM = {"VA": 15.000149, "VB": 15.034357, "SA": 15.000149, "SB": 15.034357}
GEOMETRY_METRICS = (
    "porosity_resolved",
    "matrix_fraction",
    "solid_fraction",
    "largest_cluster_fraction",
    "S_w",
    "habitat_volume_fraction",
    "wet_fraction_of_habitat",
    "spans_axes_xyz",
)


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
        raise FileExistsError(f"Refusing to overwrite Stage 22 output: {path}")
    path.write_text(json.dumps(value, indent=2, default=_json_default) + "\n", encoding="utf-8")


def _save_structure(structure_id: str, structure, report: dict) -> None:
    STRUCTURES.mkdir(parents=True, exist_ok=True)
    path = STRUCTURES / f"{structure_id}.npz"
    if path.exists():
        raise FileExistsError(f"Refusing to reuse a structure from an earlier run: {path}")
    np.savez_compressed(
        path,
        phase=P.to_np(structure.phase),
        lt=P.to_np(structure.lt),
        vox_um=np.float64(structure.vox_um),
        seed=np.int64(structure.seed),
    )
    report["structure_file"] = str(path)


def _load_structure(structure_id: str, xp):
    path = STRUCTURES / f"{structure_id}.npz"
    with np.load(path, allow_pickle=False) as data:
        structure = P.Structure()
        structure.kind = structure_id
        structure.phase = xp.asarray(data["phase"])
        structure.lt = xp.asarray(data["lt"])
        structure.n = int(data["phase"].shape[0])
        structure.vox_um = float(data["vox_um"])
        structure.seed = int(data["seed"])
        structure.info = {"structure_id": structure_id}
    return structure


def _analyse_and_save(structure_id: str, structure, xp) -> dict:
    report = P.analyse(structure, xp)
    report["structure_id"] = structure_id
    _save_structure(structure_id, structure, report)
    return report


def _ensure_structure_cache(step1: dict) -> None:
    for scan_id in SCAN_IDS:
        selected_origins = step1["positions"][scan_id]["selected"]
        structure_ids = [f"{scan_id}_c{index}" for index in range(len(selected_origins))]
        domain_id = f"{scan_id}_domain{step1['positions'][scan_id]['domain_actual_size']}"
        if not any(not (STRUCTURES / f"{structure_id}.npz").exists()
                   for structure_id in [*structure_ids, domain_id]):
            continue
        volume = scan3d.load_scan(scan_id, ROOT)
        xp = P.get_xp()
        for structure_id in structure_ids:
            path = STRUCTURES / f"{structure_id}.npz"
            if path.exists():
                continue
            report1 = step1["reports"][structure_id]
            crop = scan3d.crop_at(volume, tuple(report1["origin"]), report1["n"])
            structure, _ = scan3d.crop_to_structure(crop, scan_id, VOX_UM[scan_id], xp=xp)
            structure.phase = xp.asarray(structure.phase)
            _analyse_and_save(structure_id, structure, xp)
        path = STRUCTURES / f"{domain_id}.npz"
        if not path.exists():
            report1 = step1["reports"][domain_id]
            crop = scan3d.crop_at(volume, tuple(report1["origin"]), report1["actual_size"])
            structure, _ = scan3d.crop_to_structure(crop, scan_id, VOX_UM[scan_id], xp=xp)
            structure.phase = xp.asarray(structure.phase)
            _analyse_and_save(domain_id, structure, xp)
        del volume
    xp = P.get_xp()
    for kind in ("S", "V"):
        for seed in range(3):
            structure_id = f"{kind}15_s{seed}"
            if (STRUCTURES / f"{structure_id}.npz").exists():
                continue
            structure = P.build_structure(kind, seed_offset=seed, n=100, vox_um=15.0, xp=xp)
            _analyse_and_save(structure_id, structure, xp)


def _job_config_path(key: str) -> Path:
    return JOB_CONFIGS / f"{key}.json"


def _voxel_tag(vox_um: float) -> str:
    return f"{vox_um:.6f}".rstrip("0").rstrip(".")


def job(key: str) -> dict:
    with _job_config_path(key).open(encoding="utf-8") as file:
        config = json.load(file)
    xp = P.get_xp()
    structure = _load_structure(config["structure_id"], xp)
    hp = M.load_inherited()
    hp.update(config.get("model_overrides", {}))
    bio = load_yaml("configs/biology.yaml")
    bio.update(config.get("biology_overrides", {}))
    K_density = config.get("K_density", S21.load("step1.json")["K_density"])
    start = time.time()
    out, summary = M.run_structure(
        structure, hp, bio, xp, K_density, config["OM_density"], T=config["T"]
    )
    elapsed = time.time() - start
    if summary["mass_balance_rel"] > 1e-10:
        raise RuntimeError(f"{key} mass balance failed: {summary['mass_balance_rel']}")
    if structure.n == 100 and elapsed > 600:
        raise TimeoutError(f"{key} exceeded the 10 minute 100-cube run limit: {elapsed:.1f}s")
    B0 = hp["B0_density"] * structure.phase.size
    summary.update({
        "key": key,
        "structure_id": config["structure_id"],
        "T": config["T"],
        "OM_density": config["OM_density"],
        "K_density": K_density,
        "vox_um": structure.vox_um,
        "n": structure.n,
        "seconds": elapsed,
        "R_peak": float(np.max(P.to_np(out["R"]))),
        "R_floor": float(bio["m0"] * B0),
        "backend": xp.__name__,
    })
    RUNS.mkdir(parents=True, exist_ok=True)
    json_path = RUNS / f"{key}.json"
    npz_path = RUNS / f"{key}.npz"
    if json_path.exists() or npz_path.exists():
        raise FileExistsError(f"Refusing to overwrite Stage 22 run: {key}")
    arrays = {name: P.to_np(out[name]) for name in M.STAT_NAMES}
    np.savez_compressed(npz_path, **arrays)
    _write_json(json_path, summary)
    return summary


def run_jobs(keys: list[str], ngpu: int = 2) -> dict[str, dict]:
    RUNS.mkdir(parents=True, exist_ok=True)
    JOB_CONFIGS.mkdir(parents=True, exist_ok=True)
    running = {}
    pending = list(keys)
    results = {}
    while pending or running:
        for slot in range(ngpu):
            if slot in running or not pending:
                continue
            key = pending.pop(0)
            json_path, npz_path = RUNS / f"{key}.json", RUNS / f"{key}.npz"
            if json_path.exists() or npz_path.exists():
                raise FileExistsError(f"Refusing to reuse prior Stage 22 output: {key}")
            config_path = _job_config_path(key)
            if not config_path.exists():
                raise FileNotFoundError(f"Missing job configuration for {key}")
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
            structure_id = "_".join(key.split("_")[:-2])
            size = int(structure_id.rsplit("domain", 1)[1]) if "_domain" in structure_id else 100
            if size == 100 and process.poll() is None and time.time() - started > 600:
                process.terminate()
                process.wait()
                raise TimeoutError(f"{key} exceeded the 10 minute 100-cube run limit")
            if process.poll() is None:
                continue
            stderr = process.stderr.read() if process.stderr else ""
            if process.returncode:
                for other_slot, (_, other, _) in list(running.items()):
                    if other_slot != slot and other.poll() is None:
                        other.terminate()
                raise RuntimeError(f"Stage 22 run {key} failed: {stderr[-4000:]}")
            with (RUNS / f"{key}.json").open(encoding="utf-8") as file:
                results[key] = json.load(file)
            print(f"done {key} ({results[key]['seconds']:.1f}s)", flush=True)
            del running[slot]
    return results


def step2() -> dict:
    with (OUT / "step1.json").open(encoding="utf-8") as file:
        step1_result = json.load(file)
    _ensure_structure_cache(step1_result)
    keys = []
    for scan_id in SCAN_IDS:
        voxel_tag = _voxel_tag(VOX_UM[scan_id])
        keys.extend(f"{scan_id}_c{crop}_v{voxel_tag}um_T3000" for crop in range(4))
    for scan_id in SCAN_IDS:
        size = step1_result["positions"][scan_id]["domain_actual_size"]
        voxel_tag = _voxel_tag(VOX_UM[scan_id])
        keys.append(f"{scan_id}_domain{size}_v{voxel_tag}um_T3000")
    for kind in ("S", "V"):
        keys.extend(f"{kind}15_s{seed}_v15um_T3000" for seed in range(3))
    S21.check_biology()
    K_density = S21.load("step1.json")["K_density"]
    for key in keys:
        parts = key.split("_")
        structure_id = "_".join(parts[:-2])
        _write_json(_job_config_path(key), {
            "structure_id": structure_id,
            "T": 3000,
            "OM_density": 1.0,
            "K_density": K_density,
            "biology_overrides": {},
            "model_overrides": {},
        })
    results = run_jobs(keys, ngpu=2)
    output = {"K_density": K_density, "T": 3000, "OM_density": 1.0, "runs": results}
    _write_json(OUT / "step2.json", output)
    return output


def _aggregate(reports: dict[str, dict], ids: list[str]) -> dict:
    result = {"n": len(ids)}
    for metric in GEOMETRY_METRICS:
        values = [reports[structure_id][metric] for structure_id in ids]
        if metric == "spans_axes_xyz":
            result[metric] = [
                {"mean": float(np.mean([v[axis] for v in values])),
                 "sd": float(np.std([v[axis] for v in values], ddof=1)) if len(values) > 1 else 0.0}
                for axis in range(3)
            ]
        else:
            numbers = np.asarray(values, dtype=float)
            result[metric] = {
                "mean": float(numbers.mean()),
                "sd": float(numbers.std(ddof=1)) if len(numbers) > 1 else 0.0,
            }
    for percentile in ("d10", "d50", "d90"):
        values = [reports[structure_id]["LT_percentiles_um"][percentile] for structure_id in ids]
        result[f"LT_{percentile}_um"] = {
            "mean": float(np.mean(values)),
            "sd": float(np.std(values, ddof=1)) if len(values) > 1 else 0.0,
        }
    return result


def _write_step1_slices(reports: dict[str, dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    selected = ("VA_c0", "VB_c0", "SA_c0", "SB_c0", "S15_s0", "V15_s0")
    colors = ("#2786b8", "#d9a441", "#a77746", "#333333")
    cmap = ListedColormap(colors)
    fig, axes = plt.subplots(1, len(selected), figsize=(15, 3.8))
    for ax, structure_id in zip(axes, selected):
        structure = _load_structure(structure_id, np)
        phase = structure.phase[:, :, structure.n // 2]
        lt = structure.lt[:, :, structure.n // 2]
        image = np.full(phase.shape, 3, dtype=np.uint8)
        image[phase == P.MATRIX] = 2
        image[(phase == P.PORE) & (lt < 50.0)] = 0
        image[(phase == P.PORE) & (lt >= 50.0)] = 1
        ax.imshow(image.T, cmap=cmap, vmin=0, vmax=3, interpolation="nearest", origin="lower")
        ax.set_title(structure_id, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    legend = [
        Patch(color=colors[0], label="wet pore (<50 um)"),
        Patch(color=colors[1], label="drained pore (>=50 um)"),
        Patch(color=colors[2], label="matrix"),
        Patch(color=colors[3], label="inert solid"),
    ]
    fig.legend(handles=legend, loc="lower center", ncol=4, fontsize=8)
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    target = OUT / "step1_slices.png"
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite Stage 22 figure: {target}")
    fig.savefig(target, dpi=130)
    plt.close(fig)


def step1() -> dict:
    xp = P.get_xp()
    reports: dict[str, dict] = {}
    crops_by_scan: dict[str, list[str]] = {}
    positions: dict[str, dict] = {}

    for scan_id in SCAN_IDS:
        volume = scan3d.load_scan(scan_id, ROOT)
        gate_report = scan3d.scan_histogram(scan_id, volume)
        if not scan3d.scan_gate(scan_id, gate_report):
            raise RuntimeError(f"Input gate failed for {scan_id}: {gate_report}")
        main_origins = scan3d.grid_origins(volume.shape, 100)
        chosen, rejected = scan3d.draw_habitable_crops(
            volume, scan_id, VOX_UM[scan_id], main_origins, crop_size=100, n_crops=4, seed=22
        )
        positions[scan_id] = {
            "grid_origins_100": main_origins,
            "selected": chosen,
            "rejected_zero_habitat": rejected,
        }
        crops_by_scan[scan_id] = []
        for crop_index, origin in enumerate(chosen):
            structure_id = f"{scan_id}_c{crop_index}"
            crop = scan3d.crop_at(volume, origin, 100)
            structure, threshold = scan3d.crop_to_structure(crop, scan_id, VOX_UM[scan_id], xp=xp)
            structure.phase = xp.asarray(structure.phase)
            report = _analyse_and_save(structure_id, structure, xp)
            report["origin"] = origin
            report["solid_LT_threshold_um"] = threshold
            reports[structure_id] = report
            crops_by_scan[scan_id].append(structure_id)

        domain_origin, domain_size = scan3d.centered_cube_origin(volume.shape, 256)
        domain_id = f"{scan_id}_domain{domain_size}"
        crop = scan3d.crop_at(volume, domain_origin, domain_size)
        structure, threshold = scan3d.crop_to_structure(crop, scan_id, VOX_UM[scan_id], xp=xp)
        structure.phase = xp.asarray(structure.phase)
        report = _analyse_and_save(domain_id, structure, xp)
        report["origin"] = domain_origin
        report["requested_size"] = 256
        report["actual_size"] = domain_size
        report["solid_LT_threshold_um"] = threshold
        reports[domain_id] = report
        positions[scan_id]["domain_origin"] = domain_origin
        positions[scan_id]["domain_actual_size"] = domain_size
        del volume

    for kind in ("S", "V"):
        for seed in range(3):
            structure_id = f"{kind}15_s{seed}"
            structure = P.build_structure(kind, seed_offset=seed, n=100, vox_um=15.0, xp=xp)
            reports[structure_id] = _analyse_and_save(structure_id, structure, xp)

    aggregates = {}
    for scan_id in SCAN_IDS:
        aggregates[scan_id] = _aggregate(reports, crops_by_scan[scan_id])
    for soil, scan_ids in (("sand", ("SA", "SB")), ("vertisol", ("VA", "VB"))):
        all_ids = [structure_id for scan_id in scan_ids for structure_id in crops_by_scan[scan_id]]
        aggregates[soil] = _aggregate(reports, all_ids)

    result = {
        "input_gates": {
            scan_id: scan3d.scan_histogram(scan_id, scan3d.load_scan(scan_id, ROOT))
            for scan_id in SCAN_IDS
        },
        "positions": positions,
        "reports": reports,
        "aggregates_main_crops": aggregates,
        "controls": [f"{kind}15_s{seed}" for kind in ("S", "V") for seed in range(3)],
        "domain_crops": [f"{scan_id}_domain{positions[scan_id]['domain_actual_size']}" for scan_id in SCAN_IDS],
        "vertisol_plate_voxels_at_15um": {
            "median": float(P.STRUCTURE_PARAMS["vertisol"]["plate_median_um"]
                            / P.STRUCTURE_PARAMS["vertisol"]["plate_aspect"] / 15.0),
            "note": "Median synthetic plate thickness is the median plate diameter divided by plate aspect ratio and voxel size.",
        },
    }
    _write_json(OUT / "step1.json", result)
    _write_step1_slices(reports)
    print(json.dumps({
        "crop_origins": {scan_id: positions[scan_id]["selected"] for scan_id in SCAN_IDS},
        "zero_habitat_rejections": {scan_id: positions[scan_id]["rejected_zero_habitat"] for scan_id in SCAN_IDS},
        "aggregates_main_crops": aggregates,
        "vertisol_plate_voxels_at_15um": result["vertisol_plate_voxels_at_15um"],
    }, indent=2, default=_json_default))
    return result


if __name__ == "__main__":
    command = sys.argv[1]
    if command == "job":
        job(sys.argv[2])
    else:
        globals()[command]()
