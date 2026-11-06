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

import metrics
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


def _run_key(structure_id: str) -> str:
    with (OUT / "step1.json").open(encoding="utf-8") as file:
        report = json.load(file)["reports"][structure_id]
    return f"{structure_id}_v{_voxel_tag(report['vox_um'])}um_T3000"


def _load_run(structure_id: str) -> dict:
    key = _run_key(structure_id)
    with (RUNS / f"{key}.json").open(encoding="utf-8") as file:
        result = json.load(file)
    with np.load(RUNS / f"{key}.npz", allow_pickle=False) as arrays:
        result.update({name: arrays[name] for name in arrays.files})
    return result


def _normalized_curve(run: dict, m0: float) -> np.ndarray:
    denominator = m0 * float(run["B"][0])
    if denominator <= 0:
        raise ValueError(f"Nonpositive m0*B0 for {run['key']}")
    return np.asarray(run["R"], dtype=float) / denominator


def _distance(curve_a: np.ndarray, curve_b: np.ndarray, steps: int) -> float:
    return metrics.pairwise_distance(curve_a[:steps], curve_b[:steps])


def _mean_curve(curves: list[np.ndarray]) -> np.ndarray:
    if not curves:
        raise ValueError("Cannot average an empty curve set")
    return np.mean(np.stack(curves), axis=0)


def step3() -> dict:
    with (OUT / "step2.json").open(encoding="utf-8") as file:
        step2_result = json.load(file)
    bio = load_yaml("configs/biology.yaml")
    curve_ids = {
        "sand": [f"{scan_id}_c{index}" for scan_id in ("SA", "SB") for index in range(4)],
        "vertisol": [f"{scan_id}_c{index}" for scan_id in ("VA", "VB") for index in range(4)],
    }
    runs = {structure_id: _load_run(structure_id) for ids in curve_ids.values() for structure_id in ids}
    curves = {key: _normalized_curve(value, bio["m0"]) for key, value in runs.items()}
    means = {soil: _mean_curve([curves[key] for key in ids]) for soil, ids in curve_ids.items()}
    within = {}
    cross = {}
    for steps in (1500, 3000):
        within_by_soil = {}
        within_pairs = {}
        for soil, ids in curve_ids.items():
            groups = {"within_specimen": [], "between_specimen": []}
            all_distances = []
            for index_a in range(len(ids)):
                for index_b in range(index_a + 1, len(ids)):
                    a, b = ids[index_a], ids[index_b]
                    pair = {
                        "a": a,
                        "b": b,
                        "distance": _distance(curves[a], curves[b], steps),
                    }
                    all_distances.append(pair)
                    group = "within_specimen" if a.split("_")[0] == b.split("_")[0] else "between_specimen"
                    groups[group].append(pair)
            within_by_soil[soil] = {
                group: [pair["distance"] for pair in pairs] for group, pairs in groups.items()
            }
            within_pairs[soil] = all_distances
        cross_pairs = []
        for sand_id in curve_ids["sand"]:
            for vertisol_id in curve_ids["vertisol"]:
                cross_pairs.append({
                    "sand": sand_id,
                    "vertisol": vertisol_id,
                    "distance": _distance(curves[sand_id], curves[vertisol_id], steps),
                })
        all_within = [pair["distance"] for pairs in within_pairs.values() for pair in pairs]
        cross_distances = [pair["distance"] for pair in cross_pairs]
        within[steps] = {
            "by_soil": within_by_soil,
            "pairs": within_pairs,
            "max_all": float(max(all_within)),
            "p90_all": float(np.percentile(all_within, 90)),
        }
        cross[steps] = {
            "pairs": cross_pairs,
            "min": float(min(cross_distances)),
            "max": float(max(cross_distances)),
        }

    alive = {
        soil: {
            "count": sum(bool(runs[key]["alive_primary"]) for key in ids),
            "total": len(ids),
            "crops": {key: bool(runs[key]["alive_primary"]) for key in ids},
        }
        for soil, ids in curve_ids.items()
    }
    mean_distances = {
        str(steps): _distance(means["sand"], means["vertisol"], steps)
        for steps in (1500, 3000)
    }
    criteria = {
        "alive_at_least_6_each": all(alive[soil]["count"] >= 6 for soil in alive),
        "soil_mean_distance_01500_gt_0_3": mean_distances["1500"] > 0.3,
        "max_within_lt_half_min_cross_01500": within[1500]["max_all"] < 0.5 * cross[1500]["min"],
    }
    result = {
        "alive_primary": alive,
        "soil_mean_distance": mean_distances,
        "within_soil": {str(k): v for k, v in within.items()},
        "cross_soil": {str(k): v for k, v in cross.items()},
        "secondary_p90_within_vs_min_cross_01500": {
            "p90_within": within[1500]["p90_all"],
            "min_cross": cross[1500]["min"],
            "p90_lt_half_min_cross": within[1500]["p90_all"] < 0.5 * cross[1500]["min"],
        },
        "criteria": criteria,
        "pass": all(criteria.values()),
        "active_step2": "step2.json",
        "K_density": step2_result["K_density"],
    }
    curves_path = OUT / "step3_soil_mean_curves.npz"
    if curves_path.exists():
        raise FileExistsError(f"Refusing to overwrite Stage 22 output: {curves_path}")
    np.savez_compressed(curves_path, sand=means["sand"], vertisol=means["vertisol"])
    _write_json(OUT / "step3.json", result)
    print(json.dumps({
        "alive_primary": alive,
        "soil_mean_distance": mean_distances,
        "max_within_1500": within[1500]["max_all"],
        "p90_within_1500": within[1500]["p90_all"],
        "min_cross_1500": cross[1500]["min"],
        "criteria": criteria,
        "pass": result["pass"],
    }, indent=2))
    return result


def _control_mean(kind: str, m0: float) -> tuple[np.ndarray, list[np.ndarray]]:
    curves = [_normalized_curve(_load_run(f"{kind}15_s{seed}"), m0) for seed in range(3)]
    return _mean_curve(curves), curves


def _stage21_mean(kind: str, m0: float) -> tuple[np.ndarray, list[np.ndarray]]:
    curves = []
    for seed in range(3):
        run = S21.load_run(S21.jkey(kind, seed, 1.0))
        curves.append(np.asarray(run["R"], dtype=float) / (m0 * float(run["B"][0])))
    return _mean_curve(curves), curves


def step4() -> dict:
    with (OUT / "step2.json").open(encoding="utf-8") as file:
        step2_result = json.load(file)
    bio = load_yaml("configs/biology.yaml")
    m0 = float(bio["m0"])
    scan_ids = {
        "sand": [f"{scan_id}_c{index}" for scan_id in ("SA", "SB") for index in range(4)],
        "vertisol": [f"{scan_id}_c{index}" for scan_id in ("VA", "VB") for index in range(4)],
    }
    scan_curves = {
        soil: [_normalized_curve(_load_run(structure_id), m0) for structure_id in ids]
        for soil, ids in scan_ids.items()
    }
    scan_means = {soil: _mean_curve(curves) for soil, curves in scan_curves.items()}
    controls = {kind: _control_mean(kind, m0) for kind in ("S", "V")}
    synthetic5 = {kind: _stage21_mean(kind, m0) for kind in ("S", "V")}
    distances = {}
    for steps in (1500, 3000):
        distances[str(steps)] = {
            "S_scan_vs_S15": _distance(scan_means["sand"], controls["S"][0], steps),
            "V_scan_vs_V15": _distance(scan_means["vertisol"], controls["V"][0], steps),
            "S15_vs_V15": _distance(controls["S"][0], controls["V"][0], steps),
            "S_scan_vs_V_scan": _distance(scan_means["sand"], scan_means["vertisol"], steps),
        }
    distances["1500"].update({
        "S15_vs_S5": _distance(controls["S"][0], synthetic5["S"][0], 1500),
        "V15_vs_V5": _distance(controls["V"][0], synthetic5["V"][0], 1500),
    })
    control_noise = {
        kind: [
            _distance(curves[a], curves[b], 3000)
            for a in range(3) for b in range(a + 1, 3)
        ]
        for kind, (_, curves) in controls.items()
    }
    domain_comparisons = {}
    for scan_id in SCAN_IDS:
        specimen_curves = [
            _normalized_curve(_load_run(f"{scan_id}_c{index}"), m0) for index in range(4)
        ]
        specimen_mean = _mean_curve(specimen_curves)
        domain_curve = _normalized_curve(_load_run(f"{scan_id}_domain256"), m0)
        domain_comparisons[scan_id] = {
            "distance_to_100cube_specimen_mean_03000": _distance(domain_curve, specimen_mean, 3000),
        }
    source_effect = float(np.mean([
        distances["3000"]["S_scan_vs_S15"],
        distances["3000"]["V_scan_vs_V15"],
    ]))
    resolution_effect = float(np.mean([
        distances["1500"]["S15_vs_S5"],
        distances["1500"]["V15_vs_V5"],
    ]))
    result = {
        "distances": distances,
        "control_seed_noise_03000": control_noise,
        "domain_size_comparisons": domain_comparisons,
        "factorial_2x2": {
            "synthetic_5um_S_vs_V": _distance(synthetic5["S"][0], synthetic5["V"][0], 1500),
            "synthetic_15um_S_vs_V": distances["1500"]["S15_vs_V15"],
            "scan_15um_S_vs_V": distances["1500"]["S_scan_vs_V_scan"],
            "mean_structure_source_effect_at_15um": source_effect,
            "mean_resolution_effect_in_synthetic": resolution_effect,
            "larger_mean_curve_shift": "structure source" if source_effect > resolution_effect else "resolution",
        },
        "active_step2": "step2.json",
        "K_density": step2_result["K_density"],
    }
    curve_path = OUT / "step4_curves.png"
    if curve_path.exists():
        raise FileExistsError(f"Refusing to overwrite Stage 22 figure: {curve_path}")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    x = np.arange(3000)
    colors = {"sand": "#c45127", "vertisol": "#4b398a"}
    for soil, curves in scan_curves.items():
        for curve in curves:
            axes[0].plot(x, curve, color=colors[soil], alpha=0.16, linewidth=0.8)
        axes[0].plot(x, scan_means[soil], color=colors[soil], linewidth=2.2, label=f"{soil} scan mean")
        axes[1].plot(x, scan_means[soil], color=colors[soil], linewidth=2.2, label=f"{soil} scan")
    for kind, label in (("S", "S15 control"), ("V", "V15 control")):
        color = colors["sand" if kind == "S" else "vertisol"]
        axes[1].plot(x, controls[kind][0], color=color, linestyle="--", linewidth=1.6, label=label)
        axes[1].plot(
            np.arange(len(synthetic5[kind][0])), synthetic5[kind][0],
            color=color, linestyle=":", linewidth=1.5, label=f"{kind}5 Stage 21",
        )
    for ax in axes:
        ax.set_xlabel("model step")
        ax.set_ylabel("R/(m0 B0)")
        ax.legend(fontsize=8)
    axes[0].set_title("Scan crops and soil means")
    axes[1].set_title("Soil means and synthetic controls")
    fig.tight_layout()
    fig.savefig(curve_path, dpi=130)
    plt.close(fig)
    _write_json(OUT / "step4.json", result)
    print(json.dumps(result, indent=2))
    return result


def _load_run_key(key: str) -> dict:
    with (RUNS / f"{key}.json").open(encoding="utf-8") as file:
        result = json.load(file)
    with np.load(RUNS / f"{key}.npz", allow_pickle=False) as arrays:
        result.update({name: arrays[name] for name in arrays.files})
    return result


def _peak_order(sand: dict, vertisol: dict) -> tuple[str, str]:
    amplitude = "Sand" if sand["R_peak"] > vertisol["R_peak"] else (
        "Vertisol" if vertisol["R_peak"] > sand["R_peak"] else "tie"
    )
    sand_step = int(np.argmax(sand["R"]))
    vertisol_step = int(np.argmax(vertisol["R"]))
    timing = "Sand earlier" if sand_step < vertisol_step else (
        "Vertisol earlier" if vertisol_step < sand_step else "tie"
    )
    return amplitude, timing


def step5() -> dict:
    with (OUT / "step2.json").open(encoding="utf-8") as file:
        step2_result = json.load(file)
    with (OUT / "step3.json").open(encoding="utf-8") as file:
        step3_result = json.load(file)
    if not step3_result["pass"]:
        raise RuntimeError("Step 5 is gated on a passing Step 3")
    bio = load_yaml("configs/biology.yaml")
    soil_ids = {
        "sand": [f"{scan_id}_c{index}" for scan_id in ("SA", "SB") for index in range(4)],
        "vertisol": [f"{scan_id}_c{index}" for scan_id in ("VA", "VB") for index in range(4)],
    }
    mean_curves = {}
    with np.load(OUT / "step3_soil_mean_curves.npz", allow_pickle=False) as arrays:
        mean_curves = {soil: arrays[soil] for soil in soil_ids}
    representatives = {}
    for soil, ids in soil_ids.items():
        candidates = []
        for structure_id in ids:
            run = _load_run(structure_id)
            curve = _normalized_curve(run, bio["m0"])
            candidates.append((metrics.pairwise_distance(curve, mean_curves[soil]), structure_id, run))
        _, selected_id, selected_run = min(candidates, key=lambda value: value[0])
        representatives[soil] = {
            "structure_id": selected_id,
            "distance_to_soil_mean": float(min(candidates, key=lambda value: value[0])[0]),
            "run": selected_run,
        }

    hp = M.load_inherited()
    parameters = {
        "r_max": ("biology", float(bio["r_max"])),
        "Ks": ("biology", float(bio["Ks"])),
        "Y": ("biology", float(bio["Y"])),
        "m0": ("biology", float(bio["m0"])),
        "m_s": ("biology", float(bio["m_s"])),
        "k_dis": ("biology", float(bio["k_dis"])),
        "k_ex": ("model", float(hp["k_ex"])),
        "K_density": ("capacity", float(step2_result["K_density"])),
    }
    keys = {}
    base_key = {}
    for parameter, (target, value) in parameters.items():
        for factor in (0.5, 1.5):
            perturbation = value * factor
            for soil in ("sand", "vertisol"):
                structure_id = representatives[soil]["structure_id"]
                key = f"sens_{parameter}_x{factor}_{structure_id}_v{_voxel_tag(VOX_UM[structure_id[:2]])}um_T3000"
                base_key[(parameter, factor, soil)] = key
                biology_overrides, model_overrides = {}, {}
                K_density = step2_result["K_density"]
                if target == "biology":
                    biology_overrides[parameter] = perturbation
                elif target == "model":
                    model_overrides[parameter] = perturbation
                else:
                    K_density = perturbation
                _write_json(_job_config_path(key), {
                    "structure_id": structure_id,
                    "T": 3000,
                    "OM_density": 1.0,
                    "K_density": K_density,
                    "biology_overrides": biology_overrides,
                    "model_overrides": model_overrides,
                })
                keys[(parameter, factor, soil)] = key

    results = run_jobs(list(keys.values()), ngpu=2)
    output = {"representative_crops": {}, "perturbations": {}}
    base_orders = _peak_order(representatives["sand"]["run"], representatives["vertisol"]["run"])
    for soil in ("sand", "vertisol"):
        output["representative_crops"][soil] = {
            "structure_id": representatives[soil]["structure_id"],
            "distance_to_soil_mean": representatives[soil]["distance_to_soil_mean"],
        }
    for parameter, (_, value) in parameters.items():
        output["perturbations"][parameter] = {}
        for factor in (0.5, 1.5):
            sand_key = keys[(parameter, factor, "sand")]
            vertisol_key = keys[(parameter, factor, "vertisol")]
            sand_run, vertisol_run = results[sand_key], results[vertisol_key]
            sand_arrays, vertisol_arrays = _load_run_key(sand_key), _load_run_key(vertisol_key)
            orders = _peak_order(sand_arrays, vertisol_arrays)
            output["perturbations"][parameter][str(factor)] = {
                "base_value": value,
                "perturbed_value": value * factor,
                "alive_primary_both": bool(sand_run["alive_primary"] and vertisol_run["alive_primary"]),
                "alive_primary": {
                    "sand": bool(sand_run["alive_primary"]),
                    "vertisol": bool(vertisol_run["alive_primary"]),
                },
                "distance_01500": _distance(sand_arrays["R"], vertisol_arrays["R"], 1500),
                "peak_amplitude_order": {"base": base_orders[0], "perturbed": orders[0],
                                         "unchanged": base_orders[0] == orders[0]},
                "peak_time_order": {"base": base_orders[1], "perturbed": orders[1],
                                    "unchanged": base_orders[1] == orders[1]},
                "run_keys": {"sand": sand_key, "vertisol": vertisol_key},
            }
    output["n_runs"] = len(keys)
    output["active_step2"] = "step2.json"
    _write_json(OUT / "step5.json", output)
    print(json.dumps({
        "representative_crops": output["representative_crops"],
        "n_runs": output["n_runs"],
        "alive_both_counts": {
            parameter: sum(
                item["alive_primary_both"]
                for item in perturbations.values()
            )
            for parameter, perturbations in output["perturbations"].items()
        },
    }, indent=2))
    return output


if __name__ == "__main__":
    command = sys.argv[1]
    if command == "job":
        job(sys.argv[2])
    else:
        globals()[command]()
