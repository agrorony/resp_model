"""Stage 23 Loess-chemistry workflow. Run geometry before respiration."""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

import chem23
import model3d as M
import packing3d as P
import scan3d
import scan3d_loess
from model import load_yaml


ROOT = Path(os.environ.get("RB", r"E:\PROJECTS\Yael_Mishael\Rony\remote_computer backup"))
OUT = Path("results") / "stage23"
STRUCTURES = OUT / "structures"
RUNS = OUT / "runs"
JOB_CONFIGS = OUT / "job_configs"
SCAN_IDS = ("VA", "VB", "SA", "SB", "LA", "LB", "LC")
SOIL_FOR_SCAN = {
    "SA": "S",
    "SB": "S",
    "VA": "V",
    "VB": "V",
    "LA": "L",
    "LB": "L",
    "LC": "L",
}
VOX_UM = {
    "VA": 15.000149,
    "VB": 15.034357,
    "SA": 15.000149,
    "SB": 15.034357,
    **scan3d_loess.VOX_UM,
}
K_DENSITY = 0.16445
T_STEPS = 3000
DT = 0.05
N_CROPS = 4
CROP_SIZE = 100
R_LIT_DAYS = {"low": 0.36, "primary": 2.6, "high": 9.2}
PERIODS_DAYS = ((0.0, 3.0), (3.0, 7.0), (7.0, 10.0), (10.0, 15.0))

INPUT_PATHS = {
    "VA": Path("nnUNet_resources/bnei_reem_specA_unified/inference_concatenated/bnei_reem_specA_unified.nii.gz"),
    "VB": Path("nnUNet_resources/bnei_reem_specB_unified/inference_concatenated/bnei_reem_specB_unified.nii.gz"),
    "SA": Path("10.5/rehovot_samp_2.npy"),
    "SB": Path("10.5/Rehovot_samp3_highkV_Cu0.11_15um.npy"),
    "LA": Path("nnUNet4SoilXrayCT/05_evaluation/psd/pom_analysis_20260829_roi_expansion/pipeline_work/mishmar_native_5p85um/mishmar_native_5p85um_label_downsample.nii.gz"),
    "LB": Path("nnUNet4SoilXrayCT/05_evaluation/psd/pom_analysis_20260829_roi_expansion/pipeline_work/mishmar_second_8p8um/mishmar_second_8p8um_label_downsample.nii.gz"),
    "LC": Path("nnUNet_resources/mishmar_hanegev_Cu011_samp_2_Rec_nlm/inference_output_concat_loess_i2/mishmar_hanegev_Cu011_samp_2_Rec_nlm.nii.gz"),
}


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
        raise FileExistsError(f"Refusing to overwrite Stage 23 output: {path}")
    path.write_text(json.dumps(value, indent=2, default=_json_default) + "\n", encoding="utf-8")


def _read_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as file:
        return json.load(file)


def _load_scan(scan_id: str) -> np.ndarray:
    if scan_id in scan3d.SOIL_BY_ID:
        return scan3d.load_scan(scan_id, ROOT)
    if scan_id in scan3d_loess.LOESS_IDS:
        return scan3d_loess.load_scan(scan_id, ROOT)
    raise ValueError(f"Unknown scan ID: {scan_id}")


def _save_input_figure(scan_id: str, volume: np.ndarray) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Circle, Patch

    if volume.dtype == np.bool_:
        cmap = ListedColormap(("#333333", "#2786b8"))
        vmin, vmax = 0, 1
        legend = [Patch(color="#333333", label="solid"), Patch(color="#2786b8", label="pore")]
    else:
        cmap = ListedColormap(("#d9d9d9", "#333333", "#a77746", "#bc8f5a", "#ff00ff", "#2786b8"))
        vmin, vmax = 0, 5
        legend = [
            Patch(color="#d9d9d9", label="label 0 / matrix"),
            Patch(color="#333333", label="label 1 / solid"),
            Patch(color="#a77746", label="label 2 / POM"),
            Patch(color="#bc8f5a", label="label 3 / POM"),
            Patch(color="#2786b8", label="label 5 / pore"),
        ]
    x_mid, y_mid, z_mid = (size // 2 for size in volume.shape)
    slices = (
        volume[:, :, z_mid].T,
        volume[:, y_mid, :].T,
        volume[x_mid, :, :].T,
    )
    labels = ("central xy", "central xz", "central yz")
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, image, title in zip(axes, slices, labels):
        ax.imshow(image, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest", origin="lower")
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
    if scan_id == "LB":
        cx, cy, radius = scan3d_loess.inscribed_disk(volume.shape)
        axes[0].add_patch(Circle((cx, cy), radius, fill=False, color="white", linewidth=1.0, linestyle="--"))
        axes[0].add_patch(
            Circle((cx, cy), 0.90 * radius, fill=False, color="#00ffff", linewidth=0.8)
        )
        legend.extend(
            (
                Patch(fill=False, edgecolor="white", linestyle="--", label="inscribed disk"),
                Patch(fill=False, edgecolor="#00ffff", label="crop-safe 0.90 radius"),
            )
        )
    fig.suptitle(f"{scan_id} input labels")
    fig.legend(handles=legend, loc="lower center", ncol=min(5, len(legend)), fontsize=8)
    fig.tight_layout(rect=(0, 0.12, 1, 0.94))
    target = OUT / f"inputs_{scan_id}.png"
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite Stage 23 input figure: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=110)
    plt.close(fig)
    return target


def _verify_biology_matches_main() -> dict:
    current = Path("configs/biology.yaml").read_bytes().replace(b"\r\n", b"\n")
    head = subprocess.check_output(["git", "show", "HEAD:configs/biology.yaml"])
    main = subprocess.check_output(["git", "show", "main:configs/biology.yaml"])
    origin_main = subprocess.check_output(["git", "show", "origin/main:configs/biology.yaml"])
    if current != head or head != main or head != origin_main:
        raise RuntimeError("configs/biology.yaml is not byte-identical to main and origin/main")
    return {
        "sha256": hashlib.sha256(head).hexdigest(),
        "matches_main": True,
        "matches_origin_main": True,
    }


def _reproduce_stage22_sa_c0() -> dict:
    step1 = _read_json(Path("results/stage22_scan_structure/step1.json"))
    origin = tuple(int(value) for value in step1["positions"]["SA"]["selected"][0])
    stored_path = Path(
        "results/stage22_scan_structure/runs/SA_c0_v15.000149um_T3000.npz"
    )
    if not stored_path.is_file():
        raise FileNotFoundError(f"Missing Stage 22 reproduction curve: {stored_path}")
    stored_summary = _read_json(stored_path.with_suffix(".json"))
    if stored_summary["T"] != T_STEPS or stored_summary["OM_density"] != 1.0:
        raise ValueError(f"Unexpected Stage 22 SA_c0 reference metadata: {stored_summary}")
    stage22_k_density = float(stored_summary["K_density"])
    volume = scan3d.load_scan("SA", ROOT)
    crop = scan3d.crop_at(volume, origin, CROP_SIZE)
    if crop.shape != (CROP_SIZE, CROP_SIZE, CROP_SIZE):
        raise ValueError(f"Stage 22 SA_c0 crop has unexpected shape: {crop.shape}")
    xp = P.get_xp()
    structure, threshold = scan3d.crop_to_structure(
        crop, "SA", VOX_UM["SA"], seed=22, xp=xp
    )
    structure.phase = xp.asarray(structure.phase)
    P.analyse(structure, xp)
    hp = M.load_inherited()
    hp["B0_density"] = chem23.MODEL_REFERENCE_B0_DENSITY
    bio = load_yaml("configs/biology.yaml")
    current, summary = M.run_structure(
        structure, hp, bio, xp, stage22_k_density,
        chem23.MODEL_REFERENCE_OM_DENSITY, T=T_STEPS,
    )
    if summary["mass_balance_rel"] > 1e-10:
        raise RuntimeError(
            f"Stage 22 SA_c0 reproduction failed mass balance: {summary['mass_balance_rel']}"
        )
    with np.load(stored_path, allow_pickle=False) as stored:
        reference = np.asarray(stored["R"], dtype=float)
    reproduced = np.asarray(P.to_np(current["R"]), dtype=float)
    if reproduced.shape != reference.shape:
        raise ValueError(
            f"Stage 22 reproduction length mismatch: {reproduced.shape} vs {reference.shape}"
        )
    scale = max(float(np.max(np.abs(reference))), 1e-300)
    pointwise_scale = np.maximum(np.abs(reference), scale * 1e-12)
    max_relative = float(np.max(np.abs(reproduced - reference) / pointwise_scale))
    max_scaled_error = float(np.max(np.abs(reproduced - reference)) / scale)
    result = {
        "reference_run": str(stored_path),
        "structure_id": "SA_c0",
        "origin": origin,
        "shape": list(crop.shape),
        "fines_solid_threshold_um": threshold,
        "chemistry": "chem_eq",
        "K_density": stage22_k_density,
        "T": T_STEPS,
        "max_relative_difference": max_relative,
        "max_difference_over_peak": max_scaled_error,
        "gate_max_relative_lt_1e-8": max_relative < 1e-8,
    }
    if not result["gate_max_relative_lt_1e-8"]:
        raise RuntimeError(f"Stage 22 SA_c0 reproduction gate failed: {result}")
    return result


def _validate_registered_inputs() -> dict:
    missing = [scan_id for scan_id, rel in INPUT_PATHS.items() if not (ROOT / rel).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing registered scan inputs: {missing}")
    reports = {}
    for scan_id in SCAN_IDS:
        volume = _load_scan(scan_id)
        if scan_id in scan3d.SOIL_BY_ID:
            report = scan3d.scan_histogram(scan_id, volume)
            passed = scan3d.scan_gate(scan_id, report)
            report["gate"] = {
                "labels_and_porosity": bool(passed),
                "pass": bool(passed),
            }
        else:
            report = scan3d_loess.scan_histogram(scan_id, volume)
            passed = bool(report["gate"]["pass"])
        if not passed:
            raise RuntimeError(f"Input label/porosity gate failed for {scan_id}: {report}")
        report["input_path"] = str(ROOT / INPUT_PATHS[scan_id])
        report["input_figure"] = str(_save_input_figure(scan_id, volume))
        reports[scan_id] = report
        del volume
    return reports


def step0() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    biology = _verify_biology_matches_main()
    reproduction = _reproduce_stage22_sa_c0()
    input_gates = _validate_registered_inputs()
    chem23.validate_chemistry()
    test_command = [
        sys.executable,
        "-m",
        "unittest",
        "test_stage22_scan_structure",
        "test_stage23",
    ]
    subprocess.run(
        test_command,
        check=True,
        env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"),
    )
    result = {
        "stage": 23,
        "step": 0,
        "main_commit": subprocess.check_output(
            ["git", "rev-parse", "main"], text=True
        ).strip(),
        "biology_config": biology,
        "stage22_reproduction": reproduction,
        "chemistry": chem23.CHEMISTRY,
        "chemistry_validation": "Derived columns agree with prompt values to four significant figures.",
        "input_gates": input_gates,
        "LB_gate_amendment": (
            "The 2026-11-08 amendment gates LB on the whole-volume pore fraction; "
            "the inscribed-disk pore fraction is reported but is not a gate."
        ),
        "tests": "test_stage22_scan_structure and test_stage23 passed.",
    }
    _write_json(OUT / "preflight.json", result)
    print(json.dumps({
        "stage22_reproduction": reproduction,
        "pore_fractions": {
            key: value.get("whole_volume_pore_fraction", value.get("pore_fraction"))
            for key, value in input_gates.items()
        },
        "LB_inside_disk": input_gates["LB"]["inscribed_disk"],
        "biology_config_sha256": biology["sha256"],
        "tests": result["tests"],
    }, indent=2, default=_json_default))
    return result


GEOMETRY_METRICS = (
    "porosity_resolved",
    "matrix_fraction",
    "solid_fraction",
    "largest_cluster_fraction",
    "S_w",
    "habitat_volume_fraction",
    "wet_fraction_of_habitat",
)


def _cache_structure(structure_id: str, structure: P.Structure, report: dict) -> None:
    STRUCTURES.mkdir(parents=True, exist_ok=True)
    path = STRUCTURES / f"{structure_id}.npz"
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite Stage 23 structure: {path}")
    np.savez_compressed(
        path,
        phase=P.to_np(structure.phase),
        lt=P.to_np(structure.lt),
        vox_um=np.float64(structure.vox_um),
        seed=np.int64(structure.seed),
    )
    report["structure_file"] = str(path)


def _analyse_structure(
    structure_id: str,
    structure: P.Structure,
    xp,
    extra: dict | None = None,
) -> dict:
    report = P.analyse(structure, xp)
    report["structure_id"] = structure_id
    if extra:
        report.update(extra)
    _cache_structure(structure_id, structure, report)
    return report


def _load_structure(structure_id: str, xp) -> P.Structure:
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


def _make_structure(crop: np.ndarray, scan_id: str, xp) -> tuple[P.Structure, float | None]:
    if scan_id in scan3d.SOIL_BY_ID:
        structure, threshold = scan3d.crop_to_structure(
            crop, scan_id, VOX_UM[scan_id], seed=22, xp=xp
        )
        structure.phase = xp.asarray(structure.phase)
        return structure, threshold
    structure = scan3d_loess.crop_to_structure(crop, scan_id, xp)
    return structure, None


def _aggregate_geometry(reports: dict[str, dict], structure_ids: list[str]) -> dict:
    result = {"n": len(structure_ids)}
    for metric in GEOMETRY_METRICS:
        values = np.asarray([reports[item][metric] for item in structure_ids], dtype=float)
        result[metric] = {
            "mean": float(values.mean()),
            "sd": float(values.std(ddof=1)) if values.size > 1 else 0.0,
        }
    spans = np.asarray(
        [reports[item]["spans_axes_xyz"] for item in structure_ids], dtype=float
    )
    result["spans_axes_xyz"] = [
        {
            "mean": float(spans[:, axis].mean()),
            "sd": float(spans[:, axis].std(ddof=1)) if spans.shape[0] > 1 else 0.0,
        }
        for axis in range(3)
    ]
    for percentile in ("d10", "d50", "d90"):
        values = np.asarray(
            [reports[item]["LT_percentiles_um"][percentile] for item in structure_ids],
            dtype=float,
        )
        result[f"LT_{percentile}_um"] = {
            "mean": float(values.mean()),
            "sd": float(values.std(ddof=1)) if values.size > 1 else 0.0,
        }
    return result


def _geometry_aggregate_sets(
    crops_by_scan: dict[str, list[str]],
    arm_d_ids: list[str],
    domain_ids: list[str],
) -> dict[str, list[str]]:
    result = {scan_id: list(crops_by_scan[scan_id]) for scan_id in SCAN_IDS}
    result.update(
        {
            "S": crops_by_scan["SA"] + crops_by_scan["SB"],
            "V": crops_by_scan["VA"] + crops_by_scan["VB"],
            "L": crops_by_scan["LA"] + crops_by_scan["LB"] + crops_by_scan["LC"],
            "SA_armD": [item for item in arm_d_ids if item.startswith("SA_")],
            "SB_armD": [item for item in arm_d_ids if item.startswith("SB_")],
            "S_armD": list(arm_d_ids),
            "LA_domain": ["LA_domain256"],
            "LB_domain": ["LB_domain256"],
            "LC_domain": ["LC_domain256"],
            "L_domain": list(domain_ids),
        }
    )
    return result


def _write_step1_slices() -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    selected_ids = [f"{scan_id}_c0" for scan_id in SCAN_IDS]
    colors = ("#2786b8", "#d9a441", "#a77746", "#333333")
    cmap = ListedColormap(colors)
    fig, axes = plt.subplots(1, len(selected_ids), figsize=(17, 3.8))
    for ax, structure_id in zip(axes, selected_ids):
        structure = _load_structure(structure_id, np)
        phase = structure.phase[:, :, structure.n // 2]
        lt = structure.lt[:, :, structure.n // 2]
        image = np.full(phase.shape, 3, dtype=np.uint8)
        image[phase == P.MATRIX] = 2
        image[(phase == P.PORE) & (lt < 50.0)] = 0
        image[(phase == P.PORE) & (lt >= 50.0)] = 1
        ax.imshow(image.T, cmap=cmap, vmin=0, vmax=3, interpolation="nearest", origin="lower")
        ax.set_title(structure_id, fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.legend(
        handles=[
            Patch(color=colors[0], label="wet pore (<50 um)"),
            Patch(color=colors[1], label="drained pore (>=50 um)"),
            Patch(color=colors[2], label="matrix"),
            Patch(color=colors[3], label="inert solid"),
        ],
        loc="lower center",
        ncol=4,
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    target = OUT / "step1_slices.png"
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite Stage 23 figure: {target}")
    fig.savefig(target, dpi=130)
    plt.close(fig)
    return target


def _format_mean_sd(summary: dict, metric: str, digits: int = 3) -> str:
    values = summary[metric]
    return f"{values['mean']:.{digits}f} ± {values['sd']:.{digits}f}"


def _write_predictions(step1_result: dict, overwrite: bool = False) -> Path:
    soils = step1_result["aggregates"]
    table = []
    for soil_code, label in (("S", "Sand"), ("V", "Vertisol"), ("L", "Loess")):
        row = soils[soil_code]
        spans = "/".join(
            f"{item['mean']:.2f}" for item in row["spans_axes_xyz"]
        )
        table.append(
            f"| {label} | {_format_mean_sd(row, 'porosity_resolved')} | "
            f"{_format_mean_sd(row, 'matrix_fraction')} | "
            f"{_format_mean_sd(row, 'solid_fraction')} | "
            f"{_format_mean_sd(row, 'LT_d10_um', 1)} | "
            f"{_format_mean_sd(row, 'LT_d50_um', 1)} | "
            f"{_format_mean_sd(row, 'LT_d90_um', 1)} | "
            f"{_format_mean_sd(row, 'largest_cluster_fraction')} | "
            f"{spans} | "
            f"{_format_mean_sd(row, 'S_w')} | "
            f"{_format_mean_sd(row, 'habitat_volume_fraction')} | "
            f"{_format_mean_sd(row, 'wet_fraction_of_habitat')} |"
        )
    loess = soils["L"]
    sand = soils["S"]
    vertisol = soils["V"]
    loess_wet = loess["wet_fraction_of_habitat"]["mean"]
    sand_wet = sand["wet_fraction_of_habitat"]["mean"]
    vertisol_wet = vertisol["wet_fraction_of_habitat"]["mean"]
    if min(sand_wet, vertisol_wet) < loess_wet < max(sand_wet, vertisol_wet):
        relative_wetness = "between Sand and Vertisol"
    elif loess_wet > max(sand_wet, vertisol_wet):
        relative_wetness = "higher than both Sand and Vertisol"
    else:
        relative_wetness = "lower than both Sand and Vertisol"
    loess_regime = (
        "Loess has "
        f"{loess['porosity_resolved']['mean']:.3f} resolved porosity and "
        f"{loess_wet:.3f} wet fraction of habitat, {relative_wetness} "
        f"(Sand {sand_wet:.3f}; Vertisol {vertisol_wet:.3f}). "
    )
    if loess["wet_fraction_of_habitat"]["mean"] > 0.5:
        loess_regime += (
            "The prediction is therefore a wet-pore diffusion response with a substantial "
            "drained-film/exchange contribution, rather than a purely drained regime."
        )
    else:
        loess_regime += (
            "The prediction is a drained-film/exchange-dominated response, with the wet-pore "
            "fraction limiting direct diffusion."
        )
    text = f"""# Stage 23 predictions (pre-run)

These predictions were written after Step 1 geometry and before any Stage 23 respiration run.
They use only the fixed scan geometries and chemistry constants. No laboratory respiration
measurements were opened or used.

## Geometry informing the predictions

| Structure soil | Porosity | Matrix | Solid | LT d10 (µm) | LT d50 (µm) | LT d90 (µm) | Largest pore cluster | Spans x/y/z | S_w | Habitat | Wet habitat fraction |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(table)}

Under the inherited matric rule, pore voxels with local thickness below 50 µm are wet and
diffusive; larger pores remain habitat through the drained-film diffusivity and exchange.
{loess_regime}

## Arm A ordering prediction

Loess is predicted to have an intermediate structure-only response between Sand and Vertisol
in peak timing and curve shape: its resolved porosity and pore LT d50 are intermediate, while
its wet habitat fraction and S_w are higher than both other soils. This may make its response
earlier or more active than the porosity ordering alone suggests. The exact peak height depends
on connectivity and habitat fraction. This is a geometry-based prediction, not a claim about
measured respiration. The Loess specimens are expected to show appreciable specimen-to-specimen
spread because LA and LB are downsampled labels while LC is a native 15 µm scan.

## Structure × chemistry cross

| Output | Predicted dominant factor | Pre-run reason |
|---|---|---|
| Peak height (R_peak/floor) | Chemistry | OM and B0 are changed together; the chemistry-specific OM:B0 ratio ranges from about {chem23.CHEMISTRY['chem_S']['OM_density']/chem23.CHEMISTRY['chem_S']['B0_density']:.1f} (Sand) to {chem23.CHEMISTRY['chem_V']['OM_density']/chem23.CHEMISTRY['chem_V']['B0_density']:.1f} (Vertisol), changing substrate availability per initial biomass. |
| Peak time | Structure | Pore LT distributions, connected clusters, wet habitat and matrix exchange distances alter substrate delivery times even when chemistry is held fixed. |
| Cumulative R to day 15 | Chemistry | OM_density spans {min(chem23.CHEMISTRY[c]['OM_density'] for c in ('S','V','L')):.3f}–{max(chem23.CHEMISTRY[c]['OM_density'] for c in ('S','V','L')):.3f}; the largest expected direct effect is on total available substrate, with geometry modulating how much is accessed by day 15. |

These are directional predictions only; the fixed 3×3 cross will quantify main effects and
interaction without tuning or dropping crops.

## Arm D prediction

Filling the narrowest Sand pore voxels with matrix rather than carving fines from inert solid
is expected to reduce resolved porosity and the narrow-pore habitat fraction, while increasing
matrix volume and proximity for exchange. The change should broaden or delay the response and
may change the Sand-versus-Loess/Vertisol ordering; arm D is a model-choice check and will not
be adopted.
"""
    target = Path("docs/stage_results/STAGE23_PREDICTIONS.md")
    if target.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite Stage 23 predictions: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding="utf-8")
    return target


def step1() -> dict:
    if not (OUT / "preflight.json").is_file():
        raise FileNotFoundError("Run Stage 23 Step 0 before Step 1")
    xp = P.get_xp()
    stage22 = _read_json(Path("results/stage22_scan_structure/step1.json"))
    reports: dict[str, dict] = {}
    positions: dict[str, dict] = {}
    crops_by_scan: dict[str, list[str]] = {}
    origin_checks = {}
    arm_d_ids = []
    domain_ids = []

    for scan_id in SCAN_IDS:
        volume = _load_scan(scan_id)
        circular = scan_id == "LB"
        origins = scan3d_loess.grid_origins(volume.shape, CROP_SIZE, circular=circular)
        if scan_id in scan3d.SOIL_BY_ID:
            chosen, rejected = scan3d.draw_habitable_crops(
                volume,
                scan_id,
                VOX_UM[scan_id],
                origins,
                crop_size=CROP_SIZE,
                n_crops=N_CROPS,
                seed=22,
            )
        else:
            chosen, rejected = scan3d_loess.draw_habitable_crops(
                volume,
                scan_id,
                origins,
                crop_size=CROP_SIZE,
                n_crops=N_CROPS,
                seed=22,
            )
        if scan_id in stage22["positions"]:
            registered = [
                tuple(int(value) for value in point)
                for point in stage22["positions"][scan_id]["selected"]
            ]
            match = chosen == registered
            origin_checks[scan_id] = {
                "reproduced_stage22_origins": match,
                "selected": chosen,
                "stage22_selected": registered,
            }
            if not match:
                raise RuntimeError(
                    f"Stage 22 crop origins changed for {scan_id}: {chosen} != {registered}"
                )
        if scan_id == "LB":
            outside = [
                origin
                for origin in chosen
                if not scan3d_loess.crop_inside_disk(origin, CROP_SIZE, volume.shape)
            ]
            if outside:
                raise RuntimeError(f"LB crop reaches outside the core: {outside}")

        positions[scan_id] = {
            "grid_origins_100": origins,
            "selected": chosen,
            "rejected_zero_habitat": rejected,
            "circular_constraint": circular,
        }
        crops_by_scan[scan_id] = []
        for crop_index, origin in enumerate(chosen):
            structure_id = f"{scan_id}_c{crop_index}"
            crop = (
                scan3d.crop_at(volume, origin, CROP_SIZE)
                if scan_id in scan3d.SOIL_BY_ID
                else scan3d_loess.crop_at(volume, origin, CROP_SIZE)
            )
            structure, threshold = _make_structure(crop, scan_id, xp)
            report = _analyse_structure(
                structure_id,
                structure,
                xp,
                {
                    "origin": origin,
                    "solid_LT_threshold_um": threshold,
                    "fines_rule": "solid_lt" if scan_id in ("SA", "SB") else "not_applicable",
                },
            )
            reports[structure_id] = report
            crops_by_scan[scan_id].append(structure_id)
            if scan_id in ("SA", "SB"):
                pore_phase, n_fines = scan3d_loess.pore_fill_fines_phase(
                    crop, VOX_UM[scan_id], seed=22, xp=xp
                )
                alternative = P.Structure()
                alternative.kind = f"{structure_id}_porefill"
                alternative.n = CROP_SIZE
                alternative.vox_um = VOX_UM[scan_id]
                alternative.seed = 22
                alternative.phase = xp.asarray(pore_phase)
                alternative.info = {
                    "scan_id": scan_id,
                    "fines_rule": "stage21_pore_fill",
                    "n_fines_voxels": n_fines,
                }
                arm_d_id = f"{structure_id}_porefill"
                arm_d_ids.append(arm_d_id)
                reports[arm_d_id] = _analyse_structure(
                    arm_d_id,
                    alternative,
                    xp,
                    {
                        "origin": origin,
                        "fines_rule": "stage21_pore_fill",
                        "n_fines_voxels": n_fines,
                    },
                )

        if scan_id in scan3d_loess.LOESS_IDS:
            domain_origin, domain_size = scan3d.centered_cube_origin(volume.shape, 256)
            if domain_size != 256:
                raise ValueError(f"{scan_id} cannot supply the required 256-cube domain check")
            if scan_id == "LB" and not scan3d_loess.crop_inside_disk(
                domain_origin, domain_size, volume.shape
            ):
                raise RuntimeError(
                    f"LB 256-cube domain reaches outside the core: {domain_origin}"
                )
            crop = scan3d_loess.crop_at(volume, domain_origin, domain_size)
            domain = scan3d_loess.crop_to_structure(crop, scan_id, xp)
            domain_id = f"{scan_id}_domain256"
            reports[domain_id] = _analyse_structure(
                domain_id,
                domain,
                xp,
                {
                    "origin": domain_origin,
                    "requested_size": 256,
                    "actual_size": domain_size,
                    "fines_rule": "not_applicable",
                    "lb_inside_crop_safe_disk": (
                        scan3d_loess.crop_inside_disk(domain_origin, domain_size, volume.shape)
                        if scan_id == "LB"
                        else None
                    ),
                },
            )
            domain_ids.append(domain_id)
            positions[scan_id]["domain_origin"] = domain_origin
            positions[scan_id]["domain_actual_size"] = domain_size
        del volume

    aggregate_ids = _geometry_aggregate_sets(crops_by_scan, arm_d_ids, domain_ids)
    aggregates = {
        key: _aggregate_geometry(reports, structure_ids)
        for key, structure_ids in aggregate_ids.items()
    }
    result = {
        "stage": 23,
        "step": 1,
        "input_gates": {
            key: value["gate"]
            for key, value in _read_json(OUT / "preflight.json")["input_gates"].items()
        },
        "positions": positions,
        "origin_checks_against_stage22": origin_checks,
        "reports": reports,
        "aggregates": aggregates,
        "arm_d_structure_ids": arm_d_ids,
        "domain_structure_ids": domain_ids,
        "crop_count": sum(len(ids) for ids in crops_by_scan.values()),
        "lb_crop_constraint": {
            "grid": "centered non-overlapping grid from run_tortuosity_subvolumes.py",
            "corner_radius_fraction": 0.90,
            "selected_origins_inside_crop_safe_disk": True,
            "domain_256_inside_crop_safe_disk": True,
        },
        "respiration_runs_started": False,
    }
    slices_path = _write_step1_slices()
    prediction_path = _write_predictions(result)
    result["slice_figure"] = str(slices_path)
    result["predictions"] = str(prediction_path)
    _write_json(OUT / "step1.json", result)
    print(json.dumps({
        "crop_count": result["crop_count"],
        "origin_checks_against_stage22": origin_checks,
        "aggregates": aggregates,
        "arm_d_structures": arm_d_ids,
        "domain_structures": domain_ids,
        "slice_figure": str(slices_path),
        "predictions": str(prediction_path),
    }, indent=2, default=_json_default))
    return result


def refresh_step1_outputs() -> dict:
    """Finalize geometry summaries/predictions from this run's Step 1 outputs."""
    result = _read_json(OUT / "step1.json")
    if result.get("respiration_runs_started"):
        raise RuntimeError("Cannot refresh Step 1 predictions after respiration runs began")
    crops_by_scan = {
        scan_id: [f"{scan_id}_c{index}" for index in range(N_CROPS)]
        for scan_id in SCAN_IDS
    }
    arm_d_ids = result["arm_d_structure_ids"]
    domain_ids = result["domain_structure_ids"]
    aggregate_ids = _geometry_aggregate_sets(crops_by_scan, arm_d_ids, domain_ids)
    result["aggregates"] = {
        key: _aggregate_geometry(result["reports"], ids)
        for key, ids in aggregate_ids.items()
    }
    _write_predictions(result, overwrite=True)
    target = OUT / "step1.json"
    target.write_text(json.dumps(result, indent=2, default=_json_default) + "\n", encoding="utf-8")
    return result


def main(argv: list[str] | None = None) -> None:
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 1:
        raise SystemExit("Usage: python stage23_run.py step0|step1|step2|step3|step4|step5|step6|step7")
    steps = {
        "step0": step0,
        "step1": step1,
        "refresh_step1": refresh_step1_outputs,
    }
    if args[0] not in steps:
        raise SystemExit(f"Step {args[0]!r} is not implemented yet")
    steps[args[0]]()


if __name__ == "__main__":
    main()
