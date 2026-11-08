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
    subprocess.run(test_command, check=True)
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


def main(argv: list[str] | None = None) -> None:
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 1:
        raise SystemExit("Usage: python stage23_run.py step0|step1|step2|step3|step4|step5|step6|step7")
    steps = {
        "step0": step0,
    }
    if args[0] not in steps:
        raise SystemExit(f"Step {args[0]!r} is not implemented yet")
    steps[args[0]]()


if __name__ == "__main__":
    main()
