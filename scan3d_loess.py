"""Loess scan loading, circular crop constraints, and Stage 23 phase mapping."""
from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np

import packing3d as P
import scan3d


LOESS_IDS = ("LA", "LB", "LC")
VOX_UM = {"LA": 15.007595, "LB": 14.998099, "LC": 15.000149}
EXPECTED_PORE = {"LA": 0.2657, "LB": 0.2034, "LC": 0.2618}
EXPECTED_POM = {"LA": 0.020, "LB": 0.079, "LC": 0.001}
VALID_LABELS = (0, 1, 2, 3, 5)


def load_scan(scan_id: str, root: str | Path) -> np.ndarray:
    paths = {
        "LA": (
            "nnUNet4SoilXrayCT",
            "05_evaluation",
            "psd",
            "pom_analysis_20260829_roi_expansion",
            "pipeline_work",
            "mishmar_native_5p85um",
            "mishmar_native_5p85um_label_downsample.nii.gz",
        ),
        "LB": (
            "nnUNet4SoilXrayCT",
            "05_evaluation",
            "psd",
            "pom_analysis_20260829_roi_expansion",
            "pipeline_work",
            "mishmar_second_8p8um",
            "mishmar_second_8p8um_label_downsample.nii.gz",
        ),
        "LC": (
            "nnUNet_resources",
            "mishmar_hanegev_Cu011_samp_2_Rec_nlm",
            "inference_output_concat_loess_i2",
            "mishmar_hanegev_Cu011_samp_2_Rec_nlm.nii.gz",
        ),
    }
    if scan_id not in paths:
        raise ValueError(f"Unknown Loess scan ID: {scan_id}")
    import nibabel as nib

    path = Path(root).joinpath(*paths[scan_id])
    return np.asanyarray(nib.load(str(path)).dataobj)


def inscribed_disk(shape: Iterable[int]) -> tuple[float, float, float]:
    dims = tuple(int(value) for value in shape)
    if len(dims) != 3:
        raise ValueError(f"Expected a 3-D volume shape, got {dims}")
    cx, cy = dims[0] / 2.0, dims[1] / 2.0
    return cx, cy, min(cx, cy)


def disk_label_counts(volume: np.ndarray) -> dict[int, int]:
    cx, cy, radius = inscribed_disk(volume.shape)
    xx, yy = np.ogrid[: volume.shape[0], : volume.shape[1]]
    disk = (xx - cx) ** 2 + (yy - cy) ** 2 <= radius**2
    counts: dict[int, int] = {}
    for z in range(volume.shape[2]):
        values, numbers = np.unique(volume[:, :, z][disk], return_counts=True)
        for value, number in zip(values, numbers):
            key = int(value)
            counts[key] = counts.get(key, 0) + int(number)
    return counts


def _fractions(counts: dict[int, int]) -> dict[int, float]:
    total = sum(counts.values())
    return {value: number / total for value, number in counts.items()}


def scan_histogram(scan_id: str, volume: np.ndarray) -> dict:
    values, numbers = np.unique(volume, return_counts=True)
    counts = {int(value): int(number) for value, number in zip(values, numbers)}
    invalid = {value: number for value, number in counts.items() if value not in VALID_LABELS}
    whole_fractions = _fractions(counts)
    pom_whole = whole_fractions.get(2, 0.0) + whole_fractions.get(3, 0.0)
    result = {
        "dtype": str(volume.dtype),
        "shape": list(volume.shape),
        "whole_volume_counts": counts,
        "whole_volume_fractions": whole_fractions,
        "invalid_counts": invalid,
        "whole_volume_pore_fraction": whole_fractions.get(5, 0.0),
        "whole_volume_pom_fraction": pom_whole,
    }
    if scan_id == "LB":
        disk_counts = disk_label_counts(volume)
        disk_fractions = _fractions(disk_counts)
        result["inscribed_disk"] = {
            "radius_voxels": inscribed_disk(volume.shape)[2],
            "voxel_count": sum(disk_counts.values()),
            "counts": disk_counts,
            "fractions": disk_fractions,
            "pore_fraction": disk_fractions.get(5, 0.0),
            "pom_fraction": disk_fractions.get(2, 0.0) + disk_fractions.get(3, 0.0),
        }
    expected_pore = EXPECTED_PORE[scan_id]
    result["gate"] = {
        "valid_labels": not invalid,
        "whole_volume_pore_fraction_within_1pp": (
            abs(result["whole_volume_pore_fraction"] - expected_pore) <= 0.01
        ),
    }
    result["gate"]["pass"] = all(result["gate"].values())
    return result


def crop_inside_disk(
    origin: tuple[int, int, int], crop_size: int, shape: Iterable[int], safety: float = 0.90
) -> bool:
    dims = tuple(int(value) for value in shape)
    cx, cy, radius = inscribed_disk(dims)
    radius *= safety
    x, y, _ = origin
    corners = (
        (x + dx - cx, y + dy - cy)
        for dx in (0, crop_size)
        for dy in (0, crop_size)
    )
    return max(np.hypot(*corner) for corner in corners) <= radius


def grid_origins(
    shape: Iterable[int], crop_size: int, circular: bool = False
) -> list[tuple[int, int, int]]:
    dims = tuple(int(value) for value in shape)
    if not circular:
        return scan3d.grid_origins(dims, crop_size)
    if len(dims) != 3 or crop_size <= 0 or any(value < crop_size for value in dims):
        raise ValueError(f"Crop size {crop_size} does not fit volume shape {dims}")
    axes = []
    for length in dims:
        count = length // crop_size
        offset = (length - count * crop_size) // 2
        axes.append([offset + index * crop_size for index in range(count)])
    origins = []
    for x in axes[0]:
        for y in axes[1]:
            for z in axes[2]:
                origin = (x, y, z)
                if crop_inside_disk(origin, crop_size, dims):
                    origins.append(origin)
    return origins


def crop_at(volume: np.ndarray, origin: tuple[int, int, int], size: int) -> np.ndarray:
    x, y, z = origin
    return np.asarray(volume[x : x + size, y : y + size, z : z + size])


def crop_to_structure(crop: np.ndarray, scan_id: str, xp=np) -> P.Structure:
    if scan_id not in LOESS_IDS:
        raise ValueError(f"Unknown Loess scan ID: {scan_id}")
    if not np.all(np.isin(crop, VALID_LABELS)):
        invalid = np.unique(crop[~np.isin(crop, VALID_LABELS)]).tolist()
        raise ValueError(f"Unexpected label(s) in {scan_id} crop: {invalid}")
    phase = np.full(crop.shape, P.SOLID, dtype=np.uint8)
    phase[crop == 5] = P.PORE
    phase[np.isin(crop, (0, 2, 3))] = P.MATRIX
    structure = P.Structure()
    structure.kind = scan_id
    structure.n = int(crop.shape[0])
    structure.vox_um = VOX_UM[scan_id]
    structure.seed = 22
    structure.phase = xp.asarray(phase)
    structure.info = {"scan_id": scan_id}
    return structure


def draw_habitable_crops(
    volume: np.ndarray,
    scan_id: str,
    origins: list[tuple[int, int, int]],
    crop_size: int = 100,
    n_crops: int = 4,
    seed: int = 22,
) -> tuple[list[tuple[int, int, int]], list[tuple[int, int, int]]]:
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(origins))
    chosen: list[tuple[int, int, int]] = []
    rejected: list[tuple[int, int, int]] = []
    for index in order:
        origin = origins[int(index)]
        crop = crop_at(volume, origin, crop_size)
        if crop.shape != (crop_size, crop_size, crop_size):
            raise ValueError(f"Crop selector requires {crop_size}-cube origins")
        structure = crop_to_structure(crop, scan_id)
        lt = P.local_thickness(structure.phase == P.PORE, VOX_UM[scan_id], np)
        if np.any((structure.phase == P.PORE) & (lt >= 10.0)):
            chosen.append(origin)
            if len(chosen) == n_crops:
                return chosen, rejected
        else:
            rejected.append(origin)
    raise ValueError(
        f"Only {len(chosen)} habitable crops found for {scan_id}; rejected {rejected}"
    )


def pore_fill_fines_phase(
    crop: np.ndarray, vox_um: float, seed: int = 22, xp=np
) -> tuple[np.ndarray, int]:
    """Build the Stage 21 Sand alternative by filling pore voxels with fines."""
    pore_np = np.asarray(crop == 1)
    n_target = int(round(P.fines_volume_fraction() * pore_np.size))
    if n_target > int(pore_np.sum()):
        raise ValueError("Sand fines target exceeds the available pore volume")
    pore = xp.asarray(pore_np)
    thickness = P.local_thickness(pore, vox_um, xp)
    fines = P._fill_fines(pore, thickness, n_target, vox_um, seed, xp)
    phase = np.full(crop.shape, P.SOLID, dtype=np.uint8)
    phase[pore_np] = P.PORE
    phase[P.to_np(fines)] = P.MATRIX
    if np.any(phase[np.asarray(crop) == 0] != P.SOLID):
        raise AssertionError("Pore-fill fines changed a solid voxel")
    return phase, n_target
