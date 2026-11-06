"""Scan-backed structures for Stage 22 (Sand and Vertisol only)."""
from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np

import packing3d as P


SOIL_BY_ID = {"VA": "vertisol", "VB": "vertisol", "SA": "sand", "SB": "sand"}


def load_scan(scan_id: str, root: str | Path) -> np.ndarray:
    """Load only the segmented scan volume identified by the Stage 22 prompt."""
    root = Path(root)
    if scan_id == "VA":
        import nibabel as nib

        path = root / "nnUNet_resources" / "bnei_reem_specA_unified" / "inference_concatenated" / "bnei_reem_specA_unified.nii.gz"
        return np.asanyarray(nib.load(str(path)).dataobj)
    if scan_id == "VB":
        import nibabel as nib

        path = root / "nnUNet_resources" / "bnei_reem_specB_unified" / "inference_concatenated" / "bnei_reem_specB_unified.nii.gz"
        return np.asanyarray(nib.load(str(path)).dataobj)
    if scan_id == "SA":
        return np.load(root / "10.5" / "rehovot_samp_2.npy", mmap_mode="r")
    if scan_id == "SB":
        return np.load(root / "10.5" / "Rehovot_samp3_highkV_Cu0.11_15um.npy", mmap_mode="r")
    raise ValueError(f"Unknown scan ID: {scan_id}")


def scan_histogram(scan_id: str, volume: np.ndarray) -> dict:
    values, counts = np.unique(volume, return_counts=True)
    size = int(volume.size)
    fractions = {int(v): float(n / size) for v, n in zip(values, counts)}
    if SOIL_BY_ID[scan_id] == "vertisol":
        invalid = {int(v): int(n) for v, n in zip(values, counts) if int(v) not in {0, 1, 2, 3, 5}}
        return {
            "dtype": str(volume.dtype),
            "shape": list(volume.shape),
            "fractions": fractions,
            "invalid_counts": invalid,
            "pore_fraction": fractions.get(5, 0.0),
            "pom_fraction": fractions.get(2, 0.0) + fractions.get(3, 0.0),
        }
    return {
        "dtype": str(volume.dtype),
        "shape": list(volume.shape),
        "unique_values": [v.item() for v in values],
        "pore_fraction": float(np.count_nonzero(volume == 1) / size),
    }


def scan_gate(scan_id: str, report: dict) -> bool:
    expected_pore = {"VA": 15.17, "VB": 21.64, "SA": 30.50, "SB": 40.85}[scan_id] / 100
    pore_ok = abs(report["pore_fraction"] - expected_pore) <= 0.01
    if SOIL_BY_ID[scan_id] == "vertisol":
        expected_pom = {"VA": 0.24, "VB": 0.82}[scan_id] / 100
        return not report["invalid_counts"] and pore_ok and abs(report["pom_fraction"] - expected_pom) <= 0.005
    return set(report["unique_values"]).issubset({0, 1}) and pore_ok


def select_fines(
    solid: np.ndarray,
    local_thickness_um: np.ndarray,
    pore_distance_vox: np.ndarray,
    fraction: float,
    seed: int = 22,
) -> tuple[np.ndarray, float]:
    """Select the thinnest solid voxels, then nearest to pores, then by seeded tie-break."""
    if solid.shape != local_thickness_um.shape or solid.shape != pore_distance_vox.shape:
        raise ValueError("Solid, local-thickness, and pore-distance arrays must have matching shapes")
    if not 0.0 <= fraction <= 1.0:
        raise ValueError("Matrix fraction must be between zero and one")
    solid_flat = np.flatnonzero(solid)
    n_target = int(round(fraction * solid.size))
    if n_target > solid_flat.size:
        raise ValueError("Requested matrix fraction exceeds the available solid fraction")
    selected = np.zeros(solid.size, dtype=bool)
    if n_target:
        rng = np.random.default_rng(seed)
        tie = rng.random(solid_flat.size)
        order = np.lexsort(
            (tie, pore_distance_vox.ravel()[solid_flat], local_thickness_um.ravel()[solid_flat])
        )
        selected[solid_flat[order[:n_target]]] = True
        threshold = float(local_thickness_um.ravel()[solid_flat[order[n_target - 1]]])
    else:
        threshold = 0.0
    return selected.reshape(solid.shape), threshold


def crop_to_structure(crop: np.ndarray, scan_id: str, vox_um: float, seed: int = 22):
    """Map one scan crop to a packing3d.Structure with Stage 22 phase rules."""
    phase = np.empty(crop.shape, dtype=np.uint8)
    if SOIL_BY_ID[scan_id] == "vertisol":
        phase[crop == 5] = P.PORE
        phase[np.isin(crop, (0, 2, 3))] = P.MATRIX
        phase[crop == 1] = P.SOLID
        if not np.all(np.isin(crop, (0, 1, 2, 3, 5))):
            raise ValueError(f"Unexpected Vertisol label in {scan_id} crop")
        threshold = None
    else:
        pore = np.asarray(crop == 1)
        solid = ~pore
        thickness = P.local_thickness(solid, vox_um, np)
        pore_distance = P.edt(solid, np)
        matrix_fraction = P.fines_volume_fraction()
        matrix, threshold = select_fines(solid, thickness, pore_distance, matrix_fraction, seed)
        phase.fill(P.SOLID)
        phase[pore] = P.PORE
        phase[matrix] = P.MATRIX
    structure = P.Structure()
    structure.kind = scan_id
    structure.n = int(crop.shape[0])
    structure.vox_um = float(vox_um)
    structure.seed = int(seed)
    structure.phase = phase
    structure.info = {"scan_id": scan_id, "solid_lt_threshold_um": threshold}
    return structure, threshold


def grid_origins(shape: Iterable[int], crop_size: int) -> list[tuple[int, int, int]]:
    dims = tuple(int(v) for v in shape)
    if len(dims) != 3 or crop_size <= 0 or any(v < crop_size for v in dims):
        raise ValueError(f"Crop size {crop_size} does not fit volume shape {dims}")
    return [
        (x, y, z)
        for x in range(0, dims[0] - crop_size + 1, crop_size)
        for y in range(0, dims[1] - crop_size + 1, crop_size)
        for z in range(0, dims[2] - crop_size + 1, crop_size)
    ]


def crop_at(volume: np.ndarray, origin: tuple[int, int, int], size: int) -> np.ndarray:
    x, y, z = origin
    return np.asarray(volume[x : x + size, y : y + size, z : z + size])


def draw_habitable_crops(
    volume: np.ndarray,
    scan_id: str,
    vox_um: float,
    origins: list[tuple[int, int, int]],
    crop_size: int = 100,
    n_crops: int = 4,
    seed: int = 22,
) -> tuple[list[tuple[int, int, int]], list[tuple[int, int, int]]]:
    """Draw unique origins and reject only crops without >=10 um pore habitat."""
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(origins))
    chosen: list[tuple[int, int, int]] = []
    rejected: list[tuple[int, int, int]] = []
    for index in order:
        origin = origins[int(index)]
        crop = crop_at(volume, origin, crop_size)
        if crop.shape != (crop_size, crop_size, crop_size):
            raise ValueError(f"Crop selector requires {crop_size}-cube origins")
        structure, _ = crop_to_structure(crop, scan_id, vox_um)
        lt = P.local_thickness(structure.phase == P.PORE, vox_um, np)
        if np.any((structure.phase == P.PORE) & (lt >= 10.0)):
            chosen.append(origin)
            if len(chosen) == n_crops:
                return chosen, rejected
        else:
            rejected.append(origin)
    raise ValueError(f"Only {len(chosen)} habitable crops found for {scan_id}; rejected {rejected}")


def centered_cube_origin(shape: Iterable[int], requested_size: int) -> tuple[tuple[int, int, int], int]:
    dims = tuple(int(v) for v in shape)
    size = min(requested_size, *dims)
    if size <= 0:
        raise ValueError(f"Invalid volume shape: {dims}")
    origin = tuple((v - size) // 2 for v in dims)
    return origin, size
