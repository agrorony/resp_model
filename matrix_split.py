"""Fixed Stage 24 matrix/inert-solid partition rules."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

import packing3d as P
import scan3d


@dataclass(frozen=True)
class Variant:
    name: str
    sand_fraction: float | None = None
    fine_soil_inert_share: float = 0.0
    affects: str = "all"


BASE_SAND_FRACTION = P.fines_volume_fraction()
VARIANTS = (
    Variant("base", BASE_SAND_FRACTION, 0.0, "all"),
    Variant("S04", 0.040, 0.0, "sand"),
    Variant("S16", 0.160, 0.0, "sand"),
    Variant("S30", 0.300, 0.0, "sand"),
    Variant("VL25", BASE_SAND_FRACTION, 0.25, "fine"),
    Variant("VL50", BASE_SAND_FRACTION, 0.50, "fine"),
    Variant("CONV", 0.300, 0.50, "all"),
)


def variant_map() -> dict[str, Variant]:
    return {variant.name: variant for variant in VARIANTS}


def variant_changes(variant: Variant, scan_id: str) -> bool:
    if variant.name == "base":
        return False
    soil = scan3d.SOIL_BY_ID.get(scan_id)
    if soil == "sand":
        return variant.affects in ("sand", "all")
    return variant.affects in ("fine", "all")


def split_phase(
    crop: np.ndarray,
    scan_id: str,
    vox_um: float,
    variant: Variant,
) -> tuple[np.ndarray, dict]:
    """Return phase codes and fixed split diagnostics for one segmented crop."""
    if scan_id in scan3d.SOIL_BY_ID:
        if scan3d.SOIL_BY_ID[scan_id] == "sand":
            pore = crop == 1
            solid = ~pore
            if variant.sand_fraction is None:
                raise ValueError(f"{variant.name} has no Sand matrix target")
            thickness = P.to_np(P.local_thickness(solid, vox_um, np))
            pore_distance = P.to_np(P.edt(solid, np))
            solid_flat = np.flatnonzero(solid)
            n_target_requested = int(round(variant.sand_fraction * crop.size))
            n_target = min(n_target_requested, int(solid_flat.size))
            matrix = np.zeros(crop.shape, dtype=bool)
            if n_target:
                tie = np.random.default_rng(22).random(solid_flat.size)
                order = np.lexsort(
                    (
                        tie,
                        pore_distance.ravel()[solid_flat],
                        thickness.ravel()[solid_flat],
                    )
                )
                matrix.ravel()[solid_flat[order[:n_target]]] = True
            phase = np.full(crop.shape, P.SOLID, dtype=np.uint8)
            phase[pore] = P.PORE
            phase[matrix] = P.MATRIX
            return phase, {
                "matrix_target_fraction_of_crop": variant.sand_fraction,
                "matrix_target_voxels_requested": n_target_requested,
                "matrix_target_voxels_selected": n_target,
                "solid_voxels_available": int(solid_flat.size),
                "target_clamped_to_available_solid": n_target < n_target_requested,
            }

    if not np.all(np.isin(crop, (0, 1, 2, 3, 5))):
        raise ValueError(f"Unexpected fine-soil labels in {scan_id} crop")
    pore = crop == 5
    original_matrix = np.isin(crop, (0, 2, 3))
    phase = np.full(crop.shape, P.SOLID, dtype=np.uint8)
    phase[pore] = P.PORE
    phase[original_matrix] = P.MATRIX
    n_matrix = int(original_matrix.sum())
    n_convert = int(round(variant.fine_soil_inert_share * n_matrix))
    if n_convert:
        non_pore = ~pore
        local_thickness = P.to_np(P.local_thickness(non_pore, vox_um, np))
        pore_distance = P.to_np(P.edt(non_pore, np))
        candidates = np.flatnonzero(original_matrix)
        tie = np.random.default_rng(24).random(candidates.size)
        order = np.lexsort(
            (
                tie,
                -pore_distance.ravel()[candidates],
                -local_thickness.ravel()[candidates],
            )
        )
        phase.ravel()[candidates[order[:n_convert]]] = P.SOLID
    return phase, {
        "original_matrix_voxels": n_matrix,
        "fine_soil_inert_share_target": variant.fine_soil_inert_share,
        "matrix_voxels_converted_to_inert": n_convert,
    }


def geometry(phase: np.ndarray, vox_um: float, om_densities: dict[str, float]) -> dict:
    """Measure phase fractions, interface, distance-to-pore, and OM loading."""
    pore = phase == P.PORE
    matrix = phase == P.MATRIX
    solid = phase == P.SOLID
    faces = 0
    for axis in range(3):
        lo = [slice(None)] * 3
        hi = [slice(None)] * 3
        lo[axis] = slice(None, -1)
        hi[axis] = slice(1, None)
        faces += int(np.count_nonzero(matrix[tuple(lo)] & pore[tuple(hi)]))
        faces += int(np.count_nonzero(pore[tuple(lo)] & matrix[tuple(hi)]))
    if not matrix.any():
        raise ValueError("A Stage 24 crop contains no matrix voxels")
    distance_um = P.to_np(P.edt(~pore, np)) * vox_um
    return {
        "matrix_fraction": float(matrix.mean()),
        "solid_fraction": float(solid.mean()),
        "pore_fraction": float(pore.mean()),
        "matrix_pore_interface_faces_per_crop_voxel": faces / float(phase.size),
        "mean_matrix_distance_to_pore_um": float(distance_um[matrix].mean()),
        "matrix_voxels": int(matrix.sum()),
        "matrix_interface_faces": faces,
        "om_per_matrix_voxel": {
            arm: float(density * phase.size / matrix.sum())
            for arm, density in om_densities.items()
        },
    }
