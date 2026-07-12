"""
Emergent-distinctness check + shape auto-describer (soil_respiration_prompt_v2.md
"New success criterion").

Distinctness is measured pairwise between the three (self-normalized) R(t)
curves; shapes are *described*, never targeted. The search loop only ever
optimizes distinctness.
"""
from __future__ import annotations

import numpy as np


def normalize(R):
    R = np.asarray(R, dtype=float)
    peak = R.max()
    if peak <= 0:
        return R.copy()
    return R / peak


def pearson_corr(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a - a.mean()
    b = b - b.mean()
    denom = np.sqrt((a ** 2).sum() * (b ** 2).sum())
    if denom <= 0:
        return 0.0
    return float((a * b).sum() / denom)


def pairwise_distance(Ri, Rj, l2_weight=0.5):
    """dist = 1 - pearson_corr(norm_Ri, norm_Rj), plus a normalized L2 term
    (mean absolute difference of the normalized curves) so two curves with
    matching correlation but different amplitude/offset still register as
    distinct."""
    ni, nj = normalize(Ri), normalize(Rj)
    corr_term = 1.0 - pearson_corr(ni, nj)
    l2_term = float(np.mean(np.abs(ni - nj)))
    return corr_term + l2_weight * l2_term


def peak_count(R, prominence_frac=0.1):
    """count local maxima in R whose prominence exceeds prominence_frac*peak
    (simple, dependency-free peak counter)."""
    R = np.asarray(R, dtype=float)
    peak = R.max()
    if peak <= 0 or len(R) < 3:
        return 0
    thresh = prominence_frac * peak
    d = np.diff(R)
    sign = np.sign(d)
    count = 0
    for i in range(1, len(sign)):
        if sign[i - 1] > 0 and sign[i] <= 0:
            # local max candidate at index i
            left = R[max(0, i - 5):i + 1].min()
            right = R[i:min(len(R), i + 6)].max()
            if R[i] - min(left, right) >= thresh or R[i] >= peak - thresh:
                count += 1
    return max(count, 1) if peak > 0 else 0


def describe_shape(R, flat_tol=0.10):
    """rising / falling / flat / single-peak / multi-peak, via trend sign +
    peak count. Descriptive only -- never used to steer the search."""
    R = np.asarray(R, dtype=float)
    n = len(R)
    third = n // 3
    e = R[:third].mean()
    l = R[n - third:].mean()
    peak = R.max()
    npeaks = peak_count(R)

    if npeaks >= 2:
        shape = "multi-peak"
    elif abs(l - e) < flat_tol * peak:
        shape = "flat"
    elif l > e:
        shape = "rising" if npeaks <= 1 else "single-peak"
    else:
        shape = "falling"

    return {"e": float(e), "l": float(l), "peak": float(peak),
            "n_peaks": int(npeaks), "shape": shape}


def is_nontrivial(out):
    """biomass actually grows, total respiration > 0."""
    grids = out["grids"]
    B0_total = float(grids["B0"].sum())
    B_final_total = float(out["B_final"].sum())
    R_total = float(out["R_t"].sum())
    return (B_final_total > B0_total) and (R_total > 0.0)


def emergent_distinctness(out_A, out_B, out_C, threshold=0.3, l2_weight=0.5):
    """Success = all three pairwise distances exceed `threshold` AND every
    soil is non-trivial. Returns (success: bool, detail: dict)."""
    curves = {"A": out_A["R_t"], "B": out_B["R_t"], "C": out_C["R_t"]}
    outs = {"A": out_A, "B": out_B, "C": out_C}

    distances = {}
    for i, j in [("A", "B"), ("A", "C"), ("B", "C")]:
        distances[f"{i}{j}"] = pairwise_distance(curves[i], curves[j], l2_weight)

    nontrivial = {s: is_nontrivial(outs[s]) for s in "ABC"}
    shapes = {s: describe_shape(curves[s]) for s in "ABC"}

    dist_ok = all(d > threshold for d in distances.values())
    trivial_ok = all(nontrivial.values())
    success = dist_ok and trivial_ok

    return success, {
        "distances": distances, "threshold": threshold,
        "nontrivial": nontrivial, "shapes": shapes,
        "min_distance": min(distances.values()),
    }
