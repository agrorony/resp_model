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


def classify_outcome(out, dead_frac=0.02):
    """Stage 3 outcome axis (soil_respiration_prompt_stage3.md SS8): every
    soil's run is `thriving` / `declining` / `dead`, via final-vs-initial
    biomass and the respiration trend -- never a disqualifier for success, a
    structurally-caused die-off is a legitimate, strongly distinct pattern
    in its own right.

    - `dead`: final biomass has collapsed to <= `dead_frac` of its initial
      total (or respiration has effectively stopped in the last few steps).
    - `thriving`: biomass held or grew (final >= initial) OR respiration is
      still trending up/flat late in the window.
    - `declining`: shrank but did not collapse.
    """
    grids = out["grids"]
    B0_total = float(grids["B0"].sum())
    B_final_total = float(out["B_final"].sum())
    R = np.asarray(out["R_t"], dtype=float)
    ratio = B_final_total / B0_total if B0_total > 0 else 0.0

    tail = min(5, len(R))
    if ratio <= dead_frac or R[-tail:].sum() <= 0.0:
        return "dead"

    shape = describe_shape(R)
    if ratio >= 1.0 or shape["l"] >= shape["e"]:
        return "thriving"
    return "declining"


def label_connected(mask):
    """4-connectivity connected-component labeling of a boolean mask
    (plain BFS flood fill -- grids are <=50x50 so this is fast without
    scipy.ndimage.label; numpy-only per stage-2 invariant 8).

    Returns (labels, n_components): labels[i,j] == 0 for background
    (mask False), 1..n_components over the mask's connected components."""
    mask = np.asarray(mask, dtype=bool)
    labels = np.zeros(mask.shape, dtype=int)
    ny, nx = mask.shape
    n = 0
    for i0 in range(ny):
        for j0 in range(nx):
            if mask[i0, j0] and labels[i0, j0] == 0:
                n += 1
                stack = [(i0, j0)]
                labels[i0, j0] = n
                while stack:
                    i, j = stack.pop()
                    for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                        ii, jj = i + di, j + dj
                        if 0 <= ii < ny and 0 <= jj < nx and mask[ii, jj] and labels[ii, jj] == 0:
                            labels[ii, jj] = n
                            stack.append((ii, jj))
    return labels, n


def largest_cluster_fraction(water_mask):
    """Largest connected water cluster (4-connectivity) as a fraction of all
    water-filled cells. 0.0 if the soil is bone dry (no water-filled cells)."""
    water_mask = np.asarray(water_mask, dtype=bool)
    total = int(water_mask.sum())
    if total == 0:
        return 0.0
    labels, n = label_connected(water_mask)
    sizes = [int((labels == k).sum()) for k in range(1, n + 1)]
    return max(sizes) / total if sizes else 0.0


def isolation_ratio(water_mask):
    """Stage 3 SS9 hypothesis-test diagnostic: number of isolated water
    clusters (4-connectivity connected components of the water-filled mask)
    divided by the total number of saturated cells. Returns
    (isolation_ratio, n_clusters). 0.0 (0 clusters) if bone dry -- a dry
    soil has no water network to fragment."""
    water_mask = np.asarray(water_mask, dtype=bool)
    total = int(water_mask.sum())
    if total == 0:
        return 0.0, 0
    _, n = label_connected(water_mask)
    return n / total, n


HIGH_PERCENTILE = 75.0  # "high-OM" / "high-K" cells = top quartile of that field


def om_biomass_connectivity(water_mask, OM, B0, high_percentile=HIGH_PERCENTILE):
    """Are the high-OM (finest-pore) cells in the same connected water
    component as where biomass actually sits? Returns the fraction of
    INITIAL biomass (B0-weighted, over the WHOLE grid) that sits in a water
    component containing at least one high-OM cell. "High-OM" = top
    `high_percentile` of the OM field's values (OM is heavily skewed toward
    fine pores by construction, so its top quartile concentrates most of the
    mass)."""
    water_mask = np.asarray(water_mask, dtype=bool)
    B0_total = float(B0.sum())
    if B0_total <= 0 or not water_mask.any():
        return 0.0

    om_thresh = np.percentile(OM, high_percentile)
    high_om = water_mask & (OM >= om_thresh)
    if not high_om.any():
        return 0.0

    labels, _ = label_connected(water_mask)
    om_components = set(np.unique(labels[high_om])) - {0}
    reachable = np.isin(labels, list(om_components))
    return float(B0[reachable].sum() / B0_total)


def accessible_om_fraction(water_mask, OM, K, high_percentile=HIGH_PERCENTILE):
    """Total OM reachable through the water network from any biomass-
    favorable cell, over total OM (summed over the WHOLE grid, water and
    dry, since OM_total is a fixed constant regardless of saturation).
    "Biomass-favorable" = a water-filled cell in the top `high_percentile`
    of K (the near-r_opt habitat where b0_field concentrates the initial
    population). A dry cell has D=0 on every face (harmonic mean with a
    zero side is zero), so it is isolated from the water network -- its OM
    is only ever locally accessible, never "reachable" from elsewhere."""
    water_mask = np.asarray(water_mask, dtype=bool)
    OM_total = float(OM.sum())
    if OM_total <= 0 or not water_mask.any():
        return 0.0

    k_thresh = np.percentile(K, high_percentile)
    high_k = water_mask & (K >= k_thresh)
    if not high_k.any():
        return 0.0

    labels, _ = label_connected(water_mask)
    k_components = set(np.unique(labels[high_k])) - {0}
    reachable = np.isin(labels, list(k_components))
    return float(OM[reachable].sum() / OM_total)


def habitat_fraction(K):
    """Fraction of ALL cells that are habitat at all, i.e. K > 0 -- the pores
    inside the 30-150 um physiological window (plus its upper decay tail).
    Cells with K == 0 can hold water and conduct, but nothing lives there."""
    return float((np.asarray(K) > 0).mean())


def wet_habitat_fraction(water_mask, K):
    """Stage 3.5 SS4 key diagnostic: fraction of ALL cells that are BOTH
    water-filled AND habitat (K > 0). Because small pores fill first, a soil
    whose habitat sits entirely in coarse pores can be very wet overall and
    still have a bone-dry habitat -- its biomass then sits in cells with D=0
    on every face, cut off from the substrate dissolving out of the wet fine
    matrix. This is the variable that actually decides the outcome at a fixed
    theta, and neither theta nor the PSD alone reveals it."""
    water_mask = np.asarray(water_mask, dtype=bool)
    return float((water_mask & (np.asarray(K) > 0)).mean())


def habitat_wet_share(water_mask, K):
    """The same quantity normalized by the habitat itself: what SHARE of this
    soil's habitat cells are water-filled (0.0 if the soil has no habitat)."""
    water_mask = np.asarray(water_mask, dtype=bool)
    habitat = np.asarray(K) > 0
    if not habitat.any():
        return 0.0
    return float((water_mask & habitat).sum() / habitat.sum())


def mean_substrate_in_habitat(S, K):
    """Stage 4 killer metric (soil_respiration_prompt_stage4_exploratory.md
    SS4): mean substrate concentration S in habitat cells (K>0). Stage 3.5
    showed a wet habitat is not enough -- the frozen OM rule could dump
    substrate so far from the habitat that diffusion never delivered it, so
    this is the direct check of whether carbon actually reaches the
    microbes, not merely whether their cells are wet."""
    K = np.asarray(K)
    S = np.asarray(S)
    habitat = K > 0
    if not habitat.any():
        return 0.0
    return float(S[habitat].mean())


def connectivity_diagnostics(grids, high_percentile=HIGH_PERCENTILE):
    """Bundle the percolation/connectivity diagnostics (stage-2 SS5, stage-3
    SS9's isolation_ratio, stage-3.5's wet-habitat fraction) for one soil's
    grids dict (as returned by model.build_grids)."""
    water_mask, OM, K, B0 = grids["water_mask"], grids["OM"], grids["K"], grids["B0"]
    iso_ratio, n_clusters = isolation_ratio(water_mask)
    return {
        "theta": grids["theta"],
        "largest_cluster_frac": largest_cluster_fraction(water_mask),
        "om_biomass_connectivity": om_biomass_connectivity(water_mask, OM, B0, high_percentile),
        "accessible_om_fraction": accessible_om_fraction(water_mask, OM, K, high_percentile),
        "isolation_ratio": iso_ratio,
        "n_clusters": n_clusters,
        "n_saturated_cells": int(water_mask.sum()),
        "habitat_fraction": habitat_fraction(K),
        "wet_habitat_fraction": wet_habitat_fraction(water_mask, K),
        "habitat_wet_share": habitat_wet_share(water_mask, K),
    }


def emergent_distinctness(out_A, out_B, out_C, threshold=0.3, l2_weight=0.5,
                           require_nontrivial=True):
    """Success = all three pairwise distances exceed `threshold` (Stage 1/2:
    AND every soil is non-trivial; Stage 3, soil_respiration_prompt_stage3.md
    SS8: pass `require_nontrivial=False` -- a structurally-caused die-off is
    a legitimate, strongly distinct outcome, not a disqualifier. `outcomes`
    (thriving/declining/dead per soil) is always reported as its own axis,
    regardless of `require_nontrivial`. Returns (success: bool, detail: dict)."""
    curves = {"A": out_A["R_t"], "B": out_B["R_t"], "C": out_C["R_t"]}
    outs = {"A": out_A, "B": out_B, "C": out_C}

    distances = {}
    for i, j in [("A", "B"), ("A", "C"), ("B", "C")]:
        distances[f"{i}{j}"] = pairwise_distance(curves[i], curves[j], l2_weight)

    nontrivial = {s: is_nontrivial(outs[s]) for s in "ABC"}
    outcomes = {s: classify_outcome(outs[s]) for s in "ABC"}
    shapes = {s: describe_shape(curves[s]) for s in "ABC"}

    dist_ok = all(d > threshold for d in distances.values())
    trivial_ok = all(nontrivial.values())
    success = dist_ok and (trivial_ok if require_nontrivial else True)

    return success, {
        "distances": distances, "threshold": threshold,
        "nontrivial": nontrivial, "outcomes": outcomes, "shapes": shapes,
        "min_distance": min(distances.values()),
        "require_nontrivial": require_nontrivial,
    }
