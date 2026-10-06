"""
Pore-size field generator and pore -> parameter maps (soil_respiration_prompt_v2.md SA-B;
Stage 3 empirical-PSD mapping and physiological K/D, soil_respiration_prompt_stage3.md SS3-5).

Stage 1/2: every cell gets a pore radius r(i,j) from a spatially-correlated
lognormal field. Stage 3: every cell gets a pore DIAMETER d(i,j) in microns
by mapping the same kind of spatially-correlated field through a soil's real
measured pore-size distribution (empirical-quantile mapping) instead of a
lognormal quantile function. Either way, K and D are FUNCTIONS of that field
using constants fixed in configs/mapping.yaml -- soils differ only through
their pore field, never through K/D directly.
"""
from __future__ import annotations

import numpy as np


def _gaussian_kernel1d(sigma, truncate=4.0):
    """1D Gaussian kernel, matching scipy.ndimage.gaussian_filter's default
    (`radius = int(truncate*sigma + 0.5)`, normalized to sum 1)."""
    radius = int(truncate * sigma + 0.5)
    x = np.arange(-radius, radius + 1)
    kernel = np.exp(-0.5 * (x / sigma) ** 2)
    return kernel / kernel.sum()


def gaussian_filter_2d(a, sigma, truncate=4.0):
    """Separable 2D Gaussian smoothing (1D kernel applied along each axis)
    with reflect padding, replacing scipy.ndimage.gaussian_filter so the
    whole repo stays numpy-only. Note: scipy's mode="reflect" duplicates the
    edge value (d c b a | a b c d), which is numpy.pad's "symmetric" mode --
    numpy's own "reflect" mode does NOT duplicate the edge and would not
    match scipy's default."""
    kernel = _gaussian_kernel1d(sigma, truncate)
    radius = (len(kernel) - 1) // 2

    def _convolve_axis(arr, axis):
        padded = np.pad(arr, [(radius, radius) if ax == axis else (0, 0)
                               for ax in range(arr.ndim)], mode="symmetric")
        out = np.zeros_like(arr, dtype=float)
        for i, w in enumerate(kernel):
            sl = [slice(None)] * arr.ndim
            sl[axis] = slice(i, i + arr.shape[axis])
            out += w * padded[tuple(sl)]
        return out

    return _convolve_axis(_convolve_axis(a, axis=0), axis=1)


def _correlated_gaussian_z(n, lam, seed):
    """Spatially-correlated field standardized to mean 0 / std 1: draw white
    noise, convolve with a Gaussian kernel of width `lambda` (the
    correlation length), standardize. Shared by the lognormal (Stage 1/2)
    and empirical-PSD (Stage 3) pore-field generators below -- only the
    quantile function applied to `z` differs."""
    rng = np.random.default_rng(seed)
    noise = rng.standard_normal((n, n))
    lam = max(float(lam), 1e-6)
    smooth = gaussian_filter_2d(noise, sigma=lam)
    return (smooth - smooth.mean()) / smooth.std()


def _apply_aggregate_boundaries(z, field):
    """Threshold the high-gradient (domain-boundary) cells of `field` down
    to a fine value (its own 10th percentile), producing peds separated by
    fine-pore boundaries. `field` may hold radii (Stage 1/2) or diameters
    (Stage 3) -- the rule ("boundaries are fine") is the same either way."""
    gy, gx = np.gradient(z)
    grad_mag = np.hypot(gy, gx)
    boundary = grad_mag > np.percentile(grad_mag, 75)
    fine_val = np.percentile(field, 10)
    return np.where(boundary, np.minimum(field, fine_val), field)


def generate_pore_field(n, mu, sigma, lam, seed, aggregate=False):
    """Spatially-correlated lognormal pore-radius field r(i,j) (Stage 1/2).

    Map the standardized correlated field through the lognormal quantile
    r = exp(mu + sigma*z) -- equivalent to passing the standardized field
    through the lognormal CDF/quantile since z is already standard normal by
    construction. `aggregate` additionally carves fine-pore ped boundaries.
    """
    z = _correlated_gaussian_z(n, lam, seed)
    r = np.exp(mu + sigma * z)
    if aggregate:
        r = _apply_aggregate_boundaries(z, r)
    return r


def rank_uniform(z):
    """Empirical uniform quantile of `z` via a rank transform (numpy-only:
    argsort-based, no scipy.stats.norm.cdf/erf needed). Exact for the
    field's own empirical distribution regardless of how close it is to
    truly Gaussian -- which is what "convert to uniform quantiles" means for
    an empirical-quantile mapping (soil_respiration_prompt_stage3.md SS3)."""
    flat = z.ravel()
    order = np.argsort(flat)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(flat.size)
    u = (ranks + 0.5) / flat.size
    return u.reshape(z.shape)


def generate_diameter_field_from_psd(n, lam, seed, bin_edges_um, cdf_edges, aggregate=False):
    """Spatially-correlated pore-DIAMETER field d(i,j) in microns (Stage 3),
    via empirical-quantile mapping: same correlated-Gaussian recipe as
    `generate_pore_field`, converted to uniform quantiles (`rank_uniform`),
    then mapped through THIS soil's empirical inverse-CDF (`bin_edges_um`,
    `cdf_edges` from psd_data.load_soil_psd) by linear interpolation within
    the measured bins. The field's marginal diameter distribution then
    matches the soil's measured PSD; only the spatial *arrangement* (the
    correlated-Gaussian recipe) is synthetic and shared across soils."""
    z = _correlated_gaussian_z(n, lam, seed)
    u = rank_uniform(z)
    d = np.interp(u, cdf_edges, bin_edges_um)
    if aggregate:
        d = _apply_aggregate_boundaries(z, d)
    return d


def K_of_r(r, mapping_cfg):
    """K(r) = K_max * exp(-0.5*((ln r - ln r_opt)/w_K)^2), floored at k_min."""
    r_opt = mapping_cfg["r_opt"]
    w_K = mapping_cfg["w_K"]
    K_max = mapping_cfg["K_max"]
    k_min = mapping_cfg.get("k_min", 0.005)
    hump = np.exp(-0.5 * ((np.log(r) - np.log(r_opt)) / w_K) ** 2)
    return k_min + (K_max - k_min) * hump


def D_of_r(r, mapping_cfg):
    """D(r) = D_min + (D_max - D_min) * (r / (r + r_ref))^p -- rises with
    pore size and saturates at D_max as r -> infinity (a literal power law
    (r/r_ref)^p would be unbounded and could blow past the numerical
    stability cap for a coarse/large-mu soil)."""
    D_min = mapping_cfg["D_min"]
    D_max = mapping_cfg["D_max"]
    r_ref = mapping_cfg["r_ref"]
    p = mapping_cfg["p"]
    return D_min + (D_max - D_min) * (r / (r + r_ref)) ** p


def K_window(d, mapping_cfg):
    """Stage 3 physiological carrying capacity K(d), d = pore diameter in
    microns (soil_respiration_prompt_stage3.md SS4): hard zero below
    `K_d_low` (too tight for microbial activity), full plateau at `K_max`
    across the microbially-relevant window [`K_d_low`, `K_d_high`], smooth
    Gaussian decay above `K_d_high` toward ~0 (width `K_w_decay`), no floor.
    Cells with K=0 are non-habitat -- they still hold water and conduct
    (conduits, not homes)."""
    d_lo = mapping_cfg["K_d_low"]
    d_hi = mapping_cfg["K_d_high"]
    w = mapping_cfg["K_w_decay"]
    K_max = mapping_cfg["K_max"]
    return np.where(
        d < d_lo, 0.0,
        np.where(d <= d_hi, K_max, K_max * np.exp(-((d - d_hi) / w) ** 2)),
    )


def D_of_d(d, mapping_cfg):
    """Stage 3 D(d), d = pore diameter in microns (soil_respiration_prompt_stage3.md
    SS5): same saturating form as `D_of_r`, rescaled to a diameter reference
    `D_d_ref` (a physically sensible mid-window diameter) instead of the
    Stage 1/2 radius-based `r_ref`. Reuses the shared `D_min`/`D_max`/`p`."""
    D_min = mapping_cfg["D_min"]
    D_max = mapping_cfg["D_max"]
    d_ref = mapping_cfg["D_d_ref"]
    p = mapping_cfg["p"]
    return D_min + (D_max - D_min) * (d / (d + d_ref)) ** p


def om_field(r, total, fine_bias):
    """Rule-based OM: protected in fine pores, weight = r^-fine_bias,
    renormalized to a fixed total. Same rule/exponent for every soil, so
    layout differs only because the pore field differs."""
    weight = r ** (-fine_bias)
    return total * weight / weight.sum()


def truncate_psd_floor(psd, floor_um, n_bins=None):
    """Stage 4 Variant A (soil_respiration_prompt_stage4_exploratory.md SS2.1):
    drop all pore volume below `floor_um` and renormalize, so the finest cell
    representable on the grid is `floor_um` -- sub-floor clay is simply not on
    the grid, removing the sub-micron OM traps by construction.

    Resamples the ORIGINAL empirical CDF (`psd["bin_edges_um"]`,
    `psd["cdf_edges"]`) onto a fresh log-spaced grid from `floor_um` to the
    original max diameter, so the result is a clean PSD dict of the same shape
    regardless of whether the input was a measured or parametric PSD (no
    partial-bin bookkeeping). Returns `(new_psd, report)`, where `report` has
    the dropped volume fraction and before/after total porosity."""
    edges, cdf = psd["bin_edges_um"], psd["cdf_edges"]
    d_max = float(edges[-1])
    if n_bins is None:
        n_bins = len(edges) - 1

    if floor_um <= edges[0]:
        new_edges = edges.copy()
    else:
        new_edges = np.logspace(np.log10(floor_um), np.log10(d_max), n_bins + 1)

    raw_cdf = np.interp(new_edges, edges, cdf)
    dropped_fraction = float(raw_cdf[0])
    new_cdf = (raw_cdf - raw_cdf[0]) / (raw_cdf[-1] - raw_cdf[0])
    new_cdf = np.maximum.accumulate(new_cdf)

    centers = np.sqrt(new_edges[:-1] * new_edges[1:])
    vol_frac = np.diff(new_cdf)
    differential = vol_frac / np.diff(new_edges)

    new_psd = dict(psd)
    new_psd.update({
        "bin_centers_um": centers, "bin_edges_um": new_edges,
        "volume_count": vol_frac, "cumulative_porosity": new_cdf[1:],
        "differential_psd": differential, "cdf_edges": new_cdf,
    })

    porosity_before = psd.get("total_porosity")
    porosity_after = (porosity_before * (1.0 - dropped_fraction)
                      if porosity_before is not None else None)
    report = {
        "floor_um": floor_um,
        "dropped_volume_fraction": dropped_fraction,
        "porosity_before": porosity_before,
        "porosity_after": porosity_after,
    }
    return new_psd, report


def b0_field(K, total):
    """Rule-based initial biomass: proportional to local K, so the
    community seeds itself wherever habitat is good, not at typed-in
    centers."""
    if K.sum() <= 0:
        return np.zeros_like(K)
    return total * K / K.sum()
