"""
Literature-informed PARAMETRIC pore-size distributions (Stage 3.5,
soil_respiration_prompt_stage3_5.md SS2).

Stage 3 fed the model each soil's real *measured* (image-derived) PSD. Those
PSDs turned out to be poor: 30 um minimum resolvable diameter, imperfect Otsu
segmentation, and whole pore modes missing from the highest-porosity sample.
Stage 3.5 keeps the entire Stage-3 pipeline and swaps only the source of the
PSD: each soil's PSD is now a **mixture of lognormal modes over pore DIAMETER
(microns)**, with mode diameters / geometric SDs / volume weights taken from
the soil-physics literature for the three real soil types in the experiment
(Loess, Sand, Vertisol -- see configs/psd_literature.yaml for the numbers and
their rationale).

The output of `build_soil_psd` is deliberately the SAME dict shape that
`psd_data.load_soil_psd` returns for a measured PSD, so everything downstream
(pore_field.generate_diameter_field_from_psd's empirical-quantile mapping,
psd_data.habitable_volume_fraction, psd_data.median_diameter_um) works
unchanged on either kind of PSD.

These are literature-informed REPRESENTATIVE distributions for the three soil
TYPES, not measurements of the specific samples in the experiment.

Convention: a mode's `mode_um` is the peak of the volume distribution plotted
against LOG diameter -- i.e. the geometric mean / median exp(mu) of that
lognormal component, which is where the mode appears as a bump on the usual
log-x PSD plot. `gsd` is the geometric standard deviation exp(sigma_ln), so
sigma_ln = ln(gsd).
"""
from __future__ import annotations

import math

import numpy as np
import yaml


def _norm_cdf(z):
    """Standard normal CDF Phi(z) = 0.5*(1 + erf(z/sqrt(2))), numpy-only (no
    scipy): math.erf applied elementwise. The arrays here are bin-edge grids
    of a few hundred points, so the Python-level loop is irrelevant."""
    z = np.asarray(z, dtype=float)
    erf = np.array([math.erf(v / math.sqrt(2.0)) for v in z.ravel()]).reshape(z.shape)
    return 0.5 * (1.0 + erf)


def _lognormal_cdf(d, mode_um, gsd):
    """CDF of one lognormal mode at diameters `d` (microns): the volume
    fraction of THAT mode finer than d. `mode_um` is the median/log-peak
    diameter exp(mu), `gsd` the geometric SD exp(sigma_ln)."""
    sigma_ln = math.log(float(gsd))
    z = (np.log(np.asarray(d, dtype=float)) - math.log(float(mode_um))) / sigma_ln
    return _norm_cdf(z)


def build_mixture_psd(modes, d_min_um, d_max_um, n_bins):
    """Build a PSD dict from a list of lognormal modes.

    `modes` is a list of dicts with `mode_um`, `gsd`, `weight` (volume
    weight; weights are renormalized to sum to 1). Returns the same keys as
    psd_data.load_soil_psd, so the parametric and measured PSDs are
    interchangeable everywhere downstream.

    The mixture is evaluated on a log-spaced diameter grid and then
    renormalized over [d_min_um, d_max_um], so the truncated tails outside
    the grid do not leak volume: `cdf_edges` runs exactly 0 -> 1 and is a
    valid inverse-CDF (quantile) function for the pore-field mapping.
    """
    edges = np.logspace(math.log10(d_min_um), math.log10(d_max_um), int(n_bins) + 1)
    centers = np.sqrt(edges[:-1] * edges[1:])   # geometric bin centers

    weights = np.array([float(m["weight"]) for m in modes], dtype=float)
    weights = weights / weights.sum()

    cdf = np.zeros_like(edges)
    for w, m in zip(weights, modes):
        cdf += w * _lognormal_cdf(edges, m["mode_um"], m["gsd"])

    # renormalize over the truncated grid: cdf_edges[0] == 0, cdf_edges[-1] == 1
    cdf = (cdf - cdf[0]) / (cdf[-1] - cdf[0])
    cdf = np.maximum.accumulate(cdf)            # guard strict monotonicity for np.interp

    volume_fraction = np.diff(cdf)              # volume fraction per bin
    differential = volume_fraction / np.diff(edges)   # volume fraction per micron

    return {
        "bin_centers_um": centers,
        "bin_edges_um": edges,
        "volume_count": volume_fraction,        # "count" == volume fraction here
        "cumulative_porosity": cdf[1:],
        "differential_psd": differential,
        "cdf_edges": cdf,
        "modes": [dict(m) for m in modes],
    }


def load_literature_config(path="configs/psd_literature.yaml"):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def build_soil_psd(soil_name, config_path="configs/psd_literature.yaml"):
    """Build one soil's literature-informed PSD from configs/psd_literature.yaml."""
    cfg = load_literature_config(config_path)
    if soil_name not in cfg["soils"]:
        raise KeyError(f"no literature PSD for soil {soil_name!r}; "
                       f"have {sorted(cfg['soils'])}")
    soil = cfg["soils"][soil_name]
    grid = cfg["grid"]
    psd = build_mixture_psd(soil["modes"], grid["d_min_um"], grid["d_max_um"],
                            grid["n_bins"])
    psd["total_porosity"] = soil.get("total_porosity")
    psd["name"] = soil_name
    return psd


def volume_fraction_between(psd, d_lo, d_hi):
    """Volume fraction of pores with diameter in [d_lo, d_hi], read straight
    off the CDF (exact for the parametric mixture, unlike a bin-center sum)."""
    edges, cdf = psd["bin_edges_um"], psd["cdf_edges"]
    lo = float(np.interp(d_lo, edges, cdf))
    hi = float(np.interp(d_hi, edges, cdf))
    return hi - lo
