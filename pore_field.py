"""
Pore-size field generator and pore -> parameter maps (soil_respiration_prompt_v2.md SA-B).

Every cell gets a pore radius r(i,j) from a spatially-correlated lognormal
field. K(r) and D(r) are then FUNCTIONS of that field using constants fixed
in configs/mapping.yaml -- soils differ only through (mu, sigma, lambda,
aggregate), never through K/D directly.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter


def generate_pore_field(n, mu, sigma, lam, seed, aggregate=False):
    """Spatially-correlated lognormal pore-radius field r(i,j).

    Recipe: draw white noise, convolve with a Gaussian kernel of width
    `lambda` (the correlation length), standardize to mean 0 / std 1, then
    map through the lognormal quantile r = exp(mu + sigma*z) -- equivalent to
    passing the standardized field through the lognormal CDF/quantile since z
    is already standard normal by construction. `aggregate` additionally
    thresholds the high-gradient (domain-boundary) cells down to fine pores,
    producing peds separated by fine-pore boundaries.
    """
    rng = np.random.default_rng(seed)
    noise = rng.standard_normal((n, n))
    lam = max(float(lam), 1e-6)
    smooth = gaussian_filter(noise, sigma=lam, mode="reflect")
    z = (smooth - smooth.mean()) / smooth.std()
    r = np.exp(mu + sigma * z)

    if aggregate:
        gy, gx = np.gradient(z)
        grad_mag = np.hypot(gy, gx)
        boundary = grad_mag > np.percentile(grad_mag, 75)
        fine_r = np.percentile(r, 10)
        r = np.where(boundary, np.minimum(r, fine_r), r)

    return r


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


def om_field(r, total, fine_bias):
    """Rule-based OM: protected in fine pores, weight = r^-fine_bias,
    renormalized to a fixed total. Same rule/exponent for every soil, so
    layout differs only because the pore field differs."""
    weight = r ** (-fine_bias)
    return total * weight / weight.sum()


def b0_field(K, total):
    """Rule-based initial biomass: proportional to local K, so the
    community seeds itself wherever habitat is good, not at typed-in
    centers."""
    if K.sum() <= 0:
        return np.zeros_like(K)
    return total * K / K.sum()
