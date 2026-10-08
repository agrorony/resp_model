"""
Pore-size field generator and pore -> parameter maps (soil_respiration_prompt_v2.md SA-B,
extended in v3 to n-D grids, structure archetypes and equal-totals rules).

Every cell gets a pore radius r from a spatially-structured lognormal field.
K(r) and D(r) are then FUNCTIONS of that field using constants fixed in
configs/mapping.yaml -- soils differ only through their pore field (texture
mu/sigma + spatial arrangement), never through K/D directly.

v3 equal-totals rule (MODEL_SPEC.md SS13): total OM, total initial biomass AND
total carrying capacity sum(K) are fixed per unit volume and identical for
every soil. K(r) keeps its shape (where habitat is good), but is renormalized
to the shared total exactly like OM and B0 already were in v2.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter


# ---------------------------------------------------------------- latent field

def _standardize(z):
    return (z - z.mean()) / z.std()


def gaussian_latent(shape, lam, rng, aniso=1.0):
    """Correlated Gaussian random field: white noise smoothed by a Gaussian
    kernel of width `lambda` (correlation length). `aniso` stretches the
    kernel along axis 0 (the vertical axis in 3D), giving layered fabric."""
    noise = rng.standard_normal(shape)
    lam = max(float(lam), 1e-6)
    widths = [lam] * len(shape)
    widths[0] = lam * float(aniso)
    return _standardize(gaussian_filter(noise, sigma=widths, mode="reflect"))


def _cell_centers(shape):
    return np.stack(np.meshgrid(*[np.arange(n) + 0.5 for n in shape], indexing="ij"), -1)


def peds_latent(shape, n_peds, rng, roughness=0.3):
    """Aggregated fabric: Voronoi peds. Latent value is high on ped faces
    (inter-aggregate macropores) and low in ped interiors (intra-aggregate
    micropores), plus a little small-scale roughness."""
    pts = rng.uniform(0, 1, (int(n_peds), len(shape))) * np.array(shape)
    g = _cell_centers(shape)
    d = np.sort(np.linalg.norm(g[..., None, :] - pts, axis=-1), axis=-1)
    face_dist = d[..., 1] - d[..., 0]          # 0 exactly on a ped face
    z = -_standardize(face_dist)
    return _standardize(z + roughness * gaussian_latent(shape, 1.0, rng))


def biopore_latent(shape, n_tubes, rng, background_lam=2.0, background_w=0.5):
    """Biopore fabric: straight tubular macropores (root/earthworm channels),
    preferentially vertical (axis 0), over a correlated background matrix."""
    g = _cell_centers(shape)
    dmin = np.full(shape, np.inf)
    for _ in range(int(n_tubes)):
        p = rng.uniform(0, 1, len(shape)) * np.array(shape)
        v = rng.normal(size=len(shape))
        v[0] = abs(v[0]) * 3.0
        v /= np.linalg.norm(v)
        rel = g - p
        d = np.linalg.norm(rel - (rel @ v)[..., None] * v, axis=-1)
        dmin = np.minimum(dmin, d)
    z = -_standardize(dmin)
    return _standardize(z + background_w * gaussian_latent(shape, background_lam, rng))


def hierarchical_latent(shape, lam_small, lam_large, w_large, rng):
    """Dual-scale fabric: fine-scale texture nested in large-scale domains."""
    return _standardize(gaussian_latent(shape, lam_small, rng)
                        + w_large * gaussian_latent(shape, lam_large, rng))


def latent_field(shape, pore_cfg):
    """Dispatch on `pore.structure` (default 'gaussian')."""
    rng = np.random.default_rng(pore_cfg["seed"])
    kind = pore_cfg.get("structure", "gaussian")
    if kind == "gaussian":
        return gaussian_latent(shape, pore_cfg["lambda"], rng, pore_cfg.get("aniso", 1.0))
    if kind == "peds":
        return peds_latent(shape, pore_cfg["n_peds"], rng, pore_cfg.get("roughness", 0.3))
    if kind == "biopores":
        return biopore_latent(shape, pore_cfg["n_tubes"], rng)
    if kind == "hierarchical":
        return hierarchical_latent(shape, pore_cfg["lambda"], pore_cfg["lambda_large"],
                                   pore_cfg.get("w_large", 1.0), rng)
    raise ValueError(f"unknown pore structure: {kind}")


def generate_pore_field(shape, pore_cfg):
    """Spatially-structured lognormal pore-radius field r = exp(mu + sigma*z).

    `aggregate: true` (v2 option, kept for gaussian fabric) additionally
    thresholds the high-gradient (domain-boundary) cells down to fine pores.
    """
    if isinstance(shape, int):
        shape = (shape, shape)
    shape = tuple(int(s) for s in shape)
    z = latent_field(shape, pore_cfg)
    r = np.exp(pore_cfg["mu"] + pore_cfg["sigma"] * z)

    if pore_cfg.get("aggregate", False):
        grad_mag = np.sqrt(sum(g ** 2 for g in np.gradient(z)))
        boundary = grad_mag > np.percentile(grad_mag, 75)
        fine_r = np.percentile(r, 10)
        r = np.where(boundary, np.minimum(r, fine_r), r)

    return r


def shuffled(r, seed=12345):
    """Same pore-size multiset, spatial arrangement destroyed (S4 control)."""
    rng = np.random.default_rng(seed)
    return rng.permutation(r.ravel()).reshape(r.shape)


# ------------------------------------------------------------ pore -> params

def K_shape_of_r(r, mapping_cfg):
    """K(r) = k_min + (K_max - k_min) * exp(-0.5*((ln r - ln r_opt)/w_K)^2)."""
    r_opt = mapping_cfg["r_opt"]
    w_K = mapping_cfg["w_K"]
    K_max = mapping_cfg["K_max"]
    k_min = mapping_cfg.get("k_min", 0.005)
    hump = np.exp(-0.5 * ((np.log(r) - np.log(r_opt)) / w_K) ** 2)
    return k_min + (K_max - k_min) * hump


def K_of_r(r, mapping_cfg):
    """Carrying capacity. With `K_density` set (v3 equal-totals rule), the
    K(r) shape is renormalized so sum(K) = K_density * n_cells for every
    soil; without it (v2 configs) the raw shape is returned."""
    K = K_shape_of_r(r, mapping_cfg)
    if "K_density" in mapping_cfg:
        K = K * mapping_cfg["K_density"] * K.size / K.sum()
    return K


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


def totals(mapping_cfg, n_cells):
    """Shared OM / B0 totals. v3: per-cell densities x n_cells (so a bigger or
    3D grid keeps the same amount per unit volume); v2: fixed totals."""
    if "OM_density" in mapping_cfg:
        return mapping_cfg["OM_density"] * n_cells, mapping_cfg["B0_density"] * n_cells
    return mapping_cfg["OM_total"], mapping_cfg["B0_total"]
