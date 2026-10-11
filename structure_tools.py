"""
Reusable structure tooling ported from the v3 branch
(claude/brave-brahmagupta-edjq14: pore_field.py, model.py, scan.py) into
main -- Stage 17 (docs/stage_results/STAGE17_RESULTS.md SS6).

Ported AS TOOLING, not as main's physics: nothing here is imported by
dual_porosity.py, and no earlier stage's result depends on it. main stays the
source of truth for data (measured/literature PSDs), physics (dual porosity,
implicit transport, matric/global saturation) and success criteria (SS4).

Two adaptations to main's conventions:
  * numpy-only (MODEL_SPEC.md "Numpy-only (from Stage 2)"): v3's
    scipy.ndimage.gaussian_filter is replaced by a separable FFT-free
    Gaussian convolution with reflect padding; v3's results are therefore not
    bit-reproducible through this module (they are through the v3 branch).
  * every fabric generator returns a standardized LATENT field z; mapping z
    to pore diameters is left to the caller (e.g. through main's own
    empirical-quantile PSD pipeline), so a measured PSD can be laid out with
    any fabric archetype without changing its marginal distribution
    (`latent_to_psd_quantiles`).

Contents:
  n-D transport   harmonic_mean / face_conductances / conservative_divergence
                  (any number of dimensions; exactly conservative)
  fabrics         gaussian (lambda, anisotropy), peds (Voronoi aggregates),
                  biopores (tubular channels), hierarchical (dual-scale)
  controls        shuffled (S4: same multiset, arrangement destroyed),
                  om_hoarding (R3 statistic)
  harness         seed_robustness / shuffle_control around ANY run function
"""
from __future__ import annotations

import itertools

import numpy as np


# ------------------------------------------------------------ n-D transport

def harmonic_mean(a, b):
    """Face conductivity; 0 if either side is 0 (a true barrier)."""
    out = np.zeros_like(a, dtype=float)
    mask = (a > 0) & (b > 0)
    out[mask] = 2.0 * a[mask] * b[mask] / (a[mask] + b[mask])
    return out


def face_conductances(D):
    """Harmonic-mean face conductivity along every axis of an n-D field."""
    faces = []
    for ax in range(D.ndim):
        n = D.shape[ax]
        lo = tuple(slice(0, n - 1) if a == ax else slice(None) for a in range(D.ndim))
        hi = tuple(slice(1, n) if a == ax else slice(None) for a in range(D.ndim))
        faces.append((lo, hi, harmonic_mean(D[lo], D[hi])))
    return faces


def conservative_divergence(S, D, dx, faces=None):
    """div(D grad S), zero-flux boundaries, any dimension. Each face flux is
    added to one cell and subtracted from its neighbour, so the grid sum is
    zero to round-off. Explicit-Euler stability: max(D)*dt/dx^2 < 1/(2*ndim)."""
    if faces is None:
        faces = face_conductances(D)
    out = np.zeros_like(S, dtype=float)
    for lo, hi, Df in faces:
        F = Df * (S[hi] - S[lo]) / dx
        out[lo] += F / dx
        out[hi] -= F / dx
    return out


# ---------------------------------------------------------------- fabrics

def _standardize(z):
    sd = z.std()
    return (z - z.mean()) / (sd if sd > 0 else 1.0)


def gaussian_smooth(x, widths):
    """Separable Gaussian smoothing, reflect boundaries, numpy only."""
    out = np.asarray(x, dtype=float)
    for ax, w in enumerate(widths):
        w = float(w)
        if w <= 0:
            continue
        r = max(1, int(4.0 * w + 0.5))
        k = np.exp(-0.5 * (np.arange(-r, r + 1) / w) ** 2)
        k /= k.sum()
        pad = [(0, 0)] * out.ndim
        pad[ax] = (r, r)
        xp = np.pad(out, pad, mode="symmetric")
        out = np.apply_along_axis(lambda v: np.convolve(v, k, mode="valid"), ax, xp)
    return out


def gaussian_latent(shape, lam, rng, aniso=1.0):
    """Correlated Gaussian field, correlation length `lam`; `aniso`
    stretches it along axis 0 (elongated/columnar domains)."""
    widths = [max(float(lam), 1e-6)] * len(shape)
    widths[0] *= float(aniso)
    return _standardize(gaussian_smooth(rng.standard_normal(shape), widths))


def _cell_centers(shape):
    return np.stack(np.meshgrid(*[np.arange(n) + 0.5 for n in shape], indexing="ij"), -1)


def peds_latent(shape, n_peds, rng, roughness=0.3):
    """Voronoi aggregates: high on ped faces (inter-aggregate macropores),
    low in ped interiors (intra-aggregate micropores)."""
    pts = rng.uniform(0, 1, (int(n_peds), len(shape))) * np.array(shape)
    d = np.sort(np.linalg.norm(_cell_centers(shape)[..., None, :] - pts, axis=-1), axis=-1)
    z = -_standardize(d[..., 1] - d[..., 0])
    return _standardize(z + roughness * gaussian_latent(shape, 1.0, rng))


def biopore_latent(shape, n_tubes, rng, background_lam=2.0, background_w=0.5):
    """Straight tubular channels, preferentially along axis 0, over a
    correlated background."""
    g = _cell_centers(shape)
    dmin = np.full(shape, np.inf)
    for _ in range(int(n_tubes)):
        p = rng.uniform(0, 1, len(shape)) * np.array(shape)
        v = rng.normal(size=len(shape))
        v[0] = abs(v[0]) * 3.0
        v /= np.linalg.norm(v)
        rel = g - p
        dmin = np.minimum(dmin, np.linalg.norm(rel - (rel @ v)[..., None] * v, axis=-1))
    z = -_standardize(dmin)
    return _standardize(z + background_w * gaussian_latent(shape, background_lam, rng))


def hierarchical_latent(shape, lam_small, lam_large, w_large, rng):
    """Fine texture nested in large-scale domains."""
    return _standardize(gaussian_latent(shape, lam_small, rng)
                        + w_large * gaussian_latent(shape, lam_large, rng))


def latent_field(shape, fabric_cfg, seed):
    """Dispatch on fabric_cfg['structure'] (gaussian|peds|biopores|hierarchical)."""
    rng = np.random.default_rng(seed)
    kind = fabric_cfg.get("structure", "gaussian")
    if kind == "gaussian":
        return gaussian_latent(shape, fabric_cfg["lambda"], rng, fabric_cfg.get("aniso", 1.0))
    if kind == "peds":
        return peds_latent(shape, fabric_cfg["n_peds"], rng, fabric_cfg.get("roughness", 0.3))
    if kind == "biopores":
        return biopore_latent(shape, fabric_cfg["n_tubes"], rng)
    if kind == "hierarchical":
        return hierarchical_latent(shape, fabric_cfg["lambda"], fabric_cfg["lambda_large"],
                                   fabric_cfg.get("w_large", 1.0), rng)
    raise ValueError(f"unknown fabric: {kind}")


def latent_to_psd_quantiles(z, bin_edges_um, cdf_edges):
    """Lay a latent field out with a GIVEN pore-size distribution: rank-map
    z through the PSD's empirical quantile function (main's psd_parametric /
    psd_measured `bin_edges_um` + `cdf_edges`), so the marginal PSD is
    exactly the data's and only the arrangement comes from the fabric."""
    u = (np.argsort(np.argsort(z.ravel(), kind="stable"), kind="stable") + 0.5) / z.size
    logd = np.interp(u, np.asarray(cdf_edges, float), np.log(np.asarray(bin_edges_um, float)))
    return np.exp(logd).reshape(z.shape)


# ---------------------------------------------------------------- controls

def shuffled(field, seed=12345):
    """S4 control: identical multiset of values, spatial arrangement destroyed."""
    return np.random.default_rng(seed).permutation(np.ravel(field)).reshape(np.shape(field))


def om_hoarding(OM, support_mask=None):
    """R3 statistic: max cell OM over the mean OM of its support."""
    OM = np.asarray(OM, float)
    sup = OM[support_mask] if support_mask is not None else OM
    return float(OM.max() / sup.mean()) if sup.size and sup.mean() > 0 else float("nan")


# ---------------------------------------------------------------- harness

def seed_robustness(run_fn, distance_fn, soils, seed_sets):
    """S3: run_fn(soil, seed) -> curve. For each seed set (one seed per soil)
    compute all pairwise between-soil distances; also the within-soil
    seed-to-seed distances. Pass iff every set's min between-soil distance
    beats `threshold` (checked by the caller) and mean within < 0.5 x min
    between (returned as `separable`)."""
    curves = {(s, i): run_fn(s, seeds[k]) for i, seeds in enumerate(seed_sets)
              for k, s in enumerate(soils)}
    between = [min(distance_fn(curves[(a, i)], curves[(b, i)])
                   for a, b in itertools.combinations(soils, 2))
               for i in range(len(seed_sets))]
    within = [distance_fn(curves[(s, i)], curves[(s, j)])
              for s in soils for i, j in itertools.combinations(range(len(seed_sets)), 2)]
    w = float(np.mean(within)) if within else float("nan")
    return {"min_between_per_set": between, "within_mean": w,
            "separable": bool(w < 0.5 * min(between))}


def shuffle_control(run_fn, distance_fn, soils):
    """S4: run_fn(soil, shuffle: bool) -> curve. Returns structured vs
    shuffled pairwise distances and the drop (arrangement contribution)."""
    out = {}
    for a, b in itertools.combinations(soils, 2):
        d_s = distance_fn(run_fn(a, False), run_fn(b, False))
        d_u = distance_fn(run_fn(a, True), run_fn(b, True))
        out[f"{a}-{b}"] = {"structured": d_s, "shuffled": d_u, "drop": d_s - d_u}
    return out
