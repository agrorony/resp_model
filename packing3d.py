"""
Stage 21 -- 3D particle-packing structure (Sand and Vertisol only).

Builds a cubic voxel domain (default 256^3 at 5 um = 1.28 mm) from particles:
round quartz grains (Sand) or plate-like porous clay aggregates (Vertisol).
No measured data of any kind enters: every parameter below is fixed from the
cited literature (see STRUCTURE_PARAMS and docs/stage_results/STAGE21_RESULTS.md)
BEFORE any respiration run and is never adjusted afterwards.

Every array operation goes through the array module `xp` (cupy if a GPU and
CuPy are present, else numpy; override with env RESP_BACKEND=numpy).

Voxel phases (uint8 `phase`): 0 = resolved pore, 1 = porous sub-resolution
matrix, 2 = inert solid.
"""
from __future__ import annotations

import os

import numpy as np

PORE, MATRIX, SOLID = 0, 1, 2

# ---------------------------------------------------------------- backend


def get_xp():
    if os.environ.get("RESP_BACKEND", "").lower() != "numpy":
        try:
            import cupy as cp
            if cp.cuda.runtime.getDeviceCount() > 0:
                return cp
        except Exception:
            pass
    return np


def to_np(a):
    return a.get() if hasattr(a, "get") else np.asarray(a)


def edt(mask, xp):
    """Euclidean distance (voxel units) of every True voxel to the nearest
    False voxel. The domain edge is NOT a wall (no-flux boundary only)."""
    if xp is np:
        from scipy.ndimage import distance_transform_edt
        return distance_transform_edt(mask)
    from cupyx.scipy.ndimage import distance_transform_edt
    return distance_transform_edt(mask)


def label6(mask, xp):
    """6-connected component labels."""
    if xp is np:
        from scipy.ndimage import label
        return label(mask)
    from cupyx.scipy.ndimage import label
    return label(mask)


# ------------------------------------------------- literature parameters
# All values fixed before any respiration run. "src" keys name entries of
# the reference list in STAGE21_RESULTS.md.
STRUCTURE_PARAMS = {
    "sand": {
        # Grain sizes: USDA fine-to-medium sand (fine 100-250 um, medium
        # 250-500 um; Soil Survey Manual 2017). Lognormal in phi units with
        # median 250 um (the fine/medium boundary) and Folk & Ward (1957)
        # inclusive graphic sorting sigma_phi = 0.5 (upper limit of
        # "moderately well sorted"), truncated to the fine-to-medium range.
        "grain_median_um": 250.0,
        "grain_sigma_phi": 0.5,
        "grain_min_um": 100.0,
        "grain_max_um": 500.0,
        # Grain solid fraction of the domain. Random close packing of equal
        # spheres is 0.64 (Scott 1960; Torquato et al. 2000); a sandy soil
        # at bulk density 1.6 g/cm3 with particle density 2.65 g/cm3 has
        # solid fraction 0.60. 0.60 is used.
        "grain_solid_fraction": 0.60,
        # Fines: USDA "sand" class needs silt + 1.5*clay < 15 %. Adopted:
        # silt 4 %, clay 3 % by mass -> 4 + 1.5*3 = 8.5 % (< 15). Fines mass
        # fraction f = 0.07. Volume chain (particle density rho for all
        # minerals): fines solid volume = f/(1-f) * grain volume
        # = 0.0753 * 0.60 = 0.0452 of the domain; the fines aggregate is
        # sub-resolution porous with internal porosity 0.45, so the matrix
        # occupies 0.0452/(1-0.45) = 0.0822 of the domain voxels.
        "silt_pct": 4.0,
        "clay_pct": 3.0,
        "fines_internal_porosity": 0.45,
    },
    "vertisol": {
        # Vertisols: >= 30 % clay (Soil Survey Staff, Keys to Soil Taxonomy).
        # Plate-like aggregates (lenticular / platy peds, Soil Survey Manual
        # 2017). Resolvable aggregate class: microaggregates < 250 um
        # (Tisdall & Oades 1982). Discs: major diameter lognormal, median
        # 160 um, geometric SD 1.3, truncated [80, 250] um; aspect ratio
        # (diameter / thickness) 4.
        "plate_median_um": 160.0,
        "plate_gsd": 1.3,
        "plate_min_um": 80.0,
        "plate_max_um": 250.0,
        "plate_aspect": 4.0,
        # Resolved inter-aggregate porosity 0.25 -> matrix voxel fraction
        # 0.75. Swollen Vertisol aggregates coalesce, so porous plates are
        # placed as a penetrable (Boolean) population; the union of the
        # plates is the matrix.
        "resolved_porosity": 0.25,
    },
}

PREDICTION_SEED = {"sand": 2101, "vertisol": 2102}


# --------------------------------------------------------- particle sizes

def sample_grains(rng, params, target_volume_um3):
    """Lognormal (in phi) truncated grain diameters until the summed sphere
    volume reaches target_volume_um3; the last grain is scaled so the sum
    matches exactly (nominal solid fraction)."""
    med_phi = -np.log2(params["grain_median_um"] / 1000.0)
    d, vol = [], 0.0
    while vol < target_volume_um3:
        phi = rng.normal(med_phi, params["grain_sigma_phi"])
        dum = 1000.0 * 2.0 ** (-phi)
        if not (params["grain_min_um"] <= dum <= params["grain_max_um"]):
            continue
        d.append(dum)
        vol += np.pi / 6.0 * dum ** 3
    d = np.array(d)
    over = vol - target_volume_um3
    last_v = np.pi / 6.0 * d[-1] ** 3
    d[-1] = (6.0 / np.pi * (last_v - over)) ** (1.0 / 3.0)
    return d


def pack_spheres_periodic(diam, box_um, rng, max_sweeps=20000, tol=2e-3, log=None):
    """Random periodic packing of spheres by compression + relaxation
    (Jodrey-Tory / Lubachevsky-Stillinger-type quasi-static growth): start at
    random positions with all diameters scaled by lam << 1; at each lam push
    overlapping pairs apart until max overlap < tol * mean diameter, then
    increase lam. Stops at lam = 1 (target solid fraction reached) or when
    jammed (lam_final < 1, reported). Contacts (overlap <= tol) are allowed.
    Returns positions (um), final lam, and max relative overlap."""
    n = len(diam)
    pos = rng.uniform(0.0, box_um, size=(n, 3))
    r = diam / 2.0
    rsum = r[:, None] + r[None, :]
    mean_d = float(diam.mean())
    lam, dlam = 0.55, 0.02
    iu = np.triu_indices(n, 1)
    best = (pos.copy(), lam)
    fails = 0
    while lam < 1.0 and dlam > 1e-5:
        trial = min(lam + dlam, 1.0)
        p = best[0].copy()
        ok = False
        for sweep in range(max_sweeps):
            dv = p[:, None, :] - p[None, :, :]
            dv -= box_um * np.round(dv / box_um)
            dist = np.sqrt((dv ** 2).sum(-1)) + np.eye(n)
            ov = trial * rsum - dist
            np.fill_diagonal(ov, 0.0)
            mx = ov.max()
            if mx < tol * mean_d:
                ok = True
                break
            ov = np.clip(ov, 0.0, None)
            push = (ov / dist)[:, :, None] * dv
            p = (p + 0.5 * push.sum(1)) % box_um
        if ok:
            best = (p, trial)
            lam = trial
            dlam = min(dlam * 1.5, 0.05)
        else:
            dlam *= 0.5
    pos, lam = best
    dv = pos[:, None, :] - pos[None, :, :]
    dv -= box_um * np.round(dv / box_um)
    dist = np.sqrt((dv ** 2).sum(-1)) + np.eye(n)
    ov = lam * rsum - dist
    np.fill_diagonal(ov, 0.0)
    return pos, float(lam), float(ov.max() / mean_d)


# --------------------------------------------------------- voxelization

def _patch_indices(c_um, half_um, vox, n, xp):
    c = c_um / vox
    h = int(np.ceil(half_um / vox)) + 1
    i0 = int(np.floor(c)) - h
    idx = np.arange(i0, i0 + 2 * h + 1)
    return idx, xp.asarray(idx % n), ((idx + 0.5) * vox - c_um)


def stamp_ellipsoid(vol, c_um, semi_ab_um, semi_c_um, normal, vox, xp):
    """OR an oblate spheroid (semi-axes a=b=semi_ab in the plane normal to
    `normal`, c along it) into the periodic boolean volume `vol` in place."""
    n = vol.shape[0]
    half = max(semi_ab_um, semi_c_um)
    ix, wx, px = _patch_indices(c_um[0], half, vox, n, xp)
    iy, wy, py = _patch_indices(c_um[1], half, vox, n, xp)
    iz, wz, pz = _patch_indices(c_um[2], half, vox, n, xp)
    px = xp.asarray(px)[:, None, None]
    py = xp.asarray(py)[None, :, None]
    pz = xp.asarray(pz)[None, None, :]
    nx, ny, nz = normal
    h = px * nx + py * ny + pz * nz
    rr = px * px + py * py + pz * pz
    inside = (h * h) / semi_c_um ** 2 + (rr - h * h) / semi_ab_um ** 2 <= 1.0
    sl = (wx[:, None, None], wy[None, :, None], wz[None, None, :])
    vol[sl] |= inside


def _random_normals(rng, k):
    v = rng.normal(size=(k, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


# ------------------------------------------------------ local thickness

def local_thickness(pore, vox_um, xp, step=0.25):
    """Local thickness of the pore phase (Hildebrand & Ruegsegger 1997):
    diameter of the largest inscribed sphere that contains the voxel.
    EDT-based: for radii r from large to small, voxels within r of any centre
    whose inscribed radius is >= r receive diameter 2r (if not already set).
    Inscribed radius of a centre = EDT - 0.5 (distance from voxel centre to
    the wall face), so a 1-voxel channel has thickness 1 voxel. Radius levels
    step by `step` voxels (diameter resolution 2*step voxels).
    Returns float32 diameters in um (0 on non-pore)."""
    dt = edt(pore, xp)
    R = xp.where(pore, dt, 0.0)
    lt = xp.zeros(pore.shape, dtype=xp.float32)
    rmax = float(R.max())
    if rmax <= 0:
        return lt
    top = np.ceil(rmax / step) * step
    levels = np.arange(top, 1.0 - 1e-9, -step)
    for r in levels:
        C = R >= (r - 1e-9)
        if not bool(C.any()):
            continue
        cover = edt(~C, xp) <= r
        cover &= pore
        cover &= lt == 0
        lt[cover] = 2.0 * r - 1.0
    return lt * vox_um


# ------------------------------------------------------------ structures

class Structure:
    pass


def _fill_fines(pore, lt_um, n_target, vox_um, rng_seed, xp):
    """Fill the n_target pore voxels with the smallest local thickness first
    (ties: nearest to a wall, then random) with fines matrix."""
    dt = edt(pore, xp)
    g = xp.random.RandomState(rng_seed) if xp is not np else np.random.RandomState(rng_seed)
    tie = g.random_sample(pore.shape).astype(xp.float32)
    score = xp.where(pore, lt_um + 0.01 * dt.astype(xp.float32) + 1e-3 * tie, xp.float32(1e30))
    flat = score.ravel()
    order = xp.argsort(flat)[:n_target]
    sel = xp.zeros(flat.shape, dtype=bool)
    sel[order] = True
    return sel.reshape(pore.shape)


def build_structure(kind, seed_offset=0, n=256, vox_um=5.0, xp=None, log=print):
    """kind: S (spheres + inert solid + fines), V (plates + porous aggregates),
    Sp (plates + inert solid + fines), Vp (spheres + porous aggregates),
    Val (V with the same particles, plates aligned along z)."""
    xp = xp or get_xp()
    L = n * vox_um
    N = n ** 3
    sp, vp = STRUCTURE_PARAMS["sand"], STRUCTURE_PARAMS["vertisol"]
    material = "sand" if kind in ("S", "Sp") else "vertisol"
    seed = PREDICTION_SEED[material] + 1000 * seed_offset
    rng = np.random.default_rng(seed)
    st = Structure()
    st.kind, st.n, st.vox_um, st.seed = kind, n, vox_um, seed
    particles = {}
    build_info = {}

    if material == "sand":
        target_v = sp["grain_solid_fraction"] * L ** 3
        solid = xp.zeros((n, n, n), dtype=bool)
        if kind == "S":
            d = sample_grains(rng, sp, target_v)
            pos, lam, ovmax = pack_spheres_periodic(d, L, rng)
            build_info.update(n_grains=len(d), lam_final=lam, max_rel_overlap=ovmax,
                              nominal_solid_fraction=sp["grain_solid_fraction"] * lam ** 3)
            for c, dd in zip(pos, d):
                stamp_ellipsoid(solid, c, lam * dd / 2, lam * dd / 2, (0, 0, 1.0), vox_um, xp)
        else:
            # Sp: equal-volume plates (same grain-volume distribution), random
            # orientation, penetrable placement to the same solid fraction.
            A = vp["plate_aspect"]
            count = 0
            target_vox = sp["grain_solid_fraction"] * N
            med_phi = -np.log2(sp["grain_median_um"] / 1000.0)
            while int(solid.sum()) < target_vox:
                phi = rng.normal(med_phi, sp["grain_sigma_phi"])
                dd = 1000.0 * 2.0 ** (-phi)
                if not (sp["grain_min_um"] <= dd <= sp["grain_max_um"]):
                    continue
                D = dd * A ** (1.0 / 3.0)
                c = rng.uniform(0, L, 3)
                nrm = _random_normals(rng, 1)[0]
                stamp_ellipsoid(solid, c, D / 2, D / (2 * A), nrm, vox_um, xp)
                count += 1
            build_info.update(n_plates=count)
        pore = ~solid
        lt0 = local_thickness(pore, vox_um, xp)
        n_fines = int(round(fines_volume_fraction() * N))
        fines = _fill_fines(pore, lt0, n_fines, vox_um, seed, xp)
        phase = xp.zeros((n, n, n), dtype=xp.uint8)
        phase[solid] = SOLID
        phase[fines] = MATRIX
        build_info["fines_volume_fraction_target"] = fines_volume_fraction()
    else:
        A = vp["plate_aspect"]
        target_vox = (1.0 - vp["resolved_porosity"]) * N
        med = vp["plate_median_um"]
        sg = np.log(vp["plate_gsd"])
        plist = []
        mat_v = xp.zeros((n, n, n), dtype=bool)
        mat_s = xp.zeros((n, n, n), dtype=bool)
        mat_al = xp.zeros((n, n, n), dtype=bool)
        # Population of V is generated once; V', V_al re-use its centres/sizes.
        cnt = 0
        rr = np.random.default_rng(seed)
        while True:
            D = float(np.exp(np.log(med) + sg * rr.normal()))
            if not (vp["plate_min_um"] <= D <= vp["plate_max_um"]):
                continue
            c = rr.uniform(0, L, 3)
            nrm = _random_normals(rr, 1)[0]
            plist.append((c, D, nrm))
            stamp_ellipsoid(mat_v, c, D / 2, D / (2 * A), nrm, vox_um, xp)
            cnt += 1
            if cnt % 5 == 0 and int(mat_v.sum()) >= target_vox:
                break
        # trim to the first particle count at which the target was crossed
        build_info.update(n_plates=len(plist))
        if kind == "V":
            m = mat_v
        elif kind == "Val":
            m = mat_al
            for c, D, nrm in plist:
                stamp_ellipsoid(m, c, D / 2, D / (2 * A), (0.0, 0.0, 1.0), vox_um, xp)
        else:  # Vp
            m = mat_s
            for c, D, nrm in plist:
                de = D / A ** (1.0 / 3.0)
                stamp_ellipsoid(m, c, de / 2, de / 2, (0.0, 0.0, 1.0), vox_um, xp)
        phase = xp.zeros((n, n, n), dtype=xp.uint8)
        phase[m] = MATRIX
    st.phase = phase
    st.info = build_info
    return st


def fines_volume_fraction():
    sp = STRUCTURE_PARAMS["sand"]
    f = (sp["silt_pct"] + sp["clay_pct"]) / 100.0
    assert sp["silt_pct"] + 1.5 * sp["clay_pct"] < 15.0, "USDA sand class"
    solid_fines = f / (1.0 - f) * sp["grain_solid_fraction"]
    return solid_fines / (1.0 - sp["fines_internal_porosity"])


# --------------------------------------------------------------- reports

def geodesic_within(matrix, pore, steps, xp):
    """Matrix voxels within `steps` 6-connected voxel steps of a pore voxel,
    walking through matrix voxels only."""
    def nb(m):
        o = m.copy()
        o[1:] |= m[:-1]; o[:-1] |= m[1:]
        o[:, 1:] |= m[:, :-1]; o[:, :-1] |= m[:, 1:]
        o[:, :, 1:] |= m[:, :, :-1]; o[:, :, :-1] |= m[:, :, 1:]
        return o
    reached = nb(pore) & matrix
    for _ in range(steps - 1):
        reached |= nb(reached) & matrix
    return reached


def pore_axis_span(cluster, xp):
    out = []
    for ax in range(3):
        lo = xp.take(cluster, 0, axis=ax).any()
        hi = xp.take(cluster, -1, axis=ax).any()
        out.append(bool(lo and hi))
    return out


def analyse(st, xp, d_cut_um=50.0, hab_min_um=10.0, k_hi_um=150.0, geo_um=20.0):
    """Geometry report (no respiration)."""
    phase, vox = st.phase, st.vox_um
    N = phase.size
    pore = phase == PORE
    matrix = phase == MATRIX
    solid = phase == SOLID
    lt = local_thickness(pore, vox, xp)
    lab, nlab = label6(pore, xp)
    counts = xp.bincount(lab.ravel())
    counts[0] = 0
    big = int(counts.argmax()) if nlab > 0 else 0
    cluster = lab == big
    vals = lt[pore]
    vals_np = to_np(vals)
    q = {f"d{p}": float(np.percentile(vals_np, p)) for p in (10, 50, 90)} if vals_np.size else {}
    nonsolid = int((~solid).sum())
    wet_pore = pore & (lt < d_cut_um)
    wet = wet_pore | matrix
    S_w = float(wet.sum()) / nonsolid
    hab = pore & (lt >= hab_min_um)
    wet_hab = hab & (lt < d_cut_um)
    geo = geodesic_within(matrix, pore, int(round(geo_um / vox)), xp)
    nm = int(matrix.sum())
    hist_edges = np.array([0, 5, 10, 20, 30, 50, 75, 100, 150, 250, 1e9])
    hist = np.histogram(vals_np, bins=hist_edges)[0] / max(vals_np.size, 1)
    rep = {
        "kind": st.kind, "n": st.n, "vox_um": vox,
        "porosity_resolved": float(pore.sum()) / N,
        "matrix_fraction": nm / N,
        "solid_fraction": float(solid.sum()) / N,
        "LT_percentiles_um": q,
        "LT_mean_um": float(vals_np.mean()) if vals_np.size else 0.0,
        "LT_hist_edges_um": hist_edges[:-1].tolist() + ["inf"],
        "LT_hist_volume_frac": hist.tolist(),
        "largest_cluster_fraction": float(counts.max()) / max(int(pore.sum()), 1),
        "n_pore_clusters": int(nlab),
        "spans_axes_xyz": pore_axis_span(cluster, xp),
        "S_w": S_w,
        "habitat_volume_fraction": float(hab.sum()) / N,
        "habitat_in_spanning_cluster": float((hab & cluster).sum()) / max(int(hab.sum()), 1),
        "wet_fraction_of_habitat": float(wet_hab.sum()) / max(int(hab.sum()), 1),
        "OM_within_20um_geodesic": float(geo.sum()) / max(nm, 1),
        "build": {k: (float(v) if isinstance(v, (np.floating, float)) else v)
                  for k, v in st.info.items()},
    }
    st.lt = lt
    return rep
