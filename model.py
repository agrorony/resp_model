"""
Pore-based soil environment respiration engine (soil_respiration_prompt_v2.md).

Every cell has a pore radius r(i,j) (pore_field.generate_pore_field). Both the
local carrying capacity K and the local diffusion rate D are FUNCTIONS of that
pore field (pore_field.K_of_r / D_of_r) -- soils differ only through their
pore-distribution params, never through K/D directly.

    dB/dt  = r_max * f(S) * B * (1 - B/K_ij) - m0*B - m_s*(1 - f(S))*B
    dS/dt  = div( D_ij * grad(S) ) + k_dis * OM_ij - uptake
    dOM/dt = -k_dis * OM_ij

    f(S)    = S / (S + Ks)
    uptake  = (1/Y) * r_max * f(S) * B * (1 - B/K_ij)

Respiration:
    R_ij(t) = (1 - Y) * [substrate used for growth in cell ij] + m0 * B_ij
    R(t)    = sum_ij R_ij(t)

D is a spatial field, so substrate transport is conservative-flux finite
volume with harmonic-mean face conductivities (correct for sharp D contrasts
and for D -> 0 at unsaturated/disconnected cells), zero-flux (Neumann)
boundaries.

Water saturation gates which cells conduct at all (percolation, stage 2/3):
`saturation.mode` in a soil config is one of:
  - fully_wet    : all pores conduct (stage 1).
  - global_theta : same theta for every soil; r_cut is THIS soil's own
                    theta-quantile of r, so cells with r <= r_cut conduct
                    (stage 2).
  - retention    : same physical r_cut (configs/mapping.yaml matric_r_cut)
                    for every soil; theta emerges from each soil's own pore
                    distribution (stage 3).

No soil-specific values are hard-coded here: everything that differs between
soils comes from the pore field built from the soil config; everything
shared (biology constants, mapping constants) comes from biology.yaml /
mapping.yaml.
"""
from __future__ import annotations

import hashlib

import numpy as np
import yaml

from pore_field import (generate_pore_field, K_of_r, D_of_r, om_field,
                        b0_field, totals)

# Hard-rule R1 (MODEL_SPEC.md SS13): biology.yaml must stay byte-identical.
BIOLOGY_SHA256 = "c76d9869dc5847f830ac3caadcce4c68986a9180875f3333504c2482a4b0bf5f"

# Hard-rule R3: no single cell may hoard OM (v2 soil C had 36% of all OM in
# one cell -- an amount artifact, not structure).
OM_MAX_OVER_MEAN = 30.0


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def check_biology_frozen(path="configs/biology.yaml"):
    with open(path, "rb") as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    assert digest == BIOLOGY_SHA256, "R1 violated: configs/biology.yaml was modified"


def harmonic_mean(a, b):
    """Face conductivity between two cells; 0 if either side is 0 (a
    disconnected/unsaturated cell truly blocks transport)."""
    out = np.zeros_like(a)
    mask = (a > 0) & (b > 0)
    out[mask] = 2.0 * a[mask] * b[mask] / (a[mask] + b[mask])
    return out


def face_conductances(D):
    """Harmonic-mean face conductivity along every axis (precomputed once)."""
    faces = []
    for ax in range(D.ndim):
        n = D.shape[ax]
        lo = tuple(slice(0, n - 1) if a == ax else slice(None) for a in range(D.ndim))
        hi = tuple(slice(1, n) if a == ax else slice(None) for a in range(D.ndim))
        faces.append((lo, hi, harmonic_mean(D[lo], D[hi])))
    return faces


def conservative_divergence(S, D, dx, faces=None):
    """div(D * grad(S)) via face-centered harmonic-mean conductivities,
    zero-flux (Neumann) boundaries, any number of dimensions. Sum over the
    whole grid is exactly zero (each face flux is added to one cell and
    subtracted from its neighbour), so this term alone conserves total S."""
    if faces is None:
        faces = face_conductances(D)
    out = np.zeros_like(S)
    for lo, hi, Df in faces:
        F = Df * (S[hi] - S[lo]) / dx
        out[lo] += F / dx
        out[hi] -= F / dx
    return out


def water_mask_and_theta(r, sat_cfg, mapping_cfg):
    mode = sat_cfg.get("mode", "fully_wet")

    if mode == "fully_wet":
        return np.ones_like(r, dtype=bool), 1.0

    if mode == "global_theta":
        theta = float(sat_cfg["theta"])
        if theta >= 1.0:
            return np.ones_like(r, dtype=bool), 1.0
        if theta <= 0.0:
            return np.zeros_like(r, dtype=bool), 0.0
        r_cut = np.quantile(r, theta)
        return r <= r_cut, theta

    if mode == "retention":
        r_cut = mapping_cfg["matric_r_cut"]
        mask = r <= r_cut
        return mask, float(mask.mean())

    raise ValueError(f"unknown saturation mode: {mode}")


def grid_shape(structure_cfg):
    g = structure_cfg["grid"]
    if "shape" in g:
        return tuple(int(s) for s in g["shape"])
    return (int(g["n"]), int(g["n"]))


def build_grids(structure_cfg, mapping_cfg, r=None):
    """Build the pore field and all derived fields (K, D, OM, B0, S0) for one
    soil. This is the ONLY place per-soil freedom enters the model. `r` may
    be passed in directly (used by the S4 shuffle control)."""
    shape = grid_shape(structure_cfg)
    dx = structure_cfg["grid"]["dx"]
    assert np.prod(shape) <= 50 ** 3, "grid too large for seconds-per-run"

    if r is None:
        r = generate_pore_field(shape, structure_cfg["pore"])

    K = K_of_r(r, mapping_cfg)
    D_full = D_of_r(r, mapping_cfg)

    sat_cfg = structure_cfg.get("saturation", {"mode": "fully_wet"})
    water_mask, theta = water_mask_and_theta(r, sat_cfg, mapping_cfg)
    D = np.where(water_mask, D_full, 0.0)

    OM_total, B0_total = totals(mapping_cfg, r.size)
    OM = om_field(r, OM_total, mapping_cfg["om_fine_bias"])
    B0 = b0_field(K, B0_total)
    S0 = np.zeros(shape)

    om_hoard = float(OM.max() / OM.mean())
    if "OM_density" in mapping_cfg:
        # R2: equal totals per unit volume (OM, B0 and habitat sum(K)).
        n = r.size
        assert abs(OM.sum() / (mapping_cfg["OM_density"] * n) - 1) < 1e-9, "R2: OM total"
        assert abs(B0.sum() / (mapping_cfg["B0_density"] * n) - 1) < 1e-9, "R2: B0 total"
        assert abs(K.sum() / (mapping_cfg["K_density"] * n) - 1) < 1e-9, "R2: sum(K)"

    return {
        "n": shape[0], "shape": shape, "dx": dx, "r": r, "K": K,
        "D_full": D_full, "D": D, "water_mask": water_mask, "theta": theta,
        "OM": OM, "B0": B0, "S0": S0, "om_hoard": om_hoard,
    }


def r3_ok(grids):
    """R3: max cell OM <= OM_MAX_OVER_MEAN x mean cell OM."""
    return grids["om_hoard"] <= OM_MAX_OVER_MEAN


def simulate(biology_cfg, structure_cfg, mapping_cfg, r=None, T=None):
    """Integrate the model for one soil. Returns time series and final fields.

    Optional v3 mechanism (`mapping_cfg["enzyme"]["enabled"]`, MODEL_SPEC.md
    SS14): OM depolymerization is catalysed by extracellular enzymes that
    biomass produces and that diffuse (slowly) through water-filled pores:

        dE/dt  = div(f_E * D * grad E) + a_E * B - d_E * E
        release = k_dis * OM * E / (E + K_E)

    With it disabled, release = k_dis * OM exactly as in v1/v2.
    """
    grids = build_grids(structure_cfg, mapping_cfg, r=r)
    dx = grids["dx"]
    D, K = grids["D"], grids["K"]
    ndim = D.ndim

    r_max = biology_cfg["r_max"]
    Ks = biology_cfg["Ks"]
    Y = biology_cfg["Y"]
    m0 = biology_cfg["m0"]
    m_s = biology_cfg["m_s"]
    k_dis = biology_cfg["k_dis"]
    dt = biology_cfg["dt"]
    T = int(T if T is not None else structure_cfg["T"])

    # R8: explicit-Euler diffusion stability limit is 1/(2*ndim).
    assert D.max() * dt / dx ** 2 < 1.0 / (2 * ndim), \
        "diffusion stability violated: reduce D, dt, or increase dx"

    enz = mapping_cfg.get("enzyme", {}) or {}
    use_enz = bool(enz.get("enabled", False))
    faces_S = face_conductances(D)
    if use_enz:
        a_E, d_E, K_E, f_E = enz["a_E"], enz["d_E"], enz["K_E"], enz["f_E"]
        faces_E = face_conductances(f_E * D)
        E = np.zeros_like(D)

    B = grids["B0"].copy()
    S = grids["S0"].copy()
    OM = grids["OM"].copy()

    R_t = np.zeros(T)
    cum_co2 = np.zeros(T)
    running = 0.0

    for t in range(T):
        f_S = S / (S + Ks)
        growth = r_max * f_S * B * (1.0 - B / K)
        uptake = (1.0 / Y) * growth

        R_field = (1.0 - Y) * growth + m0 * B
        R_t[t] = R_field.sum()
        running += R_t[t] * dt
        cum_co2[t] = running

        if use_enz:
            release = k_dis * OM * E / (E + K_E)
            dE = conservative_divergence(E, None, dx, faces_E) + a_E * B - d_E * E
            E = np.clip(E + dt * dE, 0.0, None)
        else:
            release = k_dis * OM

        dB = growth - m0 * B - m_s * (1.0 - f_S) * B
        dS = conservative_divergence(S, None, dx, faces_S) + release - uptake

        B = np.clip(B + dt * dB, 0.0, None)
        S = np.clip(S + dt * dS, 0.0, None)
        OM = np.clip(OM - dt * release, 0.0, None)

    return {
        "R_t": R_t,
        "cum_co2": cum_co2,
        "B_final": B,
        "S_final": S,
        "OM_final": OM,
        "grids": grids,
        "dt": dt,
    }
