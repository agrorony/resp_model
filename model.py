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

import numpy as np
import yaml

from pore_field import generate_pore_field, K_of_r, D_of_r, om_field, b0_field


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def harmonic_mean(a, b):
    """Face conductivity between two cells; 0 if either side is 0 (a
    disconnected/unsaturated cell truly blocks transport)."""
    out = np.zeros_like(a)
    mask = (a > 0) & (b > 0)
    out[mask] = 2.0 * a[mask] * b[mask] / (a[mask] + b[mask])
    return out


def conservative_divergence(S, D, dx):
    """div(D * grad(S)) via face-centered harmonic-mean conductivities,
    zero-flux (Neumann) boundaries. Sum over the whole grid is exactly zero
    (telescoping face fluxes), so this term alone conserves total S."""
    Dx = harmonic_mean(D[:-1, :], D[1:, :])
    Fx = Dx * (S[1:, :] - S[:-1, :]) / dx
    Fx_pad = np.pad(Fx, ((1, 1), (0, 0)))
    div_x = (Fx_pad[1:, :] - Fx_pad[:-1, :]) / dx

    Dy = harmonic_mean(D[:, :-1], D[:, 1:])
    Fy = Dy * (S[:, 1:] - S[:, :-1]) / dx
    Fy_pad = np.pad(Fy, ((0, 0), (1, 1)))
    div_y = (Fy_pad[:, 1:] - Fy_pad[:, :-1]) / dx

    return div_x + div_y


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


def build_grids(structure_cfg, mapping_cfg):
    """Build the pore field and all derived fields (K, D, OM, B0, S0) for one
    soil. This is the ONLY place per-soil freedom enters the model."""
    n = structure_cfg["grid"]["n"]
    dx = structure_cfg["grid"]["dx"]
    assert n <= 50, "grid must stay <= 50x50 (v2 invariant 3)"
    pore_cfg = structure_cfg["pore"]

    r = generate_pore_field(
        n, pore_cfg["mu"], pore_cfg["sigma"], pore_cfg["lambda"],
        pore_cfg["seed"], pore_cfg.get("aggregate", False),
    )

    K = K_of_r(r, mapping_cfg)
    D_full = D_of_r(r, mapping_cfg)

    sat_cfg = structure_cfg.get("saturation", {"mode": "fully_wet"})
    water_mask, theta = water_mask_and_theta(r, sat_cfg, mapping_cfg)
    D = np.where(water_mask, D_full, 0.0)

    OM = om_field(r, mapping_cfg["OM_total"], mapping_cfg["om_fine_bias"])
    B0 = b0_field(K, mapping_cfg["B0_total"])
    S0 = np.zeros((n, n))

    return {
        "n": n, "dx": dx, "r": r, "K": K, "D_full": D_full, "D": D,
        "water_mask": water_mask, "theta": theta,
        "OM": OM, "B0": B0, "S0": S0,
    }


def simulate(biology_cfg, structure_cfg, mapping_cfg):
    """Integrate the model for one soil. Returns time series and final fields."""
    grids = build_grids(structure_cfg, mapping_cfg)
    dx = grids["dx"]
    D, K = grids["D"], grids["K"]

    r_max = biology_cfg["r_max"]
    Ks = biology_cfg["Ks"]
    Y = biology_cfg["Y"]
    m0 = biology_cfg["m0"]
    m_s = biology_cfg["m_s"]
    k_dis = biology_cfg["k_dis"]
    dt = biology_cfg["dt"]
    T = structure_cfg["T"]

    assert D.max() * dt / dx ** 2 < 0.2, "diffusion stability violated: reduce D, dt, or increase dx"

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

        dB = growth - m0 * B - m_s * (1.0 - f_S) * B
        dS = conservative_divergence(S, D, dx) + k_dis * OM - uptake
        dOM = -k_dis * OM

        B = np.clip(B + dt * dB, 0.0, None)
        S = np.clip(S + dt * dS, 0.0, None)
        OM = np.clip(OM + dt * dOM, 0.0, None)

    return {
        "R_t": R_t,
        "cum_co2": cum_co2,
        "B_final": B,
        "S_final": S,
        "OM_final": OM,
        "grids": grids,
        "dt": dt,
    }
