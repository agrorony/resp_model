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
  - matric_d     : same air-entry DIAMETER d_cut in microns for every soil
                    (a shared matric potential); theta emerges per soil from
                    its own PSD (stage 3.5 diagnostic).

No soil-specific values are hard-coded here: everything that differs between
soils comes from the pore field built from the soil config; everything
shared (biology constants, mapping constants) comes from biology.yaml /
mapping.yaml.
"""
from __future__ import annotations

import numpy as np
import yaml

from pore_field import (
    generate_pore_field, K_of_r, D_of_r,
    generate_diameter_field_from_psd, K_window, D_of_d,
    om_field, b0_field, truncate_psd_floor,
)
import psd_data
import psd_parametric


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

    if mode == "matric_d":
        # Stage 3.5 diagnostic: a shared MATRIC POTENTIAL rather than a shared
        # theta. All soils sit at the same air-entry diameter `d_cut` (microns)
        # -- every pore finer than d_cut is water-filled -- and each soil's
        # theta then EMERGES from its own PSD. This, not a shared theta, is what
        # a closed humid jar actually imposes on three different soils; it is
        # reported alongside the shared-theta headline run, never instead of it.
        d_cut = float(sat_cfg.get("d_cut", mapping_cfg["matric_d_cut"]))
        mask = r <= d_cut
        return mask, float(mask.mean())

    raise ValueError(f"unknown saturation mode: {mode}")


def build_grids(structure_cfg, mapping_cfg):
    """Build the pore field and all derived fields (K, D, OM, B0, S0) for one
    soil. This is the ONLY place per-soil freedom enters the model.

    `pore.mode` selects the pore-field generator: "lognormal" (Stage 1/2,
    default -- `pore.mu`/`pore.sigma`), "psd" (Stage 3 -- `pore.psd_dir`
    points at a soil's real measured PSD, see psd_data.load_soil_psd), or
    "psd_parametric" (Stage 3.5 -- `pore.psd_soil` names a soil in
    `pore.psd_config` = configs/psd_literature.yaml, a literature-informed
    lognormal mixture, see psd_parametric.build_soil_psd). Both PSD modes put
    pore DIAMETER in microns on the field (not radius) and share the same
    empirical-quantile mapping and the same K/D maps (K_window/D_of_d); the
    lognormal mode uses K_of_r/D_of_r.

    Stage 4 Variant A (soil_respiration_prompt_stage4_exploratory.md SS2) adds
    OPTIONAL `pore` keys in PSD modes: `psd_truncate_floor_um` (drop pore
    volume below this floor and renormalize, pore_field.truncate_psd_floor)
    and `wet_diffusion_scale` (multiply D only in water-filled cells). Neither
    key is set by Stage 1-3.5 configs, so this is a no-op for them."""
    n = structure_cfg["grid"]["n"]
    dx = structure_cfg["grid"]["dx"]
    assert n <= 200, "grid must stay <= 200x200 (stage 3 invariant, soil_respiration_prompt_stage3.md SS7)"
    pore_cfg = structure_cfg["pore"]
    mode = pore_cfg.get("mode", "lognormal")

    psd_truncation = None
    if mode in ("psd", "psd_parametric"):
        if mode == "psd":
            psd = psd_data.load_soil_psd(pore_cfg["psd_dir"])
        else:
            psd = psd_parametric.build_soil_psd(
                pore_cfg["psd_soil"],
                pore_cfg.get("psd_config", "configs/psd_literature.yaml"),
            )
        floor_um = pore_cfg.get("psd_truncate_floor_um")
        if floor_um is not None:
            psd, psd_truncation = truncate_psd_floor(psd, float(floor_um))
        r = generate_diameter_field_from_psd(
            n, pore_cfg["lambda"], pore_cfg["seed"],
            psd["bin_edges_um"], psd["cdf_edges"],
            pore_cfg.get("aggregate", False),
        )
        K = K_window(r, mapping_cfg)
        D_full = D_of_d(r, mapping_cfg)
    elif mode == "lognormal":
        r = generate_pore_field(
            n, pore_cfg["mu"], pore_cfg["sigma"], pore_cfg["lambda"],
            pore_cfg["seed"], pore_cfg.get("aggregate", False),
        )
        K = K_of_r(r, mapping_cfg)
        D_full = D_of_r(r, mapping_cfg)
    else:
        raise ValueError(f"unknown pore.mode: {mode}")

    sat_cfg = structure_cfg.get("saturation", {"mode": "fully_wet"})
    water_mask, theta = water_mask_and_theta(r, sat_cfg, mapping_cfg)
    wet_diffusion_scale = float(pore_cfg.get("wet_diffusion_scale", 1.0))
    D = np.where(water_mask, D_full * wet_diffusion_scale, 0.0)

    OM = om_field(r, mapping_cfg["OM_total"], mapping_cfg["om_fine_bias"])
    B0 = b0_field(K, mapping_cfg["B0_total"])
    S0 = np.zeros((n, n))

    return {
        "n": n, "dx": dx, "r": r, "K": K, "D_full": D_full, "D": D,
        "water_mask": water_mask, "theta": theta,
        "OM": OM, "B0": B0, "S0": S0, "psd_truncation": psd_truncation,
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

    # The Stage-3 physiological K(d) is a HARD zero below K_d_low (30 um): those
    # cells are non-habitat -- they hold water and conduct, but nothing lives
    # there. The logistic term 1 - B/K is then 0/0 there, so it must be masked,
    # not divided. (Stage 3's measured PSDs bottomed out at 30 um, so K was
    # never exactly 0 and this never fired; the Stage-3.5 literature PSDs go
    # down to clay micropores, so most of the Vertisol's grid is K == 0.)
    # K is fixed in time, so the mask is built once.
    habitable = K > 0
    K_safe = np.where(habitable, K, 1.0)

    B = grids["B0"].copy()
    S = grids["S0"].copy()
    OM = grids["OM"].copy()

    R_t = np.zeros(T)
    cum_co2 = np.zeros(T)
    S_habitat_mean_t = np.zeros(T)
    running = 0.0

    for t in range(T):
        f_S = S / (S + Ks)
        logistic = np.where(habitable, 1.0 - B / K_safe, 0.0)
        growth = r_max * f_S * B * logistic
        uptake = (1.0 / Y) * growth

        R_field = (1.0 - Y) * growth + m0 * B
        R_t[t] = R_field.sum()
        running += R_t[t] * dt
        cum_co2[t] = running
        # Stage-4 mechanism diagnostic (soil_respiration_prompt_stage4_
        # exploratory.md SS4): mean substrate concentration reaching the
        # habitat, not just whether it is wet -- the Stage-3.5 killer metric.
        S_habitat_mean_t[t] = float(S[habitable].mean()) if habitable.any() else 0.0

        dB = growth - m0 * B - m_s * (1.0 - f_S) * B
        dS = conservative_divergence(S, D, dx) + k_dis * OM - uptake
        dOM = -k_dis * OM

        B = np.clip(B + dt * dB, 0.0, None)
        S = np.clip(S + dt * dS, 0.0, None)
        OM = np.clip(OM + dt * dOM, 0.0, None)

    return {
        "R_t": R_t,
        "cum_co2": cum_co2,
        "S_habitat_mean_t": S_habitat_mean_t,
        "B_final": B,
        "S_final": S,
        "OM_final": OM,
        "grids": grids,
        "dt": dt,
    }
