"""
Stage 18 single-compartment model: uniform OM per cell, real measured PSDs
(docs/stage_results/STAGE18_RESULTS.md).

Rules (fixed, shared by every soil, no per-soil knob):
  1. OM(i,j) = OM_total / N for every cell. No pore-size weighting. Because
     each cell's diameter is drawn from the soil's real measured PSD, the
     number of cells at each diameter already follows that soil's real
     pore-volume distribution, so uniform-per-cell OM IS pore-volume-
     proportional OM, with no invented weighting.
  2. D(d) is the only pore-size-dependent transport law: smooth, continuous,
     >= D_min > 0 (water gating by the unchanged saturation physics aside).
  3. K(d) is a smooth lognormal hump, floored at k_min (never zero, no cutoff).
  4. Equal totals enforced and asserted: sum(OM) = OM_total, sum(B0) =
     B0_total, sum(K) = K_mean * N.
  5. Real measured PSDs only (pore.mode: psd_measured).
  6. Saturation physics unchanged (model.water_mask_and_theta; the measured
     configs use global_theta, theta = 0.6).
  7. No enzyme, no macro/matrix compartments, no categorical pore classes.

Biology (configs/biology.yaml) and the per-cell equations are main's
`model.simulate` unchanged:
    dB/dt  = r_max f(S) B (1 - B/K) - m0 B - m_s (1 - f(S)) B
    dS/dt  = div(D grad S) + k_dis OM - (1/Y) r_max f(S) B (1 - B/K)
    dOM/dt = -k_dis OM
    R      = sum[(1 - Y) growth + m0 B]
numpy-only; explicit conservative-flux transport (model.conservative_divergence).
"""
from __future__ import annotations

import numpy as np

import psd_measured
from model import conservative_divergence, load_yaml, water_mask_and_theta
from pore_field import b0_field, generate_diameter_field_from_psd


def D_of_d(d, m):
    return m["D_min"] + (m["D_max"] - m["D_min"]) * (d / (d + m["D_d_ref"])) ** m["p"]


def K_shape_of_d(d, m):
    hump = np.exp(-0.5 * (np.log(d / m["K_d_opt"]) / m["w_K"]) ** 2)
    return m["k_min"] + (m["K_max"] - m["k_min"]) * hump


def diameter_field(structure_cfg, shuffle_seed=None):
    """Real measured PSD laid out with main's shared correlated-Gaussian
    arrangement recipe. `shuffle_seed` permutes the field: same diameter
    multiset, arrangement destroyed (robustness control)."""
    pore = structure_cfg["pore"]
    assert pore["mode"] == "psd_measured", "Stage 18 uses real measured PSDs only"
    psd = psd_measured.build_soil_psd(pore["psd_soil"], pore.get("psd_config", "configs/psd_measured.yaml"))
    d = generate_diameter_field_from_psd(structure_cfg["grid"]["n"], pore["lambda"], pore["seed"],
                                         psd["bin_edges_um"], psd["cdf_edges"], pore.get("aggregate", False))
    if shuffle_seed is not None:
        d = np.random.default_rng(int(shuffle_seed)).permutation(d.ravel()).reshape(d.shape)
    return d


def build_grids(structure_cfg, m, shuffle_seed=None, transport=True):
    n = structure_cfg["grid"]["n"]
    N = n * n
    d = diameter_field(structure_cfg, shuffle_seed)

    K_shape = K_shape_of_d(d, m)
    K = K_shape * (m["K_mean"] * N / K_shape.sum())          # rule 4: sum(K) = K_mean*N
    OM = np.full((n, n), m["OM_total"] / N)                   # rule 1: uniform per cell
    B0 = b0_field(K, m["B0_total"])

    sat_cfg = structure_cfg.get("saturation", {"mode": "fully_wet"})
    water_mask, theta = water_mask_and_theta(d, sat_cfg, m)  # rule 6: unchanged
    D = np.where(water_mask, D_of_d(d, m), 0.0)
    if not transport:
        # Diagnostic control only (stage18_run.py): every face conductance
        # zero, i.e. no substrate moves between cells. Not a model option.
        D = np.zeros_like(D)

    for name, arr, tot in (("OM", OM, m["OM_total"]), ("B0", B0, m["B0_total"]),
                           ("K", K, m["K_mean"] * N)):
        assert abs(arr.sum() / tot - 1.0) < 1e-9, f"equal-totals rule violated for {name}"
    assert K.min() > 0, "K must never be zero"

    return {"n": n, "dx": structure_cfg["grid"]["dx"], "d": d, "K": K, "D": D,
            "water_mask": water_mask, "theta": theta, "OM": OM, "B0": B0, "S0": np.zeros((n, n))}


def simulate(biology_cfg, structure_cfg, m, shuffle_seed=None, transport=True):
    g = build_grids(structure_cfg, m, shuffle_seed, transport)
    dx, D, K = g["dx"], g["D"], g["K"]
    r_max, Ks, Y = biology_cfg["r_max"], biology_cfg["Ks"], biology_cfg["Y"]
    m0, m_s, k_dis, dt = biology_cfg["m0"], biology_cfg["m_s"], biology_cfg["k_dis"], biology_cfg["dt"]
    T = int(structure_cfg["T"])
    assert D.max() * dt / dx ** 2 < 0.2, "diffusion stability violated"
    S_scale = float(m["OM_total"])  # transport-error normalization (total substrate source)

    B, S, OM = g["B0"].copy(), g["S0"].copy(), g["OM"].copy()
    R_t, B_t = np.zeros(T), np.zeros(T)
    transport_mass_err = 0.0
    for t in range(T):
        f_S = S / (S + Ks)
        growth = r_max * f_S * B * (1.0 - B / K)
        R_t[t] = ((1.0 - Y) * growth + m0 * B).sum()
        B_t[t] = B.sum()
        div = conservative_divergence(S, D, dx)
        transport_mass_err = max(transport_mass_err, abs(float(div.sum())))
        dB = growth - m0 * B - m_s * (1.0 - f_S) * B
        dS = div + k_dis * OM - growth / Y
        B = np.clip(B + dt * dB, 0.0, None)
        S = np.clip(S + dt * dS, 0.0, None)
        OM = np.clip(OM - dt * k_dis * OM, 0.0, None)
    out = {"R_t": R_t, "B_t": B_t, "cum_co2": np.cumsum(R_t) * dt, "B_final": B, "S_final": S,
           "OM_final": OM, "grids": g, "dt": dt, "transport_mass_err_max": transport_mass_err,
           "transport_mass_err_rel": transport_mass_err / S_scale}
    return out


def is_nontrivial(out):
    """main's metrics.is_nontrivial criterion: biomass grows, respiration > 0."""
    return bool(out["B_final"].sum() > out["grids"]["B0"].sum() and out["R_t"].sum() > 0)
