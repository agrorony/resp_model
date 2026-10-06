"""
Loader for real MEASURED pore-size distributions + Track E connectivity data
(Stage 13, soil_respiration_prompt_stage13_real_data.md), replacing
psd_parametric's literature-mixture PSDs for the `pore.mode: psd_measured`
config path. Source: the companion "resarch exercise" project's CT-scan +
nnU-Net segmentation + porespy topology-metrics pipeline -- see
`configs/psd_measured.yaml` and each soil's `data/psd_measured/<soil>/
source.txt` for full provenance.

The real pipeline's psd_table.csv exports (unlike the Stage-3 `data/psd/
<soil>/` copies psd_data.py reads) ship only bin CENTERS (`Diameter_um`),
not a separate bin-edges file -- `_edges_from_log_centers` reconstructs the
n+1 edges from the n log-spaced centers by inverting psd_parametric.
build_mixture_psd's own edge->center convention (center = sqrt(edge_lo *
edge_hi)), recovering the true edges essentially exactly for a genuinely
log-spaced grid (every soil's real table here is log-spaced end to end).

`build_soil_psd`'s return dict has the SAME shape as psd_data.load_soil_psd
and psd_parametric.build_soil_psd, so build_dual_grids and everything
downstream (pore_field.generate_diameter_field_from_psd's empirical-
quantile mapping, psd_data.habitable_volume_fraction/median_diameter_um)
work unchanged regardless of which of the three PSD sources fed it.
"""
from __future__ import annotations

import csv
import os

import numpy as np
import yaml

import psd_data


def _edges_from_log_centers(centers):
    """Invert psd_parametric's center = sqrt(edge_lo*edge_hi) convention to
    recover n+1 bin edges from n log-spaced bin centers. Interior edges are
    the geometric mean of adjacent centers; the two outer edges extrapolate
    the local log-spacing of the nearest pair."""
    centers = np.asarray(centers, dtype=float)
    log_c = np.log(centers)
    interior = np.exp(0.5 * (log_c[:-1] + log_c[1:]))
    step_first = log_c[1] - log_c[0]
    step_last = log_c[-1] - log_c[-2]
    first_edge = np.exp(log_c[0] - step_first / 2.0)
    last_edge = np.exp(log_c[-1] + step_last / 2.0)
    return np.concatenate([[first_edge], interior, [last_edge]])


def load_measured_psd_table(csv_path):
    """Read a real pipeline's psd_table.csv: returns (bin_centers_um,
    volume_count, cumulative_porosity, differential_psd). Extra columns
    (Euler_Number, Connectivity_*, Tortuosity_*, etc., repeated on every
    row in the real export) are ignored here -- their values are instead
    curated, with provenance, in configs/psd_measured.yaml."""
    centers, vol, cum, diff = [], [], [], []
    with open(csv_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            centers.append(float(row["Diameter_um"]))
            vol.append(float(row["Volume_Count"]))
            cum.append(float(row["Cumulative_Porosity"]))
            diff.append(float(row["Differential_PSD"]))
    return (np.array(centers), np.array(vol), np.array(cum), np.array(diff))


def load_soil_psd_measured(csv_path, bin_edges_csv=None):
    """Bundle one soil's real measured PSD, same dict shape as
    psd_data.load_soil_psd. If `bin_edges_csv` is given (Rehovot's Stage-3
    copy ships its own, exact edges file), it is used directly instead of
    the log-center reconstruction."""
    centers, vol, cum, diff = load_measured_psd_table(csv_path)
    if bin_edges_csv and os.path.exists(bin_edges_csv):
        edges = psd_data.load_bin_edges(os.path.dirname(bin_edges_csv))
    else:
        edges = _edges_from_log_centers(centers)
    cdf_edges = psd_data.empirical_cdf_at_edges(vol)
    return {
        "bin_centers_um": centers,
        "bin_edges_um": edges,
        "volume_count": vol,
        "cumulative_porosity": cum,
        "differential_psd": diff,
        "cdf_edges": cdf_edges,
    }


def load_measured_config(path="configs/psd_measured.yaml"):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def build_soil_psd(soil_name, config_path="configs/psd_measured.yaml"):
    """Build one soil's real measured PSD from configs/psd_measured.yaml.
    Attaches `porosity`/`connectivity` (the curated real Track E scalars)
    onto the returned dict for Stage 13's reporting -- dual_porosity.py
    never reads these two keys (mirrors psd_parametric.build_soil_psd's
    unused `total_porosity`), they exist purely for stage13_run.py."""
    cfg = load_measured_config(config_path)
    if soil_name not in cfg["soils"]:
        raise KeyError(f"no measured PSD for soil {soil_name!r}; "
                       f"have {sorted(cfg['soils'])}")
    soil = cfg["soils"][soil_name]
    csv_dir = os.path.dirname(soil["psd_csv"])
    bin_edges_csv = os.path.join(csv_dir, "bin_edges_um.csv")
    psd = load_soil_psd_measured(soil["psd_csv"], bin_edges_csv)
    psd["total_porosity"] = soil.get("porosity")
    psd["connectivity"] = soil.get("connectivity")
    psd["name"] = soil_name
    return psd
