"""
Loader for the real, measured pore-size distributions (Stage 3,
soil_respiration_prompt_stage3.md SS2).

Each soil's data lives in `data/psd/<name>/` (copied once from the
diagnostics pipeline output on Z:, see `source.txt` in each folder):

- `psd_table.csv` -- one row per bin: `Diameter_um` (bin center),
  `Volume_Count` (raw voxel count in that bin, i.e. volume-weighted),
  `Cumulative_Porosity` (cumulative volume fraction through that bin's right
  edge), `Differential_PSD` (volume fraction density per micron).
- `bin_edges_um.csv` -- the 51 bin edges (50 bins) in microns, exactly as
  produced by the diagnostics pipeline (log-spaced, soil-specific -- we use
  each soil's own bins, never re-binned).

Pore sizes are DIAMETERS in microns, not radii.
"""
from __future__ import annotations

import csv
import os

import numpy as np


def load_psd_table(soil_dir):
    """Read `psd_table.csv`: returns (bin_centers_um, volume_count,
    cumulative_porosity, differential_psd), each a 1D numpy array."""
    path = os.path.join(soil_dir, "psd_table.csv")
    centers, vol, cum, diff = [], [], [], []
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            centers.append(float(row["Diameter_um"]))
            vol.append(float(row["Volume_Count"]))
            cum.append(float(row["Cumulative_Porosity"]))
            diff.append(float(row["Differential_PSD"]))
    return (np.array(centers), np.array(vol), np.array(cum), np.array(diff))


def load_bin_edges(soil_dir):
    """Read `bin_edges_um.csv`: returns the n+1 bin edges in microns."""
    path = os.path.join(soil_dir, "bin_edges_um.csv")
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    return np.array([float(r[0]) for r in rows[1:]])


def empirical_cdf_at_edges(volume_count):
    """Cumulative volume fraction at each bin EDGE (length n+1, starting at
    0 at the left-most edge), built from the bin's own volume counts -- the
    inverse of this is the quantile-mapping function used to turn a uniform
    field into a diameter field (pore_field.generate_diameter_field_from_psd)."""
    total = volume_count.sum()
    fraction = volume_count / total
    return np.concatenate([[0.0], np.cumsum(fraction)])


def load_soil_psd(soil_dir):
    """Bundle everything one soil's empirical PSD needs, in one call."""
    centers, vol, cum, diff = load_psd_table(soil_dir)
    edges = load_bin_edges(soil_dir)
    cdf_edges = empirical_cdf_at_edges(vol)
    return {
        "bin_centers_um": centers,
        "bin_edges_um": edges,
        "volume_count": vol,
        "cumulative_porosity": cum,
        "differential_psd": diff,
        "cdf_edges": cdf_edges,
    }


def habitable_volume_fraction(psd, d_lo=30.0, d_hi=150.0):
    """Volume fraction of pores with diameter in [d_lo, d_hi] (the K(d)
    plateau window) -- a PSD summary stat used as the Stage-3 SS9 baseline
    to compare against the isolation ratio."""
    centers, vol = psd["bin_centers_um"], psd["volume_count"]
    mask = (centers >= d_lo) & (centers <= d_hi)
    return float(vol[mask].sum() / vol.sum())


def median_diameter_um(psd):
    """Diameter at which the empirical CDF first reaches 0.5."""
    edges, cdf = psd["bin_edges_um"], psd["cdf_edges"]
    return float(np.interp(0.5, cdf, edges))
