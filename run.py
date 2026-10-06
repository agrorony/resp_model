"""
Run one or all soils: integrates the model, saves R(t) / cumulative CO2 CSVs,
and produces the required figures (pore-size/K/D/water-filled maps,
respiration curves, cumulative CO2).

Usage:
    python run.py                # run all three soils, write everything to results/
    python run.py --soil A       # run only soil A
"""
from __future__ import annotations

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from model import load_yaml, simulate

SOILS = ["A", "B", "C"]
COLORS = {"A": "#1f77b4", "B": "#d62728", "C": "#2ca02c"}


def run_soil(soil, biology_cfg, mapping_cfg, results_dir="results/baseline"):
    structure_cfg = load_yaml(f"configs/soil_{soil}.yaml")
    out = simulate(biology_cfg, structure_cfg, mapping_cfg)

    os.makedirs(results_dir, exist_ok=True)
    dt = out["dt"]
    t = np.arange(len(out["R_t"])) * dt

    csv_path = os.path.join(results_dir, f"soil_{soil}_timeseries.csv")
    np.savetxt(
        csv_path,
        np.column_stack([t, out["R_t"], out["cum_co2"]]),
        delimiter=",",
        header="time,R_t,cumulative_co2",
        comments="",
    )
    return out, t


def plot_structure_maps(soil, out, results_dir="results"):
    grids = out["grids"]
    has_percolation = not np.all(grids["water_mask"])
    fields = [
        ("pore radius r", grids["r"]),
        ("K(r)", grids["K"]),
        ("D used", grids["D"]),
    ]
    if has_percolation:
        fields.append(("water-filled mask", grids["water_mask"].astype(float)))

    fig, axes = plt.subplots(1, len(fields), figsize=(4 * len(fields), 4))
    if len(fields) == 1:
        axes = [axes]
    for ax, (name, field) in zip(axes, fields):
        im = ax.imshow(field, cmap="viridis", origin="lower")
        ax.set_title(f"Soil {soil}: {name}")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle(f"Soil {soil} (theta={grids['theta']:.3f})")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, f"structure_maps_{soil}.png"), dpi=150)
    plt.close(fig)


def plot_curves(all_out, results_dir="results"):
    fig, ax = plt.subplots(figsize=(7, 5))
    for soil in SOILS:
        out, t = all_out[soil]
        ax.plot(t, out["R_t"], label=f"Soil {soil}", color=COLORS[soil])
    ax.set_xlabel("time")
    ax.set_ylabel("R(t) - respiration rate")
    ax.set_title("Respiration curves")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "respiration_curves.png"), dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    for soil in SOILS:
        out, t = all_out[soil]
        ax.plot(t, out["cum_co2"], label=f"Soil {soil}", color=COLORS[soil])
    ax.set_xlabel("time")
    ax.set_ylabel("cumulative CO2")
    ax.set_title("Cumulative CO2")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "cumulative_co2.png"), dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--soil", choices=SOILS, default=None,
                         help="run a single soil instead of all three")
    parser.add_argument("--results-dir", default="results/baseline")
    args = parser.parse_args()

    biology_cfg = load_yaml("configs/biology.yaml")
    mapping_cfg = load_yaml("configs/mapping.yaml")
    soils_to_run = [args.soil] if args.soil else SOILS

    all_out = {}
    for soil in soils_to_run:
        out, t = run_soil(soil, biology_cfg, mapping_cfg, args.results_dir)
        all_out[soil] = (out, t)
        plot_structure_maps(soil, out, args.results_dir)
        print(f"soil {soil}: R peak={out['R_t'].max():.4f}, "
              f"R end={out['R_t'][-1]:.4f}, cum={out['cum_co2'][-1]:.4f}, "
              f"theta={out['grids']['theta']:.3f}")

    if len(soils_to_run) == 3:
        plot_curves(all_out, args.results_dir)


if __name__ == "__main__":
    main()
