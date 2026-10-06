# Re-runs the frozen Stage 12 ensembles (same seeds, same params) only to save the per-seed R(t) curves for the appendix figure.
import os, sys, numpy as np
sys.path.insert(0, os.getcwd())
import stage12_run as s
from model import load_yaml
bio = load_yaml("configs/biology.yaml"); mp = load_yaml("configs/mapping.yaml")
for soil in s.SOILS:
    seeds = [s.SEED_BASE[soil] + i for i in range(s.M_SEEDS)]
    rows, curves = s.run_ensemble(soil, seeds, bio, mp)
    np.save(f"results/appendix_figure/curves_{soil}.npy", curves)
    print(soil, "done", flush=True)
