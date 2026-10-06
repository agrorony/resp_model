# Re-runs frozen Stage 12 ensemble members (same seeds/params) to save per-seed R(t) curves for the appendix figure.
import os, sys, time, numpy as np
sys.path.insert(0, os.getcwd())
import stage12_run as s
from model import load_yaml
soil, i0, i1 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
bio = load_yaml("configs/biology.yaml"); mp = load_yaml("configs/mapping.yaml")
for i in range(i0, i1):
    f = f"results/appendix_figure/curve_{soil}_{i:02d}.npy"
    if os.path.exists(f): continue
    t=time.time(); _, c = s.run_one(soil, s.SEED_BASE[soil] + i, bio, mp); np.save(f, c); print(soil, i, round(time.time()-t,1), flush=True)
