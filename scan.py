"""
Candidate-library scan + triple selection + S3/S4 checks (v3 search, MODEL_SPEC.md SS13).

One scan round = simulate a library of candidate pore fields under the
shared rules of a path (grid, T, saturation mode, mechanics), drop every
candidate that violates R3 (OM hoarding) or is trivial, then choose the
three candidates whose minimum pairwise distance is largest. The winning
triple is then checked for:

  S3  seed robustness  -- the same three soil definitions with 3 seed sets
  S4  arrangement      -- each soil's pore field spatially shuffled (same
                          pore-size multiset, no arrangement): how much
                          distance is lost?

Nothing here touches biology.yaml or mapping.yaml; the path settings are
passed in explicitly and logged.
"""
from __future__ import annotations

import copy
import itertools
import json
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import metrics
import model
from pore_field import shuffled

BIO = "configs/biology.yaml"
MAP = "configs/mapping.yaml"


def soil_cfg(pore, shape, T, sat_mode, seed):
    pore = dict(pore)
    pore["seed"] = seed
    return {"grid": {"shape": list(shape), "dx": 1.0}, "T": T, "pore": pore,
            "saturation": {"mode": sat_mode}}


def _mapping(enzyme):
    mp = model.load_yaml(MAP)
    mp["enzyme"]["enabled"] = bool(enzyme)
    return mp


def _run(args):
    cfg, enzyme, shuffle = args
    bio = model.load_yaml(BIO)
    mp = _mapping(enzyme)
    r = None
    if shuffle:
        r = shuffled(model.build_grids(cfg, mp)["r"])
    out = model.simulate(bio, cfg, mp, r=r)
    g = out["grids"]
    return {"R_t": out["R_t"], "nontrivial": metrics.is_nontrivial(out),
            "r3": model.r3_ok(g), "hoard": g["om_hoard"], "theta": g["theta"],
            "cum": float(out["cum_co2"][-1])}


def run_many(cfgs, enzyme=False, shuffle=False, workers=4):
    with ProcessPoolExecutor(workers) as ex:
        return list(ex.map(_run, [(c, enzyme, shuffle) for c in cfgs]))


def dist(a, b):
    return metrics.pairwise_distance(a, b)


def best_triples(curves, k=5):
    """Top-k triples by min pairwise distance (exhaustive)."""
    names = list(curves)
    D = {}
    for a, b in itertools.combinations(names, 2):
        D[(a, b)] = D[(b, a)] = dist(curves[a], curves[b])
    scored = []
    for tri in itertools.combinations(names, 3):
        m = min(D[p] for p in itertools.combinations(tri, 2))
        scored.append((m, tri))
    scored.sort(reverse=True)
    return scored[:k], D


def library(mus, sigmas, structures):
    """structures: list of (label, dict of pore keys)."""
    lib = {}
    for mu in mus:
        for sg in sigmas:
            for label, extra in structures:
                pore = {"mu": mu, "sigma": sg}
                pore.update(extra)
                lib[f"mu{mu}_s{sg}_{label}"] = pore
    return lib


def scan(lib, shape, T, sat_mode, enzyme=False, seed=1, workers=4):
    names = list(lib)
    res = run_many([soil_cfg(lib[n], shape, T, sat_mode, seed) for n in names],
                   enzyme=enzyme, workers=workers)
    res = dict(zip(names, res))
    ok = {n: v for n, v in res.items() if v["r3"] and v["nontrivial"]}
    return res, ok


def check_triple(tri_pores, shape, T, sat_mode, enzyme=False,
                 seed_sets=((1, 2, 3), (11, 12, 13), (21, 22, 23)), workers=4):
    """S2/S3/S4 for one triple of pore dicts (A, B, C)."""
    cfgs, keys = [], []
    for si, seeds in enumerate(seed_sets):
        for soil, pore, seed in zip("ABC", tri_pores, seeds):
            cfgs.append(soil_cfg(pore, shape, T, sat_mode, seed))
            keys.append((si, soil))
    res = dict(zip(keys, run_many(cfgs, enzyme=enzyme, workers=workers)))
    shuf = dict(zip("ABC", run_many(cfgs[:3], enzyme=enzyme, shuffle=True, workers=workers)))

    per_set = []
    for si in range(len(seed_sets)):
        d = {a + b: dist(res[(si, a)]["R_t"], res[(si, b)]["R_t"])
             for a, b in itertools.combinations("ABC", 2)}
        valid = all(res[(si, s)]["nontrivial"] and res[(si, s)]["r3"] for s in "ABC")
        per_set.append({"distances": d, "valid": valid,
                        "pass": valid and min(d.values()) > 0.3})
    within = [dist(res[(i, s)]["R_t"], res[(j, s)]["R_t"])
              for s in "ABC" for i, j in itertools.combinations(range(len(seed_sets)), 2)]
    min_between = min(min(p["distances"].values()) for p in per_set)
    d_shuf = {a + b: dist(shuf[a]["R_t"], shuf[b]["R_t"]) for a, b in itertools.combinations("ABC", 2)}
    drop = {k: per_set[0]["distances"][k] - d_shuf[k] for k in d_shuf}
    return {
        "S2": per_set[0]["pass"],
        "S3a_all_seed_sets_pass": all(p["pass"] for p in per_set),
        "S3b_within_mean": float(np.mean(within)),
        "S3b_min_between": float(min_between),
        "S3b": float(np.mean(within)) < 0.5 * min_between,
        "S4_shuffled_distances": d_shuf,
        "S4_drop": drop,
        "S4": max(drop.values()) >= 0.1,
        "per_seed_set": per_set,
        "theta": {s: res[(0, s)]["theta"] for s in "ABC"},
        "hoard": {s: res[(0, s)]["hoard"] for s in "ABC"},
        "curves": {s: res[(0, s)]["R_t"] for s in "ABC"},
    }


def fmt(d):
    return json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in d.items()})
