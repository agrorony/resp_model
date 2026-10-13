"""
Stage 21 runner (3D packing; Sand and Vertisol ONLY, no Loess, no measured data).

    python stage21_run.py step0        code validation (point source, extrusion, LT)
    python stage21_run.py step1        geometry only
    python stage21_run.py step2        regime map (S, V, one seed)
    python stage21_run.py step3        distinctness + controls
    python stage21_run.py step4        seed robustness
    python stage21_run.py pilot        10 um / 128^3 pilot
"""
from __future__ import annotations

import copy
import hashlib
import json
import os
import sys
import time

import numpy as np

import model3d as M
import packing3d as P
from model import load_yaml

RESULTS_DIR = os.path.join("results", "stage21")
BIOLOGY_SHA256 = "c76d9869dc5847f830ac3caadcce4c68986a9180875f3333504c2482a4b0bf5f"
OM_DENSITIES = (0.037, 0.084, 0.192, 0.439, 1.0)
THRESHOLD = 0.3


def check_biology():
    with open("configs/biology.yaml", "rb") as f:
        assert hashlib.sha256(f.read().replace(b"\r\n", b"\n")).hexdigest() == BIOLOGY_SHA256, "biology.yaml modified"


def save(name, obj):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(os.path.join(RESULTS_DIR, name), "w") as f:
        json.dump(obj, f, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))


def load(name):
    with open(os.path.join(RESULTS_DIR, name)) as f:
        return json.load(f)


# ------------------------------------------------------------------ step 0

def _point_source(xp):
    n = 33
    D = xp.ones((n, n, n))
    dt, dx, steps = 0.05, 1.0, 200
    out = {}
    for solver in ("cg", "jacobi"):
        dif = M.Diffuser(D, dx, dt, xp, solver, iterations=400, relaxation=1.0, tol=1e-15, cg_rtol=1e-14)
        u = xp.zeros((n, n, n)); u[n // 2, n // 2, n // 2] = 1.0
        for _ in range(steps):
            u = dif.solve(u, x0=u)
        ax = np.arange(n) - n // 2
        uu = P.to_np(u)
        var = [float((uu.sum(axis=tuple(a for a in range(3) if a != k)) * ax ** 2).sum() / uu.sum())
               for k in range(3)]
        out[solver] = {"mass_rel_err": float(abs(uu.sum() - 1.0)), "var_xyz": var,
                       "isotropy_rel_spread": float((max(var) - min(var)) / np.mean(var)),
                       "var_analytic_2Dt": 2 * 1.0 * dt * steps, "min": float(uu.min())}
    out["cg_vs_jacobi_max_abs"] = float(np.abs(
        P.to_np(_run_ps(xp, "cg")) - P.to_np(_run_ps(xp, "jacobi"))).max())
    return out


def _run_ps(xp, solver):
    n = 33
    D = xp.ones((n, n, n))
    dif = M.Diffuser(D, 1.0, 0.05, xp, solver, iterations=400, relaxation=1.0, tol=1e-15)
    u = xp.zeros((n, n, n)); u[n // 2, n // 2, n // 2] = 1.0
    for _ in range(60):
        u = dif.solve(u, x0=u)
    return u


def _extrusion_fields(soil, nz=4):
    """Stage 20 A1 (F1 + W0) 2D grids at n = 40, OM 0.439, with the OLD rules,
    reassembled exactly as simulate_dual_porosity does, extruded along z."""
    import dual_porosity as dp
    import stage20_run as s20
    bio = load_yaml("configs/biology.yaml")
    cfg = s20.structure_cfg(soil, 40, "W0")
    m = s20.mapping_for(0.439, "F1")
    ref = dp.simulate_dual_porosity(bio, cfg, m, wicking_enabled=False)
    g = dp.build_dual_grids(cfg, m)
    D_scale = dp._variantc_D_scale(m, cfg["pore"].get("psd_soil", ""))
    macro, wet = g["macro_mask"], g["water_mask"].copy()
    pm = g["porous_matrix_mask"]
    dry = macro & ~wet
    Dbase = g["D_full"] * D_scale
    _, _, scale = dp.habitat_properties(macro)
    kcl = float(m.get("stage4_variantC_cluster_K_scale", 0.75))
    dcl = float(m.get("stage4_variantC_cluster_D_scale", 0.75))
    f_rel = float(g["S_w"]) ** (10.0 / 3.0)
    K_t = np.where(macro, g["K"] * (1.0 + kcl * scale), 0.0)
    K_t = K_t * float(m["stage19_K_density"]) * K_t.size / K_t.sum()
    Dm = np.where(wet, Dbase * (1.0 + dcl * scale), 0.0)
    Dm = np.where(dry, Dbase * f_rel, Dm)
    kmd = float(m.get("stage4_variantC_matrix_diffusion", 0.03)) * D_scale
    Dt = Dm + np.where(wet & pm, kmd, 0.0)
    ext = lambda a: np.repeat(np.asarray(a)[:, :, None], nz, axis=2)
    F = {"D": ext(Dt), "K": ext(K_t), "B0": ext(g["B0"]), "OM_a": ext(g["OM_coating"]),
         "OM_b": ext(g["OM_trapped"]), "k_dis_b": float(m["k_dis_slow"]),
         "pore": ext(macro), "matrix": ext(pm), "hab_wet": ext(macro & wet),
         "hab_dry": ext(dry), "k_ex": float(m["k_leak"])}
    return bio, cfg, m, ref, F, g["dx"], nz


def _extrusion(soil, xp):
    bio, cfg, m, ref, F, dx, nz = _extrusion_fields(soil)
    Fx = {k: (xp.asarray(v) if isinstance(v, np.ndarray) else v) for k, v in F.items()}
    out = M.simulate3d(Fx, bio, cfg["T"], bio["dt"], dx, xp, solver="jacobi", exchange="oneway",
                       z_diffusion=False, iterations=int(m["stage6_implicit_iterations"]),
                       relaxation=float(m["stage6_implicit_relaxation"]),
                       jacobi_tol=float(m["stage11_implicit_tol"]))
    res = {}
    for name, a, b in (("R", out["R"] / nz, ref["R_t"]), ("B", out["B"] / nz, ref["B_total_t"]),
                       ("fuel", out["fuel"] / nz, ref["fuel_remaining_t"] * ref["grids"]["om_total_init"])):
        res[name] = float(np.abs(a - b).max() / np.abs(b).max())
    res["R_rel_err_maxnorm"] = res["R"]
    res["R_peak_2d"] = float(ref["R_t"].max())
    return res


def step0():
    check_biology()
    xp = P.get_xp()
    res = {"backend": xp.__name__}
    res["point_source_3d"] = {"backend_" + xp.__name__: _point_source(xp)}
    res["extrusion"] = {}
    for soil in ("sand", "vertisol"):
        res["extrusion"][soil] = {"gpu" if xp is not np else "cpu": _extrusion(soil, xp)}
        if xp is not np:
            res["extrusion"][soil]["cpu"] = _extrusion(soil, np)
    # local thickness: analytic checks
    g = np.indices((64, 64, 64))
    sph = ((g[0] - 32) ** 2 + (g[1] - 32) ** 2 + (g[2] - 32) ** 2) <= 20 ** 2
    lt = P.to_np(P.local_thickness(xp.asarray(sph), 1.0, xp))
    slab = np.zeros((20, 20, 20), bool); slab[:, :, 5:8] = True
    lts = P.to_np(P.local_thickness(xp.asarray(slab), 1.0, xp))
    res["local_thickness"] = {"sphere_d40_range": [float(lt[sph].min()), float(lt[sph].max())],
                              "slab_3voxel": float(np.unique(lts[slab])[0])}
    hp = M.load_inherited()
    d = np.random.default_rng(0).uniform(1, 600, 5000)
    from pore_field import K_window, D_of_d
    mp = load_yaml("configs/mapping.yaml")
    mp2 = dict(mp, K_d_low=10.0)
    res["KD_match_pore_field"] = {
        "K": float(np.abs(P.to_np(M.K_window3(xp.asarray(d), hp, xp, 10.0)) - K_window(d, mp2)).max()),
        "D": float(np.abs(P.to_np(M.D_of_d3(xp.asarray(d), hp, xp)) - D_of_d(d, mp)).max())}
    save("step0.json", res)
    return res


# ------------------------------------------------------------------ step 1

KINDS = ("S", "V", "Sp", "Vp", "Val")
N_FULL, VOX_FULL = 256, 5.0


def get_structure(kind, seed_offset=0, n=N_FULL, vox=VOX_FULL, xp=None):
    """Build (or load from the untracked cache) a structure; computes its LT."""
    xp = xp or P.get_xp()
    os.makedirs(os.path.join(RESULTS_DIR, "cache"), exist_ok=True)
    fn = os.path.join(RESULTS_DIR, "cache", f"{kind}_s{seed_offset}_n{n}_v{vox:g}.npz")
    if os.path.exists(fn):
        z = np.load(fn, allow_pickle=True)
        st = P.Structure()
        st.kind, st.n, st.vox_um, st.seed = kind, n, vox, int(z["seed"])
        st.phase = xp.asarray(z["phase"]); st.info = z["info"].item()
        st.lt = xp.asarray(z["lt"])
        return st
    st = P.build_structure(kind, seed_offset, n, vox, xp)
    st.lt = P.local_thickness(st.phase == P.PORE, vox, xp)
    np.savez_compressed(fn, phase=P.to_np(st.phase), lt=P.to_np(st.lt), seed=st.seed, info=np.array(st.info, dtype=object))
    return st


def geometry_gate(rep_by_kind):
    out = {}
    for k in ("S", "V"):
        r = rep_by_kind[k]
        out[k] = {"spans_any_axis": bool(any(r["spans_axes_xyz"])),
                  "habitat_gt_1pct": bool(r["habitat_volume_fraction"] > 0.01)}
    out["pass"] = all(all(v.values()) for v in out.values() if isinstance(v, dict))
    return out


def render_slices(structs, fn):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    fig, ax = plt.subplots(2, len(structs), figsize=(3.2 * len(structs), 6.6))
    cm = ListedColormap(["white", "#c8a165", "#3a3a3a"])
    for j, st in enumerate(structs):
        z = st.n // 2
        ph = P.to_np(st.phase[:, :, z]); lt = P.to_np(st.lt[:, :, z])
        ax[0, j].imshow(ph.T, cmap=cm, vmin=0, vmax=2, origin="lower")
        ax[0, j].set_title(f"{st.kind}: phase (z mid)", fontsize=9)
        im = ax[1, j].imshow(np.ma.masked_equal(lt, 0).T, cmap="viridis", origin="lower", vmin=0, vmax=250)
        ax[1, j].set_title(f"{st.kind}: local thickness (um)", fontsize=9)
        for a in ax[:, j]:
            a.set_xticks([]); a.set_yticks([])
    fig.colorbar(im, ax=ax[1, :], shrink=0.7, label="um")
    fig.savefig(fn, dpi=110); plt.close(fig)


def step1(n=N_FULL, vox=VOX_FULL, tag=""):
    xp = P.get_xp()
    hp = M.load_inherited()
    reps, structs, times = {}, [], {}
    for k in KINDS:
        t = time.time()
        st = get_structure(k, 0, n, vox, xp)
        reps[k] = P.analyse(st, xp)
        times[k] = time.time() - t
        structs.append(st)
        print(k, f"{times[k]:.1f}s", flush=True)
    kn = {k: float(M.native_K(structs[KINDS.index(k)], hp, xp).sum()) / structs[0].phase.size
          for k in ("S", "V")}
    K_density = float(np.mean(list(kn.values())))
    res = {"reports": reps, "native_K_density": kn, "K_density": K_density,
           "gate": geometry_gate(reps), "build_seconds": times}
    save(f"step1{tag}.json", res)
    render_slices(structs, os.path.join(RESULTS_DIR, f"step1_slices{tag}.png"))
    return res

# ------------------------------------------------------------------ main

if __name__ == "__main__":
    cmd = sys.argv[1]
    t0 = time.time()
    r = globals()[cmd]()
    print(json.dumps(r, indent=1, default=str)[:6000] if r is not None else "")
    print(f"[{cmd}] {time.time() - t0:.1f} s")
