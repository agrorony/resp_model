"""
Stage 21 -- 3D simulation engine (Sand and Vertisol only).

Same biology as every earlier stage (configs/biology.yaml, unchanged):
    dB/dt = r_max f(S) B (1 - B/K) - m0 B - m_s (1 - f(S)) B
    dS/dt = div(D grad S) + k_dis OM - uptake (+ matrix<->pore exchange)
    dOM/dt = -k_dis OM,    f(S) = S/(S+Ks),  uptake = growth / Y
    R = (1 - Y) growth + m0 B

`simulate3d` is a generic engine on explicit fields. It has two modes so the
2D code can be reproduced exactly (Step 0 extrusion test) and the new physics
can run on the same code path:
    exchange = 'oneway'  matrix -> pore leak, flux k * S_matrix   (Stages 4-20)
    exchange = 'twoway'  flux k_ex * (S_matrix - S_pore) per shared face
    solver   = 'jacobi'  fixed-budget relaxed Jacobi, exactly the 2D scheme
    solver   = 'cg'      Jacobi-preconditioned conjugate gradients (SPD system
                         (I + dt/dx^2 L) u = rhs; mass conserved to the CG
                         tolerance)
Every array operation uses the array module `xp` (cupy or numpy).
"""
from __future__ import annotations

import numpy as np

import packing3d as P

# ------------------------------------------------------- constants (inherited)
D_ONLY_FLAGGED = "all values below are inherited from configs/mapping.yaml; none tuned"


def load_inherited():
    from model import load_yaml
    m = load_yaml("configs/mapping.yaml")
    return {
        "K_max": float(m["K_max"]), "K_d_high": float(m["K_d_high"]),
        "K_w_decay": float(m["K_w_decay"]),
        "D_min": float(m["D_min"]), "D_max": float(m["D_max"]),
        "D_d_ref": float(m["D_d_ref"]), "p": float(m["p"]),
        "D_scale": float(m["stage7_D_scale"]),
        "k_matrix_diff": float(m["stage4_variantC_matrix_diffusion"]),
        "k_ex": float(m["k_leak"]),               # inherited value of k_leak
        "film_K_scale": float(m["stage7_film_K_scale"]),
        "d_cut_um": float(m["matric_d_cut_um"]),
        "K_d_low_um": 10.0,                       # Stage 21: 2 voxels at 5 um
        "B0_density": 0.0125,
        "T": 1500, "voxel_ref_um": 10.0,
    }


def K_window3(d, hp, xp, d_low):
    d_hi, w = hp["K_d_high"], hp["K_w_decay"]
    return xp.where(d < d_low, 0.0,
                    xp.where(d <= d_hi, hp["K_max"], hp["K_max"] * xp.exp(-((d - d_hi) / w) ** 2)))


def D_of_d3(d, hp, xp):
    return hp["D_min"] + (hp["D_max"] - hp["D_min"]) * (d / (d + hp["D_d_ref"])) ** hp["p"]


# --------------------------------------------------------------- operators

def face_conductivities(D, xp, z_diffusion=True):
    """Harmonic-mean face conductivities along each axis (zero if either side
    is zero)."""
    out = []
    for ax in range(3):
        lo = [slice(None)] * 3; hi = [slice(None)] * 3
        lo[ax] = slice(None, -1); hi[ax] = slice(1, None)
        a, b = D[tuple(lo)], D[tuple(hi)]
        s = a + b
        f = xp.where((a > 0) & (b > 0), 2.0 * a * b / xp.where(s > 0, s, 1.0), 0.0)
        if ax == 2 and not z_diffusion:
            f = xp.zeros_like(f)
        out.append(f)
    return out


def diag_weight(faces, shape, xp):
    dw = xp.zeros(shape, dtype=faces[0].dtype)
    for ax, f in enumerate(faces):
        lo = [slice(None)] * 3; hi = [slice(None)] * 3
        lo[ax] = slice(None, -1); hi[ax] = slice(1, None)
        dw[tuple(lo)] += f
        dw[tuple(hi)] += f
    return dw


def neighbor_sum(u, faces, xp):
    out = xp.zeros_like(u)
    for ax, f in enumerate(faces):
        lo = [slice(None)] * 3; hi = [slice(None)] * 3
        lo[ax] = slice(None, -1); hi[ax] = slice(1, None)
        out[tuple(lo)] += f * u[tuple(hi)]
        out[tuple(hi)] += f * u[tuple(lo)]
    return out


class Diffuser:
    """Implicit step u - dt div(D grad u) = rhs, zero-flux boundaries."""

    def __init__(self, D, dx, dt, xp, solver="cg", z_diffusion=True,
                 iterations=90, relaxation=0.85, tol=None, cg_rtol=1e-11, cg_maxit=400):
        self.xp, self.solver = xp, solver
        self.c = dt / (dx * dx)
        self.faces = face_conductivities(D, xp, z_diffusion)
        self.dw = diag_weight(self.faces, D.shape, xp)
        self.diag = 1.0 + self.c * self.dw
        self.iterations, self.relaxation, self.tol = iterations, relaxation, tol
        self.cg_rtol, self.cg_maxit = cg_rtol, cg_maxit
        self.last_iters = 0

    def matvec(self, u):
        return self.diag * u - self.c * neighbor_sum(u, self.faces, self.xp)

    def solve(self, rhs, x0=None):
        xp = self.xp
        if self.solver == "jacobi":
            u = rhs.copy()
            for it in range(self.iterations):
                prop = (rhs + self.c * neighbor_sum(u, self.faces, xp)) / self.diag
                u_new = xp.clip(self.relaxation * prop + (1.0 - self.relaxation) * u, 0.0, None)
                if self.tol is not None:
                    delta = float(xp.abs(u_new - u).max())
                    u = u_new
                    if delta < self.tol:
                        break
                else:
                    u = u_new
            self.last_iters = it + 1
            return u
        # preconditioned CG
        x = rhs.copy() if x0 is None else x0.copy()
        r = rhs - self.matvec(x)
        bn = float(xp.sqrt((rhs * rhs).sum()))
        if bn == 0.0:
            self.last_iters = 0
            return x
        z = r / self.diag
        p = z.copy()
        rz = float((r * z).sum())
        it = 0
        while it < self.cg_maxit:
            if float(xp.sqrt((r * r).sum())) <= self.cg_rtol * bn:
                break
            Ap = self.matvec(p)
            alpha = rz / float((p * Ap).sum())
            x += alpha * p
            r -= alpha * Ap
            z = r / self.diag
            rz_new = float((r * z).sum())
            p = z + (rz_new / rz) * p
            rz = rz_new
            it += 1
        self.last_iters = it
        return x


class Exchange:
    """Matrix <-> pore exchange across shared faces."""

    def __init__(self, matrix, pore, k, mode, xp):
        self.xp, self.k, self.mode = xp, k, mode
        self.pairs = []
        for ax in range(3):
            lo = [slice(None)] * 3; hi = [slice(None)] * 3
            lo[ax] = slice(None, -1); hi[ax] = slice(1, None)
            lo, hi = tuple(lo), tuple(hi)
            m_lo_p_hi = (matrix[lo] & pore[hi])
            p_lo_m_hi = (pore[lo] & matrix[hi])
            self.pairs.append((lo, hi, m_lo_p_hi.astype(xp.float64), p_lo_m_hi.astype(xp.float64)))

    def flux(self, S):
        """Net per-voxel gain rate (per unit time) from exchange."""
        xp, k = self.xp, self.k
        ex = xp.zeros_like(S)
        for lo, hi, a, b in self.pairs:
            Slo, Shi = S[lo], S[hi]
            if self.mode == "oneway":
                f_a = k * Slo * a           # matrix(lo) -> pore(hi)
                f_b = k * Shi * b           # matrix(hi) -> pore(lo)
            else:
                f_a = k * (Slo - Shi) * a
                f_b = k * (Shi - Slo) * b
            ex[lo] += f_b - f_a
            ex[hi] += f_a - f_b
        return ex


# -------------------------------------------------------------- the engine

STAT_NAMES = ["R", "B", "growth_R", "maint_R", "fuel", "diff_wet", "diff_dry",
              "ex_wet", "ex_dry", "S_total", "clip_added", "uptake_dt", "source_dt", "iters"]


def simulate3d(F, bio, T, dt, dx, xp, solver="cg", exchange="twoway", z_diffusion=True,
               iterations=90, relaxation=0.85, jacobi_tol=None, cg_rtol=1e-11):
    """F: dict of fields (xp arrays, shape (nx, ny, nz)):
         D        per-voxel diffusivity (pore + matrix, already scaled)
         K        habitat capacity (0 = not habitat)
         B0       initial biomass
         OM_a     OM pool dissolving at k_dis
         OM_b     optional second pool dissolving at F['k_dis_b'] (legacy)
         pore, matrix  boolean masks (exchange partners)
         hab_wet, hab_dry  boolean masks (observational split)
         k_ex     exchange coefficient
       Returns dict of per-step arrays (numpy) and final fields (xp)."""
    r_max, Ks, Y, m0, m_s, k_dis = (bio[k] for k in ("r_max", "Ks", "Y", "m0", "m_s", "k_dis"))
    K, B = F["K"], F["B0"].copy()
    S = xp.zeros_like(B)
    OM_a = F["OM_a"].copy()
    OM_b = F["OM_b"].copy() if "OM_b" in F else None
    k_b = F.get("k_dis_b", k_dis)
    habitable = K > 0
    K_safe = xp.where(habitable, K, 1.0)
    dif = Diffuser(F["D"], dx, dt, xp, solver, z_diffusion, iterations, relaxation,
                   jacobi_tol, cg_rtol)
    exch = Exchange(F["matrix"], F["pore"], F["k_ex"], exchange, xp)
    assert dt * 6.0 * F["k_ex"] < 1.0, "explicit exchange stability"
    hw, hd = F["hab_wet"].astype(xp.float64), F["hab_dry"].astype(xp.float64)
    stats = xp.zeros((T, len(STAT_NAMES)))
    for t in range(T):
        f_S = S / (S + Ks)
        logistic = xp.where(habitable, 1.0 - B / K_safe, 0.0)
        growth = r_max * f_S * B * logistic
        uptake = growth / Y
        gR = ((1.0 - Y) * growth).sum()
        mR = (m0 * B).sum()
        dB = growth - m0 * B - m_s * (1.0 - f_S) * B
        src = k_dis * OM_a + (k_b * OM_b if OM_b is not None else 0.0)
        pre = S + dt * src
        S_dif = dif.solve(pre, x0=S)
        ex = dt * exch.flux(S_dif)
        S_new = S_dif - dt * uptake + ex
        S_cl = xp.clip(S_new, 0.0, None)
        row = [gR + mR, B.sum(), gR, mR, OM_a.sum() + (OM_b.sum() if OM_b is not None else 0.0),
               ((S_dif - pre) * hw).sum(), ((S_dif - pre) * hd).sum(), (ex * hw).sum(), (ex * hd).sum(),
               S_cl.sum(), (S_cl - S_new).sum(), dt * uptake.sum(), dt * src.sum(),
               float(dif.last_iters)]
        for j, v in enumerate(row):
            stats[t, j] = v
        B = xp.clip(B + dt * dB, 0.0, None)
        S = S_cl
        OM_a = xp.clip(OM_a - dt * k_dis * OM_a, 0.0, None)
        if OM_b is not None:
            OM_b = xp.clip(OM_b - dt * k_b * OM_b, 0.0, None)
    st = P.to_np(stats)
    out = {n: st[:, i] for i, n in enumerate(STAT_NAMES)}
    out["S_final"], out["B_final"] = S, B
    out["B0_sum"] = float(F["B0"].sum())
    # mass balance: S_total(T-1) = sum(source) - sum(uptake) + sum(clip added)
    bal = out["S_total"][-1] - (out["source_dt"].sum() - out["uptake_dt"].sum() + out["clip_added"].sum())
    out["mass_balance_rel"] = float(abs(bal) / max(out["source_dt"].sum(), 1e-300))
    out["clip_added_rel"] = float(out["clip_added"].sum() / max(out["source_dt"].sum(), 1e-300))
    return out


def summarize(out, bio):
    """Liveness / shape metrics (Stage 19 definitions)."""
    R = out["R"]
    pk = int(np.argmax(R))
    floor = bio["m0"] * out["B0_sum"]
    gs = out["growth_R"][pk] / max(out["growth_R"][pk] + out["maint_R"][pk], 1e-300)
    d = {
        "B_peak_over_B0": float(out["B"].max() / out["B0_sum"]),
        "R_peak_over_floor": float(R.max() / floor),
        "growth_share_at_peak": float(gs),
        "peak_frac": pk / (len(R) - 1),
        "R_end_over_peak": float(R[-1] / R.max()),
        "fuel_remaining_final": float(out["fuel"][-1] / out["fuel"][0]) if out["fuel"][0] > 0 else 0.0,
        "diffusion_into_wet_habitat": float(out["diff_wet"].sum()),
        "diffusion_into_drained_habitat": float(out["diff_dry"].sum()),
        "exchange_into_wet_habitat": float(out["ex_wet"].sum()),
        "exchange_into_drained_habitat": float(out["ex_dry"].sum()),
        "mass_balance_rel": out["mass_balance_rel"],
        "clip_added_rel": out["clip_added_rel"],
        "mean_cg_iters": float(out["iters"].mean()),
    }
    d["alive_primary"] = bool(d["B_peak_over_B0"] > 1.0 and d["R_peak_over_floor"] >= 2.0
                              and d["growth_share_at_peak"] > 0.5)
    d["alive_secondary"] = bool(d["R_peak_over_floor"] >= 2.0 and d["growth_share_at_peak"] > 0.5)
    return d


# ------------------------------------------- Stage 21 physics: field building

def native_K(st, hp, xp, lt=None):
    """Unscaled habitat capacity of a structure (film scale on drained)."""
    phase = st.phase
    pore = phase == P.PORE
    lt = st.lt if lt is None else lt
    wet = pore & (lt < hp["d_cut_um"])
    K = xp.where(pore, K_window3(lt, hp, xp, hp["K_d_low_um"]), 0.0)
    K = xp.where(pore & ~wet, K * hp["film_K_scale"], K)
    return K


def build_fields(st, hp, xp, K_density, OM_density, D_scale=None):
    """Stage 21 rules (see STAGE21_RESULTS.md). st.lt must exist (analyse())."""
    phase, lt = st.phase, st.lt
    n = phase.shape[0]
    N = float(phase.size)
    pore = phase == P.PORE
    matrix = phase == P.MATRIX
    solid = phase == P.SOLID
    wet_pore = pore & (lt < hp["d_cut_um"])
    drained = pore & ~wet_pore
    S_w = float(wet_pore.sum() + matrix.sum()) / float((~solid).sum())
    Ds = hp["D_scale"] if D_scale is None else D_scale
    Dfull = D_of_d3(lt, hp, xp) * Ds
    D = xp.where(wet_pore, Dfull, 0.0)
    D = xp.where(drained, Dfull * S_w ** (10.0 / 3.0), D)
    D = xp.where(matrix, hp["k_matrix_diff"] * Ds, D)
    Kn = native_K(st, hp, xp)
    K = Kn * (K_density * N / float(Kn.sum()))
    B0 = K * (hp["B0_density"] * N / float(K.sum()))
    OM = xp.where(matrix, OM_density * N / float(matrix.sum()), 0.0)
    habitable = K > 0
    F = {"D": D, "K": K, "B0": B0, "OM_a": OM, "pore": pore, "matrix": matrix,
         "hab_wet": habitable & wet_pore, "hab_dry": habitable & drained, "k_ex": hp["k_ex"]}
    meta = {"S_w": S_w, "N": N, "K_native_density": float(Kn.sum()) / N,
            "n_habitat_voxels": int(habitable.sum()),
            "matrix_voxels": int(matrix.sum())}
    return F, meta


def run_structure(st, hp, bio, xp, K_density, OM_density, T=None, D_scale=None, solver="cg"):
    F, meta = build_fields(st, hp, xp, K_density, OM_density, D_scale)
    N = meta["N"]
    # equal totals as per-voxel densities, asserted every run
    assert abs(float(F["OM_a"].sum()) / (OM_density * N) - 1) < 1e-9, "OM total"
    assert abs(float(F["B0"].sum()) / (hp["B0_density"] * N) - 1) < 1e-9, "B0 total"
    assert abs(float(F["K"].sum()) / (K_density * N) - 1) < 1e-9, "K total"
    dx = st.vox_um / hp["voxel_ref_um"]
    out = simulate3d(F, bio, T or hp["T"], bio["dt"], dx, xp, solver=solver, exchange="twoway")
    res = summarize(out, bio)
    res.update(meta)
    return out, res
