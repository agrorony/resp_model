from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from typing import Dict, List, Tuple

import numpy as np


Array = np.ndarray


@dataclass(frozen=True)
class Formulation:
    transport_kind: str = "local_implicit"
    d_map: str = "linear"
    n: int = 32
    l: float = 1.0
    t_final: float = 12.0
    steps: int = 120
    corr_length_cells: float = 3.2
    threshold_a: float = 0.34
    threshold_b: float = 0.56
    reservoir_cut: float = 0.18
    d_scale: float = 0.015
    d_floor: float = 0.0
    d_cut: float = 0.18
    d_center: float = 0.52
    d_width: float = 0.28
    d_power: float = 1.3
    kappa_scale: float = 1.15
    w0_scale: float = 0.85
    w0_power: float = 1.1
    down_coating_fraction: float = 0.78
    down_coating_thickness: int = 2
    lambda_decay: float = 0.14
    a0: float = 0.1
    mu: float = 0.22
    v0: float = 0.78
    implicit_iterations: int = 90
    implicit_relaxation: float = 0.85
    nonlocal_rate: float = 2.8
    active_d_threshold: float = 1e-8


@dataclass
class RunResult:
    times: Array
    observable: Array
    s: Array
    d: Array
    kappa: Array
    w0: Array
    w0_coating: Array
    w0_trapped: Array
    separation: float
    reach: float
    transport_drift: float
    reservoir_mass: float
    reaction_mass: float


def formulation_to_dict(formulation: Formulation) -> Dict[str, float | int | str]:
    return asdict(formulation)


def with_updates(formulation: Formulation, **updates: float | int | str) -> Formulation:
    return replace(formulation, **updates)


def phi(u: Array, a0: float) -> Array:
    return u / (u + a0 + 1e-12)


def normalize01(values: Array) -> Array:
    lower = float(values.min())
    upper = float(values.max())
    if upper - lower < 1e-12:
        return np.zeros_like(values)
    return (values - lower) / (upper - lower)


def filtered_noise(rng: np.random.Generator, n: int, sigma_cells: float) -> Array:
    raw = rng.standard_normal((n, n))
    fy = np.fft.fftfreq(n)[:, None]
    fx = np.fft.fftfreq(n)[None, :]
    kernel = np.exp(-2.0 * (np.pi**2) * sigma_cells**2 * (fx * fx + fy * fy))
    smooth = np.fft.ifft2(np.fft.fft2(raw) * kernel).real
    smooth -= smooth.mean()
    std = float(smooth.std())
    return smooth / (std + 1e-12)


def make_structure_field(kind: str, seed: int, formulation: Formulation) -> Array:
    rng = np.random.default_rng(seed)
    n = formulation.n
    corr = formulation.corr_length_cells
    x = np.linspace(-1.0, 1.0, n)[None, :]
    y = np.linspace(-1.0, 1.0, n)[:, None]
    if kind == "up":
        noise_a = filtered_noise(rng, n, corr)
        noise_b = filtered_noise(rng, n, max(1.2, 0.55 * corr))
        radial = np.sqrt(x * x + y * y)
        field = 0.18 + 0.46 * radial + 0.14 * noise_a + 0.05 * noise_b**3
        return np.clip(field, 0.0, 1.0)
    if kind == "down":
        shifted = np.sqrt((x + 0.45) ** 2 + (1.15 * y) ** 2)
        noise = 0.025 * filtered_noise(rng, n, max(1.0, 0.75 * corr))
        field = 0.12 + 0.68 * shifted + noise
        return np.clip(field, 0.0, 1.0)
    raise ValueError(f"Unknown structure kind: {kind}")


def normalize_weights(weights: Array) -> Array:
    total = float(weights.sum())
    if total <= 1e-12:
        return np.zeros_like(weights)
    return weights / total


def neighbor_count(mask: Array) -> Array:
    out = np.zeros_like(mask, dtype=int)
    out[1:, :] += mask[:-1, :]
    out[:-1, :] += mask[1:, :]
    out[:, 1:] += mask[:, :-1]
    out[:, :-1] += mask[:, 1:]
    return out


def bfs_distance(seed_mask: Array, allowed_mask: Array | None = None) -> Array:
    n = seed_mask.shape[0]
    dist = np.full((n, n), np.inf, dtype=float)
    frontier = np.argwhere(seed_mask)
    if len(frontier) == 0:
        return dist
    for i, j in frontier:
        dist[i, j] = 0.0
    queue: List[Tuple[int, int]] = [(int(i), int(j)) for i, j in frontier]
    head = 0
    while head < len(queue):
        ci, cj = queue[head]
        head += 1
        base = dist[ci, cj] + 1.0
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ni = ci + di
            nj = cj + dj
            if not (0 <= ni < n and 0 <= nj < n):
                continue
            if allowed_mask is not None and not allowed_mask[ni, nj]:
                continue
            if base < dist[ni, nj]:
                dist[ni, nj] = base
                queue.append((ni, nj))
    return dist


def build_up_source(s: Array, d: Array, reaction_mask: Array, formulation: Formulation) -> Tuple[Array, Array, Array]:
    active = d > max(formulation.active_d_threshold, 1e-12)
    reaction_neighbors = neighbor_count(reaction_mask) > 0
    reaction_access = active & reaction_neighbors
    reachable = np.isfinite(bfs_distance(reaction_access, active))
    geodesic = bfs_distance(reaction_access, reachable)
    source_candidates = reachable & (~reaction_mask) & (geodesic >= 2.0)
    distance_weight = np.where(np.isfinite(geodesic), geodesic, 0.0)
    structure_weight = 0.35 + 0.65 * np.clip(1.0 - s, 0.0, 1.0) ** formulation.w0_power
    source_weight = source_candidates.astype(float) * (distance_weight**2) * structure_weight
    if source_weight.sum() <= 1e-12:
        source_weight = source_candidates.astype(float)
    w0 = formulation.w0_scale * normalize_weights(source_weight)
    return w0, np.zeros_like(w0), w0.copy()


def build_down_source(s: Array, d: Array, reaction_mask: Array, formulation: Formulation) -> Tuple[Array, Array, Array]:
    active = d > max(formulation.active_d_threshold, 1e-12)
    reaction_neighbors = neighbor_count(reaction_mask) > 0
    reaction_access = active & reaction_neighbors
    reachable = np.isfinite(bfs_distance(reaction_access, active))
    geodesic = bfs_distance(reaction_access, reachable)
    coating_candidates = reachable & (~reaction_mask) & (geodesic <= float(formulation.down_coating_thickness))
    if not coating_candidates.any():
        coating_candidates = reachable & (~reaction_mask) & reaction_neighbors
    coating_weight = np.clip(float(formulation.down_coating_thickness) + 1.0 - geodesic, 0.0, None)
    coating_weight *= coating_candidates.astype(float)
    if coating_weight.sum() <= 1e-12:
        coating_weight = coating_candidates.astype(float)

    low_d = d <= max(formulation.active_d_threshold, 1e-12)
    interior = neighbor_count(low_d) >= 2
    trapped_candidates = low_d & (~reaction_mask) & interior
    if not trapped_candidates.any():
        far_cells = (~reaction_mask) & (geodesic > float(formulation.down_coating_thickness) + 1.0)
        far_cells &= np.isfinite(geodesic)
        if not far_cells.any():
            far_cells = (~reachable) & (~reaction_mask)
        trapped_candidates = far_cells
    trapped_weight = trapped_candidates.astype(float) * np.clip(formulation.threshold_a - s, 0.0, None) ** formulation.w0_power
    if trapped_weight.sum() <= 1e-12:
        finite_geodesic = np.where(np.isfinite(geodesic), geodesic, 0.0)
        trapped_weight = trapped_candidates.astype(float) * (1.0 + finite_geodesic)
    if trapped_weight.sum() <= 1e-12:
        trapped_weight = trapped_candidates.astype(float)

    coating_mass = formulation.w0_scale * formulation.down_coating_fraction
    trapped_mass = formulation.w0_scale * (1.0 - formulation.down_coating_fraction)
    coating = coating_mass * normalize_weights(coating_weight)
    trapped = trapped_mass * normalize_weights(trapped_weight)
    w0 = coating + trapped
    return w0, coating, trapped


def build_coefficients(s: Array, formulation: Formulation, instance_kind: str) -> Dict[str, Array]:
    a = formulation.threshold_a
    b = formulation.threshold_b
    if formulation.d_map == "linear":
        scale = np.clip((s - formulation.d_cut) / max(1.0 - formulation.d_cut, 1e-9), 0.0, 1.0)
        d = formulation.d_floor + formulation.d_scale * scale**formulation.d_power
    elif formulation.d_map == "inverse":
        scale = np.clip((formulation.threshold_b - s) / max(formulation.threshold_b, 1e-9), 0.0, 1.0)
        d = formulation.d_floor + formulation.d_scale * scale**formulation.d_power
    elif formulation.d_map == "band":
        span = max(formulation.d_width, 1e-9)
        scale = np.clip(1.0 - np.abs(s - formulation.d_center) / span, 0.0, 1.0)
        d = formulation.d_floor + formulation.d_scale * scale**formulation.d_power
    else:
        raise ValueError(f"Unknown D map: {formulation.d_map}")

    reaction_mask = (s >= a) & (s <= b)
    band_center = 0.5 * (a + b)
    band_span = max(0.5 * (b - a), 1e-9)
    band_shape = np.clip(1.0 - np.abs(s - band_center) / band_span, 0.0, 1.0)
    kappa = formulation.kappa_scale * reaction_mask.astype(float) * (0.35 + 0.65 * band_shape)
    if instance_kind == "down":
        w0, w0_coating, w0_trapped = build_down_source(s, d, reaction_mask, formulation)
    else:
        w0, w0_coating, w0_trapped = build_up_source(s, d, reaction_mask, formulation)
    return {"d": d, "kappa": kappa, "w0": w0, "w0_coating": w0_coating, "w0_trapped": w0_trapped}


def harmonic_mean(left: Array, right: Array) -> Array:
    out = np.zeros_like(left)
    denom = left + right
    mask = denom > 1e-14
    out[mask] = 2.0 * left[mask] * right[mask] / denom[mask]
    return out


def prepare_local_cache(d: Array, dx: float) -> Dict[str, Array | float]:
    east = harmonic_mean(d[:, :-1], d[:, 1:])
    south = harmonic_mean(d[:-1, :], d[1:, :])
    diag_weight = np.zeros_like(d)
    diag_weight[:, :-1] += east
    diag_weight[:, 1:] += east
    diag_weight[:-1, :] += south
    diag_weight[1:, :] += south
    return {"east": east, "south": south, "diag_weight": diag_weight, "dx2": dx * dx}


def apply_local_operator(u: Array, cache: Dict[str, Array | float]) -> Array:
    east = cache["east"]
    south = cache["south"]
    dx2 = float(cache["dx2"])
    out = np.zeros_like(u)
    jump_east = east * (u[:, 1:] - u[:, :-1])
    out[:, :-1] += jump_east
    out[:, 1:] -= jump_east
    jump_south = south * (u[1:, :] - u[:-1, :])
    out[:-1, :] += jump_south
    out[1:, :] -= jump_south
    return out / dx2


def local_neighbor_sum(u: Array, cache: Dict[str, Array | float]) -> Array:
    east = cache["east"]
    south = cache["south"]
    out = np.zeros_like(u)
    out[:, :-1] += east * u[:, 1:]
    out[:, 1:] += east * u[:, :-1]
    out[:-1, :] += south * u[1:, :]
    out[1:, :] += south * u[:-1, :]
    return out


def local_implicit_step(u_rhs: Array, cache: Dict[str, Array | float], dt: float, iterations: int, relaxation: float) -> Array:
    diag_weight = cache["diag_weight"]
    dx2 = float(cache["dx2"])
    diag = 1.0 + dt * diag_weight / dx2
    u = u_rhs.copy()
    for _ in range(iterations):
        rhs = u_rhs + (dt / dx2) * local_neighbor_sum(u, cache)
        proposal = rhs / diag
        u = relaxation * proposal + (1.0 - relaxation) * u
        u = np.clip(u, 0.0, None)
    return u


def label_components(mask: Array) -> List[Array]:
    n = mask.shape[0]
    visited = np.zeros_like(mask, dtype=bool)
    components: List[Array] = []
    for i in range(n):
        for j in range(n):
            if not mask[i, j] or visited[i, j]:
                continue
            stack = [(i, j)]
            visited[i, j] = True
            coords: List[Tuple[int, int]] = []
            while stack:
                ci, cj = stack.pop()
                coords.append((ci, cj))
                for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ni = ci + di
                    nj = cj + dj
                    if 0 <= ni < n and 0 <= nj < n and mask[ni, nj] and not visited[ni, nj]:
                        visited[ni, nj] = True
                        stack.append((ni, nj))
            component = np.array(coords, dtype=int)
            components.append(component)
    return components


def prepare_nonlocal_cache(d: Array, formulation: Formulation, dx: float) -> Dict[str, object]:
    active = d > max(formulation.active_d_threshold, formulation.d_floor + 1e-12)
    components = label_components(active)
    component_diameters: List[float] = []
    for component in components:
        if len(component) <= 1:
            component_diameters.append(0.0)
            continue
        dy = component[:, 0][:, None] - component[:, 0][None, :]
        dx_cells = component[:, 1][:, None] - component[:, 1][None, :]
        diameter = float(np.sqrt((dy * dy + dx_cells * dx_cells).max())) * dx
        component_diameters.append(diameter)
    return {"components": components, "diameters": component_diameters}


def apply_nonlocal_step(u: Array, cache: Dict[str, object], dt: float, rate: float) -> Array:
    relaxed = u.copy()
    decay = np.exp(-rate * dt)
    for component in cache["components"]:
        if len(component) == 0:
            continue
        values = relaxed[component[:, 0], component[:, 1]]
        mean_value = float(values.mean())
        relaxed_values = mean_value + decay * (values - mean_value)
        relaxed[component[:, 0], component[:, 1]] = relaxed_values
    return np.clip(relaxed, 0.0, None)


def pairwise_nearest_distance(mask_a: Array, mask_b: Array, dx: float) -> float:
    coords_a = np.argwhere(mask_a)
    coords_b = np.argwhere(mask_b)
    if len(coords_a) == 0 or len(coords_b) == 0:
        return float("inf")
    delta = coords_a[:, None, :] - coords_b[None, :, :]
    distances = np.sqrt((delta * delta).sum(axis=2))
    nearest = distances.min(axis=1)
    return float(np.median(nearest)) * dx


def transport_reach(formulation: Formulation, d: Array, cache: Dict[str, object]) -> float:
    if formulation.transport_kind == "nonlocal":
        diameters = cache.get("diameters", [])
        if not diameters:
            return 0.0
        return float(np.median(np.asarray(diameters)))
    return float(np.sqrt(max(float(d.max()), 0.0) * formulation.t_final))


def prepare_transport_cache(formulation: Formulation, d: Array, dx: float) -> Dict[str, object]:
    if formulation.transport_kind in {"local_explicit", "local_implicit"}:
        return prepare_local_cache(d, dx)
    if formulation.transport_kind == "nonlocal":
        return prepare_nonlocal_cache(d, formulation, dx)
    raise ValueError(f"Unknown transport kind: {formulation.transport_kind}")


def transport_step(u: Array, formulation: Formulation, cache: Dict[str, object], dt: float) -> Array:
    if formulation.transport_kind == "local_explicit":
        return np.clip(u + dt * apply_local_operator(u, cache), 0.0, None)
    if formulation.transport_kind == "local_implicit":
        return local_implicit_step(u, cache, dt, formulation.implicit_iterations, formulation.implicit_relaxation)
    if formulation.transport_kind == "nonlocal":
        return apply_nonlocal_step(u, cache, dt, formulation.nonlocal_rate)
    raise ValueError(f"Unknown transport kind: {formulation.transport_kind}")


def explicit_stability_margin(formulation: Formulation, d: Array, dx: float) -> float:
    if formulation.transport_kind != "local_explicit":
        return 0.0
    dt = formulation.t_final / formulation.steps
    return float(d.max()) * dt / (dx * dx)


def observable(u: Array, v: Array, a0: float, dx: float) -> float:
    return float((phi(u, a0) * v).sum()) * dx * dx


def simulate(formulation: Formulation, instance_kind: str, seed: int, source_enabled: bool = True) -> RunResult:
    n = formulation.n
    dx = formulation.l / formulation.n
    dt = formulation.t_final / formulation.steps
    s = make_structure_field(instance_kind, seed, formulation)
    coeffs = build_coefficients(s, formulation, instance_kind)
    d = coeffs["d"]
    kappa = coeffs["kappa"]
    w0 = coeffs["w0"] if source_enabled else np.zeros_like(coeffs["w0"])
    w0_coating = coeffs["w0_coating"] if source_enabled else np.zeros_like(coeffs["w0_coating"])
    w0_trapped = coeffs["w0_trapped"] if source_enabled else np.zeros_like(coeffs["w0_trapped"])
    cache = prepare_transport_cache(formulation, d, dx)

    u = np.zeros((n, n), dtype=float)
    v = formulation.v0 * (kappa > 0.0).astype(float)
    w = w0.copy()
    times = np.linspace(0.0, formulation.t_final, formulation.steps + 1)
    r = np.zeros(formulation.steps + 1, dtype=float)
    transport_drift = 0.0
    decay_factor = np.exp(-formulation.lambda_decay * dt)

    for step in range(formulation.steps):
        w_next = w * decay_factor
        source_increment = w - w_next
        w = w_next

        before_transport = u + source_increment
        after_transport = transport_step(before_transport, formulation, cache, dt)
        mass_before = float(before_transport.sum()) * dx * dx
        mass_after = float(after_transport.sum()) * dx * dx
        transport_drift = max(transport_drift, abs(mass_after - mass_before))

        activation = phi(after_transport, formulation.a0)
        sink = activation * v
        u = np.clip(after_transport - dt * sink, 0.0, None)
        v = np.clip(v + dt * (kappa * (activation * v * (1.0 - v) - formulation.mu * v)), 0.0, 1.0)
        r[step + 1] = observable(u, v, formulation.a0, dx)

    separation = pairwise_nearest_distance(w0 > 0.0, kappa > 0.0, dx)
    reach = transport_reach(formulation, d, cache)
    reservoir_mass = float(w0.sum()) * dx * dx
    reaction_mass = float((kappa > 0.0).sum()) * dx * dx
    return RunResult(
        times=times,
        observable=r,
        s=s,
        d=d,
        kappa=kappa,
        w0=w0,
        w0_coating=w0_coating,
        w0_trapped=w0_trapped,
        separation=separation,
        reach=reach,
        transport_drift=transport_drift,
        reservoir_mass=reservoir_mass,
        reaction_mass=reaction_mass,
    )


def normalized_curve(values: Array) -> Array:
    peak = float(values.max())
    if peak <= 1e-12:
        return np.zeros_like(values)
    return values / peak


def first_last_third_means(values: Array) -> Tuple[float, float]:
    size = len(values)
    third = max(size // 3, 1)
    early = float(values[:third].mean())
    late = float(values[-third:].mean())
    return early, late


def curve_distance(a: Array, b: Array) -> float:
    aa = normalized_curve(a)
    bb = normalized_curve(b)
    std_a = float(aa.std())
    std_b = float(bb.std())
    if std_a <= 1e-12 or std_b <= 1e-12:
        corr_term = 1.0
    else:
        corr = float(np.corrcoef(aa, bb)[0, 1])
        corr_term = 1.0 - np.clip(corr, -1.0, 1.0)
    l2_term = float(np.sqrt(np.mean((aa - bb) ** 2)))
    return corr_term + l2_term


def baseline_formulation() -> Formulation:
    return Formulation(
        transport_kind="local_implicit",
        d_map="linear",
        d_scale=0.018,
        d_floor=0.0,
        threshold_a=0.28,
        threshold_b=0.48,
        reservoir_cut=0.14,
        down_coating_fraction=0.8,
        down_coating_thickness=2,
        lambda_decay=0.14,
        mu=0.22,
        v0=0.78,
        t_final=14.0,
        steps=160,
        nonlocal_rate=2.8,
    )


def seeded_candidate_formulation() -> Formulation:
    return Formulation(
        transport_kind="local_implicit",
        d_map="linear",
        n=32,
        l=1.0,
        t_final=14.0,
        steps=170,
        corr_length_cells=3.4,
        threshold_a=0.28,
        threshold_b=0.47,
        reservoir_cut=0.14,
        d_scale=0.022,
        d_floor=0.0,
        d_cut=0.1,
        d_center=0.46,
        d_width=0.34,
        d_power=1.1,
        kappa_scale=1.1,
        w0_scale=1.3,
        w0_power=1.05,
        down_coating_fraction=0.85,
        down_coating_thickness=1,
        lambda_decay=0.18,
        a0=0.075,
        mu=0.30,
        v0=0.86,
        implicit_iterations=110,
        implicit_relaxation=0.85,
        nonlocal_rate=4.4,
        active_d_threshold=1e-6,
    )