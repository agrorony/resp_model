"""Stage 22 post-processing for registered model-time mappings."""
from __future__ import annotations

import numpy as np


R_LIT_DAYS = {"low": 0.36, "primary": 2.6, "high": 9.2}
PERIODS_DAYS = ((0.0, 3.0), (3.0, 7.0), (7.0, 10.0), (10.0, 15.0))


def day_to_step(days: float, r_lit: float, dt: float = 0.05) -> float:
    if days < 0 or r_lit <= 0 or dt <= 0:
        raise ValueError("Days must be nonnegative and r_lit and dt must be positive")
    return days * r_lit / dt


def step_to_day(step: float, r_lit: float, dt: float = 0.05) -> float:
    if step < 0 or r_lit <= 0 or dt <= 0:
        raise ValueError("Step must be nonnegative and r_lit and dt must be positive")
    return step * dt / r_lit


def period_means(curve: np.ndarray, r_lit: float, dt: float = 0.05) -> np.ndarray:
    curve = np.asarray(curve, dtype=float)
    times_days = np.arange(curve.size, dtype=float) * dt / r_lit
    means = []
    for start, end in PERIODS_DAYS:
        values = curve[(times_days >= start) & (times_days < end)]
        means.append(float(values.mean()) if values.size else float("nan"))
    return np.asarray(means)


def classify_period_shape(means: np.ndarray) -> str:
    values = np.asarray(means, dtype=float)
    if values.shape != (4,) or not np.all(np.isfinite(values)):
        raise ValueError("Four finite period means are required")
    maximum = float(values.max())
    if values[3] == maximum and values[3] > 1.1 * values[0]:
        return "rise"
    if (values[1] == maximum or values[2] == maximum) and values[3] < 0.9 * maximum:
        return "rise-then-fall"
    if values[0] == maximum and values[3] < 0.9 * values[0]:
        return "fall"
    return "flat"


def physiological_outputs(
    respiration: np.ndarray,
    biomass: np.ndarray,
    uptake_dt: np.ndarray,
    initial_om: float,
    r_lit: float,
    dt: float = 0.05,
) -> dict:
    respiration = np.asarray(respiration, dtype=float)
    biomass = np.asarray(biomass, dtype=float)
    uptake_dt = np.asarray(uptake_dt, dtype=float)
    if not (respiration.shape == biomass.shape == uptake_dt.shape):
        raise ValueError("Respiration, biomass, and uptake arrays must have matching shapes")
    if initial_om <= 0:
        raise ValueError("Initial organic matter must be positive")
    steps = np.arange(respiration.size)
    day15 = day_to_step(15.0, r_lit, dt)
    within_15_days = steps < day15
    day10_to15 = (steps * dt / r_lit >= 10.0) & (steps * dt / r_lit < 15.0)
    cumulative_respiration = float(respiration[within_15_days].sum() * dt)
    cumulative_uptake = float(uptake_dt[within_15_days].sum())
    if cumulative_uptake == 0:
        raise ValueError("Cumulative substrate uptake is zero; CUE is undefined")
    if not np.any(day10_to15):
        raise ValueError("No model steps map to days 10-15")
    qco2_model = float(np.mean(respiration[day10_to15] / biomass[day10_to15]))
    return {
        "cumulative_respiration_model_units": cumulative_respiration,
        "cumulative_respiration_fraction_initial_OM": cumulative_respiration / initial_om,
        "qCO2_per_day": qco2_model * r_lit,
        "CUE": 1.0 - cumulative_respiration / cumulative_uptake,
    }
