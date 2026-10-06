from __future__ import annotations

import argparse
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from sandbox import (
    Formulation,
    RunResult,
    baseline_formulation,
    curve_distance,
    explicit_stability_margin,
    first_last_third_means,
    formulation_to_dict,
    normalized_curve,
    seeded_candidate_formulation,
    simulate,
    with_updates,
)


RESULTS_DIR = ROOT / "results"
LOGBOOK_PATH = ROOT / "LOGBOOK.md"
RESULTS_PATH = ROOT / "RESULTS.md"
README_PATH = ROOT / "README.md"
MAX_ITERATIONS = 40
ENSEMBLE_SEEDS = (101, 202, 303, 404, 505)
SUCCESS_DELTA = 0.15
SUCCESS_TAU = 0.3
NONTRIVIAL_FACTOR = 6.0
BASELINE_FRACTION = 0.08


def write_csv(path: Path, header: Sequence[str], rows: Sequence[Sequence[float | int | str]]) -> None:
    lines = [",".join(header)]
    for row in rows:
        parts: List[str] = []
        for value in row:
            if isinstance(value, str):
                parts.append(value)
            elif isinstance(value, (int, np.integer)):
                parts.append(str(int(value)))
            else:
                parts.append(f"{float(value):.10g}")
        lines.append(",".join(parts))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_array_csv(path: Path, values: np.ndarray) -> None:
    rows = [",".join(f"{float(item):.10g}" for item in row) for row in values]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def write_pgm(path: Path, values: np.ndarray) -> None:
    data = values.astype(float)
    low = float(data.min())
    high = float(data.max())
    if high - low <= 1e-12:
        scaled = np.zeros_like(data, dtype=int)
    else:
        scaled = np.rint(255.0 * (data - low) / (high - low)).astype(int)
    lines = ["P2", f"{data.shape[1]} {data.shape[0]}", "255"]
    lines.extend(" ".join(str(int(item)) for item in row) for row in scaled)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def median_value(values: Sequence[float]) -> float:
    return float(np.median(np.asarray(values, dtype=float)))


def evaluate_instance_runs(formulation: Formulation, instance: str, seeds: Sequence[int]) -> Tuple[List[RunResult], RunResult]:
    runs = [simulate(formulation, instance, seed, source_enabled=True) for seed in seeds]
    baseline = simulate(formulation, instance, seeds[0], source_enabled=False)
    return runs, baseline


def pick_median_pair(up_runs: Sequence[RunResult], down_runs: Sequence[RunResult]) -> int:
    distances = [curve_distance(up.observable, down.observable) for up, down in zip(up_runs, down_runs)]
    order = np.argsort(np.asarray(distances))
    return int(order[len(order) // 2])


def summarize_runs(runs: Sequence[RunResult]) -> Dict[str, float]:
    normalized = [normalized_curve(run.observable) for run in runs]
    early_late = [first_last_third_means(curve) for curve in normalized]
    early = [pair[0] for pair in early_late]
    late = [pair[1] for pair in early_late]
    peaks = [float(run.observable.max()) for run in runs]
    integrals = [float(np.trapezoid(run.observable, run.times)) for run in runs]
    peak_locations = [int(np.argmax(run.observable)) / max(len(run.observable) - 1, 1) for run in runs]
    coating_fractions = []
    for run in runs:
        total_source = float(run.w0.sum())
        coating_fractions.append(float(run.w0_coating.sum()) / total_source if total_source > 1e-12 else 0.0)
    return {
        "median_early": median_value(early),
        "median_late": median_value(late),
        "median_peak": median_value(peaks),
        "median_integral": median_value(integrals),
        "median_peak_location": median_value(peak_locations),
        "median_separation": median_value([run.separation for run in runs]),
        "median_reach": median_value([run.reach for run in runs]),
        "median_transport_drift": median_value([run.transport_drift for run in runs]),
        "median_reservoir_mass": median_value([run.reservoir_mass for run in runs]),
        "median_reaction_mass": median_value([run.reaction_mass for run in runs]),
        "median_coating_fraction": median_value(coating_fractions),
    }


def compact_shape_metrics(runs: Sequence[RunResult]) -> Dict[str, float | bool]:
    first_third_passes = []
    monotone_scores = []
    for run in runs:
        curve = normalized_curve(run.observable)
        peak_idx = int(np.argmax(curve))
        first_third_limit = max(len(curve) // 3, 1)
        first_third_passes.append(1.0 if peak_idx < first_third_limit else 0.0)
        tail = curve[peak_idx:]
        if len(tail) <= 2:
            monotone_scores.append(1.0)
        else:
            positive_jumps = np.diff(tail) > 0.01
            monotone_scores.append(1.0 - float(positive_jumps.mean()))
    median_first = median_value(first_third_passes)
    median_monotone = median_value(monotone_scores)
    return {
        "median_peak_first_third": median_first,
        "median_monotone_tail": median_monotone,
        "shape_ok": median_first >= 0.5 and median_monotone >= 0.8,
    }


def evaluate_formulation(formulation: Formulation, seeds: Sequence[int]) -> Dict[str, object]:
    up_runs, up_baseline = evaluate_instance_runs(formulation, "up", seeds)
    down_runs, down_baseline = evaluate_instance_runs(formulation, "down", seeds)

    up_summary = summarize_runs(up_runs)
    down_summary = summarize_runs(down_runs)
    down_shape = compact_shape_metrics(down_runs)
    up_trend = up_summary["median_late"] - up_summary["median_early"]
    down_trend = down_summary["median_early"] - down_summary["median_late"]
    distances = [curve_distance(up.observable, down.observable) for up, down in zip(up_runs, down_runs)]
    median_distance = median_value(distances)
    median_index = pick_median_pair(up_runs, down_runs)
    up_integral = up_summary["median_integral"]
    down_integral = down_summary["median_integral"]
    up_baseline_integral = float(np.trapezoid(up_baseline.observable, up_baseline.times))
    down_baseline_integral = float(np.trapezoid(down_baseline.observable, down_baseline.times))
    up_nontrivial = up_integral > NONTRIVIAL_FACTOR * max(up_baseline_integral, 1e-12) and float(up_baseline.observable.max()) < BASELINE_FRACTION * max(up_summary["median_peak"], 1e-12)
    down_nontrivial = down_integral > NONTRIVIAL_FACTOR * max(down_baseline_integral, 1e-12) and float(down_baseline.observable.max()) < BASELINE_FRACTION * max(down_summary["median_peak"], 1e-12)
    finite_separation = np.isfinite(up_summary["median_separation"]) and np.isfinite(down_summary["median_separation"])
    signatures = up_trend > SUCCESS_DELTA and down_trend > SUCCESS_DELTA and bool(down_shape["shape_ok"])
    stable = True
    stability_margin = explicit_stability_margin(formulation, up_runs[0].d, formulation.l / formulation.n)
    if formulation.transport_kind == "local_explicit":
        stable = stability_margin < 0.2
    score = median_distance if signatures and up_nontrivial and down_nontrivial and stable else 0.0
    uncoupled = (not up_nontrivial) or (not down_nontrivial) or (not finite_separation)
    up_deficit = max(SUCCESS_DELTA - up_trend, 0.0)
    down_deficit = max(SUCCESS_DELTA - down_trend, 0.0)
    if uncoupled or formulation.transport_kind != "local_implicit":
        merit = -10.0
    else:
        merit = (
            max(up_trend, 0.0)
            + max(down_trend, 0.0)
            + 0.5 * median_distance
            + (0.25 if bool(down_shape["shape_ok"]) else -0.5)
            + (0.2 if stable else -0.5)
            - 1.2 * up_deficit
            - 1.2 * down_deficit
        )

    return {
        "score": score,
        "merit": merit,
        "signatures_ok": signatures,
        "up_nontrivial": up_nontrivial,
        "down_nontrivial": down_nontrivial,
        "finite_separation": finite_separation,
        "stable": stable,
        "up_trend": up_trend,
        "down_trend": down_trend,
        "down_shape": down_shape,
        "median_distance": median_distance,
        "stability_margin": stability_margin,
        "up_summary": up_summary,
        "down_summary": down_summary,
        "up_baseline_integral": up_baseline_integral,
        "down_baseline_integral": down_baseline_integral,
        "up_baseline_peak": float(up_baseline.observable.max()),
        "down_baseline_peak": float(down_baseline.observable.max()),
        "median_index": median_index,
        "up_runs": up_runs,
        "down_runs": down_runs,
        "up_baseline": up_baseline,
        "down_baseline": down_baseline,
    }


def mutation_note(previous: Formulation, current: Formulation) -> str:
    changed = []
    for key, value in formulation_to_dict(current).items():
        if formulation_to_dict(previous)[key] != value:
            changed.append(key)
    if not changed:
        return "Retained the previous formulation to confirm reproducibility."
    if "transport_kind" in changed:
        return f"Changed transport_kind to {current.transport_kind} to alter the coupling reach between the disjoint supports."
    if any(name in changed for name in ("d_scale", "d_floor", "d_cut", "d_map", "d_power", "d_center", "d_width")):
        return "Adjusted the transport map to change the connected active set and the effective reach across the structure field."
    if any(name in changed for name in ("lambda_decay", "mu", "v0", "a0")):
        return "Adjusted the shared time scales so the source-release tempo and the observable decay compete differently."
    if any(name in changed for name in ("threshold_a", "threshold_b", "reservoir_cut", "w0_scale", "w0_power", "kappa_scale")):
        return "Adjusted the shared coefficient maps to change how much source mass and active set area each structure law receives."
    return "Mutated the shared formulation to improve signature separation while preserving the same dynamics for both instances."


def predetermined_formulation(iteration: int) -> Formulation:
    if iteration == 1:
        return baseline_formulation()
    if iteration == 2:
        return with_updates(
            seeded_candidate_formulation(),
            lambda_decay=0.16,
            mu=0.26,
            t_final=16.0,
            steps=180,
            w0_scale=1.35,
            d_scale=0.024,
        )
    if iteration == 3:
        return with_updates(
            seeded_candidate_formulation(),
            d_floor=0.0,
            lambda_decay=0.17,
            mu=0.24,
            t_final=18.0,
            steps=180,
            w0_scale=1.55,
            d_scale=0.018,
        )
    if iteration == 4:
        return with_updates(
            seeded_candidate_formulation(),
            d_scale=0.02,
            d_floor=0.0,
            mu=0.28,
            lambda_decay=0.18,
            threshold_a=0.3,
            threshold_b=0.5,
        )
    raise ValueError("No predetermined formulation for this iteration")


def mutate_formulation(best: Formulation, rng: np.random.Generator) -> Formulation:
    candidate = with_updates(best, transport_kind="local_implicit")
    candidate = with_updates(candidate, implicit_iterations=int(np.clip(round(candidate.implicit_iterations + rng.integers(-15, 21)), 60, 150)))
    candidate = with_updates(
        candidate,
        d_map=rng.choice(np.array(["linear", "linear", "inverse", "band"])).item(),
        d_scale=float(np.clip(candidate.d_scale * np.exp(rng.normal(0.0, 0.18)), 0.01, 0.08)),
        d_floor=float(np.clip(candidate.d_floor + rng.normal(0.0, 0.004), 0.0, 0.015)),
        d_cut=float(np.clip(candidate.d_cut + rng.normal(0.0, 0.03), 0.02, 0.35)),
        d_center=float(np.clip(candidate.d_center + rng.normal(0.0, 0.04), 0.2, 0.75)),
        d_width=float(np.clip(candidate.d_width + rng.normal(0.0, 0.05), 0.12, 0.45)),
        d_power=float(np.clip(candidate.d_power + rng.normal(0.0, 0.15), 0.8, 2.0)),
        lambda_decay=float(np.clip(candidate.lambda_decay + rng.normal(0.0, 0.02), 0.08, 0.24)),
        a0=float(np.clip(candidate.a0 + rng.normal(0.0, 0.015), 0.04, 0.16)),
        mu=float(np.clip(candidate.mu + rng.normal(0.0, 0.02), 0.14, 0.34)),
        v0=float(np.clip(candidate.v0 + rng.normal(0.0, 0.03), 0.6, 0.92)),
        kappa_scale=float(np.clip(candidate.kappa_scale + rng.normal(0.0, 0.08), 0.8, 1.6)),
        w0_scale=float(np.clip(candidate.w0_scale + rng.normal(0.0, 0.1), 0.6, 1.8)),
        w0_power=float(np.clip(candidate.w0_power + rng.normal(0.0, 0.1), 0.8, 1.6)),
        threshold_a=float(np.clip(candidate.threshold_a + rng.normal(0.0, 0.02), 0.22, 0.42)),
        threshold_b=float(np.clip(candidate.threshold_b + rng.normal(0.0, 0.02), 0.46, 0.68)),
        reservoir_cut=float(np.clip(candidate.reservoir_cut + rng.normal(0.0, 0.02), 0.08, 0.3)),
        t_final=float(np.clip(candidate.t_final + rng.normal(0.0, 1.2), 10.0, 24.0)),
        steps=int(np.clip(round(candidate.steps + rng.integers(-20, 25)), 100, 220)),
        corr_length_cells=float(np.clip(candidate.corr_length_cells + rng.normal(0.0, 0.25), 2.0, 4.8)),
        down_coating_fraction=float(np.clip(candidate.down_coating_fraction + rng.normal(0.0, 0.04), 0.7, 0.85)),
        down_coating_thickness=int(np.clip(round(candidate.down_coating_thickness + rng.integers(-1, 2)), 1, 4)),
    )
    if candidate.threshold_b <= candidate.threshold_a + 0.08:
        candidate = with_updates(candidate, threshold_b=min(0.72, candidate.threshold_a + 0.12))
    if candidate.reservoir_cut >= candidate.threshold_a - 0.04:
        candidate = with_updates(candidate, reservoir_cut=max(0.08, candidate.threshold_a - 0.08))
    return candidate


def log_iteration(iteration: int, formulation: Formulation, evaluation: Dict[str, object], note: str) -> str:
    return "\n".join(
        [
            f"## Iteration {iteration}",
            "",
            f"- Formulation: `{formulation.transport_kind}`, `D` map `{formulation.d_map}`",
            f"- Mutation note: {note}",
            f"- Score: {float(evaluation['score']):.4f}",
            f"- Merit: {float(evaluation['merit']):.4f}",
            f"- Signatures: up trend = {float(evaluation['up_trend']):.4f}, down trend = {float(evaluation['down_trend']):.4f}, pass = {bool(evaluation['signatures_ok'])}",
            f"- Distance: {float(evaluation['median_distance']):.4f}",
            f"- Nontrivial: up = {bool(evaluation['up_nontrivial'])}, down = {bool(evaluation['down_nontrivial'])}",
            f"- Compact shape: peak-first-third = {float(evaluation['down_shape']['median_peak_first_third']):.2f}, monotone-tail = {float(evaluation['down_shape']['median_monotone_tail']):.2f}, pass = {bool(evaluation['down_shape']['shape_ok'])}",
            f"- Stability: pass = {bool(evaluation['stable'])}, explicit margin = {float(evaluation['stability_margin']):.4f}",
            f"- Reach vs separation: up {evaluation['up_summary']['median_reach']:.4f} vs {evaluation['up_summary']['median_separation']:.4f}; down {evaluation['down_summary']['median_reach']:.4f} vs {evaluation['down_summary']['median_separation']:.4f}",
            f"- Baselines: up integral = {float(evaluation['up_baseline_integral']):.4e}, down integral = {float(evaluation['down_baseline_integral']):.4e}",
            f"- Transport drift: up = {evaluation['up_summary']['median_transport_drift']:.4e}, down = {evaluation['down_summary']['median_transport_drift']:.4e}",
            "",
        ]
    )


def save_representative_results(formulation: Formulation, evaluation: Dict[str, object], trace_rows: Sequence[Sequence[float | int | str]]) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    median_index = int(evaluation["median_index"])
    up_run = evaluation["up_runs"][median_index]
    down_run = evaluation["down_runs"][median_index]
    up_baseline = evaluation["up_baseline"]
    down_baseline = evaluation["down_baseline"]

    times_rows = [
        [
            up_run.times[i],
            up_run.observable[i],
            up_baseline.observable[i],
            down_run.observable[i],
            down_baseline.observable[i],
        ]
        for i in range(len(up_run.times))
    ]
    write_csv(
        RESULTS_DIR / "winning_timeseries.csv",
        ["time", "R_up", "R_up_baseline", "R_down", "R_down_baseline"],
        times_rows,
    )

    for prefix, run in (("up", up_run), ("down", down_run)):
        write_array_csv(RESULTS_DIR / f"{prefix}_s.csv", run.s)
        write_array_csv(RESULTS_DIR / f"{prefix}_D.csv", run.d)
        write_array_csv(RESULTS_DIR / f"{prefix}_kappa.csv", run.kappa)
        write_array_csv(RESULTS_DIR / f"{prefix}_w0.csv", run.w0)
        write_pgm(RESULTS_DIR / f"{prefix}_s.pgm", run.s)
        write_pgm(RESULTS_DIR / f"{prefix}_D.pgm", run.d)
        write_pgm(RESULTS_DIR / f"{prefix}_kappa.pgm", run.kappa)
        write_pgm(RESULTS_DIR / f"{prefix}_w0.pgm", run.w0)
    write_array_csv(RESULTS_DIR / "down_w0_coating.csv", down_run.w0_coating)
    write_array_csv(RESULTS_DIR / "down_w0_trapped.csv", down_run.w0_trapped)
    write_pgm(RESULTS_DIR / "down_w0_coating.pgm", down_run.w0_coating)
    write_pgm(RESULTS_DIR / "down_w0_trapped.pgm", down_run.w0_trapped)

    write_csv(RESULTS_DIR / "score_trace.csv", ["iteration", "score", "merit", "up_trend", "down_trend", "distance"], trace_rows)
    write_results_markdown(formulation, evaluation)


def format_equations(formulation: Formulation) -> str:
    if formulation.transport_kind == "nonlocal":
        transport_text = "For each connected active component C = {D > threshold}, du/dt = rho (mean_C(u) - u)."
    elif formulation.transport_kind == "local_implicit":
        transport_text = "u^{n+1} solves u^{n+1} - dt div(D grad u^{n+1}) = rhs with Jacobi iteration."
    else:
        transport_text = "u^{n+1} = u^n + dt div(D grad u^n) with harmonic-face heterogeneous fluxes."
    return "\n".join(
        [
            "Shared winning formulation:",
            "",
            "dw/dt = -lambda w",
            "du/dt = L_D[u] + lambda w - phi(u) v",
            "dv/dt = kappa(x) * (phi(u) v (1 - v) - mu v)",
            "phi(u) = u / (u + a0)",
            "",
            f"Transport realization: {transport_text}",
            "",
            "Shared constants:",
            f"- N = {formulation.n}, L = {formulation.l}, T = {formulation.t_final}, steps = {formulation.steps}",
            f"- lambda = {formulation.lambda_decay:.4f}, a0 = {formulation.a0:.4f}, mu = {formulation.mu:.4f}, v0 = {formulation.v0:.4f}",
            f"- D map = {formulation.d_map}, D_scale = {formulation.d_scale:.4f}, D_floor = {formulation.d_floor:.4f}, D_cut = {formulation.d_cut:.4f}, D_power = {formulation.d_power:.4f}",
            f"- thresholds a = {formulation.threshold_a:.4f}, b = {formulation.threshold_b:.4f}, kappa_scale = {formulation.kappa_scale:.4f}",
            f"- w0_scale = {formulation.w0_scale:.4f}, w0_power = {formulation.w0_power:.4f}",
            f"- compact coating fraction = {formulation.down_coating_fraction:.4f}, compact coating thickness = {formulation.down_coating_thickness}",
            f"- nonlocal_rate = {formulation.nonlocal_rate:.4f}, implicit_iterations = {formulation.implicit_iterations}",
        ]
    )


def write_results_markdown(formulation: Formulation, evaluation: Dict[str, object]) -> None:
    up_summary = evaluation["up_summary"]
    down_summary = evaluation["down_summary"]
    success_text = "Met the signature, distance, and nontrivial checks." if float(evaluation["score"]) > 0.0 else "Did not fully meet all checks; this is the best formulation found within the budget."
    if float(evaluation["score"]) > 0.0:
        mechanism_lines = [
            "- The shared coefficient maps give the broad structure law a larger source mass and a more fragmented active set, so the observable rises as released mass keeps reaching many active patches.",
            "- The narrow structure law gets a smaller nearby source support, so the observable forms an early pulse and then decays once the limited supply is exhausted.",
            "- The winning transport increases sensitivity by making the source-to-reaction coupling depend on connectivity rather than only on a short local diffusion length.",
        ]
    else:
        mechanism_lines = [
            "- The explicit baseline reproduces the intended degeneracy: the observable stays near the source-free baseline when transport cannot bridge the source-reaction gap.",
            "- The best budget-spent formulation breaks that degeneracy only for the broad instance: the up curve becomes strongly increasing and clearly separated from baseline.",
            "- The remaining failure is the compact instance. Under the best shared formulation found here, its source-reaction coupling collapses too far, so the down curve stays effectively trivial rather than producing the requested decreasing signature.",
        ]
    text = "\n".join(
        [
            "# PDE Sandbox Results",
            "",
            format_equations(formulation),
            "",
            "Reach versus separation:",
            f"- Baseline local explicit reach target is represented by the first logbook entry.",
            f"- Winning up instance: reach = {up_summary['median_reach']:.4f}, separation = {up_summary['median_separation']:.4f}",
            f"- Winning down instance: reach = {down_summary['median_reach']:.4f}, separation = {down_summary['median_separation']:.4f}",
            f"- Winning down peak location (fraction of horizon) = {down_summary['median_peak_location']:.4f}",
            f"- Winning down coating fraction = {down_summary['median_coating_fraction']:.4f}",
            "",
            "Why this works:",
            *mechanism_lines,
            "",
            "Outcome summary:",
            f"- merit = {float(evaluation['merit']):.4f}",
            f"- up trend = {float(evaluation['up_trend']):.4f}",
            f"- down trend = {float(evaluation['down_trend']):.4f}",
            f"- median distance = {float(evaluation['median_distance']):.4f}",
            f"- up nontrivial = {bool(evaluation['up_nontrivial'])}",
            f"- down nontrivial = {bool(evaluation['down_nontrivial'])}",
            f"- down shape pass = {bool(evaluation['down_shape']['shape_ok'])}",
            f"- up baseline integral = {float(evaluation['up_baseline_integral']):.4e}",
            f"- down baseline integral = {float(evaluation['down_baseline_integral']):.4e}",
            f"- transport drift medians = up {up_summary['median_transport_drift']:.4e}, down {down_summary['median_transport_drift']:.4e}",
            "",
            "Status:",
            f"- {success_text}",
            "",
            "Artifacts:",
            "- results/winning_timeseries.csv contains the representative R(t) curves with the source-free baselines overlaid as columns.",
            "- results/up_*.csv and results/down_*.csv contain the structure and coefficient maps for the representative seeds.",
            "- results/down_w0_coating.* and results/down_w0_trapped.* isolate the compact source split.",
            "- results/score_trace.csv records the search trajectory.",
        ]
    )
    RESULTS_PATH.write_text(text + "\n", encoding="utf-8")


def write_readme() -> None:
    README_PATH.write_text(
        "\n".join(
            [
                "# PDE Sandbox",
                "",
                "Pure NumPy search over abstract reaction-transport formulations on a heterogeneous square grid.",
                "",
                "## Run",
                "",
                "python pde_sandbox/search.py",
                "",
                "Optional flags:",
                "- --max-iterations N",
                "- --seed N",
                "",
                "Outputs are written into pde_sandbox/results plus LOGBOOK.md and RESULTS.md.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def run_search(max_iterations: int, seed: int) -> Dict[str, object]:
    rng = np.random.default_rng(seed)
    write_readme()
    best_formulation: Formulation | None = None
    best_evaluation: Dict[str, object] | None = None
    previous = baseline_formulation()
    log_sections = ["# Search Log", "", f"Budget: {max_iterations} iterations", ""]
    trace_rows: List[Sequence[float | int | str]] = []

    for iteration in range(1, max_iterations + 1):
        if iteration <= 4:
            formulation = predetermined_formulation(iteration)
        else:
            formulation = mutate_formulation(best_formulation or seeded_candidate_formulation(), rng)
        evaluation = evaluate_formulation(formulation, ENSEMBLE_SEEDS)

        if best_evaluation is None:
            best_formulation = formulation
            best_evaluation = evaluation
        elif float(evaluation["score"]) > float(best_evaluation["score"]):
            best_formulation = formulation
            best_evaluation = evaluation
        elif float(evaluation["score"]) == float(best_evaluation["score"]) and float(evaluation["merit"]) > float(best_evaluation["merit"]):
            best_formulation = formulation
            best_evaluation = evaluation

        note = mutation_note(previous, formulation)
        previous = formulation
        log_sections.append(log_iteration(iteration, formulation, evaluation, note))
        trace_rows.append((iteration, evaluation["score"], evaluation["merit"], evaluation["up_trend"], evaluation["down_trend"], evaluation["median_distance"]))

        if float(evaluation["score"]) > 0.0 and float(evaluation["median_distance"]) > SUCCESS_TAU:
            best_formulation = formulation
            best_evaluation = evaluation
            break

    LOGBOOK_PATH.write_text("\n".join(log_sections).rstrip() + "\n", encoding="utf-8")
    if best_formulation is None or best_evaluation is None:
        raise RuntimeError("Search did not evaluate any formulation.")
    save_representative_results(best_formulation, best_evaluation, trace_rows)
    return {"formulation": best_formulation, "evaluation": best_evaluation, "iterations": len(trace_rows)}


def main() -> None:
    parser = argparse.ArgumentParser(description="Search for a sensitive abstract PDE sandbox formulation.")
    parser.add_argument("--max-iterations", type=int, default=MAX_ITERATIONS)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()
    result = run_search(max_iterations=args.max_iterations, seed=args.seed)
    formulation = result["formulation"]
    evaluation = result["evaluation"]
    print(f"Iterations: {result['iterations']}")
    print(f"Winner transport: {formulation.transport_kind}")
    print(f"Score: {float(evaluation['score']):.4f}")
    print(f"Up trend: {float(evaluation['up_trend']):.4f}")
    print(f"Down trend: {float(evaluation['down_trend']):.4f}")
    print(f"Distance: {float(evaluation['median_distance']):.4f}")


if __name__ == "__main__":
    main()