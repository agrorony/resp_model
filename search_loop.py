"""
Budget-capped structure search loop (soil_respiration_prompt_v2.md, New
success criterion + STAGES).

Each iteration: run soils A, B, C (biology and mapping.yaml frozen), score
emergent distinctness with metrics.emergent_distinctness, append a
LOGBOOK.md entry, and stop on success or when MAX_ITERATIONS is spent. If not
successful, nudge ONLY the pore-distribution params (mu, sigma, lambda,
aggregate) of whichever soil pair is least distinct, never biology or
mapping.yaml, and never hand-place OM/biomass.
"""
from __future__ import annotations

import datetime
import os

import yaml

from model import load_yaml, simulate
import metrics

MAX_ITERATIONS = 12
DISTINCTNESS_THRESHOLD = 0.3

CONFIG_PATH = "configs/soil_{}.yaml"
LOGBOOK_PATH = "LOGBOOK.md"


def run_all(biology_cfg, mapping_cfg):
    out = {}
    for soil in "ABC":
        st = load_yaml(CONFIG_PATH.format(soil))
        out[soil] = (simulate(biology_cfg, st, mapping_cfg), st)
    return out


def summarize_structure(st):
    p = st["pore"]
    return {
        "mu": p["mu"], "sigma": p["sigma"], "lambda": p["lambda"],
        "aggregate": p.get("aggregate", False),
        "saturation_mode": st.get("saturation", {}).get("mode", "fully_wet"),
        "T": st["T"],
    }


def _least_distinct_pair(distances):
    return min(distances, key=distances.get)


# mu is bounded so K never collapses to its floor everywhere (which would
# make growth trivial and, worse, make two collapsed soils *look* identical
# -- both just decay near zero). Bounds keep at least part of the grid
# within ~2*w_K of ln(r_opt)=1.61 (r_opt=5, w_K=0.45 in mapping.yaml) for
# every soil, so every soil retains a genuinely viable habitat fraction.
MU_BOUNDS = (0.8, 2.6)
LAMBDA_BOUNDS = (1.0, 10.0)
MU_STEP = 0.2
LAMBDA_STEP = 0.7
SIGMA_STEP = 0.08
SIGMA_BOUNDS = (0.2, 0.9)


def nudge_pair(pair, sts, iteration):
    """Take a small, bounded step that pushes the two least-distinct soils'
    pore params apart. Alternates which axis moves (texture mu, spatial
    lambda, or texture spread sigma / aggregate) across iterations so the
    search explores more than one axis instead of running one knob to its
    rails."""
    a, b = pair[0], pair[1]
    pa, pb = sts[a]["pore"], sts[b]["pore"]
    axis = ["mu", "lambda", "sigma"][iteration % 3]

    if axis == "mu":
        lo, hi = (a, b) if pa["mu"] <= pb["mu"] else (b, a)
        sts[lo]["pore"]["mu"] = max(sts[lo]["pore"]["mu"] - MU_STEP, MU_BOUNDS[0])
        sts[hi]["pore"]["mu"] = min(sts[hi]["pore"]["mu"] + MU_STEP, MU_BOUNDS[1])
        note = (f"soils {a},{b} least distinct: stepped mu apart to "
                f"({sts[lo]['pore']['mu']:.2f} vs {sts[hi]['pore']['mu']:.2f})")
    elif axis == "lambda":
        lo, hi = (a, b) if pa["lambda"] <= pb["lambda"] else (b, a)
        sts[lo]["pore"]["lambda"] = max(sts[lo]["pore"]["lambda"] - LAMBDA_STEP, LAMBDA_BOUNDS[0])
        sts[hi]["pore"]["lambda"] = min(sts[hi]["pore"]["lambda"] + LAMBDA_STEP, LAMBDA_BOUNDS[1])
        note = (f"soils {a},{b} least distinct: stepped lambda apart to "
                f"({sts[lo]['pore']['lambda']:.2f} vs {sts[hi]['pore']['lambda']:.2f})")
    else:
        lo, hi = (a, b) if pa["sigma"] <= pb["sigma"] else (b, a)
        sts[lo]["pore"]["sigma"] = max(sts[lo]["pore"]["sigma"] - SIGMA_STEP, SIGMA_BOUNDS[0])
        sts[hi]["pore"]["sigma"] = min(sts[hi]["pore"]["sigma"] + SIGMA_STEP, SIGMA_BOUNDS[1])
        # also toggle aggregate on the higher-sigma (more heterogeneous) soil
        sts[hi]["pore"]["aggregate"] = not sts[hi]["pore"].get("aggregate", False)
        note = (f"soils {a},{b} least distinct: stepped sigma apart to "
                f"({sts[lo]['pore']['sigma']:.2f} vs {sts[hi]['pore']['sigma']:.2f}), "
                f"toggled aggregate on soil {hi}")

    return sts, note


def append_logbook(iteration, results, sts, detail, success_flag, stage_label):
    lines = [f"\n## {stage_label} -- Iteration {iteration} -- {datetime.datetime.now().isoformat(timespec='seconds')}\n"]
    for soil in "ABC":
        st = sts[soil]
        struct = summarize_structure(st)
        shape = detail["shapes"][soil]
        lines.append(f"**Soil {soil}**\n")
        lines.append(f"- pore params: mu={struct['mu']:.3f}, sigma={struct['sigma']:.3f}, "
                     f"lambda={struct['lambda']:.3f}, aggregate={struct['aggregate']}, "
                     f"saturation={struct['saturation_mode']}\n")
        lines.append(f"- auto-described shape: {shape['shape']} (e={shape['e']:.3f}, l={shape['l']:.3f}, "
                     f"peak={shape['peak']:.3f}, n_peaks={shape['n_peaks']})\n")
        lines.append(f"- non-trivial: {detail['nontrivial'][soil]}\n")
    lines.append(f"- pairwise distances (threshold={detail['threshold']}): {detail['distances']}\n")
    lines.append(f"- note: {results.get('note', '')}\n")
    lines.append(f"\n**Iteration result: {'SUCCESS' if success_flag else 'not yet successful'}**\n")
    with open(LOGBOOK_PATH, "a", encoding="utf-8") as f:
        f.writelines(lines)


def run_search(stage_label="Stage 1"):
    if not os.path.exists(LOGBOOK_PATH):
        with open(LOGBOOK_PATH, "w", encoding="utf-8") as f:
            f.write("# LOGBOOK\n\nBudget-capped pore-distribution search "
                     "(soil_respiration_prompt_v2.md). MAX_ITERATIONS = "
                     f"{MAX_ITERATIONS}. Biology and configs/mapping.yaml are "
                     "frozen throughout; only configs/soil_{A,B,C}.yaml pore "
                     "params are ever modified by nudges.\n")

    success_flag = False
    detail = None
    for iteration in range(1, MAX_ITERATIONS + 1):
        biology_cfg = load_yaml("configs/biology.yaml")
        mapping_cfg = load_yaml("configs/mapping.yaml")
        run_out = run_all(biology_cfg, mapping_cfg)

        outs = {s: run_out[s][0] for s in "ABC"}
        sts = {s: run_out[s][1] for s in "ABC"}

        success_flag, detail = metrics.emergent_distinctness(
            outs["A"], outs["B"], outs["C"], threshold=DISTINCTNESS_THRESHOLD,
        )

        results = {"note": "target met, no change needed" if success_flag else ""}

        if success_flag:
            append_logbook(iteration, results, sts, detail, success_flag, stage_label)
            print(f"[{stage_label}] Success at iteration {iteration}. min_distance={detail['min_distance']:.3f}")
            break

        pair = _least_distinct_pair(detail["distances"])
        sts, note = nudge_pair(pair, sts, iteration)
        results["note"] = note
        for soil in "ABC":
            with open(CONFIG_PATH.format(soil), "w") as f:
                yaml.safe_dump(sts[soil], f, sort_keys=False)

        append_logbook(iteration, results, sts, detail, success_flag, stage_label)
        print(f"[{stage_label}] Iteration {iteration}: not yet successful "
              f"(min_distance={detail['min_distance']:.3f}), nudged {pair}.")

    if not success_flag:
        with open(LOGBOOK_PATH, "a", encoding="utf-8") as f:
            f.write(f"\n**[{stage_label}] Budget exhausted after {MAX_ITERATIONS} "
                    "iterations without full success. See the final iteration "
                    "above for the best configuration reached.**\n")
        print(f"[{stage_label}] Budget of {MAX_ITERATIONS} iterations exhausted without full success.")

    return success_flag, detail


if __name__ == "__main__":
    run_search()
