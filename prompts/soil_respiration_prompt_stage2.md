# Claude Code task: Stage 2 — water-filled pores & percolation

You are extending the EXISTING v2 project (pore-based model: `pore_field.py`,
`model.py`, `metrics.py`, `run.py`, `search_loop.py`, `configs/`, `MODEL_SPEC.md`,
`LOGBOOK.md`, `results/`). Stage 1 (fully wet) is done. Do NOT rewrite it —
evolve it, keep it runnable at every step.

## 0. Budget — read first

```
MAX_ITERATIONS = 12   # per the search loop; hard ceiling for this session
```
Stop and write up whatever you have when the budget or the Stage-2 checkpoint is
reached. Do not grind past it.

## 1. Success criterion & seeds (read before the physics)

The research question is simply: **do different soil structures produce different
respiration patterns?** So the goal is distinct curves, and a single realization
is fine — different soils SHOULD keep different seeds (they are meant to be
genuinely different structures; do NOT force a shared seed across the three
soils).

**Success** = the three `R(t)` curves are distinct — all pairwise distances
`> 0.3` — and every soil is non-trivial (`B_final_total > B0_total`,
`sum(R_t) > 0`). One realization per soil is acceptable.

Honest context from Stage 1 (state it, but don't over-engineer around it): at full
saturation only soil **B** separated clearly; A and C came out close (distance
~0.16, and only a lucky seed pushed the configured pair to 0.33). That is an
acceptable *negative* result — those two structures just don't produce clearly
different respiration when every pore conducts. The point of Stage 2 is to test
whether adding water/percolation gives genuinely different structures a strong
enough lever to produce genuinely different respiration.

Do not cherry-pick seeds or raise the threshold to fake a pass. If structures
don't separate, report it as a result.

**Later (optional, NOT this stage):** a robustness/significance pass may be added
to quantify how much the structural parameters actually matter — e.g. seed
ensembles, and testing whether different spatial *arrangements* of the same
texture statistics converge to the same respiration. Leave the door open for it,
but do not implement it now.

## 2. Invariants (unchanged)

- `configs/biology.yaml` and `configs/mapping.yaml` are **frozen** — never edited.
- The only per-soil freedom is the `pore` block and `saturation` in
  `configs/soil_*.yaml`.
- Deterministic; grid <= 50x50; seconds per run; budget-capped.

## 3. Stage-2 physics — water-filled pores & percolation

The model already supports this via `saturation.mode: global_theta`
(`model.water_mask_and_theta`): small pores fill first (capillarity), so a cell
is water-filled if `r <= r_cut`, where `r_cut` is that soil's own `theta`-quantile
of `r`. Air-filled cells get `D = 0`; harmonic-mean faces make a dry cell truly
block transport. Switch the three soils to `global_theta`.

**The mechanism to expect (state it, then test it):** as a soil dries, only the
*fine* pores stay water-filled — and fine pores have *low* D. So the conducting
network gets both slower and more fragmented as `theta` drops, and at some point
the water-filled path from OM-rich cells (finest pores) to high-K/biomass cells
(mid pores) **disconnects**. Whether and when that happens depends on the pore
distribution and spatial arrangement — a much stronger, less-tunable lever than
Stage 1's fully-wet transport-speed difference. This is your original saturation
idea ("the largest water-filled pore controls flow") made physical.

## 4. Do a theta sweep, not just one value

For each soil, sweep `theta` across a range (e.g. 1.0, 0.8, 0.6, 0.4, 0.2) and
record how `R(t)`, cumulative CO2, and the connectivity diagnostics (§5) change.
The interesting result is how the three soils *diverge as they dry*, so present
respiration vs theta, not only a single saturation.

## 5. Add percolation / connectivity diagnostics

New functions (numpy only — see §8 note on scipy) reporting, per soil per theta:
- **Largest connected water cluster** as a fraction of water-filled cells
  (4-connectivity labelling of the water mask).
- **OM-to-biomass connectivity**: are the high-OM (finest-pore) cells in the same
  connected water component as the high-K (biomass) cells? Report the fraction of
  initial biomass that sits in a water component containing appreciable OM.
- **Accessible-OM fraction**: total OM reachable through the water network from
  any biomass cell, over total OM.

Log these alongside the respiration metrics — they are the physical explanation
for whatever the curves do.

## 6. Honest framing (important)

Do NOT force three-way distinctness. Report what emerges. If only some soils
separate, say so, and explain it via the connectivity diagnostics. Compare Stage 2
to Stage 1's weak A-C result: did percolation give genuinely different structures
a strong enough lever to produce genuinely different respiration where fully-wet
transport could not?

## 7. Deliverables (add to / update the project)

- `metrics.py` — keep the emergent-distinctness check (§1); add the connectivity
  diagnostics (§5).
- A Stage-2 runner (extend `run.py` or add `stage2_run.py`) — the theta sweep,
  writing curves + diagnostics to `results/`.
- `configs/soil_*.yaml` — shared seed; `saturation.mode: global_theta` (+ the
  swept theta handled by the runner).
- `results/` — per soil: **water-filled maps** (at a chosen theta), respiration
  curves at that theta, and sweep plots: **respiration/cumulative vs theta** and
  **connectivity vs theta**.
- `LOGBOOK.md` — append every iteration and the Stage-2 checkpoint: per-soil pore
  params, pairwise distances, connectivity numbers, honest verdict.
- `STAGE2_RESULTS.md` — mirror `STAGE1_RESULTS.md`: idea, method, the numpy
  conversion, what emerged, connectivity explanation, limitations.
- `MODEL_SPEC.md` — fill in the Stage-2 details (it already scaffolds §9
  `global_theta`); record that the repo is numpy-only and that success = distinct
  curves (single realization).
- `report_prompt.md` — extend so the final report covers Stages 1–2, the seed
  caveat and fix, and interprets distinctness via percolation.

## 8. Convert the whole repo to numpy-only (required)

Remove the scipy dependency across the ENTIRE repo. Currently `pore_field.py`
uses `scipy.ndimage.gaussian_filter`; replace it with a numpy separable-Gaussian
convolution — a 1D Gaussian kernel applied along each axis with reflect padding,
`truncate=4.0` so `radius = int(4*sigma + 0.5)`, matching scipy's defaults as
closely as possible. Scan every module for any other scipy (or other non-numpy
heavy) imports and replace them too, so the project runs with `numpy` + `pyyaml`
+ `matplotlib` alone. Verify the Stage-1 curves are essentially unchanged after
the swap (tiny numeric differences from the kernel are fine; the shapes and
distinctness must hold), and update the README / any requirements accordingly.

## 9. Order of work

1. Convert the repo to numpy-only (§8); confirm the Stage-1 curves still hold.
2. Switch soils to `global_theta`; add connectivity diagnostics (§5); verify mass
   conservation still holds with `D=0` regions (it already does).
3. Theta sweep (§4); figures; distinctness per §1.
4. Checkpoint: STOP, write `STAGE2_RESULTS.md`, update `MODEL_SPEC.md`,
   `LOGBOOK.md`, `report_prompt.md`.

Remember: biology and mapping frozen; only the pore field + saturation differ;
different soils keep different seeds; success = distinct curves (single
realization fine), never targeted shapes; the repo is numpy-only; stop at the
budget and the checkpoint.
