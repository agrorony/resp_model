# Claude Code task: Stage 15 — reintegrating Loess (three-soil final demonstration)

Stage 10 explicitly dropped Loess to focus on getting Sand and Vertisol's
per-soil-tuned framework right; Stages 11–14 never brought it back. This stage
attempts a genuine three-soil version of Stage 10/11's §3 pass/fail
demonstration, now (per Stage 13) on real measured data. This is the
highest-uncertainty item in this continuation — earlier attempts at a three-soil
simultaneous fit (Stages 7–9) failed on distinctness. Budget-cap the search and
report an honest failure if it doesn't work, exactly as Stages 2–9 did — do not
force a fit by quietly loosening the success criteria.

## 1. What "success" means here

Reuse Stage 10 §3's explicit numeric pass/fail table, extended with a third row
for Loess's target shape ("low/erratic" per the measured experiment — see
`Soil_CO2_Respiration_v4.xlsx` Results sheet and `docs/stage_results/
STAGE3_5_RESULTS.md`'s framing of what "erratic" means operationally). Reuse
Stage 12's pairwise distinctness metric extended to three soils (three pairwise
distances, all > 0.3, time-normalized since windows may differ per soil as they
did for Sand/Vertisol in Stage 10).

## 2. Per-soil knobs, not a shared search

Follow Stage 10's precedent: give Loess its own `om_total`, coating fraction,
`D_scale`, and observation window `T`, tuned directly against the pass/fail table
— the same lever set already proven to work for Sand and Vertisol, not a new
mechanism. Use Loess's real PSD/connectivity data from Stage 13 as the structural
input, exactly as Sand/Vertisol now do.

## 3. Budget-capped search

Cap at `MAX_ITERATIONS=15` directed probes (Stage 7/9 style), each checked at
BOTH n=40 and n=128 as Stage 9 required (the two grids disagreed there — don't
assume they'll agree here). If Loess's shape can be hit individually but breaks
three-way distinctness against Sand and/or Vertisol (the Stage 7/9 failure mode:
soils blooming into the same shape once growth clears at all), name that
mechanism explicitly rather than chasing a 16th probe.

## 4. Honest-failure path (expected to be a live possibility)

If three-way success is not reached within budget: report which soil pair(s)
failed to separate, by how much, and the best mechanistic explanation available
(reuse Stage 9's diagnosis — Loess and Vertisol are both large, well-connected,
many-region soils whose recruitment is fast and front-loaded under a shared
transport regime — check whether that still applies now that both soils use real
rather than literature PSD data). A negative result here, clearly explained, is a
legitimate and reportable deliverable — consistent with how Stages 2–9 are
written up.

## 5. Deliverables

- `stage15_run.py`: three-soil Loess-reintegration runner, the budget-capped
  search loop and its trace, the headline three-soil comparison figure
  (`R(t)` normalized, all three), the §3-extended pass/fail table.
- `results/stage15/`: headline CSVs, comparison figures, search trace,
  distinctness table (n=128 and n=40).
- `docs/stage_results/STAGE15_RESULTS.md`: explicit success/failure verdict per
  soil and per pairwise distinctness check, with the named mechanism either way.
- Update `MODEL_SPEC.md`, `LOGBOOK.md`, `report_prompt.md`.

## 6. Order of work

1. Build Loess's real-data structure config (reuse Stage 13's loader).
2. Extend the §3 pass/fail table and distinctness metric to three soils.
3. Budget-capped search (≤15 iterations, both grids) for Loess's own knobs.
4. Three-soil headline run at whatever configuration search settles on.
5. Report success or the specific, mechanistically-named failure — do not loosen
   the 0.3 distinctness threshold or the pass/fail table to manufacture success.
6. STOP; write `STAGE15_RESULTS.md`.
