# Claude Code task: Stage 12 — statistical significance via seed ensembles

An ANALYSIS stage, not a modeling stage. Build on Stage 11: keep the two soils
(Sand, Vertisol) and their Stage-11 parameters EXACTLY (om_total, coating fraction,
D_scale, biology, T — all frozen). Do NOT re-tune the model. The goal is to show,
statistically, that the simple structural features we defined — not a lucky random
realization — drive each soil's respiration shape. numpy-only; n=128; use the
Stage-11 optimized code (fast).

## 1. Seed ensembles

For EACH soil, run an ensemble of **M random seeds** (M ≥ 30; go higher if runtime
allows — Stage 11 made single runs ~3–20 s). Each seed generates a fresh pore-field
realization with the SAME structural parameters (the soil's frozen Stage-11 config);
ONLY the pore-field generation seed changes between ensemble members. Nothing else —
not om_total, coating, D_scale, or biology — varies within a soil's ensemble. This is
the whole point: same structural class, different random layout.

Collect `R(t)` for every seed of every soil, plus per-seed structural descriptors
that vary with the realization (number of macropore regions, coating-cell count,
dry-habitat fraction, matrix fraction) and the shape summary metrics below.

## 2. Averaging and confidence bands

Sand and Vertisol use different windows (`T_sand` vs `T_vertisol`), so resample each
`R(t)` onto a common **normalized fraction-of-window axis** (as Stage 11 did for
distinctness). Per soil, compute the **mean curve** and a **95% band** (mean ± 1.96·SE,
or the 2.5/97.5 percentiles across seeds). Plot both soils' mean curves with their
bands on the same normalized axis — the headline figure.

## 3. Per-seed shape metrics (for the tests)

Per seed, per soil, record: `R_peak`, `peak_frac` (time-to-peak / window),
`R_end/R_peak` (crash vs sustained), cumulative CO2, and late-third/early-third mean
ratio. Report mean ± SD per soil for each.

## 4. Significance — between-soil vs within-soil variability

The claim to test: **the between-soil difference far exceeds the seed-to-seed
(within-soil) spread.** Show it two ways (implement the stats in numpy — no scipy):

- **Metric-level:** for the key shape metrics (`peak_frac`, `R_end/R_peak`,
  late/early ratio), test Sand vs Vertisol across the two ensembles with a
  **Mann–Whitney U** (rank-based, robust) AND report an **effect size** (Cohen's d or
  rank-biserial). Report the two ensembles' distributions (box/violin) so the
  separation is visible.
- **Curve-level:** compute normalized-curve distances WITHIN each soil (all
  seed-to-seed pairs) and BETWEEN soils (all cross pairs). Show the **between-soil
  distances ≫ within-soil distances** (report both distributions; a permutation test
  or a bootstrap CI on the gap). Non-overlapping mean curves ± bands is the visual
  version of the same claim.

## 5. Do the simple features drive it? — the interpretation, backed by numbers

State the conclusion explicitly and support it with §4: because each soil's
structural class is held fixed while only the random layout changes, and each soil
robustly reproduces its characteristic shape (Sand burst-then-crash, Vertisol slow
rise) with between-soil separation far larger than seed noise, the respiration shape
is determined by the **structural features** (Sand: few concentrated coating-fed
regions, small matrix; Vertisol: many distributed regions, large matrix), not by the
particular realization.

**Optional causal check (recommended if runtime allows, ≤ a small sweep):** to show
the FEATURE itself is the driver (not just "Sand always bursts"), vary ONE structural
feature over a few values × a few seeds and show the shape metric moves with it —
e.g. sweep the coating fraction (or the correlation length that sets region count)
and show `peak_frac` / `R_end/R_peak` shift monotonically. Keep this bounded; it is
supporting evidence, not the main deliverable.

## 6. Deliverables

- `stage12_run.py`: the seed-ensemble loop (M seeds × 2 soils, n=128), the
  resampling/averaging, the numpy stats (Mann–Whitney, effect size, within/between
  distance distributions, bootstrap/permutation), and the figures.
- `results/stage12/`: the mean±band curves figure; per-seed metric distributions
  (box/violin) per soil; the within-vs-between curve-distance figure; a CSV of every
  seed's metrics + structural descriptors; a stats summary CSV (U, p, effect size per
  metric). Optional feature-sweep figure if done.
- `docs/stage_results/STAGE12_RESULTS.md`: the ensemble mean shapes, the significance
  results (with the actual U/p/effect-size numbers), the between≫within finding, and
  the explicit conclusion that the defined structural features drive the shapes.
  Update `MODEL_SPEC.md`, `LOGBOOK.md`, `report_prompt.md`.

## 7. Order of work

1. Seed-ensemble runner (M ≥ 30 seeds × Sand, Vertisol), frozen Stage-11 params.
2. Resample, average, confidence bands; per-seed metrics + structural descriptors.
3. Between-vs-within stats (Mann–Whitney + effect size; curve-distance distributions).
4. (Optional) the bounded feature-sweep causal check.
5. STOP; write `STAGE12_RESULTS.md` with the numbers and the conclusion.

Model and per-soil parameters stay frozen — this stage only adds seeds, averaging,
and statistics. Report the total runtime.
