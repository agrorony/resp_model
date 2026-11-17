# Fix — growth respiration is not carbon-closed

You are working in the `resp_model` repo on the GPU machine that ran Stages 22–23. Read the whole prompt first. Work in order.
Stop and report at any failed gate; do not improvise around it.

## 0. The problem (verified by reading the code and by post-hoc recomputation)

In every engine the growth-associated respiration is `(1 − Y) · growth`, where `growth` is new biomass per step, but substrate
uptake is `uptake = growth / Y`. Carbon taken up splits into biomass (`Y · uptake = growth`) and CO2 (`(1 − Y) · uptake`). The correct growth
respiration is therefore `(1 − Y) · uptake = (1 − Y) · growth / Y`, which is `1/Y` (2.5x at Y = 0.40) larger than the code. Maintenance
`m0 · B` is unaffected.

Locations (verify by grep, do not trust line numbers):
- `model.py`: `R_field = (1.0 - Y) * growth + m0 * B` (docstring line ~17 says `(1 - Y) * [substrate used for growth]`, which is the correct
  form; the code multiplied biomass growth instead of substrate used).
- `dual_porosity.py`: `R_field = (1.0 - Y) * growth + m0 * B` and `growth_term_t[t] = ((1.0 - Y) * growth).sum()`.
- `model3d.py`: both the GPU fast path (`gR = (1.0 - Y) * growth.sum()`) and the generic path (`gR = ((1.0 - Y) * growth).sum()`),
  plus the docstring `R = (1 - Y) growth + m0 B`.
- `MODEL_SPEC.md` / `README.md`: wherever the R formula is written.

R is an output only. It does not enter dB, S, OM or uptake. The mass-balance check only covers substrate (source − uptake + clip), so it never saw this.
A symptom already visible in Stage 23: reported CUE = 1 − cumR/cumUptake is 0.67–0.69, above Y = 0.40, which is impossible with maintenance respiration.

Post-hoc check already done outside the repo (scratch, not committed): using `R_new = growth_R / Y + maint_R` on the stored Stage 23 traces, all
registered distinctness verdicts and Arm A/B primary-mapping soil rankings stay the same, R scales by about 2.1–2.3x, CUE becomes 0.31–0.33.
Do not rely on that scratch work. Redo it here under the rules below.

## 1. Hard rules

1. Biology unchanged: `configs/biology.yaml` byte-identical to main. No constant is tuned. No new variant, no new structure, no new chemistry.
2. Only the respiration **output** formula changes. B, S, OM, uptake, S_total, mass balance and all habitat/delivery fields must stay identical to the
   stored runs (max relative difference 0, or < 1e-12 if the GPU makes 0 impossible).
3. Keep the old behaviour available: add a keyword `legacy_growth_resp: bool = False` to the engines (`model.py`, `dual_porosity.py`,
   `simulate3d`). `True` reproduces the old formula exactly. Default is the corrected formula. Every reproduction gate of earlier stages runs with
   `legacy_growth_resp=True`.
4. Do not overwrite any stored result. Corrected outputs go to new folders (`results/stage23_rfix/`, etc.). Existing `results/stageNN/` stay as they are and
   are labelled as legacy in the report.
5. Do not touch `m_s`, the starvation term or any other part of the biology (see Step 6, report only).
6. Do not open lab/measurement files, except for Step 5 if the Figure A1 script already reads them. Then use that script's own inputs, unchanged.
7. Never reuse a stale file as if it were fresh. Every corrected curve in a gate comes from a run executed in this session (check timestamps).
8. Branch `agents/fix-growth-respiration`. Commit after each step. Push and open a PR at the end. Do not merge.

## 2. Step 0 — find every place

Grep all `.py`, `.md`, `.yaml` for `(1 - Y)`, `(1.0 - Y)`, `1 - yield`, `growth_R`, `growth_term`, `R_field`, `gR`. Write the list to
`docs/stage_results/RFIX_LOCATIONS.md` (file, line, what it computes, change needed or not). Include the stage runners
(`stage*_run.py`) and `metrics.py`, `timemap.py`, `stage23_analysis.py`: any place that recomputes R, growth share or CUE from stored traces.
Stop if any place uses the formula in a way not covered by this prompt, and report it.

## 3. Step 1 — code change

1. In the three engines, replace growth respiration by `(1 − Y) · uptake` (equivalently `(1 − Y) · growth / Y`), behind `legacy_growth_resp`.
   Keep `m0 · B` as is. Keep the output names (`R`, `growth_R`, `maint_R`) so analysis code keeps working.
2. Update docstrings, `MODEL_SPEC.md`, `README.md` so the formula is written once, explicitly: `R = (1 − Y) · uptake + m0 · B`, with
   `uptake = growth / Y`.
3. Add a carbon-closure diagnostic to `simulate3d` output (new fields, no change to existing ones): per run,
   `closure_growth = Σ uptake·dt − Σ growth·dt − Σ growth_R·dt` (must be 0 within 1e-10 relative to Σ uptake·dt under the corrected formula; the legacy
   formula must fail it, and a test must show that).
4. New tests in `test_fix_growth_resp.py`:
   - toy case, one voxel, constant S: corrected R equals `(1 − Y)/Y · growth + m0 B`;
   - legacy flag reproduces the old number exactly;
   - closure diagnostic passes for corrected and fails for legacy.
5. Run the existing tests. Commit.

## 4. Step 2 — reproduction and consistency gates (fresh GPU runs)

For Stage 23 arm A and arm B, crops `SA_c0`, `VA_c0`, `LC_c0` (6 runs, same constants, T 3000, K_density 0.16445 as Stage 23):
1. Run with `legacy_growth_resp=True`. All stored arrays (R, B, growth_R, maint_R, S_total, uptake_dt, …) must match the stored Stage 23
   arrays (max relative difference < 1e-8). If not, stop.
2. Run with the default (corrected). Then:
   - `B`, `S_total`, `uptake_dt`, `fuel`, `maint_R`, and the delivery fields identical to the legacy run;
   - `R_new == growth_R_legacy / Y + maint_R_legacy` (max relative difference < 1e-10);
   - carbon closure within 1e-10;
   - mass balance within 1e-10.
If any fails, stop.

Because R does not feed back, the other 125 Stage 23 runs are NOT rerun. Their corrected curves come from the stored traces with
`growth_R / Y + maint_R`, justified by the gate above. State this in the report.

## 5. Step 3 — corrected Stage 23 analysis from stored traces

Create `stage23_rfix_analysis.py` (new file, imports `stage23_analysis` helpers, does not edit them). It reads the stored Stage 23 traces from
`results/stage23/runs/` and builds the corrected R as in Step 2, then reruns steps 3, 4, 5 and 6 of the Stage 23 analysis into `results/stage23_rfix/`.

Important: several fields in the stored run summaries are derived from R and are **not** in the arrays. They must be recomputed from the corrected
curves, not copied: `R_peak_over_floor`, `peak_frac` (peak fraction of T), `R_end_over_peak`, `growth_share_at_peak`, `R_peak`.
Also recompute CUE (expect below Y = 0.40), qCO2, cumulative R per volume and per mass, P1–P4 means and shape classes.
The floor stays `m0 · B0`.

Deliver `docs/stage_results/STAGE23_RFIX.md` with, each as old (legacy) → new (corrected):
1. Arm A distinctness per pair (distance at 0–1500 and 0–3000, criteria, verdict), peak height and peak time order.
2. Arm B distinctness, cumulative R per mass rank, P4 rate rank, qCO2 rank, at low / primary / high mapping.
3. Shape classes per soil, arm and mapping.
4. Two-way decomposition shares (structure / chemistry / interaction / residual) for every output.
5. Soil summary table (R_peak/floor, peak fraction, growth share, CUE, qCO2, cumulative R).
6. A list **"what flipped"**: every verdict, rank or class that changed, with the numbers, and whether the old margin was already within noise.
7. A list **"numbers to update in the appendix"**: every R-based number the Stage 23 results file or the appendix text quotes, with the corrected value.

## 6. Step 4 — earlier stages and Stage 24

- Stages 22 and earlier: do not rerun. In `docs/stage_results/RFIX_LOCATIONS.md` add a table listing each stage, whether its headline conclusion depends on the
  absolute R scale, and whether it can be corrected from stored traces. Correct from stored traces only Stage 22 (the scan-structure stage) if its traces
  contain `growth_R` and `maint_R`; otherwise say so.
- Stage 24 (branch `agents/stage24-matrix-sensitivity`, results in `results/stage24/`): if it is present in this checkout or its branch can be fetched, apply the same
  post-hoc correction to its variants × outputs table into `results/stage24_rfix/` and add the same old → new comparison, with special attention to
  the "same / different from base" columns. If it is not available, say so and stop for this item. Do not merge the Stage 24 PR.

## 7. Step 5 — Figure A1 data

Find the script that produced Appendix Figure A1 (look in `results/appendix_figure/` and the repo root). Run it unchanged except that R comes from the corrected traces,
writing to `results/appendix_figure_rfix/`. Report which panels change visibly. If its comparison panels normalise R (for example to the first period),
say whether the normalisation cancels the 1/Y factor. If the script reads lab files, use its own inputs exactly as before and do nothing else with them.
If it needs a path that is missing, stop and report the path.

## 8. Step 6 — report only: starvation loss

`m_s = 0.10`. The term `m_s (1 − f) B` removes biomass in `dB` but its carbon goes neither to CO2 nor back to S. Measure and report, per soil in the
Stage 23 arm A primary runs (from stored traces where possible, otherwise from the six fresh runs): the cumulative carbon removed by this term, as a fraction of cumulative uptake and of
cumulative R. Do not change the model. State plainly that this is a second non-closed carbon flow and that the decision (starved biomass to CO2, back to S, or leave as is)
belongs to Rony.

## 9. Write-up

1. `docs/stage_results/RFIX_REPORT.md`: the problem, the code diff summary, the gates (numbers), the Stage 23 old → new tables, "what flipped",
   "numbers to update in the appendix", Stage 22 / 24 status, Figure A1 status, the starvation-loss measurement, limitations.
2. **Verdict in plain words:** (a) which Stage 23 conclusions survive unchanged, (b) which numbers or classes change, (c) whether anything in the appendix text becomes false.
3. LOGBOOK headline, run-time table, commit, push, PR.

## 10. Stop rules

Stop and report if: a stored trace lacks `growth_R` or `maint_R`; the legacy gate fails; B, S, uptake or maintenance R differ between legacy and corrected runs;
carbon closure or mass balance exceeds 1e-10; `biology.yaml` differs from main; a formula location is found that this prompt does not cover.

Final message: the gate numbers, the "what flipped" list (or "nothing flipped"), the corrected soil summary table, the starvation-loss fractions, and the PR link.
