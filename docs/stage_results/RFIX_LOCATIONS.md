# Growth-respiration formula audit

The scan covered repository `.py`, `.md`, and `.yaml` files for `(1 - Y)`,
`(1.0 - Y)`, `1 - yield`, `growth_R`, `growth_term`, `R_field`, and `gR`.
Line numbers below refer to the pre-fix source scan.

| File and location | What it computes or documents | Required action |
|---|---|---|
| `model.py:17,234` | Documents and computes growth-associated respiration from biomass growth. | Correct the computation; retain the old expression behind `legacy_growth_resp=True`. |
| `dual_porosity.py:1468,1477` | Computes total respiration and `growth_term_t` from biomass growth. | Correct both outputs; retain the old expression behind `legacy_growth_resp=True`. |
| `model3d.py:8,353,376` | Documents and computes growth respiration in the GPU fast path and generic path. | Correct both paths and document the substrate-based formula; retain the old expression behind `legacy_growth_resp=True`. |
| `stage10_run.py:109,253-254` | Reads the engine's growth/maintenance split and labels the plotted growth term as `(1-Y)*growth`. | Keep reading the engine output and update the outdated plot label. |
| `stage13_run.py:110`, `stage15_run.py:114,171`, `stage16_run.py:164`, `stage19_run.py:121` | Read `growth_term_t` from the dual-porosity engine for peak or growth-share summaries; they do not independently recompute respiration. | No formula change; corrected engine output flows through. |
| `MODEL_SPEC.md:107`, `README.md:163` | Describe the respiration observable or growth/maintenance outputs. | State the corrected equation explicitly. |
| `prompts/soil_respiration_prompt.md:74` | Historical prompt text documenting the correct substrate-based equation. | No computation; historical document is not a runtime source. |
| `LOGBOOK.md:873`, `MODEL_SPEC.md:855`, `README.md:547`, `docs/stage_results/STAGE10_RESULTS.md:177`, `report_prompt.md:489`, `prompts/soil_respiration_prompt_stage10_final.md:15,40` | Historical descriptions of growth-associated respiration or its share; no executable recomputation. | Preserve as historical records; distinguish new corrected results in the rfix report. |
| `metrics.py:94,114,304` | Summarizes or compares the supplied `R_t` arrays; it does not reconstruct growth respiration. | No formula change; use corrected `R_t` when analyzing new outputs. |
| `timemap.py:50-77` | Computes CUE from supplied respiration and uptake arrays; it does not reconstruct respiration. | No formula change; pass corrected respiration curves. |
| `stage23_analysis.py:71-80,254-260,721-742` | Loads stored `R`/`uptake_dt` arrays and derives normalized curves and physiology; stored-R-dependent summaries require recomputation from corrected curves. | Leave this historical analysis module unchanged; the new Stage 23 rfix analysis must recompute derived values. |
| `stage23_run.py:900-930`, `stage24_run.py` | Persists engine output and summarizes its `R`; no independent respiration equation was found. | No formula change; pass the legacy flag only in old-result reproduction gates. |
| Other `stage*_run.py` files | Search found no independent `(1-Y)` respiration calculation. The stage 10/13/15/16/19 consumers are listed above. | No formula change beyond the stage 10 plot label; reproduction checks must explicitly request legacy behavior. |

No separate respiration equation was found in `metrics.py`, `timemap.py`,
`stage23_analysis.py`, or the stage runners beyond the engine sites and the
stage 10 plot label listed above. The only post-hoc CUE calculation found is
`timemap.py`, which uses the supplied `R` trace.

## Stage impact and trace availability

| Stage | Headline depends on absolute R? | Stored growth/maintenance traces? | Correction status |
|---|---|---|---|
| 1 | Yes; respiration curves and integrals. | No paired component traces in the committed results. | Not corrected; do not rerun. |
| 2 | Yes; respiration curves and sweep metrics. | No paired component traces. | Not corrected; do not rerun. |
| 3 | Yes; respiration curves and comparisons. | No paired component traces. | Not corrected; do not rerun. |
| 3.5 | Yes; respiration curves and mechanism comparisons. | No paired component traces. | Not corrected; do not rerun. |
| 4 | Yes; comparative respiration outputs. | No paired component traces. | Not corrected; do not rerun. |
| 5 | Yes; R-based mechanism outputs. | No paired component traces. | Not corrected; do not rerun. |
| 6 | Yes; R-based mechanism outputs. | No paired component traces. | Not corrected; do not rerun. |
| 7 | Yes; R-based transport outputs. | No paired component traces. | Not corrected; do not rerun. |
| 8 | Yes; R-based depletion outputs. | No paired component traces. | Not corrected; do not rerun. |
| 9 | Yes; R-based shape and ranking outputs. | No paired component traces. | Not corrected; do not rerun. |
| 10 | Yes; R peak and growth-share outputs. | No paired component traces. | Not corrected; do not rerun. |
| 11 | Yes; R-based identity and headline outputs. | No paired component traces. | Not corrected; do not rerun. |
| 12 | Yes; ensemble R curves and significance metrics. | Only time-resampled R curves, without `growth_R`/`maint_R`. | Appendix Figure A1 curves were freshly rerun with corrected respiration; historical Stage 12 results remain legacy. |
| 13 | Yes; R-based diagnostics. | No paired component traces. | Not corrected; do not rerun. |
| 14 | Yes; R-based recruitment diagnostics. | No paired component traces. | Not corrected; do not rerun. |
| 15 | Yes; R-based diagnostics. | No paired component traces. | Not corrected; do not rerun. |
| 16 | Yes; includes Stage 11 R reproduction checks. | No committed paired component traces. | Reproduction calls now explicitly request legacy behavior; no stage rerun. |
| 17–18 | No Stage 17/18 run or results artifacts are present in this checkout. | None. | Not applicable. |
| 19 | Yes; R-based growth-share output. | No paired component traces. | Not corrected; do not rerun. |
| 20 | Yes; R-based film/transport outputs. | No paired component traces. | Not corrected; do not rerun. |
| 21 | Yes; R-based respiration outputs. | Yes; stored run arrays include both components. | Not corrected by this prompt's scope; do not rerun. |
| 22 equal-wetting (`results/stage22/`) | Yes; R-based sensitivity outputs. | Yes; stored run arrays include both components. | Not corrected; the requested Stage 22 target is the scan-structure stage below. |
| 22 scan-structure (`results/stage22_scan_structure/`) | Yes; normalized R curve distinctness and time-mapped physiology. | Yes; all 58 stored runs include both components. | Corrected from stored traces into `results/stage22_scan_structure_rfix/`; no reruns. |
| 23 | Yes; all respiration-based outputs. | Yes; 131 stored runs include both components. | Corrected from stored traces; six fresh GPU reproductions validate the transform. |
| 24 | Yes; variants × outputs derive from R. | Yes; 184 variant traces include both components. | Corrected post-hoc into `results/stage24_rfix/`; no new Stage 24 runs. |

Stages 1–20 have no stored per-step `growth_R`/`maint_R` arrays in this
checkout, so their historical absolute-R values cannot be post-hoc corrected
under the trace rule. Stage 21 and the separate equal-wetting Stage 22 archive
are technically correctable from their saved components, but remain legacy
because the prompt limits earlier-stage correction to Stage 22 scan-structure.
