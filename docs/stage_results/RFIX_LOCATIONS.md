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
