# Growth-respiration correction report

## Executive summary

The model previously calculated growth-associated respiration from biomass
production, `(1 - Y) * growth`, instead of substrate uptake,
`(1 - Y) * uptake` where `uptake = growth / Y`. The corrected output is
`R = (1 - Y) * uptake + m0 * B`. The old expression remains available with
`legacy_growth_resp=True`; the corrected expression is the default.

Respiration is an output and does not feed back into biomass, substrate,
uptake, or delivery. Six fresh Stage 23 GPU reproduction gates passed. The
other Stage 23 traces, all Stage 22 scan-structure traces, and all Stage 24
variant traces were corrected post-hoc into new result folders. Existing
`results/stageNN/` archives and `configs/biology.yaml` were not changed.

## Code and validation

- Corrected `model.py`, `dual_porosity.py`, and both `model3d.py` execution
  paths; retained exact legacy behavior behind `legacy_growth_resp`.
- Added `closure_growth` and `closure_growth_rel` to the 3D output, and toy
  tests covering the corrected formula, exact legacy result, and closure.
- Updated the formula documentation and Stage 10 plot label; historical
  reproduction callers explicitly request legacy behavior.
- `configs/biology.yaml` remains byte-identical to `main`.

All six Arm A/B × SA_c0/VA_c0/LC_c0 gates passed:

| Gate check | Maximum observed |
|---|---:|
| Legacy output vs stored Stage 23 trace, relative difference | 5.34e-16 |
| Corrected formula vs `growth_R_legacy / Y + maint_R`, relative difference | 5.36e-16 |
| Corrected vs legacy state/delivery arrays | 0 |
| Corrected carbon-closure relative error | 1.05e-15 |
| Mass-balance relative error | 8.04e-15 |
| Gate runtime | 427.7 s |

The legacy formula fails the new growth-closure diagnostic as expected
(relative residual 0.36); the corrected expression passes. Fresh corrected
and legacy `B_final`/`S_final` and all compared state, uptake, maintenance,
and delivery arrays are identical. The 125 remaining Stage 23 traces were
not rerun; their corrected curves use `growth_R / Y + maint_R`, justified by
these gates.

## Stage 23 results

The full old → corrected tables—including all distinctness criteria,
rankings, shape classes, decomposition shares, and per-run means ± SD—are in
[`STAGE23_RFIX.md`](./STAGE23_RFIX.md). Key corrected primary-mapping soil
means are below (each cell is legacy → corrected):

| Arm | Soil | R peak / floor | Peak fraction | Growth share at peak | CUE | qCO2/day | Cumulative R / mass |
|---|---|---:|---:|---:|---:|---:|---:|
| A | Sand | 15.6424 → 36.3379 | 0.0430143 → 0.0420974 | 0.880487 → 0.949828 | 0.656604 → 0.296604 | 0.135779 → 0.261448 | 0.0753119 → 0.154308 |
| A | Vertisol | 7.01735 → 15.2917 | 0.144923 → 0.132961 | 0.782679 → 0.903603 | 0.690865 → 0.330865 | 0.213495 → 0.455737 | 0.0425638 → 0.0921454 |
| A | Loess | 11.5518 → 25.6267 | 0.115928 → 0.106758 | 0.802166 → 0.915048 | 0.690755 → 0.330755 | 0.197495 → 0.415738 | 0.0792245 → 0.171471 |
| B | Sand | 17.2056 → 39.8524 | 0.0533511 → 0.0521007 | 0.875602 → 0.947943 | 0.667279 → 0.307279 | 0.145177 → 0.284942 | 0.0580405 → 0.120867 |
| B | Vertisol | 8.1588 → 17.7326 | 0.140047 → 0.124875 | 0.777325 → 0.902429 | 0.689933 → 0.329933 | 0.203885 → 0.431711 | 0.0699702 → 0.151228 |
| B | Loess | 10.0723 → 22.278 | 0.124152 → 0.115511 | 0.799281 → 0.912842 | 0.691118 → 0.331118 | 0.202786 → 0.428965 | 0.0620718 → 0.13452 |

All corrected primary CUE means are below `Y = 0.40`, as required. The
Arm A Sand–Vertisol registered distinctness remains True; S–L and V–L remain
False. Arm B pairwise distinctness verdicts also remain unchanged. Arm A/B
primary-mapping peak-height/time orders, Arm B primary rankings, and primary
shape classes survive unchanged.

### What flipped

- Arm B high-mapping cumulative-R/mass and P4-rate rankings each change from
  `V > S > L` to `V > L > S`. The old Sand–Loess gaps were within one
  combined crop SD.
- Arm A high Sand and Arm B high Vertisol shape classes change from
  `rise-then-fall` to `fall`. In both cases the old P1/P2 ordering was within
  one combined crop SD.
- The Arm A Stage 21 pore-fill ordering at 3000 steps changes from True to
  False: the base closest pair changes S–V → S–L, while the pore-fill
  closest pair remains S–L. No registered noise interval exists for this
  single mean-curve ordering.
- Stage 24 variant-specific flips and all changed “same/different from base”
  comparisons are listed in [`STAGE24_RFIX.md`](./STAGE24_RFIX.md).

### Appendix updates and text that is no longer accurate

Replace the Stage 23 Step 2/4 R peak/floor, growth-share, peak-fraction,
end/peak, cumulative-R, qCO2, and decomposition values; the Step 5 Arm D and
domain R summaries; and the Step 6 P1–P4, CUE, qCO2, cumulative-R, class,
rank, peak-order, and curve-distance values. The corrected values are
tabulated in [`STAGE23_RFIX.md`](./STAGE23_RFIX.md). Biomass, habitat,
delivery, and substrate-uptake results are unchanged.

The old Stage 23 CUE range (about 0.66–0.69) is not physically valid for the
stated yield: corrected primary means are 0.2966–0.3311. The Stage 23
high-mapping rank/shape statements above also change. The pore-fill
closest-pair statement at 3000 steps must be revised. Stage 22 and Stage 24
historical R-based values likewise need replacement by their corrected
reports. Other historical stage result folders are explicitly legacy where
their outputs cannot be reconstructed from stored component traces.

## Stage 22 and Stage 24

- **Stage 22 scan-structure:** all 58 stored traces were post-hoc corrected;
  no simulations were rerun. Sand–Vertisol mean distances change
  0.554377 → 0.706310 at 0–1500 and 0.203087 → 0.282005 at 0–3000. The
  registered distinctness verdict remains True. Thirty-two existing
  sensitivity traces were also corrected. Step 4 comparisons to Stage 21
  synthetic controls remain mixed-era because those controls are unchanged.
  Details: [`STAGE22_RFIX.md`](./STAGE22_RFIX.md).
- **Stage 24:** all 184 variant traces and the variants × outputs comparison
  were corrected post-hoc; no simulations were rerun and no merge was done.
  S16 Sand–Vertisol distinctness changes False → True. Other variant class
  and base-comparison changes are detailed in
  [`STAGE24_RFIX.md`](./STAGE24_RFIX.md).
- Stages 1–20 lack the paired per-step growth/maintenance traces needed for
  post-hoc correction. Stage 21 and the separate equal-wetting Stage 22
  archive remain legacy and were not in scope. See
  [`RFIX_LOCATIONS.md`](./RFIX_LOCATIONS.md) for the complete stage inventory.

## Figure A1

The Figure A1 script was run for 80 fresh simulations (40 per soil) using
corrected respiration; no lab/measurement files were read. Both Sand and
Vertisol curves change visibly. Sand mean peak changes 32.3747 → 76.0995;
Vertisol changes 2.46457 → 7.64190. The curves are time-resampled but not
divided by their initial R, so the 1/Y correction does not cancel.
Corrected peak-fraction significance is `p = 2.66e-15`. The run took
1701.25 s. Outputs and the run manifest are in
[`results/appendix_figure_rfix/`](../../results/appendix_figure_rfix/).

## Starvation loss (report only)

The term `m_s (1 - f) B`, with `m_s = 0.10`, removes biomass whose carbon is
neither emitted as CO2 nor returned to substrate. Stored-trace estimates for
Stage 23 Arm A primary runs (ranges include the bounded final step) are:

| Soil | Fraction of cumulative uptake | Fraction of corrected cumulative R |
|---|---:|---:|
| Sand | 23.231–23.245% | 30.828–30.847% |
| Vertisol | 33.326–33.333% | 49.144–49.154% |
| Loess | 32.555–32.560% | 47.639–47.646% |

The six fresh gate runs provide exact full-window c0 estimates: Sand
23.857% of uptake / 31.781% of corrected R; Vertisol 33.220% / 48.999%;
Loess 33.451% / 49.124%. This is a second non-closed carbon flow. Whether
starved biomass should become CO2, return to S, or remain as implemented is
a biological decision for Rony; the model was not changed for it.

## Limitations and verdict

The correction changes respiration outputs only; it does not repair the
separate starvation carbon flow. Earlier stages without paired component
traces, Stage 21 controls, and the equal-wetting Stage 22 archive remain
legacy. The Stage 22 Step 4 source/resolution comparison is mixed-era.

In plain words: the key Stage 23 primary distinctness, ordering, ranking,
and shape conclusions survive, but the absolute respiration scale roughly
doubles, CUE drops below yield as carbon closure requires, and two
high-mapping shape/rank conclusions plus the pore-fill ordering change.
Some Appendix text is therefore false or numerically obsolete—especially
the old CUE claim, the listed R-based values, and the affected high-mapping
and pore-fill statements—and should be updated using the linked detailed
reports.
