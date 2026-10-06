# Claude Code task: Stage 3.5 — literature-informed PSDs (Loess / Sand / Vertisol)

Extend the EXISTING numpy-only project (Stages 1–3 done). Do not rewrite it.

## 0. Why this stage exists

The Stage-3 image-derived PSDs were poor: 30 µm minimum resolvable diameter,
imperfect Otsu segmentation, and (per the sample with the highest true porosity,
the Vertisol) many pores were missed entirely. So Stage 3 could only tell the fine
soil apart from the two coarse ones. Stage 3.5 replaces the measured PSDs with
**literature-informed parametric PSDs** for the three real soil types in my
experiment, and tests whether the model can reproduce the **measured respiration
pattern**:

- **Vertisol — rising** over time (my data: ~104 → 142 → 158 → 173)
- **Sand — falling** (~116 → 94 → 90 → 74)
- **Loess — low and erratic** (~68 → 25 → 50 → 11)
  (mg CO2 kg⁻¹ day⁻¹, periods P1–P4 of 3/4/3/5 days.)

This experimental pattern is the validation target for this stage — but only
**qualitatively**: I do NOT need the model to match my numbers, only to reproduce
the *ranking and shapes* (Vertisol rising, Sand falling, Loess low/erratic). The
PSDs likewise need only be representative soil types, not my exact samples.

**Crucial:** my experiment was at **constant moisture** — the soil cups sat in
closed jars with water in the bottom, so humidity stayed high and constant for the
whole run. The P1–P4 changes are substrate/structure dynamics at FIXED saturation,
NOT drying. So this stage runs at constant moisture (see §4), not a θ sweep.

## 1. Invariants (unchanged)

- `biology.yaml` and the pore→parameter mapping frozen; numpy + pyyaml +
  matplotlib only; deterministic; a single **constant moisture** held over time
  (see §4) — matches the closed-jar experiment, which did not dry out.
- Keep the physiological K window **K=0 below 30 µm, plateau 30–150 µm, decay
  above 150 µm** — it is literature-supported (30–150 µm is repeatedly described
  as the optimal microbial-colony pore range). Diameters, not radii.

## 2. Literature-informed PSDs (replace the measured PSDs)

Implement each soil's PSD as a **mixture of lognormal modes** over pore diameter
(µm), feeding the SAME empirical-quantile pore-field pipeline built in Stage 3
(`generate_diameter_field_from_psd`). Add a `pore.mode: "psd_parametric"` path
alongside the existing `"psd"` (real data) and `"lognormal"` paths, so nothing
breaks. Use these as the starting, clearly-labeled **literature-informed
representative** distributions (mode = peak diameter µm, GSD = geometric SD,
w = volume weight; keep them in a config so they're easy to edit):

| Soil | mode 1 (clay/matrix) | mode 2 (intra/meso) | mode 3 (macro/channel) | total porosity |
|---|---|---|---|---|
| **Sand** | 3 µm, GSD 2.0, w 0.15 | — | 120 µm, GSD 2.2, w 0.85 | ~0.38 |
| **Loess** | 0.5 µm, GSD 2.2, w 0.20 | 12 µm, GSD 2.0, w 0.50 | 90 µm, GSD 2.2, w 0.30 | ~0.47 |
| **Vertisol** | 0.01 µm, GSD 2.5, w 0.35 | 0.35 µm, GSD 2.2, w 0.30 | 80 µm, GSD 2.0, w 0.35 | ~0.55 |

Rationale (encode as comments/citations in the config):
- **Sand**: macropore-dominated (>50 µm), narrow, well-connected, little fine
  volume → low water retention.
- **Loess**: bimodal — a silt-matrix micropore mode (~10–20 µm) plus a
  macropore/channel mode (wormholes, fissures, collapse >100 µm); meso+macropores
  are ~65–70% of pore volume.
- **Vertisol**: multimodal — a clay micropore mode (~0.009–0.012 µm), an
  intra-aggregate mode (~0.3–0.4 µm), and an inter-aggregate / shrink-swell-crack
  macropore mode (~60–100 µm). Highest total porosity, high water retention (most
  volume in fine pores), but a real macropore mode that lands **inside the
  30–150 µm habitat window** — that macropore mode is where its microbial habitat
  is.

Build each soil's empirical CDF by sampling the mixture on a fine diameter grid,
then reuse the existing quantile-mapping to lay diameters on the field. Same
shared spatial settings, different seeds per soil (as Stage 3).

## 3. Fix the D(d) fusion from Stage 3

Stage 3 fused the two coarse soils because `D(d)` saturated (`d_ref=90 µm`), so all
large pores looked alike to transport. Make `D(d)` keep discriminating across the
real diameter range: raise `d_ref` (try ~250 µm) and/or reduce the saturation so
that, at least across ~30–300 µm, `D` still rises appreciably with diameter.
Verify the numerical-stability assertion still passes and that two coarse-but-
different soils no longer collapse to the same curve. Keep it fixed/shared across
soils.

## 4. Constant moisture (the experiment was NOT a drying run)

Moisture was held constant in the experiment (§0). So do NOT sweep θ as the main
analysis. Instead:

- Run each soil at a **single constant moisture level held for the whole
  simulation** — a moist-but-aerated value (start with θ ≈ 0.6 water-filled, the
  SAME fixed value for all soils; respiration needs air, so not fully saturated).
  Small pores fill first (capillarity), as before, and the wet set stays fixed in
  time.
- The observable is the **temporal shape of R(t)** at that fixed θ, per soil —
  directly comparable to my P1–P4 trend.

Expected mechanism at constant moisture (this is the substrate-access story, not a
drying story):
- **Vertisol**: OM is protected in its abundant fine pores, which stay wet the
  whole run; substrate is released and diffuses slowly to the sparse habitat
  macropores → access keeps building → **rising** respiration.
- **Sand**: big, well-connected wet pores give immediate access to substrate → an
  early burst that then depletes → **falling**.
- **Loess**: intermediate/bimodal, smaller habitable fraction → lower, more
  **erratic** respiration.

Optionally, as a SECONDARY robustness check only, run 2–3 nearby fixed θ (e.g.
0.5 / 0.6 / 0.7) to confirm the ranking/shape is not knife-edge on the exact
value — but the headline is the single constant-θ time series vs the experiment.

## 5. Use the real soil names and compare to the experiment

- Relabel the three soils **Loess / Sand / Vertisol** throughout (configs, figures,
  CSVs, report). Map: Sand↔coarse, Loess↔intermediate/bimodal, Vertisol↔the
  high-porosity multimodal soil.
- Add an explicit comparison to the measured pattern (§0): does the model give
  Vertisol rising, Sand falling, Loess low/erratic? Plot the model respiration
  curves next to the experimental period means (a normalized/qualitative overlay
  is fine — the point is the *shape/ranking*, not absolute mg CO2). State plainly
  which of the three the model reproduces and which it doesn't.

## 6. Keep the Stage-3 machinery

Keep: the isolation-ratio diagnostic (computed at the fixed θ); the
emergent-distinctness metric with `require_nontrivial=False` (die-off = valid
distinct outcome); grid at n=128 (n=40 optional). Budget cap `MAX_ITERATIONS=12`;
this is a defined comparison, not a search.

## 7. Deliverables

- `pore_field.py` / a small `psd_parametric.py` — the lognormal-mixture PSD builder
  and `psd_parametric` pore mode.
- `configs/psd_literature.yaml` — the per-soil mixture parameters (editable) with
  citations in comments.
- `configs/mapping.yaml` — updated `D_d_ref` (D(d) fusion fix, §3).
- `configs/soil_{loess,sand,vertisol}.yaml` — renamed, `mode: psd_parametric`.
- A Stage-3.5 runner — constant-θ time-evolution run per soil, the experiment
  overlay (§5), and the optional 2–3-value θ-sensitivity check (§4).
- `results/stage35_*` — R(t) per soil at the fixed θ, experiment-overlay figure,
  water maps, isolation/distinctness, the optional θ-sensitivity panel, and the
  loaded PSDs plotted.
- `STAGE3_5_RESULTS.md`, `LOGBOOK.md`, `MODEL_SPEC.md`, `report_prompt.md` updated.

## 8. Order of work

1. Parametric PSD builder + `psd_literature.yaml`; plot the three PSDs; sanity-check
   habitable-fraction (30–150 µm) per soil.
2. D(d) fusion fix (§3); confirm two different coarse PSDs now separate.
3. Rename to Loess/Sand/Vertisol; constant-θ time-evolution run at n=128.
4. Experiment overlay (temporal shape) + optional θ-sensitivity; isolation/
   distinctness at the fixed θ.
5. Checkpoint: STOP; write `STAGE3_5_RESULTS.md`, update the rest.

## 9. Honest framing

Answer plainly: with literature PSDs, does the model reproduce the measured
ranking and shapes (Vertisol rising, Sand falling, Loess low/erratic)? If it
reproduces some and not others, say which and why (via water retention, habitable
fraction, and isolation). These PSDs are literature-informed *representative*
distributions, not measurements of my specific samples — state that clearly.
Total porosity differs between these soils (Vertisol highest) but the current grid
treats every cell as a pore and does not represent bulk porosity — note this as
the main remaining limitation and a candidate next stage.
