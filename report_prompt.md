# Report-generation prompt (run this later; do not run now)

You are writing a structured scientific report, in English, suitable for a
course appendix, about a soil-structure microbiome respiration simulation.

Read these inputs from the project folder before writing anything:

1. `MODEL_SPEC.md` — the model equations, frozen biology constants,
   structure knobs, and the reasoning for why each soil's structure produces
   its target respiration shape.
2. `LOGBOOK.md` — the complete, honest iteration-by-iteration record of the
   structure search: what was tried, the resulting metrics, and what was
   changed and why.
3. Everything in `results/`:
   - `respiration_curves.png` — R(t) for soils A, B, C on one axis.
   - `cumulative_co2.png` — cumulative CO2 for soils A, B, C.
   - `structure_maps_A.png`, `structure_maps_B.png`, `structure_maps_C.png`
     — D, K, and OM maps per soil.
   - `soil_A_timeseries.csv`, `soil_B_timeseries.csv`, `soil_C_timeseries.csv`
     — raw `time, R_t, cumulative_co2` per soil.

Write `REPORT.md` (or `REPORT.pdf` if a PDF pipeline is available) with the
following sections:

1. **Motivation** — why spatial structure alone, with identical biology,
   can produce qualitatively different respiration dynamics; the practical
   relevance to real soils (pore structure, substrate accessibility).
2. **Model and equations** — reproduce the equations from `MODEL_SPEC.md`
   (biomass, substrate, organic matter, respiration observable), explained
   in prose, with the frozen constants tabulated.
3. **The three soils and their structure** — for each soil, describe the
   diffusion rate, carrying-capacity layout, and organic-matter placement,
   embedding the corresponding `structure_maps_*.png` figure.
4. **Results** — embed `respiration_curves.png` and `cumulative_co2.png`;
   report the quantitative metrics (e, l, peak, roughness) for each soil
   from `LOGBOOK.md`'s final iteration.
5. **Interpretation** — explain, mechanistically, why each structure
   produces its shape (diffusion-front arrival time relative to the
   observation window; OM depletion; staggered multi-colony peaks for the
   flat/bumpy soil). Tie this explicitly back to the equations, not just the
   qualitative description.
6. **Limitations** — no calibration to real data; uniform D per soil (no
   D-field); Gaussian-only OM/K patch shapes; the roughness threshold in the
   flat classifier was picked empirically, not derived; deterministic
   dynamics with only a fixed-seed layout step.
7. **Budget note** — state plainly, from `LOGBOOK.md`, how many iterations
   the structure search actually took and whether all three target shapes
   were achieved. If the budget ran out before success, say so honestly and
   describe the best configuration reached.

Keep the tone scientific and precise. Do not invent results not present in
`LOGBOOK.md` or `results/`.
