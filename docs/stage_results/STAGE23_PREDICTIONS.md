# Stage 23 predictions (pre-run)

These predictions were written after Step 1 geometry and before any Stage 23 respiration run.
They use only the fixed scan geometries and chemistry constants. No laboratory respiration
measurements were opened or used.

## Geometry informing the predictions

| Structure soil | Porosity | Matrix | Solid | LT d10 (µm) | LT d50 (µm) | LT d90 (µm) | Largest pore cluster | Spans x/y/z | S_w | Habitat | Wet habitat fraction |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Sand | 0.320 ± 0.087 | 0.082 ± 0.000 | 0.598 ± 0.087 | 49.7 ± 8.0 | 101.4 ± 10.7 | 165.2 ± 21.0 | 0.977 ± 0.037 | 0.88/1.00/1.00 | 0.294 ± 0.051 | 0.320 ± 0.087 | 0.104 ± 0.018 |
| Vertisol | 0.164 ± 0.060 | 0.819 ± 0.059 | 0.017 ± 0.006 | 41.3 ± 17.1 | 93.9 ± 49.4 | 197.1 ± 86.7 | 0.710 ± 0.258 | 0.62/0.88/0.50 | 0.872 ± 0.072 | 0.164 ± 0.060 | 0.265 ± 0.134 |
| Loess | 0.246 ± 0.072 | 0.748 ± 0.073 | 0.006 ± 0.008 | 26.9 ± 15.2 | 85.6 ± 46.8 | 179.4 ± 106.3 | 0.874 ± 0.152 | 0.83/1.00/0.83 | 0.844 ± 0.094 | 0.246 ± 0.072 | 0.394 ± 0.255 |

Under the inherited matric rule, pore voxels with local thickness below 50 µm are wet and
diffusive; larger pores remain habitat through the drained-film diffusivity and exchange.
Loess has 0.246 resolved porosity and 0.394 wet fraction of habitat, higher than both Sand and Vertisol (Sand 0.104; Vertisol 0.265). The prediction is a drained-film/exchange-dominated response, with the wet-pore fraction limiting direct diffusion.

## Arm A ordering prediction

Loess is predicted to have an intermediate structure-only response between Sand and Vertisol
in peak timing and curve shape: its resolved porosity and pore LT d50 are intermediate, while
its wet habitat fraction and S_w are higher than both other soils. This may make its response
earlier or more active than the porosity ordering alone suggests. The exact peak height depends
on connectivity and habitat fraction. This is a geometry-based prediction, not a claim about
measured respiration. The Loess specimens are expected to show appreciable specimen-to-specimen
spread because LA and LB are downsampled labels while LC is a native 15 µm scan.

## Structure × chemistry cross

| Output | Predicted dominant factor | Pre-run reason |
|---|---|---|
| Peak height (R_peak/floor) | Chemistry | OM and B0 are changed together; the chemistry-specific OM:B0 ratio ranges from about 66.3 (Sand) to 94.6 (Vertisol), changing substrate availability per initial biomass. |
| Peak time | Structure | Pore LT distributions, connected clusters, wet habitat and matrix exchange distances alter substrate delivery times even when chemistry is held fixed. |
| Cumulative R to day 15 | Chemistry | OM_density spans 0.584–1.642; the largest expected direct effect is on total available substrate, with geometry modulating how much is accessed by day 15. |

These are directional predictions only; the fixed 3×3 cross will quantify main effects and
interaction without tuning or dropping crops.

## Arm D prediction

Filling the narrowest Sand pore voxels with matrix rather than carving fines from inert solid
is expected to reduce resolved porosity and the narrow-pore habitat fraction, while increasing
matrix volume and proximity for exchange. The change should broaden or delay the response and
may change the Sand-versus-Loess/Vertisol ordering; arm D is a model-choice check and will not
be adopted.
