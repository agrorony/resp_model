# Variant C OM Refactor Report

## Setup

This run applies soil-specific initial OM fractions in Variant C at a single condition: theta=0.6, n=128.
OM initialization is applied only to porous matrix cells (solid phase and macropore cavity cells receive zero OM).
Configured OM fractions:
- Sand: 2.0%
- Loess: 3.5%
- Vertisol: 5.0%

## Headline Metrics (theta=0.6, n=128)

| Soil | outcome | R_peak | mean_S_habitat_final | wet_habitat_fraction |
|---|---:|---:|---:|---:|
| Sand | declining | 0.400000 | 0.001008 | 0.4263 |
| Loess | declining | 0.400000 | 0.000928 | 0.1451 |
| Vertisol | declining | 0.400000 | 0.001460 | 0.1384 |

## Dynamic Analysis

The refactor combines three transport controls in Variant C: soil-specific OM loading, internal porous-matrix diffusion, and wet-matrix-driven network wicking toward macropore habitats.
Observed behavior: Sand (declining, shape=falling, R_peak=0.400000, mean_S_habitat_final=0.001008, wet_habitat_fraction=0.4263) ; Loess (declining, shape=falling, R_peak=0.400000, mean_S_habitat_final=0.000928, wet_habitat_fraction=0.1451) ; Vertisol (declining, shape=falling, R_peak=0.400000, mean_S_habitat_final=0.001460, wet_habitat_fraction=0.1384).
Interpretation: the custom OM loading increases available dissolved substrate in matrix storage, while matrix diffusion and wicking govern how much of that substrate actually reaches wet habitat cells.
Under this single condition, starvation is reduced where wet habitat connectivity and transfer are high; soils with lower wet habitat exposure continue to show declining/falling behavior even with larger OM loading.

## Pairwise Distinctness (theta=0.6, n=128)

| Pair | distance | pass (>0.3) |
|---|---:|---:|
| Loess vs Sand | 0.003174 | no |
| Loess vs Vertisol | 0.001673 | no |
| Sand vs Vertisol | 0.001895 | no |

Minimum pairwise distance: 0.001673 (threshold=0.300).
Distinctness success status: FAIL for the >0.3 criterion.

## Output Files

- variantC_OM_refactor_rates.png
- variantC_OM_refactor_cumulative.png
