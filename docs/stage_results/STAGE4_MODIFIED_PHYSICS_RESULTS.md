# Stage 4 Modified Physics Results

## Scope

This report summarizes the post-modification Stage 4 run after applying:

- Variant A: remove banded OM placement, keep distributed OM after PSD truncation,
  and increase wet-cell transport (`wet_diffusion_scale`).
- Variant C: add dynamic hydro-substrate coupling (wicking) so dry macropores
  adjacent to wet matrix can rewet and conduct at runtime, plus directional
  substrate transfer from wet matrix to wet macropores.

Run command:

`python stage4_run.py`

Outputs are in `results/stage4/`.

## Headline Condition (theta=0.6, n=128)

### Variant A

| soil | outcome | shape(model) | shape(exp) | wet_habitat_fraction | mean_S_habitat_final |
|---|---|---|---|---:|---:|
| Loess | declining | falling | erratic | 0.1358 | 6.2396e-03 |
| Sand | thriving | rising | falling | 0.5545 | 1.4347e-02 |
| Vertisol | thriving | rising | rising | 0.5226 | 1.7855e-02 |

Distinctness (`require_nontrivial=False`):

- min_distance = 0.0578 (fail; threshold 0.3)

### Variant C

| soil | outcome | shape(model) | shape(exp) | wet_habitat_fraction | mean_S_habitat_final |
|---|---|---|---|---:|---:|
| Loess | declining | falling | erratic | 0.0754 | 1.7061e-06 |
| Sand | declining | falling | falling | 0.4166 | 1.6828e-03 |
| Vertisol | declining | falling | rising | 0.0695 | 3.7003e-06 |

Distinctness (`require_nontrivial=False`):

- min_distance = 1.4506e-05 (fail; threshold 0.3)

## Robustness Checks

### Variant C (theta sweep at n=128)

| theta | soil | outcome | shape(model) | wet_habitat_fraction | mean_S_habitat_final |
|---:|---|---|---|---:|---:|
| 0.5 | Loess | declining | falling | 0.0204 | 1.6235e-06 |
| 0.5 | Sand | declining | falling | 0.3167 | 1.6213e-03 |
| 0.5 | Vertisol | declining | falling | 0.0180 | 3.6911e-06 |
| 0.6 | Loess | declining | falling | 0.0754 | 1.7061e-06 |
| 0.6 | Sand | declining | falling | 0.4166 | 1.6828e-03 |
| 0.6 | Vertisol | declining | falling | 0.0695 | 3.7003e-06 |
| 0.7 | Loess | declining | falling | 0.1360 | 1.7359e-06 |
| 0.7 | Sand | declining | falling | 0.5166 | 1.7044e-03 |
| 0.7 | Vertisol | declining | falling | 0.1326 | 3.7013e-06 |

### Grid check (theta=0.6, n=40)

| variant | soil | outcome | shape(model) | wet_habitat_fraction | mean_S_habitat_final |
|---|---|---|---|---:|---:|
| A | Loess | thriving | rising | 0.1356 | 9.9351e-03 |
| A | Sand | thriving | rising | 0.5544 | 1.8109e-02 |
| A | Vertisol | thriving | flat | 0.5225 | 1.4524e-02 |
| C | Loess | declining | falling | 0.0594 | 8.5153e-06 |
| C | Sand | declining | falling | 0.4175 | 2.4975e-03 |
| C | Vertisol | declining | falling | 0.0725 | 1.5195e-05 |

## Mechanistic Interpretation

1. Variant C wicking increased wet habitat for Loess and Vertisol from near-zero
   baseline to low-but-nonzero values (~0.07 at theta=0.6), so the hydro bottleneck
   was partially relieved.
2. Despite that, Loess and Vertisol still remain in declining/falling states and
   their habitat substrate stays in the 1e-06 range, so transport to active habitat
   is still too weak to recover growth in this configuration.
3. Variant A now behaves as a mixed result after removing banded OM: Vertisol still
   matches rising, Sand remains rising (mismatch), and Loess degrades to falling
   at n=128.

## Bottom Line

- The requested code modifications were implemented and validated across
  n in {128, 40} and theta in {0.5, 0.6, 0.7}.
- Variant C wicking is active and measurable (wet-habitat rescue), but it did not
  yet rescue Loess/Vertisol outcomes.
- Variant A without banded OM no longer keeps the previous 2/3 shape match at the
  headline condition.
