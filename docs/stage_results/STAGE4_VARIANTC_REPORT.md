# Stage 4 Variant C Report

## Scope

This report isolates what happened to Variant C in Stage 4 testing.
Variant C is the dual-porosity matrix/macropore path implemented in
`dual_porosity.py` and executed by `stage4_run.py`.

Test frame (same as Stage 4 head-to-head):

- Soils: Loess, Sand, Vertisol
- Constant moisture mode (`global_theta`) with theta in {0.5, 0.6, 0.7}
- Grid checks: n=128 (headline) and n=40 (robustness at theta=0.6)
- Comparator metric set: shape match, outcome class, wet-habitat diagnostics,
  substrate reaching habitat, pairwise distinctness

Primary sources:

- `results/stage4/stage4_variantC_sweep.csv`
- `results/stage4/stage4_distinctness.csv`
- `docs/stage_results/STAGE4_RESULTS.md`

## What Variant C was supposed to fix

Stage 3.5 found a transport failure: most OM sat in very fine pores and did not
reach the 30-150 um habitat in time.

Variant C addresses this by splitting pore space into:

- macropore cells (d >= 30 um): biomass and fast diffusion pathway
- matrix cells (d < 30 um): no biomass, fast-flow blocked
- matrix->macropore dissolved-substrate leakage at fixed rate `k_leak = 0.10`

Critical implementation detail: Variant C keeps the full PSD and does not alter
the global-theta wetting quantile logic.

## Observed behavior in tests

### 1) Headline run (theta=0.6, n=128)

From `stage4_variantC_sweep.csv`:

- Loess: declining, model shape falling (experiment is erratic),
  wet_habitat_fraction = 0.0, mean_S_habitat_final = 1.53e-6
- Sand: declining, model shape falling (experiment is falling),
  wet_habitat_fraction = 0.4166, mean_S_habitat_final = 1.49e-3
- Vertisol: declining, model shape falling (experiment is rising),
  wet_habitat_fraction = 0.0, mean_S_habitat_final = 3.14e-6

Net: only Sand matches the measured shape; Loess and Vertisol both collapse to
very similar starvation-like trajectories.

### 2) Theta robustness (0.5, 0.6, 0.7 at n=128)

Loess and Vertisol stay in the same failure mode across all tested theta values:

- outcome always declining
- period shape always falling
- wet_habitat_fraction remains 0.0 at theta=0.5 and 0.6, and only about 0.022 at
  theta=0.7 (still extremely dry in habitat terms)
- mean_S_habitat_final stays near 1e-6 to 3e-6

Sand remains the only soil with significant wet habitat and much larger habitat
substrate concentration (about 1.46e-3 to 1.50e-3).

### 3) Grid robustness (n=40 at theta=0.6)

Failure pattern is unchanged:

- Loess: declining, falling
- Vertisol: declining, falling
- Sand: declining, falling

The absolute magnitudes shift, but the qualitative conclusion does not.

### 4) Distinctness and degeneracy

From `stage4_distinctness.csv` for Variant C:

- theta=0.6, n=128: min_distance = 7.3895e-6 (threshold 0.3, fail)
- theta=0.5, n=128: min_distance = 7.3895e-6 (fail)
- theta=0.7, n=128: min_distance = 7.3284e-6 (fail)
- theta=0.6, n=40: min_distance = 4.2428e-5 (fail)

The smallest pair is consistently Loess vs Vertisol, indicating near-identical
curves in Variant C for those two soils.

## Why this happened (mechanistic explanation)

Variant C improved one pathway (exchange from matrix to adjacent macropores),
but it did not change the dominant bottleneck found in Stage 3.5:

- The habitat in Loess and Vertisol remains effectively dry under shared global
  theta wetting order.
- If habitat cells are mostly dry, the fast macropore pathway is not active
  enough to rescue growth.
- Leakage exists, but delivered substrate into habitat stays tiny for Loess and
  Vertisol (order 1e-6), far below Variant A levels (order 1e-2).

So Variant C inherits the same practical outcome: two soils starve similarly,
and emergent separation collapses.

## Variant C vs Variant A (concise contrast)

At theta=0.6, n=128:

- Variant C Loess/Vertisol wet_habitat_fraction: 0.0 / 0.0
- Variant A Loess/Vertisol wet_habitat_fraction: 0.136 / 0.523
- Variant C Loess/Vertisol mean_S_habitat_final: 1.53e-6 / 3.14e-6
- Variant A Loess/Vertisol mean_S_habitat_final: 1.09e-2 / 1.58e-2

Interpretation: Variant A accidentally fixes upstream wetting access (through
PSD truncation side-effects on quantile wetting), while Variant C leaves that
upstream condition mostly untouched.

## Verdict for Variant C in Stage 4

Variant C does not solve the Stage 3.5 failure for Loess and Vertisol in the
tested configuration.

- It reproduces Sand's falling trend, but Sand was already the least problematic
  soil in earlier stages.
- It fails shape recovery for Loess and Vertisol.
- It fails distinctness by a very large margin.
- Its mechanism is active but too weak relative to the unresolved habitat-wetting
  bottleneck.

Under the Stage 4 test protocol, Variant C is not sufficient on its own.
