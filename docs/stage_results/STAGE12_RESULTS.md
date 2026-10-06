# Stage 12 RESULTS -- statistical significance via seed ensembles

`prompts/soil_respiration_prompt_stage12_significance.md`. An ANALYSIS
stage, not a modeling stage: Stage 11's two soils (Sand, Vertisol) and
every one of their frozen parameters (`om_total`, coating fraction,
`D_scale`, biology, `T`) are read from `configs/mapping.yaml` exactly as
Stage 11 left them and **never re-tuned here**. numpy-only (no scipy);
n=128; the Stage-11-optimized code made this affordable.

## 0. Result up front

**M=40 random pore-field seeds per soil (80 simulations total, 13.0
minutes).** Every shape metric tested shows **complete rank separation**
between the two soils' ensembles (Mann-Whitney rank-biserial `r = ±1.000`
-- every single Sand seed's value is on one side of every single Vertisol
seed's value), with `p` in the `1e-14`-`1e-15` range and Cohen's `d` from
5.0 to 36.6 (an effect size of `d=0.8` is conventionally called "large";
these are one to two orders of magnitude beyond that). The between-soil
mean curve distance (2.06) is **13x Sand's own within-soil spread (0.16)
and 127x Vertisol's (0.016)**; a bootstrap 95% CI on the between-minus-
within gap is `[1.963, 1.989]`, nowhere near zero. An optional causal
check (sweeping Sand's own coating fraction) shows the crash depth moves
monotonically with the swept feature, not just "Sand always crashes."

## 1. Ensemble setup

For each soil, M=40 seeds (`SEED_BASE`: Sand 10000-10039, Vertisol
20000-20039, disjoint ranges so the two ensembles are independent draws)
generate a fresh pore-field realization -- same `psd_soil`, `lambda`,
`aggregate`, and (crucially) the same frozen `om_total`/coating
fraction/`D_scale`/`T` from Stage 11 -- with ONLY `pore.seed` varying.
Runtime: Sand 657.2s (16.4s/seed, `T=3000`), Vertisol 120.9s (3.0s/seed,
`T=1500`) -- 778.0s (13.0 min) total, using the Stage-11-optimized
`simulate_dual_porosity`. Each soil's `R(t)` was resampled onto a common
500-point fraction-of-own-window axis (`normalized_time_resample`, the
same approach Stage 10/11 used for cross-window distinctness) before any
averaging or cross-soil comparison.

**Structural descriptors vary sensibly with the seed** (`stage12_seed_
metrics.csv`), confirming the ensemble is genuinely sampling different
random layouts of the same structural class, not accidentally re-running
the same realization:

| soil | n_regions (macropore count) | coating cells | dry-habitat fraction | matrix fraction |
|---|---|---|---|---|
| Sand | 2.6 +/- 1.7 (range 1-9) | 1840 +/- 25 | 0.999 +/- 0.001 | 0.1834 (seed-invariant) |
| Vertisol | 71.8 +/- 6.3 (range 55-86) | 4347 +/- 145 | 0.927 +/- 0.012 | 0.6775 (seed-invariant) |

`matrix_fraction` is exactly seed-invariant for a given soil (SD=0.0) --
expected, not a bug: the empirical-quantile PSD mapping guarantees every
seed reproduces the SAME marginal diameter distribution (hence the same
volume fraction below the 30um matrix/macropore threshold), only the
SPATIAL ARRANGEMENT of which cells get which diameter changes with the
seed, which is exactly what drives `n_regions`/`coating_cells` variation.
This is the structural contrast the whole stage tests: Sand consistently
realizes as a FEW-region, small-matrix (18%) soil; Vertisol consistently
realizes as a MANY-region, large-matrix (68%) soil, regardless of which
random layout the seed draws.

## 2. Ensemble mean curves + 95% bands

`results/stage12/stage12_mean_bands.png` -- both soils' mean R(t) (over
all 40 seeds, on the normalized time axis) with a 95% band (mean +/-
1.96*SE) and the 2.5/97.5 percentile envelope. The two mean curves and
their bands do not visually overlap at any point on the normalized axis:
Sand's mean rises sharply to an early peak (normalized time ~0.17) then
crashes; Vertisol's mean climbs gradually to a peak essentially at the
very end (~0.99) and stays there. The bands are narrow relative to the
between-soil gap -- seed-to-seed noise is small next to the structural
difference.

## 3. Per-seed shape metrics (mean +/- SD, `stage12_metric_summary.csv`)

| metric | Sand | Vertisol |
|---|---|---|
| `R_peak` | 35.17 +/- 8.93 | 3.64 +/- 0.37 |
| `peak_frac` | 0.167 +/- 0.022 | 0.991 +/- 0.023 |
| `R_end_over_peak` | 0.021 +/- 0.122 | 1.000 +/- 0.001 |
| `cum_co2_final` | 1195.96 +/- 39.06 | 149.55 +/- 16.53 |
| `late_over_early` | 0.044 +/- 0.179 | 9.216 +/- 1.064 |

Every metric's Sand and Vertisol distributions are far apart relative to
their own spread (box plots: `results/stage12/stage12_metric_
distributions.png`). Note `R_end_over_peak`'s SD for Sand (0.122) looks
large relative to its mean (0.021) only because a handful of the 40 seeds
draw a pore field where Sand's habitat happens to fragment into more than
one region (up to 9, vs. the typical 1-3) and the crash is milder in those
realizations -- still nowhere near Vertisol's near-1.0 values (see SS4's
rank-based test, which is robust to this and still gives complete
separation).

## 4. Significance: between-soil vs. within-soil (§4)

**Metric-level (Mann-Whitney U, `stage12_significance.csv`):**

| metric | U (Sand) | z | p (two-sided) | Cohen's d | rank-biserial r |
|---|---|---|---|---|---|
| `peak_frac` | 0 | -7.98 | **1.33e-15** | -36.63 | **-1.000** |
| `R_end_over_peak` | 0 | -7.98 | **1.33e-15** | -11.37 | **-1.000** |
| `late_over_early` | 0 | -7.69 | **1.42e-14** | -12.02 | **-1.000** |
| `R_peak` | 1600 | 7.69 | **1.42e-14** | 4.99 | **1.000** |
| `cum_co2_final` | 1600 | 7.69 | **1.42e-14** | 34.89 | **1.000** |

`rank_biserial_r = +-1.000` on every metric means COMPLETE separation:
`U=0` (or its maximum, `n1*n2=1600`) is the most extreme value the
Mann-Whitney statistic can take -- every one of the 40 Sand values lies
strictly on one side of every one of the 40 Vertisol values, for all five
metrics, with no overlap at all across the full 80-seed ensemble.

**Curve-level (`stage12_within_vs_between.csv`,
`stage12_within_vs_between.png`):**

| | mean pairwise distance | n pairs |
|---|---|---|
| within Sand | 0.158 | 780 (C(40,2)) |
| within Vertisol | 0.016 | 780 |
| between Sand-Vertisol | **2.063** | 1600 (40x40) |

Between-soil distance is **13.1x** Sand's own within-soil spread and
**127x** Vertisol's. Bootstrap CI (5000 resamples, resampling within- and
between-soil distance pools independently with replacement) on
`mean(between) - mean(within, pooled)`: **mean=1.976, 95% CI=[1.963,
1.989]** -- the gap excludes zero by a wide margin; the CI is narrow
because both distributions (1600 between-pairs, 1560 pooled within-pairs)
are large samples.

## 5. Interpretation: the structural features drive the shape

Because (a) each soil's structural CLASS (its frozen `om_total`, coating
fraction, `D_scale`, PSD, and every biology constant) is held perfectly
fixed while only the random pore-field LAYOUT changes between the 40
members of each ensemble, and (b) every seed of a given soil reliably
reproduces that soil's characteristic shape (Sand: early sharp peak then
crash; Vertisol: gradual climb to a late, sustained peak) with
between-soil separation 13-127x larger than the soil's own seed-to-seed
noise and complete rank separation (`p<2e-14`, `|r|=1.0`) on every metric
tested -- **the respiration shape is determined by the defined structural
features, not by which particular random realization happened to be
drawn.** Sand's shape comes from having few (typically 1-3, up to 9)
concentrated, coating-fed macropore regions in a small (18%) matrix;
Vertisol's comes from having dozens (55-86) of distributed regions fed by
a large (68%) matrix, with a shared, frozen transport/dose regime applied
identically to both structural classes. The random SEED changes exactly
which cells end up in which region and how many regions form (SS1's
table) -- but never enough to flip either soil out of its own
characteristic shape family.

## 6. Optional causal check: the feature itself is the driver

`results/stage12/stage12_coating_sweep.png`, `.csv`: Sand's OWN coating
fraction (the fast-vs-slow OM placement lever, `stage11_coating_fraction_
sand`, a highest-priority per-soil override used only for this bounded
sweep) was swept over `{0.50, 0.70, 0.85, 0.97}` (0.97 = the frozen Stage
11 value), 5 seeds each (20 runs, ~6 minutes), all other parameters
(including `om_total`, `D_scale`, `T`) held at Stage 11's frozen values:

| coating fraction | `peak_frac` | `R_end/R_peak` (crash depth) |
|---|---|---|
| 0.50 | 0.197 +/- 0.026 | 0.055 +/- 0.007 |
| 0.70 | 0.171 +/- 0.017 | 0.022 +/- 0.002 |
| 0.85 | 0.170 +/- 0.032 | 0.006 +/- 0.002 |
| 0.97 (frozen) | 0.151 +/- 0.012 | **0.001 +/- 0.000** |

`R_end/R_peak` drops **monotonically** as the coating fraction rises (more
of the OM sits on the fast, near-habitat pool and less in the slow trapped
tail that would otherwise prop up `R(t)` after the burst) -- a clean,
directly causal relationship between the swept FEATURE and the crash
depth, not merely "Sand always crashes regardless of parameters." Peak
timing also shifts modestly earlier with more coating (0.197 -> 0.151),
consistent with faster, more concentrated delivery. This is supporting
evidence, as scoped (SS5's "optional... bounded... supporting evidence,
not the main deliverable"), not a re-tuning of the frozen Stage 11
configuration (`configs/mapping.yaml`'s checked-in `stage11_coating_
fraction_vertisol`/Sand's Stage 10/11 values are untouched by this
sweep -- it uses a throwaway local `mapping_cfg` copy).

## 7. Deliverables

`stage12_run.py`: the M=40-seed ensemble loop (2 soils, n=128, Stage-11-
frozen parameters, only `pore.seed` varying); `normalized_time_resample`
(cross-window curve alignment); numpy-only `rankdata`/`mann_whitney_u`/
`cohens_d`/`rank_biserial`/`bootstrap_gap_ci` (no scipy); the mean+band,
metric-distribution, and within-vs-between figures; the optional bounded
coating-fraction sweep. `results/stage12/`: `stage12_mean_bands.png`,
`stage12_metric_distributions.png`, `stage12_within_vs_between.png`,
`stage12_coating_sweep.png`; `stage12_seed_metrics.csv` (every seed's
metrics + structural descriptors), `stage12_metric_summary.csv`,
`stage12_significance.csv` (U/z/p/Cohen's d/rank-biserial per metric),
`stage12_within_vs_between.csv` (bootstrap gap CI), `stage12_coating_
sweep.csv`.

**[Stage 12] Checkpoint reached: M=40-seed ensembles for both soils (Stage
11 parameters frozen, confirmed never re-tuned), mean+band curves,
per-seed metric distributions, Mann-Whitney + effect size + within/between
curve-distance bootstrap CI all complete and showing complete separation
(`p<2e-14`, `|rank-biserial r|=1.0`, gap CI excludes 0), plus the optional
bounded causal sweep confirming the swept feature (coating fraction)
moves the shape metric monotonically; deliverables written (this file,
MODEL_SPEC.md, LOGBOOK.md, report_prompt.md). Stopping here per
soil_respiration_prompt_stage12_significance.md SS7 order of work.**
