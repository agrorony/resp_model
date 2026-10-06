# Stage 4 (EXPLORATORY) RESULTS -- two fixes for "OM never reaches the habitat", head-to-head

`prompts/soil_respiration_prompt_stage4_exploratory.md`. `MAX_ITERATIONS = 12`
was not exercised: this is a defined comparison (build both variants, run,
compare, write up), not a search.

## 0. The problem being fixed

Stage 3.5 diagnosed why the model reproduced only Sand-falling: at a shared
theta=0.6, the frozen `OM ~ d^-2.5` rule dumps ~99.96% of carbon onto
sub-micron clay cells, and diffusion (reach ~1-3 cells over the run) never
delivers it to the 30-150 um habitat -- Loess and Vertisol starve identically.
Two independent fixes were built and compared, without tuning either toward
the target pattern (Vertisol rising, Sand falling, Loess low/erratic):

- **Variant A** ("truncate + relocate OM"): stays on the single-continuum
  grid, but (1) drops all pore volume below a 10 um floor and renormalizes,
  and (2) replaces the unbounded fine-pore OM rule with a band centered just
  below the habitat (10-30 um).
- **Variant C** (dual-porosity matrix/macropore): keeps the full literature
  PSD, splits every cell into a fast macropore (d >= 30 um, biomass lives
  here) or a matrix cell (d < 30 um, no biomass, fast diffusion blocked), and
  lets matrix cells leak dissolved substrate into adjacent macropore cells at
  a fixed rate `k_leak = 0.10`.

Both ran at the same constant theta = 0.6 (Stage 3.5's headline condition),
n=128, plus a theta-robustness check (0.5/0.6/0.7) and an n=40 grid check.
Full numbers: `results/stage4/stage4_variantA_sweep.csv`,
`stage4_variantC_sweep.csv`, `stage4_distinctness.csv`,
`stage4_variantA_truncation.csv`.

## 1. Variant A: truncation + relocation, per soil

| soil | porosity before->after | dropped fraction | habitable frac before->after |
|---|---|---|---|
| Loess    | 0.470 -> 0.283 | 0.399 | 0.242 -> 0.404 |
| Sand     | 0.380 -> 0.325 | 0.145 | 0.487 -> 0.566 |
| Vertisol | 0.550 -> 0.192 | 0.651 | 0.258 -> 0.737 |

Dropping pore volume below 10 um removes the largest share of Vertisol's
total porosity (65%, since its clay/intra-aggregate modes sit at 0.01-0.35
um) and the smallest share of Sand's (14%, which barely has any sub-10um
volume to begin with) -- exactly the soils that most needed it.

**A side effect nobody engineered for turned out to be the mechanism that
actually fixes Stage 3.5's dry-habitat problem.** `global_theta` sets its
water-fill cutoff as *this soil's own theta-quantile of its pore field*.
Removing the sub-10um tail from that field shifts the whole quantile
distribution toward the coarse end, so at the *same* nominal theta=0.6, a
much larger share of each soil's 30-150 um habitat ends up wet:

| soil | wet-habitat frac, Stage 3.5 (untruncated) | wet-habitat frac, Variant A |
|---|---|---|
| Loess    | 0.000 | 0.136 |
| Sand     | (already wet) | 0.555 |
| Vertisol | 0.000 | 0.523 |

All three soils come out **thriving** (`B_final/B0` 1.2-2.4, `outcome` =
thriving in every theta/grid cell tested) -- the pure-starvation dead curve
that killed two of three soils in Stage 3.5 is gone.

**Shape match at the headline condition (theta=0.6, n=128):**

| soil | model shape | experiment | match? |
|---|---|---|---|
| Loess    | erratic | erratic | **yes** |
| Vertisol | rising  | rising  | **yes** |
| Sand     | rising  | falling | no |

Two of three shapes now match -- a reversal from Stage 3.5 (which matched
only Sand). **Robustness:** Loess-erratic and Vertisol-rising hold at every
theta in {0.5, 0.6, 0.7} at n=128 (`stage4_variantA_sweep.csv`); Sand stays
wrongly "rising" at every theta too (a robust miss, not a knife-edge one). The
Vertisol match is **not** grid-robust: at n=40 its shape flips to "erratic"
(same n=40 fragility Stage 3.5 already documented). Pairwise distinctness at
the headline condition does **not** clear the 0.3 threshold
(`min_distance=0.137`, Sand-vs-Vertisol -- both now thrive and both trend
upward, so they look more alike than Stage 3.5's dead-vs-thriving contrast
did) -- reproducing the qualitative *shape* for 2/3 soils and reproducing the
*emergent-distinctness* criterion are not the same thing here.

## 2. Variant C: dual-porosity matrix/macropore, per soil

Matrix fraction (share of cells with d < 30 um) tracks each soil's fine-pore
content almost exactly as intended: Loess 0.678, Vertisol 0.677, **Sand
0.183** -- "Vertisol/Loess = lots of high-internal-porosity matrix, Sand =
little" is qualitatively captured by construction.

**Variant C does NOT fix the dry-habitat problem, because it never touches
the mechanism that causes it.** `global_theta`'s water-fill quantile is still
computed on the FULL (untruncated) PSD, so Loess's and Vertisol's macropore
habitat is bone dry at theta=0.6 exactly as in Stage 3.5
(`wet_habitat_fraction = 0.000` for both, at every theta in {0.5,0.6,0.7}).
The dual-porosity split only unblocks a DIFFERENT bottleneck -- fast
macropore-macropore diffusion, which was never the blocker for these two
soils in the first place, since their macropores were never wet enough to
conduct at all. The `k_leak` exchange does still fire (it is not gated by
`water_mask`), but it is bounded by (a) how much of a soil's OM sits in a
matrix cell that is *directly adjacent* to a macropore cell -- most matrix
volume for Loess/Vertisol sits in the interior of large matrix patches,
untouched by any macropore neighbor -- and (b) the shared rate itself
(`k_leak=0.10`, same order as `k_dis=0.05`). The trickle that does arrive
(`mean_S_habitat_final` ~1.5e-6, vs ~0.01-0.02 for Variant A -- four orders
of magnitude less) cannot outrun maintenance decay: both soils'
`outcome = declining`, `R_peak` pinned at the trivial `m0*B0_total = 0.4`
(growth never turns on at all), and their `R(t)` curves come out
**functionally identical** -- `Loess_vs_Vertisol` distance = 7.4e-6, i.e. the
exact same "two soils starve to the same dead curve" failure mode Stage 3.5
reported, now reproduced under a structurally different mechanism.

Sand -- whose macropores were already wet without any help (matrix fraction
only 0.183) -- is the one soil where the leak pathway has room to matter, and
it does reproduce the correct **falling** shape robustly across every
theta/grid tested. But this is the same soil Stage 3.5 already got right
without any of this machinery; Variant C adds nothing there and fixes nothing
for the other two.

**Distinctness:** min_distance = 7.4e-6 at the headline condition -- as
degenerate a failure as `emergent_distinctness` can register.

## 3. Comparison figure and mechanism diagnostic

`results/stage4/stage4_comparison_periods.png` -- Variant A (top row) vs
Variant C (bottom row) vs measured P1-P4 shape, three soils each.
`results/stage4/stage4_substrate_reaching_habitat.png` -- the killer metric
side by side: mean substrate concentration actually sitting in habitat cells,
Variant A order-of-magnitude higher for Loess/Vertisol, both variants
comparable for Sand (the one soil that was never the problem). Per-soil
water/pore-field maps: `stage4_water_map_variantA_{soil}.png` (pore diameter,
K, D, water mask, banded-OM log-density) and
`stage4_water_map_variantC_{soil}.png` (adds the macro/matrix split panel).

## 4. Honest verdict (SS4/SS7)

**Variant A suffices better; Variant C's dual-porosity exchange, as built and
at the shared `k_leak` chosen, does not.** Variant A reproduces 2 of 3
measured shapes (Loess-erratic, Vertisol-rising) robustly across theta,
though the emergent-distinctness threshold is still not cleared and the
Vertisol match is grid-sensitive (n=40 flips it). Variant C reproduces only
Sand -- the soil that never needed fixing -- and reduces Loess/Vertisol to
the exact same degenerate dead-curve collapse Stage 3.5 already documented,
because its fix (unblocking fast intra-macropore diffusion, plus a bounded
matrix-to-macropore leak) targets a bottleneck that was never binding for
those two soils; the actual bottleneck (the shared-theta wetting order
leaving their habitat dry) is untouched by construction.

This is a mechanistically clean result, not a coincidence of tuning: Variant
A's truncation changes *what counts as fine* in the per-soil theta-quantile
calculation, incidentally rewetting the habitat as a side effect of removing
its finest pores from the distribution. Variant C's split leaves that
quantile calculation, and therefore the wetting order, completely alone --
so of course it inherits the same bone-dry habitat. **The real fix for
"OM never reaches the habitat" (Stage 3.5's diagnosis) was masking a deeper,
prior fix that was actually needed: "habitat never gets wet" under a shared
theta.** Variant A only reaches the OM-mobility fix at all because its PSD
surgery happens to also cure the wetting problem; a version of Variant C that
additionally decoupled the macropore habitat's wetting from the matrix's
(e.g. air-entry diameter applied only within the macropore fraction) would be
a fairer test of the dual-porosity idea on its own terms, but that is a
further stage, not this one -- per SS7, neither variant was tuned toward the
target, and this asymmetry is reported, not resolved.

## 5. What both variants still get wrong

- Neither variant reproduces Sand-falling AND the other two simultaneously;
  Variant A's Sand comes out rising (both Sand and Vertisol now thrive and
  trend upward, so they under-separate) and Variant C's Sand-falling is
  correct but isolated.
- Variant A's fix is partly an accounting artifact of the shared-theta rule
  interacting with a truncated PSD, not a first-principles argument that 10um
  is the "right" floor -- a different floor would shift the wet-habitat
  fractions by soil-specific amounts (Vertisol's clay-heavy PSD is far more
  sensitive to the floor choice than Sand's).
- Variant C's `k_leak=0.10` and the hard 30um macro/matrix split are
  documented choices, not calibrated; a much larger `k_leak`, or coupling the
  matrix's own wetting to the shared theta (currently the full PSD, not the
  macropore-only fraction, feeds the quantile cut) might change the verdict
  and was deliberately not tried here (this is a comparison, not a search).
- The literature PSDs remain representative of soil TYPES, not measurements
  of the specific experimental samples (carried over from Stage 3.5).

## 6. Deliverables

`pore_field.truncate_psd_floor`, `pore_field.om_field_banded`,
`model.build_grids`'s `psd_truncate_floor_um`/`om_rule` pore-config keys
(Variant A); `dual_porosity.py` (`build_dual_grids`, `leak_exchange`,
`simulate_dual_porosity` -- Variant C); `metrics.mean_substrate_in_habitat`;
`model.simulate`/`dual_porosity.simulate_dual_porosity`'s new
`S_habitat_mean_t` output; `configs/mapping.yaml`'s
`stage4_psd_truncate_floor_um`, `om_band_center_um`, `om_band_width_log`,
`k_leak`; `stage4_run.py` (the comparison runner); all figures/CSVs in
`results/stage4/`.

**[Stage 4] Checkpoint reached: both variants built, run across
theta/grid robustness, compared via the mechanism diagnostics (wet-habitat
fraction, mean substrate reaching habitat, isolation ratio, emergent
distinctness); deliverables written (this file, MODEL_SPEC.md, LOGBOOK.md,
report_prompt.md, README.md). Stopping here per
soil_respiration_prompt_stage4_exploratory.md SS6 order-of-work --
MAX_ITERATIONS=12 was not needed (a defined comparison, not a search).**
