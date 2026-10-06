# Stage 14 RESULTS -- investigating the `habitats_recruited(t)` saturation artifact

`prompts/soil_respiration_prompt_stage14_habitat_diagnostic.md`. A
DIAGNOSTIC stage, not a retuning stage: Stage 10/11/13's frozen calibrated
results are untouched, verified by an exact identity check, not just
convention.

## 0. Result up front

**The near-instant-equilibration hypothesis is CONFIRMED directly, not just
by elimination.** Of Vertisol's 70 macropore regions, 55 ever cross the
`growth > 1e-9` activation threshold at all, and every one of those 55 does
so within the first 2% of the observation window (median first-active
timestep = 2, out of `T=1500`) -- 34 of them are active by timestep 2
alone. This is exactly `count_active_habitats`'s per-step equilibration
crossing the threshold almost everywhere almost instantly, as Stage 11
suspected. A graded, growth-WEIGHTED metric (candidate: total raw growth
summed over all macropore cells each step, "active mass") tells a
completely different, and correct, story: it correlates with `R(t)`'s own
shape at `r=0.996`, against the boolean recruited-count's `r=0.072`
(effectively uncorrelated) -- confirming Stage 11's claim that `R(t)`'s
genuine gradual rise comes from already-active regions' growth RATE slowly
building up, not new regions joining, and that the boolean count is not
just incomplete but actively misleading as a rise-tracking diagnostic. Real
data (Stage 13) was not used for this run -- it collapses Vertisol to zero
growth everywhere, leaving no growth distribution to diagnose -- so this
uses Stage 11's frozen literature-PSD Vertisol configuration, exactly as
the prompt's own fallback allows.

## 1. Instrumentation (SS1) and identity check (SS3)

`dual_porosity.simulate_dual_porosity` gained one new, OFF-by-default
keyword, `debug_growth_distribution`: when `True`, it records the raw
per-cell `growth` value at every timestep, restricted to macropore cells,
as `growth_distribution_t` (shape `(T, n_macro_cells)`), alongside
`macro_mask_debug`/`region_labels_debug` for indexing it back onto
regions. This is a pure read-only side effect (one array assignment per
step, no branch that alters `B`, `S`, `OM_coating`, `OM_trapped`, or any
other state variable) -- confirmed, not just asserted:

| check | result |
|---|---|
| `R_t` exactly equal (debug vs. non-debug run) | **True** |
| `B_total_t` exactly equal | **True** |
| `habitats_recruited_t` exactly equal | **True** |
| `fuel_remaining_t` exactly equal | **True** |
| `R_peak` matches Stage 11 reference (<1%) | **True** |
| `B_peak/B0` matches Stage 11 reference (<1%) | **True** |

All arrays are BYTE-identical (`np.array_equal`, not merely close), the
strongest check available -- the flag introduces no floating-point
reordering or stochastic element, so exact equality is the correct bar
here, stronger than Stage 11's own noise-level tolerance.

## 2. Mechanism confirmation (SS1), with numbers

`stage14_region_recruitment.csv`, `stage14_mechanism_summary.csv`
(Vertisol, Stage 11's frozen config, n=128, T=1500):

| | |
|---|---|
| Total macropore regions | 70 |
| Regions that EVER cross `growth>1e-9` | 55 (matches `habitats_recruited_final=55` exactly) |
| Active by `t=0` | 0 |
| Active by `t<=2` | 34 (62% of the 55 that ever activate) |
| Active by `t<=30` (2% of the 1500-step window) | **55 (100%)** |
| Median first-active timestep | **2.0** (0.13% of the window) |

A direct cross-check ties this new instrumentation to the already-trusted
boolean diagnostic rather than treating it as a parallel, unverified
computation: reconstructing `habitats_recruited(t)` from each region's own
`first_active_t` (`count of regions with first_active_t <= t`, for every
`t`) matches the existing `habitats_recruited_t` array EXACTLY
(`np.array_equal`, confirmed in `stage14_mechanism_summary.csv`'s
`reconstructed_matches_habitats_recruited_t_exactly=True`).

**This directly confirms the hypothesis, not by elimination**: the growth
DISTRIBUTION itself (not just the boolean count derived from it) shows the
overwhelming majority of eventually-active regions crossing the threshold
within the first handful of timesteps -- `stage14_growth_distribution_
heatmap_zoomed.png` shows this visually as a near-vertical ignition front
at the very left edge of the window, immediately followed by a long,
slowly-brightening plateau (magnitude still building for hundreds of
timesteps after activation) rather than a diagonal recruitment front that
would indicate regions activating throughout the run.

## 3. A better graded recruitment metric exists (SS2)

Three candidate metrics were tested against `R(t)`'s own known-genuine
gradual rise (Stage 11: driven by per-region growth RATE building up, not
new regions joining):

| Candidate metric | correlation with normalized R(t) |
|---|---|
| `habitats_recruited(t)` (existing boolean count) | **0.072** (effectively uncorrelated) |
| growth-weighted "active mass" (sum of raw growth over all macropore cells, every step) | **0.996** |
| recruited count at a 10x-higher threshold (`0.1 * max single-cell growth` over the whole run) | **0.960** |

**The growth-weighted active-mass metric tracks `R(t)`'s real rate-building
story almost exactly** (`r=0.996`), visually confirmed in
`stage14_graded_metric_vs_R.png`: while `habitats_recruited(t)/max` jumps
to 1.0 within the first few timesteps and sits flat for the rest of the
run, `active_mass(t)/max` rises smoothly alongside `R(t)/max` across the
whole window. The raising-the-threshold approach also works
substantially better than the boolean count (`r=0.960`) but is less clean
than the fully graded active-mass metric (still throws away magnitude
information above the new, higher threshold).

**Per-region time-to-half-max** (a third graded metric, region-level
rather than a single run-wide scalar): for each of the 55 active regions,
the first timestep its OWN mean growth reaches half its OWN eventual peak.
Unlike `first_active_t` (median 2.0, 0.13% of window), this is genuinely
spread across the window: **median = timestep 659, 44.0% of the window**
(`stage14_activation_vs_halfmax_histograms.png` shows the two
distributions side by side -- boolean activation is a sharp spike at the
left edge, half-max timing is spread broadly across the run). This
confirms per-region growth RATE really does build up gradually over most
of the observation window, exactly matching Stage 11's own interpretation
of what drives `R(t)`'s rise -- now demonstrated directly at the per-region
level, not inferred only from the aggregate curve.

## 4. Calibrated results untouched (SS3)

Confirmed above (Section 1): `R_t`, `B_total_t`, `habitats_recruited_t`,
and `fuel_remaining_t` are byte-identical between the debug-instrumented
and non-instrumented runs, and match Stage 11's own recorded reference
numbers (`R_peak=3.756`, `B_peak/B0=2.337`) to well under 1%. `om_total`,
coating fraction, `D_scale`, and `T` were never touched by this stage.

## 5. Honest verdict

**The near-instant-equilibration hypothesis is confirmed, directly, not
merely by elimination**: growth values across Vertisol's macropore regions
cross the activation threshold near-simultaneously (55/55 eventually-active
regions within the first 2% of the window, most within the first 2
timesteps out of 1500) -- exactly the unconditionally-stable implicit
solver's per-step equilibration reaching essentially the whole connected
macropore network almost immediately, independent of `D_scale`. **A
better recruitment metric is recommended as an ADDITION, not a
replacement**: the existing boolean `habitats_recruited(t)` is not merely
incomplete but actively misleading for tracking `R(t)`'s rise
(`r=0.072`) -- a growth-weighted "active mass" metric (`r=0.996`) should be
used alongside it in any future stage that needs to visualize or reason
about staggered recruitment, while the boolean count remains useful for
exactly what it already correctly measures (how many regions have SOME
nonzero activity at all, a genuinely near-instantaneous, structurally
real property of this solver/threshold combination -- not itself a bug).

## 6. Deliverables

`dual_porosity.py`: `simulate_dual_porosity` gained the OFF-by-default
`debug_growth_distribution` keyword (purely observational, verified via
exact identity check). `stage14_run.py`: baseline + debug run, identity
check, per-region recruitment analysis (boolean vs. graded half-max),
candidate graded metrics + correlations, all figures/CSVs in
`results/stage14/`. This file, `LOGBOOK.md`, `report_prompt.md`, and
`MODEL_SPEC.md` (a one-paragraph note on the new documented
`debug_growth_distribution` kwarg only, per SS4's "only touch MODEL_SPEC
if a new diagnostic function is added to the documented API") updated.

**[Stage 14] Checkpoint reached: growth-distribution instrumentation added
(purely observational, exact identity check passes); equilibration
hypothesis CONFIRMED directly (55/55 active regions cross threshold within
2% of window, median t=2 of 1500); graded growth-weighted active-mass
metric recommended as an addition (r=0.996 vs. R(t), vs. boolean count's
r=0.072); per-region time-to-half-max independently confirms R(t)'s
genuine gradual rise (median 44% of window); all Stage 10/11/13 calibrated
results verified byte-identical; deliverables written. Stopping here per
soil_respiration_prompt_stage14_habitat_diagnostic.md SS5 order of work.**

---

# CORRECTION (2026-09-16, Stage 16) -- Section 5's causal attribution was WRONG

Appended, not substituted: everything above is left exactly as written, in
line with how this project has handled every other revised finding. This
section corrects ONE claim and leaves the rest of the stage standing.

## C.1 What is wrong

Sections 0 and 5 above attribute Vertisol's near-instant habitat activation
to **"the unconditionally-stable implicit solver's per-step equilibration
reaching essentially the whole connected macropore network almost
immediately, independent of `D_scale`."** That attribution is **refuted by
direct test** and must not be inherited by later work.

The reasoning error was one of elimination presented as direct
confirmation: Stage 14 observed *that* activation is near-simultaneous
(which is true, and remains true) and then named the solver as the cause
without ever running the model with transport switched off. The observation
was never in question; only the mechanism assigned to it.

## C.2 The refuting test

Vertisol, Stage 11's frozen config, n=128, T=1500 -- the identical setup
Stage 14 itself used. Re-run and reproduced independently during Stage 16
(`prompts/soil_respiration_prompt_stage16_shared_om_rule.md` SS0):

| run | regions active @ t=2 | recruitment reaches final at | R_peak |
|---|---|---|---|
| baseline (`D_scale=0.7`) | 34 | 1.33% of window | 3.756 |
| **transport OFF (`D_scale=1e-9`)** | **17** | **0.27% of window** | 0.400 |
| transport OFF + `coating_fraction=0` | 0 | never activates | 0.400 |
| transport ON + `coating_fraction=0` | **0** | **9.27% of window** | 0.400 |

With diffusive transport switched off entirely -- `D_scale` reduced by nine
orders of magnitude, so the solver moves essentially nothing between cells
-- regions still activate within 2-4 timesteps, and recruitment reaches its
final count *sooner* (0.27% of the window) than with transport on (1.33%).
If the solver's per-step equilibration were the cause, disabling transport
would have to slow or abolish the simultaneity. It does neither.

## C.3 The actual cause

Structural, and it lives in `split_om_coating_trapped`, not in
`implicit_diffusion_step`. The coating pool is defined as every
porous-matrix cell within `coating_thickness` geodesic steps of ANY
macropore cell. That construction builds an identical fuel shell around all
~70 regions at once; combined with `B0 ~ K` seeding biomass in every region
at `t=0` and a single global `k_dis`, every region is self-fuelled from
`t=0` and needs no transport at all to activate. **Simultaneity was built
into the OM placement, not produced by the numerics.**

The last row of the table is the decisive one and points at the fix:
removing the coating (`coating_fraction=0`) with transport left ON restores
genuinely staggered recruitment (final count reached at 9.27% of the
window, ~7x later than baseline, with nothing active at `t=2`). This is the
result Stage 16 builds its shared-OM-rule change on.

## C.4 What still stands

Everything except the causal claim. Specifically retained, unmodified:

- the `debug_growth_distribution` instrumentation and its exact identity
  check (Section 1) -- purely observational, still verified;
- the measured per-region recruitment numbers (Section 2) -- 55 of 70
  regions ever active, 34 by `t=2`, median first-active `t=2.0`. These are
  observations, and they are correct; only their explanation changes;
- the metric comparison (Section 3): boolean `habitats_recruited(t)` is
  actively misleading as a rise-tracking diagnostic (`r=0.072` against
  `R(t)`) while the growth-weighted "active mass" metric tracks it
  (`r=0.996`), and per-region time-to-half-max spreads genuinely across the
  window (median 44%). Stage 16 adopts this recommendation and retires the
  boolean count as a headline number;
- Section 4's confirmation that Stage 10/11/13's calibrated results were
  untouched by Stage 14.

**Corrected verdict for Section 5:** Vertisol's near-instant activation is
a consequence of the per-soil coating construction placing fuel around
every region simultaneously -- a structural artifact of the OM placement
rule -- NOT of the implicit solver. The solver is exonerated.
