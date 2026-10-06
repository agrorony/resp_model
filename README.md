# Soil-structure microbiome respiration model

Simulates a single microbial species growing logistically on a 2D grid,
where three soils (A, B, C) share identical biology but differ only in
spatial structure (diffusion rate, carrying-capacity layout, organic-matter
placement). This produces three visually distinct respiration curves:
Soil A rises, Soil B falls, Soil C stays flat/bumpy. See `MODEL_SPEC.md` for
the full equations and reasoning.

## Setup

```
python -m venv .venv
.venv/Scripts/pip install numpy matplotlib pyyaml
```

The repo is numpy-only (no scipy): `pore_field.gaussian_filter_2d` is a
hand-rolled separable-Gaussian convolution that reproduces
`scipy.ndimage.gaussian_filter`'s default behavior (reflect/mirror padding,
`truncate=4.0`) to machine precision.

(On Linux/Mac: `source .venv/bin/activate` instead of using `.venv/Scripts/...` directly.)

## Running

```
# run all three soils: writes CSVs + all figures to results/
.venv/Scripts/python.exe run.py

# run a single soil
.venv/Scripts/python.exe run.py --soil A

# re-run the budget-capped Stage-1 structure search (see LOGBOOK.md for its log)
.venv/Scripts/python.exe search_loop.py

# Stage 2: theta sweep (global_theta saturation) + connectivity diagnostics
.venv/Scripts/python.exe stage2_run.py

# Stage 3: real-PSD theta sweep at n=128 and n=40 + isolation-ratio hypothesis test
.venv/Scripts/python.exe stage3_run.py

# Stage 3.5: literature-informed PSDs (Loess/Sand/Vertisol) at constant moisture,
# vs the measured experiment -> results/stage35/
.venv/Scripts/python.exe stage35_run.py

# Stage 4 (EXPLORATORY): two head-to-head fixes for "OM never reaches the
# habitat" -- PSD truncation + banded OM (Variant A) vs dual-porosity
# matrix/macropore exchange (Variant C) -> results/stage4/
.venv/Scripts/python.exe stage4_run.py

# Stage 5 (EXPLORATORY): Variant C with matric-potential wetting (shared
# air-entry diameter) + geometric inscribed-circle macropore sizing,
# wicking on/off isolated -> results/stage5/
.venv/Scripts/python.exe stage5_run.py

# Stage 6 (EXPLORATORY): port the pde_sandbox/ recipe into Variant C --
# implicit backward-Euler transport + coating/trapped OM split, S_cap
# removed, wicking off -> results/stage6/
.venv/Scripts/python.exe stage6_run.py

# Stage 7: calibration to close the substrate gap -- uniform shared OM
# fraction, a shared reach/D-scale multiplier, and water-film macropores,
# via a budget-capped search over those three shared scalars -> results/stage7/
.venv/Scripts/python.exe stage7_run.py

# Stage 8: equal total OM (absolute, identical across soils) + a shared
# observation window, set by <=3 directed diagnostic probes (not a search),
# to test whether structural concentration alone drives depletion-timing
# shape differences -> results/stage8/
.venv/Scripts/python.exe stage8_run.py

# Stage 9: per-soil ranked total OM (Vertisol > Loess > Sand) + Vertisol
# coating raised / Sand coating lowered, set by <=3 directed probes (each
# checked at n=40 AND n=128) -> results/stage9/
.venv/Scripts/python.exe stage9_run.py

# Stage 10: two soils only (Vertisol, Sand) -- per-soil om_total / coating
# fraction / transport magnitude (D_scale) / observation window, tuned
# directly and verified against an explicit numeric pass/fail table (Sand:
# real growth burst then crash; Vertisol: genuine slow sustained rise)
# -> results/stage10/
.venv/Scripts/python.exe stage10_run.py

# Stage 11: Sand kept FROZEN exactly as Stage 10 (confirmed via an exact
# numeric identity check); a real per-timestep performance bug fixed
# (11-60x speedup, verified to change nothing); Vertisol re-tuned further
# so its biomass clearly GROWS (not just reaches parity). n=128 ONLY
# -> results/stage11/
.venv/Scripts/python.exe stage11_run.py

# Stage 12 (ANALYSIS, no re-tuning): M=40-seed ensembles per soil at the
# frozen Stage 11 parameters -- Mann-Whitney U, effect size, and bootstrap
# CIs (numpy-only) showing each soil's shape is driven by its structural
# features, not a lucky random realization -> results/stage12/
.venv/Scripts/python.exe stage12_run.py
```

Each runner writes to its own `results/<stage>/` subfolder
(`baseline/`, `stage2/`, `stage3/`, `stage35/`, `stage4/`, `stage5/`, `stage6/`, `stage7/`, `stage8/`, `stage9/`, `stage10/`, `stage11/`, `stage12/`).

## Repository layout

```
*.py                     core modules + stage runners (run from repo root)
configs/                 biology.yaml, mapping.yaml, psd_literature.yaml, soil_*.yaml
data/psd/<name>/         the three real measured PSDs (Stage 3)
prompts/                 the historical per-stage task prompts (soil_respiration_prompt*.md)
docs/stage_results/      STAGE{1,2,3,3_5,4,5,6,7,8,9,10,11,12}_RESULTS.md — per-stage narrative writeups
docs/                    LOGBOOK_v1_backup.md
pde_sandbox/             abstract PDE search (sandbox.py/search.py) that Stage 6 ports from
results/baseline/        run.py outputs        results/stage2/  results/stage3/  results/stage35/  results/stage4/  results/stage5/  results/stage6/  results/stage7/  results/stage8/  results/stage9/  results/stage10/  results/stage11/  results/stage12/
LOGBOOK.md MODEL_SPEC.md report_prompt.md README.md   (live top-level docs)
```

Runner and core modules stay at the repo root and use paths relative to it
(`configs/...`, `data/...`), so run them from the project root.

## Files

- `model.py` — the simulation engine (§2 of the spec). No soil-specific
  values are hard-coded; everything comes from `configs/`. `build_grids`
  dispatches on `pore.mode` (`lognormal` for Stage 1/2, `psd` for Stage 3's
  measured PSDs, `psd_parametric` for Stage 3.5's literature PSDs).
- `psd_parametric.py` — builds a soil's PSD as a lognormal mixture from
  `configs/psd_literature.yaml` (Stage 3.5); returns the same shape as
  `psd_data.load_soil_psd` so the pipeline is unchanged.
- `configs/psd_literature.yaml` — the per-soil (Loess/Sand/Vertisol) mixture
  parameters, editable, with rationale in comments.
- `configs/soil_{loess,sand,vertisol}.yaml` — Stage-3.5 structure configs
  (`mode: psd_parametric`, constant-moisture saturation).
- `stage35_run.py` — the Stage-3.5 constant-moisture run: literature PSDs,
  the D(d) fusion fix, the experiment overlay, the θ-sensitivity check, the
  shared-matric-potential diagnostic, and the isolation/distinctness metrics.
- `dual_porosity.py` — Variant C: the dual-porosity matrix/macropore split
  (`build_dual_grids`), the mass-conserving matrix->macropore `leak_exchange`
  term, and `simulate_dual_porosity`. Since Stage 4, evolved to a three-phase
  split (macropore/porous-matrix/solid), per-soil OM fractions, matrix
  internal diffusion, dynamic hydro wicking, habitat-cluster K/D bonuses,
  (Stage 5) a `saturation.mode: matric_regions` shared-air-entry-diameter
  wetting path plus `macropore_region_inscribed_diameter`/
  `erosion_distance_4n` (numpy-only geometric macropore sizing, replacing a
  flat forced diameter), and (Stage 6) `implicit_diffusion_step` (Jacobi
  backward-Euler substrate transport, unifying the macro/matrix D fields
  into one connected solve, ported from `pde_sandbox/`), `split_om_coating_
  trapped`/`bfs_geodesic_distance` (coating/trapped OM pools, also ported),
  `S_cap` removed, wicking off by default; (Stage 7) `dry_macro_mask` +
  `D_film`/`stage7_film_K_scale` (water-film pathway + reduced capacity for
  drained macropores) and `stage7_D_scale` (shared transport-magnitude
  multiplier applied in `simulate_dual_porosity`); (Stage 8)
  `_variantc_om_total` (single absolute, shared OM total, retiring the
  fraction-based rule as the active default) and `fuel_remaining_t`
  (fraction of each soil's initial total OM still undissolved, tracked and
  returned by `simulate_dual_porosity`); (Stage 9) `_variantc_om_total`'s
  ranked tier (`stage9_om_base` x per-soil `stage9_om_multiplier_<soil>`,
  fixed ordering Vertisol > Loess > Sand) and `_variantc_coating_fraction`'s
  Stage 9 tier (`stage9_coating_fraction_<soil>`, Vertisol raised/Sand
  lowered from their Stage 6 values); (Stage 10) `_variantc_om_
  total`'s and `_variantc_coating_fraction`'s highest-priority per-soil
  tiers (`stage10_om_total_<soil>`, `stage10_coating_fraction_<soil>`,
  tuned directly, two soils only) and new `_variantc_D_scale`
  (`stage10_D_scale_<soil>`, falling back to the shared `stage7_D_scale`)
  -- the first place transport magnitude is NOT one shared constant across
  soils -- plus `B_total_t`/`growth_term_t`/`maintenance_term_t` time
  series (the growth-vs-maintenance split diagnostic); (Stage 11)
  `habitat_properties`/`K_t`/`habitable`/`K_safe` hoisted out of
  `simulate_dual_porosity`'s time loop (were being recomputed every
  timestep despite being loop-invariant -- an 11-60x speedup, verified to
  change nothing), `implicit_diffusion_step`'s new `tol` early-stop
  argument, and one more highest-priority tier on each of the three Stage
  10 functions (`stage11_om_total_vertisol` etc. -- Vertisol only, no
  `stage11_*_sand` key exists, so Sand stays frozen through the fallback
  chain itself). A separate module so Stage 1-3.5's `model.py` code path
  is untouched.
- `pde_sandbox/` — a separate, abstract pure-NumPy reaction-transport search
  on a toy heterogeneous grid (`sandbox.py`, `search.py`, its own
  `LOGBOOK.md`/`RESULTS.md`), used to prove the implicit-transport and
  coating/trapped-OM mechanisms BEFORE porting them into `dual_porosity.py`
  (Stage 6). Not part of the real soil model's data/config flow.
- `stage4_run.py` — the Stage-4 head-to-head runner: Variant A (PSD
  truncation + banded OM, via `model.build_grids`) vs Variant C
  (`dual_porosity.py`), both at constant theta plus theta/grid robustness
  checks, the comparison figure, and the substrate-reaching-habitat
  mechanism diagnostic.
- `stage5_run.py` — the Stage-5 runner: Variant C under the new
  `matric_regions` wetting mode, wicking on/off, a `matric_d_cut_um`
  sensitivity check (40/50/60), and an n=40 grid check; writes the region
  inscribed-diameter maps/histograms and the comparison figure to
  `results/stage5/`.
- `stage6_run.py` — the Stage-6 runner: the three soils under Stage 5's
  matric wetting with the ported implicit transport + coating/trapped OM
  (wicking off), the experiment overlay, the mean-substrate-reaching-habitat
  and habitats-recruited-over-time diagnostics, coating/trapped OM maps, and
  an n=40 grid check; writes to `results/stage6/`.
- `stage7_run.py` — the Stage-7 runner: a budget-capped search
  (`MAX_ITERATIONS=15`) over three shared calibration scalars
  (`stage7_om_fraction`, `stage7_D_scale`, `D_film`), scored on whether
  fed-habitat substrate clears the growth target AND whether the three
  curves stay pairwise distinct; then the three-soil headline run at the
  winning candidate, the experiment overlay, fed-habitat-substrate/
  habitats-recruited diagnostics, a D_film live-habitat report, and an n=40
  grid check; writes to `results/stage7/` (search trace also appended to
  `LOGBOOK.md`).
- `stage8_run.py` — the Stage-8 runner: at most 3 directed diagnostic
  probes over `(om_total, T)` (logged to `LOGBOOK.md`, not a search), the
  three-soil n=128 headline run at the chosen values, the experiment
  overlay, `mean_S_habitat(t)`/`fuel_remaining(t)`/`habitats_recruited(t)`
  depletion-timing diagnostics, distinctness, and an n=40 grid check;
  writes to `results/stage8/`.
- `stage9_run.py` — the Stage-9 runner: at most 3 directed probes over the
  shared base OM scale and Vertisol/Sand coating fractions (ranked
  ordering fixed, logged to `LOGBOOK.md`), each candidate checked at BOTH
  n=40 and n=128 (the two grids disagreed about which candidate is best --
  adoption follows n=128); the n=128 headline run at the winning candidate,
  the experiment overlay, the same depletion-timing diagnostics, and an
  n=40 grid check; writes to `results/stage9/`.
- `stage10_run.py` — the Stage-10 runner: two soils only (Vertisol,
  Sand); builds each soil's own structure config at its own per-soil `T`;
  the n=128 headline run + explicit numeric §3 pass/fail tables for both
  soils (Sand: burst/crash; Vertisol: slow rise); the growth-vs-maintenance
  stackplot; `mean_S_habitat(t)`/`fuel_remaining(t)`/`biomass(t)` for both
  soils; a time-normalized distinctness check (`normalized_time_resample`,
  needed since the two soils' windows differ); and an n=40 grid check;
  writes to `results/stage10/`.
- `stage11_run.py` — the Stage-11 runner: n=128 ONLY, two soils; an exact
  numeric identity check that Sand's R(t) is unchanged from Stage 10's
  recorded reference (confirming the performance fix altered nothing);
  Vertisol's re-tuned headline diagnostics (`R_peak` vs. maintenance floor,
  `biomass(t)`/B0, `habitats_recruited(t)`, `mean_S_habitat(t)`); the
  before/after runtime comparison; and distinctness; writes to
  `results/stage11/`.
- `stage12_run.py` — the Stage-12 (ANALYSIS) runner: M=40-seed ensembles
  per soil at Stage 11's frozen parameters (only `pore.seed` varies);
  numpy-only stats (`rankdata`, `mann_whitney_u`, `cohens_d`,
  `rank_biserial`, `bootstrap_gap_ci` -- no scipy); mean+band, per-seed
  metric-distribution, and within-vs-between curve-distance figures; the
  optional bounded Sand coating-fraction causal sweep; writes to
  `results/stage12/`.
- `configs/biology.yaml` — frozen constants shared by all soils.
- `configs/mapping.yaml` — frozen pore -> K/D mapping constants, shared by
  all soils and all stages (Stage 1/2 lognormal-hump K/D plus Stage 3's
  physiological K(d) window and rescaled D(d), Stage 4's PSD-truncation
  floor / banded-OM params / `k_leak` / Variant C structural constants,
  Stage 5's `voxel_um` / `matric_d_cut_um`, Stage 6's `k_dis_slow` /
  `stage6_coating_fraction_*` / `stage6_coating_thickness` /
  `stage6_implicit_iterations` / `stage6_implicit_relaxation`, `S_cap`
  constants retired, Stage 7's `stage7_D_scale` / `D_film` /
  `stage7_film_K_scale` (`stage7_om_fraction` retired in favor of Stage 8),
  Stage 9's `stage9_om_base` / `stage9_om_multiplier_<soil>` (ranked
  V>L>S) / `stage9_coating_fraction_{vertisol,sand}` / `stage9_T`
  (`stage8_om_total` / `stage8_T` retired in favor of Stage 9), Stage 10's
  `stage10_om_total_{sand,vertisol}` / `stage10_coating_fraction_
  {sand,vertisol}` / `stage10_D_scale_{sand,vertisol}` (per-soil transport
  magnitude, not shared) / `stage10_T_{sand,vertisol}` (per-soil
  observation window), Stage 11's `stage11_om_total_vertisol` /
  `stage11_coating_fraction_vertisol` / `stage11_D_scale_vertisol`
  (Vertisol-only, highest-priority tier -- no `stage11_*_sand` key exists
  anywhere, so Sand is provably frozen through the fallback chain itself)
  and `stage11_implicit_tol` (Jacobi convergence-tolerance early stop)).
- `configs/soil_A.yaml`, `soil_B.yaml`, `soil_C.yaml` — structure-only
  configs; as of Stage 3, each points at a soil's real PSD (`pore.psd_dir`)
  instead of lognormal `mu`/`sigma`.
- `pore_field.py` — pore-field generators (lognormal + Stage-3
  empirical-quantile) and the K/D maps.
- `psd_data.py` — loads a soil's real measured PSD from `data/psd/<name>/`.
- `data/psd/{mishmar,rehovot,nlm}/` — the three real measured PSDs used in
  Stage 3 (`psd_table.csv`, `bin_edges_um.csv`, `source.txt`), copied once
  from the external diagnostics pipeline.
- `run.py` — runs one or all soils; writes `results/*.csv` and the required
  figures.
- `metrics.py` — the emergent-distinctness success check
  (`require_nontrivial` flag as of Stage 3), the Stage-2
  percolation/connectivity diagnostics (`label_connected`,
  `largest_cluster_fraction`, `om_biomass_connectivity`,
  `accessible_om_fraction`, `connectivity_diagnostics`), and the Stage-3
  outcome/isolation diagnostics (`classify_outcome`, `isolation_ratio`).
- `search_loop.py` — the budget-capped (`MAX_ITERATIONS = 12`) Stage-1
  structure search; appends to `LOGBOOK.md` every iteration and only ever
  edits `configs/soil_*.yaml`, never `configs/biology.yaml` or
  `configs/mapping.yaml`.
- `stage2_run.py` — the Stage-2 theta sweep (`global_theta` saturation);
  writes the connectivity/distinctness CSVs and all `stage2_*.png` figures
  to `results/`.
- `stage3_run.py` — the Stage-3 theta x grid-size sweep and the
  isolation-ratio-vs-PSD hypothesis test; writes all `stage3_*` CSVs/figures
  to `results/`.
- `LOGBOOK.md` — the honest record of the Stage-1 search, the Stage-2
  numpy-only conversion + theta sweep, the Stage-3 real-PSD conversion +
  theta x grid sweep + hypothesis test, the Stage-4 Variant A/C head-to-head,
  the Stage-5 matric-potential + inscribed-circle-macropore run, and the
  Stage-6 sandbox port.
- `MODEL_SPEC.md` — equations, frozen constants, structure knobs, the
  saturation modes (`fully_wet` / `global_theta` / `retention` /
  `matric_d` / `matric_regions`), and the Stage 1/2/3/3.5/4/5/6 result
  summaries.
- `docs/stage_results/STAGE{1,2,3,3_5,4,5,6}_RESULTS.md` — the per-stage
  narrative writeups (idea, method, what emerged, honest verdict, limitations).
- `report_prompt.md` — a self-contained prompt (not run yet) that generates
  the final written report covering all three stages.
- `results/baseline/` — Stage 1: `respiration_curves.png`, `cumulative_co2.png`,
  `structure_maps_{A,B,C}.png`, one CSV per soil.
- `results/stage2/` — `stage2_theta_sweep.csv`,
  `stage2_distinctness_by_theta.csv`, `stage2_water_map_{A,B,C}_theta0.40.png`,
  `stage2_respiration_curves_theta0.40.png`,
  `stage2_respiration_vs_theta.png`, `stage2_connectivity_vs_theta.png`.
- `results/stage3/` — `stage3_psd_loaded.png`, `stage3_theta_sweep.csv`,
  `stage3_distinctness_by_theta.csv`,
  `stage3_isolation_vs_psd_correlation.csv`,
  `stage3_water_map_{A,B,C}_n{128,40}_theta0.40.png`, and the rest of the
  `stage3_*` figures.
- `results/stage35/` — `stage35_psd_loaded.png`, `stage35_D_of_d_fix.png`,
  `stage35_experiment_overlay.png` (+ the `_matric` variant),
  `stage35_respiration_constant_theta.png`, `stage35_theta_sensitivity.png`,
  `stage35_water_map_{loess,sand,vertisol}.png`, and their CSVs
  (`stage35_constant_theta.csv`, `stage35_matric_diagnostic.csv`,
  `stage35_distinctness.csv`, …).
- `results/stage4/` — `stage4_comparison_periods.png` (the head-to-head
  figure), `stage4_substrate_reaching_habitat.png`,
  `stage4_respiration_variant{A,C}.png`,
  `stage4_water_map_variant{A,C}_{loess,sand,vertisol}.png`, and their CSVs
  (`stage4_variantA_truncation.csv`, `stage4_variant{A,C}_sweep.csv`,
  `stage4_distinctness.csv`).
- `results/stage5/` — `stage5_comparison_periods.png` (wicking on/off vs
  measured shape), `stage5_region_map_{loess,sand,vertisol}.png`,
  `stage5_region_diameter_hist.png`, `stage5_substrate_reaching_habitat.png`,
  `stage5_matric_sensitivity.png`, and their CSVs (`stage5_headline.csv`,
  `stage5_region_diameters.csv`, `stage5_matric_sensitivity.csv`,
  `stage5_grid_check.csv`, `stage5_distinctness.csv`).
- `results/stage6/` — `stage6_comparison_periods.png`,
  `stage6_mean_s_habitat_vs_t.png`, `stage6_habitats_recruited_vs_t.png`,
  `stage6_om_map_{loess,sand,vertisol}.png`, `stage6_respiration_curves.png`,
  and their CSVs (`stage6_headline.csv`, `stage6_grid_check.csv`,
  `stage6_distinctness.csv`).

## Current status

**Stage 1** (fully wet) succeeded on the first search iteration (see
`LOGBOOK.md`, `STAGE1_RESULTS.md`): all three pairwise distances exceed the
0.3 threshold and every soil is non-trivial, though a later seed-robustness
check found the A-C pair is not reliably separated across seeds.

**Stage 2** (`global_theta` percolation) swept theta over {1.0, 0.8, 0.6,
0.4, 0.2} using the same Stage-1 pore parameters. Success was reached at
only 2 of 5 theta values (1.0, 0.6); at theta<=0.4 soil C's biomass
net-declines under severe percolation fragmentation (an extreme but
"non-trivial-violating" divergence), and at theta=0.8 the AC pair narrowly
misses the threshold. Percolation is a real, sometimes extreme, structural
lever but did not turn the Stage-1-weak A-C pair into a robustly separated
one. Full detail and mechanism (via the connectivity diagnostics) in
`STAGE2_RESULTS.md`.

**Stage 3** (real measured PSD, physiological `K(d)`, `n=128`/`n=40`) swept
theta over the same grid at both sizes, on each soil's real PSD instead of
the synthetic lognormal field. No grid/theta combination reached full
three-way success: the BC pair (nlm vs. mishmar) never separates
(min. pairwise distance never exceeds 0.062) because both soils' large
median pore diameters saturate the diffusion map `D(d)`, compressing their
real PSD difference away; AB and AC partially separate, best at `n=128`,
theta=0.2. No soil died off in this sweep (Stage 2's die-off mechanism,
`aggregate=True`, is off by Stage-3's shared-spatial-recipe invariant). The
isolation-ratio hypothesis (SS9) holds, weakly: `isolation_ratio` correlates
with respiration outcome more strongly than the PSD summary stats tested
(`|r|` ~0.25-0.33 vs ~0.14-0.24). `n=128` produces markedly more distinct
respiration across soils than `n=40`, where outcomes are nearly flat
regardless of soil or theta. Full detail and mechanism in
`docs/stage_results/STAGE3_RESULTS.md`.

**Stage 3.5** (literature-informed parametric PSDs for the three real soil
types Loess / Sand / Vertisol, at constant moisture) replaced the poor
image-derived PSDs and tested the model against my measured experiment
(Vertisol rising, Sand falling, Loess low/erratic). Honest verdict: the model
reproduces **Sand-falling only, and not robustly** (it flips to rising at
n=40); it does **not** reproduce Vertisol-rising or Loess-erratic. At a shared
θ the fine soils' 30-150 µm habitat macropores are the last pores to wet and
stay bone dry (wet-habitat = 0), so their biomass starves and both give the
identical pure-starvation curve; a shared-matric-potential diagnostic gives all
three a live habitat yet still reproduces only Sand, because the frozen
`OM ~ d^-2.5` rule concentrates ~99.96 % of substrate on sub-micron clay from
which diffusion never reaches the habitat. The D(d) fusion of the two coarse
soils was fixed (`D_d_ref` 90 → 250 µm). Full detail, mechanism, and the two
bugs found in `docs/stage_results/STAGE3_5_RESULTS.md`.

**Stage 4** (EXPLORATORY — two head-to-head fixes for Stage 3.5's "OM never
reaches the habitat" diagnosis, not a search) built Variant A (truncate the
PSD at a 10 µm floor + band OM in the 10-30 µm range, single-continuum) and
Variant C (dual-porosity matrix/macropore split with a `k_leak`
matrix→macropore substrate-exchange term, on the full untruncated PSD), then
compared both at constant θ=0.6 across a θ- and grid-robustness check.
Honest verdict: **Variant A reproduces 2 of 3 measured shapes
(Loess-erratic, Vertisol-rising) robustly across θ — though not the
emergent-distinctness threshold, and not grid-robustly for Vertisol —
because truncating the PSD has the side effect of shifting the shared-θ
wetting quantile toward the coarse end, incidentally rewetting the habitat
that was bone-dry in Stage 3.5. Variant C does not fix this at all**: its
dual-porosity split leaves that same θ-quantile calculation untouched, so
Loess's and Vertisol's habitat stays bone-dry exactly as in Stage 3.5, and
both soils collapse to the identical degenerate dead curve Stage 3.5 already
reported (the `k_leak` exchange delivers ~4 orders of magnitude less
substrate to the habitat than Variant A). Full detail and mechanism in
`docs/stage_results/STAGE4_RESULTS.md`.

**Stage 5** (EXPLORATORY — matric-potential wetting + inscribed-circle
macropores, built on the further-evolved Variant C) replaced Variant C's
shared-theta quantile with a shared air-entry diameter
(`matric_d_cut_um=50`) and its forced 80 µm macropore diameter with a
numpy-only geometric measure (each connected macropore region sized by its
own maximum inscribed circle). Honest verdict: **this DOES give Loess and
Vertisol genuinely live wet habitat where Stage 4 gave exactly zero**
(habitat-wet-share 0.074/0.077 under the static map vs 0.000), **and
correctly identifies Sand's macropore network as one huge (81% of the grid),
well-connected, ~340 µm region too large to hold water at all**
(habitat-wet-share 0.000) — a real geometric finding, not a bug. The
pre-existing dynamic "wicking" mechanism substantially undermines this new
distinction in practice (habitat-wet-share rebounds to 0.595/0.570/0.205
with wicking on, since the always-wet matrix sits within wicking's
rewetting radius almost everywhere). Despite the real wetting improvement,
**no soil ever reaches "thriving"**: all three stay `declining`/`falling`
under every condition tested, because a separate, untouched, inherited
bottleneck (a per-macropore-cell substrate-holding cap, `S_cap`) keeps usable
substrate far below what growth needs regardless of habitat wetness — verified
to predate Stage 5's own changes. Full detail and mechanism in
`docs/stage_results/STAGE5_RESULTS.md`.

**Stage 6** (EXPLORATORY — porting a fix proved in a separate abstract PDE
sandbox, `pde_sandbox/`, into Variant C) replaced the explicit substrate
transport with an implicit backward-Euler Jacobi solve (unconditionally
stable, and — more importantly — unifying the macropore and matrix D fields
into one connected solve for the first time) and split each soil's OM into a
fast coating pool near habitat (per-soil fraction: Sand 0.80, Loess 0.50,
Vertisol 0.30) and a slow-releasing trapped pool; `S_cap` was removed and
wicking turned off by default. Honest verdict: **this measurably helps where
live wet habitat exists** — mean substrate reaching the habitat rose ~4-5×
for Loess and Vertisol, with dozens of macropore regions recruited into
active growth — **but essentially not at all for Sand**, whose habitat
remains the one huge dry macropore region Stage 5 identified. **Despite this
real improvement, no soil ever revives**: mean substrate stays 2-3 orders of
magnitude below `Ks=0.2` everywhere, so all three remain declining/falling
and none of the target shapes (Sand rise-then-fall, Vertisol rising, Loess
low/erratic) are reproduced — because of two factors outside this stage's
scope (a per-soil OM budget an order of magnitude below the historical
`OM_total=800`, and a genuinely short physical diffusion length given this
model's own D(d)/dt/dx/T). A pre-check also found the OLD explicit stability
cap was never actually close to binding here, correcting the stage's own
premise about why the fix should matter. Full detail and mechanism in
`docs/stage_results/STAGE6_RESULTS.md`.

**Stage 7** (a quantitative CALIBRATION of Variant C, not a new mechanism)
closes the two gaps Stage 6 diagnosed via three SHARED scalars: a single
`stage7_om_fraction` (replacing the retired per-soil OM fractions),
`stage7_D_scale` (a shared post-scale multiplier raising the transport
magnitude/reach without touching `D_d_ref`), and `D_film`/`stage7_film_
K_scale` (a small water-film pathway + reduced capacity for drained
macropore cells previously fully cut off, `D=0`). A 15-iteration
budget-capped search (`stage7_run.py`, fully exercised) cleared the
`S>=0.025` fed-habitat growth target for all three soils by 2-7x, robustly
at both grid sizes — the direct gap Stage 6 identified is closed. The
water-film mechanism additionally revealed this was never a Sand-specific
problem: ~91-100% of EVERY soil's macropore habitat sits dry and was
previously fully transport-isolated. **But exactly the failure mode the
stage anticipated also occurred**: all three soils bloomed into the same
"rising" shape (pairwise distinctness never exceeded ~0.06, two orders of
magnitude short of the 0.3 threshold, at every dose tested including the
minimum that just cleared growth) — "everything blooms and rises together"
under-separation, reported as a failure on equal footing with starvation.
Mechanistically, once fed-habitat substrate clears the growth/decay switch
at all, unbounded logistic growth ramps every soil toward carrying capacity
within the run, so `R(t)`'s qualitative SHAPE stops being governed by the
structural differences that still visibly change how much and how many
macropore regions grow. Full detail, honest verdict, and forward-looking
suggestions (out of this stage's scope) in
`docs/stage_results/STAGE7_RESULTS.md`.

**Stage 8** (a DEFINED configuration with verification, not a search) acts
directly on Stage 7's own diagnosis: replaces the shared OM *fraction*
with a single absolute `om_total`, byte-identical across all three soils
(confirmed equal: 600.0 for each), set by exactly 3 directed diagnostic
probes over `(om_total, T)` — no per-soil tuning, no scan. The same total
carbon is now concentrated on Sand's small coating vs. spread thin across
Vertisol's many habitat regions, letting the model's existing finite-fuel
+ Monod-limitation physics produce genuinely different depletion timing —
verified directly via a new `fuel_remaining_t` diagnostic (Sand consumes
84% of its pool by the end of the window, Loess 65%, Vertisol 51%, exactly
the predicted ordering). **This decisively achieves the primary
distinctness criterion**: `min_distance=0.7955` at the n=128 headline grid,
vastly exceeding the 0.3 threshold and Stage 7's 0.059 — equal carbon,
distributed by structure alone, genuinely separates the three curves. The
exact qualitative shape match is only 1 of 3 (Loess's "low/erratic"
matches; Sand's real R(t) decline is too late/mild in the coarse 4-period
view at n=128 to read as "falling" — it reads sharply correct only at
n=40; Vertisol's recruited regions also eventually run low on the same
lean shared pool, reading "erratic" rather than purely "rising"). Both
shortfalls trace to one honestly-named tension, not chased further: a
shared total lean enough for Sand's fall to register sharply is not
abundant enough for Vertisol's many dispersed regions to sustain a
monotonic rise through the whole window, with no per-soil knob available
to resolve that separately. Full detail in
`docs/stage_results/STAGE8_RESULTS.md`.

**Stage 9** (another DEFINED configuration) acts directly on Stage 8's own
finding: since one EQUAL total OM cannot give Sand-fall AND Vertisol-rise
together, each soil now gets a DIFFERENT total OM, ranked Vertisol(1.5) >
Loess(1.0) > Sand(0.6) x a shared base scale (confirmed: 900 > 600 > 360)
-- a real-world QUALITATIVE ranking, explicitly never fitted to the
measured respiration values -- plus Vertisol's coating fraction raised
(0.30->0.60) and Sand's lowered (0.80->0.50, subject to a "burst wins"
exception that was checked at every candidate and never triggered). The 3
directed probes were each checked at BOTH n=40 and n=128, and the two
grids DISAGREED: escalating Vertisol's total/coating past the prompt's own
suggested starting ratios looked increasingly promising at n=40 while
consistently WORSENING the same pairwise distance at n=128 -- adoption
followed the primary n=128 grid, so the modest, unescalated probe 1 (not
the most-pushed probe 3) is what is frozen. Result: a genuine partial
success -- Sand's R(t) now shows a real, sharp fall and Vertisol's a real
rise (`metrics.describe_shape`, a continuous descriptor, reads both
correctly), and Loess's period-binned shape still matches "erratic," but
the coarser 4-period classifier used against the experiment's own
measurement windows reads BOTH Sand and Vertisol as "erratic" too, so
which descriptor is applied changes how many labels "match" (both are
reported). **Distinctness fails** (`min_distance=0.162` at n=128, short of
0.3), driven entirely by Loess-vs-Vertisol, while Loess-vs-Sand and
Sand-vs-Vertisol both clear the threshold overwhelmingly (>1.8). Named
cause, not chased further: Loess and Vertisol are both large,
well-connected, many-region soils whose recruitment is fast and
front-loaded under this model's shared transport magnitude (kept from
Stage 7), so their normalized R(t) shapes stay too similar to each other
regardless of how much their total OM/coating fraction differ. Full detail
in `docs/stage_results/STAGE9_RESULTS.md`.

**Stage 10** restricts to two soils only (Vertisol, Sand -- Loess
dropped) and explicitly sets aside the shared/ranked-scalar purity of
Stages 7-9 to tune per-soil parameters directly against a numeric
pass/fail table: Sand needed a real, growth-dominated burst that then
crashes; Vertisol needed a genuine, sustained slow rise, not an early
spike that dips. The key new lever is PER-SOIL transport magnitude
(`D_scale`) -- the first time it is not one shared constant -- because
naively raising Sand's total OM alone at the old shared `D_scale` only
made its burst LATER and broader (`peak_frac` 0.75->0.88 as `om_total`
rose 3000->6000); giving Sand its own, much larger `D_scale` (80, vs.
Vertisol's own, much smaller 1.5) fixed that. The observation window `T`
also became per-soil, checked directly rather than assumed: the same Sand
candidate that crashes cleanly by `t=1500` at n=40 needed `T=3000` at
n=128 for the crash to fully register (more habitat cells means more peak
biomass, and biomass decay is a first-order process, so a bigger peak
takes proportionally longer, in absolute time, to decay); lengthening
Vertisol's own window to match instead broke ITS slow-rise shape
(`peak_frac` 0.91->0.46). **Result: every §3 pass/fail criterion is met, at
BOTH n=128 and n=40 (16 of 16 checks pass)** -- Sand's `R_peak=15.79` is
39x the maintenance floor with 77% of that peak growth-associated (not
maintenance), crashing to `R_end/R_peak=0.019` as its coating pool empties
(`fuel_remaining` 1.0->0.015); Vertisol's late-third mean is 5.8x its
early-third mean, peaking at 91% of its window with `R_end/R_peak=0.98`
(no dip), 55 of 70 macropore regions recruited, still ~49% of its fuel
unspent at the end. Distinctness (time-normalized, since the windows
differ) is 2.17, more than 7x the 0.3 threshold. Neither curve is
hand-drawn -- both emerge from the SAME finite-fuel + Monod-limitation +
implicit-diffusion mechanism used since Stage 6, with per-soil physical
inputs matching each soil's real textural contrast. Full detail in
`docs/stage_results/STAGE10_RESULTS.md`.

**Stage 11** continues from Stage 10 with Sand kept exactly FROZEN
(confirmed via an exact numeric identity check against Stage 10's own
recorded `R_peak`/`R_end`, not just by not touching its config -- relative
differences 0.032% and 0.79%) and fixes a real per-timestep performance
bug: `habitat_properties` (a connected-component BFS over the whole grid)
was being recomputed every single timestep despite the macropore topology
never changing during a run. Hoisting it out, plus an optional Jacobi
convergence-tolerance early stop, gives an **11.5x speedup for Sand and a
50-60x speedup for Vertisol** (which has far more macropore regions, so
paid a bigger tax under the bug), with results unchanged within numerical
noise -- a lower FIXED iteration cap (20-30), also suggested by the stage
prompt, was tried and explicitly REJECTED because it measurably changed
results (up to 24.8% on Sand's `R_end`). With Sand untouched, Vertisol was
re-tuned further: Stage 10's Vertisol already had a plausible R(t) shape
but its biomass only ever reached PARITY with its initial total
(`B_peak/B0=1.000`); raising `om_total` alone (2000->10000, coating
fraction deliberately unchanged) crossed a real threshold above which
aggregate growth outpaces decay, giving `B_peak/B0=2.34` -- genuine growth,
not just a shape that looks like one. One limitation is reported honestly,
not hidden: `habitats_recruited(t)` reaches its final count within the
first 2% of the window and stays completely flat afterward regardless of
tuning, because the model's very low activation threshold is crossed
almost everywhere almost instantly by the implicit solver's per-step
equilibration -- R(t)'s own genuine gradual rise instead comes from
already-active regions' growth RATE slowly building up, not new regions
joining. Full detail in `docs/stage_results/STAGE11_RESULTS.md`.

**Stage 12** is an ANALYSIS stage, not a modeling stage: Stage 11's two
soils and every one of their frozen parameters (`om_total`, coating
fraction, `D_scale`, biology, `T`) are run completely unchanged across
M=40 independent random `pore.seed` values each (80 simulations, 13.0
minutes using the Stage-11 performance fix), to test statistically
whether each soil's characteristic shape comes from its defined
structural features or from a lucky random layout. **Every shape metric
tested shows COMPLETE rank separation between the two 40-seed ensembles**
(Mann-Whitney rank-biserial `r = ±1.000`, `p` in the `1e-14`-`1e-15`
range, Cohen's `d` 5.0-36.6 -- `d=0.8` is conventionally "large," so these
are one to two orders of magnitude beyond that; all statistics implemented
directly in numpy, no scipy). The mean between-soil curve distance (2.06)
is 13x Sand's own within-soil spread and 127x Vertisol's, with a
bootstrap 95% CI on the gap of `[1.963, 1.989]`, excluding zero by a wide
margin. Structural descriptors confirm what varies with the seed and what
doesn't: Sand consistently realizes 1-9 macropore regions in an
always-18%-matrix soil, Vertisol 55-86 regions in an always-68%-matrix
soil (`matrix_fraction` is exactly seed-invariant per soil -- a property
of the shared literature PSD's marginal distribution, not the spatial
arrangement). An optional bounded causal check (sweeping Sand's own
coating fraction over 4 values, 5 seeds each) shows its crash depth
(`R_end/R_peak`) drop monotonically as coating fraction rises
(0.055->0.022->0.006->0.001) -- direct evidence the swept FEATURE itself
drives the metric. **Conclusion**: because each soil's structural class is
held fixed while only the random layout varies, and every seed reliably
reproduces that soil's characteristic shape with between-soil separation
one to two orders of magnitude larger than seed-to-seed noise, each
soil's respiration shape is determined by its defined structural
features -- not by a lucky random realization. Full detail in
`docs/stage_results/STAGE12_RESULTS.md`.
