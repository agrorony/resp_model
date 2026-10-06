# Report-generation prompt (run this later; do not run now)

You are writing a structured scientific report, in English, suitable for a
course appendix, about a soil-structure microbiome respiration simulation.
The report covers Stage 1 (fully-wet pore-based transport), Stage 2
(water-filled pores / percolation), Stage 3 (real measured PSD, physiological
carrying capacity, bigger grids), Stage 3.5 (literature-informed parametric
PSDs for the three real soil types Loess / Sand / Vertisol, tested against
the measured experiment at constant moisture), Stage 4 (EXPLORATORY: two
head-to-head fixes -- PSD truncation + banded OM vs. dual-porosity
matrix/macropore exchange -- for Stage 3.5's "OM never reaches the habitat"
diagnosis), Stage 5 (EXPLORATORY: replacing Variant C's shared-theta
wetting with a matric potential and its forced macropore diameter with a
geometric inscribed-circle size, to test whether that alone revives
Loess/Vertisol's dead habitat), Stage 6 (EXPLORATORY: porting an
implicit-transport + coating/trapped-OM recipe, proved in a separate
abstract PDE sandbox, into Variant C, to test whether that revives growth),
Stage 7 (a quantitative CALIBRATION, not a new mechanism: a single
shared OM fraction, a shared reach/diffusion-magnitude scale, and a
water-film pathway for drained macropores, budget-capped-searched to close
the ~25x substrate-budget/reach gap Stage 6 diagnosed), Stage 8 (a
DEFINED configuration, not a search: replacing Stage 7's shared OM
*fraction* with a single absolute OM total, byte-identical across soils,
set by at most 3 directed diagnostic probes, to test whether equal carbon
concentrated vs. spread by structure produces different DEPLETION TIMING
and therefore genuinely distinct respiration shapes), Stage 9 (another
DEFINED configuration: since a single EQUAL total OM could not give
Sand-fall AND Vertisol-rise together, giving each soil a different total OM
ranked Vertisol > Loess > Sand, plus Vertisol's coating fraction raised and
Sand's lowered -- all ranking-only, never fitted to the measured
respiration values, set by at most 3 directed probes), Stage 10 (two soils
only, Vertisol and Sand -- Loess dropped: per-soil `om_total`, coating
fraction, and, for the first time, per-soil transport magnitude and
observation window, tuned directly to verify a real growth-driven
burst-then-crash for Sand and a genuine sustained slow rise for Vertisol
against an explicit, numeric pass/fail table), Stage 11 (a follow-up
to Stage 10: Sand kept FROZEN exactly as Stage 10 left it, a real
per-timestep performance bug fixed in the simulation loop and verified to
change nothing, and Vertisol re-tuned further so its biomass clearly
GROWS, not merely reaches parity with its initial total), and Stage 12 (an
ANALYSIS stage, not a modeling stage: Stage 11's two soils and every one
of their frozen parameters are run unchanged across M=40 random
pore-field seeds each, to show statistically -- Mann-Whitney U, effect
size, and bootstrap confidence intervals, all implemented in numpy only --
that each soil's characteristic respiration shape is driven by its defined
structural features, not by a lucky random realization), and Stage 13 (a
DATA-SOURCING + RECALIBRATION-CHECK stage, not a new mechanism: Stage
3.5-12's literature-parametrized PSDs are replaced with real measured PSD +
pore-network connectivity data from a companion CT-scan/segmentation
research project, and Stage 10/11's frozen per-soil configuration is
re-run unchanged against the real data -- both soils' pass/fail tables
fail completely, traced to a structural cause (the real PSDs carry zero
pore volume below the model's macro/matrix split threshold) rather than a
tuning problem, confirmed rather than fixed by a budget-capped search over
the same allowed knobs), Stage 14 (a DIAGNOSTIC stage, not a retuning
stage: instruments the simulation to record the actual per-cell growth
VALUE distribution over time, rather than only the boolean "recruited"
count Stage 11 found saturates within the first 2% of the observation
window, to directly confirm or refute the suspected cause -- confirmed
directly: the implicit solver's per-step equilibration crosses the
activation threshold almost everywhere almost instantly, while a graded,
growth-weighted alternative metric tracks R(t)'s own genuine gradual rise
almost exactly, unlike the boolean count -- an attribution Stage 16 later
partly overturned, see below), Stage 16 (retiring the per-soil coating/
trapped OM split for one shared, structure-blind placement rule with
uniform dissolution: shapes and distinctness survive the change, but
recruitment timing fails structurally because the emergent fuel-habitat
distances are only 1-2 grid cells in both soils), and Stage 15 (reintegrating
Loess for a genuine three-soil demonstration on real measured data, Sand
and Vertisol deliberately kept at their last working literature-PSD
configuration since Stage 13 showed real data breaks their mechanism
entirely -- a budget-capped per-soil search for Loess achieves robust
three-soil distinctness and a genuinely low, growth-driven Loess curve,
but the measured experiment's "erratic" (non-monotonic) target is shown,
directly and mechanistically, to be unreachable by this model's smooth,
monotonically-forcing dynamics, on real or literature PSD alike).

Note on layout: the historical per-stage prompts live in `prompts/`, the
per-stage narrative writeups in `docs/stage_results/`, and each stage's figures
in `results/<stage>/` (`baseline/`, `stage2/`, `stage3/`, `stage35/`).

Read these inputs from the project folder before writing anything:

1. `MODEL_SPEC.md` — the model equations, frozen biology constants,
   structure knobs, the pore -> K/D mapping, the `global_theta` saturation
   rule (SS9), and the Stage 1/2/3 result summaries (SS11, SS14, SS15).
2. `LOGBOOK.md` — the complete, honest record of the Stage-1 structure
   search, the Stage-2 numpy-only conversion + theta sweep, AND the Stage-3
   real-PSD conversion + theta x grid sweep + hypothesis test: what was
   tried, the resulting metrics, and what was changed and why.
3. `docs/stage_results/STAGE1_RESULTS.md`, `STAGE2_RESULTS.md`,
   `STAGE3_RESULTS.md`, and `STAGE3_5_RESULTS.md` — the per-stage narrative
   writeups (idea, method, what emerged, honest verdict, limitations).
4. Stage 1 figures/data in `results/`:
   - `respiration_curves.png`, `cumulative_co2.png` — R(t) / cumulative CO2
     at full saturation, for soils A, B, C on one axis.
   - `structure_maps_A.png`, `structure_maps_B.png`, `structure_maps_C.png`
     — pore radius, K, and D maps per soil.
   - `soil_A_timeseries.csv`, `soil_B_timeseries.csv`, `soil_C_timeseries.csv`
     — raw `time, R_t, cumulative_co2` per soil at full saturation.
5. Stage 2 figures/data in `results/`:
   - `stage2_water_map_{A,B,C}_theta0.40.png` — pore radius, K, D, and the
     water-filled mask at theta=0.4.
   - `stage2_respiration_curves_theta0.40.png` — R(t) at theta=0.4.
   - `stage2_respiration_vs_theta.png` — peak R / cumulative CO2 vs theta.
   - `stage2_connectivity_vs_theta.png` — largest water cluster fraction,
     OM-biomass connectivity, and accessible-OM fraction vs theta.
   - `stage2_theta_sweep.csv`, `stage2_distinctness_by_theta.csv` — raw
     numbers behind every Stage-2 figure.
6. Stage 3 figures/data in `results/`:
   - `stage3_psd_loaded.png` — the three soils' real, loaded empirical PSDs
     (density and CDF).
   - `stage3_water_map_{A,B,C}_n{128,40}_theta0.40.png` — pore diameter,
     K(d), D, and the water-filled mask at theta=0.4, both grid sizes.
   - `stage3_respiration_curves_theta0.40_n{128,40}.png` — R(t) at
     theta=0.4, both grid sizes.
   - `stage3_respiration_vs_theta.png` — peak R / cumulative CO2 vs theta,
     both grid sizes overlaid (solid=n128, dashed=n40).
   - `stage3_connectivity_vs_theta.png`, `stage3_isolation_vs_theta.png` —
     the Stage-2 connectivity diagnostics plus the new isolation ratio, both
     grid sizes.
   - `stage3_outcome_vs_isolation_ratio.png`,
     `stage3_outcome_vs_psd_stat.png` — the SS9 hypothesis-test scatter
     plots (cumulative CO2 vs isolation ratio / vs PSD habitable fraction).
   - `stage3_theta_sweep.csv`, `stage3_distinctness_by_theta.csv`,
     `stage3_isolation_vs_psd_correlation.csv` — raw numbers behind every
     Stage-3 figure.

Write `REPORT.md` (or `REPORT.pdf` if a PDF pipeline is available) with the
following sections:

1. **Motivation** — why spatial structure alone, with identical biology,
   can produce qualitatively different respiration dynamics; the practical
   relevance to real soils (pore structure, substrate accessibility, water
   retention/percolation).
2. **Model and equations** — reproduce the equations from `MODEL_SPEC.md`
   (biomass, substrate, organic matter, respiration observable, the pore ->
   K/D maps, and the `global_theta` water-mask rule), explained in prose,
   with the frozen constants tabulated.
3. **The three soils and their structure** — for each soil, describe the
   pore-field texture/structure parameters, the resulting K and D layout,
   and organic-matter placement, embedding the corresponding
   `structure_maps_*.png` figure.
4. **Stage 1 results (fully wet)** — embed `respiration_curves.png` and
   `cumulative_co2.png`; report the pairwise distances and shapes from
   `STAGE1_RESULTS.md`; state the honest A-C caveat (robust seed-resampling
   showed A-C is not a reliably separated pair at full saturation).
5. **Stage 2 results (percolation)** — explain the `global_theta` mechanism
   (small pores fill first; a soil's own texture/arrangement determines how
   its water network fragments as theta drops). Embed
   `stage2_water_map_{A,B,C}_theta0.40.png`,
   `stage2_respiration_curves_theta0.40.png`,
   `stage2_respiration_vs_theta.png`, and `stage2_connectivity_vs_theta.png`.
   Report the theta-by-theta distinctness table and connectivity numbers
   from `STAGE2_RESULTS.md` SS3, including where success was and was not
   achieved and why (soil C's biomass die-off at low theta, the AC pair's
   persistent fragility).
6. **Stage 3 results (real PSD, physiological K, bigger grid)** — explain
   the empirical-quantile mapping (real measured PSD replaces the lognormal
   field), the physiological `K(d)` window, and why `n=128` vs `n=40`
   matters (a bigger domain lets more independent correlation lengths and
   isolated regions form; `n=40` at `lambda=2.0` is close enough to
   spatially homogeneous that cumulative CO2 comes out nearly flat across
   soil and theta). Embed `stage3_psd_loaded.png`,
   `stage3_water_map_{A,B,C}_n128_theta0.40.png`,
   `stage3_respiration_vs_theta.png`, `stage3_isolation_vs_theta.png`, and
   the two hypothesis-test scatter plots. Report the theta x grid
   distinctness table and the SS9 correlation numbers from
   `STAGE3_RESULTS.md`, including: which pairs separated and which did not
   (the BC pair never approaches the 0.3 threshold, mechanistically because
   both soils' large median diameters saturate `D(d)`); that no soil died
   off in this sweep (contrast with Stage 2's soil-C die-off, and why —
   `aggregate=True` is off by Stage-3 invariant); and that the isolation
   ratio predicts respiration outcome better than the PSD summary stats
   tested, but only weakly (`|r|` ~0.25-0.33).
7. **Interpretation** — explain, mechanistically, why each structure
   behaves as it does at full saturation, under drying, with the real PSD,
   and across grid sizes: diffusion-front arrival time relative to the
   observation window; OM depletion; percolation fragmentation vs a soil's
   spatial correlation length (`lambda`) and `aggregate` flag (Stage 1/2)
   or PSD shape and `D(d)` saturation (Stage 3). Tie this explicitly to the
   equations and the connectivity/isolation diagnostics, not just the
   qualitative description.
8. **Limitations** — combine Stage 1's caveats (no calibration to real
   data, single-realization seeds, the shape classifier's empirically
   chosen thresholds), Stage 2's (the theta grid is coarse; the
   high-OM/high-K connectivity threshold — top quartile — is a documented
   choice, not derived; no Stage-2 parameter search was run, only a theta
   sweep at the existing Stage-1 pore parameters), and Stage 3's (the real
   PSD's spatial *arrangement* is still synthetic and shared across soils;
   only one shared `(lambda, aggregate)` recipe was tried, not swept; `K`
   window and `D` reference constants are documented choices, not
   calibrated; a seed-ensemble significance pass remains deferred
   throughout).
9. **Budget note** — state plainly, from `LOGBOOK.md`, that the Stage-1
   structure search succeeded at iteration 1, that Stage 2 used a direct
   theta sweep (no nudge-loop iterations were needed), and that Stage 3
   used a direct theta x grid sweep (also no nudge-loop iterations —
   `MAX_ITERATIONS=12` was not exercised). Report the Stage-2 success rate
   across the sweep (2 of 5 theta values) and the Stage-3 result (no
   grid/theta combination reached full three-way success; the BC pair
   stayed fused throughout) honestly, without implying either pattern was
   monotonic or fully resolved.

Stage 3.5 inputs in `results/stage35/`:
   - `stage35_psd_loaded.png` — the three literature PSDs (density + CDF) with
     the 30–150 µm habitat window; `stage35_psd_summary.csv` — porosity,
     median diameter, habitable volume fraction per soil.
   - `stage35_D_of_d_fix.png` / `.csv` — the D(d) fusion fix (D_d_ref 90 → 250).
   - `stage35_experiment_overlay.png` — model vs measured temporal shape per
     soil (the headline validation figure); `stage35_constant_theta.csv`.
   - `stage35_experiment_overlay_matric.png` + `stage35_matric_diagnostic.csv`
     — the shared-matric-potential diagnostic.
   - `stage35_theta_sensitivity.png` / `.csv`, `stage35_theta_robustness.png`,
     `stage35_water_map_{loess,sand,vertisol}.png`, `stage35_distinctness.csv`.

When covering Stage 3.5, report the honest verdict: the model reproduces the
measured shape for **Sand only (falling), and not robustly** (it flips to rising
at n=40); it does NOT reproduce Vertisol-rising or Loess-erratic. Explain the
mechanism — at a shared θ the fine soils' habitat macropores are the last to
wet and stay dry, and even under a shared matric potential the frozen OM rule
concentrates substrate on sub-micron clay from which it cannot diffuse to the
habitat. State that the literature PSDs are representative of soil TYPES, not
measurements of the specific samples.

Stage 4 inputs in `results/stage4/`:
   - `stage4_variantA_truncation.csv` — per-soil porosity/habitable-fraction
     before/after the 10 µm PSD truncation.
   - `stage4_comparison_periods.png` — THE headline comparison figure: Variant
     A (top row) vs Variant C (bottom row) vs measured P1–P4 shape, three
     soils each.
   - `stage4_substrate_reaching_habitat.png` — the killer-metric bar chart:
     mean substrate reaching the habitat, both variants, all three soils.
   - `stage4_respiration_variant{A,C}.png` — raw R(t)/cumulative-CO2 curves
     per variant; `stage4_water_map_variant{A,C}_{loess,sand,vertisol}.png` —
     pore/K/D/water-mask maps (Variant C adds the macro/matrix split panel).
   - `stage4_variant{A,C}_sweep.csv` — full per-soil × theta × grid results
     (outcome, shape, habitat/wet-habitat fraction, isolation ratio, mean
     substrate reaching habitat); `stage4_distinctness.csv` — pairwise
     distances and success flag per variant/condition.

When covering Stage 4, report the honest verdict from
`docs/stage_results/STAGE4_RESULTS.md`: **Variant A (PSD truncation + banded
OM) reproduces 2 of 3 measured shapes (Loess-erratic, Vertisol-rising)
robustly across θ, though not the emergent-distinctness threshold and not
grid-robustly for Vertisol; Variant C (dual-porosity matrix/macropore
exchange) does not fix Stage 3.5's failure at all** — because Variant A's PSD
truncation has the (unplanned) side effect of shifting the shared-θ
quantile cutoff toward the coarse end, incidentally rewetting the habitat
that was bone-dry in Stage 3.5, while Variant C leaves that quantile
calculation completely untouched, so Loess's and Vertisol's habitat stays
bone-dry and both soils collapse to the same degenerate dead curve Stage 3.5
already reported. Explain that neither variant was tuned toward the target
(SS7 of the stage-4 prompt) and that this asymmetry — a PSD fix accidentally
solving a saturation problem, while an explicit dual-porosity fix does not —
is reported as-is, not resolved.

Stage 5 inputs in `results/stage5/`:
   - `stage5_comparison_periods.png` — the headline comparison figure: three
     soils' R(t) (period-binned), model with wicking ON (solid) and OFF
     (dashed) vs measured P1–P4 shape.
   - `stage5_region_map_{loess,sand,vertisol}.png` — per-soil macropore-region
     inscribed-diameter map, wet/dry classification, and final water mask.
   - `stage5_region_diameter_hist.png` — distribution of macropore-region
     inscribed diameters per soil, with the 50 µm matric cutoff marked.
   - `stage5_substrate_reaching_habitat.png` — mean substrate reaching the
     habitat, wicking ON vs OFF, all three soils.
   - `stage5_matric_sensitivity.png` / `.csv` — wet-habitat fraction and
     emergent θ vs `matric_d_cut_um` (40/50/60).
   - `stage5_headline.csv` — per-soil, per-wicking-setting: emergent θ,
     wet-habitat fraction, habitat-wet share, mean substrate reaching
     habitat, outcome, shape, region counts; `stage5_region_diameters.csv` —
     every macropore region's size and wet/dry status; `stage5_grid_check.csv`
     (n=40); `stage5_distinctness.csv`.

When covering Stage 5, report the honest verdict from
`docs/stage_results/STAGE5_RESULTS.md`: **the matric-potential wetting rule
plus geometric (inscribed-circle) macropore sizing DOES give Loess and
Vertisol genuinely live wet habitat where Stage 4's shared-theta quantile
gave exactly zero** (habitat-wet-share 0.074/0.077 under the static map,
vs 0.000), **and correctly identifies Sand's macropore network as one huge,
well-connected, ~340 µm region too large to hold water at all**
(habitat-wet-share 0.000) — a real geometric finding, not a bug, exactly the
"Sand's large drained macropores leave it too dry" outcome the stage-5 prompt
anticipated. Explain that the pre-existing dynamic "wicking" mechanism
(kept from Variant C's baseline) substantially undermines this new static
distinction by rewetting most of what it marks dry, since the always-wet
matrix sits within wicking's rewetting radius almost everywhere — habitat-wet
-share jumps back up to 0.595/0.570/0.205 with wicking on, confirming the
exact conflict the prompt asked to isolate. State plainly that despite this
real wetting improvement, **no soil ever reaches "thriving"**: all three stay
`declining`/`falling` under every condition tested, because a separate,
untouched, inherited bottleneck (the per-macropore-cell substrate-holding
cap, `S_cap`) caps usable substrate far below what growth needs regardless of
habitat wetness — verified to predate Stage 5's changes, and left untouched
per the stage-5 prompt's "keep everything else" instruction. Note
`voxel_um`, `matric_d_cut_um`, and the literature-PSD representativeness as
key assumptions, and that the inscribed-diameter measure is an L1/Manhattan
approximation ("diamond," not a true circle), per the prompt's own specified
numpy-only method.

Stage 6 inputs in `results/stage6/`:
   - `stage6_comparison_periods.png` — the headline comparison figure: model
     R(t) (period-binned) vs measured P1–P4 shape, three soils.
   - `stage6_mean_s_habitat_vs_t.png` — the killer-metric time series, with
     `Ks=0.2` marked for scale; `stage6_habitats_recruited_vs_t.png` — count
     of actively-growing macropore regions over time, per soil.
   - `stage6_om_map_{loess,sand,vertisol}.png` — macropore/matrix/solid
     phases, the coating mask, and log-density maps of `OM_coating` and
     `OM_trapped`.
   - `stage6_respiration_curves.png` — raw R(t)/cumulative-CO2 curves.
   - `stage6_headline.csv` — per-soil emergent θ, wet-habitat fraction, mean
     substrate reaching habitat (final + time-averaged), habitats-recruited
     stats, outcome, shape, coating fraction and cell counts;
     `stage6_grid_check.csv` (n=40); `stage6_distinctness.csv`.

When covering Stage 6, report the honest verdict from
`docs/stage_results/STAGE6_RESULTS.md`: two mechanisms proved in a separate
abstract PDE sandbox (`pde_sandbox/`) — implicit backward-Euler transport
and a coating/trapped OM split with per-soil coating fraction — were ported
directly into the real soil model's Variant C. State plainly, as a
correction to the stage's own premise, that the OLD explicit stability cap
was never actually close to binding in this model (measured margins
0.035–0.088, well under the 0.2 cap), so removing it could not by itself
explain any improvement; the real new mechanism is that Stage 4/5 ran
macropore and matrix transport as two non-overlapping explicit divergences
that never let those two phases diffuse into each other directly, and Stage
6 unifies them into one connected implicit solve. Report that this
measurably helped: mean substrate reaching the habitat rose ~4–5× for Loess
and Vertisol (with dozens of macropore regions recruited into active growth)
but essentially not at all for Sand (whose habitat remains the one huge dry
macropore region Stage 5 identified, leaving little wet-adjacent matrix for
its coating to sit next to despite its high 0.8 coating fraction). State
plainly that despite this real improvement, **no soil ever revives**: mean
substrate stays 2–3 orders of magnitude below `Ks=0.2` everywhere, so all
three remain declining/falling and none of the target shapes (Sand
rise-then-fall, Vertisol rising, Loess low/erratic) are reproduced. Explain
the two out-of-scope reasons why: a per-soil OM budget an order of magnitude
below the historical shared `OM_total=800`, and a genuinely short physical
diffusion length given this model's own D(d)/dt/dx/T, unlike the sandbox's
much finer relative grid resolution. Note that Sand has now been the
recurring problem soil since Stage 3.5.

Stage 7 inputs in `results/stage7/`:
   - `stage7_comparison_periods.png` — the headline comparison figure: model
     R(t) (period-binned) vs measured P1–P4 shape, three soils, with the
     calibrated shared parameters in the subtitle.
   - `stage7_mean_s_habitat_vs_t.png` — fed-habitat substrate over time, with
     the `S=0.025` growth target and `Ks=0.2` marked; `stage7_habitats_
     recruited_vs_t.png` — count of actively-growing macropore regions.
   - `stage7_respiration_curves.png` — raw R(t)/cumulative-CO2 curves.
   - `stage7_calibration_history.png` / `.csv` — the two scoring axes (min
     fed-habitat S, min pairwise distance) across all 15 search iterations.
   - `stage7_headline.csv` — per-soil `om_total`, fed-habitat S/f_S, habitats
     recruited, `dry_macro_fraction_of_habitat`, outcome, shape, coating
     fraction; `stage7_grid_check.csv` (n=40); `stage7_distinctness.csv`.

When covering Stage 7, report the honest verdict from
`docs/stage_results/STAGE7_RESULTS.md`: a 15-iteration budget-capped search
over three SHARED scalars — a single OM fraction (replacing the per-soil
fractions), a post-scale multiplier on the transport magnitude, and a
water-film diffusivity for drained (dry) macropore cells that were
previously fully cut off (`D=0`) — closed the substrate gap **robustly**:
fed-habitat substrate cleared the `S>=0.025` growth target for all three
soils by 2–7x, at both grid sizes, and the water-film mechanism revealed
that ~91–100% of EVERY soil's habitat (not only Sand's single giant dry
region) had been completely cut off from transport before this stage.
State plainly, per the stage's own SS9 instruction, that this success comes
paired with the anticipated failure: **all three soils bloomed into the
same "rising" shape** (none of Sand rise-then-fall / Vertisol rising /
Loess low-erratic is reproduced except Vertisol's rising by coincidence),
pairwise distinctness stayed two orders of magnitude short of the 0.3
threshold throughout the entire 15-iteration search (0.027–0.032) and
across an independent, more fine-grained manual sweep at doses down to the
minimum that just cleared growth — this is not an overshoot artifact, it
is the outcome as soon as growth clears at all anywhere. Explain the
mechanism: once fed-habitat substrate clears the growth/decay switch, the
model's unbounded logistic growth (uncontested by any other soil-varying
growth-RATE mechanism) ramps every soil toward its carrying capacity within
the run, so `R(t)`'s qualitative SHAPE stops being governed by the
per-soil structural differences (coating fraction, macropore topology) that
still visibly change how MUCH and how MANY regions grow. State the
forward-looking conclusion Stage 7 draws but does not attempt (out of
scope per its own SS1): reproducing the measured qualitative shapes will
need something that differentiates soils' growth DYNAMICS or TIMING, not
merely how much carbon reaches them — e.g. a depletion-sensitive coating
release, or a structure-tied (not shared) growth-rate/capacity modifier.

Stage 8 inputs in `results/stage8/`:
   - `stage8_comparison_periods.png` — the headline comparison figure: model
     R(t) (period-binned) vs measured P1–P4 shape, three soils, with the
     chosen `(om_total, T)` in the subtitle.
   - `stage8_mean_s_habitat_vs_t.png` — fed-habitat substrate over time,
     with the `S=0.025` growth target marked; `stage8_fuel_remaining_vs_t.png`
     — fraction of each soil's (equal) initial total OM still undissolved,
     the direct depletion-timing evidence; `stage8_habitats_recruited_vs_t.png`
     — count of actively-growing macropore regions over time.
   - `stage8_respiration_curves.png` — raw R(t)/cumulative-CO2 curves.
   - `stage8_headline.csv` — per-soil `om_total_init` (confirming equality),
     fed-habitat S (final and peak, with peak timing), fuel-remaining,
     habitats recruited, outcome, shape, `R_end/R_peak`; `stage8_grid_check.csv`
     (n=40); `stage8_distinctness.csv`. The 3-probe search trace (candidate
     `(om_total, T)`, all diagnostics, and reasoning notes) is in
     `LOGBOOK.md`'s "Stage 8 directed probe" entries, not a separate CSV.

When covering Stage 8, report the honest verdict from
`docs/stage_results/STAGE8_RESULTS.md`: replacing Stage 7's shared OM
*fraction* with a single ABSOLUTE total, byte-identical across all three
soils (confirmed: `om_total_init=600.0` for all three), set by exactly 3
directed diagnostic probes (not a search), lets the model's existing
finite-fuel + Monod-limitation physics produce a REAL, soil-ordered
depletion — verified directly via the new `fuel_remaining_t` diagnostic
(Sand 84% of its pool consumed by the end of the window, Loess 65%,
Vertisol 51%), not merely inferred from R(t) shape. State plainly that this
achieves the PRIMARY success criterion decisively: pairwise distinctness at
the n=128 headline grid is `min_distance=0.7955`, vastly exceeding the 0.3
threshold and Stage 7's 0.059 — equal carbon, distributed by structure,
genuinely separates the three curves. But the exact QUALITATIVE shape match
is only 1 of 3 (Loess's "low/erratic" matches; Sand shows a real R(t)
decline from its peak but too late and mild within the coarse 4-period
binning at n=128 to read as "falling" rather than "rising" — it reads
sharply correct only at the smaller n=40 grid, a genuine grid sensitivity;
Vertisol's ~55-70 recruited macropore regions also eventually run low on
the SAME lean shared pool, so its curve reads "erratic" rather than the
hoped-for purely sustained "rising"). Explain the single root tension named
in SS6 rather than chased with a further probe: a shared total lean enough
to make Sand's fall register sharply is not abundant enough for Vertisol's
many dispersed regions to sustain a monotonic rise through the whole
window — with no per-soil knob available, by this stage's own design, to
resolve that separately for each soil.

Stage 9 inputs in `results/stage9/`:
   - `stage9_comparison_periods.png` — the headline comparison figure: model
     R(t) (period-binned) vs measured P1–P4 shape, three soils, with the
     ranked OM/coating values in the subtitle.
   - `stage9_mean_s_habitat_vs_t.png`, `stage9_fuel_remaining_vs_t.png`
     (each soil's OWN, now-unequal total OM), `stage9_habitats_recruited_vs_t.png`
     — the depletion-timing diagnostics.
   - `stage9_respiration_curves.png` — raw R(t)/cumulative-CO2 curves.
   - `stage9_headline.csv` — per-soil `om_total_init` (confirming the V>L>S
     ranking), `coating_fraction`, fed-habitat S, fuel-remaining, habitats
     recruited, outcome, BOTH shape descriptors (`curve_shape` = continuous,
     `period_shape_model` = P1–P4 binned), `R_end_over_peak`;
     `stage9_grid_check.csv` (n=40); `stage9_distinctness.csv`. The 3-probe
     trace (each probe checked at BOTH n=40 and n=128) is in `LOGBOOK.md`'s
     "Stage 9 directed probe" entries.

When covering Stage 9, report the honest verdict from
`docs/stage_results/STAGE9_RESULTS.md`: since Stage 8 showed a single EQUAL
total OM cannot give Sand-fall AND Vertisol-rise together (opposite
regimes), this stage gives each soil a DIFFERENT total OM, ranked
Vertisol(1.5) > Loess(1.0) > Sand(0.6) x a shared base scale (confirmed:
900 > 600 > 360) — a real-world QUALITATIVE ranking, explicitly never
derived from the experiment's own P1–P4 numbers — plus Vertisol's coating
fraction raised (0.30→0.60, more of its now-larger carbon reaches its many
habitat regions) and Sand's lowered (0.80→0.50, smaller fast burst + a
slow trapped tail), subject to a "burst wins" exception that was checked at
every candidate and never triggered. State the METHODOLOGICAL finding
plainly: the 3 directed probes were each checked at BOTH the fast n=40 grid
and the n=128 primary grid, and the two DISAGREED — escalating Vertisol's
total/coating past the prompt's own suggested starting ratios looked
increasingly promising at n=40 (Loess-vs-Vertisol distance rising) while
consistently WORSENING the same pair at n=128; per the stage's own
instruction that adoption must reflect the grid the success criterion is
evaluated against, the modest, unescalated probe 1 — not the most-pushed
probe 3 — is the actual best n=128 candidate and what is frozen. Report the
result as a genuine partial success: Sand's R(t) now shows a real, sharp
fall and Vertisol's a real rise (`metrics.describe_shape`, a continuous
descriptor, reads both correctly) and Loess's period-binned shape still
matches "erratic" — but the coarser 4-period classifier used against the
experiment's own measurement windows reads BOTH Sand and Vertisol as
"erratic" too, so which descriptor is applied changes how many labels
"match"; report both, not just the flattering one. State plainly that
DISTINCTNESS FAILS (`min_distance=0.162` at n=128, short of 0.3), driven
entirely by Loess-vs-Vertisol, while Loess-vs-Sand and Sand-vs-Vertisol both
clear the threshold overwhelmingly (>1.8). Name the specific cause per the
stage's own instruction rather than chasing a 4th probe: Loess and Vertisol
are both large, well-connected, many-region soils whose recruitment is fast
and front-loaded under this model's shared transport magnitude (kept from
Stage 7), so their normalized R(t) shapes stay too similar to each other
regardless of how much their total OM/coating fraction differ — and pushing
that gap further made the primary-grid separation worse, not better.

Stage 10 (FINAL) inputs in `results/stage10/`:
   - `stage10_respiration_curves.png` — the headline two-soil figure: R(t)
     for Sand and Vertisol on both an absolute-time axis (own windows) and
     a normalized (0=start, 1=end of own window) axis.
   - `stage10_sand_growth_vs_maintenance.png` — the direct evidence Sand's
     peak is growth-driven: a stackplot of the growth-associated and
     maintenance respiration terms over time, against R(t) and the
     maintenance floor.
   - `stage10_mean_s_habitat_and_fuel.png`, `stage10_biomass_vs_t.png` —
     `mean_S_habitat(t)`, `fuel_remaining(t)`, and `biomass(t)`/B0 for both
     soils.
   - `stage10_headline.csv` — per-soil parameters (`om_total`, coating
     fraction, `D_scale`, `T`) and every §3 pass/fail check with its
     numeric value; `stage10_grid_check.csv` (n=40); `stage10_
     distinctness.csv` (time-normalized, since the two soils use different
     windows).

When covering Stage 10, present it as the CONSOLIDATED FINAL DEMONSTRATION,
not another exploratory stage. State plainly that Stage 10 restricts to two
soils (Vertisol, Sand; Loess dropped) and — for the first time — tunes
per-soil parameters directly against an explicit numeric pass/fail table
(§3 of the stage prompt) rather than only structural/ranking-only levers,
including a genuinely new lever: PER-SOIL transport magnitude (`D_scale`)
and observation window (`T`), breaking the "one shared constant" invariant
every prior stage held. Explain why both had to become per-soil, with the
actual numbers checked directly (not assumed): raising Sand's total OM
alone at the shared `D_scale` made its burst LATER and broader, not
sharper (`peak_frac` 0.75→0.88 as `om_total` rose 3000→6000); giving Sand
its own, much larger `D_scale` fixed that. The same Sand candidate that
crashes cleanly at n=40 by `t=1500` needed `T=3000` at n=128 to fully
register the crash (more habitat cells → more peak biomass → a
proportionally longer, but still first-order-bounded, decay time);
lengthening Vertisol's own window to match instead broke ITS slow-rise
shape (`peak_frac` 0.91→0.46). Report the result plainly: **every §3
pass/fail criterion is met, at both n=128 and n=40 (16 of 16 checks
pass)** — Sand shows a real, growth-dominated burst (`R_peak=15.79`, 39×
the maintenance floor; 77% of R at peak is the growth term, not
maintenance) that genuinely crashes (`R_end/R_peak=0.019`) as its
concentrated coating pool empties (`fuel_remaining` 1.0→0.015); Vertisol
shows a genuine, sustained, staggered rise (55 of 70 macropore regions
recruited, late-third mean 5.8× its early-third mean, `R_end/R_peak=0.98`,
no dip) with roughly half its fuel still unspent at the end. Distinctness
(computed on time-normalized curves) is 2.17, more than 7× the 0.3
threshold. Note that neither curve was hand-drawn — both emerge from the
SAME finite-fuel + Monod-limitation + implicit-diffusion mechanism used
since Stage 6, with per-soil physical inputs chosen to match each soil's
real textural contrast (Sand: coarse, low-porosity, concentrated fast
pathways; Vertisol: fine, high-porosity, distributed slow pathways).

Stage 11 inputs in `results/stage11/`:
   - `stage11_respiration_curves.png` — R(t) for both soils, own absolute
     time axis and normalized (0=start, 1=end of own window) axis.
   - `stage11_vertisol_diagnostics.png` — R(t), biomass(t)/B0,
     habitats_recruited(t), and mean_S_habitat(t) for Vertisol in one
     figure.
   - `stage11_headline.csv` — per-soil parameters and diagnostics
     including `R_peak_over_maintenance_floor`, `B_peak_over_B0`,
     `habitats_recruited_at_2pct_window` (vs. `habitats_recruited_final`
     — the direct evidence recruitment does NOT climb, see below);
     `stage11_sand_identity_check.csv` — Sand's Stage 10 reference numbers
     vs. Stage 11's, with relative differences; `stage11_distinctness.csv`.

When covering Stage 11, present it as a FOLLOW-UP to Stage 10 (not a new
exploratory direction): Sand is kept frozen — state plainly that this was
CONFIRMED via an exact numeric identity check against Stage 10's own
recorded R_peak/R_end (relative differences 0.032% and 0.79%, both well
under 1%), not merely asserted by not touching Sand's config. Report the
performance fix precisely: a per-timestep connected-component BFS
(`habitat_properties`) was being recomputed every single step despite
being loop-invariant (the macropore topology never changes during a run);
hoisting it out, plus an optional Jacobi convergence-tolerance early stop,
gave an 11.5× speedup for Sand and a 50–60× speedup for Vertisol (which
has far more macropore regions, 70 vs. 3, so paid a proportionally larger
tax under the bug) — with results unchanged within numerical noise. State
that a lower FIXED iteration cap (20–30) was tried and explicitly
REJECTED because it measurably changed results (up to 24.8% on Sand's
`R_end`), violating the "must not alter results" requirement — this is a
case where the more conservative fix (tolerance-based early stop, not a
lower fixed cap) was the correct engineering choice, checked empirically
rather than assumed. Report Vertisol's improvement: Stage 10 already gave
a plausible R(t) SHAPE but biomass only ever reached parity with its
initial total (`B_peak/B0=1.000`); raising `om_total` alone (2000→10000,
coating fraction deliberately left unchanged) crossed a real threshold
above which aggregate growth outpaces decay, and re-slowing `D_scale`
(1.5→0.7) kept the larger dose's rise gradual rather than reintroducing a
spike — result: `R_peak=3.76` (9.4× the maintenance floor),
`B_peak/B0=2.34` (genuine growth), peak at 95% of the window,
`R_end/R_peak=0.998` (no dip). Finally, report the one honest, UNRESOLVED
limitation rather than omitting it: `habitats_recruited(t)` reaches its
final count within the first 2% of the window and stays completely flat
afterward, at both the old and new Vertisol configuration, because the
model's activation threshold for "recruited" (`growth > 1e-9`) is crossed
almost everywhere almost instantly by the unconditionally-stable implicit
solver's per-step equilibration, regardless of transport speed — R(t)'s
own genuine gradual rise instead comes from already-active regions' growth
RATE slowly building up over time, not new regions joining, a distinction
worth explaining clearly since it means the specific SS3-named diagnostic
does not literally "climb" even though the underlying mechanism is sound.

Stage 12 inputs in `results/stage12/`:
   - `stage12_mean_bands.png` — the headline figure: both soils' ensemble
     mean R(t) (M=40 seeds each) with a 95% band, on the common
     normalized fraction-of-window time axis.
   - `stage12_metric_distributions.png` — box plots of `peak_frac`,
     `R_end/R_peak`, and the late/early ratio, both soils, all 40 seeds
     each (with individual seed points overlaid).
   - `stage12_within_vs_between.png` — within-Sand, within-Vertisol, and
     between-soil pairwise normalized-curve distance distributions.
   - `stage12_coating_sweep.png` — the optional causal check: Sand's own
     `R_end/R_peak` and `peak_frac` vs. its swept coating fraction.
   - `stage12_seed_metrics.csv` — every one of the 80 seeds' shape metrics
     AND structural descriptors (`n_regions`, `coating_cells`, dry-habitat
     fraction, matrix fraction); `stage12_metric_summary.csv` (mean +/- SD
     per soil per metric); `stage12_significance.csv` (U, z, p, Cohen's d,
     rank-biserial r per metric); `stage12_within_vs_between.csv`
     (bootstrap gap CI); `stage12_coating_sweep.csv`.

When covering Stage 12, present it as a STATISTICAL VALIDATION of the
project's central claim (that structure alone, not a lucky random
realization, drives each soil's respiration shape), built entirely on
Stage 11's frozen model — state explicitly that no parameter was
re-tuned. Explain the design: for each soil, M=40 independent random
`pore.seed` values generate fresh pore-field realizations of the SAME
structural class (same `om_total`, coating fraction, `D_scale`, `T`,
biology — only the random spatial layout changes), 80 simulations total
in 13.0 minutes using the Stage-11 performance fix. Report the headline
numbers precisely: every shape metric tested shows COMPLETE rank
separation between the two 40-seed ensembles (Mann-Whitney rank-biserial
`r = ±1.000`, `p` in the `1e-14`–`1e-15` range, Cohen's `d` from 5.0 to
36.6 — note that `d=0.8` is conventionally called "large," so these are
one to two orders of magnitude beyond that); the mean between-soil curve
distance (2.06) is 13× Sand's own within-soil spread and 127× Vertisol's,
with a bootstrap 95% CI on the gap of `[1.963, 1.989]`, excluding zero by
a wide margin. State plainly that all statistics (rank transform,
Mann-Whitney U with tie/continuity correction, Cohen's d, rank-biserial
correlation, bootstrap CI) were implemented directly in numpy, not scipy,
consistent with the project's invariant throughout. Report the structural
descriptors that justify the causal story: Sand consistently realizes 1–9
macropore regions (mean 2.6) in an always-18%-matrix soil; Vertisol
consistently realizes 55–86 regions (mean 71.8) in an always-68%-matrix
soil — note that `matrix_fraction` is exactly seed-invariant per soil (a
property of the shared literature PSD's marginal distribution, not the
spatial arrangement), which is itself worth explaining since it shows
precisely what does and doesn't vary with the seed. Cover the optional
bounded causal check (sweeping Sand's own coating fraction over 4 values,
5 seeds each): `R_end/R_peak` (crash depth) drops monotonically as
coating fraction rises (0.055→0.022→0.006→0.001), direct evidence the
swept FEATURE itself drives the metric, not merely "Sand always crashes."
Close by stating the conclusion the numbers support: because each soil's
structural class is held fixed while only the random layout varies, and
every seed reliably reproduces that soil's characteristic shape with
between-soil separation one to two orders of magnitude larger than
seed-to-seed noise, each soil's respiration shape is determined by its
defined structural features — not by a lucky random realization.

Stage 13 inputs in `results/stage13/`:
   - `stage13_matrix_phase_diagnostic.csv` — the root-cause number: real
     measured PSD volume fraction below `K_d_low=30um` per soil (0.0000
     for Sand and Vertisol, 0.1852 for Loess).
   - `stage13_headline_real.csv`, `stage13_vs_stage11_comparison.csv` —
     Stage 10/11's frozen config re-run on real data vs. the exact Stage
     11 literature-PSD reference numbers.
   - `stage13_search_trace_{sand,vertisol}.csv` — the budget-capped search
     trace (2 of the allowed 15 probes per soil, each showing
     `OM_total_placed=0.0` regardless of the knob swept).
   - `stage13_region_connectivity_comparison.csv` — model emergent
     region count/largest-region fraction vs. real Gamma/r* per soil, all
     three soils.
   - `stage13_psd_real_vs_literature.png`, `stage13_respiration_curves.png`.

When covering Stage 13, present it as a DATA-SOURCING + RECALIBRATION-CHECK
stage, not a new mechanism or a new modeling direction — `dual_porosity.py`'s
physics is untouched throughout. Explain the sourcing precisely: a companion
CT-scan + nnU-Net segmentation + pore-network-topology research project
(read in full before use) provided real measured pore-size distributions
and Track E connectivity metrics (Euler number, connectivity density,
connectivity probability Gamma, degree of anisotropy, tortuosity, and the
Vogel crossover radius r*) for the same three soil types this model already
uses, replacing the literature-informed lognormal-mixture PSDs used since
Stage 3.5. State plainly, as a genuine finding rather than an
implementation detail, that real Sand and Vertisol data carry EXACTLY 0%
measured pore volume below this model's `K_d_low=30um` macro/matrix split
threshold (in every one of 50 and 32 log-spaced measurement bins
respectively) — the identical CT-resolution limitation that motivated
Stage 3.5's original move away from real image-derived PSDs, now
re-encountered independently with fresh, current real data several stages
into a far more calibrated model. Explain the mechanistic consequence
precisely: because ALL organic matter in this model is placed inside
`porous_matrix_mask` (cells with diameter below `K_d_low`), a 0% sub-30um
volume fraction makes that mask empty, so `OM_total_placed=0.0` regardless
of `om_total`. Report the headline numbers: Stage 10/11's frozen
configuration, re-run unchanged except for the PSD source, collapses both
soils to `R_peak=0.400` exactly (the trivial maintenance-only floor,
`m0*B0`) and outcome `dead`, down from Sand's `R_peak=15.79` (declining/
burst-crash) and Vertisol's `R_peak=3.76` (thriving/slow-rise) under the
literature PSD. Report the budget-capped search honestly: only 2 of the
allowed 15 probes per soil were run, each confirming `OM_total_placed=0.0`
independent of which knob (`om_total`, coating fraction, `D_scale`, `T`)
was swept — stopped early not because the search failed, but because the
structural cause was directly and repeatedly confirmed, consistent with
this project's standing practice of naming a confirmed mechanism rather
than exhausting a search budget for its own sake. Cover the region-
connectivity comparison (SS2) evenhandedly: the model's emergent macropore
structure agrees QUALITATIVELY with all three real soils' own
"one dominant, highly-connected component" story (real largest-component
fraction 92.7–96.3%; model largest-region fraction 0.998–1.000), but this
agreement is not independently earned for Sand/Vertisol (a trivial
consequence of the entire grid being one macropore region once the matrix
phase vanishes) and the model's single global diameter threshold has no
mechanism to reproduce real soils' simultaneous large-scale connectivity
and small-scale fragmentation for any soil, including Loess — state this
limitation directly rather than implying closer agreement than what was
found. Close by stating the conclusion precisely: real data changed the
calibrated behavior, the cause is structural (a genuine data-resolution/
model-convention mismatch) rather than a parameter-tuning shortfall, and
fixing it is explicitly out of this stage's scope (no physics changes, no
new free parameters) — a clearly-named, honest negative result, not a
failure to find the right numbers.

Stage 14 inputs in `results/stage14/`:
   - `stage14_identity_check.csv` — proof the new debug instrumentation
     changes nothing (byte-identical `R_t`/`B_total_t`/`habitats_
     recruited_t`/`fuel_remaining_t`, debug run vs. non-debug run).
   - `stage14_mechanism_summary.csv`, `stage14_region_recruitment.csv` —
     per-region first-activation timestep and time-to-half-max.
   - `stage14_growth_distribution_heatmap_zoomed.png` — the direct visual
     confirmation: a near-vertical activation front at the window's left
     edge, then a long, slowly-brightening plateau.
   - `stage14_activation_vs_halfmax_histograms.png` — boolean activation
     (sharp spike near t=0) vs. graded half-max timing (spread across the
     window), side by side.
   - `stage14_graded_metric_vs_R.png`, `stage14_graded_metric_
     correlations.csv` — the candidate metrics plotted against R(t) and
     their Pearson correlations.

When covering Stage 14, present it as a DIAGNOSTIC stage confirming a
previously-UNCONFIRMED suspected cause, not a fix or a retuning. Restate
the Stage 11 puzzle precisely: `habitats_recruited(t)` (a boolean count of
macropore regions with any actively-growing cell) reached its final value
within the first ~2% of Vertisol's observation window and stayed exactly
flat afterward, even though `R(t)` itself kept genuinely rising — Stage 11
suspected, but never confirmed, that the implicit solver's per-step
equilibration crosses the very low activation threshold (`growth>1e-9`)
almost everywhere almost instantly, regardless of transport speed. Report
that Stage 14 confirms this DIRECTLY, not by elimination: instrumenting the
simulation to record the actual growth VALUE on every macropore cell at
every timestep (a new, off-by-default `debug_growth_distribution` flag,
verified via an exact byte-identical check to add no side effects) shows
that of Vertisol's 70 macropore regions, all 55 that ever become active do
so within the first 2% of the 1500-step window, with a median first-active
timestep of 2. State the corollary clearly: this makes the boolean count
not merely incomplete but actively MISLEADING as a rise-tracking
diagnostic — its correlation with normalized R(t) is only r=0.072
(essentially zero), whereas a graded, growth-WEIGHTED "active mass" metric
(the sum of raw growth over every macropore cell, each step) correlates at
r=0.996, and a per-region "time to half its own eventual peak growth"
metric independently shows real growth RATE building up gradually across
44% of the window on average — both directly demonstrating, at the
mechanism level, what Stage 11 could previously only infer from R(t)'s
aggregate shape. Recommend the graded metric as an ADDITION for future
stages, not a replacement: the boolean count remains correct for the
narrower thing it actually measures (near-instantaneous first activity, a
genuine, structurally real property of this solver/threshold combination,
not a bug to be fixed). Close by confirming that none of Stage 10/11/13's
calibrated numbers moved: the debug run's R(t), biomass(t), habitats_
recruited(t), and fuel_remaining(t) are byte-identical to a non-debug run
of the same configuration, and match Stage 11's own recorded reference
numbers to well under 1%.

Stage 15 inputs in `results/stage15/`:
   - `stage15_search_trace.csv` — the 5-probe (of the allowed 15) directed
     search for Loess's own knobs, including the 2 probes that directly
     reproduced Stage 7/9's "blooming into the same shape" failure before
     the adopted, smaller-dose candidate resolved it.
   - `stage15_headline_n{40,128}.csv`, `stage15_pass_fail_n{40,128}.csv`,
     `stage15_distinctness_n{40,128}.csv` — the extended three-soil
     pass/fail table and distinctness check, at both grid resolutions.
   - `stage15_three_soil_comparison_n{40,128}.png` — the headline figure.

When covering Stage 15, present it as the highest-uncertainty item in this
continuation, per its own framing, with a genuinely mixed result — do not
round it up to a clean success or down to a clean failure. State the
deliberate, logged deviation from a literal reading of the prompt plainly
and first: Sand and Vertisol are kept at Stage 10/11's frozen,
literature-PSD configuration rather than switched to real data, because
Stage 13 already established real data structurally breaks both soils'
mechanism (zero recorded pore volume below the macro/matrix split
threshold) — using real data for all three soils here would have collapsed
Sand and Vertisol to the identical trivial floor, making any three-way
comparison meaningless. Report the search precisely: only Loess searches
its own knobs (om_total, coating fraction, D_scale, T) against real
PSD/connectivity data, and 2 of the first 5 probes directly reproduced
Stage 7/9's own named failure mode — raising Loess's total/D_scale enough
for real growth pushed its curve toward Vertisol's own rising shape,
collapsing their pairwise distance to 0.05–0.22 (well under the 0.3
threshold) — resolved not by loosening the threshold but by finding a
substantially smaller, faster-clearing dose. Report the headline result
exactly: three-way distinctness succeeds robustly at BOTH n=40 (0.659) and
n=128 (0.442); Loess's "low" target (lowest R_peak of the three soils)
succeeds at n=128 but FAILS at n=40 — a genuine grid-resolution
disagreement, not smoothed over; Loess's "erratic" target was tested
across every probe and additional configuration explored (om_total
150–40000, coating 0.1–1.0, D_scale 0.3–80, T 600–3000) and never
achieved — every single curve has exactly one peak, confirmed as a
structural limitation of this model's smooth, monotonically-forcing
mechanism (no feedback loop capable of producing a genuine dip-then-
recovery), independently reconfirming Stage 3.5's own identical finding
under literature PSD, now shown to hold under real PSD data too and
across a much broader search. Mention the byproduct finding: Vertisol's
Stage-11 configuration, never tested at n=40 before (Stage 11 restricted
itself to n=128 only), turns out to also pass every one of Stage 11's own
success criteria at n=40 — a genuine, previously-unknown robustness
result surfaced only because this stage's own instructions required
checking both grids. Close by stating the verdict precisely: this is a
real, mechanistically-explained partial result — distinctness and "low"
achieved (with one honest grid caveat), "erratic" directly shown to be
outside this model's reach — not a failure to search hard enough and not
a success dressed up as complete.

Stage 16 inputs in `results/stage16/`:
   - `stage16_geodesic_distance.csv`, `stage16_geodesic_histograms.png` --
     the stage's PRIMARY structural result: the emergent distribution of
     geodesic distance from each macropore region to its nearest
     substantial OM, per soil, in grid cells and microns.
   - `stage16_geodesic_threshold_sensitivity.csv` -- how that distribution
     moves as "substantial" is made stricter; secondary, and explicitly not
     used for any pass/fail number.
   - `stage16_ladder_trace.csv` -- the 9-probe (of the allowed 15) dose and
     transport ladder, rung by rung.
   - `stage16_ss5_reconciliation.csv` -- the diagnostic showing why the
     prompt's own SS5 evidence table does not transfer to the shared rule.
   - `stage16_headline_n{40,128}.csv`, `stage16_pass_fail_n{40,128}.csv`,
     `stage16_recruitment_timing.csv`, `stage16_distinctness_n{40,128}.csv`,
     `stage16_legacy_reproducibility.csv`, `stage16_clumping_sensitivity.csv`.
   - `stage16_headline_comparison_n{40,128}.png` -- the headline figure.

When covering Stage 16, lead with what it is: a MODEL-WIDE structural
change that deletes a hand-set per-soil carbon parameter, not a
calibration. The per-soil coating/trapped OM split is retired in favor of
ONE shared, structure-blind placement rule with spatially uniform
dissolution, under the invariant that no term in the model has a rate
depending on distance to habitat; the Stage 6-15 path is preserved behind a
config flag and still reproduces Stage 11 to within 0.004%.

Report the outcome as the genuinely split result it is, without rounding it
either way. What survived removing the knob: both soils keep their Stage
10/11 target shapes at both grid resolutions (Sand 6/6 checks, Vertisol
5/5), and Sand-vs-Vertisol distinctness is essentially unchanged from the
hand-set configuration (2.082 vs 2.001 at n=128; 2.152 vs 2.217 at n=40) --
so the per-soil carbon knob was buying nothing that structure did not
already provide, which makes the project's central claim stronger, not
weaker. What failed: recruitment timing, at every rung of the ladder and
independently of dose (`frac@2%=1.000` across a 17x range of `om_total` and
a 114x range of `D_scale`). Explain WHY, because the reason is the same
finding as the structural one: under a shared rule that fills the porous
matrix, the emergent fuel-habitat distances are only 1-2 grid cells for
BOTH soils (Sand median 10 um, Vertisol 20 um), 66 of Vertisol's 70 regions
physically touch a fuel-bearing cell, and the boolean activation test needs
only `growth > 1e-9` -- so recruitment is saturated by the pore field's own
geometry before any parameter is chosen. State plainly that the two soils'
distance distributions came out far MORE similar, and far smaller, than
expected, and that they do not land in Portell et al. (2018)'s 252 +/- 111
/ 490 +/- 262 um range; note the domain-size difference (6.8 mm vs 1.28 mm)
but also that it does not explain the gap. Note that the soils do separate
strongly (Vertisol to 200-370 um) once "substantial" OM is defined more
strictly, and that this is reported as a secondary observation rather than
swapped in as the headline definition.

Cover the ladder's quantitative answer to "how much can be derived from
structure": rung 2's shared concentration per unit matrix volume accounts
for 22.2% (linear) / 46.5% (log) of the required per-soil dose difference,
and the structure-derived `D_scale` for only 2.2% / 19.6% of the transport
difference -- so no rung was adopted and a per-soil residual of ~4.5x
remains hand-set.

Report both corrections this stage made to earlier claims, since both were
verified by direct measurement rather than argued: (i) Stage 14's causal
attribution of near-instant habitat activation to the implicit solver is
REFUTED -- with transport switched off entirely regions still activate
within 2-4 timesteps and recruitment finishes sooner, and the real cause is
the coating construction fuelling all ~70 regions simultaneously; a dated
correction is appended to `STAGE14_RESULTS.md` and Stage 14's diagnostic
findings stand. (ii) The "79% selectivity" figure is the boolean recruited
count, not Portell's peak-biomass metric; on Portell's own metric Stage
11's baseline is 14.3% and the shared rule gives 21.4%, against their 9.5
+/- 4.0% -- selectivity got slightly worse, not better.

Be explicit about one tension rather than resolving it silently: the
criterion that failed is defined on the boolean recruited count, which this
same stage's instructions (following Stage 14) direct be retired as a
headline number because it is effectively uncorrelated with R(t). Measured
on the graded per-region time-to-half-max instead, recruitment under the
shared rule is MORE spread across the window than the baseline (median 53%
vs 44%). The stage counts the criterion as failed, because that is what it
says, and flags the tension.

Close by stating the verdict precisely: a clean, mechanistically-explained
negative on recruitment timing, with everything else surviving the removal
of the per-soil knob -- and, per the stage's own instructions, the signal
to attempt a physiological rather than geometric fix (dormant/active
biomass with substrate-dependent reactivation), which was explicitly out of
scope and not attempted. Do not present the clumping diagnostic as a
rescue: it restores staggering at n=128 but fails at n=40, and is reported
as a lead only.

Keep the tone scientific and precise. Do not invent results not present in
`LOGBOOK.md`,
`docs/stage_results/STAGE{1,2,3,3_5,4,5,6,7,8,9,10,11,12,13,14,15,16}_RESULTS.md`,
or `results/`.
