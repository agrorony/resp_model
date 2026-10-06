# Stage 13 RESULTS -- real measured PSD + connectivity data

`prompts/soil_respiration_prompt_stage13_real_data.md`. A DATA-SOURCING +
RECALIBRATION-CHECK stage, not a new mechanism: `dual_porosity.py`'s physics
is untouched. Replaces Stages 3.5-12's literature-parametrized PSDs
(`configs/psd_literature.yaml`) with real measured PSD + Track E
connectivity data from the companion "resarch exercise" project.

## 0. Result up front

**Real data breaks the Stage 10/11 pass/fail table completely, for a
structural reason that no amount of retuning within the four allowed knobs
(`om_total`, coating fraction, `D_scale`, `T`) can fix: Sand and Vertisol's
real measured PSDs have exactly 0% recorded pore volume below `K_d_low=30
um`, the model's macro/matrix split threshold.** Every gram of organic
matter in this model is placed inside `porous_matrix_mask` (diameter <
`K_d_low`) -- with that mask structurally empty, `OM_total_placed=0.0` for
both soils regardless of `om_total`, coating fraction, `D_scale`, or `T`,
and both runs collapse to `R(t)` sitting exactly on the trivial maintenance
floor (`R_peak=0.4000`, outcome `dead`) for the entire window. The
budget-capped search (SS3) confirms this directly and stops after 2 probes
per soil (well under the 15-iteration cap), not because the search failed,
but because the mechanism itself has nowhere to place substrate -- an
honest, structural non-fixable result, reported plainly rather than forced.
Loess (never part of Stage 10/11's frozen set) retains 18.5% of its measured
pore volume below 30 um and was not subject to this collapse. Region-
connectivity (SS2) comparison against real Track E numbers is reported for
all three soils, without forcing agreement.

## 1. Data sourced (SS0), with exact provenance

Read `PROJECT_STATUS.md` and `DATA_CATALOG.md` in the companion "resarch
exercise" project folder in full (2026-09-14). Key findings that shaped
sourcing decisions, and exactly what was used -- full detail in each
soil's `data/psd_measured/<soil>/source.txt`:

| Soil | PSD table used | Connectivity scalars used | Limitation, stated plainly |
|---|---|---|---|
| Sand (Rehovot) | `rehovot_150z` (2026-05-24, already in `data/psd/rehovot/` since Stage 3) | Current full-volume `samp_2` run (2026-08-24/25): porosity 0.3050, chi=-53,704, conn.density 57.94 mm^-3, Gamma=0.9889, DA=0.1956, tortuosity 3.666/3.398/3.148 (all converged), r*=100.8um, PSD 30-150um frac 0.7030 | The current full-volume run's own bin-level PSD table lives only on the companion project's network share (`Z:\...`), unreachable in this environment (unmounted). `rehovot_150z` is real measured data, just an earlier, resolution-truncated (30um min diameter) processing pass -- the nearest documented substitute, not an invention. |
| Loess (Mishmar HaNegev) | Native 5.85um `mishmar_hanegev_maoz_3_5p85um` full-volume run (2026-08-03), copied fresh from the companion project | Same volume, matched-resolution scalars: chi=-65,340, conn.density 326.37 mm^-3, Gamma=0.9273, DA=0.1227, tortuosity ALL 3 axes NaN (confirmed `porespy` finite-difference solver non-convergence bug, NOT "no percolating path" -- see `connectivity_validation_summary.md` Part B), r*=67.0um, PSD 30-150um frac 0.7321. Matched-~15um n=3 physical-replicate mean+/-SE also recorded (different volumes): conn.density 151.2+/-110.7, Gamma 0.9324+/-0.0018, DA 0.148+/-0.023, r*=133.3+/-42.4um | None -- this is the best, most current, full-resolution real data available for this soil. |
| Vertisol (Bnei Re'em) | "Canonical" (`bnei_reem_fresh_bnei_reem_i4`, ~Jun/Jul 2026) full-volume run, copied fresh from the companion project | Specimen B unified pipeline (2026-09-07, the closer numeric match to "canonical"): chi=+9,771, conn.density -10.4699 mm^-3, Gamma=0.859338, DA=0.096700, tortuosity axis0=NaN (genuine non-percolation)/5.9585/2.7848, r*=920.48um, PSD 30-150um frac 0.391276. n=2 specimen mean+/-SE also recorded: pore 18.40+/-3.24%, Gamma 0.807+/-0.052, DA 0.143+/-0.046 | **Specimen identity is genuinely unresolved in the source project itself** (`DATA_CATALOG.md`, 2026-09-07 entry): an EARLIER candidate second-specimen run was found to have consumed raw rotational X-ray projections instead of reconstructed slices and was retracted; the "canonical" volume used here for its PSD table is neither of the two 2026-09-07 rebuilt specimens, and its own scalars are numerically closest to (but not identical to) Specimen B's. Flagged, not silently resolved -- this is the source project's own open item, not something Stage 13 can settle. |

All three CSVs and this table's numbers are logged with exact run
names/dates in `configs/psd_measured.yaml` and each soil's `source.txt`
(the "log exactly which files/versions" requirement, SS0.2).

## 2. Matrix-phase diagnostic -- the root cause, checked directly (SS1/SS3)

Before running any simulation, the fraction of each soil's real measured PSD
volume below `K_d_low=30um` (the model's macro/matrix split threshold) was
computed directly (`stage13_matrix_phase_diagnostic.csv`):

| Soil | volume fraction below 30um | nonzero bins below 30um |
|---|---|---|
| Sand | **0.0000** | 0 |
| Loess | 0.1852 | 4 |
| Vertisol | **0.0000** | 0 |

This is not a marginal or stochastic effect -- it is exactly zero, in every
one of Sand's 50 and Vertisol's 32 log-spaced measurement bins below 30um,
deterministically for any grid size or seed. The real, CT-derived local-
thickness PSD pipeline simply cannot resolve pores below roughly this scale
at 15um voxel resolution (Sand, Vertisol) -- **the identical limitation
`psd_parametric.py`'s own docstring already names as the reason Stage 3.5
replaced Stage 3's real image-derived PSDs with literature mixtures in the
first place** ("Stage 3 fed the model each soil's real measured PSD. Those
PSDs turned out to be poor: 30 um minimum resolvable diameter..."). Stage 13
brings real data back and hits the identical wall for two of three soils --
this validates that Stage 3.5's move was well-motivated, not merely a
convenience, rather than being a new problem introduced here.

Because `porous_matrix_mask = matrix_mask & ~solid_mask` and `matrix_mask =
~macro_mask` (`dual_porosity.build_dual_grids`), a 0% sub-30um volume
fraction means `matrix_mask` (and therefore `porous_matrix_mask`) is
EMPTY for the entire grid -- confirmed directly: `macro_fraction=1.0`,
`porous_matrix_fraction=0.0` for both soils at n=128 (`stage13_headline_
real.csv`). Every gram of organic matter this model places goes through
`split_om_coating_trapped`, which only ever assigns mass inside
`porous_matrix_mask` -- with that mask empty, `OM_coating=OM_trapped=0`
everywhere, for ANY `om_total`.

## 3. Recalibration check (SS3) -- Stage 10/11's frozen config on real data

`stage13_headline_real.csv`, `stage13_vs_stage11_comparison.csv` (n=128,
`om_total`/coating/`D_scale`/`T` exactly as Stage 10/11 left them):

| Soil | `R_peak` (Stage 11, literature) | `R_peak` (Stage 13, real) | outcome (lit -> real) | all checks pass? |
|---|---|---|---|---|
| Sand | 15.790 | **0.400** | declining -> **dead** | **False** |
| Vertisol | 3.756 | **0.400** | thriving -> **dead** | **False** |

`R_peak=0.400` for both is exactly `m0*B0=0.4`, the trivial maintenance-only
floor -- confirming, numerically, that zero growth occurs at all (no
substrate was ever placed). Neither soil's shape survives.

**Budget-capped search (SS3, `MAX_ITERATIONS=15`), same per-soil knobs
Stage 10 established:**

| Soil | probe | om_total | coating | D_scale | T | `OM_total_placed` | `R_peak` |
|---|---|---|---|---|---|---|---|
| Sand | 1 | 10x (40000) | 0.97 | 80.0 | 3000 | 0.0 | 0.400 |
| Sand | 2 | 4000 | 1.0 (2x, capped) | 80.0 | 3000 | 0.0 | 0.400 |
| Vertisol | 1 | 10x (100000) | 0.30 | 0.7 | 1500 | 0.0 | 0.400 |
| Vertisol | 2 | 10000 | 0.60 (2x) | 0.7 | 1500 | 0.0 | 0.400 |

`OM_total_placed=0.0` in every single probe, independent of the knob swept
-- **because the placement mechanism itself, not the total or its split,
requires `porous_matrix_mask` cells to exist, and none do.** The search
stops after 2 of the allowed 15 probes per soil (`stage13_search_trace_
{sand,vertisol}.csv`), not because 2 probes exhausted the space, but because
the structural cause was directly confirmed and a 3rd-15th probe would
provably repeat the identical `OM_total_placed=0.0` result -- consistent
with Stage 15's own explicit instruction to name a confirmed mechanism
rather than chasing further probes for their own sake. **This is a genuine,
reportable negative result, not a search failure**: recalibration within
the allowed lever set is structurally impossible for these two soils' real
PSDs, given the model's existing macro/matrix convention.

## 4. Region-connectivity vs. Track E (SS2) -- reported, not forced

`stage13_region_connectivity_comparison.csv` (n=128, all 3 soils):

| Soil | model `n_regions` | model largest-region fraction | real Gamma | real r\* |
|---|---|---|---|---|
| Sand | 1 | 1.000 | 0.9889 | 100.8 um |
| Loess | 4 | 0.998 | 0.9273 | 67.0 um |
| Vertisol | 1 | 1.000 | 0.859338 | 920.48 um (Specimen B); Specimen A has none in 2-1000um |

**Directional agreement**: all three real soils are reported (Part A of
`connectivity_validation_summary.md`) as overwhelmingly dominated by ONE
giant connected pore cluster (largest-component fraction 92.7-96.3%)
despite tens of thousands of separate tiny fragments. The model's own
emergent macropore structure agrees with the *qualitative* "one dominant
component, very high Gamma" story for all three soils -- largest-region
fraction is 0.998-1.000 in every case here, matching the real soils' own
extremely high largest-component share.

**Honest disagreement, not smoothed over**: for Sand and Vertisol the
model's `n_regions=1` is a trivial consequence of SS1's finding
(`macro_fraction=1.0` means the entire grid is one contiguous macropore
region by construction, not a genuine emergent network result) -- the
model cannot be said to have independently reproduced the real soils'
connectivity structure here, since there is no matrix phase left to create
any structure at all. Even for Loess, where a genuine matrix phase exists
(18.5% below 30um) and the macropore field is not the whole grid, the model
still collapses to essentially one region (4 regions, 99.8% in the largest)
-- the model's macro/matrix split is a single global diameter threshold
applied to a spatially-correlated Gaussian field, with no mechanism to
produce the real soil's co-existing large-scale connectivity AND tens of
thousands of small disconnected fragments simultaneously. Per SS2's
instruction, this disagreement is reported as-is: the model's structure
generator is qualitatively simpler than real pore networks and was never
designed to reproduce multi-scale fragmentation, only a single-scale
macropore/matrix split -- no geometry was tuned to force closer agreement.

Vertisol's real Track E story (per `crossover_radius_summary.md`) is
genuinely mixed across its two real physical specimens: Specimen A has NO
`chi(r)` crossover in the full measured range (fragment-dominated at every
scale, matching the pre-2026-09-07 "canonical" characterization the Stage
13 prompt itself quotes), while Specimen B crosses at a very large
`r*=920.48um` (the largest crossover radius of any soil/specimen in the
data). The model's own emergent single-giant-region result is at least
directionally consistent with "no meaningful internal fragmentation
regime" either way, though (as above) this is a trivial consequence of
`macro_fraction=1.0`, not an independent structural finding.

## 5. Honest verdict

**Real data materially changed the calibrated behavior -- it did not
reproduce Stage 11's numbers, and no retuning within the allowed levers
fixes it.** The cause is structural, not a tuning failure: the real,
CT-resolution-limited measured PSDs for Sand and Vertisol carry zero
recorded pore volume below this model's `K_d_low=30um` macro/matrix split,
so the entire coating/trapped organic-matter mechanism this model has used
since Stage 6 has nowhere to place substrate for either soil. This is
reported plainly, exactly as Stage 3.5 plainly reported the same
30um-resolution limitation when it first replaced real image-derived PSDs
with literature mixtures for this same reason -- Stage 13 independently
re-confirms that finding with fresh, currently-canonical real data three
stages further into a much more calibrated model, rather than resolving
it. Fixing this within Stage 13's stated scope ("do not change
`dual_porosity.py`'s physics", "do not introduce new free parameters") is
not possible; a genuine fix would require either a different macro/matrix
convention (e.g. a resolution-aware or soil-relative split threshold) or
richer real PSD data reaching below 30um for Sand/Vertisol (neither of
which this stage is scoped to attempt). Loess, never part of Stage 10/11's
frozen recalibration set, is unaffected by this collapse and is where
Stage 15's later reintegration attempt has real (not literature) structural
data to build on. The region-connectivity comparison (SS2) is honestly
mixed: qualitative "one dominant, highly connected component" agreement
across all three soils, but no genuine structural correspondence for
Sand/Vertisol given SS1's collapse, and a model structure generator too
simple to reproduce real soils' simultaneous large-scale connectivity and
small-scale fragmentation for any soil, Loess included.

## 6. Deliverables

`psd_measured.py` (new loader, mirrors `psd_data.py`/`psd_parametric.py`'s
shape exactly so `build_dual_grids` needs no other change); `configs/
psd_measured.yaml` (soil -> real PSD csv + curated, provenance-logged
connectivity scalars); `configs/soil_{loess,sand,vertisol}_measured.yaml`
(`pore.mode: psd_measured`); `data/psd_measured/<soil>/{psd_table.csv,
source.txt}` (real data + full provenance narrative per soil);
`dual_porosity.py` (`build_dual_grids` now accepts `pore.mode:
psd_measured` alongside the existing `psd_parametric`, unchanged
mechanics either way). `stage13_run.py`: matrix-phase diagnostic,
frozen-config recalibration check, budget-capped search, region-
connectivity comparison, all figures/CSVs in `results/stage13/`. This
file, `MODEL_SPEC.md`, `LOGBOOK.md`, `report_prompt.md` updated.

**[Stage 13] Checkpoint reached: real data sourced with full provenance
(companion project read in full, exact files/dates logged); matrix-phase
diagnostic identifies the root cause directly (0% real pore volume below
`K_d_low=30um` for Sand/Vertisol); Stage 10/11's frozen config re-run on
real data and both soils' pass/fail tables fail completely (`R_peak=0.400`,
outcome `dead`); budget-capped search (2/15 probes/soil, stopped early on
confirmed structural cause, not search failure) confirms no combination of
the 4 allowed knobs can fix it; region-connectivity vs. Track E reported
honestly (qualitative agreement, no forced match); deliverables written.
Stopping here per soil_respiration_prompt_stage13_real_data.md SS5 order of
work.**
