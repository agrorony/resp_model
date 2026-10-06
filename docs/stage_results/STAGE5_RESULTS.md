# Stage 5 (EXPLORATORY) RESULTS -- matric-potential wetting + inscribed-circle macropores (Variant C)

`prompts/soil_respiration_prompt_stage5_matric_inscribed.md`. `MAX_ITERATIONS
= 12` not exercised: a defined run (single matric condition + two small
robustness checks), not a search. Builds on the latest Variant C
(`dual_porosity.py`, as evolved for `stage4_run.py`: macropore/porous-matrix/
solid phases, per-soil OM fractions, `k_leak`, matrix internal diffusion,
hydro wicking, habitat-cluster K/D/S bonuses, `S_cap`) and changes only the
two things below.

## 0. Why this stage

Stage 4 showed Variant C's blocker was upstream of all its transport
machinery: the water mask was still each soil's `global_theta` quantile of
the FULL PSD, so Loess's and Vertisol's habitat started bone dry
(`wet_habitat_fraction ~ 0`) regardless of leak/wicking/OM mechanics. This
stage replaces that wetting rule with a shared matric potential and gives
macropores a real geometric size, to test whether that alone can get live
wet habitat into those two soils.

## 1. Change 1 + Change 2, as built

- **`dual_porosity.macropore_region_inscribed_diameter`**: labels connected
  macropore regions (reusing `_component_labels`) and sizes each by
  `2 * max(erosion depth in cells) * voxel_um`, via a new numpy-only
  iterative 4-neighbor erosion distance transform (`erosion_distance_4n`) --
  the exact numpy-only method the prompt specifies. `voxel_um = 10`
  (`configs/mapping.yaml`), so the grid now has a real physical scale; the
  old flat `stage4_variantC_macro_diameter_um = 80` constant is retired and
  replaced everywhere K(d)/D(d) are evaluated on macropore cells.
- **`build_dual_grids`**'s new `saturation.mode: matric_regions` path: a pore
  is wet if its SIZE < `matric_d_cut_um` (50, shared by all three soils) --
  matrix cells use their own PSD diameter `d` (always < 30um < 50, so always
  wet), macropore cells use their REGION's inscribed diameter (not the
  per-cell PSD label), solid cells are never wet. `theta` then EMERGES as
  `water_mask.mean()` rather than being imposed. The old `global_theta` path
  (`model.water_mask_and_theta`) is left completely intact for anything that
  still passes a different `saturation.mode` (Stage 4's `stage4_run.py` reruns
  unchanged, just with the retired 80um constant now replaced by geometry).
- **`simulate_dual_porosity(..., wicking_enabled=True)`**: a new toggle.
  `False` freezes `wet_mask` at exactly the static Change-1 map for the whole
  run (no dynamic rewetting); `True` is the pre-existing behavior. `k_leak`
  (the separate, always-on matrix->macropore mobile-immobile term) is
  unaffected either way.

## 2. Region geometry, per soil (n=128, `results/stage5/stage5_region_diameters.csv`)

| soil | n regions | wet (<50um) | dry (>=50um) | mean diam (um) | median | max |
|---|---|---|---|---|---|---|
| Loess    | 76 | 40 | 36 | 54.2  | 40.0 | 140.0 |
| Sand     |  3 |  1 |  2 | 140.0 | 60.0 | 340.0 |
| Vertisol | 70 | 39 | 31 | 53.7  | 40.0 | 160.0 |

**Sand's macropore network is essentially ONE region**: 13,347 of its 16,384
cells (81%) belong to a single connected macropore spanning almost the whole
grid, with an inscribed diameter of 340um -- far above the 50um cutoff, so
this entire region is dry. Its only other structure is a tiny 8-cell wet
pocket (20um) and a 24-cell dry one (60um). Loess and Vertisol, by contrast,
have dozens of small-to-medium regions straddling the cutoff almost evenly
(40 wet / 36 dry, 39 wet / 31 dry) -- their macropore geometry is genuinely
patchy at this grid resolution, Sand's is not.

## 3. Wet-habitat fraction: the headline mechanism result

| soil | wicking | emergent theta | wet_habitat_fraction | habitat_wet_share |
|---|---|---|---|---|
| Loess    | ON  | 0.705 | 0.190 | 0.595 |
| Loess    | OFF | 0.705 | 0.024 | 0.074 |
| Sand     | ON  | 0.184 | 0.167 | 0.205 |
| Sand     | OFF | 0.184 | **0.000** | **0.000** |
| Vertisol | ON  | 0.706 | 0.182 | 0.570 |
| Vertisol | OFF | 0.706 | 0.025 | 0.077 |

("emergent theta" is the STATIC, build-time value from Change 1's wet/dry
rule -- unaffected by the wicking toggle by construction, since wicking is a
runtime mechanism; `wet_habitat_fraction`/`habitat_wet_share` are measured on
the run's FINAL water mask, which wicking can expand, so they are where
wicking's effect actually shows up.)

**Under the static matric map alone (wicking OFF), Loess and Vertisol do get
genuinely live wet habitat where Stage 4 had exactly zero** (`habitat_wet_
share` 0.074 and 0.077 vs Stage 4's 0.000) -- Change 1+2 work as intended for
these two soils: their many small (<50um) macropore regions, fed by an
always-wet surrounding matrix, hold water. **Sand's habitat is COMPLETELY dry
under the static map** (`habitat_wet_share = 0.000`) -- exactly the "large
drained macropores leave Sand too dry" outcome the prompt anticipated, and
not a bug: Sand's one giant 340um region is (correctly) judged too large to
hold water at a 50um air-entry diameter, and that region contains virtually
all of Sand's habitat.

**Wicking, kept ON by default from the inherited baseline, substantially
inflates wet-habitat beyond the static map for every soil** -- Loess
0.074->0.595, Vertisol 0.077->0.570, Sand 0->0.205. Since the surrounding
matrix is essentially always wet everywhere (matrix diameter is always <30um
< 50), almost every macropore cell sits within `wicking_radius=2` of some wet
matrix cell, so wicking's rewetting rule ends up reinstating most of the
"dry" macropore volume Change 1 was designed to exclude. **This is exactly
the conflict the prompt asked to isolate**: wicking, as currently built,
largely overrides Change 1's geometric wet/dry distinction rather than
supplementing it.

## 4. Does carbon reach the microbes? Substrate reaching habitat and outcome

`mean_S_habitat_final` stays tiny in every cell of the table (0.0002-0.0014,
`results/stage5/stage5_headline.csv`) -- roughly two orders of magnitude
below the biology's half-saturation constant `Ks=0.2`. **All three soils come
out `declining` and `falling` (never thriving, never matching Loess-erratic
or Vertisol-rising), under BOTH wicking settings and across every
`matric_d_cut_um` tested (40/50/60).** This is not caused by Change 1/2
failing to wet the habitat -- Loess and Vertisol clearly do get more wet
habitat than Stage 4 gave them. It is caused by a bottleneck upstream of
wetting that Stage 5 was told to keep untouched: the inherited `S_cap`
mechanism (`stage4_variantC_S_capacity_base=0.05`, scaled by habitat-cluster
size) caps how much substrate a macropore cell's growth kinetics can ever
"see" to roughly 0.05-0.09, and even that tiny ceiling is never approached --
substrate simply never accumulates fast enough in habitat cells, wet or not,
for `f_S = S/(S+Ks)` to rise above ~0.01. Biomass then just decays via
maintenance and starvation mortality in every soil, identically, regardless
of which soil's habitat is wetter. **Distinctness is accordingly degenerate**
(`min_distance` 0.0013 wicking-on, 0.0001 wicking-off,
`results/stage5/stage5_distinctness.csv`) -- confirmed NOT an artifact of
this stage's changes: the identical `declining`/`R_peak` pinned at
`m0*B0_total=0.4` pattern was already present running the CURRENT
`dual_porosity.py` under Stage 4's OLD `global_theta` saturation mode too
(checked directly before building `stage5_run.py`), so it predates and is
independent of Change 1/2.

## 5. Robustness

- **`matric_d_cut_um` sensitivity (40/50/60, wicking ON, n=128,
  `stage5_matric_sensitivity.csv`)**: Loess and Vertisol's `wet_habitat_
  fraction` ticks up slightly from 40->50 (0.184->0.190, 0.177->0.182) then
  is flat 50->60 (no regions happen to sit in that particular window); Sand
  is completely flat at 0.167 across all three cutoffs (its only borderline
  region sits at exactly 60um, which fails the strict `<` wet test at every
  tested cutoff). No soil's outcome or shape changes across the sweep --
  robust, but only because the region-diameter distribution happens to be
  sparse right where the sweep looks; a denser sensitivity grid could still
  find sharper transitions.
- **Grid check (n=40 vs n=128, wicking ON, `stage5_grid_check.csv`)**: same
  qualitative story holds -- Loess `theta` 0.733 (vs 0.705), `wet_habitat`
  0.146 (vs 0.190); Vertisol `theta` 0.712 (vs 0.706), `wet_habitat` 0.186
  (vs 0.182); Sand `theta` 0.183 (vs 0.184), `wet_habitat` 0.139 (vs 0.167).
  All three stay `declining`/`falling` at n=40 too -- reasonably grid-robust,
  though absolute wet-habitat numbers shift by up to ~25% (Loess).

## 6. Honest verdict (SS "Honest framing")

**Change 1 + Change 2 do what they were built to do: they give Loess and
Vertisol genuinely live wet habitat (`habitat_wet_share` 0.07-0.08 even with
wicking off) where Stage 4's shared-theta quantile gave exactly zero, and
they correctly identify Sand's single, huge, well-connected macropore as too
large to hold water at a realistic air-entry diameter -- a real, informative
structural difference between the soils' geometries, not an artifact.**
However, this alone does **not** get the three R(t) shapes to match the
experiment: every soil still decays to the same trivial `declining`/`falling`
curve, because of an inherited, untouched bottleneck (`S_cap`) that caps
usable substrate in macropore cells far below what would be needed for
growth to ever exceed maintenance, regardless of how wet the habitat is.
Stage 5's fix operates entirely upstream of that bottleneck and so cannot, by
itself, resolve it -- consistent with the instruction not to tune toward the
target. Separately, wicking (part of the "keep everything else" baseline)
turns out to substantially undermine Change 1's own wet/dry distinction by
rewetting most of what the static matric map marks dry, since the always-wet
matrix sits within its rewetting radius almost everywhere; the wicking-off
run is the fairer test of Change 1/2 in isolation, and it is the one where
Sand's habitat comes out completely dry (`habitat_wet_share=0.000`) and
Loess/Vertisol's comes out modestly wet (~0.07-0.08).

## 7. Key assumptions (as requested)

- `voxel_um=10` sets the grid's physical scale; a different voxel size
  changes every region's inscribed diameter in direct proportion and was not
  swept.
- `matric_d_cut_um=50` is a documented choice, not calibrated; the
  sensitivity check (Sec. 5) found it not to matter much in the 40-60 range
  tested, but that is itself a function of where this particular grid
  realization's region diameters happen to cluster, not a general claim.
- The literature PSDs remain representative of soil TYPES, not measurements
  of the specific experimental samples (carried over from Stage 3.5).
- The erosion-based inscribed-diameter measure is an L1/Manhattan
  approximation (per the prompt's own specification), not a true Euclidean
  inscribed circle; it systematically reports a "diamond," not a disk, which
  is why the exact numeric cutoff crossings in Sec. 5 should be read
  qualitatively, not as precise geometry.

## 8. Deliverables

`dual_porosity.py`: `erosion_distance_4n`, `macropore_region_inscribed_
diameter`, the `matric_regions` saturation mode in `build_dual_grids`, the
`wicking_enabled` toggle in `simulate_dual_porosity`; `configs/mapping.yaml`:
`voxel_um`, `matric_d_cut_um` (macropore diameter constant retired);
`stage5_run.py` (headline + matric-cutoff sensitivity + wicking on/off + grid
check); all figures/CSVs in `results/stage5/`.

**[Stage 5] Checkpoint reached: region labeling + inscribed-diameter geometry,
matric-potential wetting, wicking on/off isolation, matric-cutoff and grid
robustness, all diagnostics and the comparison figure complete; deliverables
written (this file, MODEL_SPEC.md, LOGBOOK.md, report_prompt.md, README.md).
Stopping here per soil_respiration_prompt_stage5_matric_inscribed.md order-of-
work -- MAX_ITERATIONS=12 was not needed.**
