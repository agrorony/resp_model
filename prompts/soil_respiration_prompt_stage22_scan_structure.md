# Stage 22 — structure from µCT scan crops (Sand and Vertisol only)

You are working in the `resp_model` repo on the GPU machine (the one that ran Stage 21: 2x RTX A6000, CuPy).
Read this whole prompt before doing anything. Work in order. Stop and report at any failed gate; do not improvise around it.

## 0. Context (read first)

- Read `docs/stage_results/STAGE21_RESULTS.md`, `model3d.py`, `packing3d.py`, `stage21_run.py`, `stage21_steps34.py`, `metrics.py`.
- Stage 21 result: synthetic 3D packings (spheres + fines for Sand, porous plates for Vertisol, 5 µm voxels, 256³) with equal
  per-voxel OM and biomass and no soil-specific rules gave two alive, distinct, seed-robust respiration curves
  (S–V distance 1.23, seed noise ≤ 0.11). The difference comes from the matrix / resolved-pore partition and the wet / drained split.
- Stage 22 changes ONE thing: the structure comes from crops of Rony's segmented µCT scans (15 µm voxels) instead of packing.
  Everything else (biology, transport, exchange, wetting rule, densities, distance metric, liveness rules) is Stage 21's.
- It also adds two things that need no new physics: (i) a mapping of model time to days from literature growth rates,
  fixed in this prompt before any run, and (ii) model CUE and qCO₂ as outputs.
- Purpose: test whether the Stage 21 result survives on real soil structure, at the scan's resolution.

## 1. Hard rules

1. Sand and Vertisol only. No Loess.
2. Measured data allowed as input: the segmented scan volumes ONLY. Do not open, read or use any respiration, biomass,
   texture, CaCO₃ or other lab file. Do not compare model output with measured respiration; Rony does that himself.
3. Do not edit `configs/biology.yaml`, `configs/*`, `model3d.py`, `packing3d.py`, `dual_porosity.py`, `model.py`,
   or any Stage ≤ 21 runner. New code goes in new files (`scan3d.py`, `timemap.py`, `stage22_run.py`).
   Parameter changes in Step 5 are in-memory overrides only.
4. Fixed, inherited, not tuned: OM_density 1.0, K_density 0.16445, B0_density 0.0125, D_scale, k_matrix_diff, k_ex 0.10,
   film_K_scale, d_cut 50 µm, habitat lower limit 10 µm (`K_d_low_um`), dt 0.05. T = 3000 steps (Stage 21 used 1500;
   the first 1500 steps are the same run, the extension only covers 15 days at the fastest time mapping).
5. Equal totals asserted every run (as `run_structure` already does): sum(OM) = OM_density·N, sum(B0) = B0_density·N,
   sum(K) = K_density·N.
6. Order: distinctness first. Robustness/sensitivity (Step 5) only if Step 3 passes.
7. No search over parameters. The only contingency allowed is the Stage 21 OM map, described in Step 3.
8. Commit predictions before the runs they predict. Commit after each step.
9. Never reuse result files from an earlier run or stage as if they were new. Every Step 3–6 number must come from a run
   made in this stage, except the Stage 21 5 µm curves, which are read from `results/stage21/runs/` and labelled as such.

## 2. Inputs (verify every path; stop if any is missing)

Root on this machine (from Rony's analysis scripts): `RB = E:\PROJECTS\Yael_Mishael\Rony\remote_computer backup`
(same as `\\HIVE3065\Yael_Mishael\Rony\remote_computer backup`).

| ID | Soil | File | Voxel (µm) | Format | Expected whole-volume pore % (gate) |
|---|---|---|---|---|---|
| VA | Vertisol, specimen A | `RB\nnUNet_resources\bnei_reem_specA_unified\inference_concatenated\bnei_reem_specA_unified.nii.gz` | 15.000149 | nnU-Net labels | 15.17 (POM 0.24) |
| VB | Vertisol, specimen B | `RB\nnUNet_resources\bnei_reem_specB_unified\inference_concatenated\bnei_reem_specB_unified.nii.gz` | 15.034357 | nnU-Net labels | 21.64 (POM 0.82) |
| SA | Sand, samp2 | `RB\10.5\rehovot_samp_2.npy` | 15.000149 | binary, 1 = pore | 30.50 |
| SB | Sand, samp3 | `RB\10.5\Rehovot_samp3_highkV_Cu0.11_15um.npy` | 15.034357 | binary, 1 = pore | 40.85 |

Vertisol labels: the training `dataset.json` says background 0, Stones 1, POM_type1 2, POM_type2 3, unused 4, Pore 5, ignore 6;
Rony's deployed analyses use pore = 5, POM = 2. The project's Methods table disagrees with both, so VERIFY:
- Report the label histogram (fraction of each value) of VA and VB.
- Gate: label-5 fraction within ±1 percentage point of the table above, and POM (2 + 3) within ±0.5 pp. If not, stop and report.
- Any voxel with value 4 or 6, or any value not in {0,1,2,3,5}: stop and report the counts.

Sand files: report dtype, unique values and pore fraction. Gate: pore fraction within ±1 pp of the table.

Inside-the-core check: Rony's tortuosity runner (`tortuosity_subvolumes_20261004/run_tortuosity_subvolumes.py`) treated these
four volumes as cropped to the core interior (`circular_container = False`). Confirm it: save the central xy, xz and yz slices
of each volume (`results/stage22/inputs_<ID>.png`) and check that no container wall or outside air is visible. If any is,
restrict crop origins to the inscribed region (reuse that script's `grid_origins` logic) and report what you did.
The nnU-Net output is y-flipped relative to the raw NIfTI; this does not matter here (geometry only). Do not flip anything.

## 3. Phase mapping (pre-registered)

Phases are the Stage 21 ones (`packing3d.PORE, MATRIX, SOLID`).

**Vertisol:** PORE = label 5; MATRIX = labels 0, 2, 3 (POM is treated as matrix: it gets the same per-voxel OM as the rest
of the matrix; real POM positions are not used); SOLID = label 1 (stones, inert).

**Sand:** PORE = 1 in the binary mask. The scan has no matrix class, but in Stage 21 Sand's OM lives in a fines matrix
(literature volume fraction `packing3d.fines_volume_fraction()`, 0.082 of the domain). Rule (decided, do not change):
the fines are the finest solid features. Compute local thickness of the SOLID phase (`packing3d.local_thickness(solid, vox)`),
and convert solid voxels to MATRIX in order of smallest solid local thickness first (ties: nearest to a pore by EDT, then
random with seed 22) until the matrix fraction of the crop equals `fines_volume_fraction()`. The rest of the solid stays inert.
This keeps the scanned visible porosity unchanged (Stage 21 added fines into the pores; here they are carved from the solid,
because the scan already images them). Report the solid-LT threshold reached in each crop.

Build a `packing3d.Structure` per crop with `phase` (uint8), `vox_um` (the specimen's voxel size from the table), `info`,
then call `packing3d.analyse` to get `lt` and the geometry report, exactly as Stage 21 does. `run_structure` then applies
Stage 21's rules unchanged (dx = vox_um / 10).

## 4. Crops and controls (pre-registered)

- **Main crops:** 100³ voxels (≈1.5 mm). List all non-overlapping 100³ grid positions inside each volume (and inside the
  inscribed region if Step 2 required it). Draw 4 per specimen at random with `numpy.random.default_rng(22)`, without
  replacement. Reject a crop only if it has zero habitat voxels; draw the next one and report every rejection.
  16 crops in total (VA, VB, SA, SB × 4). IDs: `VA_c0 … SB_c3`.
- **Domain-size check:** one 256³ crop per specimen (the grid position closest to the volume centre). If a volume is smaller
  than 256³ in any axis, use the largest cube that fits and report it.
- **Resolution control:** Stage 21's synthetic S and V rebuilt with `packing3d.build_structure(kind, seed_offset=s, n=100, vox_um=15.0)`
  for s = 0, 1, 2 (6 structures, IDs `S15_s0 … V15_s2`). This separates "real structure" from "15 µm voxels".
  Report how many voxels thick the Vertisol plates are at 15 µm.
- **Stage 21 reference:** the 5 µm curves already in `results/stage21/runs/` (S and V, seeds 0–2, OM 1.0). Read only.

## 5. Time mapping (pre-registered; post-processing only)

The biology is in model units with r_max = 1. One model time unit = 1 / r_lit days, where r_lit is a literature maximal
specific growth rate (d⁻¹). Three values, fixed now:

| Label | r_lit (d⁻¹) | Basis (Rony verifies the sources before writing) |
|---|---|---|
| low | 0.36 | upper end of realized whole-community specific growth rates in soil, 0.0037–0.015 h⁻¹ (J. Plant Nutr. Soil Sci. 149: 457–468, 1986, "Carbon assimilation and microbial activity in soil"). A maximum cannot be below this. |
| **primary** | **2.6** | substrate-induced growth response µmax ≈ 0.11 h⁻¹ ("Estimating the biomass of microbial functional groups using rates of growth-related soil respiration") |
| high | 9.2 | generation time of the actively growing fraction 1.8–2.8 h (Blagodatskaya et al., EGU 2008 / 2010); 1.8 h → 9.2 d⁻¹ |

Day t_d corresponds to step t_d · r_lit / dt. 15 days = 108 / 780 / 2760 steps. T = 3000 covers all three.
Incubation periods (as measured in the lab): P1 days 0–3, P2 3–7, P3 7–10, P4 10–15.

## 6. Steps

### Step 0 — setup and validation
1. `git pull` on main, then branch `agents/stage22-scan-structure`.
2. Check `configs/biology.yaml` is byte-identical to main (Stage 21 `check_biology`).
3. Reproduce one Stage 21 run (S, seed 0, OM 1.0, 256³, 5 µm, T 1500) and compare R(t) with
   `results/stage21/runs/S_s0_om1_T1500_Dstd_n256.npz`. Gate: max relative difference < 1e-8. Do not overwrite the stored file.
4. Load the four volumes, run all gates in §2, save the input slice figures.
5. Write `scan3d.py` (loading, phase mapping §3, crop listing/drawing §4) and `timemap.py` (§5 conversions, period means,
   qCO₂, CUE). Unit test the Sand fines rule on a small synthetic binary: matrix fraction exact, thinnest solid first.
6. Commit.

### Step 1 — geometry and predictions (no respiration)
1. Build all 16 main crops, 4 domain-size crops and 6 control structures. For each, the Stage 21 Step 1 table:
   resolved porosity, matrix fraction, solid fraction, pore LT d10/d50/d90, largest pore cluster, spans x/y/z, S_w,
   habitat fraction, wet fraction of habitat. Add per-crop mean and SD per specimen and per soil.
2. Figure: central slice of one crop per specimen and of S15_s0 / V15_s0, phases coloured, wet vs drained pores marked
   (`results/stage22/step1_slices.png`).
3. Write `docs/stage_results/STAGE22_STEP1_PREDICTIONS.md` BEFORE any respiration run:
   - expected regime per soil (drained film + exchange vs wet-pore diffusion) and why, from the geometry table;
   - expected direction relative to Stage 21 (larger / smaller S–V distance; same / different peak order) and relative to
     the 15 µm control;
   - the period-shape classifier used in Step 6, defined now: for the four period means M1..M4, "rise" if M4 is the maximum
     and M4 > 1.1·M1; "rise-then-fall" if the maximum is M2 or M3 and M4 < 0.9·max; "fall" if M1 is the maximum and
     M4 < 0.9·M1; otherwise "flat".
4. Commit.

### Step 2 — main runs (OM 1.0, T 3000)
1. Run the 16 main crops and the 6 controls (22 runs). Use Stage 21's job runner pattern (both GPUs, one result json + npz
   per run, `results/stage22/runs/`, keys that include the crop ID and voxel size).
2. Per run: Stage 21 summary (`alive_primary`, `alive_secondary`, B_peak/B0, R_peak/floor, growth share, peak_frac,
   R_end/R_peak, delivery split wet/drained × diffusion/exchange, mass balance).
3. Commit.

### Step 3 — distinctness (gate for Step 5)
Distance = Stage 21's `metrics.pairwise_distance` on R, computed on steps 0–1500 (comparable with Stage 21) and, separately,
on 0–3000.
1. Soil-level curves: mean R/(m0·B0) over the 8 crops of each soil. S–V distance between the two means.
2. Noise: all within-soil crop-pair distances (28 per soil), split into within-specimen and between-specimen pairs.
3. All 64 cross-soil crop pairs.
4. **PASS** if (a) at least 6 of 8 crops per soil are primary-alive, (b) soil-mean S–V distance > 0.3, and
   (c) max within-soil crop-pair distance < 0.5 × min cross-soil crop-pair distance (Stage 21's criterion). Report (c) also
   with the 90th percentile instead of the max, as a secondary number; the verdict uses the max.
5. **Contingency, only if (a) fails:** run the Stage 21 OM map (0.037, 0.084, 0.192, 0.439, 1.0) on crop c0 of each specimen
   (4 crops × 5 = 20 runs), select the OM level by Stage 21's rule (both primary-alive, highest distance), rerun Step 2 at that
   level and say so in the headline. No other level, no other parameter.
6. Commit.

### Step 4 — comparisons (always run; descriptive)
1. Scan vs 15 µm control: distances S_scan–S15 and V_scan–V15 (soil means), S15–V15, and the control's own seed noise.
2. 15 µm control vs Stage 21 5 µm: S15–S5 and V15–V5 (effect of resolution alone on the synthetic structures).
3. Domain size: each 256³ crop vs the mean of its specimen's 100³ crops.
4. A 2×2 summary: structure source (synthetic / scan) × resolution (5 / 15 µm, synthetic only). Say plainly which factor moves
   the curves more.
5. Figure `results/stage22/step4_curves.png`: R/(m0·B0) vs step, (a) all crops thin + soil means thick, (b) soil means vs
   15 µm control vs Stage 21 5 µm.
6. Commit.

### Step 5 — biology sensitivity (only if Step 3 PASSED)
One crop per soil: the crop whose curve is closest to its soil mean. One-at-a-time ×0.5 and ×1.5 on r_max, Ks, Y, m0, m_s,
k_dis (in-memory override of the loaded biology dict; the file stays untouched), k_ex and K_density. 8 × 2 × 2 soils = 32 runs.
Report per perturbation: both alive?, S–V distance (0–1500), whether the peak order and peak-time order of S and V are unchanged.
No perturbed value is adopted; this step only reports robustness. Commit.

### Step 6 — time mapping and physiological outputs (post-processing of Step 2 runs)
For each r_lit (low / primary / high) and each soil (crop means ± SD across the 8 crops; also per specimen):
1. Period means P1–P4 of R, normalized to R(t=0), plus the S/V ratio per period; the shape class from the Step 1 classifier.
2. Cumulative R to day 15 (model units, and as fraction of initial OM).
3. qCO₂ = mean R/B over days 10–15, converted to d⁻¹ with r_lit; also the V/S ratio (independent of the mapping).
4. CUE over days 0–15 = 1 − ΣR / Σuptake (state that R includes maintenance, so this is a whole-community value and can be low or
   negative when uptake is small).
5. Figure `results/stage22/step6_periods.png`: bars of P1–P4 per soil, one panel per r_lit, crop points overlaid.
6. Table and figure only. No comparison with measured data, and no statement about agreement with the lab.
7. Commit.

### Step 7 — write-up
1. `docs/stage_results/STAGE22_RESULTS.md` in Stage 21's format: pre-flight, Step 0 validation, input gates (label histograms,
   pore fractions), phase mapping, geometry table, predictions vs outcomes (one line each: right / wrong), Steps 2–6 tables,
   verdict, limitations, run times.
2. Limitations must include at least: the scanned cores are intact cores while the lab incubations used disturbed soil;
   15 µm voxels make every visible pore habitat (no pores between 10 and 15 µm); Sand fines rule is a model choice;
   POM treated as matrix; time mapping depends on r_lit (state how the shape classes change across the three values);
   no oxygen limitation; single soil-level OM density; 100³ crops vs the domain-size check.
3. LOGBOOK entry (headline in the same style as Stage 21's) and the run-time table.
4. Commit, push the branch, open a PR to main. Do not merge.

## 7. Budget and stop rules

- Expected compute: Step 0 reproduction ≈ 4 min; Step 2 ≈ 22 runs at 100³ (≈ 0.5–1 min each at T 3000); domain-size check 4 runs
  at 256³, T 3000 (≈ 8–11 min each); Step 5 32 runs at 100³. Stop and report if any single 100³ run exceeds 10 min.
- Stop and report (do not work around) if: a path is missing; a label or porosity gate fails; the Stage 21 reproduction differs;
  mass balance error > 1e-10 in any run; biology.yaml differs from main.
- Final message: the headline, the Step 3 verdict with the three numbers, the Step 4 2×2, the Step 6 shape classes for the
  primary r_lit, and the PR link.
