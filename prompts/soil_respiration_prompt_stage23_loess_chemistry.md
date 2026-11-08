# Stage 23 — add Loess, and give each soil its measured Corg and Cmic

You are working in the `resp_model` repo on the GPU machine that ran Stage 22 (scan structure).
Read the whole prompt before doing anything. Work in order. Stop and report at any failed gate; do not improvise around it.

## 0. Context

- Read `docs/stage_results/STAGE22_SCAN_STRUCTURE_RESULTS.md`, `stage22_scan_structure_run.py`, `scan3d.py`, `timemap.py`,
  `test_stage22_scan_structure.py`, `model3d.py`, `packing3d.py`, `metrics.py`.
- Stage 22 (scan structure): 100³ crops of Rony's segmented µCT scans (15 µm), equal OM and biomass per voxel in both soils,
  Stage 21 rules. Sand and Vertisol were alive and distinct (0.554 on steps 0–1500), robust to biology, Sand peaking first and higher.
- Stage 23 makes two changes, on purpose and with Rony's approval:
  1. **Loess is added** (three Mishmar HaNegev scans at ~15 µm).
  2. **Each soil gets its own initial OM and biomass** from Rony's lab means (Corg and Cmic), converted to per-volume values
     with the measured core bulk density. Until now the totals were equal by design.
- Because change 2 means the soils no longer differ only in structure, the stage uses a **structure × chemistry cross**:
  every soil's crops are run with every soil's chemistry. That separates the effect of structure (rows) from the effect of the
  initial OM/biomass (columns).
- Name: this is **Stage 23**. The repo has two Stage 22 documents (equal wetting, Oct; scan structure, Nov); do not rename them.

## 1. Hard rules

1. Measured data allowed as input: the segmented scan volumes, and ONLY the Corg, Cmic and core bulk-density values written in §3
   of this prompt. Do not open any lab file (respiration, biomass, organic carbon, texture, CaCO₃, bulk density). Do not compare
   model output with measured respiration. Rony makes that comparison afterwards, using the rule fixed in §7.
2. Do not edit `configs/biology.yaml`, `configs/*`, `model3d.py`, `packing3d.py`, `scan3d.py`, `timemap.py`, or any earlier
   runner. New code goes in new files (`stage23_run.py`, plus `scan3d_loess.py` / `chem23.py` if needed). Per-soil OM and B0
   are passed as arguments / in-memory copies of `hp`; the biology dict is never changed in this stage.
3. Fixed, inherited, not tuned: biology, D_scale, k_matrix_diff, k_ex 0.10, film_K_scale, d_cut 50 µm, habitat lower limit 10 µm,
   dt 0.05, T 3000, **K_density 0.16445 (common to all soils and all arms)**, crop size 100³, crop seed 22, Sand fines rule from
   Stage 22.
4. Equal-total assertions stay, with per-arm targets: sum(OM) = OM_density·N, sum(B0) = B0_density·N, sum(K) = K_density·N,
   using that run's OM_density and B0_density.
5. No parameter search. No arm or crop may be dropped or added after seeing results.
6. Commit predictions before the runs they predict. Commit after each step.
7. Never present older results as new. Stage 22 runs may be read only for the Step 0 reproduction gate.

## 2. Inputs

`RB = E:\PROJECTS\Yael_Mishael\Rony\remote_computer backup`. Paths below are the ones recorded by Rony's tortuosity run
(`tortuosity_subvolumes.csv`, column `input_path`). Stop if any is missing.

| ID | Soil | File | Voxel (µm) | Pore label | Gate: whole-volume pore % | Notes |
|---|---|---|---|---|---|---|
| VA | Vertisol | `RB\nnUNet_resources\bnei_reem_specA_unified\inference_concatenated\bnei_reem_specA_unified.nii.gz` | 15.000149 | 5 | 15.17 | as Stage 22 |
| VB | Vertisol | `RB\nnUNet_resources\bnei_reem_specB_unified\inference_concatenated\bnei_reem_specB_unified.nii.gz` | 15.034357 | 5 | 21.64 | as Stage 22 |
| SA | Sand | `RB\10.5\rehovot_samp_2.npy` | 15.000149 | 1 (binary) | 30.50 | as Stage 22 |
| SB | Sand | `RB\10.5\Rehovot_samp3_highkV_Cu0.11_15um.npy` | 15.034357 | 1 (binary) | 40.85 | as Stage 22 |
| LA | Loess | `RB\nnUNet4SoilXrayCT\05_evaluation\psd\pom_analysis_20260829_roi_expansion\pipeline_work\mishmar_native_5p85um\mishmar_native_5p85um_label_downsample.nii.gz` | 15.007595 | 5 | 26.57 | 5.85 µm scan, majority-vote label downsample, 474³ |
| LB | Loess | `RB\nnUNet4SoilXrayCT\05_evaluation\psd\pom_analysis_20260829_roi_expansion\pipeline_work\mishmar_second_8p8um\mishmar_second_8p8um_label_downsample.nii.gz` | 14.998099 | 5 | 20.34 | 8.8 µm scan, label downsample, 947³, **circular container** |
| LC | Loess | `RB\nnUNet_resources\mishmar_hanegev_Cu011_samp_2_Rec_nlm\inference_output_concat_loess_i2\mishmar_hanegev_Cu011_samp_2_Rec_nlm.nii.gz` | 15.000149 | 5 | 26.18 (catalog 26.4) | native 15 µm; POM channel invalid (irrelevant: POM → matrix) |

Loess labels (`i2_loess` `dataset.json`): background 0, Stones 1, POM_type1 2, POM_type2 3, unused 4, Pore 5, ignore 6.
Gates for LA–LC: label histogram reported; whole-volume label-5 fraction within ±1 pp of the table; any voxel with 4, 6 or a
value outside {0,1,2,3,5}: stop and report. For LB, report both whole-volume and inscribed-disk pore fractions.
**Amendment (2026-11-08):** the LB porosity gate applies to the whole-volume fraction. The expected 20.34% was computed over
the whole 947³ volume, including the region outside the core, so the higher inside-disk fraction is expected and is not gated.
Expected POM (2+3): LA ≈ 2.0 %, LB ≈ 7.9 %, LC ≈ 0.1 %.

**LB has a circular container:** voxels outside the core must never enter a crop (outside is probably label 0 and would become
matrix). Restrict crop origins to cubes whose xy corners lie inside the inscribed disk, reusing the `grid_origins` circular logic
of Rony's `run_tortuosity_subvolumes.py` (find it under `RB`). Save central slices for every input with the disk drawn
(`results/stage23/inputs_<ID>.png`) and check LA and LC for any container or outside air too.

**Phase mapping:** Vertisol and Loess: PORE = 5; MATRIX = 0, 2, 3 (POM treated as matrix); SOLID = 1. Sand: Stage 22 rule
unchanged (`scan3d.crop_to_structure`, fines carved from the thinnest solid up to `fines_volume_fraction()`).

**Crops:** 4 per specimen, 100³, drawn with `scan3d.draw_habitable_crops` logic and `default_rng(22)` per specimen.
For VA, VB, SA, SB this must reproduce Stage 22's crops exactly (check origins against Stage 22 `step1.json`). Total 28 crops:
S 8, V 8, L 12. Domain-size check: one 256³ crop per Loess specimen (centre for LA/LC; for LB the centred cube inside the disk).

## 3. Chemistry (fixed here; do not recompute or look up)

Values are Rony's lab means. Corg = 0.58 × loss-on-ignition organic matter (as in his write-up). Per-volume values use the core
(intact cylinder) bulk density, because the model domain is an intact core.

| Soil | Cmic (mg C kg⁻¹) | LOI SOM (g g⁻¹) | Corg (g C kg⁻¹) | Core BD (g cm⁻³) | Corg per volume (mg C cm⁻³) | Cmic per volume (µg C cm⁻³) |
|---|---|---|---|---|---|---|
| Sand (Rehovot) | 199.44 | 0.016732 | 9.705 | 1.3947 | 13.535 | 278.16 |
| Vertisol (Bnei Re'em) | 406.61 | 0.048627 | 28.204 | 1.3482 | 38.023 | 548.17 |
| Loess (Mishmar HaNegev) | 315.57 | 0.027252 | 15.806 | 1.1344 | 17.931 | 357.99 |

Model inputs, scaled so that the mean over the three soils equals the Stage 22 operating point (OM 1.0, B0 0.0125):

| Chemistry | OM_density | B0_density |
|---|---|---|
| chem_S | 0.5843 | 0.008807 |
| chem_V | 1.6415 | 0.017357 |
| chem_L | 0.7741 | 0.011335 |
| chem_eq (Stage 22) | 1.0 | 0.0125 |

OM stays uniform per matrix voxel within a crop (Stage 21 rule); B0 stays proportional to K. K_density is common.
Record in the results: Corg is total organic C, so this assumes the same available fraction in every soil; LOI may overestimate
SOM in the clay-rich Vertisol; the incubated soil was disturbed while the BD and scans are of intact cores.

## 4. Time mapping and shape classifier (fixed here)

Same three r_lit as Stage 22: low 0.36, **primary 2.6**, high 9.2 d⁻¹ (sources in the Stage 22 prompt). Periods: P1 days 0–3,
P2 3–7, P3 7–10, P4 10–15.

Shape classifier for period means M1..M4, applied in this order (this replaces Stage 22's, which labelled a strong rise that
levels off as "flat"):
1. "rise-then-fall" if the maximum is M2 or M3 and M4 < 0.9·max;
2. "rise" if the maximum is M3 or M4, M4 ≥ 0.9·max and max > 1.1·M1;
3. "fall" if M1 is the maximum and M4 < 0.9·M1;
4. otherwise "flat".
Also report Stage 22's classifier for continuity.

## 5. Arms

| Arm | Structures | Chemistry | Runs |
|---|---|---|---|
| A — structure only | all 28 crops | chem_eq | 28 |
| B — own chemistry | all 28 crops | each soil's own (S crops chem_S, V crops chem_V, L crops chem_L) | 28 |
| C — swapped chemistry | all 28 crops | the two other soils' chemistries | 56 |
| D — Sand fines-rule check | 8 Sand crops | chem_eq and chem_S | 16 |
| E — Loess domain size | 3 Loess 256³ crops | chem_eq | 3 |

B + C together form the full 3 × 3 cross (structure × chemistry). Arm D uses Stage 21's original fines rule (fill the narrowest
**pore** voxels first, `packing3d._fill_fines` logic, same volume fraction) instead of carving from the solid; it tests a model
choice and its result is not adopted.

All runs: T 3000, 100³ unless stated, Stage 22 job runner pattern, both GPUs, one json + npz per run in `results/stage23/runs/`,
keys include crop ID, chemistry and fines rule.

## 6. Steps

### Step 0 — setup and validation
1. `git pull` on main; branch `agents/stage23-loess-chemistry`. Check `configs/biology.yaml` is byte-identical to main.
2. Reproduce one Stage 22 run (`SA_c0`, chem_eq, T 3000) and compare R(t) with the stored Stage 22 npz. Gate: max relative
   difference < 1e-8.
3. Load all seven inputs; run the §2 gates; save input figures.
4. Write the chemistry table of §3 into code as constants and assert the two derived columns from the raw ones (Corg per volume,
   Cmic per volume, mean-scaled OM_density and B0_density) to 4 significant figures.
5. Unit test: a run with per-soil OM/B0 passes the total assertions; arm D's fines rule fills pores, not solid.
6. Commit.

### Step 1 — geometry and predictions (no respiration)
1. Stage 22's Step 1 geometry table for all 28 crops (+ arm D Sand variants, + 3 Loess 256³): porosity, matrix, solid, pore LT
   d10/d50/d90, largest pore cluster, spans, S_w, habitat fraction, wet fraction of habitat; means and SDs per specimen and soil.
2. Slice figure with one crop per specimen (`results/stage23/step1_slices.png`).
3. Write `docs/stage_results/STAGE23_PREDICTIONS.md` BEFORE any respiration run:
   - Loess regime (wet-pore diffusion vs drained film + exchange) and where its curve should fall relative to Sand and Vertisol
     in arm A, with reasons from the geometry;
   - for the 3 × 3 cross: which factor (structure or chemistry) should dominate peak height, peak time and cumulative R;
   - expected effect of the arm D fines rule on Sand.
4. Commit.

### Step 2 — runs
Arms A, B, C, D, E (131 runs). Per run the Stage 22 summary (liveness, B_peak/B0, R_peak/floor, growth share, peak_frac,
R_end/R_peak, delivery split, mass balance). Commit.

### Step 3 — arm A: structure only, three soils
Distance = Stage 22's (`metrics.pairwise_distance` on R/(m0·B0)), on steps 0–1500 and 0–3000.
1. Soil-mean curves (mean over crops). Pairwise distances S–V, S–L, V–L.
2. For each pair: max within-soil crop-pair distance (both soils) vs min cross-soil crop-pair distance.
3. A pair is **distinct** if: ≥ 75 % of each soil's crops are primary-alive, soil-mean distance > 0.3, and max within < 0.5 × min
   cross (steps 0–1500). Report the 0–3000 numbers alongside.
4. Check that the S–V numbers equal Stage 22's (same crops, same code).
5. Commit.

### Step 4 — arms B and C: the structure × chemistry cross
1. Arm B: the same distinctness analysis as Step 3, with own chemistry.
2. Cross (B + C): for each scalar below, a 3 × 3 table of soil means (rows = structure, columns = chemistry) and a two-way
   decomposition using crops as replicates (sum of squares for structure, chemistry, interaction, residual; share of the total):
   R_peak/floor, peak time (fraction of T), cumulative R to day 15 at the primary r_lit (per volume and per mass, see §7),
   R_end/R_peak, and qCO₂ days 10–15.
3. Curve-level: mean distance between soil-mean curves that share structure but differ in chemistry, versus share chemistry
   but differ in structure.
4. One plain sentence: does structure or chemistry move the curves more, and for which outputs.
5. Figure `results/stage23/step4_cross.png`: 3 × 3 grid of soil-mean curves (crops thin), same y-axis.
6. Commit.

### Step 5 — side checks
1. Arm D: Sand stage-22 rule vs Stage 21 pore-fill rule, under chem_eq and chem_S: distances, peak height and time, delivery split.
   Does the S–V and S–L ordering of arm A / arm B change if the pore-fill Sand is used? (Compute with the arm D curves in place of
   Sand; report, do not adopt.)
2. Arm E: each Loess 256³ curve vs its specimen's 100³ mean.
3. Commit.

### Step 6 — time mapping and outputs for comparison (arms A and B)
For each r_lit and each soil (crop means ± SD; also per specimen):
1. P1–P4 period means of R, normalized to R(t = 0), and absolute (model units per voxel).
2. Shape class (§4 classifier, and Stage 22's).
3. Cumulative R to day 15 per volume, and per mass = per volume / core BD of that soil (§3).
4. qCO₂ over days 10–15 = mean R/B, converted to d⁻¹ with r_lit.
5. CUE over days 0–15 = 1 − ΣR / Σuptake.
6. Rank tables for arm B at each r_lit: soils ordered by cumulative R per mass, by P4 rate per mass, and by qCO₂.
7. Figure `results/stage23/step6_periods.png`: P1–P4 bars per soil, panels = r_lit × arm (A, B).
8. Tables and figures only; no statement about agreement with laboratory data. Commit.

### Step 7 — write-up
1. `docs/stage_results/STAGE23_RESULTS.md` in Stage 22's format: validation, input gates, chemistry table, geometry, predictions
   vs outcomes, Steps 3–6, verdict, limitations, run times.
2. Limitations at least: chemistry from lab means (no uncertainty propagated); Corg as a proxy for available OM; LOI caveat;
   disturbed incubation soil vs intact cores; 15 µm voxels (all visible pores are habitat); LA and LB are majority-vote
   downsamples of finer scans while LC is a native 15 µm scan; LC's POM channel invalid (POM → matrix, so only pore validity
   matters); Sand fines rule; time mapping depends on r_lit; no oxygen; K common to all soils.
3. LOGBOOK entry with a one-line headline, and the run-time table.
4. Commit, push, open a PR to main. Do not merge.

## 7. Comparison rule for Rony (fixed now; not computed by the agent)

After the PR, Rony compares arm B with his incubation at the primary r_lit, and reports the low/high mappings as robustness:
(1) rank order of the three soils by cumulative CO₂-C per kg at day 15; (2) shape class of each soil over P1–P4;
(3) rank order of qCO₂ over days 11–15. Agreement on each item means an exact match. Arm A is reported alongside, so that any
agreement can be attributed to structure, to chemistry, or to both, using the Step 4 cross.

## 8. Budget and stop rules

- Expected: 131 runs; 100³ runs ≈ 0.5–1 min each, 256³ ≈ 3–5 min each; about 1 hour wall time on two GPUs.
  Stop and report if a single 100³ run exceeds 10 min.
- Stop and report (do not work around) if: a path is missing; a label or porosity gate fails; the Stage 22 crop origins or the
  Stage 22 reproduction do not match; mass balance > 1e-10; biology.yaml differs from main; an LB crop touches the outside of
  the core.
- Final message: headline; arm A pairwise distinctness (three pairs); the Step 4 structure-vs-chemistry shares for peak height,
  peak time and cumulative R; arm B shape classes and rank tables at the primary r_lit; arm D effect; PR link.
