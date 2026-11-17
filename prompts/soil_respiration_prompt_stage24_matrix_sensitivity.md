# Stage 24 — how much do the results depend on the matrix / inert-solid split?

You are working in the `resp_model` repo on the GPU machine that ran Stages 22–23. Read the whole prompt first. Work in order.
Stop and report at any failed gate; do not improvise around it.

## 0. Context

- Read `docs/stage_results/STAGE23_RESULTS.md`, `stage23_run.py`, `scan3d.py`, `timemap.py`, `model3d.py`, `packing3d.py`.
- In Stages 22–23 the segmentations separate only visible pore from everything else. The model then sets the non-pore part by
  soil type:
  - **Vertisol / Loess:** all non-pore voxels are porous MATRIX (holds OM), except label 1 (about 1–2 %) as inert SOLID.
  - **Sand:** the binary Otsu mask has no matrix class, so most of the solid is inert. Only the thinnest solid is MATRIX, up to
    `fines_volume_fraction()` = 0.082 of the volume.
- So the matrix fraction differs about tenfold between soils (Sand 0.082, Vertisol about 0.82, Loess about 0.75), and it was never
  measured or varied. Arm D of Stage 23 moved Sand's matrix but kept its amount.
- Stage 24 answers two questions:
  - **A (sensitivity):** do the Stage 23 conclusions survive other matrix fractions?
  - **B (grey levels):** is there evidence in the raw scans to measure the split instead (for example multi-Otsu)?
- No parameter is adopted from this stage. It reports.

## 1. Hard rules

1. Inputs: the same seven segmented volumes, crops and chemistry constants as Stage 23 (reuse its code and constants unchanged).
   For Part B, also the raw (pre-segmentation) grey-level volumes of the same scans. No lab files.
2. Do not edit `configs/*`, `model3d.py`, `packing3d.py`, `scan3d.py`, `timemap.py` or earlier runners. New code in `stage24_run.py`
   (+ `matrix_split.py` if useful). Biology unchanged and byte-identical to main.
3. Fixed from Stage 23: crops (same 28 origins, 100³), T 3000, K_density 0.16445, chemistry constants, Stage 23 shape classifier,
   r_lit 0.36 / 2.6 / 9.2. Totals stay equal per run as in Stage 23. When the matrix fraction changes, the total OM is unchanged,
   so OM per matrix voxel changes. This is intended.
4. All variants below are fixed now. No variant may be added, dropped or tuned after seeing results.
5. Commit predictions before runs, and commit after each step. Branch `agents/stage24-matrix-sensitivity`, results in
   `results/stage24/`, push and open a PR at the end (do not merge).
6. Never present Stage 23 outputs as new. Stage 23 runs are read only for the reproduction gate.

## 2. Part A — matrix-fraction sensitivity

### Variants (pre-registered)

| ID | Sand matrix fraction (of the crop) | Vertisol and Loess: share of non-pore converted to inert | Purpose |
|---|---|---|---|
| base | 0.082 | 0 | Stage 23 reproduction |
| S04 | 0.04 | 0 | less Sand matrix |
| S16 | 0.16 | 0 | more Sand matrix |
| S30 | 0.30 | 0 | much more Sand matrix |
| VL25 | 0.082 | 0.25 | some inert grains in the fine soils |
| VL50 | 0.082 | 0.50 | many inert grains in the fine soils |
| CONV | 0.30 | 0.50 | matrix fractions brought closest together (Sand 0.30, Vertisol about 0.41, Loess about 0.37) |

**How the split is made (fixed rules):**
- **Sand:** the Stage 22 rule unchanged (`scan3d` thinnest-solid-first by solid local thickness; ties: nearest to a pore, then
  random seed 22), with the target fraction from the table instead of 0.082. If the solid phase of a crop is smaller than the
  target, use all of it and report.
- **Vertisol / Loess:** convert matrix voxels (labels 0, 2, 3 after Stage 23 mapping) to inert SOLID, largest local thickness of
  the non-pore phase first (the interiors of the thickest solid regions, farthest from pores). Ties: largest distance to a pore,
  then random seed 24. Convert until the stated share of the original matrix voxels is inert. Label 1 stays inert as before.
- Pore voxels never change.

### Runs

Arms as in Stage 23:
- **A:** equal chemistry (chem_eq).
- **B:** own chemistry (chem_S / chem_V / chem_L).

| Variant | Crops run | Arms | Runs |
|---|---|---|---|
| base | one crop per specimen (SA_c0, SB_c0, VA_c0, VB_c0, LA_c0, LB_c0, LC_c0) | A, B | 14 (gate) |
| S04, S16, S30 | 8 Sand crops (other soils unchanged, reuse Stage 23 curves) | A, B | 48 |
| VL25, VL50 | 8 Vertisol + 12 Loess crops (Sand unchanged, reuse Stage 23 curves) | A, B | 80 |
| CONV | all 28 crops | A, B | 56 |

About 198 runs at 100³, T 3000 (about 0.5–1 min each on two GPUs, so 1.5–2 h). Stop if any single run exceeds 10 min.

**Reproduction gate:** the 14 base runs must match the stored Stage 23 R(t) curves (max relative difference < 1e-8). If not,
stop. For every variant, unchanged soils use the Stage 23 runs; label this in every table.

### Geometry first (Step A1, no respiration)

For every variant and crop: matrix / solid / pore fractions, matrix–pore interface area per volume (count of matrix–pore voxel
faces / N), mean distance from matrix to the nearest pore, OM per matrix voxel in arms A and B. Write
`docs/stage_results/STAGE24_PREDICTIONS.md` before any run:
- the expected direction for Sand peak height and time in S04 → S30;
- the expected direction for Vertisol and Loess in VL25 / VL50;
- whether the soils stay distinct in CONV.

### Outputs (Step A2)

Per variant, the same summaries as Stage 23:
1. **Arm A:**
   - pairwise soil-mean distances (S–V, S–L, V–L; steps 0–1500 and 0–3000) and the Stage 23 distinctness rule per pair;
   - the order of soils by peak height and by peak time.
2. **Arm B,** at each r_lit:
   - P1–P4 period means and the shape class per soil;
   - the rank of soils by cumulative R per mass to day 15;
   - qCO₂ days 10–15 and its rank.
3. **One summary table, variants × outputs,** marking each output as same / different from base:
   - Sand–Vertisol distinct?
   - order of peak height;
   - order of peak time;
   - arm B shape classes at 2.6;
   - arm B cumulative rank at 2.6.
4. **Delivery split:** diffusion vs exchange into wet and drained habitat, per soil and variant. Did Sand's reliance on drained-pore
   exchange change?
5. **Figures:**
   - `results/stage24/sens_curves.png`: arm A and arm B soil-mean curves (primary r_lit, days 0–15), one panel per variant, same
     y-axis;
   - `results/stage24/sens_summary.png`: Sand and Vertisol peak height and peak time against Sand matrix fraction (base, S04,
     S16, S30) and against Vertisol/Loess inert share (base, VL25, VL50).

## 3. Part B — is the split visible in the grey levels?

### Inputs

The raw volumes the segmentations were made from, voxel-aligned with them. Known paths (verify; `N = RB\nnUNet_resources`):
- **Sand samp2:** `N\rehovot_inference_rehovot_samp_2\nifti_predict\rehovot_samp_2_0000.nii.gz`. The Otsu mask is
  `RB\10.5\rehovot_samp_2.npy`.
- **Vertisol A:** `N\bnei_reem_specA_unified\nifti_predict\bnei_reem_specA_unified_0000.nii.gz`. Vertisol B: the same pattern for
  `bnei_reem_specB_unified`.
- **Loess:** use LC (`mishmar_hanegev_Cu011_samp_2_Rec_nlm`, native 15 µm) and find its `nifti_predict` raw volume.
- **Sand samp3:** find its raw volume. If it cannot be found, report it and continue with samp2.

The 5.85 / 8.8 µm Loess label volumes are downsampled and do not share a grid with their raw scans, so skip them here.

**Alignment gate.** nnU-Net outputs are y-flipped relative to the raw NIfTI. Rony's `make_fig7.py`
(`07_analysis_runs/missing_data_20261004/figure7/`) flips them and finds the npy axis order by correlation. Reuse that logic.
For each scan, report the correlation between the pore mask and low grey values after alignment. If the correlation is below 0.5,
stop for that scan.

### Analysis (per scan; inside the core only, excluding a 10-voxel border)

1. **Histograms:**
   - Sand: grey values inside the solid mask. Vertisol / Loess: inside the non-pore phase (labels 0–3), and separately label 1.
   - Each both raw and after a 2-voxel erosion of that phase, which removes grain-edge partial-volume voxels.
   - Also the whole-volume histogram.
2. **Two-population test on the eroded set:**
   - fit 1- and 2-component Gaussian mixtures;
   - report ΔBIC, Ashman's D, and the weight of the minor component.
   - Rule (fixed now): "two populations" if ΔBIC > 10, Ashman's D > 2 and minor weight ≥ 0.05.
3. **Multi-Otsu (3 classes) on the whole volume:** report thresholds, and the fraction of non-pore voxels in the middle class before
   and after the 2-voxel erosion. A middle class that mostly disappears after erosion means it is grain-edge blur, not material.
4. **Preview only, no model runs:**
   - where the rule in step 2 says "two populations", report the matrix fraction a grey-level split would give (darker
     population = porous matrix, brighter = dense solid);
   - save central-slice overlays of that split for each scan (`results/stage24/greysplit_<ID>.png`).
5. **Figure** `results/stage24/grey_histograms.png`: one panel per scan, raw vs eroded histograms, with the multi-Otsu thresholds
   and the GMM components drawn.
6. **Bright coatings:** Rehovot is a red sand. State whether a bright tail (possible Fe-oxide coatings) appears in the Sand solid.

## 4. Write-up

1. **Results file:** `docs/stage_results/STAGE24_RESULTS.md` with:
   - reproduction gate;
   - Part A geometry table and predictions vs outcomes;
   - the variants × outputs table;
   - delivery split;
   - Part B alignment, the test results per scan and the preview fractions;
   - limitations.
2. **Verdict:** answer in plain words.
   - (a) Which Stage 23 conclusions survive all variants, and which flip, and in which variant?
   - (b) Does any scan show a grey-level population that could measure the matrix / inert split, and what fraction would it give?
3. **LOGBOOK** headline, run-time table, commit, push, PR.

## 5. Stop rules

Stop and report if:
- a path is missing (except Sand samp3 raw, see Part B);
- the reproduction gate fails;
- mass balance exceeds 1e-10;
- biology.yaml differs from main;
- an alignment correlation is below 0.5 (stop for that scan only).

Final message:
- the variants × outputs table;
- one line on CONV (do the soils stay distinct when matrix fractions are similar?);
- the Part B result per scan;
- the PR link.
