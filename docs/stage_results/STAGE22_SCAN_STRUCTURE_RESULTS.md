# Stage 22 — scan-derived structure (Sand and Vertisol only)

The scan-structure branch of Stage 22 tests whether the Stage 21 respiration
distinction survives when geometry comes from segmented µCT crop volumes.
Biology, transport, exchange, wetting, equal-total densities, liveness rules,
and distance metric were inherited unchanged. No respiration, biomass, texture,
CaCO3, or other laboratory data were read or used for validation.

The earlier equal-wetting Stage 22 outputs remain intact. This run uses the
separate namespace `results/stage22_scan_structure/` and the runner
`stage22_scan_structure_run.py`; derived crop arrays are local-only and ignored
by Git.

## Pre-flight and Step 0

- The working branch was created from `main` at
  `72b45fe534cb3f9f294c51d2bee8eb3aa4f1a692`, identical to the synchronized
  `origin/main` ref. `configs/biology.yaml` matched `main` byte-for-byte
  (normalized SHA-256 `c76d9869dc5847f830ac3caadcce4c68986a9180875f3333504c2482a4b0bf5f`).
- The fresh Stage 21 Sand seed-0 reproduction (OM 1.0, T=1500, 256^3, 5 µm)
  matched the stored respiration curve exactly: max relative difference 0;
  mass-balance relative error `1.02e-15`.
- The two Vertisol inputs were 650^3 uint8 label volumes; the two Sand inputs
  were 650^3 boolean masks. The central xy/xz/yz slices showed no apparent core
  wall or exterior air. The volumes were not flipped.
- **Voxel metadata note:** both segmented NIfTI headers report a 1.0 mm identity
  affine, which conflicts with the scan voxel sizes fixed in the prompt. As
  instructed, the model used the registered table values (VA/SA 15.000149 µm;
  VB/SB 15.034357 µm). Physical length-scale interpretations therefore depend
  on those registered values, not the NIfTI header.

| Input | Labels / values | Pore fraction | POM (labels 2+3) | Gate |
|---|---|---:|---:|---|
| VA | 0: 82.2763%, 1: 2.3146%, 2: 0.2403%, 5: 15.1688% | 15.1688% | 0.2403% | pass |
| VB | 0: 75.4465%, 1: 2.0962%, 2: 0.8170%, 5: 21.6403% | 21.6403% | 0.8170% | pass |
| SA | boolean {False, True}; 1 = pore | 30.4983% | n/a | pass |
| SB | boolean {False, True}; 1 = pore | 40.8479% | n/a | pass |

No Vertisol labels outside `{0,1,2,3,5}` were present. Input slice figures and
the machine-readable preflight record are in
`results/stage22_scan_structure/`.

## Phase mapping and Step 1 geometry

- Vertisol: label 5 = PORE; labels 0, 2, 3 = MATRIX; label 1 = inert SOLID.
  POM positions therefore receive the same per-matrix-voxel OM as the rest of
  the matrix.
- Sand: label 1 = PORE. The model-derived 0.082111 fines fraction was carved
  from SOLID in ascending solid local thickness; distance to a pore and seeded
  randomness break ties. Thresholds for all eight crops are recorded in
  `step1.json`.
- The 100^3 crop grid had 216 non-overlapping origins per specimen. Four crops
  were selected per specimen with seed 22; no zero-habitat crops were rejected.
  The domain-size crops were all 256^3 at origin `(197,197,197)`.
- Six controls were generated as registered: S15/V15, seed offsets 0–2,
  100^3 at 15 µm. Synthetic Vertisol median plate thickness is 40 µm, or
  2.67 voxels at 15 µm.

| Group | Porosity | Matrix | Solid | Pore LT d10/d50/d90 (µm) | Largest cluster | Spans x/y/z | S_w | Habitat | Wet habitat |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| VA | 0.131 | 0.854 | 0.016 | 30.0 / 69.4 / 159.4 | 0.538 | 0.25 / 0.75 / 0.25 | 0.913 | 0.131 | 0.354 |
| VB | 0.197 | 0.785 | 0.018 | 52.6 / 118.4 / 234.9 | 0.882 | 1.00 / 1.00 / 0.75 | 0.830 | 0.197 | 0.175 |
| Vertisol mean | 0.164 | 0.819 | 0.017 | 41.3 / 93.9 / 197.1 | 0.710 | 0.625 / 0.875 / 0.50 | 0.872 | 0.164 | 0.265 |
| SA | 0.247 | 0.082 | 0.671 | 46.9 / 93.8 / 151.9 | 0.958 | 0.75 / 1.00 / 1.00 | 0.334 | 0.247 | 0.111 |
| SB | 0.392 | 0.082 | 0.526 | 52.6 / 109.0 / 178.5 | 0.997 | 1.00 / 1.00 / 1.00 | 0.255 | 0.392 | 0.098 |
| Sand mean | 0.320 | 0.082 | 0.598 | 49.7 / 101.4 / 165.2 | 0.977 | 0.875 / 1.00 / 1.00 | 0.294 | 0.320 | 0.104 |

Per-crop geometry plus means and SDs by specimen and soil are in
`results/stage22_scan_structure/step1.json`. The predictions were committed in
`STAGE22_SCAN_STRUCTURE_STEP1_PREDICTIONS.md` before respiration runs.

**Prediction check:**

| Prediction | Outcome |
|---|---|
| Sand would be drained-film/exchange dominated. | **Right:** 68.2% of reported net delivery was into drained habitat; drained exchange was the largest single delivery component. |
| Vertisol would shift away from Stage 21 wet-pore diffusion toward drained-film/exchange. | **Wrong:** despite wet habitat falling from Stage 21's 73.8% to 26.5%, wet habitat still received 63.4% of net delivery; wet-pore diffusion remained the largest component. |
| Scan S–V distance would be smaller than Stage 21; Sand would peak first. | **Right:** distance was 0.554 over steps 0–1500 versus 1.227 in Stage 21; soil-mean peaks were at steps 129 (Sand) and 434 (Vertisol). |
| Relative to 15 µm synthetic controls, scan Sand would be higher/earlier and scan Vertisol lower/later. | **Right:** the curves and distances in Step 4 show both predicted shifts. |

## Step 2 — OM 1.0, T=3000

Twenty-six runs were made: 16 main crops, six 15 µm controls, and the four
256^3 domain-size checks. The prompt's Step 2 count says 22, while its Step 4
and compute-budget sections require the four domain-size runs; these were
included so the requested Step 4 comparison could be completed.

All 16 main crops were primary- and secondary-alive. The maximum mass-balance
relative error across all 26 runs was `1.04e-14` (below the `1e-10` stop limit).

| Soil | B_peak/B0 | R_peak/floor | Growth share at peak | Peak fraction of T | R_end/R_peak |
|---|---:|---:|---:|---:|---:|
| Sand, 8 crops | 4.493 ± 0.332 | 15.643 ± 0.781 | 0.8805 ± 0.0014 | 0.0430 ± 0.0003 | 0.2888 ± 0.0419 |
| Vertisol, 8 crops | 1.630 ± 0.147 | 7.017 ± 0.631 | 0.7827 ± 0.0025 | 0.1449 ± 0.0076 | 0.3094 ± 0.0281 |

Mean net delivery totals across the eight crops of each soil, in model units:

| Soil | Diffusion to wet | Diffusion to drained | Exchange to wet | Exchange to drained |
|---|---:|---:|---:|---:|
| Sand | 179,070 | 200,804 | 101,240 | 400,104 |
| Vertisol | 213,231 | 133,591 | 106,463 | 50,740 |

Every per-run Stage 21 summary and R/B trace is in
`results/stage22_scan_structure/runs/`.

## Step 3 — distinctness

**PASS on steps 0–1500**, the Stage 21-comparable interval:

| Gate statistic | Result |
|---|---:|
| Primary-alive crops | Sand 8/8; Vertisol 8/8 |
| Soil-mean S–V distance | 0.5544 (> 0.3) |
| Max within-soil crop-pair distance | 0.1100 |
| Min cross-soil crop-pair distance | 0.4215 |
| Max within / min cross | 0.261 (< 0.5) |
| 90th percentile within-soil distance | 0.0550 |

On steps 0–3000 the mean S–V distance falls to 0.2031; max within-soil
distance is 0.0935 versus a minimum cross-soil distance of 0.1500. Thus the
registered comparable-prefix gate passes, but the distance and separation
criterion do not persist over the full 3000-step horizon. This is included as
a material qualification, not hidden by the prefix verdict.

## Step 4 — resolution and domain-size comparisons

| Distance | Steps 0–1500 | Steps 0–3000 |
|---|---:|---:|
| S scan – S15 | 0.2365 | 0.4353 |
| V scan – V15 | 0.1906 | 0.1544 |
| S15 – V15 | 1.0603 | 0.9672 |
| S scan – V scan | 0.5544 | 0.2031 |
| S15 – Stage 21 S5 | 0.4435 | n/a (Stage 21 reference ends at 1500) |
| V15 – Stage 21 V5 | 0.0384 | n/a (Stage 21 reference ends at 1500) |

The structure-source effect at 15 µm (mean scan-vs-control distance 0.2949)
exceeded the synthetic resolution effect (mean 5-vs-15 µm distance 0.2410).
The 2×2 S–V distances on steps 0–1500 were: synthetic 5 µm 1.2274,
synthetic 15 µm 1.0603, scan 15 µm 0.5544. Domain-size distances from each
256^3 crop to its specimen's 100^3 mean were VA 0.0050, VB 0.0036, SA 0.0509,
and SB 0.0435.

## Step 5 — biology sensitivity

Step 3 passed, so the registered 32 runs were performed. Representative crops
closest to their soil means were SA_c2 and VB_c3. Across the 16 paired
perturbations, peak-amplitude order and peak-time order were unchanged in all
16; both soils remained primary-alive in 15/16. Vertisol was not primary-alive
for `Y × 0.5`. S–V distances below are on steps 0–1500.

| Parameter | Distance at 0.5× | Both alive? | Distance at 1.5× | Both alive? |
|---|---:|---|---:|---|
| r_max | 0.312 | yes | 0.673 | yes |
| Ks | 0.677 | yes | 0.512 | yes |
| Y | 0.480 | **no** | 0.508 | yes |
| m0 | 0.834 | yes | 0.390 | yes |
| m_s | 0.660 | yes | 0.525 | yes |
| k_dis | 0.496 | yes | 0.524 | yes |
| k_ex | 0.500 | yes | 0.613 | yes |
| K_density | 0.675 | yes | 0.529 | yes |

No perturbed value was adopted; all overrides were in-memory.

## Step 6 — time mapping and physiological outputs

The primary `r_lit=2.6 d^-1` shape class is **flat** for both soils. Across the
three registered mappings, both soils are rise at low, flat at primary, and
rise-then-fall at high. Period means and S/V period ratios for each mapping,
plus per-specimen and crop data, are in `results/stage22_scan_structure/step6.json`.

| r_lit | Soil | Shape | P1 / P2 / P3 / P4 (R/R0) | Cumulative R (model units) | Cumulative R / initial OM | qCO2 (d^-1) | CUE |
|---|---|---|---|---:|---:|---:|---:|
| 0.36 | Sand | rise | 1.589 / 3.987 / 7.348 / 12.081 | 9,247.7 ± 646.5 | 0.00925 ± 0.00065 | 0.0695 ± 0.0016 | 0.7171 ± 0.0024 |
| 0.36 | Vertisol | rise | 1.209 / 1.967 / 2.662 / 3.381 | 3,265.9 ± 189.7 | 0.00327 ± 0.00019 | 0.0315 ± 0.0015 | 0.6345 ± 0.0096 |
| 2.6 | Sand | flat | 9.444 / 11.220 / 10.978 / 11.089 | 105,038.8 ± 8,805.8 | 0.1050 ± 0.0088 | 0.1358 ± 0.0015 | 0.6566 ± 0.0026 |
| 2.6 | Vertisol | flat | 3.022 / 6.268 / 6.989 / 6.639 | 57,384.6 ± 4,974.8 | 0.0574 ± 0.0050 | 0.2135 ± 0.0023 | 0.6909 ± 0.0012 |
| 9.2 | Sand | rise-then-fall | 10.643 / 10.776 / 9.236 / 6.736 | 313,762.9 ± 14,727.9 | 0.3138 ± 0.0147 | 0.4184 ± 0.0025 | 0.6076 ± 0.0027 |
| 9.2 | Vertisol | rise-then-fall | 5.591 / 5.804 / 4.021 / 2.841 | 152,385.6 ± 14,321.6 | 0.1524 ± 0.0143 | 0.7216 ± 0.0026 | 0.6819 ± 0.0006 |

qCO2 Vertisol/Sand ratios are 0.453 (low), 1.572 (primary), and 1.725 (high).
CUE includes maintenance respiration; it is a whole-community quantity and
can be low or negative when uptake is small.

## Verdict

The scan-backed model passes the Stage 21-comparable distinctness gate at
OM=1.0 and T=3000 with all main crops primary-alive. Scan structure moves the
curves more than the tested synthetic resolution change at 15 µm, and the
domain-size checks are close to their corresponding crop means. However, the
distinctness is not maintained across the full 3000-step interval, the
Vertisol's delivered substrate remains wet-diffusion dominated, and predicted
period shapes depend strongly on the literature time mapping. These are
model-only outcomes; no agreement with laboratory observations is claimed.

## Limitations

- The scanned cores are intact; laboratory incubations used disturbed soil.
- At 15 µm, every visible pore is habitat under the 10 µm lower threshold;
  pores between 10 and 15 µm cannot be resolved.
- Sand fines placement by local thickness is a registered model choice.
- Vertisol POM is treated as matrix, not as its own class.
- The NIfTI affine/zoom metadata conflicts with the registered voxel sizes;
  physical scale is conditional on the prompt's fixed values.
- Time mapping depends on r_lit; the shape classes change from rise (low), to
  flat (primary), to rise-then-fall (high) for both soils.
- No oxygen limitation is modeled.
- OM density is a single soil-level value (1.0).
- The main results use 100^3 crops; only one 256^3 domain-size check was run
  per specimen.
- Only two specimens per soil are represented. The Step 3 full-horizon
  distinctness criterion does not pass.

## Run times

Per-run worker times are recorded by the runner. Aggregate worker time is
reported rather than inferred wall time; the Step 2 and Step 5 jobs were
scheduled across two GPUs.

| Step | Work | Runs | Aggregate worker time | Per-run range |
|---|---|---:|---:|---:|
| 0 | Stage 21 reproduction | 1 | not instrumented | n/a |
| 1 | Crop/control geometry | 26 structures | not instrumented | n/a |
| 2 | Main crops, domain checks, controls | 26 | 1,897.2 GPU-s (31.6 GPU-min) | 33.8–278.1 s |
| 3 | Distinctness metrics | post-processing | not instrumented | n/a |
| 4 | Comparisons and curves | post-processing | not instrumented | n/a |
| 5 | Biology sensitivity | 32 | 1,461.8 GPU-s (24.4 GPU-min) | 37.4–55.7 s |
| 6 | Time mapping and physiology | post-processing | not instrumented | n/a |

Within Step 2, the 100^3 runs took at most 56.7 s (the 600 s limit was not
approached); the four 256^3 domain checks took 180.5–278.1 s each, and the
six controls took 33.8–45.9 s each.
