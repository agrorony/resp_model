# Stage 24 — matrix/inert-solid sensitivity and raw grey levels

## Pre-flight

- Biology gate: passed; normalized SHA-256 `c76d9869dc5847f830ac3caadcce4c68986a9180875f3333504c2482a4b0bf5f`, matching HEAD, main and origin/main.
- The requested initial command outputs were empty:
  - `git --no-pager diff HEAD -- configs/biology.yaml`: *(no output)*
  - `git status --short`: *(no output)*
- The working copy used CRLF line endings; after CRLF-to-LF normalization its bytes matched the committed config. `git stash push -u -m "pre-stage24-leftovers"` returned **No local changes to save**. No change was stashed or used; the apparent difference was only line-ending representation.
- The committed file history was `ce82d87 first commit`; HEAD and main both normalize to the expected hash. The branch was created from fetched `origin/main`.
- All seven Stage 23 segmentations, five raw scans and fixed Stage 23 run artifacts were present. The alignment helper copy used was `E:\PROJECTS\Yael_Mishael\Rony\remote_computer backup\nnUNet4SoilXrayCT\06_reporting\missing_data_20261004\figure7\make_fig7.py` (SHA-256 `520713cc24749d7624b5f1e06877bc98a55e3f39757a79dc56e233a7d20cdc2f`); it is byte-identical to the handoff copy of the referenced `make_fig7.py`.
- Stage 23 crops/origins, chemistry values, `K_density=0.16445`, `T=3000`, the Stage 23 shape classifier, `r_lit=0.36/2.6/9.2`, and biology were unchanged. No lab files were opened.

## Reproduction gate

The 14 base runs (A and B, c0 for each of the seven specimens) were rerun and compared pointwise with Stage 23 R(t). The gate is maximum relative difference `<1e-8`; mass balance must remain `≤1e-10`.

| Runs | Max relative difference | Max mass-balance error | Gate |
|---:|---:|---:|---|
| 14 | 0 | 8.04e-15 | PASS |

## Part A — geometry first, then respiration

Every listed geometry is a mean ± SD across the fixed Stage 23 crop set for that soil. “Solid” is inert SOLID; the Sand target is a fraction of crop volume. Vertisol/Loess inert shares are fractions of the original matrix voxels converted to solid. Unchanged soils in sensitivity variants reuse Stage 23 curves.

| Variant | Soil | Matrix fraction | Inert solid fraction | Pore fraction | Matrix–pore faces / N | Mean matrix-to-pore distance | OM/matrix voxel A | OM/matrix voxel B |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| base | Sand | 0.0821 ± 0.0000 | 0.5981 ± 0.0866 | 0.3198 ± 0.0866 | 0.0859 | 22.28 µm | 12.1786 | 7.11657 |
| base | Vertisol | 0.8193 ± 0.0588 | 0.0169 ± 0.0057 | 0.1638 ± 0.0599 | 0.1155 | 75.09 µm | 1.22669 | 2.01371 |
| base | Loess | 0.7478 ± 0.0731 | 0.0064 ± 0.0077 | 0.2458 ± 0.0721 | 0.2151 | 56.79 µm | 1.3487 | 1.044 |
| S04 | Sand | 0.0400 ± 0.0000 | 0.6402 ± 0.0866 | 0.3198 ± 0.0866 | 0.0569 | 18.59 µm | 25 | 14.6087 |
| S04 | Vertisol | 0.8193 ± 0.0588 | 0.0169 ± 0.0057 | 0.1638 ± 0.0599 | 0.1155 | 75.09 µm | 1.22669 | 2.01371 |
| S04 | Loess | 0.7478 ± 0.0731 | 0.0064 ± 0.0077 | 0.2458 ± 0.0721 | 0.2151 | 56.79 µm | 1.3487 | 1.044 |
| S16 | Sand | 0.1600 ± 0.0000 | 0.5202 ± 0.0866 | 0.3198 ± 0.0866 | 0.1206 | 26.50 µm | 6.25 | 3.65218 |
| S16 | Vertisol | 0.8193 ± 0.0588 | 0.0169 ± 0.0057 | 0.1638 ± 0.0599 | 0.1155 | 75.09 µm | 1.22669 | 2.01371 |
| S16 | Loess | 0.7478 ± 0.0731 | 0.0064 ± 0.0077 | 0.2458 ± 0.0721 | 0.2151 | 56.79 µm | 1.3487 | 1.044 |
| S30 | Sand | 0.3000 ± 0.0000 | 0.3802 ± 0.0866 | 0.3198 ± 0.0866 | 0.1543 | 32.44 µm | 3.33333 | 1.94783 |
| S30 | Vertisol | 0.8193 ± 0.0588 | 0.0169 ± 0.0057 | 0.1638 ± 0.0599 | 0.1155 | 75.09 µm | 1.22669 | 2.01371 |
| S30 | Loess | 0.7478 ± 0.0731 | 0.0064 ± 0.0077 | 0.2458 ± 0.0721 | 0.2151 | 56.79 µm | 1.3487 | 1.044 |
| VL25 | Sand | 0.0821 ± 0.0000 | 0.5981 ± 0.0866 | 0.3198 ± 0.0866 | 0.0859 | 22.28 µm | 12.1786 | 7.11657 |
| VL25 | Vertisol | 0.6145 ± 0.0441 | 0.2217 ± 0.0166 | 0.1638 ± 0.0599 | 0.1133 | 59.51 µm | 1.63559 | 2.68494 |
| VL25 | Loess | 0.5608 ± 0.0548 | 0.1933 ± 0.0185 | 0.2458 ± 0.0721 | 0.2099 | 44.15 µm | 1.79827 | 1.392 |
| VL50 | Sand | 0.0821 ± 0.0000 | 0.5981 ± 0.0866 | 0.3198 ± 0.0866 | 0.0859 | 22.28 µm | 12.1786 | 7.11657 |
| VL50 | Vertisol | 0.4097 ± 0.0294 | 0.4266 ± 0.0308 | 0.1638 ± 0.0599 | 0.1066 | 47.68 µm | 2.45338 | 4.02741 |
| VL50 | Loess | 0.3739 ± 0.0366 | 0.3803 ± 0.0360 | 0.2458 ± 0.0721 | 0.1953 | 34.38 µm | 2.6974 | 2.088 |
| CONV | Sand | 0.3000 ± 0.0000 | 0.3802 ± 0.0866 | 0.3198 ± 0.0866 | 0.1543 | 32.44 µm | 3.33333 | 1.94783 |
| CONV | Vertisol | 0.4097 ± 0.0294 | 0.4266 ± 0.0308 | 0.1638 ± 0.0599 | 0.1066 | 47.68 µm | 2.45338 | 4.02741 |
| CONV | Loess | 0.3739 ± 0.0366 | 0.3803 ± 0.0360 | 0.2458 ± 0.0721 | 0.1953 | 34.38 µm | 2.6974 | 2.088 |

### Predictions versus outcomes

| Prediction | Observed outcome |
|---|---|
| Sand peak height decreases and time becomes later from S04 to S30. | Sand Arm A height 14.02 → 16.07; peak-time fraction 0.03643 → 0.07357. |
| Vertisol/Loess peak height increases and time becomes earlier from base through VL25 to VL50. | Vertisol height [7.017351115563388, 8.922187786316162, 11.616053177615296] and peak-time fraction [0.14492330776925644, 0.1373791263754585, 0.12012337445815272]; Loess height [11.551827648535374, 14.239950333722488, 17.757872112891118] and peak-time fraction [0.11592753139935534, 0.10384016894520397, 0.08519506502167389]. |
| Sand–Vertisol remain distinct in CONV. | CONV Sand–Vertisol distinct: False. |

### Arm A: distances and distinctness

The Stage 23 rule applies at 0–1500 steps: both soils at least 75% alive, soil-mean distance >0.3, and max within-soil distance < half the min cross-soil distance.

| Variant | Pair | Mean distance 0–1500 | Max within | Min cross | Distinct? | Mean distance 0–3000 | Curve source by soil |
|---|---|---:|---:|---:|---|---:|---|
| base | S-V | 0.5544 | 0.1100 | 0.4215 | True | 0.2031 | S: Stage 23 reused; V: Stage 23 reused |
| base | S-L | 0.5509 | 0.2745 | 0.4173 | False | 0.2474 | S: Stage 23 reused; L: Stage 23 reused |
| base | V-L | 0.0950 | 0.2745 | 0.0046 | False | 0.0761 | V: Stage 23 reused; L: Stage 23 reused |
| S04 | S-V | 0.9049 | 0.1080 | 0.7186 | True | 0.4229 | S: Stage 23 reused; V: Stage 23 reused |
| S04 | S-L | 0.8786 | 0.2745 | 0.7148 | True | 0.4549 | S: Stage 23 reused; L: Stage 23 reused |
| S04 | V-L | 0.0950 | 0.2745 | 0.0046 | False | 0.0761 | V: Stage 23 reused; L: Stage 23 reused |
| S16 | S-V | 0.2308 | 0.0749 | 0.1630 | False | 0.1056 | S: Stage 23 reused; V: Stage 23 reused |
| S16 | S-L | 0.2467 | 0.2745 | 0.1614 | False | 0.1203 | S: Stage 23 reused; L: Stage 23 reused |
| S16 | V-L | 0.0950 | 0.2745 | 0.0046 | False | 0.0761 | V: Stage 23 reused; L: Stage 23 reused |
| S30 | S-V | 0.0661 | 0.0358 | 0.0454 | False | 0.0700 | S: Stage 23 reused; V: Stage 23 reused |
| S30 | S-L | 0.0908 | 0.2745 | 0.0470 | False | 0.0570 | S: Stage 23 reused; L: Stage 23 reused |
| S30 | V-L | 0.0950 | 0.2745 | 0.0046 | False | 0.0761 | V: Stage 23 reused; L: Stage 23 reused |
| VL25 | S-V | 0.5233 | 0.1100 | 0.3959 | True | 0.1985 | S: Stage 24; V: Stage 24 |
| VL25 | S-L | 0.5118 | 0.2763 | 0.3792 | False | 0.2394 | S: Stage 24; L: Stage 24 |
| VL25 | V-L | 0.0977 | 0.2763 | 0.0035 | False | 0.0757 | V: Stage 24; L: Stage 24 |
| VL50 | S-V | 0.4763 | 0.1100 | 0.3358 | True | 0.1967 | S: Stage 24; V: Stage 24 |
| VL50 | S-L | 0.4562 | 0.2973 | 0.3183 | False | 0.2348 | S: Stage 24; L: Stage 24 |
| VL50 | V-L | 0.0976 | 0.2973 | 0.0050 | False | 0.0728 | V: Stage 24; L: Stage 24 |
| CONV | S-V | 0.0423 | 0.0358 | 0.0118 | False | 0.0382 | S: Stage 24; V: Stage 24 |
| CONV | S-L | 0.1273 | 0.2973 | 0.0233 | False | 0.0683 | S: Stage 24; L: Stage 24 |
| CONV | V-L | 0.0976 | 0.2973 | 0.0050 | False | 0.0728 | V: Stage 24; L: Stage 24 |

| Variant | Peak-height order (high to low) | Peak-time order (early to late) |
|---|---|---|
| base | S > L > V | S < L < V |
| S04 | S > L > V | S < L < V |
| S16 | S > L > V | S < L < V |
| S30 | S > L > V | S < L < V |
| VL25 | S > L > V | S < L < V |
| VL50 | L > S > V | S < L < V |
| CONV | L > S > V | S < L < V |

### Arm B: periods, shape, cumulative rank, qCO₂

P1–P4 values are mean ± SD of per-voxel R, in model units. Shape classes use the unchanged Stage 23 classifier.

| Variant | r_lit | Soil | P1 | P2 | P3 | P4 | Shape | Curve source |
|---|---|---|---:|---:|---:|---:|---|---|
| base | low | Sand | 0.00024061 ± 1.01e-05, 0.00052867 ± 4.36e-05, 0.00094759 ± 7.96e-05, 0.0016734 ± 0.000113 | rise | Stage 23 reused |
| base | low | Vertisol | 0.00046224 ± 1.33e-05, 0.00083683 ± 4.72e-05, 0.0011688 ± 7.51e-05, 0.0015083 ± 0.000106 | rise | Stage 23 reused |
| base | low | Loess | 0.00027573 ± 1.7e-05, 0.00047486 ± 8.53e-05, 0.00067725 ± 0.000161, 0.00091295 ± 0.000264 | rise | Stage 23 reused |
| base | primary | Sand | 0.0014846 ± 7.84e-05, 0.0024556 ± 0.00013, 0.0021274 ± 0.000153, 0.002094 ± 0.000128 | rise-then-fall | Stage 23 reused |
| base | primary | Vertisol | 0.0013255 ± 8.96e-05, 0.0026094 ± 0.000224, 0.0028216 ± 0.000254, 0.0026815 ± 0.000248 | rise | Stage 23 reused |
| base | primary | Loess | 0.00082633 ± 0.00024, 0.0020668 ± 0.000806, 0.0022035 ± 0.000758, 0.0019457 ± 0.000597 | rise-then-fall | Stage 23 reused |
| base | high | Sand | 0.00207 ± 0.00012, 0.0019193 ± 8.88e-05, 0.0013533 ± 5.33e-05, 0.00072879 ± 8.28e-05 | fall | Stage 23 reused |
| base | high | Vertisol | 0.0023164 ± 0.000195, 0.0023641 ± 0.000225, 0.0016673 ± 0.000175, 0.0011876 ± 0.00013 | rise-then-fall | Stage 23 reused |
| base | high | Loess | 0.0017566 ± 0.000624, 0.0015882 ± 0.000438, 0.0009249 ± 0.000183, 0.00055658 ± 7.39e-05 | fall | Stage 23 reused |
| S04 | low | Sand | 0.00025541 ± 1.03e-05, 0.00060023 ± 4.53e-05, 0.0011378 ± 9.12e-05, 0.0020298 ± 0.000143 | rise | Stage 23 reused |
| S04 | low | Vertisol | 0.00046224 ± 1.33e-05, 0.00083683 ± 4.72e-05, 0.0011688 ± 7.51e-05, 0.0015083 ± 0.000106 | rise | Stage 23 reused |
| S04 | low | Loess | 0.00027573 ± 1.7e-05, 0.00047486 ± 8.53e-05, 0.00067725 ± 0.000161, 0.00091295 ± 0.000264 | rise | Stage 23 reused |
| S04 | primary | Sand | 0.0016319 ± 9.86e-05, 0.0018469 ± 0.000138, 0.0016598 ± 0.000169, 0.0016882 ± 0.000163 | flat | Stage 23 reused |
| S04 | primary | Vertisol | 0.0013255 ± 8.96e-05, 0.0026094 ± 0.000224, 0.0028216 ± 0.000254, 0.0026815 ± 0.000248 | rise | Stage 23 reused |
| S04 | primary | Loess | 0.00082633 ± 0.00024, 0.0020668 ± 0.000806, 0.0022035 ± 0.000758, 0.0019457 ± 0.000597 | rise-then-fall | Stage 23 reused |
| S04 | high | Sand | 0.0017245 ± 0.000133, 0.00161 ± 0.00014, 0.0013222 ± 6.87e-05, 0.00094709 ± 3.07e-05 | fall | Stage 23 reused |
| S04 | high | Vertisol | 0.0023164 ± 0.000195, 0.0023641 ± 0.000225, 0.0016673 ± 0.000175, 0.0011876 ± 0.00013 | rise-then-fall | Stage 23 reused |
| S04 | high | Loess | 0.0017566 ± 0.000624, 0.0015882 ± 0.000438, 0.0009249 ± 0.000183, 0.00055658 ± 7.39e-05 | fall | Stage 23 reused |
| S16 | low | Sand | 0.00022847 ± 7.99e-06, 0.00046428 ± 3.32e-05, 0.00077755 ± 5.63e-05, 0.0012895 ± 7.49e-05 | rise | Stage 23 reused |
| S16 | low | Vertisol | 0.00046224 ± 1.33e-05, 0.00083683 ± 4.72e-05, 0.0011688 ± 7.51e-05, 0.0015083 ± 0.000106 | rise | Stage 23 reused |
| S16 | low | Loess | 0.00027573 ± 1.7e-05, 0.00047486 ± 8.53e-05, 0.00067725 ± 0.000161, 0.00091295 ± 0.000264 | rise | Stage 23 reused |
| S16 | primary | Sand | 0.0012047 ± 5.41e-05, 0.002823 ± 0.000128, 0.0025548 ± 0.000107, 0.0024138 ± 8.87e-05 | rise-then-fall | Stage 23 reused |
| S16 | primary | Vertisol | 0.0013255 ± 8.96e-05, 0.0026094 ± 0.000224, 0.0028216 ± 0.000254, 0.0026815 ± 0.000248 | rise | Stage 23 reused |
| S16 | primary | Loess | 0.00082633 ± 0.00024, 0.0020668 ± 0.000806, 0.0022035 ± 0.000758, 0.0019457 ± 0.000597 | rise-then-fall | Stage 23 reused |
| S16 | high | Sand | 0.0022717 ± 9.44e-05, 0.0020799 ± 6.85e-05, 0.0011413 ± 7.88e-05, 0.00042299 ± 4.96e-05 | fall | Stage 23 reused |
| S16 | high | Vertisol | 0.0023164 ± 0.000195, 0.0023641 ± 0.000225, 0.0016673 ± 0.000175, 0.0011876 ± 0.00013 | rise-then-fall | Stage 23 reused |
| S16 | high | Loess | 0.0017566 ± 0.000624, 0.0015882 ± 0.000438, 0.0009249 ± 0.000183, 0.00055658 ± 7.39e-05 | fall | Stage 23 reused |
| S30 | low | Sand | 0.00021498 ± 5.89e-06, 0.00039214 ± 2.4e-05, 0.00060749 ± 3.86e-05, 0.00091962 ± 4.92e-05 | rise | Stage 23 reused |
| S30 | low | Vertisol | 0.00046224 ± 1.33e-05, 0.00083683 ± 4.72e-05, 0.0011688 ± 7.51e-05, 0.0015083 ± 0.000106 | rise | Stage 23 reused |
| S30 | low | Loess | 0.00027573 ± 1.7e-05, 0.00047486 ± 8.53e-05, 0.00067725 ± 0.000161, 0.00091295 ± 0.000264 | rise | Stage 23 reused |
| S30 | primary | Sand | 0.00086392 ± 3.74e-05, 0.0025689 ± 0.000187, 0.0025919 ± 0.000158, 0.0023605 ± 0.000137 | rise | Stage 23 reused |
| S30 | primary | Vertisol | 0.0013255 ± 8.96e-05, 0.0026094 ± 0.000224, 0.0028216 ± 0.000254, 0.0026815 ± 0.000248 | rise | Stage 23 reused |
| S30 | primary | Loess | 0.00082633 ± 0.00024, 0.0020668 ± 0.000806, 0.0022035 ± 0.000758, 0.0019457 ± 0.000597 | rise-then-fall | Stage 23 reused |
| S30 | high | Sand | 0.0020896 ± 0.000126, 0.001934 ± 0.000102, 0.00097298 ± 2.41e-05, 0.00041198 ± 4.61e-05 | fall | Stage 23 reused |
| S30 | high | Vertisol | 0.0023164 ± 0.000195, 0.0023641 ± 0.000225, 0.0016673 ± 0.000175, 0.0011876 ± 0.00013 | rise-then-fall | Stage 23 reused |
| S30 | high | Loess | 0.0017566 ± 0.000624, 0.0015882 ± 0.000438, 0.0009249 ± 0.000183, 0.00055658 ± 7.39e-05 | fall | Stage 23 reused |
| VL25 | low | Sand | 0.00024061 ± 1.01e-05, 0.00052867 ± 4.36e-05, 0.00094759 ± 7.96e-05, 0.0016734 ± 0.000113 | rise | Stage 24 |
| VL25 | low | Vertisol | 0.0005017 ± 1.79e-05, 0.0010006 ± 6.29e-05, 0.0014572 ± 9.9e-05, 0.0019398 ± 0.000139 | rise | Stage 24 |
| VL25 | low | Loess | 0.00029388 ± 2.17e-05, 0.00055996 ± 0.000108, 0.00083804 ± 0.000209, 0.0011801 ± 0.000353 | rise | Stage 24 |
| VL25 | primary | Sand | 0.0014846 ± 7.84e-05, 0.0024556 ± 0.00013, 0.0021274 ± 0.000153, 0.002094 ± 0.000128 | rise-then-fall | Stage 24 |
| VL25 | primary | Vertisol | 0.0016847 ± 0.000117, 0.0033184 ± 0.000276, 0.0034957 ± 0.000299, 0.003308 ± 0.00028 | rise | Stage 24 |
| VL25 | primary | Loess | 0.0010658 ± 0.000323, 0.0026616 ± 0.000974, 0.0026917 ± 0.000834, 0.0023587 ± 0.000644 | rise-then-fall | Stage 24 |
| VL25 | high | Sand | 0.00207 ± 0.00012, 0.0019193 ± 8.88e-05, 0.0013533 ± 5.33e-05, 0.00072879 ± 8.28e-05 | fall | Stage 24 |
| VL25 | high | Vertisol | 0.0029134 ± 0.000238, 0.0029227 ± 0.000246, 0.0020458 ± 0.00017, 0.0014155 ± 0.000111 | rise-then-fall | Stage 24 |
| VL25 | high | Loess | 0.0022127 ± 0.000736, 0.0019098 ± 0.000455, 0.0010615 ± 0.00015, 0.0005873 ± 7.72e-05 | fall | Stage 24 |
| VL50 | low | Sand | 0.00024061 ± 1.01e-05, 0.00052867 ± 4.36e-05, 0.00094759 ± 7.96e-05, 0.0016734 ± 0.000113 | rise | Stage 24 |
| VL50 | low | Vertisol | 0.00056688 ± 2.76e-05, 0.0012728 ± 9.61e-05, 0.0019502 ± 0.000148, 0.0026865 ± 0.000201 | rise | Stage 24 |
| VL50 | low | Loess | 0.00032378 ± 2.94e-05, 0.0007007 ± 0.000146, 0.0011125 ± 0.000288, 0.0016544 ± 0.000506 | rise | Stage 24 |
| VL50 | primary | Sand | 0.0014846 ± 7.84e-05, 0.0024556 ± 0.00013, 0.0021274 ± 0.000153, 0.002094 ± 0.000128 | rise-then-fall | Stage 24 |
| VL50 | primary | Vertisol | 0.0022875 ± 0.000168, 0.0042616 ± 0.000371, 0.0042796 ± 0.000379, 0.0040589 ± 0.000333 | rise | Stage 24 |
| VL50 | primary | Loess | 0.0014837 ± 0.000461, 0.0034558 ± 0.00113, 0.0032618 ± 0.000851, 0.0028328 ± 0.000623 | rise-then-fall | Stage 24 |
| VL50 | high | Sand | 0.00207 ± 0.00012, 0.0019193 ± 8.88e-05, 0.0013533 ± 5.33e-05, 0.00072879 ± 8.28e-05 | fall | Stage 24 |
| VL50 | high | Vertisol | 0.0037054 ± 0.000315, 0.0036066 ± 0.000268, 0.002505 ± 0.000141, 0.0016421 ± 6.01e-05 | fall | Stage 24 |
| VL50 | high | Loess | 0.0028221 ± 0.000838, 0.002248 ± 0.000387, 0.0011273 ± 0.000123, 0.00053664 ± 0.000161 | fall | Stage 24 |
| CONV | low | Sand | 0.00021498 ± 5.89e-06, 0.00039214 ± 2.4e-05, 0.00060749 ± 3.86e-05, 0.00091962 ± 4.92e-05 | rise | Stage 24 |
| CONV | low | Vertisol | 0.00056688 ± 2.76e-05, 0.0012728 ± 9.61e-05, 0.0019502 ± 0.000148, 0.0026865 ± 0.000201 | rise | Stage 24 |
| CONV | low | Loess | 0.00032378 ± 2.94e-05, 0.0007007 ± 0.000146, 0.0011125 ± 0.000288, 0.0016544 ± 0.000506 | rise | Stage 24 |
| CONV | primary | Sand | 0.00086392 ± 3.74e-05, 0.0025689 ± 0.000187, 0.0025919 ± 0.000158, 0.0023605 ± 0.000137 | rise | Stage 24 |
| CONV | primary | Vertisol | 0.0022875 ± 0.000168, 0.0042616 ± 0.000371, 0.0042796 ± 0.000379, 0.0040589 ± 0.000333 | rise | Stage 24 |
| CONV | primary | Loess | 0.0014837 ± 0.000461, 0.0034558 ± 0.00113, 0.0032618 ± 0.000851, 0.0028328 ± 0.000623 | rise-then-fall | Stage 24 |
| CONV | high | Sand | 0.0020896 ± 0.000126, 0.001934 ± 0.000102, 0.00097298 ± 2.41e-05, 0.00041198 ± 4.61e-05 | fall | Stage 24 |
| CONV | high | Vertisol | 0.0037054 ± 0.000315, 0.0036066 ± 0.000268, 0.002505 ± 0.000141, 0.0016421 ± 6.01e-05 | fall | Stage 24 |
| CONV | high | Loess | 0.0028221 ± 0.000838, 0.002248 ± 0.000387, 0.0011273 ± 0.000123, 0.00053664 ± 0.000161 | fall | Stage 24 |

| Variant | r_lit | Cumulative R/mass rank (high to low) | qCO₂ days 10–15 (S/V/L mean ± SD) | qCO₂ rank (high to low) | Curve source |
|---|---|---|---|---|---|
| base | low | V > S > L | S 0.06401 ± 0.00191, V 0.0371 ± 0.00173, L 0.03543 ± 0.00712 | S > V > L | Stage 24 affected soils; Stage 23 reused otherwise |
| base | primary | V > L > S | S 0.1452 ± 0.00342, V 0.2039 ± 0.00308, L 0.2028 ± 0.00977 | V > L > S | Stage 24 affected soils; Stage 23 reused otherwise |
| base | high | V > S > L | S 0.4684 ± 0.0122, V 0.6887 ± 0.00497, L 0.6892 ± 0.0319 | L > V > S | Stage 24 affected soils; Stage 23 reused otherwise |
| S04 | low | S > V > L | S 0.07037 ± 0.00212, V 0.0371 ± 0.00173, L 0.03543 ± 0.00712 | S > V > L | Stage 24 affected soils; Stage 23 reused otherwise |
| S04 | primary | V > L > S | S 0.1315 ± 0.00371, V 0.2039 ± 0.00308, L 0.2028 ± 0.00977 | V > L > S | Stage 24 affected soils; Stage 23 reused otherwise |
| S04 | high | V > L > S | S 0.4132 ± 0.00445, V 0.6887 ± 0.00497, L 0.6892 ± 0.0319 | L > V > S | Stage 24 affected soils; Stage 23 reused otherwise |
| S16 | low | V > L > S | S 0.05477 ± 0.00161, V 0.0371 ± 0.00173, L 0.03543 ± 0.00712 | S > V > L | Stage 24 affected soils; Stage 23 reused otherwise |
| S16 | primary | V > S > L | S 0.1654 ± 0.00434, V 0.2039 ± 0.00308, L 0.2028 ± 0.00977 | V > L > S | Stage 24 affected soils; Stage 23 reused otherwise |
| S16 | high | V > L > S | S 0.5683 ± 0.0289, V 0.6887 ± 0.00497, L 0.6892 ± 0.0319 | L > V > S | Stage 24 affected soils; Stage 23 reused otherwise |
| S30 | low | V > L > S | S 0.04385 ± 0.00138, V 0.0371 ± 0.00173, L 0.03543 ± 0.00712 | S > V > L | Stage 24 affected soils; Stage 23 reused otherwise |
| S30 | primary | V > L > S | S 0.1898 ± 0.00497, V 0.2039 ± 0.00308, L 0.2028 ± 0.00977 | V > L > S | Stage 24 affected soils; Stage 23 reused otherwise |
| S30 | high | V > L > S | S 0.6448 ± 0.0191, V 0.6887 ± 0.00497, L 0.6892 ± 0.0319 | L > V > S | Stage 24 affected soils; Stage 23 reused otherwise |
| VL25 | low | V > L > S | S 0.06401 ± 0.00191, V 0.0435 ± 0.00191, L 0.042 ± 0.00817 | S > V > L | Stage 24 affected soils; Stage 23 reused otherwise |
| VL25 | primary | V > L > S | S 0.1452 ± 0.00342, V 0.193 ± 0.00394, L 0.1962 ± 0.0112 | L > V > S | Stage 24 affected soils; Stage 23 reused otherwise |
| VL25 | high | V > L > S | S 0.4684 ± 0.0122, V 0.6527 ± 0.00835, L 0.6689 ± 0.0384 | L > V > S | Stage 24 affected soils; Stage 23 reused otherwise |
| VL50 | low | V > L > S | S 0.06401 ± 0.00191, V 0.05235 ± 0.00206, L 0.05158 ± 0.00926 | S > V > L | Stage 24 affected soils; Stage 23 reused otherwise |
| VL50 | primary | V > L > S | S 0.1452 ± 0.00342, V 0.174 ± 0.00455, L 0.1846 ± 0.0124 | L > V > S | Stage 24 affected soils; Stage 23 reused otherwise |
| VL50 | high | V > L > S | S 0.4684 ± 0.0122, V 0.5869 ± 0.0145, L 0.6261 ± 0.0528 | L > V > S | Stage 24 affected soils; Stage 23 reused otherwise |
| CONV | low | V > L > S | S 0.04385 ± 0.00138, V 0.05235 ± 0.00206, L 0.05158 ± 0.00926 | V > L > S | Stage 24 affected soils; Stage 23 reused otherwise |
| CONV | primary | V > L > S | S 0.1898 ± 0.00497, V 0.174 ± 0.00455, L 0.1846 ± 0.0124 | S > L > V | Stage 24 affected soils; Stage 23 reused otherwise |
| CONV | high | V > L > S | S 0.6448 ± 0.0191, V 0.5869 ± 0.0145, L 0.6261 ± 0.0528 | S > L > V | Stage 24 affected soils; Stage 23 reused otherwise |

### Variants × outputs

| Variant | Sand–Vertisol distinct? | Peak-height order | Peak-time order | Arm B shapes at 2.6 (S/V/L) | Arm B cumulative rank at 2.6 | Curve source |
|---|---|---|---|---|---|---|
| base | True (same) | S > L > V (same) | S < L < V (same) | rise-then-fall/rise/rise-then-fall (same) | V > L > S (same) | Stage 24 affected soils; Stage 23 reused otherwise |
| S04 | True (same) | S > L > V (same) | S < L < V (same) | flat/rise/rise-then-fall (different) | V > L > S (same) | Stage 24 affected soils; Stage 23 reused otherwise |
| S16 | False (different) | S > L > V (same) | S < L < V (same) | rise-then-fall/rise/rise-then-fall (same) | V > S > L (different) | Stage 24 affected soils; Stage 23 reused otherwise |
| S30 | False (different) | S > L > V (same) | S < L < V (same) | rise/rise/rise-then-fall (different) | V > L > S (same) | Stage 24 affected soils; Stage 23 reused otherwise |
| VL25 | True (same) | S > L > V (same) | S < L < V (same) | rise-then-fall/rise/rise-then-fall (same) | V > L > S (same) | Stage 24 affected soils; Stage 23 reused otherwise |
| VL50 | True (same) | L > S > V (different) | S < L < V (same) | rise-then-fall/rise/rise-then-fall (same) | V > L > S (same) | Stage 24 affected soils; Stage 23 reused otherwise |
| CONV | False (different) | L > S > V (different) | S < L < V (same) | rise/rise/rise-then-fall (different) | V > L > S (same) | Stage 24 affected soils; Stage 23 reused otherwise |

`same` / `different` compares the listed output with base. The primary mapping is `r_lit=2.6`.

### Delivery split (Arm B)

The first two delivery columns partition total net delivery. Wet-exchange is its share of total delivery; drained-exchange columns show Sand’s drained-pore exchange reliance explicitly.

| Variant | Soil | Diffusion wet / total | Diffusion drained / total | Exchange wet / total | Exchange drained / total | Total diffusion | Total exchange | Drained exchange / exchange | Curve source |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| base | Sand | 23.7% | 19.9% | 12.0% | 44.4% | 43.6% | 56.4% | 78.7% | Stage 23 reused |
| base | Vertisol | 42.0% | 27.0% | 20.8% | 10.2% | 68.9% | 31.1% | 32.9% | Stage 23 reused |
| base | Loess | 49.0% | 17.5% | 26.6% | 6.9% | 66.5% | 33.5% | 20.6% | Stage 23 reused |
| S04 | Sand | 25.6% | 10.4% | 13.0% | 50.9% | 36.0% | 64.0% | 79.6% | Stage 23 reused |
| S04 | Vertisol | 42.0% | 27.0% | 20.8% | 10.2% | 68.9% | 31.1% | 32.9% | Stage 23 reused |
| S04 | Loess | 49.0% | 17.5% | 26.6% | 6.9% | 66.5% | 33.5% | 20.6% | Stage 23 reused |
| S16 | Sand | 23.1% | 31.7% | 11.5% | 33.7% | 54.8% | 45.2% | 74.6% | Stage 23 reused |
| S16 | Vertisol | 42.0% | 27.0% | 20.8% | 10.2% | 68.9% | 31.1% | 32.9% | Stage 23 reused |
| S16 | Loess | 49.0% | 17.5% | 26.6% | 6.9% | 66.5% | 33.5% | 20.6% | Stage 23 reused |
| S30 | Sand | 23.4% | 40.0% | 11.4% | 25.2% | 63.4% | 36.6% | 69.0% | Stage 23 reused |
| S30 | Vertisol | 42.0% | 27.0% | 20.8% | 10.2% | 68.9% | 31.1% | 32.9% | Stage 23 reused |
| S30 | Loess | 49.0% | 17.5% | 26.6% | 6.9% | 66.5% | 33.5% | 20.6% | Stage 23 reused |
| VL25 | Sand | 23.7% | 19.9% | 12.0% | 44.4% | 43.6% | 56.4% | 78.7% | Stage 24 |
| VL25 | Vertisol | 41.3% | 27.6% | 20.4% | 10.7% | 68.8% | 31.2% | 34.5% | Stage 24 |
| VL25 | Loess | 48.6% | 17.7% | 26.3% | 7.4% | 66.3% | 33.7% | 22.0% | Stage 24 |
| VL50 | Sand | 23.7% | 19.9% | 12.0% | 44.4% | 43.6% | 56.4% | 78.7% | Stage 24 |
| VL50 | Vertisol | 39.3% | 29.2% | 19.7% | 11.9% | 68.5% | 31.5% | 37.7% | Stage 24 |
| VL50 | Loess | 47.7% | 18.2% | 25.6% | 8.6% | 65.9% | 34.1% | 25.1% | Stage 24 |
| CONV | Sand | 23.4% | 40.0% | 11.4% | 25.2% | 63.4% | 36.6% | 69.0% | Stage 24 |
| CONV | Vertisol | 39.3% | 29.2% | 19.7% | 11.9% | 68.5% | 31.5% | 37.7% | Stage 24 |
| CONV | Loess | 47.7% | 18.2% | 25.6% | 8.6% | 65.9% | 34.1% | 25.1% | Stage 24 |

Sand’s drained-exchange share is 44.4% in base versus 25.2% in S30.

## Part B — raw grey-level split evidence

Grey histograms cover the whole volume and the core after excluding a 10-voxel border; phase masks and label 1 are shown raw and after 2-voxel erosion. Multi-Otsu thresholds are fit to the whole-volume histogram. The fixed two-population rule is applied to the eroded phase sample. GMMs use deterministic 1-D NumPy EM because the installed scikit-learn binary is incompatible with NumPy in this environment. Preview fractions are candidate darker-population volume fractions, not adopted parameters.

| Scan | Corr(pore, low grey) | Two-population test (ΔBIC, Ashman D, minor wt, pass) | Multi-Otsu thresholds | Middle non-pore fraction before / after erosion | Phase voxels raw / eroded | Candidate matrix fraction of core | Overlay |
|---|---:|---|---|---|---:|---:|---|
| SA | 0.846 | ΔBIC -118.17; D 0.67; minor 0.439; two populations: False | 0.404, 0.6572 | 0.265 / 0.060 | 172971646 / 116680633 | no qualifying split | not generated |
| SB | 0.902 | ΔBIC 31948.89; D 0.25; minor 0.437; two populations: False | 0.2768, 0.5963 | 0.224 / 0.023 | 146706691 / 89808211 | no qualifying split | not generated |
| VA | 0.782 | ΔBIC 10196.56; D 0.03; minor 0.171; two populations: False | -1.774, -0.04361 | 0.138 / 0.066 | 212450293 / 168280459 | no qualifying split | not generated |
| VB | 0.859 | ΔBIC 169569.83; D 0.48; minor 0.066; two populations: False | -1.373, -0.01758 | 0.076 / 0.025 | 194498675 / 160993562 | no qualifying split | not generated |
| LC | 0.864 | ΔBIC 100533.04; D 1.48; minor 0.010; two populations: False | -1.064, 0.1907 | 0.187 / 0.081 | 233872565 / 184121873 | no qualifying split | not generated |

For Sand, the solid histograms show the bright-tail quantiles in `grey_levels.json`; a bright tail can be consistent with dense/coated grains but is not specific evidence for Fe-oxide coatings. No scan failing the 0.5 alignment threshold was analyzed further.

| Sand scan | Solid q99.9 | Top 0.1% grey range |
|---|---:|---:|
| SA | 0.991615 | 0.991615–0.999991 |
| SB | 1 | 1–1 |

A bright tail is present in both Sand solids; samp3’s top 0.1% is clipped at the observed maximum of 1.0, so it is not evidence specific to Fe-oxide coatings.

The Multi-Otsu middle-class fraction decreases after erosion in every scan, by more than half in each case. This is consistent with the class being predominantly grain-edge partial-volume blur, not a distinct material population.

## Verdict

Arm A Sand–Vertisol distinctness flips in S16, S30, CONV. Peak-height order differs from base in VL50, CONV; peak-time order differs in no variant. Arm B primary shape classes differ in S04, S30, CONV; cumulative-R/mass rank differs in S16. Sand drained-exchange delivery share changes from 44.4% in base to 25.2% in S30 (-19.1 percentage points).

CONV: Sand–Vertisol do not remain distinct even though their matrix fractions are similar.

The grey-level verdict is determined by the fixed GMM rule: SA=no qualifying two-population split, SB=no qualifying two-population split, VA=no qualifying two-population split, VB=no qualifying two-population split, LC=no qualifying two-population split. A split is only plausible where all three pre-registered GMM criteria pass; it is not used in the model.

## Limitations

- Matrix fractions are sensitivity inputs, not measured quantities; no variant is adopted.
- Crop replicates are subsamples of the seven registered scans, not independent field replicates.
- Input paths and shapes were recorded, but source-volume checksums were not; the files at those paths were assumed immutable during the workflow.
- The dark/bright interpretation assumes darker material is porous matrix; grey-level populations may also reflect composition, density, reconstruction, noise or partial volume.
- Multi-Otsu is an intensity partition, not a mineralogical or porosity assay; erosion cannot eliminate every partial-volume effect.
- Stage 23 chemistry means and the fixed time mapping were not varied. No laboratory respiration data were opened.
- Sand’s possible bright tail is not diagnostic of Fe-oxide coatings.
- Stage 24 changes matrix amount while preserving total OM, intentionally changing OM per matrix voxel.

## Runtime and delivery log

| Run group | Runs | Mean seconds | Maximum seconds | Sum seconds |
|---|---:|---:|---:|---:|
| base reproduction | 14 | 34.8 | 43.8 | 486.5 |
| S04 | 16 | 28.9 | 32.9 | 462.2 |
| S16 | 16 | 28.9 | 32.4 | 462.6 |
| S30 | 16 | 28.9 | 33.0 | 462.2 |
| VL25 | 40 | 36.5 | 44.8 | 1461.3 |
| VL50 | 40 | 34.4 | 42.2 | 1377.3 |
| CONV | 56 | 32.9 | 43.0 | 1842.3 |

Figures: [sens_curves.png](../../results/stage24/sens_curves.png), [sens_summary.png](../../results/stage24/sens_summary.png), [grey_histograms.png](../../results/stage24/grey_histograms.png). Preview overlays are saved for qualifying scans only.

Branch: `agents/stage24-matrix-sensitivity`. Commits and publication status are recorded in `LOGBOOK.md`.
