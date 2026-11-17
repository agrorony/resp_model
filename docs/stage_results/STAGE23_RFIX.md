# Stage 23 growth-respiration correction

Legacy values are from the committed Stage 23 outputs; corrected values use `growth_R_legacy / Y + maint_R`.
Only the six registered Arm A/B × SA_c0/VA_c0/LC_c0 runs were freshly rerun as gates. The remaining 125 traces were post-hoc corrected; no Stage 23 state was rerun.

## Fresh reproduction and consistency gates

See [`gates/gate.json`](../../results/stage23_rfix/gates/gate.json) and the paired fresh traces in `results/stage23_rfix/gates/`.

## Arm A distinctness and peak order

| Pair | Distance 0–1500, legacy → corrected | Distance 0–3000, legacy → corrected | Criteria, legacy → corrected | Verdict, legacy → corrected |
|---|---:|---:|---|---|
| S-V | 0.554379 → 0.706312 | 0.20309 → 0.282006 | `{'at_least_75_percent_alive_each': True, 'mean_distance_gt_0_3': True, 'max_within_lt_half_min_cross': True}` → `{'at_least_75_percent_alive_each': True, 'mean_distance_gt_0_3': True, 'max_within_lt_half_min_cross': True}` | True → True |
| S-L | 0.550876 → 0.586714 | 0.24742 → 0.251214 | `{'at_least_75_percent_alive_each': True, 'mean_distance_gt_0_3': True, 'max_within_lt_half_min_cross': False}` → `{'at_least_75_percent_alive_each': True, 'mean_distance_gt_0_3': True, 'max_within_lt_half_min_cross': False}` | False → False |
| V-L | 0.0950152 → 0.108688 | 0.0761041 → 0.0840383 | `{'at_least_75_percent_alive_each': True, 'mean_distance_gt_0_3': False, 'max_within_lt_half_min_cross': False}` → `{'at_least_75_percent_alive_each': True, 'mean_distance_gt_0_3': False, 'max_within_lt_half_min_cross': False}` | False → False |

- Peak-height order: **S > L > V → S > L > V**.
- Peak-time order (earliest to latest): **S < L < V → S < L < V**.

Corrected Stage 22 scan-structure comparison for the shared Sand–Vertisol Arm A pair:
- 0–1500: Stage 22 0.70631, Stage 23 0.706312; equal at 4 decimals: True.
- 0–3000: Stage 22 0.282005, Stage 23 0.282006; equal at 4 decimals: True.

## Arm B distinctness and rankings

| Pair | Distance 0–1500, legacy → corrected | Distance 0–3000, legacy → corrected | Verdict, legacy → corrected |
|---|---:|---:|---|
| S-V | 0.338082 → 0.434158 | 0.160331 → 0.209175 | True → True |
| S-L | 0.317774 → 0.388574 | 0.12005 → 0.152156 | False → False |
| V-L | 0.0861694 → 0.0854677 | 0.0822923 → 0.0827829 | False → False |

| Mapping | Metric | Legacy rank | Corrected rank |
|---|---|---|---|
| low | cumulative_R_per_mass | V > S > L | V > S > L |
| low | P4_rate_per_mass_per_day | S > V > L | S > V > L |
| low | qCO2_days10_15 | S > V > L | S > V > L |
| primary | cumulative_R_per_mass | V > L > S | V > L > S |
| primary | P4_rate_per_mass_per_day | V > L > S | V > L > S |
| primary | qCO2_days10_15 | V > L > S | V > L > S |
| high | cumulative_R_per_mass | V > S > L | V > L > S |
| high | P4_rate_per_mass_per_day | V > S > L | V > L > S |
| high | qCO2_days10_15 | L > V > S | L > V > S |

## Shape classes

| Arm | Mapping | Soil | Legacy | Corrected |
|---|---|---|---|---|
| A | low | S | rise | rise |
| A | low | V | rise | rise |
| A | low | L | rise | rise |
| A | primary | S | flat | flat |
| A | primary | V | rise | rise |
| A | primary | L | rise-then-fall | rise-then-fall |
| A | high | S | rise-then-fall | fall |
| A | high | V | rise-then-fall | rise-then-fall |
| A | high | L | fall | fall |
| B | low | S | rise | rise |
| B | low | V | rise | rise |
| B | low | L | rise | rise |
| B | primary | S | rise-then-fall | rise-then-fall |
| B | primary | V | rise | rise |
| B | primary | L | rise-then-fall | rise-then-fall |
| B | high | S | fall | fall |
| B | high | V | rise-then-fall | fall |
| B | high | L | fall | fall |

## Two-way decomposition shares

| Output | Component | Legacy | Corrected |
|---|---|---:|---:|
| R_peak_over_floor | structure | 0.584115 | 0.614294 |
| R_peak_over_floor | chemistry | 0.0120262 | 0.0107015 |
| R_peak_over_floor | interaction | 0.0732457 | 0.0676929 |
| R_peak_over_floor | residual | 0.330613 | 0.307311 |
| peak_time_fraction | structure | 0.853789 | 0.828265 |
| peak_time_fraction | chemistry | 0.0527121 | 0.0760927 |
| peak_time_fraction | interaction | 0.00888025 | 0.00933186 |
| peak_time_fraction | residual | 0.0846185 | 0.08631 |
| cumulative_R_per_volume | structure | 0.186414 | 0.165243 |
| cumulative_R_per_volume | chemistry | 0.545941 | 0.547245 |
| cumulative_R_per_volume | interaction | 0.0367827 | 0.0456443 |
| cumulative_R_per_volume | residual | 0.230863 | 0.241867 |
| cumulative_R_per_mass | structure | 0.199208 | 0.195758 |
| cumulative_R_per_mass | chemistry | 0.499238 | 0.489926 |
| cumulative_R_per_mass | interaction | 0.0597268 | 0.0671644 |
| cumulative_R_per_mass | residual | 0.241827 | 0.247152 |
| R_end_over_peak | structure | 0.160259 | 0.234149 |
| R_end_over_peak | chemistry | 0.262211 | 0.208529 |
| R_end_over_peak | interaction | 0.351235 | 0.261229 |
| R_end_over_peak | residual | 0.226295 | 0.296093 |
| qCO2_days10_15 | structure | 0.878916 | 0.878916 |
| qCO2_days10_15 | chemistry | 0.0679096 | 0.0679096 |
| qCO2_days10_15 | interaction | 0.00369532 | 0.00369532 |
| qCO2_days10_15 | residual | 0.0494789 | 0.0494789 |

## Soil summary, Arm A/B primary mapping

| Arm | Soil | Metric | Legacy mean ± SD | Corrected mean ± SD |
|---|---|---|---:|---:|
| A | S | R_peak_over_floor | 15.6424 ± 0.781 | 36.3379 ± 1.83 |
| A | S | peak_frac | 0.0430143 ± 0.000252 | 0.0420974 ± 0.000236 |
| A | S | growth_share_at_peak | 0.880487 ± 0.00144 | 0.949828 ± 0.000542 |
| A | S | R_end_over_peak | 0.28883 ± 0.0419 | 0.230292 ± 0.0324 |
| A | S | CUE | 0.656604 ± 0.00261 | 0.296604 ± 0.00261 |
| A | S | qCO2 | 0.135779 ± 0.00148 | 0.261448 ± 0.0037 |
| A | S | cumulative_R_per_volume | 0.105038 ± 0.00881 | 0.215213 ± 0.0188 |
| A | S | cumulative_R_per_mass | 0.0753119 ± 0.00631 | 0.154308 ± 0.0134 |
| A | V | R_peak_over_floor | 7.01735 ± 0.631 | 15.2917 ± 1.38 |
| A | V | peak_frac | 0.144923 ± 0.00762 | 0.132961 ± 0.00727 |
| A | V | growth_share_at_peak | 0.782679 ± 0.00248 | 0.903603 ± 0.00168 |
| A | V | R_end_over_peak | 0.309421 ± 0.0281 | 0.30172 ± 0.028 |
| A | V | CUE | 0.690865 ± 0.00117 | 0.330865 ± 0.00117 |
| A | V | qCO2 | 0.213495 ± 0.00233 | 0.455737 ± 0.00584 |
| A | V | cumulative_R_per_volume | 0.0573845 ± 0.00497 | 0.12423 ± 0.011 |
| A | V | cumulative_R_per_mass | 0.0425638 ± 0.00369 | 0.0921454 ± 0.00816 |
| A | L | R_peak_over_floor | 11.5518 ± 4.13 | 25.6267 ± 9.43 |
| A | L | peak_frac | 0.115928 ± 0.0195 | 0.106758 ± 0.0179 |
| A | L | growth_share_at_peak | 0.802166 ± 0.0154 | 0.915048 ± 0.00848 |
| A | L | R_end_over_peak | 0.198079 ± 0.0794 | 0.189596 ± 0.0794 |
| A | L | CUE | 0.690755 ± 0.00111 | 0.330755 ± 0.00111 |
| A | L | qCO2 | 0.197495 ± 0.0113 | 0.415738 ± 0.0283 |
| A | L | cumulative_R_per_volume | 0.0898722 ± 0.0293 | 0.194516 ± 0.0634 |
| A | L | cumulative_R_per_mass | 0.0792245 ± 0.0258 | 0.171471 ± 0.0559 |
| B | S | R_peak_over_floor | 17.2056 ± 0.749 | 39.8524 ± 1.76 |
| B | S | peak_frac | 0.0533511 ± 0.000309 | 0.0521007 ± 0.000463 |
| B | S | growth_share_at_peak | 0.875602 ± 0.00209 | 0.947943 ± 0.000801 |
| B | S | R_end_over_peak | 0.105002 ± 0.032 | 0.0877568 ± 0.0258 |
| B | S | CUE | 0.667279 ± 0.003 | 0.307279 ± 0.003 |
| B | S | qCO2 | 0.145177 ± 0.00342 | 0.284942 ± 0.00856 |
| B | S | cumulative_R_per_volume | 0.0809491 ± 0.00472 | 0.168573 ± 0.0104 |
| B | S | cumulative_R_per_mass | 0.0580405 ± 0.00338 | 0.120867 ± 0.00748 |
| B | V | R_peak_over_floor | 8.1588 ± 0.733 | 17.7326 ± 1.59 |
| B | V | peak_frac | 0.140047 ± 0.00697 | 0.124875 ± 0.00784 |
| B | V | growth_share_at_peak | 0.777325 ± 0.00205 | 0.902429 ± 0.00161 |
| B | V | R_end_over_peak | 0.321877 ± 0.0281 | 0.31251 ± 0.0277 |
| B | V | CUE | 0.689933 ± 0.00118 | 0.329933 ± 0.00118 |
| B | V | qCO2 | 0.203885 ± 0.00308 | 0.431711 ± 0.00771 |
| B | V | cumulative_R_per_volume | 0.0943339 ± 0.00817 | 0.203885 ± 0.018 |
| B | V | cumulative_R_per_mass | 0.0699702 ± 0.00606 | 0.151228 ± 0.0133 |
| B | L | R_peak_over_floor | 10.0723 ± 3.73 | 22.278 ± 8.44 |
| B | L | peak_frac | 0.124152 ± 0.0174 | 0.115511 ± 0.0164 |
| B | L | growth_share_at_peak | 0.799281 ± 0.0131 | 0.912842 ± 0.00749 |
| B | L | R_end_over_peak | 0.189851 ± 0.0833 | 0.182808 ± 0.0829 |
| B | L | CUE | 0.691118 ± 0.00189 | 0.331118 ± 0.00189 |
| B | L | qCO2 | 0.202786 ± 0.00977 | 0.428965 ± 0.0244 |
| B | L | cumulative_R_per_volume | 0.0704142 ± 0.0239 | 0.152599 ± 0.052 |
| B | L | cumulative_R_per_mass | 0.0620718 ± 0.021 | 0.13452 ± 0.0458 |

## Complete Arm B/C structure × chemistry outputs

| Structure × chemistry | Metric | Legacy mean ± SD | Corrected mean ± SD |
|---|---|---:|---:|
| S×chem_S | R_peak_over_floor | 17.2056 ± 0.749 | 39.8524 ± 1.76 |
| S×chem_S | peak_time_fraction | 0.0533511 ± 0.000309 | 0.0521007 ± 0.000463 |
| S×chem_S | cumulative_R_per_volume | 0.0809491 ± 0.00472 | 0.168573 ± 0.0104 |
| S×chem_S | cumulative_R_per_mass | 0.0580405 ± 0.00338 | 0.120867 ± 0.00748 |
| S×chem_S | R_end_over_peak | 0.105002 ± 0.032 | 0.0877568 ± 0.0258 |
| S×chem_S | qCO2_days10_15 | 0.145177 ± 0.00342 | 0.284942 ± 0.00856 |
| S×chem_V | R_peak_over_floor | 13.4999 ± 0.866 | 31.3504 ± 2.05 |
| S×chem_V | peak_time_fraction | 0.0352201 ± 0.000173 | 0.0346365 ± 0.000214 |
| S×chem_V | cumulative_R_per_volume | 0.126706 ± 0.0129 | 0.256062 ± 0.0266 |
| S×chem_V | cumulative_R_per_mass | 0.0908482 ± 0.00922 | 0.183597 ± 0.0191 |
| S×chem_V | R_end_over_peak | 0.491863 ± 0.0323 | 0.378762 ± 0.0261 |
| S×chem_V | qCO2_days10_15 | 0.127825 ± 0.000671 | 0.241562 ± 0.00168 |
| S×chem_L | R_peak_over_floor | 15.3291 ± 0.683 | 35.5194 ± 1.6 |
| S×chem_L | peak_time_fraction | 0.0471407 ± 0.000173 | 0.0461821 ± 0.000356 |
| S×chem_L | cumulative_R_per_volume | 0.0936578 ± 0.00676 | 0.193188 ± 0.0147 |
| S×chem_L | cumulative_R_per_mass | 0.0671526 ± 0.00485 | 0.138516 ± 0.0105 |
| S×chem_L | R_end_over_peak | 0.192243 ± 0.0388 | 0.156944 ± 0.0303 |
| S×chem_L | qCO2_days10_15 | 0.140122 ± 0.00232 | 0.272305 ± 0.00581 |
| V×chem_S | R_peak_over_floor | 5.80265 ± 0.533 | 12.6615 ± 1.17 |
| V×chem_S | peak_time_fraction | 0.153551 ± 0.00567 | 0.144048 ± 0.00509 |
| V×chem_S | cumulative_R_per_volume | 0.0325924 ± 0.00284 | 0.0704038 ± 0.0063 |
| V×chem_S | cumulative_R_per_mass | 0.0241747 ± 0.00211 | 0.0522206 ± 0.00468 |
| V×chem_S | R_end_over_peak | 0.304703 ± 0.0283 | 0.297409 ± 0.0282 |
| V×chem_S | qCO2_days10_15 | 0.218517 ± 0.0018 | 0.468292 ± 0.0045 |
| V×chem_V | R_peak_over_floor | 8.1588 ± 0.733 | 17.7326 ± 1.59 |
| V×chem_V | peak_time_fraction | 0.140047 ± 0.00697 | 0.124875 ± 0.00784 |
| V×chem_V | cumulative_R_per_volume | 0.0943339 ± 0.00817 | 0.203885 ± 0.018 |
| V×chem_V | cumulative_R_per_mass | 0.0699702 ± 0.00606 | 0.151228 ± 0.0133 |
| V×chem_V | R_end_over_peak | 0.321877 ± 0.0281 | 0.31251 ± 0.0277 |
| V×chem_V | qCO2_days10_15 | 0.203885 ± 0.00308 | 0.431711 ± 0.00771 |
| V×chem_L | R_peak_over_floor | 5.98597 ± 0.542 | 13.0454 ± 1.19 |
| V×chem_L | peak_time_fraction | 0.148133 ± 0.00684 | 0.137588 ± 0.00641 |
| V×chem_L | cumulative_R_per_volume | 0.0441162 ± 0.00382 | 0.0952873 ± 0.00845 |
| V×chem_L | cumulative_R_per_mass | 0.0327223 ± 0.00283 | 0.0706774 ± 0.00627 |
| V×chem_L | R_end_over_peak | 0.306991 ± 0.0282 | 0.299755 ± 0.0281 |
| V×chem_L | qCO2_days10_15 | 0.215955 ± 0.00209 | 0.461887 ± 0.00521 |
| L×chem_S | R_peak_over_floor | 9.95207 ± 3.81 | 22.001 ± 8.6 |
| L×chem_S | peak_time_fraction | 0.133683 ± 0.0159 | 0.125597 ± 0.015 |
| L×chem_S | cumulative_R_per_volume | 0.0528331 ± 0.0185 | 0.114766 ± 0.0406 |
| L×chem_S | cumulative_R_per_mass | 0.0465736 ± 0.0163 | 0.101169 ± 0.0358 |
| L×chem_S | R_end_over_peak | 0.183924 ± 0.0862 | 0.177529 ± 0.0855 |
| L×chem_S | qCO2_days10_15 | 0.207581 ± 0.00863 | 0.440952 ± 0.0216 |
| L×chem_V | R_peak_over_floor | 12.5129 ± 4.06 | 27.913 ± 9.43 |
| L×chem_V | peak_time_fraction | 0.0990608 ± 0.0241 | 0.0893909 ± 0.0212 |
| L×chem_V | cumulative_R_per_volume | 0.138614 ± 0.04 | 0.297428 ± 0.0848 |
| L×chem_V | cumulative_R_per_mass | 0.122191 ± 0.0353 | 0.262189 ± 0.0748 |
| L×chem_V | R_end_over_peak | 0.232249 ± 0.0674 | 0.217627 ± 0.0691 |
| L×chem_V | qCO2_days10_15 | 0.182905 ± 0.0144 | 0.379262 ± 0.036 |
| L×chem_L | R_peak_over_floor | 10.0723 ± 3.73 | 22.278 ± 8.44 |
| L×chem_L | peak_time_fraction | 0.124152 ± 0.0174 | 0.115511 ± 0.0164 |
| L×chem_L | cumulative_R_per_volume | 0.0704142 ± 0.0239 | 0.152599 ± 0.052 |
| L×chem_L | cumulative_R_per_mass | 0.0620718 ± 0.021 | 0.13452 ± 0.0458 |
| L×chem_L | R_end_over_peak | 0.189851 ± 0.0833 | 0.182808 ± 0.0829 |
| L×chem_L | qCO2_days10_15 | 0.202786 ± 0.00977 | 0.428965 ± 0.0244 |

## Curve-level structure/chemistry distances

| Horizon | Comparison | Legacy mean ± SD | Corrected mean ± SD |
|---|---|---:|---:|
| 0–1500 | same_structure_different_chemistry | 0.114989 ± 0.157 | 0.105985 ± 0.132 |
| 0–1500 | same_chemistry_different_structure | 0.39439 ± 0.251 | 0.458485 ± 0.282 |
| 0–3000 | same_structure_different_chemistry | 0.082826 ± 0.116 | 0.0695774 ± 0.091 |
| 0–3000 | same_chemistry_different_structure | 0.197472 ± 0.136 | 0.217793 ± 0.124 |

## Arm D and Loess domain side checks

| Check | Legacy | Corrected |
|---|---:|---:|
| Arm D chem_eq mean peak/floor, solid-fines → pore-fill | 15.6424 → 13.7026 | 36.3379 → 31.7395 |
| Arm D chem_eq mean peak fraction, solid-fines → pore-fill | 0.0430143 → 0.042139 | 0.0420974 → 0.0413471 |
| Arm D chem_S mean peak/floor, solid-fines → pore-fill | 17.2056 → 14.6657 | 39.8524 → 33.8445 |
| Arm D chem_S mean peak fraction, solid-fines → pore-fill | 0.0533511 → 0.0527259 | 0.0521007 → 0.051142 |
| Loess LA domain peak/floor | 16.9935 | 38.1649 |
| Loess LA domain peak fraction | 0.0910303 | 0.0853618 |
| Loess LB domain peak/floor | 12.0674 | 26.7726 |
| Loess LB domain peak fraction | 0.107703 | 0.099033 |
| Loess LC domain peak/floor | 9.13705 | 20.1233 |
| Loess LC domain peak fraction | 0.120707 | 0.109703 |

## P1–P4 means by arm, mapping, and soil

The first table is absolute R per voxel; the second is normalized to R(0).

| Arm | Mapping | Soil | Legacy P1 / P2 / P3 / P4 | Corrected P1 / P2 / P3 / P4 |
|---|---|---|---|---|
| A | low | S | 0.000397359 / 0.000996832 / 0.00183707 / 0.00302013 | 0.000638604 / 0.00215755 / 0.00423119 / 0.00708257 |
| A | low | V | 0.000302208 / 0.000491774 / 0.000665548 / 0.0008453 | 0.000402181 / 0.000914586 / 0.00136912 / 0.0018235 |
| A | low | L | 0.000321912 / 0.000602442 / 0.000889456 / 0.00123298 | 0.000451158 / 0.00118694 / 0.00191507 / 0.00275979 |
| A | primary | S | 0.00236107 / 0.002805 / 0.0027445 / 0.00277224 | 0.00541413 / 0.00587599 / 0.00544415 / 0.00533745 |
| A | primary | V | 0.000755495 / 0.00156698 / 0.00174726 / 0.00165965 | 0.0015798 / 0.0034862 / 0.00379608 / 0.00354289 |
| A | primary | L | 0.00110897 / 0.00268849 / 0.00275304 / 0.00244549 | 0.00242925 / 0.00602739 / 0.00591879 / 0.00513164 |
| A | high | S | 0.0026608 / 0.00269397 / 0.00230886 / 0.00168392 | 0.00559771 / 0.00508637 / 0.00425803 / 0.00309831 |
| A | high | V | 0.00139765 / 0.00145102 / 0.00100532 / 0.00071013 | 0.00304791 / 0.00307347 / 0.00211956 / 0.00150346 |
| A | high | L | 0.00225661 / 0.00202167 / 0.00121206 / 0.000742262 | 0.00495317 / 0.00420228 / 0.00250747 / 0.00154841 |
| B | low | S | 0.000240608 / 0.000528672 / 0.000947595 / 0.0016734 | 0.000352164 / 0.0010933 / 0.00213655 / 0.00390366 |
| B | low | V | 0.000462244 / 0.000836831 / 0.00116879 / 0.0015083 | 0.00066426 / 0.0016477 / 0.00249413 / 0.00333253 |
| B | low | L | 0.000275731 / 0.000474859 / 0.00067725 / 0.000912955 | 0.000368916 / 0.000900782 / 0.00142206 / 0.0020088 |
| B | primary | S | 0.00148455 / 0.00245557 / 0.00212737 / 0.00209398 | 0.00340466 / 0.00530535 / 0.00427677 / 0.00411006 |
| B | primary | V | 0.00132552 / 0.00260939 / 0.00282159 / 0.00268149 | 0.00285255 / 0.0057978 / 0.00609517 / 0.00567774 |
| B | primary | L | 0.000826329 / 0.00206681 / 0.00220351 / 0.00194566 | 0.00177327 / 0.00463945 / 0.00476346 / 0.00410528 |
| B | high | S | 0.00206995 / 0.00191925 / 0.00135327 / 0.00072879 | 0.00441663 / 0.00370798 / 0.00256874 / 0.00138939 |
| B | high | V | 0.00231644 / 0.00236414 / 0.00166726 / 0.00118758 | 0.00505798 / 0.00495744 / 0.00347659 / 0.0024923 |
| B | high | L | 0.0017566 / 0.00158823 / 0.000924904 / 0.000556581 | 0.00385477 / 0.00332094 / 0.00192532 / 0.00116722 |

| Arm | Mapping | Soil | Legacy normalized P1 / P2 / P3 / P4 | Corrected normalized P1 / P2 / P3 / P4 |
|---|---|---|---|---|
| A | low | S | 1.58944 / 3.98733 / 7.34828 / 12.0805 | 2.55442 / 8.63019 / 16.9248 / 28.3303 |
| A | low | V | 1.20883 / 1.96709 / 2.66219 / 3.3812 | 1.60872 / 3.65835 / 5.47648 / 7.294 |
| A | low | L | 1.28765 / 2.40977 / 3.55782 / 4.93191 | 1.80463 / 4.74778 / 7.6603 / 11.0391 |
| A | primary | S | 9.44429 / 11.22 / 10.978 / 11.0889 | 21.6565 / 23.5039 / 21.7766 / 21.3498 |
| A | primary | V | 3.02198 / 6.26792 / 6.98906 / 6.63859 | 6.3192 / 13.9448 / 15.1843 / 14.1715 |
| A | primary | L | 4.4359 / 10.754 / 11.0122 / 9.78198 | 9.71702 / 24.1095 / 23.6752 / 20.5266 |
| A | high | S | 10.6432 / 10.7759 / 9.23546 / 6.73567 | 22.3909 / 20.3455 / 17.0321 / 12.3933 |
| A | high | V | 5.59061 / 5.80407 / 4.0213 / 2.84052 | 12.1916 / 12.2939 / 8.47823 / 6.01386 |
| A | high | L | 9.02643 / 8.08668 / 4.84825 / 2.96905 | 19.8127 / 16.8091 / 10.0299 / 6.19363 |
| B | low | S | 1.36594 / 3.00128 / 5.37951 / 9.49991 | 1.99924 / 6.20669 / 12.1292 / 22.1611 |
| B | low | V | 1.33153 / 2.41056 / 3.36679 / 4.34479 | 1.91346 / 4.74633 / 7.18455 / 9.59964 |
| B | low | L | 1.21629 / 2.09467 / 2.98744 / 4.02717 | 1.62734 / 3.97348 / 6.27292 / 8.86111 |
| B | primary | S | 8.42781 / 13.9403 / 12.0771 / 11.8875 | 19.3283 / 30.1185 / 24.2792 / 23.3328 |
| B | primary | V | 3.81827 / 7.51656 / 8.12782 / 7.72426 | 8.21701 / 16.701 / 17.5576 / 16.3552 |
| B | primary | L | 3.64505 / 9.11696 / 9.71998 / 8.58257 | 7.82215 / 20.4653 / 21.0123 / 18.109 |
| B | high | S | 11.7511 / 10.8956 / 7.68254 / 4.13735 | 25.0732 / 21.0502 / 14.5828 / 7.88757 |
| B | high | V | 6.67268 / 6.81011 / 4.80267 / 3.42093 | 14.5699 / 14.2803 / 10.0146 / 7.17929 |
| B | high | L | 7.74859 / 7.00589 / 4.07988 / 2.45516 | 17.0039 / 14.6491 / 8.49287 / 5.14878 |

## What flipped

- Arm B high cumulative_R_per_mass rank: V > S > L → V > L > S; legacy S−L mean gap=0.00308799, combined crop SD=0.0350587 (within one combined SD).
- Arm B high P4_rate_per_mass_per_day rank: V > S > L → V > L > S; legacy S−L mean gap=0.000293511, combined crop SD=0.000810809 (within one combined SD).
- Arm A high S shape: rise-then-fall → fall (legacy P2−P1=3.31764e-05, combined crop SD=0.000301155; old P1/P2 ordering within one SD; corrected P1/P2 means=0.00559771/0.00508637).
- Arm B high V shape: rise-then-fall → fall (legacy P2−P1=4.77077e-05, combined crop SD=0.000298008; old P1/P2 ordering within one SD; corrected P1/P2 means=0.00505798/0.00495744).
- Arm A Stage 21 pore-fill ordering at 3000 steps: True → False; base closest pair S-V → S-L, pore-fill closest pair S-L → S-L; legacy S–V/S–L distance gap 0.0134378, corrected gap 0.0884881. No registered noise interval was calculated for this single mean-curve ordering.

## Numbers to update in the appendix

Replace the Stage 23 Step 2 and Step 4 peak/floor, growth-share, peak-fraction, end/peak, cumulative-R, qCO2, and decomposition numbers with the complete old-to-new tables above. Replace the Step 5 Arm D and domain peak/floor and peak-time numbers, the Step 6 P1–P4 means, shape classes, CUE, qCO2 and cumulative-R/rank tables, and the peak-order/curve-distance figures. Biomass, habitat, delivery, and substrate-uptake results are unchanged.
