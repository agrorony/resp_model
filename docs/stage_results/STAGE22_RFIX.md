# Stage 22 scan-structure respiration correction

All 58 stored Stage 22 scan-structure traces were post-hoc corrected from `growth_R` and `maint_R`; no Stage 22 or earlier simulation was rerun. The Stage 22 normalized-curve analysis was rerun from corrected traces.

| Metric | Legacy | Corrected |
|---|---:|---:|
| Sand–Vertisol mean distance, 0–1500 | 0.554377 | 0.70631 |
| Sand–Vertisol mean distance, 0–3000 | 0.203087 | 0.282005 |
| Maximum within-soil distance, 0–1500 | 0.110019 | 0.0989515 |
| Minimum cross-soil distance, 0–1500 | 0.421508 | 0.550747 |
| Registered distinctness verdict | True | True |

The Step 4 comparison against Stage 21 synthetic controls remains mixed-era: the Stage 22 scan curves are corrected, while Stage 21 control traces are retained unchanged, as the prompt permits correction of Stage 22 only. Its distances are not a same-formula comparison.

## Step 4 curve-distance comparisons

| Horizon | Comparison | Legacy | Corrected |
|---|---|---:|---:|
| 0–1500 | S_scan_vs_S15 | 0.236518 | 0.187916 |
| 0–1500 | V_scan_vs_V15 | 0.190595 | 0.218947 |
| 0–1500 | S15_vs_V15 | 1.0603 | 1.01129 |
| 0–1500 | S_scan_vs_V_scan | 0.554377 | 0.70631 |
| 0–1500 | S15_vs_S5 | 0.443479 | 0.622112 |
| 0–1500 | V15_vs_V5 | 0.0384361 | 0.0817483 |
| 0–3000 | S_scan_vs_S15 | 0.435321 | 0.338451 |
| 0–3000 | V_scan_vs_V15 | 0.154426 | 0.167562 |
| 0–3000 | S15_vs_V15 | 0.967214 | 0.870859 |
| 0–3000 | S_scan_vs_V_scan | 0.203087 | 0.282005 |

Mean scan-source effect at 15 µm: 0.213557 → 0.203431.
Mean resolution effect in synthetic controls: 0.240958 → 0.35193.
The larger shift remains **resolution → resolution**.

## Time-mapped summaries

| Mapping | Soil | Shape | CUE | qCO2/day | Period means R/R0 (P1–P4) |
|---|---|---|---:|---:|---|
| low | sand | rise → rise | 0.7171 → 0.3571 | 0.0694768 → 0.162892 | 1.5894 / 3.9873 / 7.3483 / 12.081 → 2.5544 / 8.6302 / 16.925 / 28.331 |
| low | vertisol | rise → rise | 0.634464 → 0.274464 | 0.031476 → 0.0678899 | 1.2088 / 1.9671 / 2.6622 / 3.3812 → 1.6087 / 3.6583 / 5.4765 / 7.294 |
| primary | sand | flat → flat | 0.656604 → 0.296604 | 0.135779 → 0.261449 | 9.4444 / 11.22 / 10.978 / 11.089 → 21.657 / 23.504 / 21.777 / 21.35 |
| primary | vertisol | flat → flat | 0.690865 → 0.330865 | 0.213495 → 0.455738 | 3.022 / 6.2679 / 6.9891 / 6.6386 → 6.3192 / 13.945 / 15.184 / 14.172 |
| high | sand | rise-then-fall → fall | 0.607641 → 0.247641 | 0.418415 → 0.770038 | 10.643 / 10.776 / 9.2355 / 6.7356 → 22.391 / 20.346 / 17.032 / 12.393 |
| high | vertisol | rise-then-fall → rise-then-fall | 0.681939 → 0.321939 | 0.721632 → 1.52808 | 5.5906 / 5.8041 / 4.0213 / 2.8405 → 12.192 / 12.294 / 8.4782 / 6.0139 |

## Sensitivity checks from stored Stage 22 traces

The 32 existing sensitivity runs were post-hoc corrected using each job's configured Y; no sensitivity run was repeated.

| Parameter | Factor | Legacy distance | Corrected distance | Peak amplitude order unchanged, legacy → corrected | Peak time order unchanged, legacy → corrected |
|---|---:|---:|---:|---|---|
| r_max | 0.5 | 0.31162 | 0.428442 | True → True | True → True |
| r_max | 1.5 | 0.672672 | 0.823125 | True → True | True → True |
| Ks | 0.5 | 0.67715 | 0.820651 | True → True | True → True |
| Ks | 1.5 | 0.511894 | 0.683186 | True → True | True → True |
| Y | 0.5 | 0.479666 | 0.647459 | True → True | True → True |
| Y | 1.5 | 0.508102 | 0.624024 | True → True | True → True |
| m0 | 0.5 | 0.834047 | 0.911726 | True → True | True → True |
| m0 | 1.5 | 0.390313 | 0.574201 | True → True | True → True |
| m_s | 0.5 | 0.659711 | 0.840844 | True → True | True → True |
| m_s | 1.5 | 0.524691 | 0.683871 | True → True | True → True |
| k_dis | 0.5 | 0.495862 | 0.759365 | True → True | True → True |
| k_dis | 1.5 | 0.524085 | 0.598428 | True → True | True → True |
| k_ex | 0.5 | 0.500164 | 0.666811 | True → True | True → True |
| k_ex | 1.5 | 0.613423 | 0.771359 | True → True | True → True |
| K_density | 0.5 | 0.675004 | 0.811871 | True → True | True → True |
| K_density | 1.5 | 0.529089 | 0.683472 | True → True | True → True |
