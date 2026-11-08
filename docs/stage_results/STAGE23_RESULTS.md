# Stage 23 — Loess and soil-specific chemistry

Stage 23 adds three Loess scan specimens and replaces the equal initial OM/B0
setting with the fixed soil-specific chemistry values in the prompt. Arms B+C
cross each structure soil with all three chemistries. Biology, transport,
wetting, `K_density`, time maps, crop selection, and the Sand solid-fines rule
remain fixed. No respiration, biomass, texture, carbonate, or other laboratory
files were opened; model output was not compared with measurements.

## Pre-flight

The Stage 23 branch was created after fetching `origin/main` at
`eecc7b217aaa0388ab41ad441efe10cf8c58c2fa`. The 2026-11-08 LB amendment is
recorded in the Stage 23 prompt on this branch: the LB porosity gate uses the
whole-volume value; its inside-disk value is reported but is not gated.
`configs/biology.yaml` is byte-identical to both `main` and `origin/main`
(normalized SHA-256 `c76d9869dc5847f830ac3caadcce4c68986a9180875f3333504c2482a4b0bf5f`).

The Stage 22 `SA_c0` run was rebuilt from the scan at origin `(400,300,100)`
and compared with its stored Stage 22 curve at `T=3000`. Maximum relative
difference was `0` (gate `<1e-8`). For this reproduction only, the recorded
Stage 22 `K_density` (`0.16445337235927582`) was used; Stage 23 production
runs use the registered common value `0.16445`.

All seven registered paths existed and passed their label/porosity gates.
The label histograms below are whole-volume percentages. The allowed Loess
label set is `{0,1,2,3,5}`; no `4`, `6`, or other values were observed.

| Input | Whole-volume label fractions (%) | Pore % | POM % | Gate |
|---|---|---:|---:|---|
| VA | 0: 82.2763; 1: 2.3146; 2: 0.2403; 5: 15.1688 | 15.1688 | 0.2403 | Pass |
| VB | 0: 75.4465; 1: 2.0962; 2: 0.8170; 5: 21.6403 | 21.6403 | 0.8170 | Pass |
| SA | binary `{0,1}`; 1 is pore | 30.4983 | n/a | Pass |
| SB | binary `{0,1}`; 1 is pore | 40.8479 | n/a | Pass |
| LA | 0: 70.4685; 1: 0.9322; 2: 2.0304; 5: 26.5689 | 26.5689 | 2.0304 | Pass |
| LB | 0: 71.1399; 1: 0.6579; 2: 7.8593; 5: 20.3430 | 20.3430 | 7.8593 | Pass |
| LC | 0: 72.6253; 1: 0.8730; 2: 0.1196; 5: 26.3822 | 26.3822 | 0.1196 | Pass |

**LB amendment detail:** the full-volume fraction `20.3430%` passes against
`20.34% ±1 pp`. The actual inscribed disk (radius `473.5` voxels; 666,956,948
voxels) has `24.7125%` pore and `2.0508%` POM. As amended, this higher
inside-disk value is reported, not gated. Independently, every selected LB
100-cube and the centered 256-cube has all xy corners inside the registered
crop-safe radius (`0.90 × 473.5` voxels).

Central xy/xz/yz slices were saved for all inputs:
[`inputs_VA.png`](../../results/stage23/inputs_VA.png),
[`inputs_VB.png`](../../results/stage23/inputs_VB.png),
[`inputs_SA.png`](../../results/stage23/inputs_SA.png),
[`inputs_SB.png`](../../results/stage23/inputs_SB.png),
[`inputs_LA.png`](../../results/stage23/inputs_LA.png),
[`inputs_LB.png`](../../results/stage23/inputs_LB.png), and
[`inputs_LC.png`](../../results/stage23/inputs_LC.png). LA and LC show no
apparent container wall or exterior-air region in the checked central slices.
The LB image displays both the inscribed-disk and crop-safe boundaries.

### Fixed chemistry

The supplied Corg, Cmic, LOI and intact-core bulk-density means were used
without looking up or opening any laboratory file. Derived volume
concentrations and scaled model densities agree with the prompt values to four
significant figures.

| Soil | Cmic (mg C/kg) | LOI SOM (g/g) | Corg (g C/kg) | Core BD (g/cm³) | Corg (mg C/cm³) | Cmic (µg C/cm³) | OM_density | B0_density |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Sand | 199.44 | 0.016732 | 9.705 | 1.3947 | 13.5356 | 278.1590 | 0.584349 | 0.00880746 |
| Vertisol | 406.61 | 0.048627 | 28.204 | 1.3482 | 38.0246 | 548.1916 | 1.641575 | 0.01735760 |
| Loess | 315.57 | 0.027252 | 15.806 | 1.1344 | 17.9303 | 357.9826 | 0.774076 | 0.01133494 |
| `chem_eq` | — | — | — | — | — | — | 1.0 | 0.0125 |

`K_density=0.16445` was common to every Stage 23 arm. OM remains uniform per
matrix voxel; initial B0 remains proportional to K. Corg is total organic C,
so equal availability is assumed across soils. LOI may overestimate SOM in the
clay-rich Vertisol; the incubated soil was disturbed, whereas the scan and
bulk-density measurements are from intact cores.

## Step 1 — geometry and predictions

The 28 primary 100³ crops comprise Sand 8, Vertisol 8 and Loess 12. The
selected VA/VB/SA/SB origins reproduce Stage 22 exactly; no zero-habitat crops
were rejected. Loess candidate counts were LA 64, LB 333 (circular constraint),
and LC 288; four crops were selected per specimen, with none rejected for zero
habitat. The eight Arm D Sand variants and three Loess 256³ domains were also
measured before any Stage 23 respiration run.

| Structure soil | Porosity | Matrix | Solid | LT d10/d50/d90 (µm) | Largest cluster | Spans x/y/z | S_w | Habitat | Wet habitat fraction |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Sand | 0.320 ± 0.087 | 0.082 ± 0.000 | 0.598 ± 0.087 | 49.7 / 101.4 / 165.2 | 0.977 ± 0.037 | 0.88 / 1.00 / 1.00 | 0.294 ± 0.051 | 0.320 ± 0.087 | 0.104 ± 0.018 |
| Vertisol | 0.164 ± 0.060 | 0.819 ± 0.059 | 0.017 ± 0.006 | 41.3 / 93.9 / 197.1 | 0.710 ± 0.258 | 0.62 / 0.88 / 0.50 | 0.872 ± 0.072 | 0.164 ± 0.060 | 0.265 ± 0.134 |
| Loess | 0.246 ± 0.072 | 0.748 ± 0.073 | 0.006 ± 0.008 | 26.9 / 85.6 / 179.4 | 0.874 ± 0.152 | 0.83 / 1.00 / 0.83 | 0.844 ± 0.094 | 0.246 ± 0.072 | 0.394 ± 0.255 |

Full per-crop, per-specimen, soil-mean and SD metrics (including all Arm D and
domain structures) are in [`step1.json`](../../results/stage23/step1.json).
The selected crop slices are in
[`step1_slices.png`](../../results/stage23/step1_slices.png); the pre-run
prediction record is [`STAGE23_PREDICTIONS.md`](./STAGE23_PREDICTIONS.md).

The Arm D pore-fill geometry has mean resolved porosity `0.2377` versus
`0.3198` for the standard Sand crops; wet habitat fraction falls from `0.1043`
to `0.0360`. Matrix fraction remains `0.0821`, but the added matrix occupies
pore rather than inert-solid voxels. Mean porosities of the centered Loess
256³ domains are LA `0.2957`, LB `0.2706`, and LC `0.2825`.

### Predictions versus outcomes

| Pre-run prediction | Outcome |
|---|---|
| Loess would be drained-film/exchange dominated because its wet-habitat fraction was below one half. | **Wrong:** Arm A Loess net delivery was 75.2% into wet habitat; wet-pore diffusion was the largest single component. |
| Loess would be intermediate between Sand and Vertisol in the Arm A mean response. | **Partly right:** Loess mean peak ratio and peak time were intermediate, but it was not distinct from either soil under the registered pairwise criteria. |
| Chemistry would dominate peak height and cumulative R; structure would dominate peak timing. | **Mixed:** structure dominated R_peak/floor (58.4% vs chemistry 1.2%) and peak time (85.4% vs 5.3%). Chemistry dominated cumulative R per volume (54.6% vs 18.6%) and per mass (49.9% vs 19.9%). |
| Pore-filling fines would broaden or delay Sand response. | **Partly wrong:** the alternative lowered peak ratio but shifted the mean peak fraction slightly earlier, not later. Diffusion share fell and exchange share rose. |

## Step 2 — fixed respiration arms

All 131 registered runs completed: A 28, B 28, C 56, D 16 and E 3.
All 131 were both primary-alive and secondary-alive. Maximum mass-balance relative error was
`1.23e-14` (gate `1e-10`). The longest 100³ run was `177.6 s`, below the
10-minute stop limit. Every run used the CuPy backend.

| Arm / structure soil | n | B_peak/B0 | R_peak/floor | Growth share at peak | Peak fraction of T | R_end/R_peak |
|---|---:|---:|---:|---:|---:|---:|
| A Sand | 8 | 4.493 ± 0.332 | 15.64 ± 0.78 | 0.8805 ± 0.0014 | 0.0430 ± 0.0003 | 0.2888 ± 0.0419 |
| A Vertisol | 8 | 1.630 ± 0.147 | 7.017 ± 0.631 | 0.7827 ± 0.0025 | 0.1449 ± 0.0076 | 0.3094 ± 0.0281 |
| A Loess | 12 | 2.679 ± 0.94 | 11.55 ± 4.1 | 0.8022 ± 0.015 | 0.1159 ± 0.019 | 0.1981 ± 0.079 |
| B Sand | 8 | 4.315 ± 0.19 | 17.21 ± 0.75 | 0.8756 ± 0.0021 | 0.0534 ± 0.0003 | 0.1050 ± 0.032 |
| B Vertisol | 8 | 1.984 ± 0.18 | 8.159 ± 0.733 | 0.7773 ± 0.0020 | 0.1400 ± 0.0070 | 0.3219 ± 0.0281 |
| B Loess | 12 | 2.292 ± 0.83 | 10.07 ± 3.73 | 0.7993 ± 0.013 | 0.1242 ± 0.017 | 0.1899 ± 0.083 |

Mean net delivery fractions (diffusion / exchange; model units are retained in
the per-run JSON) were:

| Arm | Sand | Vertisol | Loess |
|---|---:|---:|---:|
| A | 43.1% / 56.9% | 68.8% / 31.2% | 66.6% / 33.4% |
| B | 43.6% / 56.4% | 68.9% / 31.1% | 66.5% / 33.5% |

Every registered Stage 22-format per-run summary field and every R/B trace
is stored in
[`results/stage23/runs/`](../../results/stage23/runs/). The run manifest and
aggregate summaries are [`run_plan.json`](../../results/stage23/run_plan.json)
and [`step2.json`](../../results/stage23/step2.json).

## Step 3 — Arm A structure-only distinctness

All crops in each compared soil were primary-alive. The registered distinctness
criteria apply to steps 0–1500: at least 75% alive in both soils, soil-mean
distance greater than 0.3, and max within-soil distance less than half the min
cross-soil distance.

| Pair | Mean distance 0–1500 | Max within | Min cross | Distinct? | Mean distance 0–3000 |
|---|---:|---:|---:|---|---:|
| Sand–Vertisol | 0.5544 | 0.1100 | 0.4215 | **Yes** | 0.2031 |
| Sand–Loess | 0.5509 | 0.2745 | 0.4173 | No: within/cross separation | 0.2474 |
| Vertisol–Loess | 0.0950 | 0.2745 | 0.0046 | No: distance and separation | 0.0761 |

The Sand–Vertisol Stage 23 Arm A values match the Stage 22 values at four
decimal places. At full precision, the differences are approximately
`2e-6`–`4e-6`, because the Stage 23 fixed `K_density=0.16445` differs from the
Stage 22 recorded `0.16445337235927582`; the Stage 22 reproduction gate itself
matched exactly using the stored Stage 22 value.

## Step 4 — Arm B and structure × chemistry cross

Arm B distinctness: Sand–Vertisol passes on steps 0–1500 (distance `0.3381`,
max within `0.0855`, min cross `0.2283`); Sand–Loess fails the within/cross
separation criterion (`0.2688` vs `0.2221`); Vertisol–Loess fails both the
mean-distance and separation criteria (`0.0862`, `0.2688` vs `0.0143`).
All crops remain primary-alive.

The table gives crop mean ± SD in each structure-soil × chemistry cell.
Cumulative R per mass divides per-voxel cumulative R by the **structure
soil's** core BD. qCO2 uses the primary `r_lit` mapping.

| Structure × chemistry | n | R_peak/floor | Peak fraction | Cum. R/volume | Cum. R/mass | R_end/R_peak | qCO2 days 10–15 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Sand × chem_S | 8 | 17.21 ± 0.75 | 0.05335 ± 0.00031 | 0.08095 ± 0.00472 | 0.05804 ± 0.00338 | 0.1050 ± 0.0320 | 0.1452 ± 0.0034 |
| Sand × chem_V | 8 | 13.50 ± 0.87 | 0.03522 ± 0.00017 | 0.1267 ± 0.0129 | 0.09085 ± 0.00922 | 0.4919 ± 0.0323 | 0.1278 ± 0.0007 |
| Sand × chem_L | 8 | 15.33 ± 0.68 | 0.04714 ± 0.00017 | 0.09366 ± 0.00676 | 0.06715 ± 0.00485 | 0.1922 ± 0.0388 | 0.1401 ± 0.0023 |
| Vertisol × chem_S | 8 | 5.803 ± 0.533 | 0.1536 ± 0.0057 | 0.03259 ± 0.00284 | 0.02417 ± 0.00211 | 0.3047 ± 0.0283 | 0.2185 ± 0.0018 |
| Vertisol × chem_V | 8 | 8.159 ± 0.733 | 0.1400 ± 0.0070 | 0.09433 ± 0.00817 | 0.06997 ± 0.00606 | 0.3219 ± 0.0281 | 0.2039 ± 0.0031 |
| Vertisol × chem_L | 8 | 5.986 ± 0.542 | 0.1481 ± 0.0068 | 0.04412 ± 0.00382 | 0.03272 ± 0.00283 | 0.3070 ± 0.0282 | 0.2160 ± 0.0021 |
| Loess × chem_S | 12 | 9.952 ± 3.81 | 0.1337 ± 0.0159 | 0.05283 ± 0.0185 | 0.04657 ± 0.0163 | 0.1839 ± 0.0862 | 0.2076 ± 0.0086 |
| Loess × chem_V | 12 | 12.51 ± 4.06 | 0.09906 ± 0.0241 | 0.1386 ± 0.0400 | 0.1222 ± 0.0353 | 0.2322 ± 0.0674 | 0.1829 ± 0.0144 |
| Loess × chem_L | 12 | 10.07 ± 3.73 | 0.1242 ± 0.0174 | 0.07041 ± 0.0239 | 0.06207 ± 0.0210 | 0.1899 ± 0.0833 | 0.2028 ± 0.0098 |

Two-way sums-of-squares shares of total variation (structure, chemistry,
interaction, residual):

| Output | Structure | Chemistry | Interaction | Residual |
|---|---:|---:|---:|---:|
| R_peak/floor | 58.4% | 1.2% | 7.3% | 33.1% |
| Peak fraction of T | 85.4% | 5.3% | 0.9% | 8.5% |
| Cumulative R/volume | 18.6% | 54.6% | 3.7% | 23.1% |
| Cumulative R/mass | 19.9% | 49.9% | 6.0% | 24.2% |
| R_end/R_peak | 16.0% | 26.2% | 35.1% | 22.6% |
| qCO2 days 10–15 | 87.9% | 6.8% | 0.4% | 4.9% |

On steps 0–1500, the mean distance for curves sharing structure but differing
in chemistry was `0.115 ± 0.157`; for curves sharing chemistry but differing
in structure it was `0.394 ± 0.251`. On steps 0–3000 these means were
`0.083 ± 0.116` and `0.197 ± 0.136`, respectively. **Structure moves the
normalized curves more; chemistry explains more variation in cumulative R.**
The cross curves are in
[`step4_cross.png`](../../results/stage23/step4_cross.png).

## Step 5 — side checks

| Sand chemistry | Mean curve distance 0–1500 / 0–3000 | R_peak/floor, solid-fines → pore-fill | Peak fraction, solid-fines → pore-fill | Diffusion share, solid-fines → pore-fill |
|---|---:|---:|---:|---:|
| chem_eq | 0.1189 / 0.0857 | 15.64 → 13.70 | 0.04301 → 0.04214 | 43.0% → 39.2% |
| chem_S | 0.0611 / 0.0441 | 17.21 → 14.67 | 0.05335 → 0.05273 | 43.6% → 41.3% |

The Stage 21 pore-fill rule increases the exchange share correspondingly.
Under Arm A it changes which Sand comparison is closer at 0–3000 (S–V to
S–Loess), but not at 0–1500. Under Arm B it changes neither interval's
S–Vertisol vs S–Loess ordering. The rule is a sensitivity check only and is
not adopted.

| Loess specimen | 256³ domain vs 100³ crop-mean distance, 0–1500 / 0–3000 |
|---|---:|
| LA | 0.01894 / 0.01535 |
| LB | 0.00975 / 0.00978 |
| LC | 0.00638 / 0.00784 |

## Step 6 — time mapping and physiology

The table gives the new Stage 23 classifier followed in parentheses by the
Stage 22 classifier for continuity.

| Arm | Soil | Low | Primary | High |
|---|---|---|---|---|
| A | Sand | rise (rise) | flat (flat) | rise-then-fall (rise-then-fall) |
| A | Vertisol | rise (rise) | rise (flat) | rise-then-fall (rise-then-fall) |
| A | Loess | rise (rise) | rise-then-fall (rise-then-fall) | fall (fall) |
| B | Sand | rise (rise) | rise-then-fall (rise-then-fall) | fall (fall) |
| B | Vertisol | rise (rise) | rise (flat) | rise-then-fall (rise-then-fall) |
| B | Loess | rise (rise) | rise-then-fall (rise-then-fall) | fall (fall) |

Arm B primary-mapping physiological summaries (crop mean ± SD):

| Soil | Cumulative R/volume | Cumulative R/mass | P4 rate/mass/day | qCO2 days 10–15 (d⁻¹) | CUE |
|---|---:|---:|---:|---:|---:|
| Sand | 0.08095 ± 0.00472 | 0.05804 ± 0.00338 | 0.003904 ± 0.000238 | 0.1452 ± 0.0034 | 0.6673 ± 0.0030 |
| Vertisol | 0.09433 ± 0.00817 | 0.06997 ± 0.00606 | 0.005171 ± 0.000478 | 0.2039 ± 0.0031 | 0.6899 ± 0.0012 |
| Loess | 0.07041 ± 0.02386 | 0.06207 ± 0.02103 | 0.004459 ± 0.001369 | 0.2028 ± 0.0098 | 0.6911 ± 0.0019 |

Arm B ranks at the primary mapping are Vertisol > Loess > Sand for cumulative
R/mass, P4 rate/mass/day, and qCO2. Arm A and B crop means/SDs, specimen
summaries, all P1–P4 absolute and R0-normalized means, CUE, and the low/high
rank tables are in [`step6.json`](../../results/stage23/step6.json). The
period figure is [`step6_periods.png`](../../results/stage23/step6_periods.png).
These are model outputs only; no agreement with laboratory respiration is
claimed.

## Verdict

Arm A retains the Stage 22 Sand–Vertisol distinction on the registered
0–1500 interval, but Loess does not form a third distinct structure-only
group under the same criterion. Arm B is also distinct for Sand–Vertisol
only. In the 3×3 cross, structure is the larger main effect for peak timing
and normalized curve movement; chemistry is larger for cumulative R. The
pre-run prediction that chemistry would dominate peak ratio was not supported.
Pore-filling Sand fines reduce peak ratio and wet habitat but do not produce
the predicted delay. The Arm E domain checks remain close to the corresponding
100³ specimen means.

## Limitations

- Chemistry uses fixed lab means; no chemistry uncertainty was propagated.
- Corg is treated as available OM at the same fraction in every soil; LOI may
  overestimate SOM in clay-rich Vertisol.
- The incubated soil was disturbed, whereas bulk density and scans represent
  intact cores.
- At 15 µm, all resolved pores meeting the inherited 10 µm habitat threshold
  are treated as habitat; smaller unresolved pores are not modeled.
- LA and LB are majority-vote label downsampled scans; LC is a native 15 µm
  scan. LC's POM channel is invalid, so labels 2/3 are mapped to matrix and
  only pore-label validity is used.
- Sand fines placement follows the Stage 22 solid-LT rule; the Stage 21
  pore-fill alternative is not adopted.
- Physiological interpretation depends on the fixed `r_lit` time mapping.
- Oxygen limitation is absent. `K_density` is common across all soils and
  arms.
- Crop replicates are subsamples of the registered scans, not independent
  field replicates. Arm B ranks are model rankings only.

## Run times

Per-run elapsed times; Arm E used one GPU sequentially for the larger domains.
All runs used CuPy, and Arms A–D were scheduled on both GPUs.

| Arm | Runs | Mean seconds | Maximum seconds | Sum seconds |
|---|---:|---:|---:|---:|
| A | 28 | 62.8 | 177.6 | 1,757 |
| B | 28 | 60.8 | 147.7 | 1,701 |
| C | 56 | 42.4 | 115.4 | 2,375 |
| D | 16 | 31.7 | 39.2 | 508 |
| E | 3 | 348.6 | 349.9 | 1,046 |
| **Total** | **131** | — | — | **7,387** |

The longest 100³ run was below the 600-second stop limit. The 256³ runs were
about 5.8 minutes each, slightly above the planning estimate but below any
registered stop threshold. Full per-run timings and mass-balance values are
recorded in the run JSON files.
