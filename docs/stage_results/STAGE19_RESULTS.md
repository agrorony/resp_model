# Stage 19 RESULTS -- Stage 8 revived (equal OM + per-soil coating fraction): alive, but NOT distinct

`stage19_run.py`, `configs/soil_{sand,vertisol}_stage19.yaml`, `results/stage19/`.
**Two soils only: Sand and Vertisol.** Loess was not run in any form.

## 0. Verdict

**FAIL at Step 3. Both soils are alive at n=128 with identical per-cell
densities, but their curves are not distinct: min pairwise distance 0.238
against the 0.3 bar. Step 4 (robustness) was not run.**

- **Liveness is fixed, and the fix is grid-independent.** With totals as
  densities, the selected density gives the same regime at n=40 and n=128
  (distance 0.223 vs 0.238; Sand B_peak/B0 3.02 vs 2.89; Vertisol 1.50 vs
  1.48). Stage 8's n=40 / n=128 disagreement is gone.
- **Liveness and distinctness trade off along the only tuning axis
  allowed.** Distance falls monotonically as OM_density rises: 0.733 -> 0.365
  -> 0.223 -> 0.101 over the four densities where anything grows. Both soils
  are alive only at the two highest densities, and there the curves are
  already converging.
- **What distinctness remains is carried by the per-soil coating fraction,
  not by structure.** Setting both soils to coating_fraction 0.5 cuts the
  distance from 0.238 to 0.070 and reverses which soil peaks first.

## 1. What was run (rules R1-R6, fixed in `stage19_run.py` before any run)

- **Soils:** Sand and Vertisol, literature PSDs (`psd_parametric`, lambda
  2.0, seeds 2 / 3), `matric_regions` wetting (d_cut 50 um).
- **Biology** byte-identical (SHA-256 asserted on every step). **One
  dissolution rate:** `k_dis_slow = k_dis = 0.05`. No enzyme.
  **`om_fine_bias = 0`**: Stage 8 also weighted OM *within* each pool by
  d^-2.5. R3 forbids that, so OM is uniform within the coating pool and
  within the trapped pool. This is a deliberate difference from the
  original Stage 8.
- **OM placement:** `om_mode: coating_trapped`, coating = porous-matrix
  cells within 2 geodesic steps of a macropore. `coating_fraction` is an
  explicit parameter in each soil's new config: Sand 0.8, Vertisol 0.3
  (Stage 8's values). It is in separate `*_stage19.yaml` files so that
  earlier stages reading `configs/soil_<s>.yaml` are unaffected.
- **Equal densities, asserted on every run:**
  - sum(OM) = OM_density x N
  - sum(B0) = 0.0125 x N
  - sum(K_t) = K_density x N, on the effective capacity after the cluster
    multiplier (Stage 17's `stage17_K_density` logic, ported as
    `stage19_K_density`)

  K_density = 0.1096 was fixed by rule: the mean of the two soils' native
  effective densities at n=128 (Sand 0.018, Vertisol 0.201). The same value
  is used at both grids.
- **Shared and inherited from Stage 8:** `D_scale` 12.2582, `D_film`
  0.3831, film K scale 0.45, T = 1500 (75 time units), dt 0.05. Every
  per-soil Stage 9-16 key is stripped from the mapping, so no per-soil
  D_scale, T or dose can leak in. `stage11_implicit_tol` (a solver
  tolerance, <1% effect per Stage 11) is kept.
- **Model code:** the new switches in `dual_porosity.py` (`stage19_*`
  densities, `pore.coating_fraction`, wet/dry delivery series) are off by
  default. With them unset, output is byte-identical to `main` (Stage 16
  probe 9 at n=40 and the legacy `coating_trapped` path, `np.array_equal`).

## 2. Step 0 -- budget from the grids alone (`stage19_step0_budget.json`)

Regime predictions were committed before any simulation
(`results/stage19/stage19_step0_predictions.md`, commit c9f82c4).

| | Sand n=40 | Sand n=128 | Vertisol n=40 | Vertisol n=128 |
|---|---|---|---|---|
| macropore regions | 1 | 3 | 9 | 70 |
| habitat cells / N | 0.82 | 0.82 | 0.32 | 0.32 |
| porous-matrix cells / N | 0.12 | 0.12 | 0.44 | 0.44 |
| **dry fraction of habitat** | **1.000** | **0.999** | **0.893** | **0.913** |
| coating / trapped cells | 166 / 24 | 1,844 / 109 | 407 / 297 | 4,269 / 2,946 |
| OM per coating cell (x OM_density) | 7.7 | 7.1 | 1.18 | 1.15 |
| OM per trapped cell (x OM_density) | 13.3 | 30.1 | 3.8 | 3.9 |
| median geodesic distance trapped OM -> habitat | 4 cells (40 um) | 4 (40 um) | 6 (60 um) | 5 (50 um; 22 cells unreachable) |
| trapped diffusion time L^2/D_matrix | 44 | 44 | 98 | 68 |

| OM_density | OM per unit habitat capacity (both soils, by R5) | coating supply / maintenance, Sand | Vertisol |
|---|---|---|---|
| 0.037 | 0.34 | 5.9 | 2.2 |
| 0.084 | 0.77 | 13.5 | 5.1 |
| 0.192 | 1.76 | 30.8 | 11.5 |
| 0.439 | 4.00 | 70.2 | 26.3 |
| 1.000 | 9.13 | 160 | 60 |

Supply/maintenance = k_dis x coating OM / (m0 x B0) at t=0.

Timescales:
- max growth 1/r_max = 1
- coating release 1/k_dis = 20
- matrix-to-macro leak 1/k_leak = 10
- starvation 1/(m0 + m_s) = 8.3
- window 75

Growth threshold S* = 0.0245 (f* = 0.109). Film diffusivity in dry habitat
is 4.70; matrix diffusivity is 0.37.

**Pre-registered predictions and how they fared:**

| prediction (before Step 1) | outcome |
|---|---|
| x = 0.037 starves, both soils | **correct**: both dead |
| x >= 0.084 alive, both soils | **wrong**: Sand alive only from 0.192, Vertisol only from 0.439 |
| x = 0.44 and 1.0 converge toward a saturated bloom | **correct direction**: distance 0.223 then 0.101 |
| supply-limited window ~0.084-0.19 is where distinct, alive curves live | **half right**: curves are most distinct there, but Vertisol is not alive there |

The threshold estimate was too optimistic because it ignored the transport
step. Coating OM dissolves in the porous matrix, not in habitat, and must
leak or diffuse into habitat (~10+ time units) while the seeded biomass
starves (8.3). That lag costs Vertisol most: it has the lower coating
fraction and its coating OM is spread thinly (1.15x mean density).

## 3. Step 1 -- regime map, n=40, one seed (`stage19_step1_map.csv`, `stage19_step1_map.png`)

Liveness (pre-registered, both soils required): B_peak/B0 > 1 AND
R_peak >= 2 m0 B0 AND growth term > 50% of R at the peak.

| OM_density | soil | alive | B_peak/B0 | R_peak/(m0 B0) | growth share at peak | peak @ | fuel left | distance |
|---|---|---|---|---|---|---|---|---|
| 0.037 | Sand | no | 1.000 | 1.00 | 0.00 | 0.00 | 0.023 | 0.140 |
| 0.037 | Vertisol | no | 1.000 | 1.00 | 0.00 | 0.00 | 0.023 | |
| 0.084 | Sand | no | 1.000 | 1.87 | 0.78 | 0.35 | 0.023 | **0.733** |
| 0.084 | Vertisol | no | 1.000 | 1.00 | 0.00 | 0.00 | 0.023 | |
| 0.192 | Sand | **yes** | 1.250 | 5.63 | 0.80 | 0.28 | 0.023 | 0.365 |
| 0.192 | Vertisol | no | 1.000 | 2.69 | 0.78 | 0.41 | 0.023 | |
| **0.439** | Sand | **yes** | 3.017 | 12.80 | 0.80 | 0.24 | 0.023 | **0.223** |
| **0.439** | Vertisol | **yes** | 1.497 | 6.62 | 0.79 | 0.31 | 0.023 | |
| 1.000 | Sand | yes | 6.315 | 21.43 | 0.82 | 0.18 | 0.023 | 0.101 |
| 1.000 | Vertisol | yes | 3.482 | 14.07 | 0.79 | 0.26 | 0.023 | |

- **Alive at both soils: 0.439 and 1.0. Selected 0.439** (higher distance
  of the two).
- **The closest miss is OM_density 0.192**, distance 0.365. Vertisol there
  has a real, growth-dominated respiration rise (2.69x floor, 78% growth),
  but its total biomass never exceeds its t=0 value (B_peak/B0 = 1.000:
  biomass first decays, then recovers to just under B0). By the
  pre-registered criterion it is not alive, and the criterion was not
  relaxed.
- **`fuel_remaining_final` is 0.023 at every point.** With one uniform
  dissolution rate, every OM pool dissolves as exp(-k_dis t); after 75 time
  units 2.3% is left regardless of placement or dose. Placement changes
  *where* substrate appears, not *how much* or *when* it dissolves.

## 4. Step 2 -- confirmation at n=128 (`stage19_step2_n128.json`)

At OM_density 0.439:

| soil | alive | B_peak/B0 | R_peak/(m0 B0) | growth share | peak @ | R_end/R_peak | shape |
|---|---|---|---|---|---|---|---|
| Sand | yes | 2.892 | 11.97 | 0.80 | 0.24 | 0.17 | falling |
| Vertisol | yes | 1.475 | 6.38 | 0.79 | 0.31 | 0.46 | flat |

Distance 0.238 (n=40: 0.223). The same regime at both grids, as R5
intended.

## 5. Step 3 -- distinctness at n=128 and controls (`stage19_step3_controls.json`, `stage19_step3_controls.png`)

**Distinctness: 0.238 < 0.3 -> FAIL** (both soils alive).

| control | distance | both alive | notes |
|---|---|---|---|
| baseline (cf 0.8 / 0.3) | 0.238 | yes | Sand peaks at 0.24 and falls to 0.17 of peak; Vertisol peaks at 0.31 and holds 0.46 |
| **(a) equal coating_fraction 0.5** | **0.070** | yes | Sand peak 0.28, Vertisol peak **0.245 (now earlier than Sand)**; R changes up to 25% (Sand) / 35% (Vertisol) |
| **(b) transport off** (D_scale 1e-9) | 0.496 | **no** | B_peak/B0 = 1.000 both; R changes 83% (Sand) / 79% (Vertisol) |

- **(a) The per-soil coating fraction carries most of the distinctness, and
  its direction.** At equal coating fractions the two real pore structures
  differ by 0.070, and Vertisol, not Sand, is the faster soil. Structure
  alone supplies about a quarter of the baseline distance; the hand-set
  per-soil parameter supplies the rest, and also sets which soil is the
  "fast" one.
- **(b) Transport is essential for life here, unlike Stage 18.** With
  diffusion off, neither soil's biomass ever exceeds B0 (Sand's R still
  reaches 3.7x the floor from coating cells adjacent to habitat; Vertisol
  1.7x). The larger distance (0.496) is two soils failing in different
  ways, the same artifact that inflated Stage 8.
- **(c) Delivery runs almost entirely through dry water films.**

  | soil | wet fraction of habitat | delivered to habitat via wet pores | via dry film |
  |---|---|---|---|
  | Sand | 0.06% | 5.1 (0.07%) | 6,762 (99.9%) |
  | Vertisol | 8.7% | 918 (20%) | 3,580 (80%) |

  Wet fraction is constant in time (wicking off, static matric map).
  Delivery totals are summed over the window. So "transport" in this model
  is mostly `D_film` x `D_scale`, one shared film diffusivity acting on
  dry habitat. That is a physically thin pathway carrying most of the
  load, and it is the same for both soils.

## 6. Honest status

- **Achieved:** both soils alive, with equal OM, B0 and habitat-capacity
  densities, no per-soil dose, D_scale or T, one dissolution rate, and the
  same regime at n=40 and n=128. This fixes the two defects Stage 8 had
  (dead soils, grid dependence).
- **Not achieved:** distinctness. At the only density where both soils
  live and that the procedure selects, the distance is 0.238. What there
  is comes mostly from the per-soil coating fraction (0.238 -> 0.070 at
  equal fractions). The number that "passed" in Stage 8 (0.7955) is
  reproduced in kind at OM_density 0.084 (0.733), but only because one
  soil is dead there.
- **Not run:** Step 4 (seeds, shuffle control, mass conservation), per the
  stop-early rule.
- **No rescue attempted:** OM_density was the only quantity varied, over
  the five pre-registered values. No other parameter was changed.

## 7. Deliverables

- `dual_porosity.py`: `stage19_om_density`, `stage19_B0_density`,
  `stage19_K_density` (on K_t), `pore.coating_fraction`, and observational
  `delivered_wet_t` / `delivered_dry_t` / `wet_habitat_frac_t` /
  `K_effective_sum`. All off or read-only by default; byte-identical to
  `main` when unset.
- `configs/soil_sand_stage19.yaml`, `configs/soil_vertisol_stage19.yaml`.
- `stage19_run.py` (`step0` ... `step3`; `step4` refuses to run unless
  Step 3 passed).
- `results/stage19/`: step-0 budget and committed predictions, step-1 map
  (CSV / JSON / PNG), step-2 n=128 JSON, step-3 controls JSON / PNG.
