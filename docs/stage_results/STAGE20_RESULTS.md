# Stage 20 RESULTS -- film transport x wetting rule, 2x2 factorial (Sand, Vertisol)

`stage20_run.py`, `results/stage20/`. Branched from Stage 19
(`claude/stage19-alive-stage8`). **Two soils only; Loess was not run in any
form.**

## 0. Verdict

**Arm A1 (physical Millington-Quirk film + matric wetting) passes Steps 1-3:
both soils alive and distinct at n=128 (distance 0.689). At Step 4 it passes
seed robustness and mass conservation, but FAILS the pre-registered shuffle
criterion at n=128.** It is not a full pass. A3 fails Step 3.

| | A0 = Stage 19 (F0+W0) | **A1 (F1+W0)** | A3 (F1+W1) |
|---|---|---|---|
| n=128 distance (OM 0.439) | 0.238 | **0.689** | 0.170 |
| both soils alive (primary) | yes | **yes** | yes |
| (a) equal coating_fraction 0.5 | 0.070 | **0.639** | 0.051 |
| (b) transport off | 0.496, both dead | 0.496, both dead | 0.161, Vertisol dead |
| Step 3 (alive AND > 0.3) | fail | **PASS** | fail |
| S3a 3 seed sets, n=128 / n=40 | -- | **pass / pass** (min 0.607 / 0.870) | -- |
| S3b seed noise < 0.5 x between | -- | **pass** (0.014 vs 0.607) | -- |
| S4 shuffle lowers distance by >= 0.1 | -- | **FAIL** at n=128 (0.689 -> 1.901); pass at n=40 (drop 0.299) | -- |
| mass conservation (< 1e-9) | -- | **pass** (5.1e-12 / 9.8e-13) | -- |

Plainly:

- **Making the film physical is what separates the soils.** Switching F0
  -> F1 at matric wetting raises distance by +0.26 to +0.72 at every
  density, and the separation survives equal coating fractions (0.639 at
  n=128). In Stage 19, equal coating fractions cut the distance to 0.070.
  For the first time in this series, the distinctness is carried by
  structure acting through transport, not by the per-soil coating
  fraction.
- **The shuffle control fails as pre-registered at n=128, for a reason
  worth stating.** Shuffling the pore field does not make the soils
  converge. It kills Vertisol: its macropores break into 2,088 tiny regions
  (median inscribed diameter 20 um, below the 30 um habitat threshold).
  Distance rises (0.689 -> 1.901) because one soil is dead, the same artifact
  as the transport-off control. Arrangement clearly matters, but the S4
  metric as I pre-registered it cannot show that when the shuffle destroys
  the habitat itself. It is recorded as FAIL, not reinterpreted.
- **Equal saturation (W1) does not help.** It removes Sand's drained-film
  problem but drains all of Vertisol's habitat (its fine matrix takes all
  the water). It lowers distance at every density in both film settings.

## 1. The two factors (as specified; nothing else changed from Stage 19)

- **F0:** constant drained-pore D = D_film x D_scale = 0.3831 x 12.2582 =
  4.70 (Stage 19).
- **F1:** `stage20_film_mode: millington_quirk`:
  D_film(i,j) = D(d_ij) x D_scale x S_w^(10/3).
  - D(d) is the existing `D_of_d` (same constants), evaluated on the same
    diameter the model already uses for macropore D (the region-inscribed
    `d_struct`), with no cluster bonus.
  - S_w = wet non-solid cells / non-solid cells, computed once per soil,
    grid and wetting rule.
  - Applying a bulk-soil relation (Millington & Quirk 1961) to individual
    drained cells is an approximation. Its role here is to guarantee
    f_rel <= 1, so a film is never faster than a full pore of the same size.
- **W0:** `matric_regions`, d_cut 50 um (Stage 19).
- **W1:** `saturation.mode: equal_saturation`, S_w = 0.5. Non-solid cells
  fill in ascending order of their own per-cell PSD diameter (stable sort,
  ties by cell index) until round(0.5 x n_nonsolid) are wet; solid is never
  wet. This is a new path, not `global_theta`.
- **Unchanged from Stage 19:**
  - literature PSDs, biology (SHA-256 asserted), `k_dis_slow = k_dis`,
    `om_fine_bias = 0`, no enzyme
  - `coating_trapped`, thickness 2 (coating needs a wet matrix path, so it
    follows the active wetting rule), coating_fraction 0.8 / 0.3
  - densities OM / B0 0.0125 / K_t 0.1096, asserted every run
  - D_scale 12.2582, k_leak, T 1500, dt 0.05, seeds
- **Byte-identity:** with the switches unset, Stage 19's n=128 curves at OM
  0.439 reproduce exactly (`np.array_equal`, both soils). A0 is Stage 19's
  results, reused.

## 2. Step 0 -- before any run (`stage20_step0_budget.json`; predictions committed in 2973dde, `stage20_step0_predictions.md`)

| soil, n | wetting | S_w | f_rel | wet habitat | drained D F0 | drained D F1 (median) | drained / full-pore D, F0 -> F1 | coating / trapped cells | native K_t density |
|---|---|---|---|---|---|---|---|---|---|
| Sand 40 | W0 | 0.127 | 0.0010 | 0.0% | 4.70 | 0.009 | -- (no wet habitat) | 166 / 24 | 0.286 |
| Sand 40 | W1 | 0.500 | 0.099 | 42.7% | 4.70 | 0.89 | 0.30 -> 0.057 | 166 / 24 | 0.436 |
| Sand 128 | W0 | 0.128 | 0.0011 | 0.06% | 4.70 | 0.013 | 10.5 -> 0.029 | 1,844 / 109 | 0.018 |
| Sand 128 | W1 | 0.500 | 0.099 | 42.7% | 4.70 | 1.23 | 0.22 -> 0.057 | 1,844 / 109 | 0.028 |
| Vertisol 40 | W0 | 0.622 | 0.206 | 10.7% | 4.70 | 0.84 | 4.9 -> 0.87 | 407 / 297 | 0.227 |
| Vertisol 40 | W1 | 0.500 | 0.099 | **0%** | 4.70 | 0.41 | -- (no wet habitat) | 262 / 442 | 0.208 |
| Vertisol 128 | W0 | 0.614 | 0.197 | 8.7% | 4.70 | 0.63 | 4.9 -> 0.67 | 4,269 / 2,946 | 0.201 |
| Vertisol 128 | W1 | 0.500 | 0.099 | **0%** | 4.70 | 0.32 | -- (no wet habitat) | 2,688 / 4,527 | 0.187 |

- **Under W1, Vertisol's habitat is fully drained.** Its matrix takes all
  the water (87% of matrix cells wet), and its coating shrinks because 13%
  of the matrix is dry. Reported, not fixed.
- **K_t's native density shifts under W1** (Sand 0.018 -> 0.028 at n=128);
  0.1096 is kept as specified.

**Predictions vs outcome:**

| arm | predicted | outcome |
|---|---|---|
| A1 | Sand's in-habitat transport collapses (film 0.013), leak-fed only, liveness marginal; few or no both-alive densities; any distance from Sand starving | **partly wrong.** Sand is leak-fed (correct) but alive at 0.439 and 1.0, and both-alive distances are the highest in the series (0.943 / 0.819, n=40). The distance is structural (equal-cf control 0.924 / 0.930), not starvation. |
| A2 | Sand alive at lower density; distance similar to or above A0 | Sand alive from 0.192 (same as A0); **distance lower** than A0 at every density |
| A3 | Vertisol dead except at 1.0; no both-alive > 0.3 | Vertisol alive (primary) at 0.439 and 1.0; no both-alive > 0.3 -- **correct on distance, wrong on Vertisol's death** |
| overall | physical films make the trade-off worse for Vertisol | wrong for A1: physical films separate the soils, and Vertisol stays alive |

## 3. Step 1 -- regime map, n=40, one seed (`stage20_step1_map.csv/.json`, `stage20_step1_map.png`, `stage20_step1_log.txt`)

Distance per arm and the factorial effects:

| OM | A0 | A1 | A2 | A3 | F: A0->A1 | F: A2->A3 | W: A0->A2 | W: A1->A3 |
|---|---|---|---|---|---|---|---|---|
| 0.037 | 0.140 | 0.768 | 0.093 | 0.148 | +0.628 | +0.055 | -0.047 | -0.619 |
| 0.084 | 0.733 | 0.995 | 0.204 | 0.210 | +0.262 | +0.006 | -0.529 | -0.785 |
| 0.192 | 0.365 | 0.771 | 0.151 | 0.194 | +0.406 | +0.042 | -0.214 | -0.578 |
| **0.439** | 0.223 | **0.943** | 0.135 | 0.192 | +0.720 | +0.057 | -0.088 | -0.751 |
| 1.000 | 0.101 | 0.819 | 0.080 | 0.058 | +0.718 | -0.022 | -0.022 | -0.761 |

Liveness per soil, Sand / Vertisol. P = passes PRIMARY (verdict); s =
passes SECONDARY only (respiration-only: R_peak >= 2 m0 B0 and growth >
50% at peak); - = neither.

| OM | A0 | A1 | A2 | A3 |
|---|---|---|---|---|
| 0.037 | - / - | - / - | - / - | - / - |
| 0.084 | - / - | - / - | - / - | - / - |
| 0.192 | P / s | s / s | P / s | P / s |
| 0.439 | P / P | P / P | P / P | P / P |
| 1.000 | P / P | P / P | P / P | P / P |

- **Both soils PRIMARY-alive only at 0.439 and 1.0 in every arm.** Under
  the SECONDARY (respiration-only) definition, 0.192 would also count in
  all four arms. Reported for the user's decision; it was not used.
- **F effect:** large and positive under matric wetting (+0.26 to +0.72),
  near zero under W1 (Vertisol's habitat is all drained there and Sand's
  film is already moderate).
- **W effect:** negative everywhere (-0.02 to -0.79).
- **Distances where a soil is dead** (A1 at 0.037-0.192: 0.768-0.995) are
  shown for completeness only; they measure two soils failing differently.

**Control (a) (equal coating_fraction 0.5) at every candidate density with
distance > 0.3** (only A1 qualifies; A3 never exceeds 0.21):

| A1, OM | baseline | equal cf 0.5 | both alive (equal cf) |
|---|---|---|---|
| 0.037 | 0.768 | 0.710 | no |
| 0.084 | 0.995 | 0.530 | no |
| 0.192 | 0.771 | 0.483 | no |
| **0.439** | **0.943** | **0.924** | **yes** |
| 1.000 | 0.819 | 0.930 | yes |

At the both-alive densities the distance survives equal coating fractions
(0.924 / 0.930). **Structure, not the per-soil coating fraction, carries
it.**

**Where substrate enters habitat** (A1 / A3 at 0.439, n=40, summed over the
window). These are **net, signed** totals over each habitat set:
- "diffusion" = net transport by the D field across the set's boundary
  (negative = net export out of the set)
- "leak" = matrix->macro exchange (one-directional, so >= 0)

Diffusion + leak = net transport gain exactly.

| arm, soil | diffusion -> wet habitat | diffusion -> drained habitat | leak -> wet | leak -> drained | net -> drained |
|---|---|---|---|---|---|
| A1 Sand | 0 | **-215** | 0 | **740** | 525 |
| A1 Vertisol | 57 | 234 | 23 | 112 | 346 |
| A3 Sand | 357 | 171 | 137 | 1 | 172 |
| A3 Vertisol | 0 | 172 | 0 | 199 | 371 |

**In A1, Sand gets no net film delivery at all.** Leak deposits substrate
at matrix-interface habitat cells, and diffusion carries part of it back out
of the drained habitat into the wet matrix (net -215). Sand in A1 is fed
entirely by leak.

**Correction to Stage 19's wording.** Stage 19's "delivered via film"
counter summed everything entering drained habitat, including matrix leak.
With the signed split (added here, read-only), Stage 19's Sand at n=128
receives 6,761 net into drained habitat: 4,930 (**73%**) by diffusion and
1,831 (27%) by leak; Vertisol 75% / 25%. Stage 19's conclusion (mostly
film) holds, but "99.9% via film" overstated it.

*(An earlier version of this split clipped the diffusion and leak gains
separately per cell, so they did not add up to the net. It was fixed
after review; all numbers here are the signed version.)*

**Selection** (both PRIMARY-alive, highest distance, per candidate arm):
**A1 at 0.439** (0.943) and **A3 at 0.439** (0.192).

## 4. Steps 2-3 -- n=128 (`stage20_step2_n128.json`, `stage20_step3_controls.json`, `stage20_n128_curves.png`)

| | Stage 19 (A0) | A1 | A3 |
|---|---|---|---|
| distance | 0.238 | **0.689** | 0.170 |
| Sand: B_peak/B0, R_peak/floor, peak @, R_end/R_peak | 2.89, 11.97, 0.24, 0.17 | 1.34, 4.01, **0.12**, 0.85 | 2.97, 12.52, 0.21, 0.13 |
| Vertisol: B_peak/B0, R_peak/floor, peak @, R_end/R_peak | 1.48, 6.38, 0.31, 0.46 | 1.39, 5.71, 0.32, 0.54 | 1.20, 5.02, 0.24, 0.48 |
| (a) equal cf 0.5 | 0.070 | **0.639** (both alive) | 0.051 |
| (b) transport off | 0.496 (both dead) | 0.496 (both dead) | 0.161 (Vertisol dead) |
| (c) Sand: drained fraction; net diffusion / leak into drained habitat | 0.999; 4,930 / 1,831 | 0.999; **-2,424 / 7,805** | 0.573; 1,880 / 3 (wet: 3,383 / 1,546) |
| (c) Vertisol: drained fraction; net diffusion / leak into drained habitat | 0.913; 2,683 / 897 | 0.913; 2,094 / 1,240 (wet: 449 / 486) | 1.000; 1,483 / 2,130 |

**A1's curves** (`stage20_n128_curves.png`):
- **Sand** bursts early (peak at 0.12 of the window, t ~ 9), dips, then
  holds a sustained plateau (R_end/R_peak 0.85). Its habitat interior
  barely conducts (film 0.013), so after the first burst it is fed entirely
  by leak at matrix interfaces (net diffusion into its drained habitat is
  negative, -2,424).
- **Vertisol** rises smoothly to a later peak (0.32) and declines slowly.

The real-soil shape targets of earlier stages (Sand burst-and-crash,
Vertisol sustained rise) are not part of this stage's criteria. Measured
against them, A1 gets Sand's early burst but not its crash, and Vertisol's
slow rise but not its sustain.

## 5. Step 4 -- robustness of A1 (`stage20_step4_robustness.json`, `stage20_step4_log.txt`)

Criteria pre-registered in `stage20_run.py` before the run:
- **S3a** all 3 seed sets both alive and > 0.3
- **S3b** seed noise < 0.5 x min between-soil distance
- **S4** shuffle lowers distance by >= 0.1
- **MASS** transport step relative error < 1e-9

| | n=128 | n=40 |
|---|---|---|
| distance, seed sets +0 / +10 / +20 | 0.689 / 0.677 / 0.607 | 0.943 / 0.959 / 0.870 |
| both alive in all three | yes | yes |
| **S3a** | **pass** | **pass** |
| seed noise vs min between | 0.014 vs 0.607 | 0.064 vs 0.870 |
| **S3b** | **pass** | **pass** |
| shuffled distance (both alive?) | 1.901 (**Vertisol dead**) | 0.644 (**Vertisol dead**) |
| **S4** (drop >= 0.1) | **FAIL** (drop -1.212) | pass (drop 0.299) |
| **MASS** | **pass** (Sand 5.1e-12, Vertisol 9.8e-13) | -- |

**Why the shuffle kills Vertisol:**

| | macropore regions | median region inscribed diameter | macropore cells with K > 0 | wet habitat |
|---|---|---|---|---|
| Vertisol, original | 70 | 100 um | 99% | 8.7% |
| Vertisol, shuffled | 2,088 | **20 um** | **10%** | 100% |
| Sand, original | 3 | 340 um | 100% | 0.1% |
| Sand, shuffled | 24 | 100 um | 100% | 0.2% |

Habitat size comes from each connected region's inscribed diameter, so
scattering Vertisol's pores destroys 90% of its habitat. The shuffle
control therefore does not hold "the same pore sizes" constant in the sense
that matters here. It changes the amount of habitat, not just its
arrangement. The outcome shows arrangement matters a great deal, but the
pre-registered metric, which expected convergence, fails at n=128 and
passes at n=40 only because both distances are dominated by Vertisol's
death.

**One documented convention in the shuffle runs only.** In shuffled Sand,
every porous-matrix cell lies within the coating distance of a macropore,
so the trapped pool is empty, and `split_om_coating_trapped` would silently
drop its 20% of OM (placed/target 0.80, caught by the R5 assertion). With
`stage20_conserve_om_pools` (set only for shuffle runs, off by default), an
empty pool's share moves to the other pool. Shuffled Sand therefore has an
effective coating fraction of 1.0. No rate changes; unshuffled runs never
have an empty pool.

## 6. Honest status

- **Passes:** A1 is alive and distinct at n=128 (0.689); robust across
  three seed sets at both domain sizes; mass-conserving; and its
  distinctness survives equal coating fractions. This means structure,
  through physically bounded film transport, carries it. That is the first
  time in this series.
- **Fails:** the pre-registered shuffle criterion at n=128. The shuffle
  destroys Vertisol's habitat, so the control does not isolate
  arrangement.
- **Does not help:** the equal-saturation rule (W1). It lowers distance
  everywhere and drains Vertisol's habitat completely.
- **Caveats:**
  - Sand in A1 is a leak-fed soil. Almost all of its drained habitat
    interior is transport-dead at S_w = 0.13, and its life comes from the
    matrix interfaces. Whether that is physically right depends on the
    matric wetting rule draining Sand's whole connected region (the second
    defect named in the brief), which W1 was meant to test and which did
    not work out.
  - n=40 vs n=128 are two domain sizes at a 10 um voxel, not two
    resolutions.
  - The SECONDARY liveness definition would admit OM 0.192 in all arms;
    the user decides whether it should replace PRIMARY.

## 7. Deliverables

- `dual_porosity.py`, all off by default, byte-identical to Stage 19 when
  unset:
  - `stage20_film_mode: millington_quirk` (F1)
  - `saturation.mode: equal_saturation` (W1)
  - `S_w` and `film_f_rel` outputs
  - `pore.shuffle_seed`
  - `stage20_conserve_om_pools` (shuffle runs only)
  - read-only, signed diffusion/leak delivery split (sums exactly to the
    net)
- `stage19_run.py`: reports the delivery split (no behaviour change).
- `stage20_run.py`: `step0` ... `step4`.
- `results/stage20/`: step-0 budget and committed predictions; step-1 map
  (CSV / JSON / PNG / log); step-2 and step-3 JSON; `stage20_n128_curves.png`;
  step-4 JSON / log.
