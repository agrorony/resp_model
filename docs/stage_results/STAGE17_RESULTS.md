# Stage 17 RESULTS -- reconciling v3's equal-totals machinery with main's real-soil battery

`stage17_run.py` (runtime 18.3 min, 258 simulations), `results/stage17/`.
A SYNTHESIS stage: the v3 branch (`claude/brave-brahmagupta-edjq14`, branched
from the first commit, before Stages 1-16) is not merged; its rules, controls
and enzyme mechanism are ported as switches/tooling and run against main's
data, physics and SS4 battery, which remain the source of truth.

## 0. Verdict up front

**Adopted configuration: NONE. No configuration that satisfies v3's
equal-totals rules passes main's SS4 battery, at either grid. SS4(b),
recruitment timing, fails in every one of the 80 configurations tested
(72 R2 probes + 8 diagnostics), enzyme on or off.**

| SS4 sub-criterion | best R2 result, n=128 | best R2 result, n=40 |
|---|---|---|
| (a) Sand burst-crash (6 checks) | PASS in several arms -- but never together with (a) Vertisol | PASS in several arms |
| (a) Vertisol slow rise (5 checks) | PASS (`U_on`, robust on 3/3 seeds; `F_on` seed 1 only) | PASS (`F_off`, `F_on`, `U_on`) |
| (a) BOTH soils in one configuration | **FAIL** -- best is 11/11 minus one check (Sand peak at 0.45 of window vs < 1/3), in `F_on` and `U_on` | **PASS** in `F_off` (seed 1 only) and `F_on` (seed 1 only); fails on the other two seed sets |
| (b) recruitment timing (frac@2% < 0.90 AND t_final > 0.10) | **FAIL everywhere** | **FAIL everywhere** |
| (c) hotspot selectivity (reported) | Vertisol 21.4% (`U_on`) vs Portell 9.5 +/- 4.0% | 33-78% |
| (d) distinctness > 0.3 | PASS, robust (`U_on` 1.40 / 1.45 / 1.42 over 3 seeds) | PASS (0.42-1.78) |

And on the rule side, stated explicitly as asked:

- **`om_mode: shared_fine_pore` cannot satisfy R3 at any dose.** The
  hoarding ratio (max cell OM / mean porous-matrix OM) is a property of the
  `d^-2.5` rule and the pore field, not of `om_total`: 54-288 (Sand-Vertisol,
  n=128), 30-165 (n=40), worst case 746 over 3 soils x 2 grids x 3 seeds,
  against R3's cap of 30. The top 1% of cells hold 60-73% of all OM. Stage
  16's own best result (2.082/2.152) is therefore R3-noncompliant, and also
  R2-noncompliant by design (per-soil `om_total` 3,000 vs 50,000, `D_scale`
  80 vs 0.7, `T` 3000 vs 1500, and a 15x difference in total habitat
  sum(K): 173 vs 2,570 at n=128).
- **Neither side was relaxed.** Arms `F_*` keep main's rule (R3 fails,
  reported); arms `U_*` use the largest shared exponent that satisfies R3
  (bias 0.75, worst-case ratio 21.1), derived from R3 alone before any SS4
  run. Both are judged by the same unchanged SS4.

## 1. Premise corrections (verified, user-confirmed before running)

The task described main as "real measured PSDs, `global_theta`, three
soils". Main's operative Stage 16 configuration is different, and the
difference is load-bearing:

| | task premise | main's Stage 16 (`stage16_run.py`) |
|---|---|---|
| PSD | real measured (`psd_measured`) | **literature** (`psd_parametric`, `configs/soil_<s>.yaml`) |
| saturation | `global_theta` | **`matric_regions`** (d_cut 50 um; theta emerges) |
| soils under SS4 | Sand, Vertisol, Loess | **Sand, Vertisol** (Loess dropped at Stage 10, reintegrated only in Stage 15 on the retired coating path) |

The reason main does not use the measured PSDs is Stage 13's finding: the
measured Sand and Vertisol tables carry 0% volume below the 30 um
macro/matrix split, so the porous matrix is empty and **no OM can be
placed**. With the user's agreement both tracks were run:

- **Track R (as asked): measured PSDs + `global_theta` (theta 0.6), three
  soils, R2 midpoint** (om 12,247, D 7.5, T 3000, equal K). Reconfirmed
  exactly: Sand and Vertisol `om_total_placed = 0`, `R_peak = 0.400`
  (the `m0*B0` floor), dead, at both grids; SS4(a) 3/6 and 0/5 (Sand's 3
  "passes" are vacuous: no growth, no fuel). Loess places its full 12,247
  and bursts (`R_peak` 86.3 / 30.3) but has no SS4 criterion of its own.
  **On real measured data the SS4 battery cannot be passed by any shared OM
  rule that stores carbon in the sub-30 um matrix -- independent of dose,
  transport, enzymes or equal totals.**
- **Track L (primary): main's literature-PSD / `matric_regions` setup**,
  Sand + Vertisol under SS4, Loess run alongside and reported descriptively.

## 2. What was ported from v3, and how (SS6 of the task)

All switches are OFF by default; with them unset `dual_porosity.py` is
**byte-identical** to main (`np.array_equal` on Stage 16 probe 9 at n=40 --
Sand `R_peak` 32.541, Vertisol 4.415 -- and on the `coating_trapped` legacy
path).

| v3 item | where it lives on this branch |
|---|---|
| R1 biology frozen | `stage17_run.BIOLOGY_SHA256`, asserted at start (`configs/biology.yaml` untouched) |
| R2 equal totals | shared `om_total`, `D_scale`, `T` (=3000) via Stage 16's rung-1 keys; **equal habitat** via new `stage17_K_density`, applied to the *effective* capacity `K_t` (after the soil-specific cluster multiplier -- fixed after review, f1e3b68); shared density = mean native sum(K_t)/N (0.139 at n=128, 0.247 at n=40). Verified equal in all 74 R2 configurations. |
| R3 OM hoarding cap | `grids["om_hoard"]` in `build_dual_grids`; cap 30 x mean porous-matrix OM |
| S3 seed robustness | 3 seed sets (+0/+10/+20 on each soil's main seed) |
| S4 shuffle control | `pore.shuffle_seed` permutes the diameter field (same PSD multiset, arrangement destroyed) |
| S5 mass conservation | measured for main's implicit Jacobi solver: relative error 1.8e-16 (fixed iterations), 1.9e-11 (tol 1e-7 early stop) |
| P3 enzyme depolymerization | `mapping.enzyme` in `simulate_dual_porosity`: E produced by B (a_E), decays (d_E), diffuses at f_E x D_total through the same water-gated field; every OM pool releases at `rate * E/(E+K_E)` |
| n-D transport, fabric archetypes, robustness harness | `structure_tools.py` -- numpy-only port (main's Stage 2 rule), not wired into `dual_porosity` (2D implicit solver/BFS); conservation verified (sum div ~1e-14 in 2D and 3D) |

**Enzyme constants.** v3's a_E=1, d_E=0.1, f_E=0.1. K_E is ported as v3's
*dimensionless* ratio kappa = K_E/E* = 0.4, E* = a_E x (B0_total/N)/d_E,
because main's per-cell biomass at n=128 is 10x below v3's and the literal
v3 K_E=0.05 is mis-scaled there. Both were run (Part D): the literal value
kills both soils at n=128 (`R_peak` 0.40, distinctness 0.056).

## 3. The R2 ladder (Part C): 4 arms x 9 probes x 2 grids

Probes: `om_total` in {3,000, 12,247, 50,000} x `D_scale` in {0.7, 7.5, 80}
(Stage 16's own anchors), T=3000 shared, K equalized. Arms: OM rule (main's
bias 2.5 / R3-compliant 0.75) x enzyme (off/on). Selection fixed in
advance: most SS4 criteria at n=128, then at n=40, then distinctness.
Full table: `stage17_ss4_scores.csv`.

| arm | selected probe | n=128: a_S / a_V / b / d | n=40: a_S / a_V / b / d | R3 |
|---|---|---|---|---|
| `F_off` main rule, no enzyme | om 3,000, D 7.5 | P / **F** (peak 0.36, dips to 0.44) / F / P (0.874) | P / P / F / P (1.166) | **F** (288) |
| `F_on` main rule + enzyme | om 3,000, D 80 | **F** (peak 0.45) / P / F / P (1.448) | P / P / F / P (1.778) | **F** (288) |
| `U_off` R3 rule, no enzyme | om 3,000, D 0.7 | P / **F** (B_peak/B0 = 1.00, no growth) / F / P (1.656) | F / F / F / P (0.352) | P (13) |
| `U_on` R3 rule + enzyme | om 12,247, D 7.5 | **F** (peak 0.45) / P / F / P (1.403) | **F** (no crash, 0.89) / P / F / P (0.421) | P (13) |

Rung-1 behavior from Stage 16 is reproduced under equal habitat: no single
shared dose/transport gives both soils their shapes at n=128. The
enzyme arms come closest (10 of 11 shape checks), and in both the one miss
is the same: **the enzyme bootstrap delays Sand's burst** (enzyme has to
accumulate before any OM dissolves, visible as a flat start to ~0.2 of the
window in `stage17_selected_curves.png`), pushing its peak to 0.45 against
the < 1/3 bar. The window T=3000 was fixed in advance and not moved to
rescue this.

## 4. Robustness of the selected probes (Parts E-F)

`stage17_summary.json` -> `robustness`.

| arm, grid | criteria passed per seed set | seed noise vs min between-soil | shuffle: distance / drop | Vertisol graded half-max (median, p25-p75) |
|---|---|---|---|---|
| `F_off` 128 | 2, 2, 2 | 0.103 vs 0.335 | 0.321 / 0.553 | 0.26 (0.17-0.40) |
| `F_off` 40 | 3, 2, 1 | 0.337 vs 0.273 (**not separable**) | 0.266 / 0.900 | 0.09 (0.06-0.18) |
| `F_on` 128 | 2, 1, 0 | 0.513 vs 0.210 (**not separable**) | 0.461 / 0.986 | 0.11 (0.07-0.68) |
| `F_on` 40 | 3, 2, 2 | 0.322 vs 0.642 | 0.296 / 1.482 | 0.12 (0.09-0.20) |
| `U_off` 128 | 2, 2, 2 | 0.013 vs 1.549 | 0.274 / 1.383 | 0.31 (0.26-0.37) |
| `U_off` 40 | 1, 1, 1 | 0.094 vs 0.352 | 0.171 / 0.181 | 0.10 (0.09-0.11) |
| **`U_on` 128** | **2, 2, 2** (Vertisol shape + distinctness on every seed) | **0.016 vs 1.403** | 0.353 / 1.050 | **0.72 (0.05-0.85)** |
| `U_on` 40 | 2, 2, 1 | 0.117 vs 0.421 | 0.772 / **-0.351** | 0.27 (0.24-0.29) |

- **The n=40 "both shapes pass" results are not seed-robust.** `F_off` and
  `F_on` pass (a) for both soils on seed set 1 only; on the other two
  Vertisol's shape fails. Reported, not adopted.
- **S4 holds where it is informative.** Shuffling the pore field (same PSD)
  destroys Vertisol's slow-rise shape in every arm (a_V False under shuffle)
  and removes 0.55-1.48 of distinctness. At n=128 the soils' differences are
  carried by arrangement, not only by the PSD histogram. The exception is
  `U_on` n=40, where shuffling *raises* distinctness (drop -0.351).
- **Enzyme S6 (+-50% on kappa and f_E)**: at `F_on` n=40 the both-shape
  result survives kappa x0.5, x1.5 and f_E x1.5 but not f_E x0.5; at n=128
  and in `U_on` no perturbation recovers Sand's early peak. (b) fails under
  every perturbation.
- **Three soils** (Loess descriptive): Vertisol-Loess is the weak pair in
  three of the four arms at n=128 (0.024-0.073, against 0.3), and only
  `F_on` separates it (0.910); at n=40 it is 0.12-0.49. Stage 9's "two
  large, many-region soils look alike" finding persists under equal
  totals.

## 5. Q3 answered: does enzyme-mediated depolymerization fix SS4(b)?

**No.** Not at Stage 16's own working point (Part D), and not under R2
(Part C, arms `*_on`):

- At rung 3 (Stage 16 probe 9, per-soil dose/transport) the enzyme
  **kills Vertisol**: `fuel_end = 1.000` (nothing ever dissolves),
  `R_peak` at the floor, both grids. Mechanism: all OM sits in the porous
  matrix and all biomass in macropores; enzyme enters the matrix with
  diffusivity 0.03 x D_scale x f_E = 0.002 at Vertisol's D_scale 0.7, so it
  cannot reach the carbon before the seeded biomass starves
  (m0 + m_s = 0.12 per time unit). Its boolean recruitment does stagger
  (frac@2% = 0.43 at n=128 and 0.875 at n=40 -- the only two values below
  0.90 in the whole stage) -- but only because the regions that recruit are
  dying, not growing. Neither is a pass of (b) (t_final 0.035 / 0.021 <
  0.10), and both cost (a).
- Under R2, in every enzyme run where Vertisol grows at all, the boolean
  count saturates immediately again: frac@2% >= 1.0 in all of them (values
  above 1 mean more regions were active at 2% of the window than at the
  end -- early activation, then decay) and t_final <= 0.015. Where it does
  not grow (e.g. om 3,000 / D 0.7) the count is empty and (b) is undefined,
  i.e. failed.

**Why the enzyme cannot fix (b) as defined.** Stage 16 traced (b)'s
failure to 66/70 Vertisol regions touching fuel, with activation needing
only `growth > 1e-9`. The enzyme changes *how much* OM turns over -- and
it does change the dynamics: the graded per-region half-max median moves
to **0.72 of the window** in `U_on` n=128 (vs 0.534 in Stage 16 and 0.440
in Stage 11). Most regions reach half their own peak growth late, the most
staggered recruitment recorded in this project. But every region is seeded
with biomass at t=0 (`B0 ~ K`), so every region produces enzyme from t=0,
and a few steps later every fuel-adjacent region crosses `growth > 1e-9`.
The boolean criterion cannot see the staggering; Stage 16 already flagged
this tension between (b) and (c). It is reported here again and not
resolved in this stage's favor: **(b) as written FAILS**.

## 6. Honest status

What passes, robustly, under v3's equal-totals rules on main's literature
soils:

- (d) distinctness at both grids, and seed-robust at n=128
  (`U_on`/`U_off`: seed noise 0.013-0.016 vs 1.40-1.55).
- (a) Vertisol's slow rise, seed-robust at n=128 under `U_on` (R3-compliant +
  enzyme), and arrangement-dependent (lost when the field is shuffled).

What fails:

- (a) for both soils at once at n=128, in every R2 configuration. The
  closest (`F_on`, `U_on`) miss one check (Sand peak 0.45 vs < 0.33),
  caused by the enzyme lag.
- (b) recruitment timing, everywhere, for the structural reason Stage 16
  named, which enzyme gating does not remove.
- R3 for main's `shared_fine_pore` rule, at any dose.
- Everything on Track R (real measured PSDs): Sand/Vertisol cannot hold OM.

Next step, not attempted here (consistent with Stage 16 SS6): a
dormant/active biomass split with substrate-dependent reactivation, so
that activation depends on how much substrate arrives rather than whether
any arrives. That is the degree of freedom the boolean (b) test needs, and
neither equal totals nor enzyme gating supplies it.

## 7. Deliverables

- `dual_porosity.py`: `stage17_K_density` (on `K_t`), `enzyme`,
  `pore.shuffle_seed`, `om_hoard`, `K_effective_sum`, `enzyme_mean_t` --
  all off by default, byte-identical to main when unset.
- `stage17_run.py`: the whole battery; rules, ladders, arms, enzyme
  constants and selection rule fixed in the file before running.
- `structure_tools.py`: v3's n-D transport, fabric archetypes
  (gaussian/anisotropic, Voronoi peds, biopores, hierarchical), PSD-quantile
  layout, shuffle/hoarding statistics and seed/shuffle harness, numpy-only.
- `configs/mapping.yaml`: Stage 17 switches documented (comments only;
  active configuration unchanged).
- `results/stage17/`: `stage17_all_runs.csv`, `stage17_ss4_scores.csv`,
  `stage17_robustness_runs.csv`, `stage17_summary.json`,
  `stage17_selected_curves.png`.
