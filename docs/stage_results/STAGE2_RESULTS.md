# Stage 2 Results — Water-Filled Pores & Percolation

This document summarizes **Stage 2** (`soil_respiration_prompt_stage2.md`):
switching the three soils from Stage 1's fully-wet transport to
`saturation.mode: global_theta`, adding percolation/connectivity
diagnostics, and sweeping `theta` to see whether the stronger, less-tunable
percolation lever gives genuinely different pore structures genuinely
different respiration where Stage 1's fully-wet transport-speed difference
could not (Stage 1's honest A–C negative, `MODEL_SPEC.md` §11).

## 1. The numpy-only conversion (§8, done first)

`pore_field.py`'s `scipy.ndimage.gaussian_filter` was replaced with a
hand-rolled separable-Gaussian convolution (`gaussian_filter_2d`): a 1D
Gaussian kernel (`radius = int(4*sigma + 0.5)`, matching scipy's
`truncate=4.0` default) applied along each axis. The one subtlety: scipy's
`mode="reflect"` duplicates the edge value (`d c b a | a b c d`), which is
`numpy.pad`'s `"symmetric"` mode — **not** numpy's own `"reflect"` mode
(which does not duplicate the edge). Using `numpy.pad(..., mode="reflect")`
first gave large discrepancies (soil C's smoothed field differed from
scipy's by up to 0.45 at `lambda=1.0`, cascading into an R_peak difference
of 2.58 → 4.15 through the nonlinear `aggregate` thresholding step).
Switching to `mode="symmetric"` reproduces scipy's output to machine
precision (`max abs diff ~1e-16`) across `lambda` = 1, 2, 4, 10, and the full
Stage-1 three-soil run is numerically identical to the scipy version
(`R_peak`/`cum_co2` match to the digits reported in `STAGE1_RESULTS.md`).
No other module imported scipy. The repo now runs on `numpy` + `pyyaml` +
`matplotlib` alone (`README.md` updated accordingly).

## 2. Mechanism (§3) and mass conservation

`model.water_mask_and_theta`'s `global_theta` mode fills the smallest pores
first: `r_cut` is *that soil's own* `theta`-quantile of `r`, so all three
soils share the same nominal saturation but reach it through different
percolation topologies. Air-filled cells get `D=0`, and harmonic-mean face
conductivities (`model.harmonic_mean`) make a dry cell block transport on
every face — fully isolating it from its neighbors, not merely slowing it.
Re-verified numerically with `theta=0.4` on soil A's actual pore field (960
of 1600 cells dry): `mass_after - mass_before = 0.0` exactly, confirming
disconnection doesn't leak or invent mass.

## 3. Theta sweep results

Swept `theta` ∈ {1.0, 0.8, 0.6, 0.4, 0.2} with each soil's **unchanged**
Stage-1 pore parameters (`mu, sigma, lambda, aggregate` — only the
saturation mode changed). Full numbers: `results/stage2_theta_sweep.csv`,
`results/stage2_distinctness_by_theta.csv`.

| theta | AB | AC | BC | min | success | A nontrivial | B nontrivial | C nontrivial |
|---|---|---|---|---|---|---|---|---|
| 1.0 | 0.823 | 0.333 | 1.454 | 0.333 | **True** | True | True | True |
| 0.8 | 0.904 | 0.293 | 1.501 | 0.293 | False | True | True | True |
| 0.6 | 0.831 | 0.445 | 1.577 | 0.445 | **True** | True | True | True |
| 0.4 | 0.739 | 0.306 | 1.306 | 0.306 | False | True | True | **False** |
| 0.2 | 0.765 | 0.384 | 1.319 | 0.384 | False | True | True | **False** |

(threshold = 0.3 throughout, unchanged from Stage 1.)

Connectivity diagnostics per soil per theta (`metrics.connectivity_diagnostics`):
largest connected water cluster as a fraction of water-filled cells,
OM-biomass connectivity (fraction of *initial* biomass sitting in a water
component that also contains a high-OM cell, "high" = top quartile), and
accessible-OM fraction (total OM reachable from a high-K cell through the
water network, over total OM). See `results/stage2_connectivity_vs_theta.png`.

- **Soil B** (coarse, `mu=2.3`, salt-and-pepper): stays robustly connected
  even at `theta=0.2` — `largest_cluster_frac` only drops to ~0.36,
  `accessible_om_fraction` stays at ~0.60. Its coarse, spatially
  uncorrelated pore field percolates easily because there is no single
  bottleneck length scale for the water network to fragment along.
- **Soil A** (fine, `mu=0.6`, `lambda=4.0`, correlated): `largest_cluster_frac`
  stays near 1.0 (its large-scale spatial correlation means the
  water-filled fine pores stay contiguous even as coverage shrinks — theta
  is low but the *shape* of the remaining wet region is still one big
  patch), but `om_biomass_connectivity` and `accessible_om_fraction` both
  collapse fast: `accessible_om_fraction` hits **0.0** by `theta=0.6`. A
  large connected water blob is not the same as OM and biomass sharing it.
- **Soil C** (aggregated, `lambda=1.0`): fragments catastrophically —
  `largest_cluster_frac` collapses from 0.83 (theta=0.6) to 0.10
  (theta=0.4), and `accessible_om_fraction` hits 0.0 by theta=0.6 as well.
  Unlike A, this soil's *pore network itself* falls apart, not just the
  OM/K overlap, because `aggregate=True` deliberately carves fine-pore
  boundaries around each ped — exactly the domain-boundary cells that
  disconnect first as theta drops.

## 4. Honest verdict (§6)

**Percolation is a real, and in soil C's case an extreme, lever — but it did
not turn Stage 1's weak A–C pair into a robustly separated one.** Success
(`all pairwise distances > 0.3`, all soils non-trivial) is achieved at only
2 of the 5 sampled theta values (1.0 and 0.6), and fails at 0.8, 0.4, and
0.2 — non-monotonically, not via a clean "wetter succeeds, drier fails"
trend. The mechanism is visible directly in the connectivity diagnostics:

- At `theta=0.8`, AC (0.293) misses the threshold by a hair — both A and C
  are still nearly fully connected structurally but their curves happen to
  be numerically close at that particular saturation.
- At `theta=0.4` and `0.2`, soil C's biomass net-declines
  (`B_final_total < B0_total`, i.e. **`nontrivial=False`**) because its
  aggregated fine-pore boundaries fragment the water network so badly that
  most of its high-K cells lose access to any OM at all
  (`accessible_om_fraction=0.0`). This is arguably the *most* distinct
  outcome physically possible (soil C effectively dies while A and B keep
  respiring) — but the success metric explicitly requires every soil to be
  non-trivial, so a die-off does not count as "success" even though it is a
  genuine, structurally-caused divergence. This is a real limitation of the
  binary success criterion, not a bug: "distinct" and "all alive" are
  different axes, and percolation can satisfy the first while violating the
  second.
- At `theta=0.6`, both A and C have already lost most of their OM-biomass
  connectivity (0.072 and 0.187 respectively) but neither has fully died
  yet, and their curves separate more than at full saturation (AC=0.445 vs
  0.333) — the one point in the sweep where percolation clearly helped.

So the answer to the Stage-2 research question is a genuine mixed result,
not a clean win: percolation gives soil C (and to a lesser extent soil A) a
much stronger structural lever than Stage 1's fully-wet transport ever did —
strong enough, in fact, to kill soil C outright at low theta — but that
strength does not translate into a *reliably* wider AB/AC/BC margin across
the whole drying range. The AC pair, weak in Stage 1, remains the fragile
pair in Stage 2, just for a different underlying reason: instead of both
soils looking similar because transport was fast everywhere (Stage 1), they
now sometimes look similar because both are becoming substrate-starved
together as theta drops (Stage 2). Soil B, by contrast, is robust in both
stages, for the same reason both times: its coarse, uncorrelated pore field
has no natural bottleneck to fragment along.

## 5. Deliverables produced

- `results/stage2_theta_sweep.csv` — R_peak, R_end, cumulative CO2,
  B_final_total, non-triviality, and all three connectivity diagnostics,
  per soil per theta.
- `results/stage2_distinctness_by_theta.csv` — pairwise distances and
  success flag per theta.
- `results/stage2_water_map_{A,B,C}_theta0.40.png` — pore radius, K, D, and
  water-filled mask at the chosen `theta=0.4`.
- `results/stage2_respiration_curves_theta0.40.png` — R(t) for all three
  soils at `theta=0.4`.
- `results/stage2_respiration_vs_theta.png` — peak R and cumulative CO2 vs
  theta, per soil.
- `results/stage2_connectivity_vs_theta.png` — the three connectivity
  diagnostics vs theta, per soil.

## 6. Seed note

Per `soil_respiration_prompt_stage2.md` §1 ("different soils SHOULD keep
different seeds... do NOT force a shared seed across the three soils"), the
seeds were **not** changed from Stage 1 (A=1, B=2, C=3). §7's "shared seed"
bullet is read here as "keep the seed shared with Stage 1" (i.e. don't
introduce a new one), not "use one seed across A/B/C" — the two would
otherwise directly contradict §1, which is the more explicit and repeated
instruction (also restated in §9's closing summary).

## 7. Limitations (stage-2 specific, in addition to Stage 1's)

- Only the pre-existing Stage-1 pore parameters were tested; no Stage-2
  parameter search/nudge loop was run, since the sweep alone already
  answers the research question (does percolation help, and how) without
  needing to hunt for a lucky configuration — consistent with §1's
  instruction not to cherry-pick.
- The "high-OM"/"high-K" percentile threshold (top quartile) used by the
  connectivity diagnostics is a documented, reasonable choice, not derived
  from data; a different percentile would shift the absolute connectivity
  numbers (though not the qualitative A/B/C ranking, which is driven by
  the underlying pore topology).
- The theta grid (1.0, 0.8, 0.6, 0.4, 0.2) is coarse; the non-monotonic
  success pattern could look different at finer resolution, but the
  qualitative finding (AC stays fragile, C is the most percolation-fragile
  structure, B is the most robust) is unlikely to change since it tracks
  each soil's spatial correlation length directly.
- A robustness/significance pass (seed ensembles per soil, testing whether
  different arrangements of the same texture statistics converge) remains
  explicitly deferred, as in Stage 1.

## 8. Relevant files

- `pore_field.py` — numpy-only pore field + K/D maps (scipy removed).
- `model.py` — `water_mask_and_theta` (`global_theta` mode), mass-conserving
  diffusion with `D=0` percolation gating.
- `metrics.py` — `label_connected`, `largest_cluster_fraction`,
  `om_biomass_connectivity`, `accessible_om_fraction`,
  `connectivity_diagnostics`; unchanged `emergent_distinctness`.
- `stage2_run.py` — the theta sweep runner and all Stage-2 figures/CSVs.
- `configs/soil_{A,B,C}.yaml` — `saturation.mode: global_theta` (unchanged
  pore params, unchanged seeds).
- `MODEL_SPEC.md`, `LOGBOOK.md` — updated with Stage-2 detail.
- `results/` — all Stage-2 figures and CSVs (see §5).
