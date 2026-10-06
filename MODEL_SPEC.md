# MODEL_SPEC (v2 -- pore-based soil environment)

## 1. Purpose

A single microbial species grows logistically on a 2D grid. The biology is
identical everywhere and frozen across all soils, exactly as in v1. Unlike
v1, the three soils (A, B, C) no longer differ through hand-picked scalars
(a single D, hand-placed OM/colony centers). Instead, every cell has a
**pore radius** `r(i,j)`, and both the local carrying capacity `K` and the
local diffusion rate `D` are FUNCTIONS of that pore field. The three soils
are genuinely different porous media (different pore-size distributions and
spatial arrangements); their respiration differences must **emerge** from
that, not be scripted toward a target shape.

This spec covers **Stage 1** (the pore field feeds K and D, soil fully wet,
every cell conducts), **Stage 2** (water-filled pores via
`saturation.mode: global_theta`, percolation/connectivity diagnostics --
SS9, SS13-14), and **Stage 3** (real measured pore-size distributions,
physiological carrying capacity, bigger grids -- SS15).

## 2. Pore field (`pore_field.py`, SA)

Each `configs/soil_{A,B,C}.yaml` specifies a pore field via:

- **Texture**: `mu`, `sigma` -- the underlying lognormal distribution of
  pore radius (mean/spread of `ln r`).
- **Structure**: `lambda` -- the spatial correlation length of a Gaussian
  random field (small = salt-and-pepper, large = big homogeneous domains).
- **Aggregate** (`true`/`false`): if true, thresholds the field's
  high-gradient (domain-boundary) cells down to fine pores, producing peds
  separated by fine-pore boundaries.

Generation: draw white noise, Gaussian-smooth with kernel width `lambda`,
standardize to mean 0 / std 1, map through `r = exp(mu + sigma*z)`. Same
fixed seed per soil for reproducibility. Soils differ ONLY in
`(mu, sigma, lambda, aggregate, seed)` -- nothing else.

## 3. Pore size -> local parameters (`pore_field.py`, SB)

```
K(r) = k_min + (K_max - k_min) * exp( -0.5 * ((ln r - ln r_opt) / w_K)^2 )
D(r) = D_min + (D_max - D_min) * (r / (r + r_ref))^p
```

- `K(r)` peaks at `r_opt` (mid-size pores best) and is floored at `k_min` for
  extreme pore sizes -- both too-small (clogged) and too-big (dry/poorly
  retentive) cells are near-inhospitable.
- `D(r)` rises with pore size and **saturates** at `D_max` as `r -> infinity`
  (a literal `(r/r_ref)^p` power law would be unbounded and could push a
  coarse/large-`mu` soil's diffusion above the numerical stability limit).
- All seven constants (`r_opt, w_K, K_max, k_min, D_min, D_max, r_ref, p`)
  live in `configs/mapping.yaml`, fixed and shared by every soil and every
  stage -- only the pore field differs between soils, so K and D maps differ
  only because the pore fields differ.

`w_K` matters a lot: it controls how *selective* habitat is. A broad hump
(K nonzero almost everywhere) lets every cell hold some biomass and OM
locally, so no transport is required for basic respiration and structure
barely matters (see the calibration note in `LOGBOOK.md`). Stage 1 uses
`w_K = 0.45`, narrow enough that habitat genuinely concentrates.

## 4. Rule-based OM and biomass (`pore_field.py`, SC)

- **Initial biomass**: `B0(i,j) = B0_total * K(i,j) / sum(K)` -- the
  community seeds itself wherever habitat is good (proportional to local K),
  not at typed-in centers.
- **Initial OM**: `OM(i,j) = OM_total * r(i,j)^-om_fine_bias / sum(...)` --
  organic matter is protected in fine pores (weight rises as `r` shrinks),
  independent of where biomass is. Same rule and exponent (`om_fine_bias`,
  fixed in `mapping.yaml`) for all soils, so OM layout differs only because
  pore fields differ.
- `OM_total` and `B0_total` are fixed shared totals in `mapping.yaml`
  (analogous to the biology constants: invariant across soils, not
  per-soil knobs).

## 5. Variable-D diffusion (`model.py`, SD)

D is a spatial field, so substrate transport is conservative-flux finite
volume, not a plain Laplacian:

```
dS/dt = div( D * grad(S) ) + k_dis * OM - uptake
```

Face conductivities use the **harmonic mean** of the two neighboring cells'
D (`model.harmonic_mean`) -- correct for sharp D contrasts and for D -> 0
(percolation, stages 2/3). Zero-flux (Neumann) boundaries via zero-padded
face fluxes. Because face fluxes telescope, `sum(div(D*grad(S)))` over the
whole grid is exactly the (zero) boundary flux -- diffusion alone can never
create or destroy mass. Verified numerically in `LOGBOOK.md` (`mass_after -
mass_before = 0.0` exactly, including with disconnected cells).

Full per-cell equations (biology unchanged from v1):

```
dB/dt  = r_max * f(S) * B * (1 - B/K_ij) - m0*B - m_s*(1 - f(S))*B
dS/dt  = div( D_ij * grad(S) ) + k_dis * OM_ij - uptake
dOM/dt = -k_dis * OM_ij

f(S)    = S / (S + Ks)
uptake  = (1/Y) * r_max * f(S) * B * (1 - B/K_ij)
```

Respiration observable (unchanged):

```
R_ij(t) = (1 - Y) * [substrate used for growth in cell ij] + m0 * B_ij
R(t)    = sum over all cells of R_ij(t)
```

Integrated with explicit Euler; states clipped at >= 0. Numerical stability:
`max(D) * dt / dx^2 < 0.2` (enforced by an assertion in `model.simulate`;
`D_max=3.0` and `dt=0.05` give `0.15 < 0.2` with headroom).

## 6. Frozen biology constants (never tuned)

From `configs/biology.yaml`, shared by all three soils, all stages:

```
r_max = 1.00     Ks    = 0.20     Y     = 0.40
m0    = 0.02     m_s   = 0.10     k_dis = 0.05
dt    = 0.05
```

## 7. Fixed mapping constants (never tuned per soil)

From `configs/mapping.yaml`, shared by all three soils, all stages:

```
r_opt = 5.0    w_K = 0.45   K_max = 1.0   k_min = 0.005
D_min = 0.02   D_max = 3.0  r_ref = 5.0   p = 2.0
matric_r_cut = 5.0 (stage 3 retention threshold, see below)
OM_total = 800.0   B0_total = 20.0   om_fine_bias = 2.5
```

## 8. Structure knobs (the ONLY things that differ between soils)

Each `configs/soil_{A,B,C}.yaml` sets:

- **`pore`**: `mu, sigma, lambda, aggregate, seed` -- see SS2.
- **`saturation.mode`**: `fully_wet` (Stage 1), `global_theta` (Stage 2),
  or `retention` (Stage 3) -- see the per-stage sections below.
- **`grid.n`** (<=50), **`grid.dx`**, **`T`** -- kept equal across soils
  within a stage (40x40, dx=1.0, T=600 for Stage 1).

## 9. Saturation / water-filled pores

`model.water_mask_and_theta(r, saturation_cfg, mapping_cfg)`:

- **`fully_wet`** (Stage 1): every cell conducts (`theta=1.0`).
- **`global_theta`** (Stage 2): same nominal `theta` for every soil;
  `r_cut` is set per soil as THAT soil's own `theta`-quantile of its pore
  field (small pores fill first), so cells with `r <= r_cut` are
  water-filled. Same saturation, different percolation topology, because
  the pore *arrangement* differs.
- **`retention`** (Stage 3): the SAME physical `r_cut` (`matric_r_cut` in
  `mapping.yaml`, standing in for a shared matric potential) for every
  soil; `theta` then EMERGES per soil as the fraction of that soil's own
  pore field with `r <= matric_r_cut` -- fine soils (small `mu`) retain more
  water at the same suction, exactly as in real soil-water retention curves.

In all modes, `D = D(r)` where water-filled, `0` elsewhere; harmonic-mean
faces mean a `D=0` cell truly blocks transport through it (percolation).
Stage 3 keeps `global_theta` (see SS15) -- `retention` remains unused so far.

## 10. New success criterion -- emergent distinctness (`metrics.py`)

Not "A rising, B falling, C flat." Instead:

1. Run the three soils, get `R(t)` each.
2. Normalize each by its own peak; pairwise distance
   `dist = (1 - pearson_corr(norm_Ri, norm_Rj)) + 0.5 * mean(|norm_Ri - norm_Rj|)`.
3. **Success** = all three pairwise distances `> 0.3` AND every soil is
   non-trivial (`B_final_total > B0_total` and `sum(R_t) > 0`).
4. Each curve's shape is separately *described* (rising / falling / flat /
   single-peak / multi-peak, via trend sign + peak count) and logged, but
   never used to steer the search.

## 11. Stage 1 result

Achieved at iteration 1 of the budget-capped search (see `LOGBOOK.md`):
`AB=0.82, AC=0.33, BC=1.45` (all `> 0.3`), all three soils non-trivial.
Soil A (fine, `mu=0.6`, moderately correlated) rises slowly across the
window; soil B (coarse, `mu=2.3`, salt-and-pepper) bursts early and decays
as its accessible OM depletes; soil C (intermediate, aggregated, short
correlation length) dips then rises as its patchy high-K network only
gradually accumulates substrate. See `results/structure_maps_{A,B,C}.png`
for the pore/K/D fields and `results/respiration_curves.png` /
`cumulative_co2.png` for the resulting curves.

**Robustness caveat (added after independent review).** Each soil uses a
*different* seed, which confounds soil identity with the particular random
field that got drawn. Resampling the A-C pair over 64 seed combinations gives
a mean distance ~0.16 and a maximum ~0.29 -- i.e. A and C fall *below* the 0.3
threshold for essentially every seed pair, while two realizations of the *same*
soil differ by only ~0.03. The configured `AC=0.33` is an upper-tail outlier
picked out by the offline seed/parameter scan, not a robust structural
difference. Soil B, by contrast, is robustly distinct (high mean pore size ->
fast transport -> early burst then decay). So at Stage 1 only **B** separates
reliably; A and C came out similar. For the project's aim -- do different
structures give different respiration? -- that is simply an honest *negative* for
the A-C pair at full saturation. Whether the closeness reflects the structures
themselves or just this realization is a question an optional later significance
analysis (seed ensembles) could settle. Stage 2 adds the stronger, less-tunable
percolation lever to see whether more structures then diverge.

## 12. Assumptions / simplifications (Stage 1)

- Fully wet: no percolation gating yet (all cells conduct) -- Stage 2/3
  address this.
- The pore-field generator's Gaussian-smoothing + lognormal-quantile recipe
  is one reasonable choice among several valid ways to synthesize a
  spatially-correlated lognormal field; not calibrated against a real soil
  dataset.
- `w_K`, `om_fine_bias`, and the `D(r)` saturation form were calibrated once
  (documented honestly in `LOGBOOK.md`) so that pore structure has a real
  causal lever on respiration at all -- this is model design, done before
  and independent of the per-soil pore-parameter search.
- The starting `(mu, sigma, lambda, aggregate)` triples were found by a
  small offline grid scan over the viable-habitat `mu` band (see
  `LOGBOOK.md`), not purely by the iterative nudge loop -- the nudge loop
  is fully implemented and will keep iterating if a future stage's starting
  point does not already clear the threshold.
- Deterministic throughout: the only randomness is the fixed-seed pore-field
  generation per soil.
- **Seeds are intentionally different per soil.** The three soils are meant to
  be genuinely different structures, so they use different seeds, and the aim is
  only whether different structures give different respiration -- single
  realizations are acceptable. A single-realization distance does mix structure
  with the particular draw; an optional later significance analysis (seed
  ensembles; testing whether different arrangements of the same texture
  statistics converge to the same respiration) could quantify that, but it is not
  required for the current aim and is deliberately deferred.

## 13. Numpy-only (from Stage 2)

The repo depends only on `numpy` + `pyyaml` + `matplotlib` -- no scipy.
`pore_field.gaussian_filter_2d` is a hand-rolled separable-Gaussian
convolution (1D kernel along each axis, `radius = int(4*sigma + 0.5)`)
that reproduces `scipy.ndimage.gaussian_filter`'s default behavior
(`mode="reflect"`, `truncate=4.0`) to machine precision -- the one
subtlety is that scipy's `mode="reflect"` duplicates the edge value, which
is `numpy.pad`'s `"symmetric"` mode, not numpy's own (non-duplicating)
`"reflect"` mode. See `STAGE2_RESULTS.md` SS1 for the verification.

## 14. Stage 2 -- global_theta result

Stage 2 (`saturation.mode: global_theta`, `model.water_mask_and_theta`)
swept `theta` in {1.0, 0.8, 0.6, 0.4, 0.2} for all three soils, keeping
each soil's Stage-1 pore parameters and seed unchanged. Success (all
pairwise distances `> 0.3`, all soils non-trivial -- still a single
realization per soil/theta, per SS10) was reached at only 2 of the 5 theta
values (1.0, 0.6); it failed at 0.8 (AC=0.293, just under threshold) and at
0.4/0.2 (soil C's biomass net-declines under severe percolation
fragmentation, failing the non-trivial requirement even though its curve is
maximally different from A and B). The percolation lever is real and, for
soil C, extreme -- but it does not turn Stage 1's weak A-C pair into a
robustly separated one; the pair remains fragile for a different reason
(both soils becoming substrate-starved together as theta drops, rather than
both having fast-enough fully-wet transport to look alike). Soil B stays
robustly distinct at every theta, as it was in Stage 1, because its coarse,
spatially uncorrelated pore field has no bottleneck length scale to
fragment along. Full detail, the theta-by-theta table, and the connectivity
diagnostics (`metrics.connectivity_diagnostics`: largest water cluster
fraction, OM-biomass connectivity, accessible-OM fraction) that explain the
mechanism are in `STAGE2_RESULTS.md`.

## 15. Stage 3 -- real PSD, physiological K, bigger grid

Three changes replace the synthetic, hand-tuned Stage 1/2 pore field with a
real, measured one, and revise the success metric accordingly.

**Pore field (`pore_field.py`, `psd_data.py`).** Each soil's real measured
pore-size distribution (`data/psd/<name>/psd_table.csv` + `bin_edges_um.csv`,
copied once from an external diagnostics pipeline; own bins, never re-binned;
volume-weighted; diameters in microns, not radii) replaces the lognormal
quantile function. The spatially-correlated Gaussian field generation is
unchanged and factored out (`pore_field._correlated_gaussian_z`, shared by
both paths); it is converted to uniform quantiles via an argsort rank
transform (`rank_uniform`, numpy-only, exact for the field's own empirical
distribution) and mapped through the soil's empirical inverse-CDF
(`generate_diameter_field_from_psd`). Only the pore-size *distribution* is
real; the spatial *arrangement* recipe (`lambda=2.0`, `aggregate=false`) is
synthetic and, per invariant, **shared by all three soils** -- the only
systematic difference between soils in Stage 3 is their real PSD. Seeds
stay different per soil (realization noise, not a lever).

**Physiological K(d) (`pore_field.K_window`).**

```
K(d) = 0                                        for d < 30
K(d) = K_max                                    for 30 <= d <= 150
K(d) = K_max * exp( -((d - 150) / w_decay)^2 )    for d > 150   (no floor)
```

replacing the lognormal hump: a hard-zero non-habitat floor below 30 um
(too tight for microbial activity), a full plateau across the
microbially-relevant 30-150 um window, and a smooth Gaussian decay above
150 um toward (not floored at) zero. `K=0` cells still hold water and
conduct -- conduits, not homes. `K_d_low, K_d_high, K_w_decay` are fixed,
shared constants in `configs/mapping.yaml`.

**Rescaled D(d) (`pore_field.D_of_d`).** Same saturating functional form as
Stage 1/2's `D_of_r`, `D(d) = D_min + (D_max-D_min)*(d/(d+D_d_ref))^p`, with
a diameter-scale reference `D_d_ref=90` um (mid-window) instead of the
radius-scale `r_ref=5.0`. `D_min`, `D_max`, `p` are reused unchanged.
`max(D)*dt/dx^2 < 0.2` still holds (0.15, same as Stage 1/2) because `D(d)`
saturates at the same `D_max` regardless of the diameter scale involved --
**this saturation is the mechanistic reason two of the three real soils
turn out to look alike in Stage 3 (result below), since both have median
diameters far above `D_d_ref` and are therefore both already near `D_max`
almost everywhere.**

**Grid size.** The `n<=50` assertion (`model.build_grids`) was relaxed to
`n<=200` to allow `n=128` alongside the Stage 1/2 `n=40`, swept side by
side.

**Metric revision (`metrics.py`).** `classify_outcome` labels each run
`thriving`/`declining`/`dead` from final-vs-initial biomass and the
late-window respiration trend. `emergent_distinctness` gained a
`require_nontrivial` flag (default `True`, preserving Stage 1/2 behavior
exactly); Stage 3 calls it with `require_nontrivial=False` so a
structurally-caused die-off is reported, never treated as a disqualifier.
`isolation_ratio` (`n_isolated_water_clusters / n_saturated_cells`, via the
existing `label_connected` BFS) is a new connectivity diagnostic, added to
`connectivity_diagnostics`.

**Stage 3 result.** Full sweep, mechanism, and honest verdict in
`STAGE3_RESULTS.md`. Headline: at the shared spatial recipe used, the AB and
AC pairs partially separate (best at `n=128`, `theta=0.2`: AB=0.430,
AC=0.256 -- AC still short of 0.3) but the BC pair never approaches the
threshold (0.048-0.062 at any theta/grid) because both soils' large median
diameters put their `D(d)` deep in its saturated regime, compressing their
real PSD difference away. No soil died off in this sweep (Stage 2's
die-off mechanism, `aggregate=True`, is off by Stage-3 invariant). The
isolation ratio correlates with both cumulative CO2 and peak R more
strongly than the PSD summary stats tested (`|r|` ~0.25-0.33 vs ~0.14-0.24),
a real but modest edge -- SS9's hypothesis holds, weakly. `n=128` produces
substantially more distinct respiration across soils than `n=40`, where
cumulative CO2 is nearly flat regardless of soil or theta (the small domain
is close enough to spatially homogeneous, at `lambda=2.0`, that OM depletes
at nearly the same rate everywhere).

## Stage 3.5 -- literature-informed PSDs at constant moisture

Stage 3.5 keeps the entire Stage-3 machinery and swaps only the PSD SOURCE:
each soil's PSD is a **mixture of lognormal modes over pore diameter (um)**
built by `psd_parametric.build_soil_psd` from `configs/psd_literature.yaml`,
selected with `pore.mode: psd_parametric` (`pore.psd_soil` names the soil).
The builder returns the same dict shape as `psd_data.load_soil_psd`, so the
empirical-quantile pore-field mapping, `K(d)`, `D(d)`, and all PSD summary
stats are unchanged. Soils relabelled **Loess / Sand / Vertisol**.

Because the experiment was at CONSTANT moisture (closed humid jars), Stage 3.5
does NOT sweep theta: it runs each soil at a single fixed theta held for the
whole simulation and reads the temporal shape of R(t). Two saturation modes
are relevant:
- `global_theta` -- the SAME theta for all soils (the headline run, theta=0.6);
- `matric_d` -- the same air-entry DIAMETER `d_cut` (um) for all soils (a
  shared matric potential, `mapping.yaml matric_d_cut`); theta then EMERGES per
  soil from its own PSD. This is the closed-jar-faithful diagnostic.

**D(d) fusion fix (SS3):** `mapping.yaml D_d_ref` raised 90 -> 250 um so D(d)
keeps discriminating across ~30-300 um (D(300)/D(30): 8.6x -> 16.7x) and the
two coarse soils Stage 3 fused no longer collapse. `max(D)*dt/dx^2 = 0.126 <
0.2` still holds.

**Numerical fix:** the physiological `K(d)` is a hard zero below 30 um; those
cells are non-habitat. `model.simulate` now masks the logistic term `1 - B/K`
where `K == 0` instead of dividing (Stage 3's measured PSDs never hit exactly
K=0; the literature PSDs do).

**New metrics** (`metrics.py`): `habitat_fraction` (K>0), `wet_habitat_fraction`
(wet AND K>0 -- the variable that actually decides the outcome at fixed theta),
and `habitat_wet_share`.

**Stage 3.5 result.** Full detail, mechanism, and honest verdict in
`docs/stage_results/STAGE3_5_RESULTS.md`. Headline: the model reproduces the
measured shape for **Sand only (falling), and not robustly** (it flips to
rising at n=40); it does NOT reproduce Vertisol-rising or Loess-erratic. At a
shared theta small pores fill first, so the Vertisol's and Loess's habitat
macropores are the last to wet -- bone dry at theta=0.6 (wet-habitat=0.000), so
their biomass starves and both give the identical pure-starvation curve. The
shared-matric-potential diagnostic gives all three a live wet habitat yet still
reproduces only Sand: the frozen `OM ~ d^-2.5` rule, across the literature
PSDs' 7-order-of-magnitude range, concentrates ~99.96% of OM on sub-micron clay
cells and diffusion (reach ~1-3 cells) never delivers substrate to the habitat.
Main limitation: substrate never reaches the habitat (OM-placement rule +
every-cell-is-a-pore, no bulk porosity).

## Stage 4 (EXPLORATORY) -- two fixes for "OM never reaches the habitat", head-to-head

`prompts/soil_respiration_prompt_stage4_exploratory.md`. A defined comparison
(`MAX_ITERATIONS=12` not exercised), not a search: build two independent
fixes for Stage 3.5's diagnosis, run both, compare against the measured
pattern, without tuning either toward the target.

**Variant A ("truncate + relocate OM", stays single-continuum).**
`pore_field.truncate_psd_floor(psd, floor_um)` resamples a PSD's empirical
CDF onto a fresh log grid from `floor_um` up to the original max diameter,
renormalizing so all volume below `floor_um` is dropped -- numerically this
is soil-agnostic (works on measured or parametric PSDs alike) and returns a
report (`dropped_volume_fraction`, before/after porosity). `pore_field.
om_field_banded(d, total, band_center_um, band_width_log)` replaces the
unbounded `OM = total * d^-om_fine_bias / sum(...)` rule with a Gaussian band
in LOG-diameter space, `weight = exp(-0.5*((ln d - ln center)/width)^2)`,
peaking at `om_band_center_um` (10-30um, just below the habitat). Both are
wired into `model.build_grids` via two new OPTIONAL `pore` config keys,
`psd_truncate_floor_um` and `om_rule: banded` (default `fine_bias`, i.e. a
no-op for every Stage 1-3.5 config, which sets neither key).

**Variant C (dual-porosity matrix/macropore, new module `dual_porosity.py`).**
Keeps the FULL, untruncated literature PSD. Every cell is classified by the
SAME physiological threshold already used for `K_window` (`K_d_low=30um`,
`mapping_cfg["K_d_low"]`, no new constant): `macro_mask = (d >= 30um)` (fast
`D(d)`, biomass lives and respires here, exactly Stage 1-3.5's biology) or
`matrix_mask = ~macro_mask` (`K=0` here already, by construction of
`K_window`; its fast diffusion pathway is explicitly zeroed,
`D_full = where(macro_mask, D_of_d(d), 0)`, so a matrix cell can never
participate in the ordinary conservative-flux diffusion term at all). Matrix
"internal fine porosity" is NOT a separate scalar field -- the existing
fine-pore-biased `om_field(d, OM_total, om_fine_bias)`, evaluated over the
FULL (untruncated) diameter field, already concentrates most OM on matrix
cells (matrix == fine pores by definition), so no new OM-placement rule is
needed for this variant; `grids["matrix_fraction"]` (mean of `matrix_mask`)
is reported per soil as the bulk-porosity proxy the prompt asks for.

The only new physics is a mobile-immobile exchange term added to `dS/dt`:
`dual_porosity.leak_exchange(S, matrix_mask, macro_mask, k_leak)`. Built as a
face-based flux exactly like `model.conservative_divergence` (each matrix
cell's outgoing flux to an adjacent macropore cell is exactly that
macropore cell's incoming flux, so `sum(leak_exchange(...)) == 0` always,
mass-conserving regardless of grid shape or mask): for every face where one
side is `matrix` and the other `macro`, `flux = k_leak * S_matrix_side`
leaves the matrix cell and enters the macropore cell; matrix-matrix and
macro-macro faces contribute nothing. This is NOT gated by `water_mask` --
a matrix cell can leak into a macropore cell even if that macropore cell is
currently dry, representing capillary-scale matric flow that the coarse
`D=0`-blocks-everything convention would otherwise forbid. `k_leak=0.10`
(`configs/mapping.yaml`, same order as `k_dis=0.05`) is a fixed, shared
constant, not a per-soil knob. Stability: `dt*4*k_leak < 1` (worst case, a
matrix cell with 4 macropore neighbors), checked alongside the existing
`D.max()*dt/dx^2 < 0.2` diffusion bound.

Both variants also add `S_habitat_mean_t` (mean substrate concentration
across habitat cells, `K>0`, recorded every timestep) to `simulate`'s output,
and `metrics.mean_substrate_in_habitat(S, K)` reports it at any single time
-- the Stage-4 "killer metric" (SS4): whether carbon actually reaches the
microbes, not merely whether the habitat is wet (Stage 3.5's
`wet_habitat_fraction` answers a necessary but not sufficient question).

**Stage 4 result.** Full detail, mechanism, and honest verdict in
`docs/stage_results/STAGE4_RESULTS.md`; all figures/CSVs in
`results/stage4/`. Headline: **Variant A suffices better; Variant C does
not.** Variant A's PSD truncation has an unplanned but decisive side effect
-- it shifts `global_theta`'s per-soil quantile cutoff toward the coarse end
(the fine tail is gone from the distribution), which rewets the 30-150um
habitat at the SAME nominal theta=0.6 that left it bone-dry in Stage 3.5
(wet-habitat: Loess 0.000->0.136, Vertisol 0.000->0.523). All three soils
come out thriving; 2/3 shapes match the measured pattern (Loess-erratic,
Vertisol-rising) robustly across theta in {0.5,0.6,0.7} at n=128, though
Sand robustly comes out wrongly "rising" (not falling) and the Vertisol
match is not grid-robust (flips to "erratic" at n=40). Variant C leaves the
theta-quantile calculation (and therefore the wetting order) completely
untouched -- Loess's and Vertisol's habitat stays bone dry
(`wet_habitat_fraction=0.000`) exactly as in Stage 3.5, its `k_leak`
exchange delivers ~4 orders of magnitude less substrate to the habitat than
Variant A (bounded by how much matrix OM sits directly adjacent to a
macropore cell), and both soils reduce to the exact same degenerate
dead-curve collapse Stage 3.5 already reported (`Loess_vs_Vertisol`
distance = 7.4e-6). Sand -- the one soil Stage 3.5 already reproduced
unaided, whose macropores are wet regardless -- is the only soil where
Variant C's leak pathway has room to matter, and it does. The dual-porosity
idea's fix (unblocking fast intra-macropore diffusion) targets a bottleneck
that was never binding for the two soils that actually needed help; the
real, prior blocker (shared-theta wetting order) is untouched by
construction and Variant A only clears it as an accident of its PSD
surgery, not by design.

**Note on Variant C's later evolution.** Since this write-up, `dual_porosity.py`
was substantially extended beyond the leak-only version described above: a
three-phase split (macropore / porous-matrix / inert solid, `_assign_solid_
phase`), per-soil initial-OM fractions (`stage4_variantC_om_fraction_*`), a
low-rate internal matrix diffusion (`stage4_variantC_matrix_diffusion`),
dynamic hydro "wicking" (`apply_wicking`, `wicking_substrate_exchange` --
dry macropores can rewet if wet porous matrix is nearby, within
`stage4_variantC_wicking_radius`), habitat-cluster-size K/D/S bonuses
(`_cluster_scale_field`, the `stage4_variantC_cluster_*_scale` constants),
and a per-macropore-cell substrate holding capacity (`S_cap`,
`stage4_variantC_S_capacity_base`). This is the baseline Stage 5 (below)
builds on; see `dual_porosity.py`'s own docstrings for the current mechanics
in full, and `STAGE5_RESULTS.md` SS4 for why `S_cap` turned out to be a load-
bearing bottleneck independent of anything Stage 5 changed.

## Stage 5 (EXPLORATORY) -- matric-potential wetting + inscribed-circle macropores

`prompts/soil_respiration_prompt_stage5_matric_inscribed.md`. Builds on the
latest Variant C (above) and changes exactly two things, both in
`dual_porosity.py`.

**Change 1 -- matric-potential wetting (`build_dual_grids`, new `saturation.
mode: matric_regions`).** Replaces the per-soil `global_theta` quantile with
a single SHARED air-entry diameter, `matric_d_cut_um` (`configs/mapping.yaml`,
default 50 -- distinct from the pre-existing `matric_d_cut`, Stage 3.5's
single-continuum `model.py` "matric_d" diagnostic, value 60, a different
mechanism kept separate). A pore is wet if its SIZE is below the cutoff:
matrix cells use their own PSD diameter `d` (always < `K_d_low`=30um < 50, so
always wet -- the fine matrix stays saturated by construction); macropore
cells use their REGION's inscribed diameter (Change 2), not the per-cell PSD
label; solid cells are never wet regardless of size. `theta` then EMERGES
(`water_mask.mean()`) instead of being imposed. The pre-existing
`global_theta` path (`model.water_mask_and_theta`) is untouched for any
`saturation.mode` other than `matric_regions`, so Stage 4's runs still work.

**Change 2 -- geometric macropore sizing (`dual_porosity.macropore_region_
inscribed_diameter`, `dual_porosity.erosion_distance_4n`).** Retires the flat
`stage4_variantC_macro_diameter_um=80` constant that forced every macropore
cell to the same diameter. `erosion_distance_4n` is a numpy-only iterative
4-neighbor ("+"-shaped structuring element) erosion distance transform:
`dist[i,j]` = the number of erosion rounds cell `(i,j)` survives, which for a
"+"-shaped element is exactly the Manhattan (L1) distance to the nearest
non-mask cell or the grid edge (out-of-grid = background via zero-padding).
`macropore_region_inscribed_diameter` labels connected macropore regions
(reusing `dual_porosity._component_labels`) and sizes region `k` as
`2 * max(dist over region k) * voxel_um` (`voxel_um=10` um/cell,
`configs/mapping.yaml` -- the grid's physical scale). This per-cell "region
diameter" field feeds BOTH the wet/dry decision (Change 1) and `K(d)`/`D(d)`
on macropore cells (`d_struct = where(macro_mask, region_diam_field, d)`), so
a small water-holding pore and a large drained channel are no longer treated
identically.

Also added: `simulate_dual_porosity(..., wicking_enabled=True)`. Setting it
`False` freezes `wet_mask` at exactly Change 1's static map for the entire
run (no dynamic rewetting), isolating wicking's own contribution --
necessary because wicking's rewetting radius reaches almost every macropore
cell once the surrounding matrix is (by Change 1's own construction) almost
always wet, which can blur or override the static wet/dry distinction Change
1 introduces. `k_leak` (the separate, always-on matrix->macropore exchange)
is unaffected by this toggle either way.

**Stage 5 result.** Full detail, mechanism, and honest verdict in
`docs/stage_results/STAGE5_RESULTS.md`; all figures/CSVs in
`results/stage5/`. Headline: **Change 1+2 do give Loess and Vertisol
genuinely live wet habitat where Stage 4 gave exactly zero** (under the
static map alone, `habitat_wet_share` 0.074/0.077 vs Stage 4's 0.000), **and
correctly identify Sand's macropore network as one huge (81% of the grid),
well-connected, 340um-inscribed-diameter region that is too large to hold
water at all** (`habitat_wet_share=0.000` under the static map) -- a real
structural finding, not a bug. Wicking (part of the inherited baseline)
substantially undermines this distinction in practice, rewetting most of what
the static map marks dry (`habitat_wet_share` rises to 0.595/0.570/0.205 for
Loess/Vertisol/Sand with wicking on) -- confirming the exact conflict the
prompt anticipated. Despite the wet-habitat improvement, **no soil reaches
"thriving"**: all three come out `declining`/`falling` under every condition
tested (both wicking settings, `matric_d_cut_um` in {40,50,60}, n=128 and
n=40), because an inherited, untouched bottleneck -- the per-macropore-cell
`S_cap` (Sec. above) -- caps usable substrate far below what growth needs
regardless of habitat wetness; this was verified to predate Stage 5 (the
identical decay pattern appears running the current `dual_porosity.py` under
Stage 4's old `global_theta` mode too). Stage 5's fix operates entirely
upstream of that bottleneck and so cannot, by itself, resolve it.

## Stage 6 (EXPLORATORY) -- port the sandbox fix into Variant C

`prompts/soil_respiration_prompt_stage6_portback.md`. A separate exploration
track, `pde_sandbox/` (an abstract, pure-NumPy reaction-transport search on a
toy heterogeneous grid -- its own `sandbox.py`, `search.py`, `LOGBOOK.md`,
`RESULTS.md`), proved two things there that this stage ports directly into
`dual_porosity.py` (Variant C, as left after Stage 5) -- reusing the
sandbox's own implementations, not re-deriving them. Everything else from
Stage 5 (matric-potential wetting, inscribed-circle macropore geometry, the
three phases, frozen biology/K(d)) is unchanged.

**Change 1 -- implicit transport (`dual_porosity.implicit_diffusion_step`,
`_prepare_implicit_cache`, `_implicit_neighbor_sum`).** Direct ports of
`pde_sandbox.sandbox.local_implicit_step`/`prepare_local_cache`: solves
`u - dt*div(D*grad(u)) = rhs` by Jacobi iteration
(`stage6_implicit_iterations=90` rounds, `stage6_implicit_relaxation=0.85`),
harmonic-mean face conductivities (`model.harmonic_mean`) so a dry/solid
cell (D=0) is a true barrier. Unconditionally stable regardless of D, dt, or
dx -- no more `max(D)*dt/dx^2<0.2` cap. Replaces BOTH of Stage 4/5's
separate `D_macro_t`/`D_matrix_t` explicit divergences with ONE solve over
`D_total = D_macro_t + D_matrix_t` (disjoint supports, so the sum is exactly
the per-cell D) -- this also, for the first time, lets a macropore cell and
its adjacent matrix cell diffuse into each other directly (previously only
the special-cased `k_leak` term crossed that boundary at all).

**Important correction, found before running anything**: the OLD explicit
stability cap was never actually close to binding here --
`D_macro_base.max() * 1.75 * dt / dx^2` (worst case, full habitat-cluster D
bonus) measured 0.035-0.088 across the three soils, comfortably under 0.2.
So removing the cap could not, by itself, be why substrate stayed frozen
~1-3 cells in Stage 3.5/4/5 -- that reach is the genuine physical diffusion
length `sqrt(D*T_total)/dx` given this model's D(d) magnitudes and total
simulated time, unaffected by solver choice for an already-stable dt. The
real new mechanism is the topological one above (unifying macro+matrix D
into one connected solve), verified with a standalone point-source unit test
(mass conserved to 1e-13, field visibly spreads over repeated calls) before
running any soil.

**Change 2/3 -- coating/trapped OM split (`dual_porosity.split_om_coating_
trapped`, `bfs_geodesic_distance`).** Ported from `pde_sandbox.sandbox.
build_down_source`'s geodesic construction: a soil's OM is split into a
COATING pool (porous-matrix cells within `stage6_coating_thickness=2`
geodesic steps of a macropore cell, distance measured along a WET connected
matrix path -- the destination macropore cell itself need not be wet, since
"wet" gates the path, matching the model's existing convention that
`water_mask` is a transport gate, not a growth gate) and a TRAPPED pool
(everything else in the porous matrix). Coating dissolves at the normal
`k_dis`; trapped dissolves at a new `k_dis_slow=0.005` (10x slower -- a tail,
not a permanent lock). The coating FRACTION is per soil and is the one
deliberately hypothesis-motivated, tuned input (`configs/mapping.yaml`
`stage6_coating_fraction_{sand,loess,vertisol}` = 0.80/0.50/0.30): Sand
coating-dominated (fast flush hypothesis), Vertisol distributed (progressive
recruitment hypothesis), Loess intermediate.

**Change 4 -- throttles removed.** `S_cap` (Stage 4's per-macropore-cell
substrate-holding ceiling, verified in Stage 5 to be the compounding
blocker) is removed entirely -- growth now sees `S` directly, no ceiling on
"usable" substrate. Wicking defaults to OFF (`simulate_dual_porosity(...,
wicking_enabled=False)`), matching the sandbox's static structure; `k_leak`
and matrix internal diffusion (`k_matrix_diff`) remain as supporting
pathways.

New diagnostics: `simulate_dual_porosity` returns `habitats_recruited_t`
(`dual_porosity.count_active_habitats`, per-timestep count of distinct
Stage-5 macropore regions with at least one actively-growing cell) alongside
the existing `S_habitat_mean_t`.

**Stage 6 result.** Full detail, mechanism, and honest verdict in
`docs/stage_results/STAGE6_RESULTS.md`; all figures/CSVs in
`results/stage6/`. Headline: **the ported mechanisms measurably help where a
live wet macropore network exists, but do not revive growth anywhere.**
`mean_S_habitat_final` (vs. Stage 5's wicking-off baseline) rose ~4-5x for
Loess (0.00017->0.00083) and Vertisol (0.00026->0.00098), with 58/76 and
55/70 macropore regions respectively recruited into active growth by the
run's end -- genuinely more habitat gets fed than in Stage 5. Sand showed
essentially no improvement (0.00033->0.00034, only 1 of 3 regions recruits),
because its habitat is still almost entirely the one huge ~340um dry region
Stage 5 identified, leaving barely any wet-adjacent porous matrix for the
coating construction to place OM near regardless of its high (0.8)
configured coating fraction. Despite the real improvement for two of three
soils, `mean_S_habitat` stays 2-3 orders of magnitude below `Ks=0.2`
everywhere, so growth never exceeds maintenance decay and all three soils
remain `declining`/`falling`, pinned near the trivial `R_peak~m0*B0_total=
0.4` -- none of the three target shapes (Sand rise-then-fall, Vertisol
rising, Loess low/erratic) are reproduced. Two factors outside this stage's
scope explain why: (1) the per-soil `om_fraction`-based OM totals
(inherited from Stage 4) are an order of magnitude below the historical
shared `OM_total=800`; (2) this model's actual `D(d)`/`dt`/`dx`/`T`
combination gives a genuinely short physical diffusion length regardless of
solver, unlike the sandbox's much finer relative grid resolution.

## Stage 7 -- calibration to close the substrate gap (uniform OM %)

`prompts/soil_respiration_prompt_stage7_calibration.md`. A quantitative
CALIBRATION of Variant C (as left after Stage 6), not a new mechanism --
`MAX_ITERATIONS=15` for a genuine budget-capped search over three SHARED
scalars, fully exercised (15/15 iterations spent). Stage 6's two
identified gaps close here: the per-soil OM totals (~10x below the
historical `OM_total=800`) and the short physical diffusion reach relative
to matrix->habitat distances.

**Change 1 -- uniform OM percent (`dual_porosity._variantc_om_fraction`,
`configs/mapping.yaml` `stage7_om_fraction`).** A single SHARED OM fraction
of each soil's `porous_matrix_mask` cell count replaces the retired
per-soil `stage4_variantC_om_fraction_{sand,loess,vertisol}`. `om_total`
still differs between soils, but ONLY because porous-matrix cell count is
structural (Sand ~12% of the grid, Loess/Vertisol ~44%) -- the percent
itself is identical.

**Change 2 -- reach/D-scale calibration (`stage7_D_scale`,
`simulate_dual_porosity`).** A shared post-scale multiplier on
`D_macro_base`, matrix internal diffusion (`k_matrix_diff`), and `D_film`
(Change 3) -- raises the diffusion magnitude, and hence `sqrt(D*T)` reach,
without touching `D_d_ref` (Stage 3.5's soil-discriminating SHAPE of D(d),
left alone).

**Change 3 -- water-film macropores (`D_film`, `stage7_film_K_scale`,
`build_dual_grids`, `dual_porosity.dry_macro_mask`).** A drained (dry, under
the Stage 5 matric rule) macropore cell no longer has `D=0` (fully cut off)
-- it gets a small shared film diffusivity `D_film` (also scaled by
`stage7_D_scale`) and hosts activity at reduced capacity
(`K *= stage7_film_K_scale`) instead of zero. `K(d)` itself was already
independent of wet/dry status (only `D` was ever gated by water) -- this
change only needed to scale `K` down on dry cells and give them a nonzero
`D`, not build a new habitat rule.

**Calibration search result (SS5).** Geometric escalation from a low
starting dose cleared the growth target (`S>=0.025` fed-habitat substrate,
the `f_S>~0.11` threshold restated as a direct substrate value) at
iteration 8; the remaining budget probed `D_film` alone for distinctness
without success (`min_distance` stayed 0.027-0.032 across every probe,
nowhere near the 0.3 threshold). The winning shared candidate
(`om_fraction=2.45, D_scale=12.26, D_film=0.38`) is frozen in
`configs/mapping.yaml`.

**Headline result.** Fed-habitat substrate clears `S>=0.025` for all three
soils by 2-7x (Loess 0.176, Sand 0.123, Vertisol 0.050) -- the direct gap
Stage 6 identified is closed, robustly, at both n=128 and n=40. Change 3
reveals this was never a Sand-specific problem: ~91-100% of EVERY soil's
macropore habitat sits dry and was previously fully transport-isolated
(Loess 91.7%, Sand 99.9%, Vertisol 91.3% of habitat cells). But **all three
soils bloom into the same "rising" shape** (`B_final/B0` = 21.6/5.5/20.2),
and pairwise distinctness never exceeds ~0.06 -- exactly the "everything
blooms and rises together" under-separation failure SS5/SS9 anticipated,
confirmed not to be an overshoot artifact (the calibration history and an
independent manual sweep both show the shape collapse happens as soon as
growth clears at all, at every dose tested). Mechanistically: once
fed-habitat substrate clears the growth/decay switch anywhere, unbounded
logistic growth (`r_max=1.0`, no other soil-varying growth-RATE mechanism)
ramps every soil toward carrying capacity within the run, so `R(t)`'s
qualitative shape stops being governed by the structural differences
(`coating_fraction`, macropore topology) that still visibly change HOW MUCH
and HOW MANY regions grow (Stage 6's diagnostics, still present: Loess
58/76 and Vertisol 55/70 regions recruit, Sand only 1/3 since its network
is topologically one region regardless of how much of it gets a film
pathway). Full detail in `docs/stage_results/STAGE7_RESULTS.md`; all
figures/CSVs in `results/stage7/`.

## Stage 8 -- equal total OM, depletion-driven shapes

`prompts/soil_respiration_prompt_stage8_depletion.md`. A DEFINED
configuration with verification, NOT a search: `(om_total, T)` set by at
most 3 directed diagnostic probes, no per-soil tuning. Stage 7's own
diagnosis motivates the one change: with abundant, unequal carbon, every
soil bloomed into the same rising curve because carbon never depleted --
the shape difference must come from the DEPLETION TIMING of the fuel, not
from how much arrives.

**Change 1 -- equal absolute `om_total`** (`dual_porosity._variantc_om_
total`, `configs/mapping.yaml` `stage8_om_total=600.0`). A single ABSOLUTE
total, byte-identical across all three soils, replaces Stage 7's shared
*fraction* (which still gave unequal totals since `porous_matrix_mask`
size is structural). The same total carbon is now CONCENTRATED on Sand's
small coating (few habitat regions) vs. SPREAD THIN across Vertisol's many
habitats -- exactly the lever meant to let the model's existing finite-fuel
+ Monod-limitation physics (present since Stage 6) produce different
depletion timing per soil. The coating/trapped split still follows each
soil's own geometry.

**Verification** (`simulate_dual_porosity`'s new `fuel_remaining_t`,
fraction of the initial total OM still undissolved): confirmed a real,
soil-ordered depletion at the n=128 headline -- Sand 84% consumed by `t=T`,
Loess 65%, Vertisol 51% -- exactly the predicted ordering, not imposed.

**Result.** Distinctness succeeds decisively: `min_distance=0.7955` at
n=128 (vs. Stage 7's 0.059), all three pairwise distances individually
clear the 0.3 threshold by a wide margin -- confirming the depletion-timing
hypothesis for the PRIMARY "three curves distinct" criterion. But the exact
qualitative shape match is only 1 of 3 (Loess's "erratic" matches; Sand's
real decline, R_end/R_peak=0.935, is too late/mild in the 4-period view at
n=128 to read as "falling" rather than "rising" -- it reads sharply
correct only at the smaller n=40 grid; Vertisol's ~55-70 recruited regions
also eventually run low on the SAME lean shared pool, R_end/R_peak=0.599,
reading "erratic" rather than purely "rising"). Both shortfalls trace to
one tension: a shared total lean enough for Sand's fall to register
sharply is not abundant enough for Vertisol's many regions to sustain a
monotonic rise through the whole window, and vice versa -- with no
per-soil knob available (by this stage's own design) to resolve that
tension separately. Full detail in `docs/stage_results/STAGE8_RESULTS.md`;
all figures/CSVs in `results/stage8/`.

## Stage 9 -- ranked per-soil OM + placement

`prompts/soil_respiration_prompt_stage9_ranked_om.md`. Another DEFINED
configuration with verification: Stage 8 showed a single EQUAL total OM
cannot give Sand-fall AND Vertisol-rise together (opposite regimes). Three
ranking-only changes (never fitted to the measured respiration values),
`(base scale, coating fractions)` set by <=3 directed probes:

1. **Ranked `om_total`** (`_variantc_om_total`'s Stage 9 tier,
   `stage9_om_base=600.0` x per-soil `stage9_om_multiplier_<soil>`, fixed
   ordering Vertisol(1.5) > Loess(1.0) > Sand(0.6)) -- confirmed:
   Vertisol 900 > Loess 600 > Sand 360.
2. **Vertisol's coating fraction raised** (`stage9_coating_fraction_
   vertisol`, 0.30 -> 0.60).
3. **Sand's coating fraction lowered** (`stage9_coating_fraction_sand`,
   0.80 -> 0.50), subject to the "burst wins" exception (never triggered --
   Sand's rise-then-fall stayed genuine at every candidate tested).

**Directed probes (SS4).** All 3 spent, each checked at BOTH n=40 (fast)
and n=128 (primary) -- a critical methodological finding: escalating
Vertisol's total/coating (probes 2-3) looked increasingly promising at
n=40 (Loess-vs-Vertisol distance rising 0.096->0.190) but consistently
WORSENED the same pair at n=128 (0.162->0.051). Per SS4/SS6, the n=128
grid decides adoption -- probe 1 (the prompt's own suggested starting
ratios) is the actual best n=128 candidate and is what is frozen.

**Result.** Two of three shapes are reproduced: Sand's R(t) shows a real,
sharp fall (R_end/R_peak=0.212, `metrics.describe_shape` reads it
correctly as "falling") and Vertisol's reads "rising" by the same
continuous descriptor; Loess's period-binned shape ("erratic") matches its
target as it did in Stage 8. But the coarser 4-period (P1-P4) classifier
used against the experiment's own measurement windows reads BOTH Sand and
Vertisol as "erratic" (non-monotonic across the 4 bins), so which
descriptor is applied changes how many labels "match." **Distinctness
fails** (`min_distance=0.162` at n=128, short of 0.3), driven entirely by
Loess-vs-Vertisol (0.162) while Loess-vs-Sand (1.874) and Sand-vs-Vertisol
(2.007) both clear the threshold overwhelmingly. Named cause (SS6, not
chased further): Loess and Vertisol are both large, well-connected,
many-region soils (76/70 macropore regions) whose recruitment is fast and
front-loaded under this model's shared transport magnitude (`D_scale`,
kept from Stage 7) -- their normalized R(t) shapes stay too similar to
each other regardless of how much their total OM / coating fraction
differ, and pushing the gap further hurts rather than helps at the primary
grid. Full detail in `docs/stage_results/STAGE9_RESULTS.md`; all
figures/CSVs in `results/stage9/`.

## Stage 10 (FINAL) -- Sand real burst-then-crash, Vertisol slow rise

`prompts/soil_respiration_prompt_stage10_final.md`. The last stage of the
project. Two soils only -- Vertisol and Sand (Loess dropped). Stage 9's
Sand never produced a genuine growth burst (`R_peak` sat AT the trivial
maintenance floor `m0*B0=0.4`) and its Vertisol's "rise" was really an
early spike that then dipped, both because the ONE shared, fast `D_scale`
(needed since Stage 7 just to get growth to fire at all) recruits a
soil's whole habitat almost at once. Stage 10 explicitly sets aside the
shared/ranked-scalar purity of Stages 7-9 for this final two-soil
deliverable and tunes PER-SOIL physical parameters directly (never
arbitrary curve-drawing -- the mechanism itself, finite-fuel + Monod-
limitation + implicit transport, is unchanged since Stage 6):

- **Sand** (`stage10_om_total_sand=4000`, `stage10_coating_fraction_
  sand=0.97`, `stage10_D_scale_sand=80`, `stage10_T_sand=3000`):
  concentrated, fast, finite. A large total, almost entirely on the fast
  coating, with Sand's OWN transport magnitude raised well above the old
  shared value so the burst actually arrives early -- naively raising
  `om_total` alone at the shared `D_scale` was checked directly and made
  the burst LATER and broader, not sharper (`peak_frac` 0.75 -> 0.88 as
  `om_total` rose from 3000 to 6000 at the shared `D_scale`), confirming
  `D_scale` itself, not just dose, had to move.
- **Vertisol** (`stage10_om_total_vertisol=2000`, `stage10_coating_
  fraction_vertisol=0.30`, `stage10_D_scale_vertisol=1.5`,
  `stage10_T_vertisol=1500`): distributed, slow, sustained. A large total,
  mostly in the slow trapped pool, and its own transport magnitude SLOWED
  well below the old shared value (12.2582) so its ~70 dispersed macropore
  regions ignite gradually across the whole window instead of together.

`_variantc_D_scale` (`dual_porosity.py`) is the first place `D_scale` is
NOT one shared constant across soils. `T` also became per-soil, checked
directly rather than assumed: the same Sand candidate that crashes
cleanly by `t=1500` at n=40 needs `T=3000` at n=128 to fully register the
crash (a measured, not assumed, grid-scale timescale effect -- more
habitat cells means more peak biomass, and biomass decay is a first-order
process bounded by `m0+m_s`, so a larger peak takes proportionally longer
in absolute time to decay); symmetrically, lengthening Vertisol's own
window to `T=3000` was checked and made its rise complete and dip before
`t=T` too (`peak_frac` 0.91 -> 0.46). Distinctness is computed on the two
curves resampled onto a common fraction-of-own-window time axis, since
their raw `R_t` arrays have different lengths.

**Result: every §3 pass/fail criterion is met, at BOTH n=128 and n=40 (16
of 16 checks pass).** Sand: `R_peak=15.79` (39x the maintenance floor),
77% of R at peak is growth-associated (not maintenance), `B_peak/B0=
10.71`, fuel drops from 1.0 to 0.015, `R_end/R_peak=0.019` (a real crash),
peak at 18% of the window. Vertisol: late-third mean 5.8x its early-third
mean, peak at 91% of the window, `R_end/R_peak=0.98` (no dip), `R_peak`
well under the burst threshold, 55 of 70 regions recruited, still ~49% of
its fuel unspent at `t=T`. Distinctness `=2.17`, more than 7x the 0.3
threshold. Full detail, every diagnostic number, and the honest framing of
what per-soil `D_scale`/`T` mean for the project's structure-only
invariant in `docs/stage_results/STAGE10_RESULTS.md`; all figures/CSVs in
`results/stage10/`.

## Stage 11 -- freeze Sand, fix Vertisol, optimize the loop

`prompts/soil_respiration_prompt_stage11_vertisol.md`. Continues from
Stage 10 with Sand FROZEN (no Sand parameter touched, confirmed by an
exact numeric identity check, not just convention) and work restricted to
a real performance bug and Vertisol only. n=128 ONLY.

**Performance fix.** `habitat_properties(macro_mask)` (a pure-Python
connected-component BFS over the whole grid) and the `K_t`/`habitable`/
`K_safe` fields derived from it were being recomputed EVERY timestep in
`simulate_dual_porosity`'s loop despite `macro_mask` being fixed for the
whole run -- hoisted out, computed once. `implicit_diffusion_step` gained
an optional Jacobi convergence-tolerance early stop (`stage11_implicit_
tol=1e-7`; `stage6_implicit_iterations=90` becomes a safety CAP rather
than always spent in full). A lower FIXED iteration cap (20-30, as the
prompt also suggested) was tried and rejected: it changed Sand's `R_end`
by up to 24.8%, failing "must not alter results"; the tolerance-only
approach changes `R_peak`/`R_end` by <1%. Verified speedup: Sand ~11.5x
(214s -> 18.6s), Vertisol ~50-60x (155-198s -> ~3.2s, benefiting more
since it has 70 macropore regions vs. Sand's 3 -- the retired BFS scaled
with connectivity complexity, not just grid size).

**Vertisol re-tuned, Sand untouched.** Stage 10's Vertisol already showed
a plausible R(t) slow-rise shape but its biomass only ever reached PARITY
with its initial total (`B_peak/B0=1.000`, never actual growth). Raising
`om_total` alone (2000 -> 10000, `stage11_om_total_vertisol`; coating
fraction kept at Stage 10's low 0.30, deliberately unchanged) crossed a
real threshold above which aggregate biomass growth outpaces decay
(swept directly: 2000/4000 -> 1.000, 8000 -> 1.813, 10000 -> 2.337);
`D_scale` was independently slowed further (1.5 -> 0.7,
`stage11_D_scale_vertisol`) to keep the larger-dose rise's SHAPE gradual
(late peak, no dip) rather than reintroducing an early spike.

**Result**: `R_peak=3.756` (9.4x the maintenance floor `m0*B0=0.4`),
`B_peak/B0=2.337` (genuine growth, not parity), late-third mean 11.5x its
early-third mean, peak at 95% of the window, `R_end/R_peak=0.998` (no
dip), `fuel_remaining_final=0.488` (still sustained, not exhausted).
Distinctness (time-normalized, since `T_sand=3000 != T_vertisol=1500`
remain frozen from Stage 10) `=2.00`. One limitation is reported honestly,
not engineered around: `habitats_recruited(t)` reaches its final value (55
of 70 regions) within the first 2% of the window and stays exactly flat
thereafter, at both the Stage 10 and Stage 11 Vertisol configuration --
`count_active_habitats`'s very low activation threshold (`growth>1e-9`)
is crossed almost everywhere almost immediately by the unconditionally-
stable implicit solver's per-step equilibration, regardless of `D_scale`;
R(t)'s own genuine gradual rise comes from already-"active" regions'
growth RATE slowly building up, not from new regions joining over time.
Full detail in `docs/stage_results/STAGE11_RESULTS.md`; all figures/CSVs
in `results/stage11/`.

## Stage 12 -- statistical significance via seed ensembles

`prompts/soil_respiration_prompt_stage12_significance.md`. An ANALYSIS
stage, not a modeling stage: Sand and Vertisol's Stage 11 parameters
(`om_total`, coating fraction, `D_scale`, biology, `T`) are read exactly
as frozen and never re-tuned. For each soil, M=40 random `pore.seed`
values generate fresh pore-field realizations of the SAME structural
class (everything else held fixed) -- 80 simulations total, 13.0 minutes
using the Stage-11-optimized code.

**Result**: every shape metric tested (`peak_frac`, `R_end/R_peak`,
late/early ratio, `R_peak`, `cum_co2_final`) shows COMPLETE rank
separation between the two 40-seed ensembles (Mann-Whitney rank-biserial
`r = +-1.000`, `p` in `1e-14`-`1e-15`, Cohen's `d` 5.0-36.6 -- numpy-only
implementations, no scipy). Curve-level: mean within-soil pairwise
distance is 0.158 (Sand) / 0.016 (Vertisol); mean between-soil distance is
2.063 -- 13x and 127x larger respectively; a bootstrap 95% CI on the
between-minus-within gap is `[1.963, 1.989]`. An optional bounded causal
sweep (Sand's own coating fraction, 4 values x 5 seeds) shows
`R_end/R_peak` (crash depth) drops monotonically as coating fraction
rises (0.055 -> 0.022 -> 0.006 -> 0.001), confirming the swept STRUCTURAL
FEATURE itself drives the shape metric, not merely "Sand always crashes."

**Conclusion**: because each soil's structural class is held fixed while
only the random pore-field layout varies, and every seed reliably
reproduces that soil's characteristic shape with between-soil separation
one to two orders of magnitude larger than seed-to-seed noise, each
soil's respiration shape is determined by its defined structural features
(Sand: few concentrated coating-fed macropore regions in a small, 18%
matrix; Vertisol: dozens of distributed regions in a large, 68% matrix,
both seed-invariant properties of the shared literature PSD's marginal
distribution) -- not by a lucky random realization. Full detail, every
number, and the numpy-only stats implementations in
`docs/stage_results/STAGE12_RESULTS.md`; all figures/CSVs in
`results/stage12/`.

## Stage 13 -- real measured PSD + connectivity data

`prompts/soil_respiration_prompt_stage13_real_data.md`. A DATA-SOURCING +
RECALIBRATION-CHECK stage, not a new mechanism: `dual_porosity.py`'s
physics is untouched. Replaces Stages 3.5-12's literature-parametrized
PSDs (`configs/psd_literature.yaml`, `pore.mode: psd_parametric`) with
real measured PSD + Track E connectivity data from a companion CT-scan +
nnU-Net segmentation + pore-network-topology research project, via a new
`pore.mode: psd_measured` config path (`psd_measured.py`, `configs/
psd_measured.yaml`, `data/psd_measured/<soil>/`) that returns the SAME
dict shape `psd_parametric.build_soil_psd` does, so `build_dual_grids`
needed only a mode dispatch, no other change. Real PSD tables are sourced
per soil with logged provenance (exact run name/date, see each soil's
`source.txt`): Sand reuses the Stage-3-era `rehovot_150z` table (the
current full-volume run's own bin-level table lives only on an
unreachable network share -- flagged, not silently worked around); Loess
and Vertisol use freshly-copied current full-volume real tables.

**Result: real data breaks Stage 10/11's pass/fail table completely, for
a structural reason no retuning can fix.** Sand and Vertisol's real
measured PSDs carry EXACTLY 0% recorded pore volume below `K_d_low=30um`
(the model's macro/matrix split threshold), in every one of 50 and 32
log-spaced measurement bins respectively -- the same CT-resolution
limitation `psd_parametric.py`'s own docstring already names as Stage
3.5's reason for replacing real image-derived PSDs with literature
mixtures in the first place, now independently reconfirmed with fresh,
current real data. Because ALL organic matter is placed inside
`porous_matrix_mask` (diameter < `K_d_low`), a 0% sub-30um fraction leaves
that mask empty, so `OM_total_placed=0.0` for any `om_total`. Stage 10/11's
frozen configuration, re-run unchanged on real data, collapses both soils
to `R_peak=0.400` exactly (the trivial `m0*B0` maintenance floor) and
outcome `dead` -- down from Sand's literature-PSD `R_peak=15.79`
(declining/burst-crash) and Vertisol's `R_peak=3.76` (thriving/slow-rise).
A budget-capped search (`MAX_ITERATIONS=15`) over the same 4 allowed knobs
(`om_total`, coating fraction, `D_scale`, `T`) stops after only 2 probes
per soil, each confirming `OM_total_placed=0.0` independent of the knob
swept -- an honest structural negative result, not a search failure.
Loess (never part of Stage 10/11's frozen set) retains 18.5% of its real
measured volume below 30um and is unaffected by this collapse.

Region-connectivity (SS2) is compared against real Track E numbers for all
three soils without forcing agreement: the model's emergent macropore
structure agrees qualitatively with all three real soils' own
"one dominant, highly-connected component" story (real largest-component
fraction 92.7-96.3%; model largest-region fraction 0.998-1.000), but this
is not independently earned for Sand/Vertisol (a trivial consequence of
the entire grid becoming one macropore region once the matrix phase
vanishes), and the model's single global diameter threshold has no
mechanism to reproduce real soils' simultaneous large-scale connectivity
and small-scale fragmentation, for any soil. Full detail, every number,
and the full data-provenance narrative in `docs/stage_results/
STAGE13_RESULTS.md`; all figures/CSVs in `results/stage13/`.

## Stage 14 -- investigating the `habitats_recruited(t)` saturation artifact

`prompts/soil_respiration_prompt_stage14_habitat_diagnostic.md`. A
DIAGNOSTIC stage, not a retuning stage: Stage 10/11/13's frozen calibrated
results are unchanged, confirmed by an exact (byte-identical, not merely
noise-level) identity check. `simulate_dual_porosity` gained one new,
OFF-by-default keyword, `debug_growth_distribution` -- when `True`, it
additionally returns `growth_distribution_t` (per-timestep raw `growth`
value on every macropore cell, shape `(T, n_macro_cells)`) plus
`macro_mask_debug`/`region_labels_debug` for indexing it back onto
regions; this is a pure read-only side effect, verified to change no
other output.

**Result: Stage 11's suspected equilibration mechanism is CONFIRMED
directly.** Of Vertisol's 70 macropore regions, the 55 that ever cross
`growth > 1e-9` all do so within the first 2% of the 1500-step window
(median first-active timestep = 2) -- the implicit solver's per-step
equilibration reaching almost the whole connected network almost
instantly, independent of `D_scale`, exactly as Stage 11 suspected but
never confirmed. A graded, growth-WEIGHTED "active mass" metric (raw
growth summed over all macropore cells, every step) tracks `R(t)`'s own
genuine gradual rise almost exactly (`r=0.996`), while the existing
boolean `habitats_recruited(t)` count is effectively uncorrelated with it
(`r=0.072`) -- not merely incomplete, actively misleading for that
purpose. A per-region "time to half its own eventual peak growth" metric
independently confirms real per-region growth RATE builds up gradually
across most of the window (median 44% of window), directly demonstrating
at the per-region level what Stage 11 could previously only infer from the
aggregate `R(t)` curve. Recommended as an addition alongside the existing
boolean count, not a replacement -- the boolean count remains correct for
what it actually measures (near-instantaneous first activity, a real
solver/threshold property). Full detail in `docs/stage_results/
STAGE14_RESULTS.md`; all figures/CSVs in `results/stage14/`.

## Stage 15 -- reintegrating Loess (three-soil final demonstration)

`prompts/soil_respiration_prompt_stage15_loess_reintegration.md`. Attempts
a genuine three-soil version of Stage 10/11's SS3 pass/fail demonstration,
Loess reintegrated on real measured data (Stage 13). **Deliberate, logged
deviation from a literal prompt reading**: SS2 implies Sand/Vertisol are
themselves on real data by this stage, but Stage 13 found real data
structurally breaks both soils' mechanism entirely (`OM_total_placed=0.0`)
-- re-running them on real data here would collapse both to the identical
trivial floor, defeating the point of a three-way comparison. Sand and
Vertisol are therefore kept at Stage 10/11's frozen, WORKING
literature-PSD configuration; only Loess uses real PSD/connectivity data,
via its own budget-capped search (`<=15` directed probes, checked at both
n=40 and n=128) over the same 4 per-soil knobs (`om_total`, coating
fraction, `D_scale`, `T`) Stage 10 established.

**Result: a genuinely mixed, honestly-reported outcome.** Three-way
distinctness (the headline success criterion) succeeds robustly at both
grids (min pairwise distance 0.442 at n=128, 0.659 at n=40, both well
above 0.3) -- reached only after 2 of the search's 5 probes directly
reproduced Stage 7/9's own named "blooming into the same shape" failure
(raising Loess's total/`D_scale` enough to clear real growth pushed its
curve toward Vertisol's own shape, `Vertisol_vs_Loess` distance falling to
0.05-0.22), resolved by a smaller, faster-clearing dose (`om_total=750`,
coating=0.5, `D_scale=9.0`, `T=1000`) rather than by loosening the
threshold. Loess's "low" target (lowest `R_peak` of the three soils)
succeeds at n=128 but FAILS at n=40 -- a genuine grid-resolution
disagreement, exactly the kind this stage's own instructions warned might
recur. Loess's "erratic" target (non-monotonic, `metrics.describe_shape`'s
`n_peaks>=2`) was tested across every probe and configuration tried
(`om_total` 150-40000, coating 0.1-1.0, `D_scale` 0.3-80, `T` 600-3000) and
NEVER achieved -- `n_peaks=1` in every case, confirmed as a structural
property of this model's smooth, monotonically-forcing mechanism (no
feedback loop capable of a genuine dip-then-recovery), independently
reconfirming Stage 3.5's own identical finding under literature PSD, now
shown to hold under real PSD data too. A byproduct finding: Vertisol's
Stage-11 configuration (never tested at n=40 before, since Stage 11
restricted itself to n=128 only) turns out to also pass every one of
Stage 11's own success criteria at n=40 -- a genuine, previously-unknown
robustness result. Full detail, the search trace, and the extended
pass/fail table in `docs/stage_results/STAGE15_RESULTS.md`; all
figures/CSVs in `results/stage15/`.

## Stage 16 -- the shared OM rule restored; structure sets fuel-habitat geometry

`prompts/soil_respiration_prompt_stage16_shared_om_rule.md`. A MODEL-WIDE
STRUCTURAL change, not a calibration stage: the per-soil coating/trapped OM
split introduced in Stage 6 is retired as the active default and ONE
shared, structure-blind OM placement rule takes its place, so that every
difference -- between soils, and between macropore regions within a soil --
must emerge from the pore field rather than from a hand-set carbon
parameter.

**Change 1 -- `om_mode` (`dual_porosity._om_mode`, `configs/mapping.yaml`).**
The new active default `shared_fine_pore` (`place_om_shared_fine_pore`)
places ALL organic matter by the model's ORIGINAL SS4 rule -- `OM ~
d^-om_fine_bias` within the porous matrix, normalized to that soil's total
-- with spatially UNIFORM dissolution at a single `k_dis`. No coating pool,
no trapped pool, no per-soil `coating_fraction` (retired; reported as
`nan`), and no `k_dis_slow`: the slow pool is identically zero and
`simulate_dual_porosity` additionally ties `k_dis_slow` to `k_dis` under
this mode, enforcing the stage's invariant that **no term anywhere in the
model has a rate that depends on distance to habitat**. Distance enters
only through diffusion of the dissolved substrate, which the model already
does. The Stage 6-15 path stays reachable behind `om_mode: coating_trapped`
(the same preservation Stage 13 gave `psd_parametric` alongside
`psd_measured`); every Stage 4-15 runner pins it explicitly on load, and it
still reproduces Stage 11 to within 0.004%. A third mode,
`shared_lognormal_clumps`, exists only for SS7's clumping question.

**Change 2 -- fuel-habitat geometry as a measured quantity**
(`substantial_om_mask`, `analyze_fuel_habitat_distance`): the distribution,
over macropore regions, of geodesic distance to the nearest SUBSTANTIAL OM
("substantial" = matrix OM at or above the median of nonzero matrix OM, a
scale-free and therefore dose-invariant definition). This is the emergent
structural quantity the stage's argument rests on.

**Change 3 -- observational outputs** (`active_mass_t` always,
`region_biomass_t` opt-in via `track_region_biomass`): Stage 14's
recommended growth-weighted metric, and per-region peak biomass for
Portell et al. (2018)'s hotspot metric. Both verified byte-identical
against non-instrumented runs, exactly as Stage 14 verified its own.

**Change 4 -- Stage 14's causal attribution CORRECTED.** Stage 14 blamed
Vertisol's near-instant habitat activation on the implicit solver's
per-step equilibration. Direct test refutes this: with transport switched
off entirely (`D_scale=1e-9`) regions still activate within 2-4 timesteps
and recruitment reaches its final count SOONER (0.27% of window) than with
transport on (1.33%). The cause was `split_om_coating_trapped` itself --
the coating shell fuels all ~70 regions simultaneously by construction. A
dated correction is appended to `STAGE14_RESULTS.md`; its diagnostic
findings stand unchanged.

**Result: shapes and distinctness survive the change; recruitment timing
fails structurally.** Both soils keep their Stage 10/11 target shapes at
both grids under the shared rule (Sand 6/6 checks, Vertisol 5/5, n=128 and
n=40), and Sand-vs-Vertisol distinctness is essentially unchanged from
Stage 15's hand-set configuration (2.082 vs 2.001 at n=128; 2.152 vs 2.217
at n=40) -- so the per-soil carbon knob was buying nothing structure did
not already provide. The SS3 dose ladder used 9 of 15 allowed probes and
**adopted no rung**: rung 2's shared concentration per unit matrix volume
accounts for only 22.2% (linear) / 46.5% (log) of the required per-soil
dose difference, and the analogous structure-derived `D_scale` (per-soil
macropore fraction) only 2.2% / 19.6% of the transport difference. SS4(b),
recruitment timing, FAILS at every rung and is dose-independent
(`frac@2%=1.000` across `om_total` 3,000-50,000 and `D_scale` 0.7-80),
because the emergent fuel-habitat distances the shared rule produces are
**1-2 grid cells for both soils** (Sand median 10 um, Vertisol 20 um) --
66 of Vertisol's 70 regions physically touch a fuel-bearing cell, so the
boolean `growth > 1e-9` activation test is saturated by the pore field's
own geometry. Two of the prompt's premises were corrected against
measurement: its SS5 evidence table was produced by a configuration that
still punches a geodesic exclusion zone around every region (which is the
entire source of the staggering it reports), and its "79% selectivity"
baseline is the boolean recruited count, not Portell's peak-biomass metric
(on which Stage 11's baseline is 14.3% and the shared rule gives 21.4%,
against Portell's 9.5 +/- 4.0%). Per SS6 this clean negative is the signal
to attempt a physiological rather than geometric fix (dormant/active
biomass with substrate-dependent reactivation), explicitly out of scope
here. Full detail in `docs/stage_results/STAGE16_RESULTS.md`; all
figures/CSVs in `results/stage16/`.
