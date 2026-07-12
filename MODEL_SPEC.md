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

This spec covers **Stage 1**: the pore field feeds K and D, and the soil is
fully wet (every cell conducts). Stages 2/3 (percolation, per-soil
saturation) are documented separately once implemented.

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
