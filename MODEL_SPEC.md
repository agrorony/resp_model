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

---

# v3 -- structure-driven respiration with equal amounts

## 13. Equal-totals rules and success criterion

**Problem addressed.** In v2, distinct curves were partly bought with
*amounts*: total habitat sum(K) differed by up to 1.6x between soils (K(r)
was never normalized), and soil C held 36% of all OM in one cell (lognormal
tail x `r^-2.5` weighting). v3 forbids both.

Hard rules (enforced in code):

- **R1** `configs/biology.yaml` byte-identical (`model.check_biology_frozen`, SHA-256).
- **R2** equal totals *per unit volume* for every soil, asserted in
  `model.build_grids`:
  `sum(OM) = OM_density*N`, `sum(B0) = B0_density*N`, **`sum(K) = K_density*N`**
  (`OM_density=0.5, B0_density=0.0125, K_density=0.3` in `mapping.yaml`;
  the first two equal the v2 40x40 totals 800/20). K(r) keeps its hump shape
  -- *where* habitat is good -- and is renormalized to the shared total, the
  same rule OM and B0 already followed. Grid, dx, T, saturation mode and all
  mapping constants are shared.
- **R3** no OM hoarding: `max(OM) <= 30 x mean(OM)` (`model.r3_ok`); a
  candidate that violates it is rejected, never rescaled.
- **R4** per-soil freedom = the pore field only: texture (`mu`, `sigma`) and
  arrangement (`structure`, `lambda`, `aniso`, `aggregate`, archetype params, `seed`).
- **R5** mechanics changes only with a written literature rationale, applied
  to all soils, constants fixed a priori, mass conservation re-checked, and
  the +-50% sensitivity test (S6).
- **R6** distance metric, 0.3 threshold and non-trivial test unchanged (SS10).
- **R7** `MAX_ITERATIONS = 12` per path, everything logged.
- **R8** stability: `max(D)*dt/dx^2 < 1/(2*ndim)`.

Success (all required): **S1** R1-R3 hold; **S2** min pairwise distance > 0.3
and all non-trivial; **S3a** S2 holds for three seed sets; **S3b** mean
seed-to-seed distance < 0.5 x min between-soil distance; **S4** shuffling
each soil's pore field (same pore-size multiset, arrangement destroyed)
removes >= 0.1 of at least one pairwise distance -- i.e. arrangement, not
only texture, carries the result; **S5** exact mass conservation; **S6**
(only with a new mechanism) +-50% sensitivity.

Paths (stop at first success): P0 re-baseline v2 soils; P1 2D search,
current mechanics; P2 3D grid (percolation thresholds differ: ~0.59 in 2D vs
~0.31 in 3D); P3 enzyme-mediated depolymerization (SS14); P4 archetypes.

## 14. Code extensions (available to every soil)

- **n-D grids**: `grid.shape: [nx, ny]` or `[nx, ny, nz]` (legacy `grid.n`
  still read). `conservative_divergence` loops over axes; face conductivities
  are precomputed once.
- **Fabric archetypes** (`pore.structure`): `gaussian` (`lambda`, optional
  `aniso` stretching axis 0), `peds` (Voronoi aggregates: macropores on ped
  faces, micropores inside), `biopores` (tubular, preferentially axis-0
  channels over a correlated matrix), `hierarchical` (fine texture nested in
  large domains). All map the latent field through `r = exp(mu + sigma*z)`.
- **Enzyme option** (`mapping.enzyme.enabled`, off in the final result):
  `dE/dt = div(f_E*D*grad E) + a_E*B - d_E*E`, release
  `k_dis*OM*E/(E+K_E)` -- microbial access controls OM turnover (Schmidt et
  al. 2011; Dungait et al. 2012; Allison 2005). Not needed for the v3 result.

## 15. v3 result and interpretation

Reached at path **P1** (2D 40x40, T=600, fully wet, unchanged mechanics),
iteration 1 of the `scan.py` library scan (120 candidates):

| soil | mu | sigma | fabric |
|---|---|---|---|
| A | 0.6 | 0.3 | gaussian, lambda=8 (large domains) |
| B | 0.6 | 0.7 | gaussian, lambda=1.5, aniso=5 (elongated along axis 0) |
| C | 2.3 | 0.7 | gaussian, lambda=8 |

Distances AB=0.711, AC=0.718, BC=1.534; seed sets min 0.711/0.652/0.612;
seed noise 0.030; shuffled-field distances 0.18/0.15/0.50 (drops
0.53/0.57/1.03). Every soil has sum(OM)=800, sum(B0)=20, sum(K)=480; max cell
OM/mean = 3.5/20.0/12.7. The same triple also passes S1-S5 in `retention` mode.

**Why it works.** With first-order release, the OM->S supply curve is
identical in every soil (22% of OM remains at t=30 in all three). Structure
acts only on *delivery*, but under equal amounts delivery differs enough:
C (coarse, a third of its OM already inside habitat) bursts and declines; A
(large fine and mid-pore domains ~12 cells apart) peaks mid-window; B has OM
nearer to habitat but behind fine, elongated low-D zones -- with harmonic-mean
faces the lowest-conductance cells on a path are the bottleneck -- so it
rises slowly all window. Shuffling erases the bottlenecks and the three
cumulative-CO2 totals converge (176/158/181 vs 108/47/174).

**Limit found (strict tier).** If all three soils must also share the SAME
pore-size histogram, the best triple reaches only ~0.09 here (and ~0.2 in
exploratory 3D/retention runs). Reason: an identical supply curve passed
through different transport "filters" yields delayed/smoothed copies of one
shape. Breaking that limit requires structure to change *how much* OM is
turned over, i.e. microbial-access-controlled depolymerization (SS14, P3),
which is the documented next step.
