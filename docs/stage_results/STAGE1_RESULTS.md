# Stage 1 Results — Pore-Based Soil Environment (v2)

This document summarizes the results of **Stage 1** of the model
(`soil_respiration_prompt_v2.md`): a fully saturated pore field (every cell
conducts), from which K and D are derived functionally, with no hand
placement of organic matter (OM) or biomass.

## 1. The idea

Unlike v1, where each soil was built from hand-picked scalars (a uniform D,
"colonies" and OM at fixed distances), here every grid cell has a pore
radius `r(i,j)`, from which the following are derived:

```
K(r) = k_min + (K_max - k_min) * exp( -0.5 * ((ln r - ln r_opt) / w_K)^2 )
D(r) = D_min + (D_max - D_min) * (r / (r + r_ref))^p
```

All seven constants (`r_opt, w_K, K_max, k_min, D_min, D_max, r_ref, p`) are
fixed and identical for all three soils (`configs/mapping.yaml`) — the only
difference between soils is **the pore field itself**: `mu` (mean log
radius), `sigma` (spread), `lambda` (spatial correlation length), and
`aggregate` (whether the field is thresholded into "peds" separated by
fine-pore boundaries).

Initial biomass and organic matter are also derived from the field, by one
rule shared across all soils:

```
B0(i,j) = B0_total * K(i,j) / sum(K)         # settles wherever conditions are good
OM(i,j) = OM_total * r(i,j)^-om_fine_bias / sum(...)   # protected in fine pores
```

## 2. Success criterion — "emergent distinctness"

No specific shape (rising/falling/flat) was targeted for any soil. Instead:

1. Each `R(t)` curve is normalized by its own peak.
2. A pairwise distance is computed between every two soils:
   `dist = (1 - pearson_corr) + 0.5 * mean(|normalized difference|)`
3. **Success** = all three pairwise distances `> 0.3`, and every soil is
   "non-trivial" (biomass actually grows, and total respiration is
   positive).

## 3. Result achieved

The run documented in `LOGBOOK.md` (Stage 1, iteration 1) achieved **full
success**:

| Soil pair | Distance (threshold = 0.3) |
|---|---|
| A – B | 0.826 |
| A – C | 0.333 |
| B – C | 1.454 |

All three soils are "non-trivial" (B_final > B0, total R(t) positive).

### Pore parameters obtained per soil

| Soil | mu | sigma | lambda | aggregate | Texture description |
|---|---|---|---|---|---|
| A | 0.6 | 0.5 | 4.0 | no | relatively fine, medium-high spatial correlation |
| B | 2.3 | 0.4 | 2.0 | no | coarse, "salt and pepper" (short correlation) |
| C | 1.2 | 0.9 | 1.0 | yes | intermediate, aggregated (peds separated by fine boundaries) |

### Shapes obtained (auto-described, not targeted)

- **Soil A** (fine, correlated): slow rise across the whole window —
  `shape=rising`, peak ~2.42 around t≈15, then nearly flat.
- **Soil B** (coarse, scattered): sharp early burst then decay —
  `shape=falling`, peak ~8.77 around t≈7, dropping to ~3.5 by the end of
  the window. Its cumulative CO2 (~181) is significantly higher than the
  other two.
- **Soil C** (aggregated): slight dip at first, then a rise —
  `shape=rising`, peak ~2.58 toward the end of the window, because its
  fractured (aggregated) network of high-K cells accumulates substrate
  gradually.

Corresponding figures: `results/respiration_curves.png`,
`results/cumulative_co2.png`, and the structure maps
`results/structure_maps_{A,B,C}.png` (pore radius, K, D for each soil).

## 4. Mass-conservation check

`model.conservative_divergence` uses harmonic-mean face conductivities in a
conservative finite-volume form, with zero-flux (Neumann) boundaries.
Verified numerically: with a random D field (including fully disconnected
cells, D=0) and a random S field, `mass_after - mass_before = 0.0` exactly —
i.e. the diffusion term alone never creates or destroys mass, even with
sharp D contrasts.

## 5. An important insight that came up during calibration (documented in LOGBOOK.md)

A first attempt with a broad K distribution (`w_K=0.6`) and a mild fine-pore
weighting for OM showed that every cell already contained a bit of biomass
*and* a bit of OM from t=0 — so no cell actually depended on transport
(diffusion) to start respiring, and all three soils converged to nearly
identical curves (distances ~0.02–0.09) regardless of pore structure.
Narrowing the K bell curve (`w_K=0.45`, a genuinely selective habitat) and
sharpening the OM fine-pore bias (`om_fine_bias=2.5`) restored a real
structural lever: whether a soil's high-K regions (close to `r_opt`) sit
near or far from its OM-rich regions (the finest pores) — which depends on
`mu/sigma/lambda` — materially affects how fast substrate reaches biomass.

## 6. Caveats and limitations of this stage

- **Stage 1 = fully wet** — there is no percolation gate yet (every cell
  conducts); stages 2–3 will add saturation (theta) and the connectivity
  question.
- The mapping constants (`w_K`, `om_fine_bias`, the saturating form of D)
  were calibrated once, before the per-soil search, so that pore structure
  would have any causal effect on respiration at all — this is model
  design, not per-soil tuning.
- The starting points for `(mu, sigma, lambda, aggregate)` were found via a
  small grid scan outside the iterative nudge loop (`search_loop.py`),
  which is fully implemented and will keep running if a future stage
  doesn't meet the threshold from the start.
- Fully deterministic — the only randomness is the fixed seed used to
  generate each soil's pore field.

## 7. Relevant files

- `pore_field.py`, `configs/mapping.yaml` — pore field generation and K/D mappings
- `model.py` — mass-conserving diffusion with spatially varying D
- `metrics.py`, `search_loop.py` — the emergent-distinctness criterion and the budget-capped search
- `MODEL_SPEC.md` — full detail of the equations and constants
- `LOGBOOK.md` — full log of the search iterations and constant calibration
- `results/` — all figures and raw data (CSV) from Stage 1
