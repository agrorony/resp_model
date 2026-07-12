# Claude Code task v2: pore-based soil environment (staged)

You are extending the EXISTING project (`model.py`, `configs/`, `run.py`,
`metrics.py`, `search_loop.py`, `MODEL_SPEC.md`, `LOGBOOK.md`, `results/`).
Do not start from scratch; evolve what is there and keep it runnable at every step.

## Why this v2 exists (read first)

In v1 the three soils differed only in *when* substrate reached a hand-placed
colony (a hand-picked scalar D + hand-placed OM at chosen distances). That is a
timing trick, not soil. v2 replaces the hand-set knobs with a real **pore-size
field**: every cell has a pore size, and BOTH the local carrying capacity K and
the local diffusion rate D are *derived from that pore size*. The three soils
become genuinely different porous media, and their respiration differences must
**emerge**, not be scripted.

## Invariants (unchanged from v1)

1. **Biology is frozen.** The constants in `configs/biology.yaml` are shared by
   all soils and must NEVER be edited by any loop. Only structure may differ.
2. **Structure = the pore field now.** The only per-soil freedom is the
   pore-size distribution and its spatial arrangement (§A). K and D are FUNCTIONS
   of the pore field, not independent dials.
3. **Deterministic.** Randomness only via a fixed seed used to generate the pore
   field. Grid <= 50x50, seconds per run.
4. **Budget cap.** `MAX_ITERATIONS = 12` per stage. Stop and document at the cap.

## New success criterion — EMERGENT DISTINCTNESS (this replaces the v1 targets)

Do NOT aim for "A rising, B falling, C flat." Instead:

- Run the three soils, get each `R(t)`.
- Normalize each curve (e.g. divide by its own max) and compute pairwise
  dissimilarity between all three: use `dist = 1 - pearson_corr(norm_Ri, norm_Rj)`
  (optionally add a normalized L2 term). 
- **Success = all three pairwise distances exceed a threshold** (start at 0.3)
  AND every soil is non-trivial (biomass actually grows, total respiration > 0).
- Separately, *describe* each emergent shape automatically (rising / falling /
  flat / single-peak / multi-peak, via trend sign + peak count) and LOG it — but
  never steer toward a chosen label.

The loop varies only the pore-distribution parameters (§A) to make the three
soils come out distinct. It must not touch biology, and must not hand-place OM or
biomass.

---

## A. Pore field generation (the only per-soil structure)

Each `configs/soil_{A,B,C}.yaml` specifies a pore field via:

- **Texture** — the pore-size distribution: lognormal with `mu` (mean log
  radius) and `sigma` (spread). Fine soil = small `mu`; coarse = large `mu`.
- **Structure** — the spatial arrangement: a spatially-correlated Gaussian
  random field with correlation length `lambda` (small `lambda` = salt-and-pepper;
  large = big homogeneous domains). Optional `aggregate: true` thresholds the
  field into peds/clusters separated by fine-pore boundaries.

Generation recipe: draw white noise, convolve with a Gaussian kernel of width
`lambda`, standardize, then map through the lognormal CDF/quantile to get a pore
radius `r(i,j)`. Same fixed seed for reproducibility. Soils differ in
`(mu, sigma, lambda, aggregate)` — nothing else.

## B. Pore size -> local parameters

```
K(r)  = K_max * exp( -0.5 * ((ln r - ln r_opt) / w_K)^2 )      # hump: peaks at mid pore size r_opt
D(r)  = D_min + (D_max - D_min) * (r / r_ref)^p                # rises with pore size (p ~ 2)
```

- `K(r)` peaks at `r_opt` (mid-size pores best; too small clogged, too big dry),
  with a small floor so extreme cells are near-inhospitable.
- `r_opt, w_K, K_max, D_min, D_max, r_ref, p` are FIXED model constants shared by
  all soils (put them in `biology.yaml` or a new `configs/mapping.yaml`). They are
  NOT per-soil knobs — only the pore field differs between soils, so K and D maps
  differ only because the pore fields differ.

## C. Rule-based OM and biomass (no hand placement)

- **Initial biomass** `B0(i,j) proportional to K(i,j)` (small total, equal across
  soils): the community seeds itself wherever habitat is good, not at typed-in
  centers.
- **Initial OM**: distribute the fixed total OM by a physical rule tied to pore
  size — default: OM density higher in FINE-pore regions (protected organic
  matter), independent of where biomass is. State the rule in the spec; it is the
  same rule for all soils, so OM layout differs only because pore fields differ.

## D. Variable-D diffusion (IMPORTANT — do this correctly)

D is now a spatial field, so the substrate transport must be in conservative
flux form, or mass will not be conserved:

```
dS/dt = div( D * grad(S) ) + release - uptake
```

Discretize with face-centered conductivities using the **harmonic mean** of the
two neighboring cells' D (correct for sharp contrasts and for D -> 0 in later
stages). Zero-flux (Neumann) boundaries. Re-check stability with `max(D)`.
Verify mass conservation with a quick test (total S + consumed should balance
release) and note it in the LOGBOOK.

---

## STAGES (this is the scaling plan — do them in order, checkpoint between each)

### Stage 1 — pore field -> K and D (fully wet / connected)

Implement §A–§D. All pores conduct (no water yet). Run the three soils, apply the
emergent-distinctness check, let the loop vary `(mu, sigma, lambda)` within the
budget to reach distinctness. Produce Stage-1 outputs and STOP for review.

### Stage 2 — water-filled pores + connectivity (percolation)

Add a global saturation `theta` (same for all soils). Capillary rule: small pores
fill first, so a cell is **water-filled** if `r <= r_cut(theta)`, where `r_cut` is
set so the water-filled volume fraction equals `theta`. Substrate diffuses ONLY
through water-filled cells: set `D = 0` (or tiny) in air-filled cells; use the
harmonic-mean faces so disconnected regions truly stop transport. Now whether OM
can reach the biomass is a **percolation** question set by the pore distribution
and `theta`, not a distance. Re-run, re-check distinctness, checkpoint, STOP.

### Stage 3 — per-soil saturation / water retention

Give each soil its own water content `theta` (or a simple water-retention curve
mapping matric potential -> theta from its pore distribution). This adds a second,
physically-grounded axis of difference. Re-run, re-check, checkpoint.

At every stage: biology frozen, structure-only, budget-capped, emergent
distinctness (never targeted shapes).

---

## Deliverables (add to / update the existing project)

- `pore_field.py` — pore-field generator (§A) and the `K(r)`, `D(r)` maps (§B).
- Update `model.py` — variable-D conservative diffusion (§D); rule-based B0/OM (§C).
- Update `configs/soil_*.yaml` — now `(mu, sigma, lambda, aggregate)` + saturation.
- New `configs/mapping.yaml` — the fixed pore->parameter constants (§B).
- Update `metrics.py` — the emergent-distinctness check + shape auto-describer.
- Update `search_loop.py` — budget-capped search over pore-distribution params
  toward distinctness; never touches biology; logs every iteration.
- New figures in `results/` per soil and per stage: **pore-size map, K map, D map**
  (and for stage 2+, the **water-filled map**), plus the combined
  `respiration_curves.png` and `cumulative_co2.png`.
- Update `MODEL_SPEC.md` — new equations, the pore-field model, the mapping
  functions, the emergent-distinctness criterion, and a clear per-stage section.
- Append to `LOGBOOK.md` every iteration and at each stage checkpoint: params
  tried, pairwise distances, auto-described shapes, mass-conservation check,
  and an honest note on what emerged.
- Update `report_prompt.md` so the final report covers the pore-based model, the
  three stages, and interprets *why* the emergent differences arise from pore
  structure and water connectivity.

## Order of work

1. `pore_field.py` + variable-D diffusion in `model.py`; verify mass conservation
   on one soil and get a sane respiration curve.
2. Stage 1: three soils, emergent-distinctness loop, figures, checkpoint, STOP.
3. Stage 2: water-filled + percolation; re-run, checkpoint, STOP.
4. Stage 3: per-soil saturation; re-run, checkpoint.
5. Update `MODEL_SPEC.md`, `LOGBOOK.md`, `report_prompt.md`.

Remember: biology frozen; pore field is the only structure; patterns must EMERGE
and be mutually distinct, never designed; stop at the budget and at each stage.
