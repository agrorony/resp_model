# Stage 18 RESULTS -- uniform OM per cell on real measured PSDs: FAILS at step 1

`uniform_om_model.py`, `configs/mapping_stage18.yaml`, `stage18_run.py`,
`results/stage18/`. A clean single-compartment redesign of OM placement,
specified in conversation with the user. It replaces both the v3 branch's
synthetic-fabric approach and Stage 17's enzyme experiment; neither is
extended. (Numbered 18 because 17 is taken by the open Stage 17 PR.)

## 0. Verdict

**FAIL at step 1 (distinctness, single seed). Per the procedure, step 2
(seed robustness, shuffle control, mass conservation battery) was NOT run,
and no parameter of D(d), K(d) or the totals was tuned.**

| grid | Sand-Vertisol | Sand-Loess | Vertisol-Loess | threshold | soils non-trivial |
|---|---|---|---|---|---|
| **n=128 (configured, verdict)** | **0.018** | **0.0004** | **0.017** | 0.3 | **none** (B_final/B0 = 0.64 / 0.60 / 0.64) |
| n=40 (reported) | 0.046 | 0.002 | 0.043 | 0.3 | all three (2.76 / 2.69 / 2.76) |

The three real soils produce essentially the same normalized respiration
curve. The closest pair (Sand-Loess) is ~750x below the bar at n=128. At
n=128 every soil also fails non-triviality: biomass declines and R(t) is
still creeping upward at the end of the window (peak at 1.00, R_peak
1.22-1.34).

## 1. The architecture as built (exactly the specified rule set)

1. **OM(i,j) = OM_total / N for every cell**, every soil. No pore-size
   weighting. Each cell's diameter is drawn from the soil's real measured
   PSD (`pore.mode: psd_measured`, `configs/soil_<s>_measured.yaml`), so
   uniform-per-cell OM is pore-volume-proportional OM.
2. **D(d) is the only pore-size-dependent transport law**:
   `D_min + (D_max - D_min)*(d/(d+D_d_ref))^p`, main's constants (0.02, 3.0,
   250 um, 2). Smooth, continuous, >= D_min. Water gating by the unchanged
   saturation physics still sets D = 0 in air-filled cells.
3. **K(d): smooth lognormal hump, floored at k_min** (never zero, no cutoff).
4. **Equal totals, enforced and asserted every run**: sum(OM) = 800,
   sum(B0) = 20, sum(K) = K_mean x N, identical for every soil. No per-soil
   `om_total`, `D_scale` or `T`.
5. **Real measured PSDs** for Sand, Vertisol and Loess (`data/psd_measured/`).
6. **Saturation unchanged**: `global_theta`, theta = 0.6, as in the measured
   configs (`model.water_mask_and_theta`, small pores fill first).
7. **No enzyme, no macro/matrix compartments, no categorical pore classes.**
   `om_fine_bias`, `K_window`/`K_d_low`/`K_d_high` and every dual-porosity
   construct are absent from this code path and from its mapping file.
   They remain in `main`'s legacy modules only so that Stages 1-16 still
   reproduce; nothing here imports them.

Biology (`configs/biology.yaml`) and the per-cell equations are `main`'s
`model.simulate`, unchanged. The grid, T = 600 and lambda = 2.0 come from the
measured configs. Each soil keeps its configured seed.

### Two choices the spec left open, fixed before running and stated here

- **K(d)'s constants had to be put in diameter units.** `K_of_r`'s own
  constants (r_opt = 5, w_K = 0.45) are v2's arbitrary *radius* units. On
  measured *diameters* of 30-1000 um they put every cell more than 6 sigma
  off-peak, i.e. K = k_min everywhere. After equal-total renormalization
  that would be a uniform K, and the hump would do nothing. The same
  functional form was therefore rescaled to diameters, the way Stage 3
  rescaled D (`r_ref` -> `D_d_ref`), using constants `main` already had
  rather than new numbers: centre d_opt = sqrt(30 x 150) = 67.1 um and
  width w = ln(150/30)/2 = 0.805, so +-1 sigma spans `main`'s microbial
  window [30, 150] um. These are used only as the hump's centre and width;
  there is no cutoff. K_max = 1.0 and k_min = 0.005 are reused.
- **sum(K)'s shared level**: K_mean = mean over the three soils of each
  soil's own native mean K(d) at its configured seed (0.736, 0.513, 0.752
  -> 0.6673). This is a stated rule, asserted by `stage18_run.py`, not tuned.

## 2. Step 1 results in full (`stage18_summary.json`, `stage18_step1_curves.png`)

| soil | grid | non-trivial | R_peak | peak @ | R_end/R_peak | B_final/B0 | shape |
|---|---|---|---|---|---|---|---|
| Sand | 128 | **no** | 1.344 | 1.00 | 1.000 | 0.642 | rising |
| Vertisol | 128 | **no** | 1.219 | 1.00 | 1.000 | 0.598 | rising |
| Loess | 128 | **no** | 1.340 | 1.00 | 1.000 | 0.641 | rising |
| Sand | 40 | yes | 10.930 | 0.35 | 0.303 | 2.755 | falling |
| Vertisol | 40 | yes | 8.672 | 0.37 | 0.374 | 2.690 | flat |
| Loess | 40 | yes | 10.888 | 0.35 | 0.305 | 2.763 | falling |

The shape labels at n=40 differ (falling / flat / falling), but the
normalized curves the distance metric compares are nearly identical: the
largest pairwise distance is 0.046. A label difference here is a threshold
effect of `describe_shape`, not a separation.

## 3. Why it fails (diagnostics; nothing was changed, `stage18_diagnostics.txt`)

**(a) Transport is inert under uniform OM.** Switching diffusion off
entirely (D = 0 in every cell) changes R(t) by at most **1.0%** of R_peak
(Vertisol; Sand 0.09%, Loess 0.27%). Every pairwise distance moves by less
than 0.002:

| pair | n=128, D on | n=128, D = 0 | n=40, D on | n=40, D = 0 |
|---|---|---|---|---|
| Sand-Vertisol | 0.0178 | 0.0196 | 0.0460 | 0.0457 |
| Sand-Loess | 0.0004 | 0.0008 | 0.0020 | 0.0021 |
| Vertisol-Loess | 0.0174 | 0.0187 | 0.0425 | 0.0418 |

This is structural, not a calibration miss. When every cell holds the same
OM and dissolves it locally at the same `k_dis`, every cell is fed by its
own source. There is no separation between fuel and habitat for diffusion
to bridge, and no source-driven concentration gradient for D to act on.
D(d) really does differ between the soils (mean wet-cell D: Sand 0.247,
Vertisol 0.302, Loess 0.080, a 4x spread), but the model never needs to
move substrate anywhere. **The one mechanism this architecture allows to
make structure matter has nothing to transport.**

**(b) What little separation exists comes from K heterogeneity, not
transport.** With sum(K) equal, the soils differ only in how unevenly K is
spread: the coefficient of variation of K across cells is 0.71 for
Vertisol versus 0.29 (Sand) and 0.32 (Loess). Uneven local capacity makes
local growth curves differ slightly, which explains why Vertisol is the
(slightly) odd one out in both pairs it is in, and why Sand-Loess, whose K
spreads are almost the same, is the closest pair (0.0004).

**(c) n=128 starves for a supply reason shared by all soils.** `main`'s
totals (OM 800, B0 20) are absolute, not per volume. At n=128 each cell
holds OM 0.049 and releases 0.0024 per time unit; mean S at the end is
0.033, so f(S) = S/(S + Ks) is only about 0.14. Starvation death
(m_s(1 - f)B) outpaces growth, and biomass declines in all three soils
alike. At n=40 the same totals give 10x more per cell, and all three grow.
This is a property of the shared totals and the grid, identical for every
soil, so it cannot be a source of distinctness either way.

## 4. What this establishes

- Under a single compartment with uniform OM and spatially uniform
  first-order release, real measured PSDs do **not** produce distinct
  respiration dynamics: pairwise distances are 0.0004-0.046 across both
  grids, against a 0.3 bar.
- The failure is not a tuning problem. D(d) could be made arbitrarily
  more pore-size-sensitive, and the result would barely move, because
  removing transport altogether changes R(t) by at most 1%. Tuning D(d) was
  therefore not attempted, as instructed.
- For structure to act through D, some quantity has to be spatially
  separated from where it is consumed. In this rule set nothing is: OM,
  its release and its consumers share every cell. Any next step has to
  decide what, if anything, is allowed to create that separation. This
  stage does not propose or test one, in keeping with the instruction not
  to re-add discarded mechanisms.

## 5. Not run

Step 2 (3 seed sets, shuffle control, mass-conservation battery) was not
run, because step 1 failed at the verdict grid. One related fact is
already measured: the per-step transport mass error,
max |sum div(D grad S)|, is reported per run in `stage18_summary.json`.
It is the conservative-flux telescoping sum and is at round-off level. It
is not a substitute for the step-2 battery.

## 6. Deliverables

- `uniform_om_model.py`: the Stage 18 single-compartment model. It imports
  only `model.conservative_divergence`, `model.water_mask_and_theta`,
  `pore_field.generate_diameter_field_from_psd`, `pore_field.b0_field` and
  `psd_measured`; none of the discarded mechanisms is reachable from it.
- `configs/mapping_stage18.yaml`: the complete constant set for this
  architecture, with the removed keys listed as absent.
- `stage18_run.py`: the step-1 / step-2 procedure with stop-early logic.
- `results/stage18/`: `stage18_summary.json`, `stage18_step1_curves.png`,
  `stage18_step1_timeseries.csv`, `stage18_diagnostics.txt`.
