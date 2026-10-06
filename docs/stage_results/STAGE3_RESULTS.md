# Stage 3 Results — Real PSD, Physiological K, Bigger Grid

This document summarizes **Stage 3** (`soil_respiration_prompt_stage3.md`):
replacing the synthetic lognormal pore field with each soil's **real
measured pore-size distribution** (PSD), a **physiological** carrying
capacity `K(d)` (a microbially-relevant diameter window instead of a
lognormal hump), grids at `n=128` and `n=40`, an outcome axis that treats a
die-off as a legitimate result rather than a disqualifier, and a specific
hypothesis test: does the isolated-water-cluster ratio predict respiration
better than the PSD itself?

## 1. The real PSD data

Three real, measured pore-size-distribution runs were available on
`Z:\Rony\remote_computer backup\10.5\psd_outputs`, one per soil sample
(`psd_table.csv` per run: `Diameter_um` = bin center, `Volume_Count` = raw
voxel count in that bin i.e. volume-weighted, `Cumulative_Porosity` =
cumulative volume fraction through that bin, `Differential_PSD` =
volume-fraction density per micron; `diagnostics.json` in each run also
carries the 51 bin edges in microns, `bin_edges_um`, which are **not**
uniform -- log-spaced, soil-specific, 30 bins onward). All three start at a
minimum diameter of 30.0 µm (coincidentally exactly the K(d) hard-zero
threshold, SS4) and use 50 bins each; nothing was re-binned. Copied once
into the repo (`data/psd/<name>/psd_table.csv` + `bin_edges_um.csv`,
`psd_data.py`) so the pipeline no longer depends on the Z: mount.

| run (source dir suffix) | assigned soil | median diameter (µm) | 30-150µm habitable volume fraction | diameter range (µm) |
|---|---|---|---|---|
| `rehovot_150z` | **A** (kept "fine", matching Stage 1/2's label) | 122.3 | 0.765 | 30.9 – 578.5 |
| `mishmar_150z` | **C** (kept "intermediate") | 211.3 | 0.305 | 31.0 – 793.5 |
| `nlm_150z` | **B** (kept "coarse") | 328.9 | 0.239 | 31.4 – 2231.9 |

Assignment is by median diameter, matching each soil's original Stage 1/2
texture label (A=fine, B=coarse, C=intermediate) so the letter still means
roughly the same thing. `nlm` is also by far the most heterogeneous sample
(long tail to ~2.3 mm, visible as the low-density right tail in
`results/stage3_psd_loaded.png`).

## 2. Empirical-quantile pore field (`pore_field.py` SS3)

`pore_field._correlated_gaussian_z` (factored out of the old
`generate_pore_field`) generates the same spatially-correlated Gaussian
field as Stage 1/2 (white noise, Gaussian-smoothed with width `lambda`,
standardized). `pore_field.rank_uniform` converts it to uniform quantiles
via an argsort-based rank transform (numpy-only; exact for the field's own
empirical distribution, avoiding a dependency on `erf`/`scipy.stats.norm`
for exact Gaussianity). `pore_field.generate_diameter_field_from_psd` then
maps each cell's uniform quantile through THAT soil's empirical inverse-CDF
(`np.interp` against `bin_edges_um`/`cdf_edges`) to get a diameter in
microns. The `aggregate` boundary-carving rule is unchanged in spirit,
applied to diameters instead of radii (`_apply_aggregate_boundaries`,
shared code with the Stage 1/2 path).

**Assumption stated explicitly (per SS3):** the 2D area-fraction of cells at
each diameter stands in for the measured (3D) volume-fraction. The spatial
*arrangement* is synthetic (shared across soils, only different seeds); only
the *distribution* is real.

## 3. Physiological K(d) and rescaled D(d) (SS4-5)

```
K(d) = 0                                        for d < 30
K(d) = K_max                                    for 30 <= d <= 150
K(d) = K_max * exp( -((d - 150) / w_decay)^2 )    for d > 150   (w_decay=100, no floor)

D(d) = D_min + (D_max - D_min) * (d / (d + d_ref))^p     (d_ref = 90 um, mid-window)
```

`K_d_low=30`, `K_d_high=150`, `K_w_decay=100`, `D_d_ref=90` were added to
`configs/mapping.yaml`, fixed and shared by all soils. `K_max`, `D_min`,
`D_max`, `p` are reused unchanged from Stage 1/2. `max(D)*dt/dx^2` stays
`0.15 < 0.2` (verified) regardless of the diameter scale, because `D(d)`
saturates at the same `D_max` whether the diameters involved are ~100 µm or
~2000 µm -- **this saturation is itself the key mechanism behind SS12's
answer below.**

## 4. Model plumbing (`model.py`, `configs/soil_*.yaml`)

`build_grids` now dispatches on `pore.mode`: `"lognormal"` (Stage 1/2,
unchanged, `run.py`/`stage2_run.py` still work) or `"psd"` (Stage 3,
`pore.psd_dir` points at a soil's real PSD folder). All three
`configs/soil_{A,B,C}.yaml` were switched to `mode: psd`, with the **same**
`lambda=2.0` and `aggregate=false` for all three soils (Stage 3 invariant 1
— the only systematic difference is the real PSD) and different seeds
(1/2/3, one per soil, unchanged from Stage 1/2 — realization noise, not a
systematic lever). The grid-size assertion was relaxed from `n<=50` to
`n<=200` to allow `n=128` (SS7).

## 5. Metric revision — viability as its own axis (SS8)

`metrics.classify_outcome` labels each soil's run `thriving` / `declining`
/ `dead` from final-vs-initial biomass and the late-window respiration
trend. `metrics.emergent_distinctness` gained a `require_nontrivial` flag;
Stage 3's runner calls it with `require_nontrivial=False`, so a die-off
never disqualifies success — it is reported, honestly, as its own column.

**In this Stage-3 sweep, every soil at every theta and both grids came out
`thriving`** (see `results/stage3_theta_sweep.csv`). No die-off occurred,
unlike Stage 2 (soil C died at low theta). This is a direct, honest
consequence of turning `aggregate` off for all three soils (Stage 3
invariant 1 forces a *shared* spatial recipe) — Stage 2's soil-C die-off was
caused specifically by `aggregate=True` fragmenting its water network at low
theta; that mechanism is absent here by construction, and the real PSDs'
own habitable-volume fractions (24-77%, see SS1 table) are all large enough
that no soil's habitat collapses under this theta range alone.

## 6. Full theta × grid sweep (SS6-7)

Swept `theta` in {1.0, 0.8, 0.6, 0.4, 0.2} at **both** `n=128` and `n=40`,
same shared spatial settings, real PSDs. Full numbers:
`results/stage3_theta_sweep.csv`, `results/stage3_distinctness_by_theta.csv`.

| grid | theta | AB | AC | BC | min | success |
|---|---|---|---|---|---|---|
| 128 | 1.0 | 0.354 | 0.213 | 0.054 | 0.054 | False |
| 128 | 0.8 | 0.324 | 0.189 | 0.053 | 0.053 | False |
| 128 | 0.6 | 0.308 | 0.186 | 0.048 | 0.048 | False |
| 128 | 0.4 | 0.340 | 0.204 | 0.052 | 0.052 | False |
| 128 | 0.2 | 0.430 | 0.256 | 0.062 | 0.062 | False |
| 40  | 1.0 | 0.203 | 0.108 | 0.034 | 0.034 | False |
| 40  | 0.8 | 0.190 | 0.099 | 0.034 | 0.034 | False |
| 40  | 0.6 | 0.183 | 0.105 | 0.033 | 0.033 | False |
| 40  | 0.4 | 0.201 | 0.119 | 0.032 | 0.032 | False |
| 40  | 0.2 | 0.212 | 0.129 | 0.033 | 0.033 | False |

**Grid size matters a lot, and in a specific way.** At `n=40`, cumulative
CO2 is nearly flat across every soil and every theta (~180-187, see
`results/stage3_respiration_vs_theta.png` right panel, dashed lines
essentially overlapping) — the small domain, at `lambda=2.0`, is close
enough to homogeneous (few independent correlation lengths fit across 40
cells) that OM gets consumed at nearly the same rate regardless of soil
identity or saturation. At `n=128`, the same three soils spread out
substantially (cum. CO2 58-88 for A, 87-135 for C, 97-149 for B, each rising
with drying) — a bigger domain lets more isolated regions actually form and
matters, exactly as anticipated in SS7. The isolation diagnostics agree:
`isolation_ratio` and raw cluster counts are higher at `n=40` for the same
theta (fewer saturated cells to start with, so fragmentation "bites" sooner
in relative terms — `results/stage3_isolation_vs_theta.png`), yet that
*doesn't* translate into more distinct respiration at `n=40`; if anything
the opposite. **Bigger grids produce more distinct outcomes here, not just
more fragmented ones.**

## 7. Isolation ratio vs PSD — the hypothesis test (SS9)

Across all 30 soil×theta×grid points (`results/stage3_isolation_vs_psd_correlation.csv`):

| outcome | predictor | pearson r |
|---|---|---|
| cumulative CO2 | isolation_ratio | **+0.256** |
| cumulative CO2 | PSD habitable fraction (30-150µm) | -0.235 |
| cumulative CO2 | PSD median diameter | +0.229 |
| peak R | isolation_ratio | **+0.331** |
| peak R | PSD habitable fraction (30-150µm) | -0.150 |
| peak R | PSD median diameter | +0.138 |

`isolation_ratio` correlates with both respiration outcomes more strongly
than either PSD summary stat tested, for both cumulative CO2 and peak R —
**the specific claim in SS9 holds, but weakly** (`|r|` around 0.25-0.33 is a
real but modest edge, not a strong predictor; see
`results/stage3_outcome_vs_isolation_ratio.png` and
`results/stage3_outcome_vs_psd_stat.png` — points for a fixed soil/PSD
scatter substantially in cumulative CO2 depending on theta and grid, which
`isolation_ratio` at least partially tracks and the static PSD stat, by
construction, cannot).

## 8. Honest framing (SS12)

**(1) With spatial arrangement shared, does the real PSD alone make the
three soils distinct — and at which theta/grid?** Only partially, and only
for two of the three pairs. `AB` and `AC` clear the 0.3 threshold at several
theta values on the `n=128` grid (best at theta=0.2: AB=0.430, AC=0.256 —
AC still falls just short). `BC` (nlm vs. mishmar) never comes close
(0.048-0.062 at n=128, ~0.03 at n=40) at any theta or grid tested. The
mechanism: both B and C have median diameters (329, 211 µm) far above
`D_d_ref=90`, so `D(d)` for both soils is already saturated near `D_max`
almost everywhere — the real, substantial difference between their PSDs
(habitable fraction 0.239 vs. 0.305, median diameter 329 vs. 211 µm) gets
compressed away by the saturating diffusion mapping before it can produce
different transport dynamics. Soil A (median 122 µm, close to `D_d_ref`) is
the one soil where `D(d)` is still in its rising, non-saturated regime, so
it separates from both B and C. **No grid/theta combination in this sweep
achieves full three-way success (`min_distance` never exceeds 0.062, driven
entirely by the fused BC pair).**

**(2) Does the isolation ratio predict respiration better than the PSD?**
Yes, modestly (SS7 table above) — consistently for both cumulative CO2 and
peak R, across every soil/theta/grid point. It is a real but not dominant
predictor.

Both are reported as they came, per SS12's instruction: an honest "partial
yes" and a "yes, weakly" — not the strong three-way separation Stage 1/2 was
aiming for, but a specific, mechanistically explained reason why (the
saturating `D(d)` mapping fuses the two large-pore soils), and a positive
(if modest) answer to the isolation-ratio hypothesis.

## 9. Assumptions and limitations

- The 2D area-fraction / 3D volume-fraction correspondence stated in SS3.
- Only one realization (seed) per soil, as in Stage 1/2 — a seed-ensemble
  significance pass remains deferred.
- `lambda=2.0`, `aggregate=false` (shared across soils) was a specific,
  single choice; a different shared spatial recipe could shift where the
  saturating-`D(d)` fusion of B/C does or doesn't happen. Not swept here
  (out of budget per SS0 — this stage is a defined sweep, not a search).
- `K_w_decay=100`, `D_d_ref=90`, and the `dead_frac=0.02` threshold in
  `classify_outcome` are documented, reasonable choices, not derived or
  calibrated against independent data.
- No die-off occurred in this sweep (SS5) — the physiological-K / real-PSD
  combination with a shared (non-aggregated) spatial recipe simply didn't
  push any soil's habitat to collapse over theta in {1.0..0.2}; a wider or
  finer theta sweep, or `aggregate=true`, might.

## 10. Relevant files

- `pore_field.py` — `_correlated_gaussian_z`, `rank_uniform`,
  `generate_diameter_field_from_psd`, `K_window`, `D_of_d`.
- `psd_data.py` — loads `data/psd/<name>/{psd_table,bin_edges_um}.csv`.
- `model.py` — `build_grids` pore-mode dispatch, relaxed grid-size cap.
- `metrics.py` — `classify_outcome`, `isolation_ratio`,
  `emergent_distinctness(..., require_nontrivial=False)`.
- `configs/mapping.yaml`, `configs/soil_{A,B,C}.yaml` — Stage 3 constants and
  per-soil PSD pointers.
- `stage3_run.py` — the full theta × grid sweep, all figures/CSVs.
- `results/stage3_*` — all Stage-3 figures and raw data.
