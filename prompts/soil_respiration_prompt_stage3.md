# Claude Code task: Stage 3 — real PSD, physiological K, bigger grid

Extend the EXISTING numpy-only project (Stages 1–2 done). Do not rewrite it —
evolve it, keep it runnable. This stage replaces the synthetic lognormal pore
field with each soil's **real measured pore-size distribution**, makes carrying
capacity **physiological** (a microbially-relevant pore-size window), enlarges the
grid, revises the success metric so a die-off counts as a real outcome, and tests
one specific hypothesis about what predicts respiration.

## 0. Budget & checkpoint

`MAX_ITERATIONS = 12` ceiling on any search. This stage is mostly a defined sweep,
not a hunt — stop and write up at the checkpoint. Don't grind.

## 1. Invariants

- `configs/biology.yaml` frozen. The pore→parameter **mapping** (`mapping.yaml`,
  incl. the new K window and D params below) is fixed and shared by all soils.
- **The only systematic difference between soils is now their real PSD.** Use the
  SAME spatial settings (correlation length / aggregation) for all three soils —
  we are explicitly testing "is the PSD alone enough to produce different
  respiration?". Keep **different seeds** per soil (do NOT build them from one
  shared field — realization noise is fine; the systematic lever is the PSD).
- Deterministic; repo stays numpy + pyyaml + matplotlib only.

## 2. Load the real PSD data

The measured PSDs live on this machine at:
`Z:\Rony\remote_computer backup\10.5\psd_outputs`

- Discover the format first (one file per soil? bin edges + fraction columns?).
  Report what you found: which file is which soil, the **bin edges in µm**, and
  whether the fraction is volume-weighted. **Pore sizes are DIAMETERS (µm).**
- Use the experiment's **own bins** — do not re-bin. Build each soil's empirical
  PSD as a binned, volume-weighted distribution over pore diameter, and its
  empirical CDF (cumulative volume fraction vs diameter).

## 3. Pore field from the empirical PSD

Replace the lognormal mapping in `pore_field.py` with **empirical-quantile
mapping**:
1. Generate a spatially-correlated Gaussian field as now (same `lambda`/aggregate
   settings for all soils; different seed per soil), standardized to N(0,1).
2. Convert to uniform quantiles, then map each cell through THAT soil's empirical
   PSD inverse-CDF to get a pore **diameter** in µm (interpolate within the
   measured bins). The field's marginal diameter distribution then matches the
   measured PSD, binned to the experiment's bins.
- Assumption to state: the 2D area-fraction of cells at each diameter stands in
  for the measured (3D) volume-fraction. The spatial *arrangement* is synthetic
  (shared across soils); only the *distribution* is real.

## 4. Physiological carrying capacity K(d)

Replace the lognormal hump with a microbially-relevant window (d = diameter, µm):
```
K(d) = 0                                      for d < 30
K(d) = K_max                                  for 30 <= d <= 150
K(d) = K_max * exp( -((d - 150) / w_decay)^2 )  for d > 150   # smooth decay to ~0, no floor
```
- Hard zero below 30 µm (too tight for microbial activity), full plateau in the
  30–150 µm ideal window, smooth Gaussian decay above 150 µm toward ~0.
- `K_max`, `w_decay` (a fixed width, e.g. ~100 µm) and the 30/150 thresholds live
  in `mapping.yaml`, fixed and shared across soils. Note: cells with K=0 are
  non-habitat (biomass can't live there) but still hold water and conduct — they
  are conduits, not homes.

## 5. Diffusion mapping D(d)

Keep D rising with pore size, rescaled to microns and kept numerically stable:
```
D(d) = D_min + (D_max - D_min) * (d / (d + d_ref))^p
```
Set `d_ref` to a physically sensible diameter (e.g. ~90 µm, mid-window) and keep
`D_max` such that `max(D)*dt/dx^2 < 0.2` (the existing assertion must still pass).
Fixed and shared across soils. (`dx` stays a nominal lattice spacing; diameters
drive K and D, not the lattice geometry.)

## 6. Saturation — keep the equal-θ sweep (matches the experiment)

Unchanged from Stage 2: `global_theta`, small pores fill first (a cell is
water-filled if its diameter `d <= d_cut`, where `d_cut` is that soil's own
θ-quantile of d), air cells get D=0, harmonic-mean faces. Sweep the same θ grid
(1.0, 0.8, 0.6, 0.4, 0.2). Do NOT switch to per-soil suction — the experiment
imposed a common water content.

## 7. Grid size — run BOTH 128 and 40

Run the full sweep at **grid n=128** and at **n=40**, and compare. The larger grid
should let more isolated regions form and disperse; report how grid size changes
the isolation diagnostics (§9) and the respiration outcomes. Keep T and dx as in
Stage 1/2 unless stability forces a change (document if so).

## 8. Metric revision (approved) — a die-off is a distinct outcome, not a failure

- Distinctness (pairwise curve distance) is unchanged.
- **Remove the "all soils must be non-trivial" requirement from success.** Instead,
  classify each soil's outcome as `thriving / declining / dead` (via final-vs-
  initial biomass and the respiration trend), and treat a structurally-caused
  die-off as a *legitimate, strongly distinct* pattern. Report viability as its
  own axis alongside distinctness — never as a disqualifier. (Stage 2's soil-C
  death was arguably the most distinct result and should have counted.)

## 9. Hypothesis test — isolation ratio vs PSD

Test the specific claim: **respiration outcome tracks the isolated-regions /
saturated-area ratio better than it tracks the PSD.**
- Per soil, per θ, per grid: count the **number of isolated water clusters**
  (connected components of the water-filled mask) and compute
  `isolation_ratio = n_isolated_clusters / total_saturated_cells` (also report the
  raw counts and largest-cluster fraction already in `metrics`).
- Pick a respiration outcome (cumulative CO2, and also peak R) and, across all
  soil×θ×grid points, compute the correlation of the outcome with (a) the
  isolation ratio and (b) a PSD summary stat (e.g. median diameter or the
  30–150 µm habitable-volume fraction). Report both correlations and say which
  predicts better. Plot outcome vs isolation_ratio and vs the PSD stat.

## 10. Deliverables (add to / update the project)

- `pore_field.py` — empirical-PSD quantile mapping (§3); physiological `K(d)` (§4);
  rescaled `D(d)` (§5).
- `configs/mapping.yaml` — new K window (30/150/`w_decay`) and D `d_ref`, fixed.
- `configs/soil_*.yaml` — point each soil at its PSD file; shared spatial settings;
  different seeds; grid handled by the runner for both 128 and 40.
- `metrics.py` — outcome classification + viability axis (§8); `isolation_ratio`
  and the correlation analysis (§9).
- A Stage-3 runner (extend `stage2_run.py` or add `stage3_run.py`) — the θ sweep
  at both grids, all figures/CSVs.
- `results/` — respiration vs θ per soil per grid; water maps; isolation_ratio vs
  θ; the two correlation plots (outcome vs isolation_ratio, outcome vs PSD stat);
  the loaded empirical PSDs plotted per soil.
- `LOGBOOK.md`, `STAGE3_RESULTS.md`, `MODEL_SPEC.md` (add a Stage-3 section: real
  PSD, physiological K, metric change), `report_prompt.md` (extend to Stages 1–3).

## 11. Order of work

1. Load + validate the real PSDs (§2); plot them; confirm diameter + bins.
2. Empirical-quantile pore field (§3) + physiological K (§4) + D rescale (§5);
   sanity-run one soil at n=40 (respiration sane, stability assertion passes).
3. Metric revision (§8) + isolation diagnostic (§9).
4. Full θ sweep at n=128 and n=40 (§6, §7); figures; correlation test.
5. Checkpoint: STOP; write `STAGE3_RESULTS.md`, update `MODEL_SPEC.md`,
   `LOGBOOK.md`, `report_prompt.md`.

## 12. Honest framing

Answer two questions plainly: (1) with spatial arrangement shared, does the real
PSD alone make the three soils' respiration distinct — and at which θ/grid? (2)
Does the isolation ratio predict respiration better than the PSD? Report both
outcomes as they come, including "no". Biology + mapping frozen; only the real PSD
differs; die-offs are valid outcomes; stop at the budget and the checkpoint.
