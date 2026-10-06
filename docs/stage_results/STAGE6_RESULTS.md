# Stage 6 (EXPLORATORY) RESULTS -- porting the sandbox fix into Variant C

`prompts/soil_respiration_prompt_stage6_portback.md`. `MAX_ITERATIONS = 12`
not exercised: a defined port + run, not a search. Ports two conclusions
proved in an abstract PDE sandbox (`pde_sandbox/`, its own `RESULTS.md`) into
`dual_porosity.py` as left after Stage 5, changing nothing else (matric-
potential wetting, inscribed-circle macropore geometry, the three phases,
frozen biology/K(d) all kept, per SS1).

## 0. What was ported

1. **Implicit backward-Euler transport** (`dual_porosity.implicit_diffusion_
   step`, `_prepare_implicit_cache`, `_implicit_neighbor_sum` -- direct
   ports of `pde_sandbox.sandbox.local_implicit_step`/`prepare_local_cache`):
   solves `u - dt*div(D*grad(u)) = rhs` by Jacobi iteration
   (`stage6_implicit_iterations=90`, `stage6_implicit_relaxation=0.85`),
   harmonic-mean face conductivities so a dry/solid cell (D=0) is a true
   barrier. Unconditionally stable -- no `max(D)*dt/dx^2<0.2` cap.
2. **Coating/trapped OM split** (`dual_porosity.split_om_coating_trapped`,
   `bfs_geodesic_distance` -- ported from `pde_sandbox.sandbox.
   build_down_source`'s geodesic construction): a soil's OM is split into a
   COATING pool (porous-matrix cells within `stage6_coating_thickness=2`
   geodesic steps of a macropore cell, the distance measured along a wet
   matrix path) and a TRAPPED pool (everything else in the porous matrix).
   Per-soil coating fraction (`stage6_coating_fraction_{sand,loess,
   vertisol}` = 0.80 / 0.50 / 0.30, `configs/mapping.yaml`) is the one
   deliberately hypothesis-motivated input.
3. **Slow trapped release**: trapped OM dissolves at `k_dis_slow=0.005`
   (10x slower than the coating's `k_dis=0.05`) -- a tail, not a lock.
4. **Removed `S_cap` entirely**; **wicking off by default** (both retired/
   flipped in `dual_porosity.py`, `configs/mapping.yaml`'s
   `stage4_variantC_S_capacity_base`/`_cluster_S_scale` marked retired).
   `k_leak` and matrix internal diffusion (`k_matrix_diff`) stay as
   supporting pathways, now folded into the SAME implicit solve as the
   macropore D field rather than a separate zero-overlap explicit term.

A one-cell point-source unit test (`implicit_diffusion_step` on a 20x20
uniform-D grid) confirmed the port: mass conserved to `1e-13`, and the field
visibly spreads (not frozen) over repeated calls.

## 1. An important correction to the premise, found before running the port

**The old explicit stability cap was never actually close to binding in this
model.** Checked directly: `D_macro_base.max() * 1.75 * dt / dx^2` (the worst
case, with the full habitat-cluster D bonus) came out to 0.088 (Sand), 0.035
(Loess), 0.042 (Vertisol) -- all comfortably under the 0.2 cap, not near it.
So removing the cap, by itself, could not have been why substrate stayed
"frozen ~1-3 cells" in Stage 3.5/4/5 -- that reach was governed by the
genuine physical diffusion length `sqrt(D*T_total)/dx` given this model's
D(d) magnitudes and the run's total simulated time, which implicit vs.
explicit stepping does not change for a fixed, already-stable dt.

**What actually changed the picture is a real topological fix bundled into
Change 1, not the solver swap itself**: Stage 4/5 ran substrate transport as
TWO separate explicit divergences (`D_macro_t`, `D_matrix_t`), each zero
wherever the other was nonzero -- so a macropore cell and its adjacent
matrix cell could NEVER diffuse into each other directly; the only pathway
crossing that boundary was the special-cased `k_leak` term. Stage 6 combines
them into `D_total = D_macro_t + D_matrix_t` and solves ONE connected
implicit system, so real spatial diffusion between matrix and macropore
cells exists for the first time. This is the actual mechanism behind
whatever improvement follows below -- worth stating plainly since it is not
what SS0 of the prompt attributes the fix to.

## 2. Does carbon reach the habitat now? (the killer metric)

| soil | Stage 5 (wicking off, static map) | Stage 6 | ratio |
|---|---|---|---|
| Loess    | 0.00017 | 0.00083 | **4.9x** |
| Sand     | 0.00033 | 0.00034 | 1.0x (flat) |
| Vertisol | 0.00026 | 0.00098 | **3.8x** |

(`mean_S_habitat_final`, `results/stage6/stage6_headline.csv` vs.
`results/stage5/stage5_headline.csv`'s wicking-off rows -- the fair
comparison, since Stage 6 also runs with wicking off.)

**Loess and Vertisol do show a real, meaningful increase** (~4-5x) in mean
substrate reaching the habitat, and their `habitats_recruited` diagnostic
confirms this is genuinely distributed: 58/76 and 55/70 macropore regions
show active growth by the end of the run (`stage6_habitats_recruited_vs_t.
png`), vs Stage 5 where the equivalent notion barely existed. **Sand shows
essentially no improvement** (1.0x): only 1 of its 3 macropore regions
(the tiny 8-cell one) ever recruits, because Sand's habitat is still almost
entirely the same one huge, dry, 340um region Stage 5 identified -- there is
barely any wet-adjacent porous matrix for the coating construction to place
OM near in the first place (1,844 of Sand's ~1,953 porous-matrix cells at
n=128 end up "coating" only via the geodesic distance-2 rule, but that
coating sits mostly around the tiny wet region, not the huge dry one that
holds Sand's actual habitat capacity).

**Despite this real, mechanistically-explained improvement, `mean_S_habitat`
stays two-to-three orders of magnitude below `Ks=0.2`** (the biology's own
half-saturation constant) for every soil (`stage6_mean_s_habitat_vs_t.png`).
`f_S = S/(S+Ks)` therefore never rises above roughly 0.004 anywhere, growth
never exceeds maintenance decay, and **all three soils remain `declining`,
pinned near the trivial `R_peak ~ m0*B0_total = 0.4`** exactly as in Stage 5.
Emergent distinctness is accordingly still degenerate
(`min_distance=0.0035`, up from Stage 5's 0.0001-0.0013 but nowhere near the
0.3 threshold).

## 3. Target shapes: not reached

| soil | target | Stage 6 model shape | match? |
|---|---|---|---|
| Sand | sharp rise then fall | falling (no rise at all) | no |
| Vertisol | rising | falling | no |
| Loess | low/erratic | falling | no (only accidentally matches Sand's target, not its own) |

None of the three qualitative targets are reproduced. Sand's hypothesized
"coating flush" never fires because there is no meaningful live habitat for
a flush to feed into; Vertisol's hypothesized "progressive recruitment"
DOES happen structurally (55 of 70 regions recruit, `habitats_recruited_t`
rises steadily) but the absolute magnitude of recruited growth is too small
to ever overcome maintenance/starvation decay, so `R(t)` still only falls.

## 4. Grid check (n=40, `stage6_grid_check.csv`)

Same qualitative story at n=40: all three `declining`/`falling`; wet-habitat
fractions track the n=128 values reasonably (Loess 0.039 vs 0.024, Vertisol
0.032 vs 0.025, Sand 0.000 at both grids); `mean_S_habitat_final` stays in
the same 1e-4 to 1e-3 range. The finding is grid-robust, in the sense that
it is not an n=128-specific artifact.

## 5. Honest verdict (SS9)

**The two sandbox-proven mechanisms were ported faithfully and DO measurably
help where a live wet macropore network exists (Loess, Vertisol: ~4-5x more
substrate reaching habitat, dozens of regions recruited) -- but this is not
enough to revive growth, and Sand barely benefits at all.** Two reasons,
both outside what Stage 6 was asked to change:

1. **A modest total-OM budget.** The per-soil `om_fraction`-based
   initialization (inherited from Stage 4) gives each soil only
   `om_total_init` = 39 (Sand), 253 (Loess), 361 (Vertisol) units of total
   organic carbon at n=128 -- an order of magnitude below the single shared
   `OM_total=800` constant used through Stage 3.5. Even perfect delivery of
   all of it directly to habitat, with no dilution at all, would not
   obviously clear `Ks=0.2` once divided across thousands of habitat cells.
2. **A genuinely short physical diffusion length** given this model's D(d)
   magnitudes, `dt=0.05`, `dx=1.0`, and total simulated time (`T*dt=30`):
   `sqrt(D_max*T)` for the realistic region diameters seen here is on the
   order of a few grid cells, regardless of solver. The sandbox's own
   parameters used a much finer relative grid resolution (`dx=1/32` on a
   unit domain) that gave it far greater dimensionless reach for a
   comparable `D*T` -- a regime difference this port did not (and was not
   asked to) reconcile.

**Sand specifically remains the recurring problem soil across every stage**:
its macropore geometry (one huge, well-connected, ~340um region under the
Stage-5 matric rule) leaves it with almost no live wet habitat for either
the wetting rule (Stage 5) or the coating placement (Stage 6) to exploit,
regardless of its very high (0.8) coating fraction -- the fraction cannot
matter if there is barely any wet-adjacent habitat to place a coating next
to. This is reported as the mechanism, not tuned around, per SS9.

## 6. Deliverables

`dual_porosity.py`: `implicit_diffusion_step` (+ `_prepare_implicit_cache`,
`_implicit_neighbor_sum`), `bfs_geodesic_distance`,
`split_om_coating_trapped`, `count_active_habitats`, `_variantc_coating_
fraction`; `simulate_dual_porosity` rewritten for two OM pools + implicit
transport + no `S_cap` + wicking off by default. `configs/mapping.yaml`:
`stage6_implicit_iterations`, `stage6_implicit_relaxation`, `k_dis_slow`,
`stage6_coating_thickness`, `stage6_coating_fraction_{sand,loess,vertisol}`
(S_cap constants retired). `stage6_run.py`: the three-soil constant-moisture
run, experiment overlay, substrate-reaching-habitat and habitat-recruitment
diagnostics, coating/trapped OM maps, n=40 check. All figures/CSVs in
`results/stage6/`.

**[Stage 6] Checkpoint reached: implicit transport ported and unit-verified,
coating/trapped OM split with per-soil fractions and slow release,
`S_cap` removed, wicking off, three-soil run + n=40 check + all diagnostics
complete; deliverables written (this file, MODEL_SPEC.md, LOGBOOK.md,
report_prompt.md, README.md). Stopping here per soil_respiration_prompt_
stage6_portback.md order-of-work -- MAX_ITERATIONS=12 was not needed.**
