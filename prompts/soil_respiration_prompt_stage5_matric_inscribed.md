# Claude Code task: Stage 5 (EXPLORATORY) — matric-potential wetting + inscribed-circle macropores (Variant C)

Extend the EXISTING numpy-only project. Build on the **latest Variant C**
(`dual_porosity.py`, run by `stage4_run.py`) — keep all of its current mechanics
(three phases macropore/porous-matrix/solid, per-soil OM fractions, `k_leak`
matrix→macropore leak, matrix internal diffusion, wicking, habitat-cluster K/D/S
bonuses, S capacity). Change ONLY the two things below, run the three soils, and
compare to the measured pattern. This is exploratory: budget-capped
(`MAX_ITERATIONS=12`), a defined run, not a search — checkpoint and stop.

Validation target (qualitative, constant moisture): **Vertisol rising, Sand
falling, Loess low/erratic.**

## Why this stage

Stage 4 showed Variant C's blocker is upstream of all its transport machinery:
the water mask is still the per-soil `global_theta` quantile of the FULL PSD, so
Loess/Vertisol habitat starts bone dry (wet-habitat ~0) and nothing the leak /
wicking / OM terms do can feed it. This stage replaces that wetting rule with a
**shared matric potential**, and gives macropores a real geometric size.

## Change 1 — matric-potential wetting (shared across soils)

Replace the `global_theta` quantile wetting with a **fixed air-entry diameter,
identical for all three soils** (this is the shared matric potential):

- **A pore is saturated (wet) if its size < 50 µm, and dry (air-filled) if its
  size ≥ 50 µm.** Put `matric_d_cut_um = 50` in `configs/mapping.yaml`.
- For **matrix / fine cells** the size is the cell's own PSD diameter `d` (all
  < 30 µm, so all wet — the fine matrix stays saturated and keeps holding water +
  OM, as intended).
- For **macropore cells** the size is the region inscribed diameter from Change 2
  (NOT the per-cell PSD label, NOT the forced 80 µm).
- **Solid cells are never wet.**
- θ then *emerges* per soil (fraction of pore volume that is wet) rather than being
  imposed — report each soil's emergent θ. Do not sweep θ; run this single matric
  condition. (Keep a small 2–3-value sensitivity on `matric_d_cut_um`, e.g.
  40/50/60 µm, as a secondary robustness check only.)

## Change 2 — reduce macropore cell size and size each macropore region by its maximum inscribed circle

Stop forcing every macropore cell to 80 µm. Instead give the grid a physical
resolution and measure macropore size geometrically:

- Add `voxel_um` to `mapping.yaml` (the physical size of one cell; start at
  **10 µm/cell** — i.e. "reduce the macropore cell" from 80 µm to a small voxel so
  a macropore is many cells, not one). Document that the grid now has a real
  physical scale.
- **Macropore cells** are the connected void where `d ≥ 30 µm` (`K_d_low`), as now.
  Label the connected macropore **regions** (4-connectivity — reuse the existing
  component-labeling in `dual_porosity`).
- For each region, compute its **maximum inscribed circle diameter**: the largest
  disk that fits inside the region. Use a numpy-only distance transform
  (e.g. iterative 4-neighbour erosion count = distance-to-boundary in cells; region
  max × 2 × `voxel_um` = inscribed diameter in µm). If the repo already has a
  local-thickness / inscribed-circle helper (the PSD code may), reuse it; otherwise
  implement the erosion version, numpy-only.
- **Use this region inscribed diameter as the region's pore size** for: (a) the
  wet/dry decision in Change 1 (region < 50 µm → the whole region is wet; ≥ 50 µm →
  dry), and (b — recommended, since it fixes the 80 µm flattening Stage 4 flagged)
  the region's `K(d)` and `D(d)`, so a 40 µm water-holding habitat pore and a
  500 µm drained channel are no longer treated identically.

Net physical picture this creates: the fine matrix is always wet (holds water +
OM); a **small macropore region (< 50 µm inscribed) stays water-filled** and is
live, wet habitat fed by the surrounding matrix; a **large macropore region
(≥ 50 µm) drains and goes dry**. Which soils keep wet habitat now depends on their
macropore *geometry*, not a per-soil quantile — exactly the lever we want to test.

## Keep everything else

All other updated-Variant-C mechanics stay: per-soil OM fractions
(`stage4_variantC_om_fraction_*`), `k_leak`, matrix internal diffusion, wicking
(rewetting + directional transfer), habitat-cluster bonuses, `S_cap`. Biology
frozen; numpy-only; grid n=128 (n=40 as a grid check). If wicking's dynamic
rewetting conflicts with the new static wet map, run once with wicking ON and once
OFF so its contribution is isolated and reported.

## Diagnostics & comparison

Per soil report: emergent θ, **wet-habitat fraction**, the distribution of
macropore-region inscribed diameters (how many regions land < 50 vs ≥ 50 µm),
**mean substrate reaching the habitat** (the Stage-3.5/4 killer metric), outcome
class, and R(t) shape vs the measured shape. One comparison figure (three soils'
R(t) vs the P1–P4 pattern) and a per-soil map panel showing macropore regions
coloured wet/dry by inscribed size. Report pairwise distinctness too, but shape/
ranking match is the headline.

## Deliverables

- `dual_porosity.py`: new `matric_d_cut_um` wetting path; a
  `macropore_region_inscribed_diameter` function; region-size-driven K/D/wet.
- `configs/mapping.yaml`: `matric_d_cut_um=50`, `voxel_um=10`; retire/replace the
  forced `stage4_variantC_macro_diameter_um`.
- `stage5_run.py` (or a Variant-C mode switch in `stage4_run.py`): the run + the
  `matric_d_cut_um` sensitivity + wicking on/off.
- `results/stage5/`: comparison figure, per-soil wet/dry inscribed-region maps,
  CSVs (per soil: emergent θ, wet-habitat frac, inscribed-diameter stats,
  mean_S_habitat, outcome, shape).
- `docs/stage_results/STAGE5_RESULTS.md`; update `MODEL_SPEC.md`, `LOGBOOK.md`,
  `report_prompt.md`.

## Order of work

1. Add `voxel_um`; region labeling → max inscribed circle per macropore region;
   sanity-plot the inscribed-diameter distribution per soil.
2. Swap wetting to the fixed 50 µm matric cutoff (region size for macropores,
   per-cell `d` for matrix); wire region size into K/D.
3. Run three soils at n=128, wicking on and off; the `matric_d_cut_um` sensitivity;
   n=40 check.
4. Diagnostics + comparison figure.
5. Checkpoint: STOP; write `STAGE5_RESULTS.md`, update the rest.

## Honest framing

Report plainly whether the matric-potential wetting + geometric macropore sizing
gets live wet habitat into Loess/Vertisol and whether the three R(t) shapes now
match the experiment — and if Sand's large drained macropores leave it too dry to
respire, say so (that is a real, informative outcome of this wetting rule, not a
bug). Don't tune toward the target; report the mechanism via the wet-habitat and
mean-substrate-reaching-habitat diagnostics. Note `voxel_um`, `matric_d_cut_um`,
and the literature-PSD representativeness as the key assumptions.
