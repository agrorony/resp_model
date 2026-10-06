# Claude Code task: Stage 4 (EXPLORATORY) — two fixes for "OM never reaches the habitat", head-to-head

Extend the EXISTING numpy-only project. This is an **exploratory** stage: build
**two** variants that attempt to fix the Stage-3.5 failure, run both, and
**compare** them against my measured pattern. Do not commit to one in the code —
keep both runnable and let the comparison decide.

## 0. Budget & framing

`MAX_ITERATIONS = 12` ceiling. This is a defined comparison, not a search — build
both variants, run, compare, write up. Stop at the checkpoint.

Stage 3.5 diagnosis (the thing both variants fix): the frozen `OM ∝ d^-2.5` rule
dumped ~99.96% of the carbon onto sub-micron clay cells, from which diffusion
(reach ~1–3 cells/run) could never deliver it to the 30–150 µm habitat, so two of
three soils just starved to the same dead curve. Validation target is my
**constant-moisture** experiment (qualitative shape/ranking only):
Vertisol **rising**, Sand **falling**, Loess **low/erratic**.

## 1. Shared setup (identical for both variants — hold everything else fixed)

- Constant moisture: single fixed θ held for the whole run (start θ ≈ 0.6, same
  for all soils), plus the 2–3-value θ robustness check as in Stage 3.5.
- Physiological `K(d)` (0 below 30 µm, plateau 30–150, decay above 150); D(d) with
  `D_d_ref = 250` (the Stage-3.5 fusion fix); literature PSDs
  (`configs/psd_literature.yaml`); real soil names (Loess/Sand/Vertisol);
  grid n=128 (n=40 as the robustness check). Biology frozen; numpy-only.

## 2. Variant A — "truncate + relocate OM" (the lightweight fix to 3.5)

Keep the current single-continuum pore model, change only two things:

1. **Truncate the PSD at a 10 µm floor.** Drop all pore volume below 10 µm and
   renormalize, so the finest cell represented is 10 µm. Sub-10 µm clay is simply
   *not on the grid* (unresolved) — this removes the sub-micron OM traps by
   construction. Report each soil's post-truncation porosity, habitable fraction,
   and how much volume was dropped.
2. **Relocate OM to fine-but-reachable pores.** Replace the unbounded `d^-2.5`
   rule with a **banded** placement: OM density peaks in the ~10–30 µm band (fine,
   water-retaining, *just below* the habitat window) and decays away from it, so
   the carbon sits where release is **slow but the substrate can actually diffuse
   into the adjacent 30–150 µm habitat within the run**. Make the band center/width
   an explicit tunable in `mapping.yaml`. The intent: slow arrival, not *no*
   arrival.

Everything else exactly as Stage 3.5. This tests whether just removing the
unreachable traps and seating OM next to the habitat is enough.

## 3. Variant C — dual-porosity (matrix + macropore)

Two cell types instead of one continuum (the "sponge + straws" model):

- **Macropore cells** — pores with `d ≥ 30 µm` (habitat + channels): fast D(d),
  biomass lives and respires here, substrate diffuses between them as now.
- **Matrix cells** — everything finer (the `< 30 µm` fraction, incl. the sub-10 µm
  that Variant A drops): **no biomass**, block fast inter-macropore flow, but each
  matrix cell carries an internal fine-porosity (water + stored OM) and **leaks
  dissolved substrate into any adjacent macropore cell** at rate `k_leak`
  (a mobile–immobile / dual-permeability exchange term). Put most OM in the matrix.
- Represent the sub-resolution fine fraction as the matrix cells' internal porosity
  (a scalar), not as separate cells — so bulk porosity is partly captured (Vertisol
  = lots of high-internal-porosity matrix; Sand = little).

Keep the biology inside macropores identical to every other stage. Add `k_leak`
(and the matrix internal-porosity mapping) as documented, shared constants — not
per-soil knobs; only the PSD (hence the macropore/matrix split) differs per soil.

## 4. Comparison protocol (the point of this stage)

Run **both variants** for all three soils at constant θ (n=128), plus the θ- and
grid-robustness checks. Produce ONE comparison figure: model R(t) for Variant A
and Variant C, three soils each, next to the experimental P1–P4 shapes. For each
variant report:
- Per-soil model shape (rising/falling/erratic) vs measured, and the ranking.
- Pairwise distinctness + outcome (thriving/declining/dead) with
  `require_nontrivial=False`.
- The mechanism diagnostics: wet-habitat fraction, and crucially the **mean
  substrate concentration reaching the habitat** (the Stage-3.5 killer metric) —
  show whether each variant actually gets carbon to the microbes.
- Isolation ratio as before.

Then state plainly: **does Variant A suffice, or is the explicit dual-porosity
exchange (C) needed** to reproduce Vertisol-rising / Sand-falling / Loess-erratic?

## 5. Deliverables

- Variant A: PSD-truncation + banded-OM code (a `pore.mode` flag or a runner
  switch); `mapping.yaml` band params.
- Variant C: the dual-porosity model (matrix/macropore split + `k_leak` exchange +
  matrix internal porosity), as a clearly separated code path so Stage 1–3.5 still
  run.
- A Stage-4 comparison runner — both variants, the comparison figure, all CSVs.
- `results/stage4/` — comparison figure, per-variant R(t), substrate-reaching-
  habitat and wet-habitat diagnostics, θ/grid robustness.
- `docs/stage_results/STAGE4_RESULTS.md` — honest head-to-head verdict (which
  variant reproduces the pattern, which mechanism explains it, what each still
  gets wrong). Update `MODEL_SPEC.md`, `LOGBOOK.md`, `report_prompt.md`.

## 6. Order of work

1. Variant A: truncate PSD ≥10 µm + banded OM; run three soils at constant θ;
   check the substrate-reaching-habitat metric.
2. Variant C: dual-porosity matrix/macropore + `k_leak`; run three soils.
3. Comparison runner + figure + robustness (θ, n=40).
4. Checkpoint: STOP; write `STAGE4_RESULTS.md`, update the rest.

## 7. Honest framing

Report which variant (if either) reproduces the measured ranking and shapes, and
explain via the mechanism diagnostics — not by tuning to the target. These PSDs
are literature-informed representative types, not my exact samples; the point is
whether a reachable-OM continuum (A) or an explicit dual-porosity exchange (C) is
what the model needs. If neither works, say so and point to why.
