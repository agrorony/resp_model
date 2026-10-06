# Claude Code task: Stage 6 — port the sandbox fix into Variant C (implicit transport + coating/trapped OM, with slow trapped release)

Bring the two conclusions the PDE sandbox proved (`pde_sandbox/`, RESULTS.md) into
the real dual-porosity soil model (`dual_porosity.py`, as left after Stage 5).
Biology frozen; numpy-only; grid n=128; constant moisture. This is a defined port +
run, not a search (`MAX_ITERATIONS = 12` ceiling if any tuning is needed).

## 0. What the sandbox proved (the recipe to port)

1. **Local implicit backward-Euler diffusion** (Jacobi-iterated, harmonic-mean
   faces) removes the explicit stability cap that kept dissolved substrate frozen
   ~1–3 cells from where it dissolved — so substrate can actually travel to the
   habitat. This is the transport operator in `pde_sandbox/sandbox.py`.
2. A **coating-dominated source with a finite connected path to the habitat** makes
   the compact/Sand case bloom fast then fall; a **distributed source** makes the
   broad/Vertisol case recruit more habitat over time and rise. The coating/trapped
   split and its geodesic-based placement are in `pde_sandbox/sandbox.py`.

Reuse those implementations directly; do not re-derive them.

## 1. Keep from Stage 5 (unchanged)

Matric-potential wetting (50 µm air-entry) + inscribed-circle macropore geometry;
the three phases (macropore / porous-matrix / solid); per-soil total OM; frozen
biology and K(d); numpy-only.

## 2. Change 1 — implicit transport for the substrate field

Replace the explicit diffusion stepping of the dissolved-substrate field with the
sandbox's **local implicit backward-Euler solve** (`u − dt·div(D∇u) = rhs` via
Jacobi iteration, harmonic-mean face coefficients so dry/solid cells are true
barriers). Port it from `pde_sandbox/sandbox.py`. D on each cell is the wet-cell
pore-size diffusivity from Stage 5 (0 on dry/solid). Verify: nonnegativity holds,
mass is conserved up to source/uptake/boundary, and the stability cap is gone.

## 3. Change 2 — coating / trapped OM split, per soil

Split each soil's total OM into two pools placed by their position relative to the
habitat (reuse the sandbox's geodesic-to-reaction coating construction):

- **Coating** — OM on porous-matrix cells within a short geodesic distance of a
  habitat (macropore) cell, along a wet connected path. Delivers fast; guarantees a
  finite connected habitat←source path.
- **Trapped** — OM in interior matrix cells far from habitat / behind barriers.

The **coating fraction is per-soil**, reflecting the organic-matter-distribution
hypothesis:
- **Sand**: coating-dominated (~0.8) — most OM coats the matrix surface next to
  habitat → fast flush.
- **Vertisol**: distributed (coating ~0.3, most OM in the bulk matrix at varied
  distances) → progressive recruitment.
- **Loess**: intermediate (~0.5).

Put these in `configs/mapping.yaml` (or per-soil configs) as documented parameters.

## 4. Change 3 — trapped OM releases SLOWLY (not locked)

Trapped OM must NOT be permanently undeliverable. Give it a **slow release**: the
trapped pool dissolves at a rate `k_dis_slow` ≪ `k_dis` (the coating/normal rate),
so it trickles substrate out through the interior over the whole run — a slow tail
after the coating burst. Expose `k_dis_slow`. Two OM pools per cell now (coating,
trapped); both feed the same dissolved-substrate field, at different rates.

## 5. Change 4 — remove the throttles that kept growth dead

- **Remove `S_cap`** entirely (no ceiling on usable substrate) — the sandbox had
  none and growth fired; it was the compounding blocker.
- **Turn wicking OFF** — use the static Stage-5 matric wetting map (the sandbox
  used a static structure). Wicking previously overrode the wetting geometry.
- Keep `k_leak` / matrix internal diffusion as the matrix→habitat exchange; now
  that implicit diffusion + the coating actually deliver substrate, these are
  supporting, not the sole pathway.

## 6. Run and compare

Run the three soils at constant moisture (n=128), plus an n=40 grid check.
Diagnostics per soil: **mean substrate reaching the habitat over time** (this MUST
now be far above Stage 5's ~1e-4 — it is the direct check that carbon arrives),
**habitats-recruited-over-time** (count of active/respiring habitats vs t), outcome
class (thriving/declining/dead), R(t) shape, and pairwise distinctness. Headline
figure: R(t) for the three soils vs the measured P1–P4 pattern.

Target (qualitative): **Sand — sharp rise then fall; Vertisol — rising; Loess —
low/erratic.** Sand's rise-then-fall should come from the coating flush depleting,
with the slow trapped release giving a low tail; Vertisol's rise from distributed
OM recruiting more habitat over time.

## 7. Deliverables

- `dual_porosity.py`: implicit BE substrate solve (ported); coating/trapped OM
  split with per-soil coating fraction; slow trapped release (`k_dis_slow`); `S_cap`
  removed; wicking off by default.
- `configs/`: per-soil coating fractions, `k_dis_slow`, coating thickness.
- `stage6_run.py`: the three-soil constant-moisture run + experiment overlay +
  the substrate-reaching-habitat and habitat-recruitment diagnostics + n=40 check.
- `results/stage6/`: R(t) per soil, experiment overlay, mean-S-habitat vs t,
  habitats-recruited vs t, coating/trapped OM maps per soil, distinctness CSV.
- `docs/stage_results/STAGE6_RESULTS.md`; update `MODEL_SPEC.md`, `LOGBOOK.md`,
  `report_prompt.md`.

## 8. Order of work

1. Port the implicit BE substrate solve; confirm mean-S-habitat is now well above
   Stage 5 (carbon actually arrives) on one soil.
2. Add the coating/trapped split (per-soil fraction) + slow trapped release.
3. Remove `S_cap`; set wicking off.
4. Run three soils, compare to the experiment; n=40 check.
5. Checkpoint: STOP; write `STAGE6_RESULTS.md` honestly.

## 9. Honest framing

Report plainly whether the ported fix makes carbon reach the habitat (the killer
metric) and whether Sand now shows rise-then-fall and Vertisol rise. Treat Loess as
it comes. Do not tune toward the target beyond the per-soil coating fractions (the
hypothesis-motivated inputs). These PSDs remain representative soil types, not the
exact experimental samples.
