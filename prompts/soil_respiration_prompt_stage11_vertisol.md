# Claude Code task: Stage 11 — freeze Sand, fix Vertisol, optimize the loop (n=128 only)

Continue from Stage 10 (`dual_porosity.py`, `stage10_run.py`). Stage 10 gave Sand a
strong early burst-then-crash (good — that is DONE) but left Vertisol at almost
nothing. This stage: **keep Sand exactly as it is, work only on Vertisol, fix the
per-timestep performance bug, and run at n=128 only** (no n=40 anywhere). numpy-only.

## 1. Freeze Sand

Do NOT change any Sand parameter or re-tune Sand. Its burst-then-crash from Stage 10
is final. After every change below, confirm Sand's curve is unchanged.

## 2. Performance / bug fix (do this FIRST)

- **Precompute the static habitat regions ONCE, before the time loop.** `macro_mask`
  is fixed and wicking is off, so the connected macropore regions, their inscribed
  sizes, and the habitat-cluster scale never change during a run — but the simulate
  loop currently recomputes `habitat_properties` (pure-Python connected-component
  labeling) EVERY timestep. Compute it once before the loop and reuse it. This is the
  main waste.
- Hoist any other per-timestep recomputation of static quantities out of the loop.
- **Reduce the implicit-solver work per step:** cap the Jacobi iterations (e.g. ~20–30)
  and/or stop early on a convergence tolerance, instead of a fixed ~90.
- **Verify the optimization changes nothing:** rerun Sand before vs after and confirm
  its R(t) is identical within numerical noise (report R_peak, R_end for both). The
  speedup must not alter results.

## 3. Fix Vertisol only

Vertisol currently barely respires (near-zero). Raise its accessible carbon so growth
actually fires across its many habitat regions, tuned to a **slow, gradual rise** —
late-third mean `>` early-third mean, still rising / peaking late, clearly NOT a burst
and NOT flat/near-zero. Vertisol-only levers (do not touch Sand):

- total OM (raise it so growth clears across its regions),
- coating fraction / availability (enough of its carbon reachable),
- effective delivery — if it front-loads into an early spike, SLOW Vertisol's
  delivery (lower its per-soil transport magnitude and/or place carbon farther) so
  recruitment staggers over the window → a gradual rise rather than a spike.

Target: growth genuinely fires (Vertisol `R_peak` clearly above the maintenance floor
`m0·B0`, biomass grows, `habitats_recruited(t)` climbs across the window), and the
curve rises slowly through the run.

## 4. Success

- **Sand:** unchanged (burst-then-crash identical to Stage 10).
- **Vertisol:** a clear slow rise (growth fires; not near-zero, not a burst).
- Two distinct curves. n=128 ONLY.

## 5. Deliverables

- `dual_porosity.py`: static habitat regions precomputed once; reduced implicit
  iterations; Sand params frozen; Vertisol params updated.
- `stage11_run.py`: two-soil (Vertisol, Sand) run at **n=128 only**; the Sand
  before/after optimization identity check; Vertisol slow-rise diagnostics
  (`R_peak` vs maintenance floor, biomass(t), `habitats_recruited(t)`,
  `mean_S_habitat(t)`).
- `results/stage11/`: the two-soil `R(t)` figure vs the measured shapes (Sand
  falling, Vertisol rising), Vertisol diagnostics, distinctness, chosen Vertisol
  parameters, and the before/after runtime + Sand-identity confirmation.
- `docs/stage_results/STAGE11_RESULTS.md`; update `MODEL_SPEC.md`, `LOGBOOK.md`,
  `report_prompt.md`.

## 6. Order of work

1. Performance fix (§2); verify Sand's curve and speed both improve/unchanged.
2. Tune Vertisol only (§3) to a slow rise; leave Sand frozen.
3. Final two-soil run at n=128; comparison figure + diagnostics.
4. STOP; write `STAGE11_RESULTS.md`.
