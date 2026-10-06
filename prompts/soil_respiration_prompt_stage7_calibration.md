# Claude Code task: Stage 7 — calibration to close the substrate gap (uniform OM %)

Build on Stage 6 (`dual_porosity.py`). Keep its mechanism intact; this stage is a
**quantitative calibration**, not a new mechanism. Biology frozen; numpy-only;
grid n=128 (n=40 check). Budget: `MAX_ITERATIONS = 15` for the shared calibration.

## 0. Why this stage

Stage 6 ported the transport fix and the coating/trapped source. The mechanism
works *structurally* (Loess/Vertisol recruit 55–58 habitat regions), but the
substrate reaching the habitat stays ~0.001 — about **25× too low** for growth to
beat maintenance. The explicit threshold: with `r_max=1, m0=0.02, m_s=0.1, Ks=0.2`,
growth exceeds death only when **f_S > ~0.11**, i.e. **habitat substrate S ≳ 0.025**.
Stage 6 identified two quantitative gaps and one structural one:

1. the OM budget is ~10× too small (per-soil totals 39/253/361 vs the 800 used
   earlier),
2. the physical diffusion reach √(D·T) is only a few cells in this model's units,
3. Sand's habitat is one big drained (dry) region with no live wet habitat.

This stage closes all three, and returns every soil to the SAME OM percent.

## 1. Keep from Stage 6 (unchanged)

Unified implicit transport (matrix + macropore solved as ONE connected system);
coating/trapped OM split with per-soil coating *placement* from geometry; slow
trapped release (`k_dis_slow`); `S_cap` removed; wicking off; matric-potential
wetting + inscribed-circle macropore geometry; frozen biology / K(d).

## 2. Change 1 — uniform OM percent across all soils

Replace the per-soil OM fractions (Sand 0.02 / Loess 0.035 / Vertisol 0.05) with a
SINGLE shared `om_fraction` applied to every soil (same percent of porous matrix).
Now the only per-soil difference is **structure** — pore geometry and the coating
placement that follows from it — not the amount of carbon. Report each soil's
resulting `om_total`. (The coating vs trapped *split* still follows each soil's own
geometry; only the total percent is unified.)

## 3. Change 2 — close the substrate budget + reach gap (shared calibration)

The target: fed-habitat substrate must clear **S ≳ 0.025 (f_S > 0.11)** in the soils
that should grow. Tune, SHARED across all soils (never per-soil), within the budget:

- the shared `om_fraction` (raise it toward whatever clears the threshold);
- the transport reach — raise the D(d) magnitude and/or reduce the effective `dx`
  (or shrink the domain scale) so √(D·T) spans the matrix→habitat distances,
  matching the dimensionless reach the sandbox had; optionally lengthen T.

Keep every calibrated value shared across the three soils, so structure remains the
only cause of between-soil differences. Re-verify numerical stability of the
implicit solve after changing D/dx/T.

## 4. Change 3 — Sand's habitat water films

At constant high humidity a drained macropore (≥50 µm, currently D=0, dead habitat)
retains a thin wall **water film**. Give dry macropore cells a small nonzero film
diffusivity `D_film` and let them host activity at reduced capacity, instead of
being fully cut off. This gives Sand's big dry habitat a **fast but low-capacity**
supply from its adjacent coating → an early flush that depletes → fall. `D_film` is
a shared physical constant (small relative to the wet D), applied to every soil's
drained macropores; report how much live habitat it opens up per soil.

## 5. Calibration loop and success

Small budget-capped search over the SHARED scalars (`om_fraction`, D-scale,
`D_film`, optionally `dx`/`T`) — biology and per-soil structure fixed. Score each
candidate on the three-soil run by:

- **Growth cleared:** fed-habitat f_S exceeds the growth threshold where structure
  allows it (report per soil).
- **Correct distinct shapes:** Sand = sharp rise then fall, Vertisol = rising,
  Loess = low/erratic — and the three curves are pairwise distinct (> 0.3). Penalize
  the degenerate "everything blooms and rises together" outcome (over-fed
  under-separation) as much as the "everything starves" one.

Do not tune per-soil to force shapes; only the shared scalars move.

## 6. Run, compare, diagnose

Three soils, constant moisture, n=128 (+ n=40 check). Report per soil: `om_total`,
fed-habitat f_S and mean_S_habitat over time (must clear ~0.025 where growth is
expected), habitats-recruited-over-time, `D_film` live-habitat gained, outcome
class, R(t) shape, pairwise distinctness. Headline figure: R(t) for the three soils
vs the measured P1–P4 pattern.

## 7. Deliverables

- `dual_porosity.py`: single shared `om_fraction`; the reach/D-scale calibration
  hooks; `D_film` for drained macropores.
- `configs/mapping.yaml`: shared `om_fraction`, D-scale, `D_film` (per-soil OM
  fractions retired), plus whatever `dx`/`T` the calibration settles on.
- `stage7_run.py`: the shared calibration loop + the three-soil run + experiment
  overlay + the f_S / substrate-reaching-habitat / recruitment diagnostics + n=40.
- `results/stage7/`: R(t) per soil, experiment overlay, f_S-vs-t and
  mean_S_habitat-vs-t, habitats-recruited-vs-t, distinctness CSV, the winning
  shared parameters.
- `docs/stage_results/STAGE7_RESULTS.md`; update `MODEL_SPEC.md`, `LOGBOOK.md`,
  `report_prompt.md`.

## 8. Order of work

1. Uniform `om_fraction`; confirm each soil's `om_total`.
2. Add `D_film` for drained macropores; confirm Sand gains live habitat.
3. Run the shared calibration (§5) to push fed-habitat f_S past ~0.11 while keeping
   the three shapes distinct.
4. Three-soil run + experiment overlay + n=40 check.
5. Checkpoint: STOP; write `STAGE7_RESULTS.md` honestly.

## 9. Honest framing

Report the fed-habitat f_S actually reached per soil (the direct check that the gap
closed), and whether Sand now shows rise-then-fall, Vertisol rising, Loess
low/erratic, with all three distinct. If closing the substrate gap makes everything
bloom-and-rise together (under-separation), say so — that is as much a failure as
starvation, and it tells us the shapes need the geometry, not just more carbon. Only
the shared scalars were tuned; per-soil differences remain structural. PSDs stay
representative soil types, not the exact samples.
