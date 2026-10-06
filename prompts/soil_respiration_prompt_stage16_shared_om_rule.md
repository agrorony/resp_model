# Claude Code task: Stage 16 — restore the shared OM rule; let structure set fuel–habitat geometry

A MODEL-WIDE structural change, not a calibration stage and not a search.
Delete the per-soil coating/trapped OM split introduced in Stage 6 and restore
ONE shared, structure-blind organic-matter placement rule, with spatially
UNIFORM dissolution. Both soils (Sand, Vertisol) are re-run under it. The
point: every difference — between soils, and between macropore regions within
a soil — must emerge from the pore field, not from a hand-set carbon
parameter. n=128 primary, n=40 check. numpy-only, as always.

## 0. Read this first — Stage 14's causal attribution was WRONG

`docs/stage_results/STAGE14_RESULTS.md` concludes that Vertisol's
near-instant habitat activation is caused by "the unconditionally-stable
implicit solver's per-step equilibration." **That attribution has been
refuted by direct test.** Do not inherit it.

The test (Vertisol, Stage 11 frozen config, n=128, T=1500):

| run | regions active @ t=2 | recruitment reaches final at | R_peak |
|---|---|---|---|
| baseline (`D_scale=0.7`) | 34 | 1.33% of window | 3.756 |
| **transport OFF (`D_scale=1e-9`)** | **17** | **0.27% of window** | 0.400 |
| transport OFF + `coating_fraction=0` | 0 | never activates | 0.400 |
| transport ON + `coating_fraction=0` | **0** | **9.27% of window** | 0.400 |

With transport switched off entirely, regions still activate within 2–4
timesteps. The solver is not the cause. The actual cause is structural and
lives in `split_om_coating_trapped`: the coating pool is defined as every
porous-matrix cell within `coating_thickness` geodesic steps of ANY
macropore cell, which builds an identical fuel shell around all ~70 regions
simultaneously; combined with `B0 ∝ K` seeding biomass in every region at
t=0 and a single global `k_dis`, every region is self-fuelled from t=0 and
needs no transport to activate. Simultaneity is built into the
construction.

Note the last row: removing the coating restores genuinely staggered
recruitment. That is what this stage builds on.

Stage 14's *diagnostic* findings stand and remain useful — the
`debug_growth_distribution` instrumentation, the per-region
time-to-half-max analysis, and the finding that the boolean
`habitats_recruited(t)` count is actively misleading (r=0.072 against R(t))
while a growth-weighted "active mass" metric tracks it (r=0.996). Only the
causal attribution is wrong.

**Deliverable obligation:** append a dated correction section to
`STAGE14_RESULTS.md` recording this. Do not delete or rewrite its original
text — correct it in place by addition, consistent with how this project
has handled every other revised finding.

## 1. The model change

Remove the coating/trapped split as the active default:

- ALL organic matter is placed by the SAME shared rule for every soil: the
  original fine-pore weighting `OM(i,j) ∝ d(i,j)^-om_fine_bias` within the
  porous matrix (`MODEL_SPEC.md` §4), normalized to that soil's total.
- Dissolution is spatially UNIFORM: a single `k_dis` everywhere. There is no
  `k_dis_slow`, no coating pool, no trapped pool, and — this is the
  invariant this stage enforces — **no term anywhere in the model whose rate
  depends on distance to habitat.** Decomposition is a local process
  (water contact, OM chemistry); distance may enter ONLY through diffusion
  of the dissolved substrate, which the model already does.
- `coating_fraction` (per-soil, currently Sand 0.50 / Vertisol 0.30) is
  retired as an active parameter. Removing it is part of the point: it was a
  hand-set per-soil carbon knob, and the project's central claim is stronger
  without it.

Keep the Stage 6–15 coating code path reachable behind a config flag (e.g.
`om.mode: coating_trapped` vs the new `om.mode: shared_fine_pore`), exactly
as Stage 13 preserved `psd_parametric` alongside `psd_measured`, so every
earlier stage stays reproducible. Do not delete it.

## 2. The structural variable to measure — this stage's primary result

For each soil, compute and report the **distribution of geodesic distance
from each macropore region to its nearest substantial OM**, using the
existing `bfs_geodesic_distance`. "Substantial" needs a stated, defensible
definition (e.g. the nearest matrix cell whose OM mass is above the median
of nonzero matrix OM) — state whichever you choose and use it consistently.

Report the full distribution per soil, not just a mean: min / p25 / median /
p75 / p95 / max, plus a histogram figure, in grid cells AND in µm (via
`voxel_um`).

This is the emergent quantity that carries the stage's argument. Sand is one
giant macropore region in an ~18%-matrix soil; Vertisol is ~70 regions in an
~68%-matrix soil. Under one shared OM rule those two pore fields must
produce different fuel–habitat distance distributions, and that difference —
not any per-soil parameter — is the mechanism by which structure controls
recruitment. Report it as such, and report honestly if the two distributions
come out more similar than expected.

For external calibration: Portell et al. (2018, *Front. Microbiol.* 9:1583),
a pore-scale individual-based model on CT-derived soil, found geodesic
distance through the liquid phase gates colony growth, with operative
distances of 252 ± 111 µm and 490 ± 262 µm in their two particulate
scenarios. Their domain was ~6.8 mm; this model's is 1.28 mm at n=128,
`voxel_um=10`. Note whether this model's emergent distances land in a
comparable range, and note the domain-size difference plainly.

## 3. Dose and transport — a ladder, adopt the highest rung that passes

The per-soil `om_total` and `D_scale` established in Stages 10/11 are
hand-set. Test how much of that can be derived from structure instead.
Work down this ladder and adopt the FIRST rung that passes §4:

1. **One shared `om_total`** and one shared `D_scale` for both soils.
   (Evidence below says this will probably fail; run it anyway and report
   the numbers, because "maximum purity was tried and failed, here is by how
   much" is a result.)
2. **Shared OM concentration per unit matrix volume**:
   `om_total_soil = rho_shared × matrix_volume_soil`. This derives the
   per-soil total from structure rather than setting it. Matrix fractions
   are ~0.68 (Vertisol) and ~0.18 (Sand), a 3.8× ratio; the evidence below
   suggests ~17× is needed, so expect this to close part but not all of the
   gap. **Report explicitly what fraction of the required per-soil
   difference this rung accounts for** — that number is itself a finding.
   Consider and report the analogous structure-derived option for `D_scale`.
3. **Per-soil `om_total` / `D_scale`**, as Stage 10 established. If you land
   here, state plainly in the write-up how much of the between-soil
   difference is structure-derived and how much remains hand-set.

Cap the total at 15 directed probes across the whole ladder (Stage 7/9
style, each logged to `LOGBOOK.md` with its reasoning). Do not turn this
into an open search.

## 4. Success criteria

All of the following, at n=128 (primary) and n=40:

**(a) Stage 10 §3 pass/fail table, unchanged.** Sand: a real
growth-dominated burst that genuinely crashes. Vertisol: a genuine,
sustained slow rise with no dip. Use Stage 11's updated Vertisol criteria,
not Stage 10's stale upper bound on `R_peak` (see `STAGE15_RESULTS.md` §3,
which already documents why).

**(b) Recruitment no longer saturates instantly.** Report, per soil:
fraction of finally-recruited regions already active at 2% of the window,
t50, and t_final (all as fractions of the window). Vertisol's baseline is
1.00 / 0.001 / 0.013 — that is the thing being fixed.

**(c) Hotspot selectivity.** Report the percentage of macropore regions
whose peak biomass exceeds 10% of the maximum peak biomass across regions in
that soil. This is Portell et al.'s own metric, so the number is directly
comparable: they reported 9.5 ± 4.0% of colony sites. Vertisol currently
sits near 79% (55 of 70 regions activate). Report it alongside the
growth-weighted "active mass" metric Stage 14 recommended; retire the
boolean `growth > 1e-9` count as a headline number, keeping it only as a
secondary diagnostic.

**(d) Distinctness** between the two soils, time-normalized as in Stages
10–12, still above 0.3.

Do NOT loosen any threshold, and do not re-tune toward a target R(t) shape.

## 5. Evidence already in hand — do not rediscover it

Measured directly on the Stage 11 code, Vertisol n=128 T=1500 and Sand
n=128 T=3000, with `coating_fraction=0` and `k_dis_slow` set equal to
`k_dis=0.05` (i.e. the §1 shared-rule configuration):

| soil | om_total | D_scale | R_peak | B_peak/B0 | frac@2% | t_final | notes |
|---|---|---|---|---|---|---|---|
| Vertisol | 10,000 | 0.7 | 0.97 | 1.00 | 0.48 | 0.381 | staggered, no real growth |
| **Vertisol** | **50,000** | **0.7** | **6.39** | **3.74** | **0.78** | **0.295** | growth AND staggering |
| Vertisol | 200,000 | 0.7 | 19.98 | 14.32 | 0.91 | 0.245 | |
| Vertisol | 50,000 | 5.0 | 42.73 | 37.27 | 0.98 | 0.044 | fast transport kills staggering |
| Vertisol | 50,000 | 20.0 | 80.70 | 69.51 | 1.00 | 0.011 | |
| **Sand** | **3,000** | 80 | **11.94** | **8.07** | — | — | `peak_frac=0.322`, `R_end/R_peak=0.027` — burst-crash PASSES |
| Sand | 20,000 | 80 | 16.81 | 13.59 | — | — | `R_end/R_peak=0.680` — fails to crash |
| Sand | 100,000 | 80 | 19.93 | 13.88 | — | — | `R_end/R_peak=0.472` — fails to crash |

Stage 11 baselines for comparison: Vertisol `R_peak=3.756`,
`B_peak/B0=2.337`, frac@2%=1.00, t_final=0.013. Sand `R_peak=15.790`,
`peak_frac=0.181`, `R_end/R_peak=0.019`.

Two things this evidence establishes, which you should treat as known:

- **The shared rule is viable for both soils.** Sand's burst-crash survives
  removal of the coating at `om_total≈3,000`; Vertisol grows and staggers at
  `om_total≈50,000`.
- **Selectivity and growth trade off against each other**, and slow
  transport is what preserves distance-sensitivity. Across
  `coating_fraction`, `coating_thickness`, `k_dis_slow` and `D_scale`, every
  setting that raises `R_peak` pushes frac@2% toward 1.00. Raising
  `D_scale` destroys staggering fastest.

These are starting points for the §3 ladder, not answers, and they were
measured on single seeds.

## 6. Honest-failure path

If no rung of the §3 ladder satisfies §4 — in particular if the shared rule
gives selectivity but cannot give both soils their §3(a) shapes at any dose
— report that plainly and quantified: which criterion failed, for which
soil, by how much, and at which rung. Name the trade-off in §5 as the
mechanism. **Do not loosen thresholds, do not reintroduce a
distance-dependent decomposition rate, and do not chase a 16th probe.** A
clean negative here is the signal to attempt a physiological rather than
geometric fix (a dormant/active biomass split with substrate-dependent
reactivation, following Blagodatsky & Richter 1998; Wang et al. 2014, *ISME
J*; Chakrawal/Manzoni et al. 2020, *Soil Biol. Biochem.*), which is
explicitly OUT of scope for this stage.

## 7. Optional sensitivity check — explicitly NOT data-driven

If §4 passes and time allows: re-run with OM placed as discrete clumps whose
sizes are drawn from a SHARED lognormal, identical parameters for both
soils, instead of the smooth fine-pore field. Question: does clumping per se
sharpen selectivity beyond what the shared smooth rule already gives?

Keep this bounded and secondary. Do **not** parameterise clump sizes from
the companion project's measured POM data — that data is resolution-
confounded between soils (Mishmar 5.85 µm vs 15 µm elsewhere; see the
companion project's `PROJECT_STATUS.md`, which retired the POM clustering
work as not robust) and does not exist at all for Rehovot/Sand. The model's
structural inputs stay the pore-size distribution and its spatial
arrangement.

## 8. Deliverables

- `dual_porosity.py`: the new `om.mode: shared_fine_pore` placement path and
  uniform-`k_dis` dissolution, with the Stage 6–15 `coating_trapped` path
  preserved behind the flag. Retired keys documented, not deleted.
- `configs/mapping.yaml`: new shared parameters; `stage*_coating_fraction_*`
  and `k_dis_slow` marked RETIRED with a comment pointing here, in the same
  style as the earlier retirements.
- `stage16_run.py`: the §2 geodesic-distance analysis, the §3 ladder with
  its probe trace, the §4 criteria at both grids, and the headline
  comparison figure (both soils' R(t), absolute and normalized).
- `results/stage16/`: geodesic distance distributions + histogram figure per
  soil; `stage16_headline_n{40,128}.csv`; `stage16_pass_fail_n{40,128}.csv`;
  `stage16_recruitment_timing.csv` (frac@2%, t50, t_final, selectivity %);
  `stage16_ladder_trace.csv`; `stage16_distinctness_n{40,128}.csv`.
- `docs/stage_results/STAGE16_RESULTS.md`: the verdict, the geodesic
  distance distributions as the structural finding, which ladder rung was
  adopted and what fraction of the per-soil difference became
  structure-derived, and the selectivity numbers against Portell's 9.5 ± 4.0%.
- A dated correction section appended to `STAGE14_RESULTS.md` (§0 above).
- `MODEL_SPEC.md`, `LOGBOOK.md`, `report_prompt.md` updated as usual.

## 9. Order of work

1. Append the §0 correction to `STAGE14_RESULTS.md`. Do this first, so the
   record is right regardless of how the rest of the stage turns out.
2. Implement the shared-rule OM path; verify the old coating path still
   reproduces Stage 11's recorded numbers behind the flag.
3. §2 geodesic-distance analysis for both soils, before any tuning — it is a
   property of the structure and the shared rule, not of the dose.
4. §3 ladder, rung by rung, ≤15 probes total, each logged.
5. §4 criteria at n=128 and n=40 for the adopted rung.
6. (Optional) §7 shared-lognormal clumping sensitivity check.
7. STOP; write `STAGE16_RESULTS.md` with the numbers and the honest verdict.

The biology stays frozen and identical across soils throughout. Report total
runtime.
