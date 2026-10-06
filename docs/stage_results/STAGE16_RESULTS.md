# Stage 16 RESULTS -- the shared OM rule restored; structure sets fuel-habitat geometry

`prompts/soil_respiration_prompt_stage16_shared_om_rule.md`. A MODEL-WIDE
STRUCTURAL change, not a calibration stage and not a search. Runtime 4.6
min (`stage16_run.py`, n=128 primary + n=40, 9 of the allowed 15 probes).

## 0. Result up front

**The shared rule works for everything except one criterion, and the one it
fails, it fails structurally rather than for want of tuning.**

Under ONE shared, structure-blind OM placement rule with spatially uniform
dissolution -- no coating pool, no trapped pool, no per-soil
`coating_fraction`, and no term anywhere whose rate depends on distance to
habitat -- both soils still produce their target shapes, and stay as
distinct from each other as they were with the hand-set per-soil knob:

| SS4 criterion | n=128 | n=40 |
|---|---|---|
| (a) Sand burst-crash, Stage 10's table unchanged | **PASS** (6/6 checks) | **PASS** (6/6) |
| (a) Vertisol slow sustained rise, Stage 11's criteria | **PASS** (5/5 checks) | **PASS** (5/5) |
| (b) recruitment no longer saturates instantly | **FAIL** (frac@2% = 1.000) | **FAIL** (1.000) |
| (c) hotspot selectivity (reported, no threshold set) | Vertisol 21.4% | 44.4% |
| (d) distinctness > 0.3 | **PASS** (2.082) | **PASS** (2.152) |

So SS6's honest-failure path applies: **criterion (b) failed, for
Vertisol, at every rung of the ladder, and the failure is dose- and
transport-independent.** `frac@2% = 1.000` at all nine probes spanning
`om_total` 3,000-50,000 (17x) and `D_scale` 0.7-80 (114x). No threshold was
loosened, no distance-dependent decomposition rate was reintroduced, and
the probe budget was not chased -- 9 of 15 probes were used, and the tenth
would have told us nothing the first nine did not.

**Why (b) cannot be fixed by dose is the same finding as SS2's**, seen from
the other side. Under a shared rule that fills the whole porous matrix, 66
of Vertisol's 70 macropore regions are *physically adjacent* to a
fuel-bearing cell (median geodesic distance 20 um = 2 cells). Activation in
`count_active_habitats` requires only `growth > 1e-9`; every region has fuel
touching it, so every region crosses that threshold immediately. Dose scales
how *much* a region grows, never *whether* it starts. The boolean count is
saturated by construction of the pore field itself.

**Two corrections to the prompt's own premises, both load-bearing and both
verified directly:**

1. **SS5's evidence does not transfer to SS1's rule** (Section 5 below).
   SS5's table was measured with `coating_fraction=0`, described there as
   "i.e. the SS1 shared-rule configuration". It is not: with
   `coating_fraction=0` the coating POOL is empty but the geodesic
   construction still runs, placing all OM on `trapped_mask` -- the matrix
   *minus* the near-habitat shell. That punches a fuel-free exclusion zone
   around every region (0 of 70 regions adjacent to fuel, versus 66 of 70
   under the real shared rule), and that exclusion zone is the entire source
   of the staggering SS5 reports. Reproduced exactly: SS5's `R_peak=6.39 /
   frac@2%=0.78 / t_final=0.295` comes back as `6.386 / 0.778 / 0.2949`.
   SS5's proxy is itself a distance-dependent OM placement -- the very thing
   SS1 set out to remove.
2. **The "79%" selectivity baseline SS4(c) quotes is a different metric
   from the one it asks for.** 79% (55 of 70) is the BOOLEAN recruited
   count. On Portell's actual metric -- peak biomass above 10% of the
   largest region's peak -- Stage 11's baseline Vertisol is **14.3%**, and
   the shared rule gives **21.4%**. Both sit near Portell's 9.5 +/- 4.0%;
   selectivity did not improve, it got modestly worse.

## 1. The model change, and the preserved path (SS1, SS9.2)

`dual_porosity.py` gained `om_mode`, now `shared_fine_pore` by default in
`configs/mapping.yaml`:

- `place_om_shared_fine_pore` places ALL organic matter by the model's
  ORIGINAL rule (`MODEL_SPEC.md` SS4): `OM ~ d^-om_fine_bias` within the
  porous matrix, normalized to that soil's total. It never measures a
  cell's distance to a habitat.
- Dissolution is spatially uniform. The slow pool is identically zero under
  this mode, and `simulate_dual_porosity` additionally ties `k_dis_slow` to
  `k_dis` there, so a second spatially-distinct decomposition rate cannot
  re-enter through that key by a later edit.
- `coating_fraction` is retired as an active parameter and reported as
  `nan` under the shared mode rather than as a number that no longer means
  anything.

The Stage 6-15 path is preserved, not deleted, behind `om_mode:
coating_trapped` -- the same treatment Stage 13 gave `psd_parametric`
alongside `psd_measured`. Every Stage 4-15 runner pins that mode explicitly
on load, so no earlier stage silently inherits the new rule. Verified, not
assumed (`stage16_legacy_reproducibility.csv`):

| soil | R_peak now | Stage 11 reference | deviation |
|---|---|---|---|
| Sand | 15.7899 | 15.790 | **0.001%** |
| Vertisol | 3.7561 | 3.756 | **0.004%** |

Stage 16's two new observational outputs (`active_mass_t` always;
`region_biomass_t` opt-in) were checked the way Stage 14 checked its own:
`R_t`, `B_total_t`, `habitats_recruited_t` and `fuel_remaining_t` are
byte-identical (`np.array_equal`) with the instrumentation on and off.

## 2. THE STRUCTURAL RESULT: emergent fuel-habitat distance (SS2)

Definition used, stated once and applied consistently to every headline
number: **"substantial" OM = a porous-matrix cell whose OM mass is at or
above the MEDIAN of all nonzero matrix OM in that soil** (SS2's own
suggested definition). It is scale-free -- the whole field scales by one
factor when `om_total` changes, so the mask and every distance below are
invariant to the dose. Verified directly, not assumed: identical medians at
`om_total` = 1,000 / 50,000 / 1,000,000
(`stage16_geodesic_distance.csv`, `dose_invariant_median=True`).

Distance is measured by `bfs_geodesic_distance` outward from substantial OM
through everything that is not solid, and each region takes the minimum over
its own cells. `stage16_geodesic_distance.csv`, in um (voxel = 10 um):

| soil | grid | regions | unreachable | min | p25 | median | p75 | p95 | max | mean +/- sd |
|---|---|---|---|---|---|---|---|---|---|---|
| Sand | n=128 | 3 | 2 | 10 | 10 | **10** | 10 | 10 | 10 | 10.0 +/- 0.0 |
| Sand | n=40 | 1 | 0 | 10 | 10 | **10** | 10 | 10 | 10 | 10.0 +/- 0.0 |
| Vertisol | n=128 | 70 | 4 | 10 | 20 | **20** | 20 | 40 | 50 | 22.4 +/- 7.6 |
| Vertisol | n=40 | 9 | 0 | 10 | 20 | **20** | 20 | 26 | 30 | 20.0 +/- 4.7 |

**Reported honestly, as SS2 requires: the two distributions come out far
more similar than expected, and far smaller.** Sand is 1 cell, flat, with
zero spread; Vertisol is 2 cells with a thin tail to 5. One giant macropore
region in an 18%-matrix soil and 70 dispersed regions in a 68%-matrix soil
produce fuel-habitat distances that differ by a single grid cell in the
median. This is the stage's central negative result, and it is the direct
cause of the SS4(b) failure in Section 4.

**Against Portell et al. (2018)**, whose operative geodesic distances were
252 +/- 111 um and 490 +/- 262 um: this model's emergent distances are
**10-50 um, roughly 5-25x smaller**, and do not land in a comparable range.
The domain-size difference is real and stated plainly -- Portell's domain
was ~6.8 mm against this model's 1.28 mm at n=128 -- but it does not explain
the gap: 252 um would be 25 cells here, well inside a 128-cell domain, so
this model *could* have produced distances in their range and did not.

### 2b. Why the median cut is so permissive (secondary, not a redefinition)

The median is this stage's stated definition and every headline number above
uses it. It is also, under this model's `d^-2.5` rule, a very weak
discriminator: Stage 3.5 measured ~99.96% of carbon landing on sub-micron
cells, and here the top half of matrix cells by OM hold **100.0%** of
Vertisol's mass and **97.0%** of Sand's. "Nearest substantial OM" therefore
degenerates to "nearest matrix cell".

`stage16_geodesic_threshold_sensitivity.csv` sweeps how strict "substantial"
is (n=128, median distance in um). This is reported because it is where the
two soils genuinely separate -- it is NOT used for any pass/fail number, and
no threshold was chosen to match an external reference:

| percentile cut | Sand: share of OM | Sand median | Vertisol: share of OM | Vertisol median | Vertisol p95 |
|---|---|---|---|---|---|
| 50 (**stated definition**) | 97.0% | **10** | 100.0% | **20** | 40 |
| 75 | 85.7% | 10 | 98.4% | 30 | 116 |
| 90 | 64.3% | 20 | 90.4% | 40 | 150 |
| 95 | 47.6% | 20 | 80.4% | 60 | 276 |
| 99 | 18.8% | 20 | 52.0% | 200 | 478 |
| 99.9 | 4.1% | 40 | 18.9% | 370 | 750 |

Sand saturates at 20-40 um no matter how strict the cut -- its single giant
region touches everything, so it is always a couple of cells from some fuel.
Vertisol climbs to 200-370 um median and 478-750 um p95, i.e. squarely into
Portell's range, because its 70 dispersed regions really do sit at varied
distances from the large carbon concentrations. **The structural difference
between the soils is real and is visible in the tail of the OM distribution;
it is invisible at the median cut, which is the cut this stage committed
to.** Both statements are reported rather than one being chosen.

## 3. The ladder (SS3): 9 of 15 probes, no rung adopted

`stage16_ladder_trace.csv`. Every probe evaluated both soils at BOTH grids.
`ok` = that soil's full SS4(a) check set.

| probe | rung | Sand om / D | Sand ok | Vertisol om / D | Vert ok | frac@2% | t_final | distinct (n128) |
|---|---|---|---|---|---|---|---|---|
| 1 | 1 | 3,000 / 80 | **True** | 3,000 / 80 | False | 1.00 | 0.003 | 0.553 |
| 2 | 1 | 50,000 / 0.7 | False | 50,000 / 0.7 | **True** | 1.00 | 0.015 | 0.434 |
| 3 | 1 | 12,247 / 7.5 | False | 12,247 / 7.5 | **True** | 1.00 | 0.012 | 0.398 |
| 4 | 1 | 12,247 / 80 | False | 12,247 / 80 | False | 1.00 | 0.001 | 1.021 |
| 5 | 2 | 3,000 / 80 | **True** | 11,083 / 80 | False | 1.00 | 0.001 | 0.382 |
| 6 | 2 | 13,534 / 0.7 | False | 50,000 / 0.7 | **True** | 1.00 | 0.015 | 0.305 |
| 7 | 2 | 6,372 / 7.5 | False | 23,540 / 7.5 | False | 1.00 | 0.010 | 0.242 |
| 8 | 2 | 6,372 / 10.8 | False | 23,540 / 4.25 | **True** | 1.00 | 0.015 | 0.564 |
| 9 | 3 | 3,000 / 80 | **True** | 50,000 / 0.7 | **True** | 1.00 | 0.015 | **2.082** |

**Rung 1 (one shared `om_total`, one shared `D_scale`) fails, and the
prompt's expectation that it would is confirmed with numbers.** No single
shared setting gives both soils their shapes: at Sand's working point
(probe 1) Vertisol is starved; at Vertisol's (probe 2) Sand neither crashes
nor bursts; at the geometric midpoint (probe 3) Sand still fails. Probe 4
isolates the axes and confirms SS5's warning -- raising `D_scale` to 80 at
the midpoint dose breaks *both* soils and drives `t_final` to 0.001, the
fastest saturation in the whole ladder.

**Rung 2 (shared concentration per unit matrix volume) closes part of the
gap, and here is exactly how much.** Measured porous-matrix cell counts at
n=128 are 1,953 (Sand) and 7,215 (Vertisol) -- a **3.694x** ratio (identical
to the full-matrix ratio 0.6775/0.1834, since the solid share is a fixed
fraction of the matrix for every soil, so the choice of which matrix volume
to use redefines rho's units but never the rung's behavior). The per-soil
dose difference actually required is 50,000/3,000 = **16.67x**. So deriving
the total from structure accounts for:

- **22.2% of the required ratio** on a linear reading (3.694 / 16.67), or
- **46.5% of it in orders of magnitude** (log 3.694 / log 16.67) -- it
  closes a little under half the gap on a log scale.

Either way the remaining factor of **4.5x** is not structural and stays
hand-set. Probes 5-7 confirm this behaviorally: pinning rho so one soil
lands on its working dose leaves the other clearly off it.

**The structure-derived `D_scale` option SS3 asks to consider, reported.**
Probe 8 sets one shared coefficient and gives each soil `D_scale = coeff x
its own macropore fraction` (Sand 0.8166, Vertisol 0.3225), so transport too
is derived rather than set; the coefficient is fixed so the two soils' mean
matches probe 7's shared `D_scale`, making probes 7 and 8 differ *only* in
whether transport is structure-derived. The macroporosity ratio is
**2.53x** against the **114x** (80/0.7) actually required -- **2.2%
linearly, 19.6% in orders of magnitude.** It is the weaker of the two
structural levers by a wide margin. It does help: probe 8 recovers
Vertisol's shape where probe 7 lost it, and raises distinctness from 0.242
to 0.564. It does not close the gap.

**Rung 3 (per-soil, Stage 10's position) gives both shapes and by far the
best distinctness (2.082), but still fails (b)** -- so the ladder never
reaches a rung that passes, and **nothing is adopted**. Probe 9 is reported
below as the closest attempt, by a selection rule fixed in advance (most
SS4 criteria met, ties broken by distinctness) and reusing the ladder's own
run rather than spending a new probe.

## 4. SS4 criteria in full, for the closest attempt (probe 9)

`stage16_headline_n{40,128}.csv`, `stage16_pass_fail_n{40,128}.csv`,
`stage16_recruitment_timing.csv`, `stage16_distinctness_n{40,128}.csv`.

**(a) Shapes -- both soils, both grids, all checks pass, no threshold moved.**

| soil | grid | R_peak | R_end/R_peak | peak @ | B_peak/B0 | checks |
|---|---|---|---|---|---|---|
| Sand | 128 | 13.994 | 0.008 | 0.234 | 9.076 | **6/6** |
| Sand | 40 | 32.541 | 0.001 | 0.065 | 17.338 | **6/6** |
| Vertisol | 128 | 9.129 | 1.000 | 1.000 | 5.756 | **5/5** |
| Vertisol | 40 | 4.415 | 1.000 | 1.000 | 3.936 | **5/5** |

Sand bursts and genuinely crashes; Vertisol rises through the whole window
with no dip. Both clear the maintenance floor with growth-dominated
respiration, and both grids agree -- which Stage 15 found was not
guaranteed.

**(b) Recruitment timing -- THE FAILURE.** Operationalized in advance (SS4b
states no numeric bar): `frac@2% < 0.90` and `t_final > 0.10`.

| soil | grid | frac@2% | t50 | t_final | passes? |
|---|---|---|---|---|---|
| Vertisol | 128 | **1.000** | 0.0013 | **0.0147** | **No** |
| Vertisol | 40 | **1.000** | 0.0007 | 0.0007 | **No** |
| Sand | 128 | 1.000 | 0.0003 | 0.0003 | No (3 regions; metric near-meaningless) |
| Sand | 40 | 1.000 | 0.0003 | 0.0003 | No (1 region; meaningless) |

Vertisol's baseline was 1.00 / 0.001 / 0.013. It is now 1.000 / 0.0013 /
0.0147. **Essentially unchanged -- the thing this stage set out to fix did
not move**, and Section 0 explains why it cannot: 66 of 70 regions touch
fuel directly.

**But the metric that failed is the one Stage 14 showed is not to be
trusted, and the graded metric tells the opposite story.** Per-region
time-to-half-max (Stage 14's own graded metric, reused unchanged), as a
fraction of the window:

| configuration | p25 | median | p75 | max |
|---|---|---|---|---|
| Vertisol, Stage 11 baseline (coating) | 0.335 | 0.440 | 0.506 | 0.917 |
| **Vertisol, shared rule (n=128)** | **0.482** | **0.534** | **0.690** | **0.967** |
| Vertisol, shared rule (n=40) | 0.279 | 0.387 | 0.407 | 0.697 |

Measured gradually rather than as a boolean, recruitment under the shared
rule is *more* spread across the window than the baseline, not less. There
is a genuine tension inside SS4 here, reported rather than resolved in this
stage's favor: **(b) is defined on the boolean count that (c) instructs us
to retire as a headline number.** On the boolean count the stage fails; on
the graded instrument Stage 14 recommended it does not. This write-up
counts (b) as FAILED, because that is what the criterion as written says,
and flags the tension rather than picking the flattering reading.

**(c) Hotspot selectivity (reported; SS4c sets no threshold).** Portell et
al.'s own metric -- share of regions whose peak biomass exceeds 10% of the
largest region's peak:

| soil | grid | hotspot % | regions | Portell | secondary: boolean recruited % |
|---|---|---|---|---|---|
| Vertisol | 128 | **21.4%** | 15/70 | 9.5 +/- 4.0% | 78.6% |
| Vertisol | 40 | 44.4% | 4/9 | | 88.9% |
| Sand | 128 | 33.3% | 1/3 | | 33.3% |
| Sand | 40 | 100.0% | 1/1 | | 100.0% |

Stage 11's baseline Vertisol on this same metric is **14.3%**. So the shared
rule moved selectivity from 14.3% to 21.4% -- **slightly further from
Portell's 9.5 +/- 4.0%, not closer.** The prompt's stated baseline of "near
79%" is the boolean recruited count (78.6%), a different quantity; comparing
it against Portell's 9.5% is not like-for-like, which is exactly the
confusion Stage 14 warned the boolean count creates. Sand's numbers are
reported for completeness but carry no information: with 1-3 regions, and
`B0 ~ K` seeding biomass in proportion to region size, the metric is
dominated by region count and size rather than by selectivity.

**(d) Distinctness: PASS, comfortably, at both grids** -- 2.082 (n=128) and
2.152 (n=40) against the 0.3 threshold, nearly 7x the bar. Against Stage
15's frozen coating configuration for the same pair (2.001 at n=128, 2.217
at n=40), the shared rule is slightly MORE distinct at n=128 and slightly
LESS at n=40. So the honest statement is that retiring the per-soil carbon
knob left between-soil distinctness essentially unchanged -- it did not
improve it, and, importantly for this stage's argument, it did not cost
anything either.

## 5. Why SS5's evidence did not transfer (`stage16_ss5_reconciliation.csv`)

A diagnostic, not a probe: it adopts nothing and proposes nothing. Vertisol,
n=128, T=1500, `om_total=50,000`, `D_scale=0.7` -- the same point in both
cases, differing only in how OM is placed:

| | SS5's proxy (`coating_fraction=0`) | SS1's actual shared rule |
|---|---|---|
| OM-bearing matrix cells | 2,946 / 7,215 | **7,215 / 7,215** |
| regions adjacent to fuel | **0 / 70** | **66 / 70** |
| median geodesic distance | 40 um | 20 um |
| R_peak | 6.386 (SS5 records 6.39) | 9.129 |
| B_peak/B0 | 3.742 (SS5 records 3.74) | 5.756 |
| frac@2% | **0.778** (SS5 records 0.78) | **1.000** |
| t_final | **0.2949** (SS5 records 0.295) | **0.0147** |

SS5's numbers are real and reproduce to three decimals -- but they were
measured on a configuration that still runs the geodesic construction and
leaves 4,269 near-habitat matrix cells carrying **zero** OM. That exclusion
zone is a distance-dependent OM placement, and it is the whole of the
staggering. SS5's conclusion that "the shared rule is viable for both soils"
holds for growth and for Sand's burst-crash; its implication that Vertisol
would also *stagger* under the shared rule does not, and was never tested by
the configuration that produced the table.

## 6. Honest verdict (SS6)

**Which criterion failed, for which soil, by how much, at which rung:**
SS4(b), recruitment timing, for Vertisol (and trivially for Sand), at
**every** rung -- rungs 1, 2 and 3 alike. `frac@2%` is **1.000 against a
0.90 bar** and `t_final` is **0.0147 against a 0.10 bar**, i.e. off by a
factor of ~7 on the criterion that has any resolution, and identical to the
1.00 / 0.013 baseline the stage was asked to fix.

**The mechanism, named.** Not the trade-off SS5 anticipated. SS5 predicted
selectivity and growth would trade off, with slow transport preserving
distance-sensitivity; that trade-off is visible in the ladder (probe 4's
fast shared transport breaks both soils) but it is not what defeats (b).
What defeats (b) is SS2's own result: **under any shared rule that fills the
porous matrix, fuel is adjacent to essentially every macropore region (66 of
70; median distance 2 cells), and boolean activation needs only `growth >
1e-9`.** Recruitment is saturated by the geometry of the pore field before
any dose or transport parameter is chosen. This is why it is a structural
result and not a search-budget failure -- six probes were left unspent
because no further probe on these axes could change it.

**What the stage did establish.** Removing the per-soil carbon knob costs
nothing on the criteria that measure the project's central claim: both
soils keep their shapes at both grids under one shared, structure-blind
rule, and between-soil distinctness is essentially unchanged from the
hand-set per-soil coating (2.082 vs 2.001 at n=128; 2.152 vs 2.217 at
n=40). The per-soil carbon knob was buying nothing that structure was not
already providing, so the project's central claim is genuinely stronger
without `coating_fraction`. What the shared rule does not deliver is a
fuel-habitat distance distribution that differs between soils at the stated
median cut -- and, consequently, staggered boolean recruitment.

**Per SS6, this is the signal to attempt a physiological rather than
geometric fix** -- a dormant/active biomass split with substrate-dependent
reactivation (Blagodatsky & Richter 1998; Wang et al. 2014, *ISME J*;
Chakrawal/Manzoni et al. 2020, *Soil Biol. Biochem.*), which is explicitly
out of scope here and was not attempted. A dormancy threshold would make
activation depend on how much substrate arrives rather than on whether any
arrives, which is precisely the degree of freedom the boolean `growth >
1e-9` test currently lacks.

## 7. Clumping, run as an SS6 mechanism diagnostic (not SS7's check)

SS7 gates the clumping sensitivity check on SS4 passing, and SS4 did not
pass. It was run anyway, strictly to answer the question SS6 requires --
does the failure belong to *shared rules*, or only to *smooth* shared
rules? Clumping is the one placement change that stays shared and
structure-blind (one `n_clumps`, one lognormal, identical for both soils;
deliberately NOT parameterised from the companion project's retired,
resolution-confounded POM data) while still opening fuel-free gaps.
**Nothing here is adopted, no threshold moved, and it is not counted against
the SS3 probe cap.** `stage16_clumping_sensitivity.csv`:

| soil | grid | frac@2% smooth -> clumped | t_final smooth -> clumped | hotspot % | shape still ok? |
|---|---|---|---|---|---|
| Vertisol | 128 | 1.000 -> **0.739** | 0.0147 -> **0.3582** | 21.4 -> 20.0 | **True** |
| Vertisol | 40 | 1.000 -> 1.000 | 0.0007 -> 0.0007 | 44.4 -> 44.4 | **False** |
| Sand | 128 | 1.000 -> 1.000 | 0.0003 -> 0.0003 | 33.3 -> 33.3 | True |
| Sand | 40 | 1.000 -> 1.000 | 0.0003 -> 0.0003 | 100.0 -> 100.0 | True |

**The answer is that the failure belongs to smooth shared rules, not to
shared rules as such** -- at n=128, clumping alone restores staggering
(`frac@2%` 0.739 and `t_final` 0.358 would both clear the bars) while
keeping Vertisol's shape, purely by concentrating carbon so that gaps open
around some regions. It does **not** reproduce at n=40, where Vertisol's
shape check fails outright and timing does not move; with 9 regions and 704
matrix cells the clump field is too coarse to open gaps. Answering SS7's
actual question -- does clumping sharpen selectivity? -- **no**: 21.4% ->
20.0%, a negligible change.

This is reported as a lead, not a result. A one-grid effect that fails at
the other grid is exactly the kind of finding Stage 15 was burned by, and
it is not offered as a rescue of SS4.

## 8. Deliverables

- `dual_porosity.py`: `om_mode` (`shared_fine_pore` default /
  `coating_trapped` preserved / `shared_lognormal_clumps` for the SS7
  question), `place_om_shared_fine_pore`, `place_om_shared_lognormal_clumps`,
  `substantial_om_mask`, `analyze_fuel_habitat_distance`, the Stage 16 ladder
  tiers in `_variantc_om_total`/`_variantc_D_scale`, uniform-`k_dis`
  enforcement, and the observational `active_mass_t` / `region_biomass_t`
  outputs.
- `configs/mapping.yaml`: `om_mode` as the new active default, the SS3 ladder
  keys documented, and `k_dis_slow` plus all eight `*_coating_fraction_*`
  keys marked RETIRED (retained so the preserved path still runs).
- Stage 4-15 runners: each pins `om_mode: coating_trapped` on load.
- `stage16_run.py`, `results/stage16/` (13 CSVs + 3 figures), this file, and
  the dated correction appended to `STAGE14_RESULTS.md`.

**[Stage 16] Checkpoint reached: the per-soil coating/trapped OM split is
retired in favor of one shared, structure-blind rule with uniform
dissolution; the Stage 6-15 path is preserved and reproduces Stage 11 to
within 0.004%. Stage 14's causal attribution refuted and corrected in
place. SS4(a) shapes and SS4(d) distinctness PASS at both grids (2.082 /
2.152, the best recorded for this pair); SS4(b) recruitment timing FAILS at
every rung, structurally and dose-independently, because a shared rule
leaves 66 of 70 regions adjacent to fuel; distinctness is essentially
unchanged from Stage 15's per-soil coating, so retiring the knob cost
nothing. No rung of the SS3 ladder adopted;
9 of 15 probes used. Rung 2 accounts for 22.2% (linear) / 46.5% (log) of the
required per-soil dose difference and 2.2% / 19.6% of the transport
difference. Thresholds not loosened, no distance-dependent decomposition
rate reintroduced, no 16th probe chased. Stopping here per
soil_respiration_prompt_stage16_shared_om_rule.md SS6/SS9.**
