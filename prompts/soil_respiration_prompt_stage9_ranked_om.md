# Claude Code task: Stage 9 — ranked per-soil OM + placement (DEFINED calibration, NOT a search)

Build on Stage 8 (`dual_porosity.py`). Defined configuration with verification, like
Stage 8: NO open-ended search, `(scale + per-soil coating fractions)` set by at most
3 directed diagnostic probes. Biology frozen; numpy-only; grid n=128 (n=40 check).
Keep everything from Stage 8 (depletion mechanism, implicit unified transport,
coating/trapped OM with slow trapped release `k_dis_slow`, `D_film`, `D_scale`,
matric wetting + inscribed-circle geometry, finite/depleting fuel).

## 0. Why this stage

Stage 8 proved the depletion mechanism gives genuinely distinct curves
(`min_distance=0.80`) and reproduced Loess-erratic, but showed that a single EQUAL
total OM cannot give Sand-fall AND Vertisol-rise together: Sand needs a lean pool
(so it depletes → falls), Vertisol needs more accessible carbon (so its many
regions sustain → rises) — opposite regimes. This stage resolves that with three
targeted changes you specified. **The OM amounts are set ONLY by a qualitative
ranking, deliberately NOT fitted to the measured respiration values.**

## 1. Change 1 — per-soil total OM, ranked Vertisol > Loess > Sand

Give each soil a different total OM, ordered **Vertisol > Loess > Sand**. Implement
as a shared base × a per-soil multiplier that enforces this ranking (e.g. start
Vertisol 1.5 / Loess 1.0 / Sand 0.6 × base). The **ordering is fixed**; only the
base scale (and, if needed, the spread) is set by the directed probes (§4). These
amounts are ranked-only — do NOT derive them from the experiment's P1–P4 numbers or
any measured OM content; the ranking is the only real-world input. Report each
soil's `om_total` and confirm the ranking holds.

## 2. Change 2 — make Vertisol's OM more available

Raise Vertisol's **coating fraction** so more of its carbon sits on the fast,
habitat-adjacent coating (reachable) rather than trapped in the interior — i.e. more
of Vertisol's (now larger) carbon actually reaches its many habitat regions, letting
their staggered recruitment sustain a rise. (Availability = coating fraction + being
placed on the connected coating path.) Vertisol coating fraction goes UP from its
Stage-8 value.

## 3. Change 3 — Sand: less coating, more matrix

Lower Sand's **coating fraction**; move the remainder into the trapped/matrix pool
(which releases slowly via `k_dis_slow`). So Sand gets a smaller fast burst plus a
slow low tail, and — combined with its lowest total (Change 1) — its accessible pool
depletes within the window → a clean rise-then-fall. Sand coating fraction goes DOWN
from its Stage-8 0.8.

**Exception — the burst wins.** Sand's rise-then-fall (a clear early peak that then
declines) is the PRIORITY; the coating reduction is allowed only as far as it
preserves that burst. If lowering Sand's coating fraction (or its total OM) BREAKS
the burst — no clear early peak, or Sand degrades to a slow monotonic rise or a low
flat curve with no peak-then-decline — then **restore OM to Sand's coating**: raise
its coating fraction back up (and, if needed, its total) until the burst reappears.
Push "less coating" only up to the point where the rise-then-fall still holds, and
back off the moment it breaks. Report the coating fraction at which the burst broke
and the value finally adopted.

(Loess stays intermediate in both total and coating fraction.)

## 4. Directed calibration (≤3 probes, NOT a search)

Two things to set: the shared base OM scale and the per-soil coating fractions
(within the directions above). Set them by at most 3 directed probes, each justified
by named diagnostics — NOT a parameter scan, NOT per-soil shape-fitting beyond the
specified directions:

- **Vertisol check:** its fuel is NOT fully exhausted by `t=T` (fuel_remaining stays
  well above 0), its `habitats_recruited(t)` keeps climbing, and `R(t)` rises through
  the window.
- **Sand check:** its accessible (coating) pool depletes, `mean_S_habitat(t)` peaks
  then falls, and `R(t)` shows a real rise-then-fall. If the burst is broken (no
  early peak), apply the §3 Exception and restore coating before spending another
  probe — preserving Sand's burst takes precedence over minimizing its coating.
- **Loess check:** intermediate / erratic.
- **Distinctness** stays > 0.3.

Keep the ranking V > L > S fixed throughout. If 3 directed probes do not land the
three shapes, STOP and report which quantity fell short (§6) — do NOT escalate to a
scan.

## 5. Run, compare, diagnose

Three soils, n=128, plus the n=40 check. Report per soil: `om_total` (confirm
ranking), coating fraction, `mean_S_habitat(t)`, `fuel_remaining(t)`,
`habitats_recruited(t)`, outcome, R(t) shape and R_end/R_peak, and pairwise
distinctness. Headline figure: R(t) for the three soils vs the measured P1–P4 shapes.

**Grid note:** Stage 8 found Sand's fall reads much more sharply at n=40 than n=128
(aggregate-averaging effect). Report Sand's shape at BOTH grids and state plainly
which grid the rise-then-fall is clear at; do not tune to hide the grid sensitivity.

## 6. Honest framing and stop condition

Report whether the three shapes now match (Sand rise-then-fall, Vertisol rising,
Loess erratic) AND stay distinct. State explicitly that the OM amounts were set by
ranking only, not fitted to the experimental values — so a shape match is a
structural/qualitative result, not a fit. If a shape still fails, name the specific
quantity (e.g. Vertisol's recruitment still runs out of fuel; Sand's burst too small
to register a fall at n=128) — do not chase with more probes.

## 7. Deliverables

- `dual_porosity.py`: per-soil ranked `om_total` (base × ranked multiplier); per-soil
  coating fractions (Vertisol up, Sand down); keep fuel_remaining / mean_S_habitat /
  recruitment tracking.
- `configs/mapping.yaml`: the base OM scale, per-soil OM multipliers (ranked), per-soil
  coating fractions, chosen `T`.
- `stage9_run.py`: the ≤3 directed probes (logged with diagnostics), the three-soil
  n=128 run + experiment overlay + the §5 diagnostics + n=40 check.
- `results/stage9/`: R(t) per soil, experiment overlay, mean_S_habitat(t),
  fuel_remaining(t), habitats_recruited(t), distinctness CSV, chosen parameters.
- `docs/stage_results/STAGE9_RESULTS.md`; update `MODEL_SPEC.md`, `LOGBOOK.md`,
  `report_prompt.md`.

## 8. Order of work

1. Ranked per-soil `om_total` (V>L>S) + Vertisol coating up + Sand coating down;
   confirm totals and ranking.
2. ≤3 directed probes on base scale + coating fractions, each justified by the §4
   Vertisol/Sand/Loess diagnostics.
3. Three-soil n=128 run + experiment overlay + n=40 check.
4. Checkpoint: STOP; write `STAGE9_RESULTS.md` — the shape match with depletion/
   recruitment evidence, or an honest negative naming the quantity that fell short.
