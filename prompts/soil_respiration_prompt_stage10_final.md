# Claude Code task: Stage 10 (FINAL) — Sand real burst-then-crash, Vertisol slow rise

This is the LAST stage. Simplify to **TWO soils only: Vertisol and Sand.** Drop
Loess entirely (do not run or report it). Build on Stage 9 (`dual_porosity.py`),
keep its mechanism. numpy-only; grid n=128 (n=40 check). Whole task: produce the two
target respiration curves below and verify them. Per-soil parameters MAY be tuned to
hit these two shapes (the structure-only purity of earlier stages is set aside for
this final deliverable) — but the mechanism must stay the physical one described in
§2, not arbitrary curve-drawing.

## 1. The two required outcomes (this IS the task)

- **Sand — a REAL early burst, then a crash.** Growth must genuinely fire: `R(t)`
  rises to an early peak that is clearly ABOVE the maintenance floor (`m0·B0 = 0.4`),
  driven by growth-associated respiration (not just maintenance of the seeded
  biomass), then **crashes** back down as the fuel runs out. A monotonic decay from
  the maintenance floor (what Stage 9 produced: `R_peak ≈ 0.4`, biomass crashing to
  5%) does NOT count — that is starvation, not a burst.
- **Vertisol — a slow rise.** `R(t)` climbs gradually across the window and is still
  rising / peaks late. Not a front-loaded early spike that then dips.

## 2. How to get there — the physical contrast (per-soil allowed)

- **Sand = concentrated, fast, finite.** Put most of its carbon on the coating (high
  coating fraction) with enough TOTAL that growth clears strongly and early; keep the
  slow trapped release small so that once the coating pool depletes, respiration
  crashes rather than being propped up. Fast delivery (coating adjacent + adequate
  transport). The burst comes from growth on the coating carbon; the crash from that
  finite pool depleting.
- **Vertisol = distributed, slow, sustained.** Spread its carbon (lower coating
  fraction, more in the matrix), delivered gradually so recruitment is staggered over
  the window → a slow rise. If the shared fast transport (`D_scale`) front-loads
  Vertisol into an early spike, **slow Vertisol's effective delivery** (a lower
  per-soil transport magnitude and/or carbon placed farther from habitat) so its rise
  is genuinely gradual. Give it enough total to keep rising through the window.

## 3. Pass/fail criteria (must verify, report the numbers)

- **Sand burst is real:** `R_peak ≥ 2 × maintenance_floor` (i.e. ≥ ~0.8, and clearly
  higher is better), the growth-associated term is the MAJORITY of `R` at the peak
  (report growth-vs-maintenance split at peak), biomass grows (`B_peak/B0 > 1`), and
  the coating fuel visibly depletes.
- **Sand crash:** `R_end / R_peak < 0.3`, with the peak in the first ~third.
- **Vertisol slow rise:** late-third mean `>` early-third mean, the curve is
  rising/peaking in the LATE part of the window (not an early spike + dip), and it is
  clearly not a burst.
- **Distinct:** the two normalized curves are distinct (> 0.3).

If any criterion is missed, adjust the per-soil knobs toward it (this is the final
deliverable — make a genuine effort to hit both), and if one still cannot be met,
STOP and report exactly which criterion and the quantity that blocked it.

## 4. Deliverables

- `dual_porosity.py` / `stage10_run.py`: the two-soil (Vertisol, Sand) run; per-soil
  parameters that achieve §1; the growth-vs-maintenance split diagnostic at Sand's
  peak; `mean_S_habitat(t)`, `fuel_remaining(t)`, biomass(t) for both soils.
- `results/stage10/`: the two-soil `R(t)` figure vs the measured shapes (Sand
  falling, Vertisol rising), the Sand burst-proof diagnostics, the Vertisol slow-rise
  diagnostics, distinctness, n=40 check, chosen per-soil parameters.
- `docs/stage_results/STAGE10_RESULTS.md` — the two curves, the pass/fail table
  against §3, and the final per-soil parameters. Update `MODEL_SPEC.md`, `LOGBOOK.md`,
  `report_prompt.md`. Note this is the final stage: the report should read as the
  consolidated two-soil demonstration.

## 5. Order of work

1. Restrict to Vertisol + Sand.
2. Set Sand for a real growth burst then crash (§2, §3); verify `R_peak` well above
   the maintenance floor with a growth-dominated peak, then a crash.
3. Set Vertisol for a slow, gradual rise (§2, §3); slow its delivery if it
   front-loads.
4. Two-soil run n=128 + n=40 check; the §3 pass/fail table; the comparison figure.
5. STOP; write `STAGE10_RESULTS.md`.

The entire task is these two curves, verified: Sand real-burst-then-crash, Vertisol
slow rise.
