# Claude Code task: Stage 8 — equal total OM, depletion-driven shapes (DEFINED configuration, NOT a search)

Build on Stage 7 (`dual_porosity.py`). This stage is a **defined configuration with
verification**, not an exploration. There is NO search loop hunting for a shape
match, and NO per-soil parameter tuning. Biology frozen; numpy-only; grid n=128
(n=40 check).

## 0. The principle this stage is built on (from Stage 7's diagnosis)

Stage 7 closed the substrate gap — growth now fires — but with abundant carbon every
soil bloomed into the same rising curve, because carbon never depleted. Stage 7's
own conclusion: **the shape difference must come from the depletion TIMING of the
fuel, not from how much carbon arrives.** Equal *total* OM is exactly the lever that
creates that timing difference:

- **Sand** — small porous matrix, one large habitat region: the same total carbon is
  CONCENTRATED on a coating feeding few habitats → high local substrate → a fast
  burst that consumes the local pool faster than the slow trapped release refills it
  → **rise then fall**.
- **Vertisol** — large matrix, ~55–70 habitat regions: the same total carbon is
  SPREAD THIN across many habitats at staggered distances → lower local substrate,
  slower delivery, new regions igniting over time → **sustained rise**.
- **Loess** — intermediate → in between.

The mechanism that produces this (finite local fuel + Monod substrate limitation +
slow refill) is ALREADY in the model. Stage 7 drowned it in excess carbon. This
stage sets the carbon finite and EQUAL so the existing depletion physics can act.

## 1. Change 1 — equal TOTAL OM (absolute, identical across soils)

Replace the fraction-based OM initialization (both the old per-soil fractions and
Stage 7's shared *fraction*, which gave unequal totals because matrix sizes differ)
with a single absolute **`om_total`, identical for all three soils**. Report and
confirm each soil's `om_total` is now exactly equal. The coating/trapped SPLIT still
follows each soil's own geometry (that is the structural difference); only the total
mass is equalised.

## 2. Change 2 — the fuel must be finite and must depend on local depletion

Confirm the OM pool is finite and genuinely depletes where consumed: local substrate
`S` is fed by `k_dis·OM_coat + k_dis_slow·OM_trap` minus uptake, and both OM pools
decrease as they dissolve. `S_cap` stays removed. The fall in Sand must arise
because its concentrated coating pool empties and `f(S)=S/(S+Ks)` drops — verify this
directly (see §4), do not impose it.

## 3. Setting the dose and window — by reasoning, at most THREE directed probes (NOT a scan)

Two numbers must be set: the shared `om_total` and the observation window `T`. Set
them by directed reasoning, each justified by a named diagnostic — NOT by a
parameter sweep:

1. **`om_total`: the LEANEST total that clears growth.** Growth requires initial
   fed-habitat `f_S > ~0.11` (`S ≳ 0.025`). Stage 7 cleared growth at large totals;
   here, lower `om_total` until growth *just* clears in the soils that should grow
   (marginal, not abundant). Lean fuel is what allows depletion to occur within the
   window. Allowed: at most 2–3 probes to bracket this threshold, each reporting the
   fed-habitat `f_S` — not a 15-iteration loop.
2. **`T`: long enough that Sand's concentrated pool depletes, short enough that
   Vertisol is still recruiting.** Check Sand's fuel-remaining and `mean_S_habitat`
   at `t=T`: Sand should have largely spent its accessible coating pool (its
   `mean_S_habitat` has peaked and is falling) while Vertisol is still recruiting new
   regions. If the window doesn't separate them, adjust `T` (the shared experiment
   window) — never per-soil parameters.

Keep everything shared across soils. If three directed probes do not bracket a
working `(om_total, T)`, STOP and report the mechanism (see §5) — do not escalate to
a blind search.

## 4. Keep from Stage 7 (unchanged)

Implicit unified transport; coating/trapped OM placement from geometry; slow trapped
release (`k_dis_slow`); `D_film` water-film macropores; `D_scale`; matric wetting +
inscribed-circle geometry; frozen biology.

## 5. Verification diagnostics (this is what makes the stage principled)

Per soil, plot vs time — these are the direct evidence the mechanism is (or isn't)
acting, not just the final R(t):

- **`mean_S_habitat(t)`** — Sand must RISE THEN FALL (local depletion); Vertisol
  sustained/rising.
- **fuel-remaining(t)** — fraction of each soil's (equal) total OM still undissolved
  — Sand depletes fast, Vertisol slowly.
- **habitats-recruited(t)** — Sand: few, ignite early, no new ones; Vertisol: count
  rises over time.
- **R(t)** vs the measured P1–P4 pattern, and pairwise distinctness.

Success = Sand rise-then-fall, Vertisol rising, Loess intermediate/low, three curves
distinct — achieved with EQUAL total OM, driven by the depletion-timing diagnostics
above (not by any per-soil knob).

## 6. Honest framing — and the explicit stop condition

This is NOT an exploration. If, with equal total OM and a leanest-growth-clearing
dose, the shapes still collapse to "all rising," report it as an honest negative and
name the specific structural quantity that failed to separate the depletion
timescales (e.g. Sand's coating pool did not deplete faster than its refill; or
Vertisol's recruitment did not outlast Sand's depletion). Do NOT respond to a
collapse by scanning more scalars. The point of this stage is to test one clean
hypothesis — equal carbon + structural concentration ⇒ different depletion timing ⇒
different shapes — and to report whether it holds, with the diagnostics as evidence.

## 7. Deliverables

- `dual_porosity.py`: single absolute `om_total` (equal per soil); confirm finite,
  depleting fuel; fuel-remaining and `mean_S_habitat(t)` tracking.
- `configs/mapping.yaml`: `stage8_om_total`, the chosen `T` (fraction-based OM
  retired).
- `stage8_run.py`: the ≤3 directed dose/window probes (logged with their
  diagnostics), the three-soil run, experiment overlay, and the §5 time-series
  diagnostics + n=40 check.
- `results/stage8/`: R(t) per soil, experiment overlay, `mean_S_habitat(t)`,
  fuel-remaining(t), habitats-recruited(t), distinctness CSV, the chosen
  `(om_total, T)` and their justifying diagnostics.
- `docs/stage_results/STAGE8_RESULTS.md`; update `MODEL_SPEC.md`, `LOGBOOK.md`,
  `report_prompt.md`.

## 8. Order of work

1. Set equal absolute `om_total`; confirm per-soil totals identical.
2. Directed probe (≤3) for the leanest growth-clearing `om_total` and a window `T`
   that lets Sand deplete while Vertisol still recruits — each probe justified by
   `f_S`, Sand fuel-remaining, and recruitment diagnostics.
3. Three-soil run at n=128; the §5 diagnostics; n=40 check.
4. Checkpoint: STOP; write `STAGE8_RESULTS.md` — success with the depletion evidence,
   or an honest negative naming the structural quantity that failed.
