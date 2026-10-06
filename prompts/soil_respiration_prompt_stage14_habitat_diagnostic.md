# Claude Code task: Stage 14 — investigating the habitats_recruited(t) saturation artifact

A DIAGNOSTIC stage, not a retuning stage. Stage 11 reported, as an explicitly
unresolved limitation: `habitats_recruited(t)` reaches its final count within the
first ~2% of the observation window and stays completely flat afterward, for
Vertisol, regardless of tuning — even though `R(t)` itself keeps genuinely rising.
The suspected cause (never confirmed) is that the implicit backward-Euler solver's
per-step equilibration crosses the "recruited" activation threshold
(`growth > 1e-9`) almost everywhere almost instantly, independent of transport
speed. This stage confirms or refutes that mechanism and, only if a low-risk fix
exists, proposes one — without changing Stage 10/11/13's frozen calibrated results.

## 1. Confirm the mechanism

Instrument (temporarily, or via a debug flag) `simulate_dual_porosity` to record,
per timestep, the actual `growth` value distribution across macropore-region
cells — not just the boolean "active" count. Plot this for Vertisol's Stage 11 (or
Stage 13, if real data changed its config) run:
- If growth values cross the 1e-9 threshold near-simultaneously across most
  regions within the first few timesteps, that confirms the equilibration
  hypothesis directly (not just by elimination).
- If instead some regions cross late and the metric only LOOKS flat because a few
  already-large regions dominate a simple count, that's a different, more
  fixable, explanation — report whichever is actually true.

## 2. Assess whether a better recruitment metric exists

Without changing the underlying transport or growth physics, test whether a
graded metric (e.g. total growth-weighted "active mass," or a higher activation
threshold, or time-to-half-max-growth per region) tracks meaningful new
information the boolean count misses. Compare a candidate metric against `R(t)`'s
own known-genuine gradual rise (Stage 11 established this is real, driven by
per-region growth RATE, not new regions joining) — a good metric should visibly
track that rate-building-up story instead of saturating instantly.

## 3. Do not touch calibrated results

This stage must not change `om_total`, coating fraction, `D_scale`, `T`, or any
other Stage 10/11/13 frozen parameter. If a new diagnostic metric is added, it is
purely observational — confirm via an exact numeric identity check (as Stage 11
did for Sand) that `R(t)`, `B(t)`, and all existing headline numbers are byte- or
noise-level identical before/after this stage's changes.

## 4. Deliverables

- `stage14_run.py`: per-timestep growth-distribution instrumentation, the
  mechanism-confirmation figure, and (if useful) the candidate graded metric
  computed alongside the existing boolean one, both plotted against `R(t)`.
- `results/stage14/`: the growth-distribution-over-time figure; old-vs-new
  recruitment metric comparison; the identity-check confirming nothing else moved.
- `docs/stage_results/STAGE14_RESULTS.md`: state plainly whether the
  near-instant-equilibration hypothesis was confirmed, and whether a better
  metric is recommended for future stages (as an addition, not a replacement,
  unless the old one is shown to be actively misleading).
- Update `LOGBOOK.md`, `report_prompt.md`. Only touch `MODEL_SPEC.md` if a new
  diagnostic function is added to the documented API.

## 5. Order of work

1. Instrument and plot the per-timestep growth distribution (Vertisol).
2. Confirm or refute the equilibration hypothesis directly from that plot.
3. If refuted, identify the actual mechanism instead — do not force-fit the
   original hypothesis.
4. Trial a graded recruitment metric; compare against R(t)'s known rate-building
   story.
5. Verify byte/noise-level identity of all existing calibrated results.
6. STOP; write `STAGE14_RESULTS.md`.
