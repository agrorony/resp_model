# Claude Code task: Stage 13 — real measured PSD + connectivity data (replacing literature parameters)

Stages 3.5–12 used literature-parametrized pore-size distributions
(`configs/psd_literature.yaml`, via `psd_parametric.py`) for Loess/Sand/Vertisol,
and a purely geometric (inscribed-circle) macropore-sizing rule with no direct tie
to any measured connectivity number. This stage replaces both with real
measurements from the companion project, so the model stops being "informed by
soil type" and starts being "informed by this specific project's own data."

This is a DATA-SOURCING + RECALIBRATION-CHECK stage, not a new mechanism. Do not
change `dual_porosity.py`'s physics.

## 0. Find the current canonical files (do not assume stale paths)

The companion project lives in a separate connected folder ("resarch exercise").
Its data layout has changed since `data/psd/{mishmar,rehovot,nlm}/` was copied into
this repo back at Stage 3 — do not reuse those old copies uncritically. Instead:

1. Read `PROJECT_STATUS.md` and `DATA_CATALOG.md` in the "resarch exercise" folder
   first, in full, to find the CURRENT canonical location of:
   - Each soil's real measured pore-size distribution (PSD) — note Mishmar
     (Loess) has multiple resolutions (native 5.85 µm vs matched ~15 µm); prefer
     whichever resolution the project itself uses for cross-soil comparison
     (see "A2" in `PROJECT_STATUS.md`).
   - Track E connectivity metrics per soil: Euler characteristic (χ), connectivity
     density, Γ (percolation/loop-vs-fragment index), χ(r) crossover radius r*,
     tortuosity where available. These live under `Topology_Metrics_Aug2026V2`
     and/or `Statistical_Analysis_Aug2026` per the catalog — confirm by reading it,
     don't guess the path.
   - Bulk density / porosity per soil (`Table 1`/`Table 2` equivalents), needed to
     convert relative pore fractions into this model's absolute units correctly.
2. Log exactly which files/versions you used (path + date) in this stage's
   `LOGBOOK.md` entry — this matters because that project's data has been revised
   multiple times (see its "Document state" history) and a future reader needs to
   know which snapshot fed this model.
3. If a needed number genuinely isn't available (e.g. a soil's tortuosity is NaN
   in the source due to the documented `porespy` solver bug), say so plainly and
   fall back to the nearest documented substitute — do not silently interpolate
   or invent a value.

## 1. Replace the PSD source

Write a new loader (or extend `psd_data.py`) that builds each soil's PSD from the
real measured table located in step 0, in the same shape `psd_parametric.py`
currently returns, so the rest of the pipeline (`build_dual_grids` etc.) is
unchanged. Keep `psd_literature.yaml`'s mixture-fit path available behind a config
flag (`pore.mode: psd_measured` vs `psd_parametric`) rather than deleting it — this
preserves Stages 3.5–12's reproducibility.

## 2. Tie macropore/connectivity structure to the real numbers

Stage 5's inscribed-circle macropore sizing and region-connectivity behavior
should now be checked against the real Γ / χ(r) crossover radius per soil, not
just used as an internal geometric measure:
- Report, per soil, how the model's own emergent region count / largest-region
  fraction compares to the real Track E numbers (e.g. Vertisol's real network is
  reported as fragment-dominated with no crossover in range — does the model's
  Vertisol also come out as many small disconnected regions, or does it disagree?).
- Do NOT force a match by tuning geometry to hit the real numbers — report
  agreement or disagreement honestly, the same way Stage 3.5/4/5 reported
  mechanism mismatches.

## 3. Recalibration check (frozen physics, new inputs)

Re-run the Stage 10/11 per-soil calibrated configuration (Sand, Vertisol) with the
real PSD substituted for the literature one, parameters otherwise unchanged first.
Then, ONLY if the real data breaks Stage 10/11's §3 pass/fail table (burst-crash
for Sand, sustained rise for Vertisol), re-run a budget-capped search
(`MAX_ITERATIONS=15`, same style as Stage 7) over the same per-soil knobs Stage 10
already established (`om_total`, coating fraction, `D_scale`, `T`) — do not
introduce new free parameters. Report whether recalibration was needed at all;
if the real data reproduces Stage 11's numbers closely with no retuning, that
itself is a result worth stating plainly.

## 4. Deliverables

- `psd_data.py` (or a new `psd_measured.py`) + a `pore.mode: psd_measured` config
  path; new `configs/soil_{loess,sand,vertisol}_measured.yaml`.
- `stage13_run.py`: loads real PSD + connectivity data, runs the recalibration
  check for Sand + Vertisol, writes comparison figures (model vs Stage-11
  literature-PSD baseline, model region-count vs real Γ/crossover-radius).
- `results/stage13/`: headline CSV, comparison figures, a data-provenance table
  (which file/version/date each number came from).
- `docs/stage_results/STAGE13_RESULTS.md`: honest verdict — did real data change
  the calibrated behavior, and does the model's emergent connectivity agree or
  disagree with Track E's measured connectivity story.
- Update `MODEL_SPEC.md`, `LOGBOOK.md`, `report_prompt.md` (append a Stage 13
  section to the existing prompt, same style as Stages 1–12).

## 5. Order of work

1. Read `PROJECT_STATUS.md` + `DATA_CATALOG.md` in "resarch exercise"; log sources.
2. Build the real-PSD loader; verify it loads without changing `dual_porosity.py`.
3. Run Stage 10/11's frozen config on real data; check against §3 pass/fail table.
4. If needed, budget-capped recalibration (≤15 iterations) on existing knobs only.
5. Region-connectivity vs Track E comparison (report, don't force agreement).
6. STOP; write `STAGE13_RESULTS.md` with the numbers and the honest verdict.
