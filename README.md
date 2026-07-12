# Soil-structure microbiome respiration model

Simulates a single microbial species growing logistically on a 2D grid,
where three soils (A, B, C) share identical biology but differ only in
spatial structure (diffusion rate, carrying-capacity layout, organic-matter
placement). This produces three visually distinct respiration curves:
Soil A rises, Soil B falls, Soil C stays flat/bumpy. See `MODEL_SPEC.md` for
the full equations and reasoning.

## Setup

```
python -m venv .venv
.venv/Scripts/pip install numpy matplotlib pyyaml scipy   # scipy unused but optional per spec
```

(On Linux/Mac: `source .venv/bin/activate` instead of using `.venv/Scripts/...` directly.)

## Running

```
# run all three soils: writes CSVs + all figures to results/
.venv/Scripts/python.exe run.py

# run a single soil
.venv/Scripts/python.exe run.py --soil A

# re-run the budget-capped structure search (see LOGBOOK.md for its log)
.venv/Scripts/python.exe search_loop.py
```

## Files

- `model.py` — the simulation engine (§2 of the spec). No soil-specific
  values are hard-coded; everything comes from `configs/`.
- `configs/biology.yaml` — frozen constants shared by all soils.
- `configs/soil_A.yaml`, `soil_B.yaml`, `soil_C.yaml` — structure-only
  configs (diffusion, carrying capacity, OM placement, horizon).
- `run.py` — runs one or all soils; writes `results/*.csv` and the required
  figures.
- `metrics.py` — the rising/falling/flat classifier and success check.
- `search_loop.py` — the budget-capped (`MAX_ITERATIONS = 12`) structure
  search described in the original prompt; appends to `LOGBOOK.md` every
  iteration and only ever edits `configs/soil_*.yaml`, never
  `configs/biology.yaml`.
- `LOGBOOK.md` — the honest per-iteration record of what was tried and why.
- `MODEL_SPEC.md` — equations, frozen constants, structure knobs, and the
  reasoning behind each soil's configuration.
- `report_prompt.md` — a self-contained prompt (not run yet) that generates
  the final written report from `MODEL_SPEC.md`, `LOGBOOK.md`, and
  `results/`.
- `results/` — `respiration_curves.png`, `cumulative_co2.png`,
  `structure_maps_{A,B,C}.png`, and one CSV per soil with
  `time, R_t, cumulative_co2`.

## Current status

The structure search succeeded on the first iteration (see `LOGBOOK.md`):
all three soils meet their target verdict (rising / falling / flat) under
the automated `metrics.py` check, and the plotted curves in
`results/respiration_curves.png` are visually distinct.
