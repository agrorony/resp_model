# Soil-structure microbiome respiration model

Simulates a single microbial species growing logistically on a 2D (or 3D)
pore-structured grid. Three soils (A, B, C) share identical biology AND
identical amounts (organic matter, initial biomass, total habitat sum(K) per
unit volume); they differ only in their pore field (pore-size texture and
spatial arrangement). Their respiration curves are nevertheless clearly
distinct. See `V3_RESULTS.md` for the short result, `MODEL_SPEC.md` for the
equations and rules (SS13-15 for v3), and `LOGBOOK.md` for the full record.

## Setup

```
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Running

```
python run.py --results-dir results/v3   # all three soils: CSVs, figures, amounts.csv
python run.py --soil A                   # a single soil
python search_loop.py "label"            # S1-S5 check of the current configs (+ nudges if failing)
```

`scan.py` holds the library scan / triple selection / seed-robustness and
shuffle-control checks used by the search.

## Files

- `model.py` — the simulation engine (§2 of the spec). No soil-specific
  values are hard-coded; everything comes from `configs/`.
- `configs/biology.yaml` — frozen constants shared by all soils.
- `configs/soil_A.yaml`, `soil_B.yaml`, `soil_C.yaml` — structure-only
  configs (diffusion, carrying capacity, OM placement, horizon).
- `run.py` — runs one or all soils; writes `results/*.csv` and the required
  figures.
- `metrics.py` — pairwise distance, shape describer, emergent-distinctness check.
- `pore_field.py` — n-D pore-field generator, fabric archetypes, K/D/OM/B0 rules.
- `configs/mapping.yaml` — shared pore->parameter constants and per-volume amounts.
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

v3 succeeds under the equal-totals rules (see `V3_RESULTS.md`). The v2
configs and figures are kept in `configs/v2/` and `results/` for comparison.
