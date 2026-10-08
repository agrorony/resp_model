# v3 Results -- distinct respiration from structure, with equal amounts

**Goal.** Get three soils with clearly different respiration patterns where
the soils do *not* differ in the amount of organic matter, initial biomass
or total habitat -- only in their pore structure.

**Result: success** (all of S1-S5, MODEL_SPEC.md SS13) at path P1 -- 2D,
fully wet, no change to biology or mechanics.

| | Soil A | Soil B | Soil C |
|---|---|---|---|
| texture (mu, sigma) | 0.6, 0.3 | 0.6, 0.7 | 2.3, 0.7 |
| arrangement | large domains (lambda 8) | short, elongated (lambda 1.5, aniso 5) | large domains (lambda 8) |
| total OM / B0 / habitat sum(K) | 800 / 20 / 480 | 800 / 20 / 480 | 800 / 20 / 480 |
| shape | peak mid-window | slow rise | early burst, decline |
| cumulative CO2 (structured / shuffled) | 108 / 176 | 47 / 158 | 174 / 181 |

- Pairwise distance (threshold 0.3): AB 0.71, AC 0.72, BC 1.53.
- Robust: three seed sets give min distance 0.71 / 0.65 / 0.61; seed noise 0.03.
- Arrangement matters: shuffling each soil's pores (same sizes, random
  positions) drops the distances to 0.18 / 0.15 / 0.50.
- The old v2 soil C is rejected under v3: 577x OM hoarding in one cell.
- Limit: with an *identical* pore-size histogram for all three soils the best
  separation is ~0.09 -- see MODEL_SPEC.md SS15 for why and the next step.

Figures: `results/v3/respiration_curves.png` (dashed = shuffled control),
`results/v3/cumulative_co2.png`, `results/v3/structure_maps_{A,B,C}.png`
(r, K, D, log OM), `results/v3/amounts.csv`. Full record: `LOGBOOK.md` (v3 section).
