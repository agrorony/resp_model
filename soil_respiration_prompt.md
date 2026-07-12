# Claude Code task: soil-structure microbiome respiration model

You are building a small, well-documented simulation project. A single microbial
species grows logistically on a 2D grid. **The biology is identical everywhere and
must NEVER be changed.** The three soils differ ONLY in physical structure. The goal
is that the three soils produce visually distinct soil-respiration curves:

- **Soil A** — respiration rate **rises** over the run
- **Soil B** — respiration rate **falls** over the run
- **Soil C** — respiration rate stays **flat / irregular (chaotic-looking)**

There is no calibration to real data. Shapes within the observation window are what matter.

---

## 0. Budget — READ FIRST (this is the stop condition)

Do NOT loop until "perfect." Work within a fixed budget and stop cleanly when it runs out.

```
MAX_ITERATIONS = 12      # one iteration = evaluate all 3 soils + one structure nudge
```

- This is my session-budget knob. Treat it as a hard ceiling.
- When you reach MAX_ITERATIONS **or** the success check passes (whichever comes first),
  STOP tuning and go straight to the write-up phase using the best configuration so far.
- If the three behaviors were not all achieved when the budget runs out, that is fine —
  document honestly what was and wasn't achieved. Do not silently keep going.

---

## 1. Hard constraints

1. **Biology is frozen.** The constants in §3 are shared by all three soils and must not
   be edited by the tuning loop, ever. If you are tempted to change a biology constant to
   make a soil behave, STOP — that is forbidden. Only structure (§4) may change.
2. **Totals are held equal across soils** so differences come purely from spatial structure:
   total initial OM, total initial biomass, and grid size are identical for A, B, C.
   The loop may only vary the *spatial arrangement* of structure and the diffusion rate D.
3. Deterministic only (no randomness in the dynamics; a fixed seed may be used to lay out
   patches). Grid ≤ 50×50. Each full run should take seconds.
4. Python with numpy + matplotlib (+ optionally scipy, pyyaml). Keep it simple and readable.

---

## 2. Model (implement exactly this; anchored on logistic growth)

Per grid cell (i, j), three state fields evolve: biomass `B`, mobile dissolved substrate
`S`, and a solid organic-matter pool `OM`.

```
dB/dt = r_max * f(S) * B * (1 - B / K_ij) - m0 * B - m_s * (1 - f(S)) * B
dS/dt = D_ij * laplacian(S) + k_dis * OM_ij - uptake
dOM/dt = -k_dis * OM_ij
```

with

```
f(S)   = S / (S + Ks)                         # saturating substrate response (Monod-type)
uptake = (1/Y) * r_max * f(S) * B * (1 - B/K_ij)   # substrate consumed to build biomass
```

- The logistic term `(1 - B/K_ij)` is the hard, fixed physical space limit per cell.
- `m0 * B` is constant maintenance death; `m_s * (1 - f(S)) * B` is starvation death that
  grows as local substrate runs out.
- `laplacian(S)` is the 5-point discrete Laplacian; `D_ij` is the local diffusion/flow rate.
  Respect numerical stability: `max(D) * dt / dx^2 < 0.2`. Use zero-flux (Neumann) boundaries.
- Integrate with explicit Euler (or RK4 if needed for stability). Clip states at ≥ 0.

**Respiration observable (the CO₂ curve):**

```
R_ij(t) = (1 - Y) * [substrate used for growth in cell ij] + m0 * B_ij    # growth + maintenance respiration
R(t)    = sum over all cells of R_ij(t)
```

Save `R(t)` (respiration rate) and its running integral (cumulative CO₂) for each soil.

---

## 3. Frozen biology constants (DO NOT TUNE)

Start with these (dimensionless units). If a soil misbehaves, fix it with structure, not these.

```
r_max = 1.00     # max growth rate
Ks    = 0.20     # half-saturation substrate
Y     = 0.40     # yield (biomass per substrate)
m0    = 0.02     # constant maintenance death
m_s   = 0.10     # max starvation death
k_dis = 0.05     # OM dissolution rate (solid OM -> mobile S)
dt    = 0.05     # time step
T     = enough steps to see the full shape (tune the horizon, NOT the biology)
```

Also equal across soils: total initial OM, total initial biomass, grid size, dx.

---

## 4. Structure — the ONLY things the loop may vary

Give each soil three structure maps plus a scalar. These are the free knobs:

- **D (pore / flow):** a base diffusion rate per soil, and optionally a spatial D-field.
  Low D = slow substrate delivery; high D = fast.
- **K_ij (habitable sites):** carrying-capacity map. Can be uniform or patchy; patch scale
  and contrast are free. (Keep the *mean* K equal across soils; vary only the pattern.)
- **OM layout:** the fixed total OM is scattered into patches. Free knobs = number/size of
  patches and their **distance from the high-biomass / high-K zones**.

**Intended mechanism (use this to steer your nudges):**

- Soil B (**falling**): high D, OM patches adjacent to bacteria → fast early burst, then the
  finite pool depletes → respiration declines.
- Soil A (**rising**): low D, OM patches far from bacteria → substrate only trickles in;
  access keeps growing through the window → respiration climbs.
- Soil C (**flat / chaotic**): patchy K + scattered OM at intermediate D → local patches peak
  at different times; early ones fall while late ones rise → the grid total is a bumpy plateau.

The controlling idea is the ratio of substrate-supply speed to consumption speed, set entirely
by structure. Nudge D and patch distance first; K-patchiness mainly for soil C's roughness.

---

## 5. Success check (automated proxy for "visually distinct")

Sample `R(t)` over the window and split into thirds. Let `e = mean(first third)`,
`l = mean(last third)`, `peak = max(R)`.

```
rising  : (l - e) >  0.15 * peak
falling : (e - l) >  0.15 * peak
flat    : |l - e| <  0.10 * peak  AND  roughness > threshold   # roughness = std of detrended R, or count of local direction changes
```

Success = Soil A passes `rising`, Soil B passes `falling`, Soil C passes `flat`, and the three
curves are clearly separable when plotted together. This proxy only decides when to stop the
loop; the final judgement is the plotted figure.

---

## 6. The loop

```
for iteration in 1..MAX_ITERATIONS:
    run soils A, B, C   (biology frozen; only structure differs)
    compute R(t) and the §5 metrics for each
    append a LOGBOOK entry (below)
    if success (§5): break
    else: nudge STRUCTURE ONLY toward the §4 mechanism (never biology), and continue
stop (budget or success) -> go to §7 write-up with the best config so far
```

Each LOGBOOK entry records: iteration number, the structure parameters used per soil
(D, K-pattern, OM layout/distance), the resulting metrics (e, l, peak, roughness, verdict per
soil), and a one-line note on what you changed and why. This log is what lets a later prompt
reconstruct the whole process — keep it complete and honest.

---

## 7. Deliverables (create all of these in the project folder)

```
model.py            # the engine in §2, no soil-specific values hard-coded
configs/
  biology.yaml      # the frozen constants (§3) — shared
  soil_A.yaml       # structure only
  soil_B.yaml
  soil_C.yaml
run.py              # runs one or all soils, saves R(t), cumulative CO2, and structure maps
metrics.py          # the §5 classifier + success check
search_loop.py      # the §6 budget-capped loop; updates LOGBOOK.md
results/
  respiration_curves.png     # R(t) for A, B, C on one axis
  cumulative_co2.png         # cumulative CO2 for A, B, C
  structure_maps_A/B/C.png   # D, K, OM maps per soil
  *.csv                      # raw R(t), cumulative, per soil
LOGBOOK.md          # appended every iteration (§6)
MODEL_SPEC.md       # the equations (§2), frozen constants (§3), structure knobs (§4), assumptions
README.md           # how to run everything
report_prompt.md    # see §8
```

Figures required: **respiration curves, cumulative CO₂, and the structure maps.** (No animation.)

---

## 8. Final report prompt (create this file, do not run it)

Create `report_prompt.md` containing a self-contained instruction that, when run later, will:
read `MODEL_SPEC.md`, `LOGBOOK.md`, and everything in `results/`, and write a structured
scientific report in **English** — motivation, model and equations, the three soils and their
structure, results (embedding the figures), interpretation of why each shape arises from
structure alone, limitations, and an honest note on what the budget did/didn't achieve. The
report should be written so it can be lifted into a course appendix.

---

## 9. Order of work

1. Build `model.py` and get ONE soil to produce a sane respiration curve (rise, peak, decay).
2. Add the three structure configs; check each target behavior by hand once.
3. Wrap in `search_loop.py` with metrics + LOGBOOK, capped by MAX_ITERATIONS.
4. Generate the final figures and CSVs from the best config.
5. Write `MODEL_SPEC.md`, `README.md`, and `report_prompt.md`.

Remember: biology is frozen, structure does all the work, and the budget in §0 is the stop.
