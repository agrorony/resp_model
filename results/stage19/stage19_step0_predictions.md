# Stage 19 Step 0 -- regime predictions, written BEFORE any Step 1 run

Budget from the grids alone (`stage19_step0_budget.json`). Frozen biology:
f* = (m0+m_s)/(r_max+m_s) = 0.109, so S* = 0.0245. B0 starvation decay
timescale 1/(m0+m_s) = 8.3. Dissolution 1/k_dis = 20, leak 1/k_leak = 10,
max growth 1/r_max = 1, window 75 time units.

R5 shared densities: B0 0.0125, K_t 0.1096 (mean of native Sand 0.018 and
Vertisol 0.201 at n=128).

Habitat is overwhelmingly DRY in both soils (Sand 99.9-100%, Vertisol
89-91%), so delivery into habitat is mostly via the water film
(D_film x D_scale = 4.70) plus matrix->macro leak.

1. **Threshold (death vs life).** Time for habitat S to reach S* from
   coating release alone, at OM_density x and ignoring transport lag and
   uptake: Sand 0.0245*0.82/(0.05*0.8*x) = 0.50/x; Vertisol
   0.0245*0.32/(0.05*0.3*x) = 0.52/x. At x=0.037 that is ~14 time units,
   longer than the starvation timescale (8.3), plus ~10 of leak lag, so
   **x=0.037: starvation, at best partial recovery, both soils**. At
   x>=0.084 it is <=6 time units, so **alive**.
2. **Saturation (bloom).** Maximum biomass from all OM is Y*x*N = 0.4x per
   cell against the habitat ceiling 0.11 per cell, so OM can fill K once
   x >~ 0.27. **x = 0.44 and 1.0: K-saturated, Stage-7-like bloom**
   (shapes set by logistic saturation, expected to converge).
3. **Supply-limited window (the target): x ~ 0.084-0.19**, where growth is
   limited by delivery. Delivery timescales (10-20 for coating, 44-98 for
   trapped OM by matrix diffusion) exceed growth (~1) by >= 10x, so R(t)
   should follow delivery kinetics.
4. **Which soil differs how.** Sand's coating OM is concentrated (7.1-7.7x
   the mean OM density per coating cell, 80% of OM) next to one giant dry
   region; its trapped 20% sits in only 109 cells (30x density) ~4 cells
   from habitat. Vertisol spreads 30% coating thinly (1.15x) over 70
   regions and keeps 70% trapped (3.9x density) ~5 cells deep (22 trapped
   cells unreachable). Expected: Sand earlier and sharper, Vertisol slower
   and more sustained, if the supply-limited regime is reached.
