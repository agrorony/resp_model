# Stage 22 — prediction (committed after the wetting map and OFF-switch check, before any equal-wetting respiration run)

Wetting map (results/stage22/step0.json): equal-50% wetting gives Sand wet habitat 50.0% (was 1.2%), f_rel 0.099 (was 0.0060, x17);
Vertisol wet habitat 44.3% (was 73.8%), f_rel 0.099 (was 0.817, x0.12). OFF switch reproduces Stage 21 S and V bit-for-bit (rel. err 0).
Native K density shifts: Sand 0.137 -> 0.222, Vertisol 0.192 -> 0.156 (K_density stays 0.16445; the film-K 0.45 factor now applies to a different habitat share).

Predictions (written before the runs):
- Sand: much larger diffusive access (film factor x17, half of habitat wet) -> more R, faster peak, less dependence on exchange.
- Vertisol: film transport in half of its pores is now ~8x weaker -> lower/later R peak, more exchange dependence.
- The two soils move toward each other. I expect the distance to fall well below Stage 21's 1.229 but stay >= 0.3 (matrix fraction 0.08 vs 0.75, pore size scale and OM placement still differ),
  with both PRIMARY-alive. A collapse (< 0.3) is possible and would be reported as such.
