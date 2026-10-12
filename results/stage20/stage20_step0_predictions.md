# Stage 20 Step 0 -- predictions, written BEFORE any Step 1 run

Budget: `stage20_step0_budget.json` (grids only). Kept from Stage 19: K_density
0.1096, B0 0.0125, coating fractions 0.8 / 0.3, D_scale 12.2582, T 1500.

Key numbers (n=128; n=40 similar):
- Sand W0: S_w 0.128, f_rel 0.0011, 0.06% of habitat wet. Drained D: F0
  4.70 -> F1 0.013 (360x lower).
- Sand W1: S_w 0.50, f_rel 0.099, 42.7% of habitat wet (water-filled D
  ~21.7 with cluster bonus). Drained D: F1 1.23.
- Vertisol W0: S_w 0.61, f_rel 0.197, 8.7% of habitat wet. Drained/full-pore
  D ratio: F0 4.9 -> F1 0.67 (the inversion is removed).
- Vertisol W1: S_w 0.50, but **0% of habitat wet** -- the fine matrix takes
  all the water (87% of matrix cells wet). Drained D under F1 0.32. Coating
  cells fall 4,269 -> 2,688 (coating needs a wet matrix path), trapped rise
  2,946 -> 4,527. Reported, not fixed.
- Native K_t density under W1: Sand 0.028 (W0 0.018), Vertisol 0.187 (W0
  0.201); 0.1096 kept as specified.

Predictions (Stage 19 / A0: both alive only at OM >= 0.439, distance 0.223
there):

- **A1 (F1+W0).** Sand's in-habitat transport effectively disappears (film
  0.013), so substrate enters habitat only through the matrix->macro leak at
  matrix interfaces. Interfacial cells can grow to K while the interior of
  the one giant region starves. Sand liveness becomes marginal: expected
  only at the highest densities, possibly not at all by the B_peak/B0
  criterion. Vertisol's film drops 4.70 -> 0.63 (7x slower), so it needs at
  least as much OM as in A0 (alive at 1.0, maybe not at 0.439). Expected:
  few or no both-alive densities. Any distance is likely to come from Sand
  being transport-starved, not from a healthy contrast.
- **A2 (F0+W1, diagnostic).** Sand gains fast wet habitat (43% wet, D ~22)
  plus the 4.70 film, so it comes alive at lower density than in A0.
  Vertisol's habitat drains completely (film 4.70 still) and its coating
  shrinks, so it is similar to or worse than A0. Expected: Sand alive at
  more densities; both alive at 0.439 and 1.0; distance similar to or larger
  than A0.
- **A3 (F1+W1).** Sand: 43% wet fast habitat plus film 1.23, so alive from
  low-mid densities. Vertisol: fully drained habitat with film 0.32 (15x
  below A0) and less coating OM, so likely dead except at 1.0. Expected: the
  highest distances again where Vertisol is dead; both alive at most at 1.0,
  where curves tend to converge. Reaching both-alive AND distance > 0.3 is
  unlikely.
- Overall expectation: physically consistent films make the
  liveness-distinctness trade-off WORSE for Vertisol, because Stage 19's
  liveness leaned on the inflated film.
