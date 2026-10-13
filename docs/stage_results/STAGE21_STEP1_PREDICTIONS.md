# Stage 21 -- Step 1 predictions (committed BEFORE the 256^3 geometry run and before any respiration run)

Honesty note: the structure parameters were fixed first (`packing3d.STRUCTURE_PARAMS`).
A 10 um / 128^3 *geometry-only* check of the five builders was run to debug the code
(no respiration); its headline numbers inform the predictions below. No parameter was
adjusted after seeing it.

Expected geometry
- **S** (spheres + quartz + fines): resolved porosity ~0.32, matrix 0.082, solid 0.60.
  Pore network spans all axes, one dominant cluster (>99%). Pores large (LT median ~100 um),
  so S_w is low (~0.3) and most habitat (~90%) is drained (>= 50 um).
- **V** (plates + porous aggregates): resolved porosity 0.25, matrix 0.75, no solid.
  Slit-like pores, LT median ~30 um, S_w ~0.95; most habitat wet.
- **Sp** (plate-shaped solid): like S but larger LT (plates leave larger voids), even more
  drained. **Vp** (sphere aggregates): like V with somewhat wider pores. **Val**
  (aligned plates): like V in porosity; LT slightly narrower and pore network anisotropic.

Expected regimes (what the physics should do, to be compared with Step 2/3)
- **S**: most habitat is drained, so diffusion through habitat is film-limited
  (D x S_w^(10/3) ~ 0.015 x the wet value). Delivery of OM from the fines is by exchange
  into the adjacent pore voxels; local, quickly depleted. Expect an early, sharp R peak.
- **V**: habitat mostly wet, OM everywhere in a 75% matrix, within 20 um of pore for ~50%.
  Expect a broader, later and lower peak, more diffusion-coupled; transport-off should hurt
  V less than S (V's OM sits next to habitat), or the reverse if matrix->pore delivery is
  exchange-limited.
- Factorial: expect the **material** (OM location relative to habitat, S_w) effect to dominate
  the **shape** effect. Orientation (V_al vs V): small effect (distance << 0.3).
- Possible failure: V's habitat volume fraction is near its porosity (0.25) so equal-K scaling
  gives V a low per-habitat-voxel K; Sand's K is spread over a bigger habitat.
