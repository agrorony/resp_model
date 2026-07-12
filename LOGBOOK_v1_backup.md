# LOGBOOK

Budget-capped structure search (soil_respiration_prompt.md §6).
MAX_ITERATIONS = 12. Biology is frozen throughout; only configs/soil_{A,B,C}.yaml are ever modified by nudges.

## Iteration 1 -- 2026-07-11T13:29:49
**Soil A** (target: rising)
- structure: D=0.05, K_pattern=single_cluster, OM_placement=far, n_patches=3, T=600
- metrics: e=1.555, l=6.395, peak=7.311, roughness=1.382, direction_changes=3, verdict=rising
- note: target met, no change needed
**Soil B** (target: falling)
- structure: D=3.0, K_pattern=single_cluster, OM_placement=near, n_patches=3, T=600
- metrics: e=7.568, l=4.040, peak=8.420, roughness=0.937, direction_changes=1, verdict=falling
- note: target met, no change needed
**Soil C** (target: flat)
- structure: D=0.3, K_pattern=multi_cluster, OM_placement=explicit, n_patches=5, T=600
- metrics: e=4.143, l=4.403, peak=6.668, roughness=1.354, direction_changes=1, verdict=flat
- note: target met, no change needed

**Iteration result: SUCCESS**
