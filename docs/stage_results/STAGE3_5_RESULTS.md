# Stage 3.5 results — literature-informed PSDs, constant moisture

*Runner: `stage35_run.py`. Outputs: `results/stage35/`. Soils relabelled
**Loess / Sand / Vertisol** (Loess ↔ the old intermediate/bimodal soil,
Sand ↔ coarse, Vertisol ↔ the high-porosity multimodal soil).*

## 1. The idea

Stage 3 fed the model each soil's real, image-derived PSD and could only tell
the fine soil apart from the two coarse ones. Those PSDs were poor: 30 µm
minimum resolvable diameter, imperfect Otsu segmentation, and (for the
highest-porosity Vertisol) whole pore modes missing. Stage 3.5 keeps the entire
Stage-3 pipeline and swaps **only** the source of the PSD: each soil's PSD is
now a **mixture of lognormal modes over pore diameter (µm)** with mode
diameters, geometric SDs and volume weights taken from the soil-physics
literature for the three real soil types (`configs/psd_literature.yaml`,
`psd_parametric.py`, `pore.mode: psd_parametric`).

The validation target is my measured respiration pattern, **qualitatively
only** — ranking and shape, not absolute mg CO₂. My experiment was at
**constant moisture** (closed humid jars that did not dry out), so this stage
runs at a single fixed θ and reads the temporal shape of R(t), rather than
sweeping θ.

| Soil | measured P1→P4 (mg CO₂ kg⁻¹ d⁻¹) | shape |
|---|---|---|
| Vertisol | 104 → 142 → 158 → 173 | **rising** |
| Sand | 116 → 94 → 90 → 74 | **falling** |
| Loess | 68 → 25 → 50 → 11 | **low / erratic** |

## 2. What was built (matches the prompt's deliverables)

- **`psd_parametric.py`** — a lognormal-mixture PSD builder that returns the
  *same dict shape* `psd_data.load_soil_psd` returns, so everything downstream
  (the empirical-quantile pore-field mapping, the habitable-fraction and
  median-diameter stats) is untouched.
- **`configs/psd_literature.yaml`** — the per-soil mixture parameters, with the
  rationale/citations as comments. Editable in one place.
- **`configs/mapping.yaml`** — `D_d_ref` raised 90 → **250 µm** (the D(d)
  fusion fix, §4 below); added `matric_d_cut` for the shared-suction diagnostic.
- **`configs/soil_{loess,sand,vertisol}.yaml`** — renamed configs,
  `mode: psd_parametric`, constant-moisture saturation.
- **`stage35_run.py`** — the constant-θ run at n=128, the experiment overlay,
  the θ-sensitivity check, the shared-matric-potential diagnostic, and the
  isolation/distinctness diagnostics.

**Loaded PSDs (sanity check, `stage35_psd_summary.csv` / `stage35_psd_loaded.png`):**

| Soil | porosity | median µm | vol <30 µm | **habitable 30–150 µm** | vol >150 µm |
|---|---|---|---|---|---|
| Loess | 0.47 | 14.2 | 0.678 | 0.244 | 0.078 |
| Sand | 0.38 | 100.7 | 0.183 | **0.486** | 0.330 |
| Vertisol | 0.55 | 0.35 | 0.678 | 0.259 | 0.064 |

These are literature-informed **representative** distributions for the three
soil *types*, not measurements of my specific samples.

## 3. Two bugs the literature PSDs exposed (both fixed)

1. **`K == 0` division (NaN).** Stage-3's measured PSDs bottomed out at 30 µm,
   so the physiological `K(d)` (hard-zero below 30 µm) was never *exactly* zero
   and `1 − B/K` was always finite. The literature PSDs reach down to clay
   micropores, so most of the Vertisol's grid has `K == 0`; `1 − B/K` became
   `0/0` and NaNs poisoned every curve. Fixed in `model.simulate` by masking
   the logistic term where `K == 0` (those cells are non-habitat — they hold
   water and conduct, but nothing lives there), instead of dividing.

2. **D(d) fusion (§3 of the prompt).** Stage 3 fused the two coarse soils
   because `D(d)` saturated at `D_d_ref = 90 µm` — its knee sat *inside* the
   30–150 µm window, so every large pore looked alike to transport. Raising the
   reference to **250 µm** moves the knee above the real diameter range:

   | `D_d_ref` | D(30) | D(150) | D(300) | **D(300)/D(30)** | coarse-pair separation |
   |---|---|---|---|---|---|
   | 90 (Stage 3) | 0.206 | 1.184 | 1.783 | 8.6× | 1.24× |
   | **250 (Stage 3.5)** | 0.054 | 0.439 | 0.907 | **16.7×** | **1.62×** |

   D now keeps rising appreciably across 30–300 µm, and the two measured
   coarse soils that Stage 3 fused (nlm, mishmar) separate by 1.62× in mean D
   instead of 1.24×. The stability bound still holds: `max(D)·dt/dx² = 0.126 <
   0.2`. (`stage35_D_of_d_fix.png` / `.csv`.)

## 4. Headline result — constant θ = 0.6, n = 128

**The model reproduces the measured shape for Sand only, and not robustly. It
does NOT reproduce Vertisol-rising or Loess-erratic.**
(`stage35_experiment_overlay.png`, `stage35_constant_theta.csv`.)

| Soil | habitat frac (K>0) | **wet-habitat frac** | outcome | model shape | measured shape | match |
|---|---|---|---|---|---|---|
| Loess | 0.322 | **0.000** | declining | falling | erratic | ✗ |
| Sand | 0.817 | **0.417** | declining | falling | falling | ✓ |
| Vertisol | 0.323 | **0.000** | declining | falling | rising | ✗ |

All three curves *decline*. Loess and Vertisol are **bit-for-bit identical**
(both `B_final/B0 = 0.0270`, exactly the pure-starvation floor
`exp(−(m0+m_s)·T)`), so their pairwise distinctness is ~10⁻⁶ — the model gives
them literally the same curve.

### Why — the mechanism (this is the real finding)

Small pores fill first (capillarity), so a shared θ wets a *different diameter
range* in each soil. At θ = 0.6 the wetting cutoff diameter is:

- **Sand → 127 µm** — its 30–150 µm habitat is *below* the cutoff, so the
  habitat is wet (wet-habitat 0.42). Sand is the only soil that can respire.
- **Loess → 20.6 µm** and **Vertisol → 0.75 µm** — their habitat macropores
  are *above* the cutoff, so the habitat is **bone dry**. A dry cell has D = 0
  on every face, so the biomass that seeds there (b0 ∝ K) is cut off from all
  substrate and simply starves. Two soils with dry habitat give the identical
  starvation curve, regardless of how different their PSDs are.

There is a **second, deeper bottleneck** that limits even Sand. The OM-placement
rule (`OM ∝ d^-2.5`, frozen in `mapping.yaml`) was calibrated for the Stage-3
diameter range (30–2000 µm). Over the literature PSDs' 7-order-of-magnitude
range (down to ~10⁻⁴ µm clay), `d^-2.5` puts **~99.96 % of the OM on
sub-micron, non-habitat cells**. The substrate has to diffuse from there to the
habitat, but the diffusion reach over the whole run is only ~1–3 cells
(`√(D·T)` for D ≈ 0.05–0.33). So even in Sand the mean substrate concentration
that reaches the habitat is ~0.002 (vs `Ks = 0.2`, i.e. `f(S) ≈ 0.01`): Sand's
"falling" is the depletion of a tiny initial substrate pool, not a genuine
early burst. **Sand matches by having the least-starved habitat, not by
reproducing the mechanism.**

## 5. Robustness — the one match is grid-fragile

At **n = 40** the same constant-θ run flips Sand to *thriving / rising*
(`B_final/B0 = 0.23`), because a coarser grid places the wet fine-pore OM and
the wet habitat cells closer together, so substrate bridges the gap. At n = 128
Sand declines. So the single reproduced shape (Sand-falling) is a
grid-resolution artifact of that ~1–3-cell diffusion bottleneck, not a stable
model prediction. (`stage35_constant_theta_n40.csv`,
`stage35_respiration_constant_theta_n40.png`.)

The θ-sensitivity sweep (`stage35_theta_sensitivity.png/.csv`) confirms the
ranking is not knife-edge on θ = 0.6 — it is *robustly wrong* in the same way
across θ ∈ {0.4…0.9}: Loess and Vertisol only get *any* wet habitat above
θ ≈ 0.7, and even then it is a thin coarse-pore rim, never enough to invert the
ranking.

## 6. Diagnostic — a shared *matric potential* instead of a shared θ

A closed humid jar imposes the same **matric potential** on all three soils,
and three soils at one potential sit at three *different* θ (that is what a
water-retention curve is). The `matric_d` saturation mode holds the air-entry
diameter fixed (`d_cut = 60 µm`) and lets θ emerge per soil
(`stage35_matric_diagnostic.csv`):

| Soil | emergent θ | wet-habitat frac | outcome | model shape |
|---|---|---|---|---|
| Loess | 0.786 | 0.108 | declining | falling |
| Sand | 0.311 | 0.128 | declining | falling |
| Vertisol | 0.769 | 0.091 | declining | falling |

This **corrects the water-retention artifact** — now all three soils have a
live (non-zero) wet habitat — yet the model *still* reproduces only Sand and
*still* declines everywhere, because the OM/diffusion bottleneck (§4) is
independent of it. So the failure to get Vertisol-rising is not merely the
shared-θ assumption; the substrate cannot reach the habitat in any of these
configurations.

## 7. Honest verdict

With literature-informed PSDs at constant moisture, the model **does not**
reproduce the measured ranking and shapes:

- **Sand — falling: reproduced**, but for the wrong reason (least-starved
  habitat) and not robustly (flips to rising at n = 40).
- **Vertisol — rising: not reproduced.** The model gives it a dead-flat
  starvation decline. The intended story (OM protected in wet fine pores,
  slowly released to the sparse habitat macropores, access building over time →
  rising) requires substrate to actually traverse from the fine matrix to the
  habitat; under the frozen OM rule and D(d), it never arrives.
- **Loess — erratic: not reproduced.** Monotonic starvation decline; the
  bimodal structure produces no non-monotonicity here.

## 8. Main remaining limitation (candidate next stage)

The blocker is **not** the PSDs — they are correct and the D(d) fusion is
fixed. It is that **the substrate never reaches the habitat**, for two coupled
reasons:

1. **OM hyper-concentration.** The frozen `OM ∝ d^-2.5` rule, applied across
   the literature PSDs' enormous dynamic range, dumps essentially all OM onto a
   handful of sub-micron clay cells, from which diffusion cannot deliver it to
   the 30–150 µm habitat within the run. The intended "OM protected in fine
   pores, slowly released" mechanism needs an OM rule (or a dissolution/
   diffusion balance) that keeps substrate *mobile toward* the habitat, not
   trapped 7 orders of magnitude away in diameter. This is the first thing to
   revisit.
2. **Every grid cell is a pore.** The grid does not represent bulk porosity, so
   a soil that is 68 % sub-micron clay *by volume* (Vertisol) devotes 68 % of
   its **cells** to isolated dry OM traps. Total porosity differs across these
   soils (Vertisol highest) but the model cannot express it. Representing bulk
   porosity — so fine-pore volume is packing/retention, not a wall of dead
   cells — is the natural next stage.
