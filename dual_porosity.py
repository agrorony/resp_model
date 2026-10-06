"""
Stage 4 Variant C -- dual-porosity (matrix + macropore) transport
(prompts/soil_respiration_prompt_stage4_exploratory.md SS3).

Stage 3.5's diagnosis: the frozen single-continuum `OM ~ d^-2.5` rule dumps
~99.96% of carbon onto sub-micron clay cells, from which diffusion (reach
~1-3 cells/run) can never deliver it to the 30-150 um habitat. Variant A
(model.build_grids' `psd_truncate_floor_um` + wet-diffusion scaling) keeps
distributed OM on the truncated grid and speeds wet-cell transport. Variant C instead
keeps the FULL literature PSD (down to sub-micron clay) but splits every cell
into one of two kinds -- the "sponge + straws" model:

- MACROPORE cells (d >= mapping.yaml K_d_low, i.e. the same 30 um
  physiological threshold used everywhere else): fast D(d), biomass lives and
  respires here exactly as in every other stage, substrate diffuses between
  them via the usual conservative-flux transport.
- MATRIX cells (d < 30 um, including the sub-10 um fraction Variant A drops):
  K=0 (no biomass -- K_window is already a hard zero there, so this falls out
  of the existing physiological rule, no new code needed), and their fast
  diffusion pathway is switched off (`D=0`) so they cannot move substrate
  through the regular conservative-flux term -- they "block fast
  inter-macropore flow". Their OM is simply the standard fine-pore-biased
  `om_field` evaluated on the FULL (untruncated) diameter field: because fine
  pores are exactly the matrix cells, this rule already concentrates most OM
  there without any new placement logic. A matrix cell's own dissolved
  substrate (from local `k_dis * OM`) is not stranded, though: it LEAKS into
  any adjacent macropore cell at a fixed rate `k_leak` (`leak_exchange`
  below), a mobile-immobile / dual-permeability exchange term -- this is the
  only pathway matrix-stored carbon has to reach the habitat.

"Represent the sub-resolution fine fraction as the matrix cells' internal
porosity, not as separate cells" (SS3): no separate scalar field is added for
this -- the matrix cell's local OM density (from the existing fine-pore-bias
rule, now evaluated over the FULL PSD instead of a truncated one) already
stands in for how much fine/internal porosity that cell holds. Bulk porosity
is captured only qualitatively this way (a soil whose PSD puts more volume
below 30 um has a larger matrix_fraction and more OM concentrated within it --
reported per soil, see `build_dual_grids`'s `matrix_fraction`), exactly as the
prompt asks (Vertisol = lots of high-internal-porosity matrix; Sand = little).

Everything else -- biology, K(d), D(d), the empirical-quantile PSD mapping,
water saturation -- is identical to Stage 3.5; transport differs from
`model.simulate` by the added matrix->macro exchange terms and dynamic
hydro-substrate wicking (macro cells adjacent to wet matrix cells can rewet and
conduct during the run). Kept as a separate module so Stage 1-3.5 code paths
are untouched.

Stage 6 (prompts/soil_respiration_prompt_stage6_portback.md) ports two
conclusions proved in an abstract PDE sandbox (`pde_sandbox/`, its own
RESULTS.md) into this module:

1. **Implicit transport** (`implicit_diffusion_step`, ported from
   `pde_sandbox.sandbox.local_implicit_step`): the explicit conservative-flux
   update used through Stage 5 has a hard `max(D)*dt/dx^2 < 0.2` stability
   cap, which is exactly why dissolved substrate stayed frozen ~1-3 cells
   from where it dissolved (the recurring Stage 3.5/4/5 diagnosis). Solving
   `u - dt*div(D*grad(u)) = rhs` by Jacobi iteration instead removes that cap
   entirely -- unconditionally stable regardless of D, dt, or dx.
2. **Coating/trapped OM split** (`split_om_coating_trapped`, ported from
   `pde_sandbox.sandbox.build_down_source`'s geodesic construction): OM is no
   longer one pool per cell. COATING sits on porous-matrix cells within a
   short geodesic distance of a habitat cell along a wet connected path (fast
   release, a guaranteed finite habitat<-source path); TRAPPED sits in
   interior matrix cells far from habitat or behind barriers, releasing at a
   much slower rate `k_dis_slow` (a tail, not a permanent lock). The coating
   FRACTION is set per soil, the one deliberately hypothesis-motivated input.

Also removed (Stage 6 SS "Change 4"): the `S_cap` substrate-holding ceiling
(verified to be the compounding blocker keeping growth pinned regardless of
habitat wetness) and wicking is off by default (the sandbox used a static
structure; wicking previously overrode the Stage-5 wetting geometry).
`k_leak` and matrix internal diffusion stay as supporting pathways.

Stage 7 (prompts/soil_respiration_prompt_stage7_calibration.md) is a
quantitative CALIBRATION, not a new mechanism: Stage 6 proved the ported
transport + coating/trapped OM mechanisms work structurally but left
fed-habitat substrate ~25x too low for growth to beat maintenance decay
everywhere. Three changes, all shared/uniform across soils:

1. **Uniform OM percent** (`_variantc_om_fraction`): a single shared
   `stage7_om_fraction` (`configs/mapping.yaml`) replaces the old per-soil
   `stage4_variantC_om_fraction_<soil>` constants -- only each soil's
   porous_matrix_mask SIZE (structure) still makes `om_total` differ.
2. **Reach/D-scale calibration** (`stage7_D_scale`, `simulate_dual_porosity`):
   a shared post-scale multiplier on `D_macro_base`/matrix internal
   diffusion/`D_film`, raising diffusion's absolute magnitude (and hence
   `sqrt(D*T)` reach) without touching `D_d_ref` (Stage 3.5's soil-
   discriminating SHAPE of D(d), left alone).
3. **Water-film macropores** (`D_film`, `stage7_film_K_scale`,
   `build_dual_grids`): a DRAINED macropore cell (dry under the Stage 5
   matric rule) is no longer fully cut off (D=0) -- it keeps a small,
   shared film diffusivity and hosts activity at reduced capacity instead,
   giving soils like Sand (whose habitat is almost entirely one huge dry
   region) a fast-but-low-capacity supply pathway.

Stage 8 (prompts/soil_respiration_prompt_stage8_depletion.md) is a DEFINED
configuration with verification, not a search: Stage 7 closed the substrate
gap but, with abundant carbon, every soil bloomed into the same rising
curve because carbon never depleted. Stage 8's one change is Change 1 --
a single ABSOLUTE `stage8_om_total`, IDENTICAL for every soil (`build_
dual_grids`), replacing Stage 7's shared *fraction* (which still gave
unequal totals because `porous_matrix_mask` size is structural). The
coating/trapped split still follows each soil's own geometry -- only the
total mass is equalised. This is meant to let the model's EXISTING finite-
fuel + Monod-limitation physics (already present since Stage 6) actually
act: the same total carbon concentrated on Sand's small coating vs. spread
across Vertisol's many habitats should deplete on different timescales.
`simulate_dual_porosity` additionally tracks `fuel_remaining_t` (fraction
of a soil's initial total OM still undissolved) as direct depletion
evidence, alongside the existing `S_habitat_mean_t`/`habitats_recruited_t`.

Stage 9 (prompts/soil_respiration_prompt_stage9_ranked_om.md) is again a
DEFINED configuration with verification: Stage 8 proved the depletion
mechanism gives genuinely distinct curves and reproduced Loess-erratic, but
a single EQUAL total OM could not give Sand-fall AND Vertisol-rise
together (opposite regimes: Sand needs a lean pool to deplete, Vertisol
needs abundant, ACCESSIBLE carbon to sustain many regions). Three targeted
changes, ranking-only (never fitted to the measured respiration values):

1. **Ranked per-soil `om_total`** (`_variantc_om_total`, `stage9_om_base` x
   per-soil `stage9_om_multiplier_<soil>`, fixed ordering Vertisol > Loess >
   Sand): a real-world qualitative input (more carbon in a heavy, porous
   Vertisol than a coarse, low-porosity Sand), not a fitted amount.
2. **Vertisol's coating fraction raised** (`_variantc_coating_fraction`,
   `stage9_coating_fraction_vertisol`, up from Stage 6-8's 0.30): more of
   its now-larger carbon actually reaches its many habitat regions.
3. **Sand's coating fraction lowered** (`stage9_coating_fraction_sand`,
   down from 0.80): a smaller fast burst plus a slow trapped tail --
   subject to the SS3 "the burst wins" exception: never lowered past the
   point where Sand's rise-then-fall itself breaks.

Stage 10 (prompts/soil_respiration_prompt_stage10_final.md) is the FINAL
stage: two soils only (Vertisol, Sand -- Loess dropped). Stage 9 got the
qualitative direction right but its Sand never produced a genuine growth
burst (`R_peak` sat AT the trivial maintenance floor `m0*B0`, i.e. pure
decay of the seeded biomass, not growth-driven respiration) and its
Vertisol's "rise" was really an early spike that then dipped, because the
shared, fast `D_scale` (needed since Stage 7 just to get growth to fire at
all) recruits a soil's whole habitat almost at once rather than gradually.
Stage 10 explicitly sets aside the "shared/ranked scalars only" purity of
Stages 7-9 for this final two-soil deliverable and allows PER-SOIL tuning
of the physical parameters (never arbitrary curve-drawing):

1. **Per-soil `om_total`** (`_variantc_om_total`'s new highest-priority
   tier, `stage10_om_total_<soil>`) and **per-soil coating fraction**
   (`_variantc_coating_fraction`'s `stage10_coating_fraction_<soil>`):
   Sand gets a strongly total- and coating-concentrated dose so growth
   genuinely fires and then crashes as that finite pool empties; Vertisol
   gets a large, more matrix-distributed total so its many regions keep
   recruiting through the whole window.
2. **Per-soil transport magnitude** (`_variantc_D_scale`, `stage10_D_
   scale_<soil>`) -- the first time `D_scale` is NOT one shared constant:
   Vertisol's own reach is slowed relative to Sand's so its dispersed
   habitat regions ignite gradually over the window instead of together.

Stage 11 (prompts/soil_respiration_prompt_stage11_vertisol.md) keeps Sand's
Stage 10 configuration FROZEN and works only on Vertisol, plus a
performance fix:

1. **Performance**: `habitat_properties(macro_mask)` (a pure-Python
   connected-component BFS over the whole grid) and the `K_t`/`habitable`/
   `K_safe` fields derived from it were being recomputed EVERY timestep
   inside `simulate_dual_porosity`'s loop even though `macro_mask` is fixed
   for the whole run (wicking, even when enabled, only ever touches
   `wet_mask`) -- hoisted out of the loop, computed once. Combined with an
   optional Jacobi convergence-tolerance early stop
   (`implicit_diffusion_step`'s new `tol` argument, `stage11_implicit_tol`
   in `configs/mapping.yaml`, `stage6_implicit_iterations` kept as a safety
   CAP rather than lowered), this gave a ~7-60x speedup (soil- and
   region-count-dependent) with results identical to the pre-Stage-11
   version to within numerical noise (verified directly, not assumed --
   see `docs/stage_results/STAGE11_RESULTS.md`).
2. **Vertisol only** (`stage11_om_total_vertisol`, `stage11_coating_
   fraction_vertisol`, `stage11_D_scale_vertisol` -- new highest-priority
   tiers in `_variantc_om_total`/`_variantc_coating_fraction`/
   `_variantc_D_scale`; Sand has no Stage 11 key for any of the three, so
   it falls through to its Stage 10 values unchanged): total OM raised
   further and transport slowed further so biomass clearly grows above its
   initial total (not just clearing the trivial maintenance-floor
   comparison), while the rise stays genuinely gradual (peaks very late,
   no dip).

Stage 16 (prompts/soil_respiration_prompt_stage16_shared_om_rule.md) is a
MODEL-WIDE STRUCTURAL change, not a calibration stage: it retires the
per-soil coating/trapped OM split introduced in Stage 6 and restores ONE
shared, structure-blind organic-matter placement rule with spatially
UNIFORM dissolution, so that every difference -- between soils, and between
macropore regions within a soil -- must emerge from the pore field rather
than from a hand-set per-soil carbon parameter.

The change is motivated by a direct REFUTATION of Stage 14's causal
attribution (see the dated correction appended to
`docs/stage_results/STAGE14_RESULTS.md`). Stage 14 blamed Vertisol's
near-instant habitat activation on `implicit_diffusion_step`'s per-step
equilibration; re-running with transport switched off entirely
(`D_scale=1e-9`) leaves regions activating within 2-4 timesteps and
recruitment reaching its final count SOONER (0.27% of the window) than
with transport on (1.33%). The solver is exonerated. The real cause is
`split_om_coating_trapped` itself: defining the coating pool as every
porous-matrix cell within `coating_thickness` geodesic steps of ANY
macropore cell builds an identical fuel shell around all ~70 regions at
once, and with `B0 ~ K` seeding biomass everywhere at t=0 plus a single
global `k_dis`, every region is self-fuelled from t=0 and needs no
transport to activate. Simultaneity was built into the construction.

1. **Shared OM placement** (`place_om_shared_fine_pore`, selected by
   `om_mode: shared_fine_pore`, the new active default in
   `configs/mapping.yaml`): ALL organic matter is placed by the SAME rule
   for every soil -- the original fine-pore weighting
   `OM(i,j) ~ d(i,j)^-om_fine_bias` within the porous matrix
   (`MODEL_SPEC.md` SS4), normalized to that soil's total. No geodesic
   construction, no per-soil share.
2. **Uniform dissolution**: a single `k_dis` everywhere. Under this mode
   the slow pool is identically zero, so `k_dis_slow` has no effect and no
   coating/trapped pools exist. This enforces the stage's invariant: **no
   term anywhere in the model has a rate that depends on distance to
   habitat.** Decomposition is local (water contact, OM chemistry);
   distance enters ONLY through diffusion of the dissolved substrate,
   which the model already does.
3. **`coating_fraction` retired as an active parameter** -- it was a
   hand-set per-soil carbon knob, and the project's central claim (that
   structure alone drives the difference) is stronger without it.

The Stage 6-15 `split_om_coating_trapped` path is NOT deleted: it stays
reachable behind `om_mode: coating_trapped`, exactly as Stage 13 preserved
`psd_parametric` alongside `psd_measured`, and every Stage 6-15 runner
pins that mode explicitly on load so each earlier stage stays
independently reproducible (verified, not assumed -- see
`docs/stage_results/STAGE16_RESULTS.md`).
"""
from __future__ import annotations

import numpy as np

from model import load_yaml, harmonic_mean, water_mask_and_theta
from pore_field import (
    generate_diameter_field_from_psd,
    K_window,
    D_of_d,
    om_field,
    b0_field,
    gaussian_filter_2d,
)
import psd_parametric
import psd_measured


OM_MODE_SHARED = "shared_fine_pore"
OM_MODE_COATING = "coating_trapped"
OM_MODE_CLUMPS = "shared_lognormal_clumps"


def _om_mode(mapping_cfg):
    """Stage 16 (soil_respiration_prompt_stage16_shared_om_rule.md SS1):
    which organic-matter placement rule is active.

    - `shared_fine_pore` (the NEW ACTIVE DEFAULT, set in
      `configs/mapping.yaml`): ONE shared, structure-blind rule for every
      soil -- `place_om_shared_fine_pore` -- with spatially uniform
      dissolution at a single `k_dis`. No coating pool, no trapped pool,
      no per-soil `coating_fraction`, and no term anywhere whose rate
      depends on distance to habitat.
    - `coating_trapped`: the Stage 6-15 path (`split_om_coating_trapped`),
      PRESERVED, not deleted, exactly as Stage 13 preserved
      `psd_parametric` alongside `psd_measured`. Every Stage 6-15 runner
      pins this explicitly right after loading `mapping.yaml`, so each
      earlier stage stays independently reproducible even though the
      file's own default has moved on.

    Defaults to `coating_trapped` when the key is absent entirely, so an
    old config file that predates Stage 16 still reproduces its own
    behavior rather than silently switching rule."""
    mode = str(mapping_cfg.get("om_mode", OM_MODE_COATING)).strip().lower()
    valid = (OM_MODE_SHARED, OM_MODE_COATING, OM_MODE_CLUMPS)
    if mode not in valid:
        raise ValueError(f"unknown om_mode {mode!r}; expected one of {valid}")
    return mode


def _variantc_om_fraction(mapping_cfg, soil_name):
    """Stage 7 Change 1 (soil_respiration_prompt_stage7_calibration.md SS2):
    a single SHARED `stage7_om_fraction` for every soil, retiring the old
    per-soil `stage4_variantC_om_fraction_<soil>` constants -- the only thing
    that still differs between soils is the porous_matrix_mask SIZE
    (structure), never the percent. Falls back to a per-soil legacy key only
    if the shared one is absent (keeps older configs runnable)."""
    if "stage7_om_fraction" in mapping_cfg:
        return float(mapping_cfg["stage7_om_fraction"])
    soil_key = str(soil_name).strip().lower()
    frac_key = f"stage4_variantC_om_fraction_{soil_key}"
    if frac_key in mapping_cfg:
        return float(mapping_cfg[frac_key])
    return None


def _variantc_om_total(mapping_cfg, soil_name, porous_matrix_cells):
    """Stage 10 (soil_respiration_prompt_stage10_final.md, the final stage):
    when a `stage10_om_total_<soil>` is configured, it is used DIRECTLY,
    per soil -- Stage 10 explicitly sets aside the "shared/ranked scalars
    only" purity of Stages 7-9 for its two-soil (Vertisol, Sand) final
    deliverable, in favor of whatever total each soil needs to hit its
    target shape (Sand: strongly clears growth so a real burst fires;
    Vertisol: enough total to sustain a rise through the whole window).

    Stage 11 (soil_respiration_prompt_stage11_vertisol.md SS3) adds one
    more, highest-priority tier, `stage11_om_total_<soil>` -- Sand is
    FROZEN (SS1: no Stage 11 key for `sand`, so it keeps using its Stage 10
    value unchanged) while Vertisol's total is raised further.

    Falls back, in order, to: Stage 10's direct per-soil `stage10_om_
    total_<soil>`; Stage 9's ranked `stage9_om_base` x per-soil
    `stage9_om_multiplier_<soil>`; Stage 8's single ABSOLUTE `stage8_om_
    total` (identical for every soil); Stage 7's shared *fraction*
    (`_variantc_om_fraction`) -- so earlier stages stay independently
    reproducible."""
    soil_key = str(soil_name).strip().lower()
    # Stage 16 (soil_respiration_prompt_stage16_shared_om_rule.md SS3): the
    # dose LADDER, highest-priority tiers, tried in the order the prompt
    # sets out -- the point being to derive as much of the per-soil total
    # from STRUCTURE as possible rather than hand-setting it.
    #   rung 3  `stage16_om_total_<soil>`  -- per-soil, fully hand-set
    #           (Stage 10's position; the fallback, not the goal).
    #   rung 2  `stage16_om_rho_shared`    -- ONE shared OM concentration
    #           per unit matrix volume; the per-soil total is DERIVED,
    #           `om_total = rho * porous_matrix_cells`. `porous_matrix_
    #           cells` (not the full matrix) is used because that is where
    #           OM is physically placed; since `stage4_variantC_solid_
    #           fraction` is a fixed share of the matrix for every soil,
    #           the Vertisol:Sand ratio is identical either way (3.69x),
    #           so this choice only redefines rho's units, never the rung's
    #           behavior.
    #   rung 1  `stage16_om_total_shared`  -- ONE shared absolute total for
    #           both soils; maximum purity.
    # Listed highest-priority first so the runner adopts a rung simply by
    # setting its key; unset rungs fall through to the Stage 6-15 tiers
    # below, leaving every earlier stage independently reproducible.
    stage16_key = f"stage16_om_total_{soil_key}"
    if stage16_key in mapping_cfg:
        return float(mapping_cfg[stage16_key])
    if "stage16_om_rho_shared" in mapping_cfg:
        return float(mapping_cfg["stage16_om_rho_shared"]) * float(porous_matrix_cells)
    if "stage16_om_total_shared" in mapping_cfg:
        return float(mapping_cfg["stage16_om_total_shared"])
    stage11_key = f"stage11_om_total_{soil_key}"
    if stage11_key in mapping_cfg:
        return float(mapping_cfg[stage11_key])
    stage10_key = f"stage10_om_total_{soil_key}"
    if stage10_key in mapping_cfg:
        return float(mapping_cfg[stage10_key])
    if "stage9_om_base" in mapping_cfg:
        mult_key = f"stage9_om_multiplier_{soil_key}"
        multiplier = float(mapping_cfg.get(mult_key, 1.0))
        return float(mapping_cfg["stage9_om_base"]) * multiplier
    if "stage8_om_total" in mapping_cfg:
        return float(mapping_cfg["stage8_om_total"])
    om_frac = _variantc_om_fraction(mapping_cfg, soil_name)
    if om_frac is None:
        return float(mapping_cfg["OM_total"])
    return float(om_frac) * float(porous_matrix_cells)


def _variantc_D_scale(mapping_cfg, soil_name):
    """Stage 10 (soil_respiration_prompt_stage10_final.md SS2): a PER-SOIL
    transport-magnitude override, `stage10_D_scale_<soil>` -- breaking, for
    the final two-soil (Vertisol, Sand) deliverable only, Stages 7-9's
    invariant that `D_scale` is one SHARED constant. Vertisol's own reach
    needs to be SLOWED relative to Sand's so its many, dispersed habitat
    regions recruit gradually across the whole window instead of almost
    all at once (the fast, shared `D_scale` that Stage 7 needed just to get
    growth to fire at all otherwise front-loads every soil's recruitment,
    as Stage 9 found). Stage 11 (soil_respiration_prompt_stage11_vertisol.md
    SS3) adds `stage11_D_scale_<soil>` as a higher-priority tier still --
    Sand is FROZEN (no Stage 11 key for `sand`, so it falls through to its
    Stage 10 value unchanged) while Vertisol's is slowed further. Falls
    back, in order, to Stage 10's `stage10_D_scale_<soil>`, then the shared
    `stage7_D_scale`, so Stages 7-10 stay independently reproducible."""
    soil_key = str(soil_name).strip().lower()
    # Stage 16 SS3: the transport half of the dose ladder, same priority
    # convention as `_variantc_om_total` above -- `stage16_D_scale_shared`
    # is the ONE shared transport magnitude both soils get at rungs 1-2;
    # `stage16_D_scale_<soil>` is rung 3's per-soil fallback.
    stage16_key = f"stage16_D_scale_{soil_key}"
    if stage16_key in mapping_cfg:
        return float(mapping_cfg[stage16_key])
    if "stage16_D_scale_shared" in mapping_cfg:
        return float(mapping_cfg["stage16_D_scale_shared"])
    stage11_key = f"stage11_D_scale_{soil_key}"
    if stage11_key in mapping_cfg:
        return float(mapping_cfg[stage11_key])
    stage10_key = f"stage10_D_scale_{soil_key}"
    if stage10_key in mapping_cfg:
        return float(mapping_cfg[stage10_key])
    return float(mapping_cfg.get("stage7_D_scale", 1.0))


def _variantc_coating_fraction(mapping_cfg, soil_name, default=0.5):
    """Per-soil coating fraction -- the hypothesis-motivated share of a
    soil's OM placed in the fast coating pool rather than the slow-releasing
    trapped interior. Stage 9 (soil_respiration_prompt_stage9_ranked_om.md
    SS2-3) raises Vertisol's (more of its now-larger carbon actually reaches
    its many habitat regions, so their staggered recruitment can sustain a
    rise) and lowers Sand's (a smaller fast burst plus a slow low tail,
    combined with Sand's lowest total, so its accessible pool depletes
    within the window -- subject to the SS3 Exception: never lowered past
    the point where Sand's rise-then-fall burst itself breaks) from their
    Stage 6 values (Sand 0.80, Loess 0.50, Vertisol 0.30,
    `stage6_coating_fraction_<soil>`). Stage 10 (final stage, two soils
    only) adds one more, highest-priority tier, `stage10_coating_fraction_
    <soil>`, tuned directly per soil rather than ranking-only. Prefers
    Stage 10's key, then Stage 9's, else falls back to the Stage 6 key, so
    every earlier stage's configs stay independently reproducible. Stage 11
    (soil_respiration_prompt_stage11_vertisol.md SS3) adds one more,
    highest-priority tier, `stage11_coating_fraction_<soil>` -- Sand is
    FROZEN (no Stage 11 key for `sand`) while Vertisol's is raised further
    (more of its much-larger Stage 11 total reaches its habitat)."""
    soil_key = str(soil_name).strip().lower()
    stage11_key = f"stage11_coating_fraction_{soil_key}"
    if stage11_key in mapping_cfg:
        return float(mapping_cfg[stage11_key])
    stage10_key = f"stage10_coating_fraction_{soil_key}"
    if stage10_key in mapping_cfg:
        return float(mapping_cfg[stage10_key])
    stage9_key = f"stage9_coating_fraction_{soil_key}"
    if stage9_key in mapping_cfg:
        return float(mapping_cfg[stage9_key])
    frac_key = f"stage6_coating_fraction_{soil_key}"
    if frac_key in mapping_cfg:
        return float(mapping_cfg[frac_key])
    return float(default)


def leak_exchange(S, matrix_mask, macro_mask, k_leak):
    """Mobile-immobile exchange term (SS3): each matrix cell leaks dissolved
    substrate into any adjacent (4-connectivity) macropore cell at rate
    `k_leak`. Matrix-matrix and macropore-macropore faces contribute nothing.
    Built as a face-based flux (like model.conservative_divergence) so every
    unit leaving a matrix cell enters exactly one neighboring macropore cell
    -- exactly mass-conserving regardless of grid shape or masks."""
    exch = np.zeros_like(S)

    def _apply(axis):
        lo = (slice(None, -1), slice(None)) if axis == 0 else (slice(None), slice(None, -1))
        hi = (slice(1, None), slice(None)) if axis == 0 else (slice(None), slice(1, None))

        m_lo, mc_lo, S_lo = matrix_mask[lo], macro_mask[lo], S[lo]
        m_hi, mc_hi, S_hi = matrix_mask[hi], macro_mask[hi], S[hi]

        flux_lo_to_hi = np.where(m_lo & mc_hi, k_leak * S_lo, 0.0)  # matrix(lo) -> macro(hi)
        flux_hi_to_lo = np.where(m_hi & mc_lo, k_leak * S_hi, 0.0)  # matrix(hi) -> macro(lo)

        exch[lo] += flux_hi_to_lo - flux_lo_to_hi
        exch[hi] += flux_lo_to_hi - flux_hi_to_lo

    _apply(axis=0)
    _apply(axis=1)
    return exch


def _adjacent4(mask):
    """Boolean 4-neighbor dilation used to detect matrix-macro interfaces."""
    adj = np.zeros_like(mask, dtype=bool)
    adj[:-1, :] |= mask[1:, :]
    adj[1:, :] |= mask[:-1, :]
    adj[:, :-1] |= mask[:, 1:]
    adj[:, 1:] |= mask[:, :-1]
    return adj


def _dilate_steps(mask, steps):
    out = np.asarray(mask, dtype=bool)
    for _ in range(max(int(steps), 0)):
        out = out | _adjacent4(out)
    return out


def _component_labels(mask):
    """4-connected-component labeling (numpy-only BFS)."""
    mask = np.asarray(mask, dtype=bool)
    labels = np.zeros(mask.shape, dtype=int)
    ny, nx = mask.shape
    n = 0
    sizes = []
    for i0 in range(ny):
        for j0 in range(nx):
            if not mask[i0, j0] or labels[i0, j0] != 0:
                continue
            n += 1
            stack = [(i0, j0)]
            labels[i0, j0] = n
            size = 0
            while stack:
                i, j = stack.pop()
                size += 1
                for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ii, jj = i + di, j + dj
                    if 0 <= ii < ny and 0 <= jj < nx and mask[ii, jj] and labels[ii, jj] == 0:
                        labels[ii, jj] = n
                        stack.append((ii, jj))
            sizes.append(size)
    return labels, np.asarray(sizes, dtype=float)


def _cluster_scale_field(labels, sizes):
    """Map per-cluster sizes to [0, 1] and broadcast back to the grid."""
    scale = np.zeros_like(labels, dtype=float)
    if sizes.size == 0:
        return scale
    if sizes.size == 1 or float(sizes.max()) == float(sizes.min()):
        norm = np.ones_like(sizes, dtype=float)
    else:
        norm = (sizes - sizes.min()) / (sizes.max() - sizes.min())
    for k, v in enumerate(norm, start=1):
        scale[labels == k] = float(v)
    return scale


def _assign_solid_phase(matrix_mask, seed, solid_fraction, sigma):
    """Create contiguous solid mineral patches inside the matrix phase."""
    matrix_mask = np.asarray(matrix_mask, dtype=bool)
    if solid_fraction <= 0.0 or not matrix_mask.any():
        return np.zeros_like(matrix_mask, dtype=bool)
    if solid_fraction >= 1.0:
        return matrix_mask.copy()

    rng = np.random.default_rng(int(seed) + 7919)
    raw = rng.standard_normal(matrix_mask.shape)
    smooth = gaussian_filter_2d(raw, sigma=max(float(sigma), 1e-6))
    mat_vals = smooth[matrix_mask]
    thr = np.quantile(mat_vals, 1.0 - float(solid_fraction))
    solid = matrix_mask & (smooth >= thr)
    return solid


def _iter_offset_slices(shape, di, dj):
    """Return overlapping source/target slices for a grid shift (di, dj)."""
    ny, nx = shape
    if di >= 0:
        src_i = slice(0, ny - di)
        dst_i = slice(di, ny)
    else:
        src_i = slice(-di, ny)
        dst_i = slice(0, ny + di)
    if dj >= 0:
        src_j = slice(0, nx - dj)
        dst_j = slice(dj, nx)
    else:
        src_j = slice(-dj, nx)
        dst_j = slice(0, nx + dj)
    return (src_i, src_j), (dst_i, dst_j)


def apply_wicking(wet_mask, porous_matrix_mask, macro_mask, radius):
    """Hydro coupling: a dry macropore rewets if wet porous matrix is nearby."""
    wet_matrix = porous_matrix_mask & wet_mask
    influence = _dilate_steps(wet_matrix, radius)
    wickable_macro = macro_mask & ~wet_mask & influence
    return wet_mask | wickable_macro, wickable_macro


def wicking_substrate_exchange(S, porous_matrix_mask, macro_mask, wet_mask, k_wick, radius):
    """Directional transfer from wet porous matrix to wet macro within radius."""
    exch = np.zeros_like(S)
    radius = int(max(radius, 1))
    for di in range(-radius, radius + 1):
        for dj in range(-radius, radius + 1):
            dist = abs(di) + abs(dj)
            if dist == 0 or dist > radius:
                continue
            src, dst = _iter_offset_slices(S.shape, di, dj)
            m_src = porous_matrix_mask[src]
            mc_dst = macro_mask[dst]
            w_src = wet_mask[src]
            w_dst = wet_mask[dst]
            S_src = S[src]
            flux = np.where(m_src & mc_dst & w_src & w_dst, (k_wick / dist) * S_src, 0.0)
            exch[src] -= flux
            exch[dst] += flux
    return exch


def habitat_properties(macro_mask):
    """Connected macropore cavities define habitat networks."""
    labels, sizes = _component_labels(macro_mask)
    scale = _cluster_scale_field(labels, sizes)
    return labels, sizes, scale


def erosion_distance_4n(mask):
    """Stage 5 (prompts/soil_respiration_prompt_stage5_matric_inscribed.md
    Change 2): numpy-only iterative 4-neighbor erosion distance transform.
    `dist[i,j]` is the number of erosion rounds cell `(i,j)` survives before
    disappearing -- for a "+"-shaped (4-neighbor) structuring element this is
    exactly the Manhattan (L1) distance to the nearest non-mask cell OR the
    grid edge (out-of-grid is treated as background, via zero-padding), the
    standard numpy-only distance-to-boundary measure the prompt asks for."""
    mask = np.asarray(mask, dtype=bool)
    dist = np.zeros(mask.shape, dtype=int)
    current = mask.copy()
    depth = 0
    while current.any():
        depth += 1
        dist = np.where(current, depth, dist)
        padded = np.pad(current, 1, mode="constant", constant_values=False)
        current = (
            padded[1:-1, 1:-1] & padded[:-2, 1:-1] & padded[2:, 1:-1]
            & padded[1:-1, :-2] & padded[1:-1, 2:]
        )
    return dist


def macropore_region_inscribed_diameter(macro_mask, voxel_um):
    """Stage 5 Change 2: label connected macropore regions (4-connectivity,
    reusing `_component_labels`) and size EACH REGION by its maximum
    inscribed circle: `2 * max(erosion depth in cells) * voxel_um`. A large,
    well-connected region (a drained channel) gets a large diameter; a small,
    compact region (a water-holding pore) gets a small one -- replacing the
    Stage-4 flat 80 um constant that treated every macropore cell alike.

    Returns `(labels, region_diam_um, diam_field)`: `labels` as from
    `_component_labels`; `region_diam_um[k-1]` is region `k`'s inscribed
    diameter (1-indexed regions, 0-indexed array); `diam_field` broadcasts
    each macro cell's OWN region's diameter back onto the grid (0 elsewhere,
    i.e. on matrix/solid cells)."""
    labels, sizes = _component_labels(macro_mask)
    dist = erosion_distance_4n(macro_mask)
    n_regions = sizes.size
    region_diam = np.zeros(n_regions, dtype=float)
    for k in range(1, n_regions + 1):
        region_diam[k - 1] = 2.0 * float(dist[labels == k].max()) * voxel_um
    diam_field = np.zeros(macro_mask.shape, dtype=float)
    for k in range(1, n_regions + 1):
        diam_field[labels == k] = region_diam[k - 1]
    return labels, region_diam, diam_field


def bfs_geodesic_distance(seed_mask, allowed_mask):
    """Numpy-array BFS geodesic distance (4-connectivity), ported from
    `pde_sandbox.sandbox.bfs_distance`: distance (in grid steps) from the
    nearest seed cell, traveling only through `allowed_mask` cells. Returns
    `np.inf` where unreachable."""
    seed_mask = np.asarray(seed_mask, dtype=bool)
    allowed_mask = np.asarray(allowed_mask, dtype=bool)
    ny, nx = seed_mask.shape
    dist = np.full((ny, nx), np.inf, dtype=float)
    frontier = np.argwhere(seed_mask & allowed_mask)
    if len(frontier) == 0:
        return dist
    for i, j in frontier:
        dist[i, j] = 0.0
    queue = [(int(i), int(j)) for i, j in frontier]
    head = 0
    while head < len(queue):
        ci, cj = queue[head]
        head += 1
        base = dist[ci, cj] + 1.0
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ni, nj = ci + di, cj + dj
            if 0 <= ni < ny and 0 <= nj < nx and allowed_mask[ni, nj] and base < dist[ni, nj]:
                dist[ni, nj] = base
                queue.append((ni, nj))
    return dist


def place_om_shared_fine_pore(d_struct, porous_matrix_mask, om_fine_bias, om_total):
    """Stage 16 Change 1 (soil_respiration_prompt_stage16_shared_om_rule.md
    SS1): ONE shared, structure-blind organic-matter placement rule, used
    for EVERY soil, replacing Stage 6-15's per-soil coating/trapped split.

    This is the model's ORIGINAL rule (`MODEL_SPEC.md` SS4,
    `pore_field.om_field`): weight each cell by the fine-pore bias
    `d^-om_fine_bias` and normalize to that soil's total -- restricted here
    to the porous matrix, since OM is physically stored there and not in
    solids or open macropore voids (the same restriction
    `split_om_coating_trapped` already applied to both of its pools).

    What makes this the stage's central change is what it does NOT do: it
    never measures a cell's distance to a habitat. The Stage 6 coating
    construction did -- it defined the fast pool geodesically, which built
    an identical fuel shell around every macropore region at once and made
    all ~70 of Vertisol's regions self-fuelled from t=0 (the refuted
    Stage 14 attribution; see this module's docstring). Here the only thing
    that decides where carbon sits is the pore field itself, so the
    fuel-habitat geometry that results is emergent rather than imposed, and
    `analyze_fuel_habitat_distance` below measures it.

    Dissolution of this field is spatially UNIFORM at the single `k_dis`
    (there is no second rate), which is the other half of the invariant
    this stage enforces: no term anywhere in the model has a rate that
    depends on distance to habitat.

    Returns a single OM density field on the full grid (zero off the porous
    matrix)."""
    if not np.asarray(porous_matrix_mask, dtype=bool).any() or om_total <= 0.0:
        return np.zeros_like(d_struct)
    weight = np.where(porous_matrix_mask, d_struct ** (-om_fine_bias), 0.0)
    wsum = float(weight.sum())
    if wsum <= 0.0:
        # Degenerate pore field (no positive weight anywhere in the matrix):
        # fall back to a flat spread over the matrix rather than silently
        # placing zero carbon.
        return np.where(porous_matrix_mask, om_total / float(porous_matrix_mask.sum()), 0.0)
    return om_total * weight / wsum


def place_om_shared_lognormal_clumps(porous_matrix_mask, om_total, seed,
                                      n_clumps, log_mu, log_sigma):
    """Stage 16 SS7 -- the OPTIONAL, explicitly secondary sensitivity check:
    OM placed as discrete clumps whose radii are drawn from a SHARED
    lognormal with IDENTICAL parameters for both soils, instead of the
    smooth fine-pore field. The question it answers is narrow: does
    clumping PER SE sharpen hotspot selectivity beyond what the shared
    smooth rule already gives?

    Everything about the clump-size distribution is shared between soils --
    `n_clumps`, `log_mu`, `log_sigma` are single global parameters, so no
    per-soil carbon knob returns by the back door. Only WHERE the clumps
    can land differs, because clump centres are drawn from that soil's own
    porous-matrix cells; the pore field remains the sole structural input.

    Deliberately NOT parameterised from the companion project's measured
    POM data (SS7): that data is resolution-confounded between soils
    (Mishmar 5.85 um vs 15 um elsewhere) and its clustering work was
    retired as not robust in the companion project's own PROJECT_STATUS.md,
    and it does not exist at all for Rehovot/Sand. The model's structural
    inputs stay the pore-size distribution and its spatial arrangement.

    The clump seed is derived from the soil's own pore-field seed so the
    realization tracks that soil's structure reproducibly."""
    porous_matrix_mask = np.asarray(porous_matrix_mask, dtype=bool)
    OM = np.zeros(porous_matrix_mask.shape, dtype=float)
    candidates = np.argwhere(porous_matrix_mask)
    if candidates.shape[0] == 0 or om_total <= 0.0 or n_clumps <= 0:
        return OM

    rng = np.random.default_rng(int(seed) + 16000)
    picks = rng.integers(0, candidates.shape[0], size=int(n_clumps))
    radii = rng.lognormal(mean=log_mu, sigma=log_sigma, size=int(n_clumps))

    ny, nx = porous_matrix_mask.shape
    yy, xx = np.ogrid[:ny, :nx]
    for idx, radius in zip(picks, radii):
        ci, cj = candidates[idx]
        rad = max(float(radius), 0.5)
        disc = ((yy - ci) ** 2 + (xx - cj) ** 2) <= rad * rad
        cells = disc & porous_matrix_mask
        if not cells.any():
            # Clump centre sits in an isolated matrix cell: keep its mass
            # on that one cell rather than discarding it, so the soil's
            # total is conserved.
            cells = np.zeros_like(porous_matrix_mask)
            cells[ci, cj] = True
        # Mass proportional to the clump's own area, spread evenly inside
        # it -- a clump is a lump of particulate carbon, not a point.
        OM[cells] += (rad * rad) / float(cells.sum())

    total = float(OM.sum())
    if total <= 0.0:
        return np.where(porous_matrix_mask, om_total / float(porous_matrix_mask.sum()), 0.0)
    return OM * (om_total / total)


def substantial_om_mask(OM, porous_matrix_mask):
    """Stage 16 SS2's stated definition of "substantial" organic matter: a
    porous-matrix cell whose OM mass is at or above the MEDIAN of all
    NONZERO matrix OM in that soil.

    Stated plainly because the choice matters and the prompt asks for it to
    be defensible and used consistently: the median of the nonzero
    population is (a) scale-free -- it does not move when a ladder rung
    changes `om_total`, since the whole field scales by one factor, so the
    distance distributions it produces are a property of STRUCTURE alone
    and are directly comparable between rungs and between soils; and (b)
    robust to the extreme right skew the `d^-om_fine_bias` rule produces
    (Stage 3.5 measured ~99.96% of carbon landing on sub-micron cells under
    the single-continuum version of this rule), where a mean-based
    threshold would sit far out in the tail and select almost nothing."""
    OM = np.asarray(OM, dtype=float)
    nonzero = OM[porous_matrix_mask & (OM > 0.0)]
    if nonzero.size == 0:
        return np.zeros(OM.shape, dtype=bool), float("nan")
    threshold = float(np.median(nonzero))
    return (porous_matrix_mask & (OM >= threshold)), threshold


def analyze_fuel_habitat_distance(grids, om_field_override=None):
    """Stage 16 SS2 -- THIS STAGE'S PRIMARY STRUCTURAL RESULT: the
    distribution, over macropore REGIONS, of the geodesic distance from a
    region to its nearest SUBSTANTIAL organic matter.

    Method, per soil: mark substantial OM (`substantial_om_mask`), then run
    `bfs_geodesic_distance` OUT FROM that OM through everything that is not
    solid (porous matrix + macropore cells -- the phases substrate can
    actually move through; solids are barriers). Each macropore region's
    distance is the MINIMUM over its own cells, i.e. how far that region's
    nearest fuel is. Regions with no reachable substantial OM are reported
    separately as unreachable rather than being folded in as a huge finite
    number, which would silently bias every percentile.

    Distances are returned both in grid cells and in microns (via
    `voxel_um`), since the external comparison this stage is measured
    against (Portell et al. 2018, Front. Microbiol. 9:1583 -- operative
    geodesic distances of 252 +/- 111 um and 490 +/- 262 um) is in microns.

    Returns a dict with the full distribution (min/p25/median/p75/p95/max,
    mean, std) in both units, the per-region distances, and the counts."""
    porous_matrix_mask = grids["porous_matrix_mask"]
    macro_mask = grids["macro_mask"]
    solid_mask = grids["solid_mask"]
    region_labels = grids["region_labels"]
    voxel_um = float(grids["voxel_um"])
    OM = grids["OM"] if om_field_override is None else om_field_override

    sub_mask, threshold = substantial_om_mask(OM, porous_matrix_mask)
    # Substrate travels through the pore space of both phases; only the
    # solid phase is a barrier. This matches the transport the model
    # actually solves (matrix internal diffusion + macropore D + k_leak),
    # rather than inventing a different connectivity for the diagnostic.
    allowed = (porous_matrix_mask | macro_mask) & ~solid_mask
    dist = bfs_geodesic_distance(sub_mask, allowed)

    n_regions = int(region_labels.max())
    per_region = np.full(n_regions, np.inf, dtype=float)
    for k in range(1, n_regions + 1):
        cells = region_labels == k
        if not cells.any():
            continue
        d_cells = dist[cells]
        finite = d_cells[np.isfinite(d_cells)]
        if finite.size:
            per_region[k - 1] = float(finite.min())

    finite_mask = np.isfinite(per_region)
    finite_vals = per_region[finite_mask]

    def _pct(q):
        return float(np.percentile(finite_vals, q)) if finite_vals.size else float("nan")

    stats_cells = {
        "min": float(finite_vals.min()) if finite_vals.size else float("nan"),
        "p25": _pct(25), "median": _pct(50), "p75": _pct(75), "p95": _pct(95),
        "max": float(finite_vals.max()) if finite_vals.size else float("nan"),
        "mean": float(finite_vals.mean()) if finite_vals.size else float("nan"),
        "std": float(finite_vals.std()) if finite_vals.size else float("nan"),
    }
    return {
        "per_region_cells": per_region,
        "per_region_um": per_region * voxel_um,
        "stats_cells": stats_cells,
        "stats_um": {k: v * voxel_um for k, v in stats_cells.items()},
        "substantial_threshold": threshold,
        "substantial_cells": int(sub_mask.sum()),
        "porous_matrix_cells": int(porous_matrix_mask.sum()),
        "n_regions": n_regions,
        "n_regions_reachable": int(finite_mask.sum()),
        "n_regions_unreachable": int((~finite_mask).sum()),
        "voxel_um": voxel_um,
    }


def split_om_coating_trapped(d_struct, porous_matrix_mask, macro_mask, water_mask,
                              om_fine_bias, om_total, coating_fraction, coating_thickness):
    """Stage 6 Change 2 (soil_respiration_prompt_stage6_portback.md), ported
    from `pde_sandbox.sandbox.build_down_source`'s geodesic source
    construction: split a soil's total OM into a COATING pool (porous-matrix
    cells within `coating_thickness` geodesic steps of a habitat/macropore
    cell, the distance measured along a WET connected matrix path -- fast,
    a guaranteed finite connected habitat<-source path) and a TRAPPED pool
    (everything else in the porous matrix -- interior cells far from habitat
    or behind dry/solid barriers). `coating_fraction` is the per-soil,
    hypothesis-motivated share of OM in the coating. Both pools are weighted
    internally by the same fine-pore bias rule as `pore_field.om_field`, just
    renormalized within each pool.

    The macropore cell itself need not be WET for a neighboring matrix cell
    to count as coating -- "along a wet connected path" gates the PATH
    (the matrix, which Stage 5's matric rule keeps essentially always wet),
    not the destination habitat cell. Requiring the destination to be wet
    too would make a soil whose entire macropore network is dry under the
    static matric map (Sand, per Stage 5) have literally zero coating
    regardless of its configured coating fraction, which would silently
    defeat the coating hypothesis for exactly the soil it is meant to test --
    the model already treats "wet" purely as a transport gate (K/growth
    never check water_mask directly), so this is consistent with that
    convention, not a special case.

    Returns `(OM_coating, OM_trapped, coating_mask, trapped_mask)`."""
    wet_matrix = porous_matrix_mask & water_mask
    macro_adjacent = _adjacent4(macro_mask)
    seed = wet_matrix & macro_adjacent
    geodesic = bfs_geodesic_distance(seed, wet_matrix)

    coating_mask = porous_matrix_mask & np.isfinite(geodesic) & (geodesic <= coating_thickness)
    if not coating_mask.any():
        # No wet connected path this realization (e.g. the whole matrix is
        # dry) -- fall back to cells directly touching a macropore, so
        # "coating" still means something even when the wet network is this
        # fragmented, rather than silently placing zero coating OM.
        coating_mask = porous_matrix_mask & macro_adjacent
    trapped_mask = porous_matrix_mask & ~coating_mask

    def _weighted(mask, total):
        if not mask.any() or total <= 0.0:
            return np.zeros_like(d_struct)
        weight = np.where(mask, d_struct ** (-om_fine_bias), 0.0)
        wsum = float(weight.sum())
        if wsum <= 0.0:
            return np.where(mask, total / float(mask.sum()), 0.0)
        return total * weight / wsum

    om_coating_total = om_total * coating_fraction
    om_trapped_total = om_total * (1.0 - coating_fraction)
    OM_coating = _weighted(coating_mask, om_coating_total)
    OM_trapped = _weighted(trapped_mask, om_trapped_total)
    return OM_coating, OM_trapped, coating_mask, trapped_mask


def _prepare_implicit_cache(D, dx):
    """Face conductivities (harmonic mean) + per-cell diagonal weight for
    `implicit_diffusion_step` below -- ported from
    `pde_sandbox.sandbox.prepare_local_cache`, using the repo's own
    `model.harmonic_mean`/dx convention."""
    east = harmonic_mean(D[:, :-1], D[:, 1:])
    south = harmonic_mean(D[:-1, :], D[1:, :])
    diag_weight = np.zeros_like(D)
    diag_weight[:, :-1] += east
    diag_weight[:, 1:] += east
    diag_weight[:-1, :] += south
    diag_weight[1:, :] += south
    return east, south, diag_weight, dx * dx


def _implicit_neighbor_sum(u, east, south):
    out = np.zeros_like(u)
    out[:, :-1] += east * u[:, 1:]
    out[:, 1:] += east * u[:, :-1]
    out[:-1, :] += south * u[1:, :]
    out[1:, :] += south * u[:-1, :]
    return out


def implicit_diffusion_step(u_rhs, D, dx, dt, iterations, relaxation, tol=None):
    """Stage 6 Change 1 (soil_respiration_prompt_stage6_portback.md), ported
    directly from `pde_sandbox.sandbox.local_implicit_step`: solves
    `u - dt*div(D*grad(u)) = u_rhs` by Jacobi iteration with harmonic-mean
    face conductivities (so a dry/solid cell, D=0, is a true barrier -- its
    face conductivity is exactly 0). Unlike the explicit conservative-flux
    update used through Stage 5, this has NO `max(D)*dt/dx^2 < 0.2` stability
    cap -- that cap was exactly what kept dissolved substrate frozen ~1-3
    cells from where it dissolved (the recurring Stage 3.5/4/5 diagnosis).

    Stage 11 (soil_respiration_prompt_stage11_vertisol.md SS2) adds an
    optional convergence-tolerance early stop (`tol`, max absolute change
    between successive iterates): the fixed `iterations=90` budget was
    conservative headroom, not a requirement -- most steps converge far
    sooner, so stopping as soon as the change per iteration drops below
    `tol` (the `iterations` argument becomes a safety CAP, not a fixed
    count) cuts wasted Jacobi work without changing the converged answer.
    `tol=None` reproduces the exact old fixed-iteration behavior."""
    east, south, diag_weight, dx2 = _prepare_implicit_cache(D, dx)
    diag = 1.0 + dt * diag_weight / dx2
    u = u_rhs.copy()
    for _ in range(iterations):
        rhs = u_rhs + (dt / dx2) * _implicit_neighbor_sum(u, east, south)
        proposal = rhs / diag
        u_new = relaxation * proposal + (1.0 - relaxation) * u
        u_new = np.clip(u_new, 0.0, None)
        if tol is not None:
            delta = float(np.abs(u_new - u).max())
            u = u_new
            if delta < tol:
                break
        else:
            u = u_new
    return u


def build_dual_grids(structure_cfg, mapping_cfg):
    """Build the macro/matrix split and all derived fields for one soil.
    `structure_cfg["pore"]` must be a `psd_parametric` config (the FULL,
    untruncated literature PSD -- Variant C's whole point is that the
    sub-resolution fraction stays on the grid, just cut off from fast
    transport, rather than being dropped as Variant A does)."""
    n = structure_cfg["grid"]["n"]
    dx = structure_cfg["grid"]["dx"]
    pore_cfg = structure_cfg["pore"]
    pore_mode = pore_cfg.get("mode", "lognormal")
    assert pore_mode in ("psd_parametric", "psd_measured"), \
        "dual_porosity requires the full literature PSD (pore.mode: psd_parametric) " \
        "or the real measured PSD (pore.mode: psd_measured, Stage 13)"

    if pore_mode == "psd_measured":
        # Stage 13 (soil_respiration_prompt_stage13_real_data.md): real,
        # measured PSD from the companion project's CT/segmentation
        # pipeline, replacing the literature mixture -- SAME downstream
        # shape (psd_measured.build_soil_psd mirrors psd_parametric.
        # build_soil_psd's return dict exactly), so nothing below this
        # point needs to know which source built `psd`.
        psd = psd_measured.build_soil_psd(
            pore_cfg["psd_soil"], pore_cfg.get("psd_config", "configs/psd_measured.yaml"))
    else:
        psd = psd_parametric.build_soil_psd(
            pore_cfg["psd_soil"], pore_cfg.get("psd_config", "configs/psd_literature.yaml"))
    d = generate_diameter_field_from_psd(
        n, pore_cfg["lambda"], pore_cfg["seed"],
        psd["bin_edges_um"], psd["cdf_edges"], pore_cfg.get("aggregate", False))

    d_split = mapping_cfg["K_d_low"]  # 30 um macro/matrix split threshold
    macro_mask = d >= d_split
    matrix_mask = ~macro_mask

    solid_fraction = float(mapping_cfg.get("stage4_variantC_solid_fraction", 0.35))
    solid_sigma = float(mapping_cfg.get("stage4_variantC_solid_patch_sigma", 2.0))
    solid_mask = _assign_solid_phase(matrix_mask, pore_cfg["seed"], solid_fraction, solid_sigma)
    porous_matrix_mask = matrix_mask & ~solid_mask

    # Stage 5 Change 2 (soil_respiration_prompt_stage5_matric_inscribed.md):
    # size each connected macropore REGION by its own maximum inscribed
    # circle instead of forcing every macropore cell to one flat diameter --
    # a 40 um water-holding pore and a 500 um drained channel are no longer
    # treated identically by K(d)/D(d).
    voxel_um = float(mapping_cfg.get("voxel_um", 10.0))
    region_labels, region_diam_um, region_diam_field = macropore_region_inscribed_diameter(
        macro_mask, voxel_um)
    d_struct = np.where(macro_mask, region_diam_field, d)

    K = K_window(d_struct, mapping_cfg)
    K = np.where(macro_mask, K, 0.0)
    D_macro_base = D_of_d(d_struct, mapping_cfg)
    D_macro_base = np.where(macro_mask, D_macro_base, 0.0)

    sat_cfg = structure_cfg.get("saturation", {"mode": "fully_wet"})
    if sat_cfg.get("mode") == "matric_regions":
        # Stage 5 Change 1: a SHARED air-entry diameter for every soil (the
        # matric potential). Matrix cells use their own PSD diameter d (all
        # < K_d_low=30um < the cutoff, so always wet -- the fine matrix stays
        # saturated); macropore cells use their REGION's inscribed diameter
        # (NOT the per-cell PSD label, NOT a forced constant) -- a small
        # region stays water-filled, a large one drains and goes dry. Solid
        # cells are never wet. theta then EMERGES (not imposed).
        d_cut = float(sat_cfg.get("matric_d_cut_um", mapping_cfg.get("matric_d_cut_um", 50.0)))
        size_for_wetting = np.where(macro_mask, region_diam_field, d)
        water_mask = size_for_wetting < d_cut
        theta = float(water_mask.mean())
    else:
        water_mask, theta = water_mask_and_theta(d, sat_cfg, mapping_cfg)
    water_mask = water_mask & ~solid_mask

    # Stage 7 Change 3 (soil_respiration_prompt_stage7_calibration.md SS4):
    # a DRAINED macropore (macro cell the matric rule marked dry) is not
    # fully cut off at constant high humidity -- it retains a thin wall
    # water film. Give it a small nonzero film diffusivity (`D_film`,
    # applied in simulate_dual_porosity's per-step D_macro_t so the shared
    # `stage7_D_scale` calibration lever also reaches it) and let it host
    # activity at REDUCED capacity (`stage7_film_K_scale`, e.g. Sand's one
    # huge dry region goes from a total habitat dead zone to a fast-but-low-
    # capacity fringe). K(d) itself is unaffected by wet/dry status
    # (already the case before Stage 7 -- only D was ever gated by water),
    # so this only needs to scale K down on dry macro cells, not build a new
    # habitat rule.
    dry_macro_mask = macro_mask & ~water_mask
    film_k_scale = float(mapping_cfg.get("stage7_film_K_scale", 1.0))
    K = np.where(dry_macro_mask, K * film_k_scale, K)

    D_film = float(mapping_cfg.get("D_film", 0.0))
    D = np.where(water_mask, D_macro_base, np.where(dry_macro_mask, D_film, 0.0))

    # Fine-pore-biased OM over the FULL PSD, physically stored only inside the
    # porous matrix phase (not in solids and not in open macropore voids).
    # Stage 8 Change 1: `_variantc_om_total` returns a single ABSOLUTE total,
    # IDENTICAL for every soil, when `stage8_om_total` is configured -- the
    # deliberate lever this stage tests: the SAME total carbon, concentrated
    # on Sand's small coating vs. spread across Vertisol's many habitats,
    # should deplete on very different timescales. The coating/trapped
    # SPLIT (below) still follows each soil's own geometry -- only the total
    # mass is equal.
    soil_name = pore_cfg.get("psd_soil", "")
    om_total = _variantc_om_total(mapping_cfg, soil_name, porous_matrix_mask.sum())

    # Stage 6 Change 2/3 (soil_respiration_prompt_stage6_portback.md): split
    # that total OM into a fast COATING pool (near habitat, along a wet path)
    # and a slow-releasing TRAPPED pool (interior/behind barriers), per-soil
    # coating fraction. Needs `water_mask`, so this comes after Change 1.
    # Stage 16 Change 1 (soil_respiration_prompt_stage16_shared_om_rule.md
    # SS1) makes this a branch: `shared_fine_pore` (the new active default)
    # places ALL of a soil's OM by the one shared fine-pore rule with NO
    # geodesic construction and NO per-soil share, while `coating_trapped`
    # keeps the Stage 6-15 path reachable and byte-identical.
    om_mode = _om_mode(mapping_cfg)
    coating_thickness = int(mapping_cfg.get("stage6_coating_thickness", 2))
    if om_mode in (OM_MODE_SHARED, OM_MODE_CLUMPS):
        # The whole field dissolves at the single, uniform `k_dis`. It is
        # carried in the `OM_coating` slot purely because that is the slot
        # `simulate_dual_porosity`'s source term already releases at
        # `k_dis`; the legacy NAME is kept so nothing downstream has to
        # branch, but there is no coating here in any physical sense --
        # `OM_trapped` is identically zero, so the `k_dis_slow` term
        # vanishes from the equations entirely and no rate anywhere depends
        # on distance to habitat.
        if om_mode == OM_MODE_CLUMPS:
            OM_coating = place_om_shared_lognormal_clumps(
                porous_matrix_mask, max(om_total, 0.0), pore_cfg["seed"],
                int(mapping_cfg.get("stage16_clump_count", 60)),
                float(mapping_cfg.get("stage16_clump_log_mu", 0.7)),
                float(mapping_cfg.get("stage16_clump_log_sigma", 0.6)))
        else:
            OM_coating = place_om_shared_fine_pore(
                d_struct, porous_matrix_mask, mapping_cfg["om_fine_bias"], max(om_total, 0.0))
        OM_trapped = np.zeros_like(OM_coating)
        coating_mask = porous_matrix_mask.copy()
        trapped_mask = np.zeros_like(porous_matrix_mask, dtype=bool)
        # NaN, not 1.0 or 0.0: `coating_fraction` is RETIRED under this
        # mode, and a NaN makes any downstream code that still reads it
        # visibly not-applicable rather than quietly reporting a number
        # that no longer means anything.
        coating_fraction = float("nan")
    else:
        coating_fraction = _variantc_coating_fraction(mapping_cfg, soil_name)
        OM_coating, OM_trapped, coating_mask, trapped_mask = split_om_coating_trapped(
            d_struct, porous_matrix_mask, macro_mask, water_mask,
            mapping_cfg["om_fine_bias"], max(om_total, 0.0),
            coating_fraction, coating_thickness,
        )
    OM = OM_coating + OM_trapped  # combined field, kept for metrics.py's
                                  # connectivity diagnostics (om_biomass_
                                  # connectivity / accessible_om_fraction),
                                  # which only need a single OM density field

    B0 = b0_field(K, mapping_cfg["B0_total"])
    S0 = np.zeros((n, n))

    labels, sizes, scale = habitat_properties(macro_mask)

    return {
        "n": n, "dx": dx, "r": d, "d": d, "K": K,
        "D_full": D_macro_base, "D": D,
        "water_mask": water_mask, "theta": theta,
        "macro_mask": macro_mask,
        "matrix_mask": matrix_mask,
        "porous_matrix_mask": porous_matrix_mask,
        "solid_mask": solid_mask,
        "dry_macro_mask": dry_macro_mask,
        "dry_macro_cells": int(dry_macro_mask.sum()),
        "dry_macro_fraction_of_habitat": (
            float(dry_macro_mask.sum()) / float(macro_mask.sum()) if macro_mask.any() else 0.0
        ),
        "D_film": D_film,
        "matrix_fraction": float(matrix_mask.mean()),
        "porous_matrix_fraction": float(porous_matrix_mask.mean()),
        "solid_fraction": float(solid_mask.mean()),
        "om_fraction_init": (
            float(mapping_cfg["stage7_om_fraction"]) if "stage7_om_fraction" in mapping_cfg
            and "stage8_om_total" not in mapping_cfg else None
        ),
        "om_total_init": float(om_total),
        "om_mode": om_mode,
        "om_total_placed": float(OM.sum()),
        "coating_fraction": float(coating_fraction),
        "coating_mask": coating_mask,
        "trapped_mask": trapped_mask,
        "habitat_labels": labels,
        "habitat_scale": scale,
        "habitat_count": int(sizes.size),
        "mean_habitat_size": float(sizes.mean()) if sizes.size else 0.0,
        "voxel_um": voxel_um,
        "region_labels": region_labels,
        "region_diam_um": region_diam_um,
        "region_diam_field": region_diam_field,
        "OM": OM, "OM_coating": OM_coating, "OM_trapped": OM_trapped,
        "B0": B0, "S0": S0, "psd_truncation": None,
    }


def count_active_habitats(growth, region_labels, growth_threshold=1e-9):
    """Stage 6 SS6 'habitats recruited over time': count of distinct
    macropore REGIONS (Stage 5's inscribed-diameter regions, `region_labels`)
    with at least one actively-growing cell this step."""
    active_cells = growth > growth_threshold
    if not active_cells.any():
        return 0
    active_labels = np.unique(region_labels[active_cells])
    return int((active_labels > 0).sum())


def simulate_dual_porosity(biology_cfg, structure_cfg, mapping_cfg, wicking_enabled=False,
                            debug_growth_distribution=False, track_region_biomass=False):
    """Integrate the dual-porosity model for one soil. Same per-cell biology
    and respiration observable as model.simulate.

    Stage 6 (soil_respiration_prompt_stage6_portback.md) replaces the
    explicit conservative-flux substrate transport with the sandbox-ported
    `implicit_diffusion_step` (Change 1: no more `max(D)*dt/dx^2 < 0.2` cap),
    splits OM into two pools -- `OM_coating` (dissolves at `k_dis`) and
    `OM_trapped` (dissolves at the much slower `k_dis_slow`, Change 3) --
    and removes the `S_cap` substrate ceiling entirely (Change 4). `k_leak`
    (matrix->macropore mobile-immobile exchange) and matrix internal
    diffusion (folded into the SAME implicit solve as the macropore D field,
    since they're now solved together rather than as two separate blocked-
    off explicit divergences) remain as supporting pathways.

    `wicking_enabled` (Stage 5 mechanism, default OFF as of Stage 6 SS5 --
    the sandbox this stage ports from used a static structure, and wicking
    previously overrode the Stage-5 wetting geometry): the dynamic
    hydro-wicking mechanism (`apply_wicking` + `wicking_substrate_exchange`)
    can rewet a macropore region the static matric map marked dry. Set True
    to re-enable it (e.g. for a direct Stage-5-vs-6 comparison).

    `debug_growth_distribution` (Stage 14, soil_respiration_prompt_
    stage14_habitat_diagnostic.md SS1, default OFF): when True, records the
    raw per-cell `growth` value (not just the `count_active_habitats`
    boolean) at every timestep, restricted to macropore cells
    (`macro_mask`), as `growth_distribution_t` (shape `(T, n_macro_cells)`)
    in the returned dict, plus `macro_mask` and `region_labels` for
    indexing it back onto regions. PURELY OBSERVATIONAL -- this flag adds
    one extra assignment per step and one extra output key; it does not
    read, branch on, or alter any state the physics depends on, so it
    cannot change `R_t`/`B_total_t`/any other existing output (verified
    directly in stage14_run.py's identity check, not just asserted here).

    Stage 16 (soil_respiration_prompt_stage16_shared_om_rule.md SS4c) adds
    two more OBSERVATIONAL outputs, needed for the hotspot-selectivity and
    recruitment-timing criteria:

    - `active_mass_t` (always recorded): the growth-weighted "active mass"
      metric Stage 14 recommended -- total raw `growth` summed over all
      macropore cells each step. Stage 14 measured this at `r=0.996`
      against `R(t)`'s own shape, versus the boolean
      `habitats_recruited_t` count's `r=0.072`, so Stage 16 reports it as
      the headline recruitment diagnostic and keeps the boolean count only
      as a secondary one.
    - `region_biomass_t` (opt-in via `track_region_biomass`, shape
      `(T, n_regions)`): per-macropore-region total biomass each step, via
      one `np.bincount` over `region_labels` per step. Needed for the
      per-region PEAK biomass that Portell et al. (2018)'s hotspot metric
      is defined on (share of regions whose peak exceeds 10% of the largest
      region's peak).

    Both are pure read-only side effects of state the loop already
    computes, exactly like `debug_growth_distribution`: neither is read
    back by the physics, so neither can alter `R_t`/`B_total_t`/any
    existing output (verified by stage16_run.py's own identity check, not
    merely asserted)."""
    grids = build_dual_grids(structure_cfg, mapping_cfg)
    dx = grids["dx"]
    # Stage 7 Change 2 (soil_respiration_prompt_stage7_calibration.md SS3): a
    # post-scale multiplier on the transport magnitude -- raises D(d)'s
    # absolute reach (sqrt(D*T)) without touching D_d_ref (which stays
    # Stage 3.5's soil-discriminating SHAPE of D(d)). Stage 7-9 kept this
    # ONE constant shared by every soil; Stage 10 (final stage)
    # `_variantc_D_scale` allows a PER-SOIL override so Vertisol's reach can
    # be deliberately slowed relative to Sand's, for the two-soil final
    # deliverable only.
    soil_name = structure_cfg["pore"].get("psd_soil", "")
    D_scale = _variantc_D_scale(mapping_cfg, soil_name)
    D_macro_base, K_base = grids["D_full"] * D_scale, grids["K"]
    dry_macro_mask = grids["dry_macro_mask"]
    D_film = float(mapping_cfg.get("D_film", 0.0)) * D_scale
    wet_mask = grids["water_mask"].copy()
    macro_mask = grids["macro_mask"]
    porous_matrix_mask = grids["porous_matrix_mask"]
    region_labels = grids["region_labels"]
    k_leak = mapping_cfg["k_leak"]
    k_wick = float(mapping_cfg.get("k_wick", k_leak))
    wick_radius = int(mapping_cfg.get("stage4_variantC_wicking_radius", 2))
    k_matrix_diff = float(mapping_cfg.get("stage4_variantC_matrix_diffusion", 0.03)) * D_scale

    d_cluster_scale = float(mapping_cfg.get("stage4_variantC_cluster_D_scale", 0.75))
    k_cluster_scale = float(mapping_cfg.get("stage4_variantC_cluster_K_scale", 0.75))

    implicit_iterations = int(mapping_cfg.get("stage6_implicit_iterations", 90))
    implicit_relaxation = float(mapping_cfg.get("stage6_implicit_relaxation", 0.85))
    # Stage 11 SS2: `stage6_implicit_iterations` becomes a safety CAP once a
    # convergence tolerance is set; the solve stops as soon as consecutive
    # iterates change by less than `stage11_implicit_tol`, instead of always
    # spending the full fixed budget. `None` (the default if unset)
    # reproduces the exact pre-Stage-11 fixed-iteration behavior.
    implicit_tol = mapping_cfg.get("stage11_implicit_tol")
    implicit_tol = float(implicit_tol) if implicit_tol is not None else None
    k_dis_slow = float(mapping_cfg.get("k_dis_slow", 0.005))
    # Stage 16 SS1: under the shared rule, dissolution is spatially UNIFORM
    # -- one `k_dis` everywhere, no second rate. The slow pool is already
    # identically zero under this mode (`build_dual_grids`), so the
    # `k_dis_slow * OM_trapped` term contributes nothing either way; tying
    # `k_dis_slow` to `k_dis` here makes that an ENFORCED invariant of the
    # mode rather than a property that merely happens to hold because an
    # array is empty, so no future edit can reintroduce a second,
    # spatially-distinct decomposition rate through this path unnoticed.
    if _om_mode(mapping_cfg) in (OM_MODE_SHARED, OM_MODE_CLUMPS):
        k_dis_slow = float(biology_cfg["k_dis"])

    r_max = biology_cfg["r_max"]
    Ks = biology_cfg["Ks"]
    Y = biology_cfg["Y"]
    m0 = biology_cfg["m0"]
    m_s = biology_cfg["m_s"]
    k_dis = biology_cfg["k_dis"]
    dt = biology_cfg["dt"]
    T = structure_cfg["T"]

    # No explicit-diffusion stability cap to check (Change 1) -- the implicit
    # Jacobi solve is unconditionally stable regardless of D, dt, dx. Only
    # the remaining EXPLICIT exchange terms (leak, and wicking if enabled)
    # still need their own bound.
    assert dt * 4.0 * (k_leak + (k_wick if wicking_enabled else 0.0)) < 1.0, \
        "exchange stability violated: reduce k_leak/k_wick or dt"

    B = grids["B0"].copy()
    S = grids["S0"].copy()
    OM_coating = grids["OM_coating"].copy()
    OM_trapped = grids["OM_trapped"].copy()

    om_total_init = float(grids["om_total_init"])

    R_t = np.zeros(T)
    cum_co2 = np.zeros(T)
    S_habitat_mean_t = np.zeros(T)
    habitat_count_t = np.zeros(T)
    habitats_recruited_t = np.zeros(T)
    fuel_remaining_t = np.zeros(T)
    B_total_t = np.zeros(T)
    growth_term_t = np.zeros(T)
    maintenance_term_t = np.zeros(T)
    running = 0.0

    if debug_growth_distribution:
        n_macro_cells = int(macro_mask.sum())
        growth_distribution_t = np.zeros((T, n_macro_cells), dtype=float)

    # Stage 16 SS4c observational outputs (see docstring).
    active_mass_t = np.zeros(T)
    n_regions_total = int(region_labels.max())
    if track_region_biomass:
        region_biomass_t = np.zeros((T, n_regions_total), dtype=float)
    region_labels_flat = region_labels.ravel()

    # Stage 11 SS2 (soil_respiration_prompt_stage11_vertisol.md): the
    # connected macropore regions, their inscribed sizes, and the
    # habitat-cluster scale depend ONLY on `macro_mask`, which is fixed for
    # the whole run (wicking, even when enabled, only ever touches
    # `wet_mask`, never `macro_mask`) -- so `habitat_properties` was pure
    # wasted work recomputed every single timestep (a pure-Python
    # connected-component BFS over the whole n x n grid, T times). Computed
    # ONCE here and reused for every step; this was the dominant per-step
    # cost.
    habitat_labels, habitat_sizes, habitat_scale = habitat_properties(macro_mask)
    habitat_count_value = float(habitat_sizes.size)
    habitat_labels_last = habitat_labels
    habitat_scale_last = habitat_scale

    # `K_t`/`habitable`/`K_safe` depend only on `macro_mask`/`K_base`/
    # `habitat_scale`, all loop-invariant (Stage 11 SS2: hoisted alongside
    # `habitat_properties` for the same reason -- nothing here changes
    # per-timestep since wicking never touches `macro_mask`).
    K_t = np.where(macro_mask, K_base * (1.0 + k_cluster_scale * habitat_scale), 0.0)
    habitable = K_t > 0
    K_safe = np.where(habitable, K_t, 1.0)

    for t in range(T):
        habitat_count_t[t] = habitat_count_value

        # Dynamic hydro coupling: OFF by default (Stage 6 SS5) so the static
        # matric-potential wet/dry map stays exactly as build_dual_grids made
        # it, matching the sandbox's static structure.
        if wicking_enabled:
            wet_mask, _ = apply_wicking(wet_mask, porous_matrix_mask, macro_mask, wick_radius)

        # Stage 7 Change 3: a DRAINED macropore cell (macro_mask & ~wet_mask,
        # recomputed each step so wicking, if enabled, still tracks it
        # correctly) is not fully cut off (D=0) any more -- it keeps a small
        # water-film pathway D_film instead, so it can receive/exchange
        # substrate at reduced rate instead of being wholly isolated.
        dry_macro_t = macro_mask & ~wet_mask
        D_macro_t = np.where(wet_mask, D_macro_base * (1.0 + d_cluster_scale * habitat_scale), 0.0)
        D_macro_t = np.where(dry_macro_t, D_film, D_macro_t)
        D_matrix_t = np.where(wet_mask & porous_matrix_mask, k_matrix_diff, 0.0)
        D_total = D_macro_t + D_matrix_t  # disjoint supports (macro vs matrix) -> a plain sum is the per-cell D

        # Change 4: no S_cap -- growth sees the substrate concentration
        # directly, with no ceiling on how much of it is "usable".
        f_S = S / (S + Ks)
        logistic = np.where(habitable, 1.0 - B / K_safe, 0.0)
        growth = r_max * f_S * B * logistic
        uptake = (1.0 / Y) * growth

        R_field = (1.0 - Y) * growth + m0 * B
        R_t[t] = R_field.sum()
        running += R_t[t] * dt
        cum_co2[t] = running
        # Stage 10 (soil_respiration_prompt_stage10_final.md SS3): the
        # growth-vs-maintenance split, so a real burst (growth-associated
        # respiration the MAJORITY of R at the peak) can be told apart from
        # pure decay of the seeded biomass (R pinned at m0*B0).
        B_total_t[t] = float(B.sum())
        growth_term_t[t] = float(((1.0 - Y) * growth).sum())
        maintenance_term_t[t] = float((m0 * B).sum())
        S_habitat_mean_t[t] = float(S[habitable].mean()) if habitable.any() else 0.0
        habitats_recruited_t[t] = count_active_habitats(growth, region_labels)
        if debug_growth_distribution:
            growth_distribution_t[t] = growth[macro_mask]
        # Stage 16 SS4c: growth-weighted "active mass" (Stage 14's
        # recommended graded replacement for the boolean recruited count),
        # and per-region biomass for the Portell-style hotspot metric.
        active_mass_t[t] = float(growth[macro_mask].sum())
        if track_region_biomass and n_regions_total > 0:
            region_biomass_t[t] = np.bincount(
                region_labels_flat, weights=B.ravel(), minlength=n_regions_total + 1)[1:]
        # Stage 8 Change 2 (soil_respiration_prompt_stage8_depletion.md SS5):
        # fraction of this soil's (equal) total OM still undissolved -- the
        # direct evidence for whether a soil's fuel pool is actually
        # depleting (Sand should fall fast) or barely touched (Vertisol
        # should fall slowly), not inferred only from R(t)'s shape.
        fuel_remaining_t[t] = (
            float(OM_coating.sum() + OM_trapped.sum()) / om_total_init if om_total_init > 0 else 0.0)

        dB = growth - m0 * B - m_s * (1.0 - f_S) * B

        # Change 1: implicit transport. Ported sandbox pattern -- add the
        # (explicit) OM-release source BEFORE the implicit diffusion solve,
        # apply the (explicit) uptake sink and matrix<->macro exchange terms
        # AFTER it.
        source_increment = dt * (k_dis * OM_coating + k_dis_slow * OM_trapped)
        S_diffused = implicit_diffusion_step(
            S + source_increment, D_total, dx, dt, implicit_iterations, implicit_relaxation,
            tol=implicit_tol)
        S_new = S_diffused - dt * uptake + dt * leak_exchange(S_diffused, porous_matrix_mask, macro_mask, k_leak)
        if wicking_enabled:
            S_new = S_new + dt * wicking_substrate_exchange(
                S_diffused, porous_matrix_mask, macro_mask, wet_mask, k_wick, wick_radius)

        dOM_coating = -k_dis * OM_coating
        dOM_trapped = -k_dis_slow * OM_trapped

        B = np.clip(B + dt * dB, 0.0, None)
        S = np.clip(S_new, 0.0, None)
        OM_coating = np.clip(OM_coating + dt * dOM_coating, 0.0, None)
        OM_trapped = np.clip(OM_trapped + dt * dOM_trapped, 0.0, None)

    grids_out = dict(grids)
    grids_out["water_mask"] = wet_mask
    dry_macro_last = macro_mask & ~wet_mask
    grids_out["D"] = np.where(
        wet_mask, D_macro_base * (1.0 + d_cluster_scale * habitat_scale_last),
        np.where(dry_macro_last, D_film, 0.0))
    grids_out["dry_macro_mask"] = dry_macro_last
    grids_out["dry_macro_cells"] = int(dry_macro_last.sum())
    grids_out["dry_macro_fraction_of_habitat"] = (
        float(dry_macro_last.sum()) / float(macro_mask.sum()) if macro_mask.any() else 0.0)
    grids_out["habitat_labels"] = habitat_labels_last
    grids_out["habitat_scale"] = habitat_scale_last
    grids_out["habitat_count"] = int(habitat_count_t[-1]) if habitat_count_t.size else 0
    grids_out["mean_habitat_count_t"] = float(habitat_count_t.mean()) if habitat_count_t.size else 0.0
    grids_out["OM"] = OM_coating + OM_trapped
    grids_out["OM_coating"] = OM_coating
    grids_out["OM_trapped"] = OM_trapped

    result = {
        "R_t": R_t,
        "cum_co2": cum_co2,
        "S_habitat_mean_t": S_habitat_mean_t,
        "habitat_count_t": habitat_count_t,
        "habitats_recruited_t": habitats_recruited_t,
        "fuel_remaining_t": fuel_remaining_t,
        "B_total_t": B_total_t,
        "growth_term_t": growth_term_t,
        "maintenance_term_t": maintenance_term_t,
        "active_mass_t": active_mass_t,
        "B_final": B,
        "S_final": S,
        "OM_final": OM_coating + OM_trapped,
        "OM_coating_final": OM_coating,
        "OM_trapped_final": OM_trapped,
        "grids": grids_out,
        "dt": dt,
    }
    if track_region_biomass:
        result["region_biomass_t"] = region_biomass_t
    if debug_growth_distribution:
        result["growth_distribution_t"] = growth_distribution_t
        result["macro_mask_debug"] = macro_mask
        result["region_labels_debug"] = region_labels
    return result
