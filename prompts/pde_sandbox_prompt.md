# Claude Code task: PDE sandbox — search for a transport formulation with sensitive output

Build a SELF-CONTAINED numpy-only sandbox in a NEW folder `pde_sandbox/`. Do not
import, modify, or depend on any existing project code. This is a pure
mathematical problem: a reaction–diffusion system on a heterogeneous domain whose
scalar output functional currently degenerates (becomes insensitive to its
coefficient fields). Find, by a budget-capped search over the formulation, a
version whose output is sensitive and produces two distinct, robust temporal
signatures. **No physical, biological, or material interpretation anywhere in the
code, comments, or report — this is abstract math only.**

## 0. Budget

`MAX_ITERATIONS = 40` (cheap 2D runs; a search, not a single build). Stop when the
success criterion (§5) holds under reseeding, or the budget is spent — then report
the best formulation found.

## 1. Fields and evolution (the minimal system)

Domain Ω = [0,L]² on an N×N grid, Neumann (zero-flux) boundaries, t ∈ [0,T].
Three nonnegative scalar fields u, v, w:

```
∂ₜ w = −λ w
∂ₜ u = 𝓛_D[u] + λ w − φ(u) v
∂ₜ v = κ(x) · ( φ(u) v (1 − v) − μ v )

φ(u) = u / (u + a₀)          # saturating, φ(0)=0
```

Initial data: u(·,0)=0, w(·,0)=w₀(x), v(·,0)=v₀·𝟙_{Ω_r}.
`𝓛_D` is the transport operator (a search choice, §4). λ, a₀, μ, v₀ are constants.

**Observable (the only output that matters):**
```
R(t) = ∫_Ω φ(u) v dx        (discrete: sum over cells)
```

## 2. Coefficient fields from a structure field

A scalar structure field s(x) > 0 is drawn from a spatially-correlated random
field (correlation length ℓ_c). Three fixed maps produce the coefficients:

```
D(x)  = 𝒟(s)   ≥ 0          transport coefficient; ALLOWED to be 0 on a subset
κ(x)  = 𝒦(s)   ≥ 0          reaction weight, support Ω_r = { s ∈ [a,b] }
w₀(x) = 𝒲(s)   ≥ 0          initial reservoir, support ⊂ { s < a }   (DISJOINT from Ω_r)
```

The essential structure — keep it in every formulation: **supp(w₀) and Ω_r are
disjoint, coupled only through D(x).** This separation is the whole problem.

## 3. Two instances (this stage: TWO, not three)

Two fixed structure-field laws generate two coefficient triples via the SAME maps
and SAME dynamics:

- **Instance ↑** — a broad s-distribution: Ω_r appears as many small patches
  scattered at a wide range of distances from a large, spread-out reservoir
  support (large total ∫w₀). Target: **R↑(t) increasing** on [0,T].
- **Instance ↓** — a narrow s-distribution: Ω_r is a compact region adjacent to a
  small reservoir support (small ∫w₀). Target: **R↓(t) decreasing** on [0,T]
  (after any brief initial transient).

The two instances differ ONLY in the structure-field law P (and seed). Everything
else — maps 𝒟,𝒦,𝒲, constants, operator, scheme, T, L, N — is SHARED and is what
the search varies. Success requires the ONE shared formulation to turn structure ↑
into an increasing output and structure ↓ into a decreasing output. Do not give
the two instances different dynamics or different constants — that would defeat the
point (structure must cause the signature).

## 4. What the search may vary (the formulation)

Anything below, shared across both instances:

- **Time-integration scheme for u** — explicit Euler (stability cap
  `max D·dt/dx² < 0.2`, which forces the diffusion length small) vs an **implicit /
  Crank–Nicolson** solve (unconditionally stable, lets the effective diffusion
  length span the domain). *Try implicit first — it is the leading hypothesis.*
- **Transport operator 𝓛_D** — standard heterogeneous divergence-form
  `∇·(D∇u)` (harmonic-mean faces) vs a **nonlocal / graph-Laplacian** relaxation
  that equilibrates u quickly within each connected {D>0} component (so the tempo
  is set by λ, not slow diffusion).
- Magnitudes and maps: the scale of 𝒟, the constants λ, a₀, μ, v₀.
- The horizon T and domain scale L/N (i.e. the separation Δ measured in cells).
- The maps' thresholds a, b and the reservoir scale in 𝒲.

## 5. Success criterion ("works")

Normalize R̂ = R / max R. Split [0,T] in thirds; e, l = mean of first, last third.
Over an ensemble of ≥ 5 seeds per instance, require robustly:

- **Correct signatures:** median(l−e) > +δ for Instance ↑, and median(e−l) > +δ for
  Instance ↓ (start δ = 0.15·peak).
- **Distinct:** ‖R̂↑ − R̂↓‖ (1 − Pearson corr, plus a normalized L2 term) > τ
  (start τ = 0.3), at the median seed.
- **Nontrivial:** ∫₀ᵀ ∫_Ω φ(u) v dx dt > θ · ∫₀ᵀ R dt with θ near 1 by construction,
  AND the output is NOT reproducible by the source-free baseline — i.e. run each
  instance once with w₀ ≡ 0 and confirm R collapses to ~0; the real run must be
  well above it. This rules out the degenerate "v just decays" solution.

Report the score `Score = 𝟙[both signatures] · (‖R̂↑−R̂↓‖) · 𝟙[both nontrivial]`.

## 6. Hard constraints

u, v, w ≥ 0 preserved every step (clip); the u-transport conserves mass up to the
source (λw), sink (φ(u)v), and boundary flux (verify the pure-transport step sums
to the boundary flux, ~0 for Neumann); the scheme is stable at the chosen (dt,dx);
fully deterministic given seed. numpy only (implement any implicit solve with a
numpy banded/sparse-free iterative method — e.g. Jacobi/Gauss–Seidel/conjugate
gradient hand-rolled — or an ADX/operator-splitting scheme; do NOT add scipy).

## 7. The loop

```
for iter in 1..MAX_ITERATIONS:
    pick/mutate a shared formulation (§4)
    build coeffs for Instance ↑ and Instance ↓ from their fixed structure laws
    integrate both to T; compute R↑(t), R↓(t) and the w₀≡0 baselines
    evaluate §5; append a LOGBOOK entry (formulation, signatures, distance, nontriv, score)
    if success holds over the seed ensemble: break
keep the best formulation; write it out explicitly
```

Each LOGBOOK entry records the exact formulation tried and why it was mutated —
enough to reconstruct the search. Do not tune the two structure laws toward the
target; only the shared formulation is searched.

## 8. Deliverables (all inside `pde_sandbox/`)

- `sandbox.py` — fields, coefficient construction, the selectable operator/scheme,
  the integrator, the observable.
- `search.py` — the §7 loop + §5 scoring, writing `LOGBOOK.md`.
- `results/` — R↑(t) and R↓(t) for the winning formulation (with the w₀≡0
  baselines overlaid), the s / D / κ / w₀ field maps for both instances, and a
  score-vs-iteration trace.
- `LOGBOOK.md` — every iteration.
- `RESULTS.md` — the winning formulation written out as explicit equations +
  constants, the diffusion-length-vs-separation numbers (ℓ(T) vs Δ) before/after,
  why it works, and any honest failure to fully meet §5.
- `README.md` — how to run.

## 9. Order of work

1. Build the minimal system (§1–2) with explicit Euler + `∇·(D∇u)`; reproduce the
   degeneracy (R↑ ≈ R↓, both ≈ the w₀≡0 baseline) and record ℓ(T) vs Δ.
2. Add the implicit/Crank–Nicolson u-solver and/or the nonlocal operator; rerun.
3. Run the §7 search over the remaining levers; converge on a formulation meeting §5.
4. Checkpoint: STOP; write `RESULTS.md` with the explicit winning formulation.

Keep it purely mathematical throughout. The single goal: a shared transport
formulation under which the disjoint source→reaction coupling makes the output R(t)
sensitive to structure, yielding one increasing and one decreasing signature,
robustly and nontrivially.
