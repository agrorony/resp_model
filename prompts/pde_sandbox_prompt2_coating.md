# Claude Code task: PDE sandbox, continue — lock local-implicit transport, engineer the sharp-rise-then-fall via a coating/trapped source

Continue working IN the existing `pde_sandbox/` (do not restart; keep `sandbox.py`,
`search.py`, `LOGBOOK.md`, the harness and the good engineering already there —
the fixed reaction logistic, the import fix, the honest fallback/RESULTS, the
seed-ensemble scoring, the full-budget run). Purely mathematical system as before;
`MAX_ITERATIONS = 40`.

## 0. Where we are

Two earlier passes exist. The FIRST used **local implicit backward-Euler diffusion**
(M-matrix ⇒ nonnegativity at any dt, harmonic-mean faces as true barriers) and got
BOTH signatures robust. The SECOND regressed by letting the search adopt a
**nonlocal within-component** operator (`du/dt = ρ(mean_C(u) − u)`): it aces the
broad instance but leaves the compact instance **disconnected** (`separation = inf`,
`down_trend = 0.0`, `down nontrivial = False`, strict score 0). The saved
formulation is currently the worse nonlocal one — recover the local-implicit
configuration from `LOGBOOK.md` if needed.

This stage merges the two: **the first agent's local-implicit transport** as the
locked operator, **the second agent's harness/discipline**, and a source-placement
fix that makes the compact instance produce a clean **sharp rise then fall**.

## 1. Lock the transport operator

Make **local implicit backward-Euler diffusion with harmonic-mean face coefficients**
the DEFAULT and primary operator (it bridges a finite gap with a real, graded delay
— exactly what a rise-then-fall needs). The nonlocal within-component operator may
remain in the code only as an optional comparison; it must NOT be the operator the
search ships, and the search must never select a formulation in which either
instance is uncoupled (see §4).

## 2. The focus target: sharp rise then fall (the compact / "down" instance)

This stage's priority is getting the compact instance to a robust **sharp-rise-then-
fall** signature: R(t) climbs quickly, peaks in the FIRST third, then declines
monotonically; nontrivial; stable across ≥5 seeds. Keep the broad / "up" instance's
increasing signature (already robust). Success = both signatures present, distinct
(‖R̂↑ − R̂↓‖ > τ), both nontrivial, robust over the seed ensemble.

## 3. Source construction — coating vs trapped (the key change)

Partition the compact instance's source field w₀ into two placements relative to
the transport network:

- **Coating (majority fraction, ~0.7–0.85 of the source mass):** a thin layer on
  the interface between the connected transport-accessible region ({D > 0}) and the
  reaction region Ω_r — i.e. immediately adjacent to Ω_r along a **finite, connected**
  path. This guarantees `separation < ∞` for the compact instance (fixing the
  Stage-2 `separation = inf` failure) and gives fast early delivery → the sharp rise.
- **Trapped (minority fraction, ~0.15–0.3):** placed in interior / low-D pockets
  (behind {D = 0} barriers), effectively undeliverable within T.

Because the coating source is finite and depletes (∂ₜw = −λw plus consumption), the
compact instance blooms fast off the coating and then declines as the coating runs
out and only trapped source remains — the rise-then-fall. Expose the coating
fraction and the coating thickness as search parameters.

*(Physical correspondence, for context only — keep the code abstract: this mirrors
the research hypothesis that in a sandy medium most organic matter sits as a coating
on the matrix surface with a smaller trapped fraction. Use it to motivate the
placement, not to add any material/biological terms to the math.)*

For the broad / "up" instance, keep its distributed source at many staggered
distances (progressive recruitment → increasing) — do not force the coating there.

## 4. Fix the search objective

Penalize triviality DURING mutation, not only at final scoring: any candidate in
which EITHER instance is uncoupled (`nontrivial = False` or `separation = inf`) gets
a hard-zero / large-negative merit, so the search cannot "win" by acing one instance
while the other dies. Keep the mean-over-≥5-seeds scoring so fragile configs are
penalized. The shipped formulation must satisfy both instances or be reported as an
honest partial result.

## 5. Keep intact

The abstract (u, v, w) system, φ(u)=u/(u+a₀), the observable R=∫φ(u)v, Neumann
boundaries, nonnegativity, determinism-by-seed, numpy-only, and the w₀≡0 baseline
check (real run must sit well above the source-free baseline). Do not tune the two
structure laws toward the target beyond the coating/trapped placement rule; the
shared formulation is what the search varies.

## 6. Deliverables (update in place)

- `sandbox.py`: local-implicit BE as default operator; the coating/trapped source
  construction for the compact instance (with coating fraction + thickness params).
- `search.py`: the §4 triviality penalty in the objective.
- `results/`: R↑(t) and R↓(t) for the winning formulation with w₀≡0 baselines
  overlaid; the compact instance's coating/trapped source map and its D/κ maps;
  score-vs-iteration trace.
- `RESULTS.md`: the explicit winning formulation + constants, the compact instance's
  `separation` (must be finite) and peak-location (first third), ℓ(T) vs Δ, and an
  honest statement of what met/failed the success gate.
- `LOGBOOK.md`: every iteration.

## 7. Order of work

1. Restore local-implicit BE diffusion as the default; confirm it reproduces the
   earlier both-signatures-coupled behavior (not the nonlocal regression).
2. Implement the coating/trapped source split for the compact instance; verify its
   `separation` is now finite and it produces a rise-then-fall on a single seed.
3. Add the §4 triviality penalty to the objective.
4. Run the 40-iteration budget focused on making the compact rise-then-fall robust
   while keeping the broad instance increasing; score over the seed ensemble.
5. Checkpoint: STOP; regenerate `RESULTS.md` honestly (success or best partial).

The single goal of this stage: with the locked local-implicit operator and a
coating-dominated compact source, find a shared formulation where the compact
instance robustly shows a **sharp rise then fall** and the broad instance stays
increasing — two distinct, nontrivial, robust signatures.
