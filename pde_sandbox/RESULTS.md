# PDE Sandbox Results

Shared winning formulation:

dw/dt = -lambda w
du/dt = L_D[u] + lambda w - phi(u) v
dv/dt = kappa(x) * (phi(u) v (1 - v) - mu v)
phi(u) = u / (u + a0)

Transport realization: u^{n+1} solves u^{n+1} - dt div(D grad u^{n+1}) = rhs with Jacobi iteration.

Shared constants:
- N = 32, L = 1.0, T = 14.0, steps = 160
- lambda = 0.1400, a0 = 0.1000, mu = 0.2200, v0 = 0.7800
- D map = linear, D_scale = 0.0180, D_floor = 0.0000, D_cut = 0.1800, D_power = 1.3000
- thresholds a = 0.2800, b = 0.4800, kappa_scale = 1.1500
- w0_scale = 0.8500, w0_power = 1.1000
- compact coating fraction = 0.8000, compact coating thickness = 2
- nonlocal_rate = 2.8000, implicit_iterations = 90

Reach versus separation:
- Baseline local explicit reach target is represented by the first logbook entry.
- Winning up instance: reach = 0.5020, separation = 0.1127
- Winning down instance: reach = 0.5020, separation = 0.0442
- Winning down peak location (fraction of horizon) = 0.3125
- Winning down coating fraction = 0.8000

Why this works:
- The shared coefficient maps give the broad structure law a larger source mass and a more fragmented active set, so the observable rises as released mass keeps reaching many active patches.
- The narrow structure law gets a smaller nearby source support, so the observable forms an early pulse and then decays once the limited supply is exhausted.
- The winning transport increases sensitivity by making the source-to-reaction coupling depend on connectivity rather than only on a short local diffusion length.

Outcome summary:
- merit = 1.1719
- up trend = 0.2573
- down trend = 0.1714
- median distance = 0.5863
- up nontrivial = True
- down nontrivial = True
- down shape pass = True
- up baseline integral = 0.0000e+00
- down baseline integral = 0.0000e+00
- transport drift medians = up 5.7244e-13, down 3.1545e-12

Status:
- Met the signature, distance, and nontrivial checks.

Artifacts:
- results/winning_timeseries.csv contains the representative R(t) curves with the source-free baselines overlaid as columns.
- results/up_*.csv and results/down_*.csv contain the structure and coefficient maps for the representative seeds.
- results/down_w0_coating.* and results/down_w0_trapped.* isolate the compact source split.
- results/score_trace.csv records the search trajectory.
