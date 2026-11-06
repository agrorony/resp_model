# Stage 22 — scan-structure geometry and pre-run predictions

Committed before any Stage 22 scan-structure respiration run. This is a separate
scan-structure branch of Stage 22; the earlier equal-wetting Stage 22 artifacts
remain unchanged.

## Geometry

All four inputs passed the registered label and porosity gates. Main crops are
100³ voxels, selected from the 216 non-overlapping grid positions per specimen
with seed 22. No zero-habitat crops were rejected. The four domain-size checks
are 256³ voxels at origin `(197, 197, 197)`. Full per-crop reports, all candidate
origins, specimen/soil means and SDs, and structure files are recorded in
`results/stage22_scan_structure/step1.json`.

| Group | Resolved porosity | Matrix | Solid | Pore LT d10 / d50 / d90 (µm) | Largest pore cluster | Spans x / y / z | S_w | Habitat fraction | Wet fraction of habitat |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| VA | 0.131 | 0.854 | 0.016 | 30.0 / 69.4 / 159.4 | 0.538 | 0.25 / 0.75 / 0.25 | 0.913 | 0.131 | 0.354 |
| VB | 0.197 | 0.785 | 0.018 | 52.6 / 118.4 / 234.9 | 0.882 | 1.00 / 1.00 / 0.75 | 0.830 | 0.197 | 0.175 |
| Vertisol mean | 0.164 | 0.819 | 0.017 | 41.3 / 93.9 / 197.1 | 0.710 | 0.625 / 0.875 / 0.50 | 0.872 | 0.164 | 0.265 |
| SA | 0.247 | 0.082 | 0.671 | 46.9 / 93.8 / 151.9 | 0.958 | 0.75 / 1.00 / 1.00 | 0.334 | 0.247 | 0.111 |
| SB | 0.392 | 0.082 | 0.526 | 52.6 / 109.0 / 178.5 | 0.997 | 1.00 / 1.00 / 1.00 | 0.255 | 0.392 | 0.098 |
| Sand mean | 0.320 | 0.082 | 0.598 | 49.7 / 101.4 / 165.2 | 0.977 | 0.875 / 1.00 / 1.00 | 0.294 | 0.320 | 0.104 |

The Sand fines mapping reached these solid local-thickness thresholds:

| Crop | SA_c0 | SA_c1 | SA_c2 | SA_c3 | SB_c0 | SB_c1 | SB_c2 | SB_c3 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Threshold (µm) | 135.0 | 142.5 | 135.0 | 135.0 | 112.8 | 127.8 | 105.2 | 105.2 |

At 15 µm, the synthetic Vertisol control's median plate thickness is 40 µm,
or about 2.67 voxels. Its geometry and the Sand control's geometry are in the
same Step 1 JSON.

## Pre-registered expectations

- **Expected transport regime:** Sand is expected to be dominated by drained
  pore films plus matrix/pore exchange: only 10.4% of habitat is wet and
  `S_w = 0.294`. Vertisol is also expected to shift toward drained-film
  transport and exchange rather than the Stage 21 wet-pore-diffusion regime:
  its wet habitat fraction is 26.5%, far below Stage 21's 73.8%, although its
  `S_w = 0.872` should make drained films less strongly suppressed than in Sand.
- **Versus Stage 21 at 5 µm:** scan Sand has a similar porosity and matrix
  fraction to synthetic Sand, but a smaller median pore LT (101 versus 125 µm)
  and more wet habitat (10.4% versus 1.2%); expect a higher/earlier Sand
  response. Scan Vertisol has more matrix and less resolved pore volume than
  synthetic Vertisol, with a larger median pore LT (94 versus 35 µm) and less
  wet habitat (26.5% versus 73.8%); expect a lower/later Vertisol response.
  Together these changes are predicted to make the S–V respiration distance
  smaller than Stage 21's 1.23. The peak-time order is predicted to remain
  Sand before Vertisol.
- **Versus the 15 µm synthetic controls:** scan Sand has comparable porosity
  and matrix fraction to S15, but substantially more wet habitat (10.4% versus
  about 1.3%); expect scan Sand to be higher/earlier than S15. Scan Vertisol
  has less porosity and more matrix than V15, and far less wet habitat (26.5%
  versus about 89.7%); expect scan Vertisol to be lower/later than V15.

The period-shape classifier for Step 6 is fixed as follows for period means
`M1..M4`: **rise** if `M4` is a maximum and `M4 > 1.1*M1`; **rise-then-fall**
if `M2` or `M3` is a maximum and `M4 < 0.9*max(M1..M4)`; **fall** if `M1` is
a maximum and `M4 < 0.9*M1`; otherwise **flat**.
