# Stage 24 predictions (pre-run)

This record is based only on the fixed Stage 23 crops, labels, chemistry and the pre-registered
matrix-split variants. It was written before any Stage 24 respiration runs. All totals remain
fixed within each run; changing the matrix voxel count changes OM loading per matrix voxel.
Unchanged soils in a variant will reuse Stage 23 curves.

## Pre-flight worktree check

`git diff HEAD -- configs/biology.yaml` and `git status --short` were empty. The working file
contained CRLF line endings; after the required CRLF-to-LF normalization its SHA-256 is
`c76d9869dc5847f830ac3caadcce4c68986a9180875f3333504c2482a4b0bf5f`, matching HEAD, main and
origin/main. `git stash push -u -m "pre-stage24-leftovers"` reported **No local changes to save**;
therefore there was no edit to stash or use.

## Geometry (crop mean ± SD; interface is matrix–pore voxel faces per crop voxel)

| Variant | Soil | Matrix fraction | Inert solid fraction | Interface / volume | Mean matrix distance (µm) | OM/matrix voxel A | OM/matrix voxel B |
|---|---|---:|---:|---:|---:|---:|---:|
| base | Sand | 0.0821 ± 0.0000 | 0.5981 ± 0.0866 | 0.0859 | 22.28 | 12.1786 | 7.11657 |
| base | Vertisol | 0.8193 ± 0.0588 | 0.0169 ± 0.0057 | 0.1155 | 75.09 | 1.22669 | 2.01371 |
| base | Loess | 0.7478 ± 0.0731 | 0.0064 ± 0.0077 | 0.2151 | 56.79 | 1.3487 | 1.044 |
| S04 | Sand | 0.0400 ± 0.0000 | 0.6402 ± 0.0866 | 0.0569 | 18.59 | 25 | 14.6087 |
| S04 | Vertisol | 0.8193 ± 0.0588 | 0.0169 ± 0.0057 | 0.1155 | 75.09 | 1.22669 | 2.01371 |
| S04 | Loess | 0.7478 ± 0.0731 | 0.0064 ± 0.0077 | 0.2151 | 56.79 | 1.3487 | 1.044 |
| S16 | Sand | 0.1600 ± 0.0000 | 0.5202 ± 0.0866 | 0.1206 | 26.50 | 6.25 | 3.65218 |
| S16 | Vertisol | 0.8193 ± 0.0588 | 0.0169 ± 0.0057 | 0.1155 | 75.09 | 1.22669 | 2.01371 |
| S16 | Loess | 0.7478 ± 0.0731 | 0.0064 ± 0.0077 | 0.2151 | 56.79 | 1.3487 | 1.044 |
| S30 | Sand | 0.3000 ± 0.0000 | 0.3802 ± 0.0866 | 0.1543 | 32.44 | 3.33333 | 1.94783 |
| S30 | Vertisol | 0.8193 ± 0.0588 | 0.0169 ± 0.0057 | 0.1155 | 75.09 | 1.22669 | 2.01371 |
| S30 | Loess | 0.7478 ± 0.0731 | 0.0064 ± 0.0077 | 0.2151 | 56.79 | 1.3487 | 1.044 |
| VL25 | Sand | 0.0821 ± 0.0000 | 0.5981 ± 0.0866 | 0.0859 | 22.28 | 12.1786 | 7.11657 |
| VL25 | Vertisol | 0.6145 ± 0.0441 | 0.2217 ± 0.0166 | 0.1133 | 59.51 | 1.63559 | 2.68494 |
| VL25 | Loess | 0.5608 ± 0.0548 | 0.1933 ± 0.0185 | 0.2099 | 44.15 | 1.79827 | 1.392 |
| VL50 | Sand | 0.0821 ± 0.0000 | 0.5981 ± 0.0866 | 0.0859 | 22.28 | 12.1786 | 7.11657 |
| VL50 | Vertisol | 0.4097 ± 0.0294 | 0.4266 ± 0.0308 | 0.1066 | 47.68 | 2.45338 | 4.02741 |
| VL50 | Loess | 0.3739 ± 0.0366 | 0.3803 ± 0.0360 | 0.1953 | 34.38 | 2.6974 | 2.088 |
| CONV | Sand | 0.3000 ± 0.0000 | 0.3802 ± 0.0866 | 0.1543 | 32.44 | 3.33333 | 1.94783 |
| CONV | Vertisol | 0.4097 ± 0.0294 | 0.4266 ± 0.0308 | 0.1066 | 47.68 | 2.45338 | 4.02741 |
| CONV | Loess | 0.3739 ± 0.0366 | 0.3803 ± 0.0360 | 0.1953 | 34.38 | 2.6974 | 2.088 |

## Directional predictions

- **Sand S04 → S30:** increasing matrix amount extends the Stage 22 thin-solid-first selection
  into thicker, more pore-distant solids. I predict lower peak height and later peak time as
  substrate is spread over more matrix and less of the additional matrix is directly coupled to
  pore habitat.
- **Vertisol/Loess base → VL25 → VL50:** converting the thickest, most pore-distant matrix first
  leaves less matrix closer to pores and raises OM per remaining matrix voxel. I predict higher
  peak height and earlier peak time as more of the fixed OM pool is concentrated at coupled
  matrix locations.
- **CONV:** the matrix fractions are substantially closer than in base, but pore topology,
  habitat fraction and the spatial geometry of selected matrix remain soil-specific. I predict
  that Sand–Vertisol remains distinct under the Stage 23 Arm A rule, although the mean curves
  may move closer.

These are directional predictions, not fitted parameters. No variant will be changed in response
to the outcomes.
