# Stage 24 post-hoc growth-respiration correction

Stage 24 was present on the checked-out `agents/stage24-matrix-sensitivity` base. Its 184 variant traces were corrected post-hoc from their stored `growth_R` and `maint_R`; no Stage 24 simulation was rerun and the Stage 24 branch was not merged.
The base variant was validated against the corrected Stage 23 analysis before comparing variants.

| Variant | Output | Legacy value | Corrected value | Same/different from base, legacy | Same/different from base, corrected |
|---|---|---|---|---|---|
| base | Arm A Sand–Vertisol distinctness | True | True | same | same |
| base | Arm A peak-height order | S / L / V | S / L / V | same | same |
| base | Arm A peak-time order | S / L / V | S / L / V | same | same |
| base | Arm B primary shapes (S/V/L) | rise-then-fall / rise / rise-then-fall | rise-then-fall / rise / rise-then-fall | same | same |
| base | Arm B primary cumulative-R/mass rank | V / L / S | V / L / S | same | same |
| S04 | Arm A Sand–Vertisol distinctness | True | True | same | same |
| S04 | Arm A peak-height order | S / L / V | S / L / V | same | same |
| S04 | Arm A peak-time order | S / L / V | S / L / V | same | same |
| S04 | Arm B primary shapes (S/V/L) | flat / rise / rise-then-fall | rise-then-fall / rise / rise-then-fall | different | same |
| S04 | Arm B primary cumulative-R/mass rank | V / L / S | V / L / S | same | same |
| S16 | Arm A Sand–Vertisol distinctness | False | True | different | same |
| S16 | Arm A peak-height order | S / L / V | S / L / V | same | same |
| S16 | Arm A peak-time order | S / L / V | S / L / V | same | same |
| S16 | Arm B primary shapes (S/V/L) | rise-then-fall / rise / rise-then-fall | rise-then-fall / rise / rise-then-fall | same | same |
| S16 | Arm B primary cumulative-R/mass rank | V / S / L | V / S / L | different | different |
| S30 | Arm A Sand–Vertisol distinctness | False | False | different | different |
| S30 | Arm A peak-height order | S / L / V | S / L / V | same | same |
| S30 | Arm A peak-time order | S / L / V | S / L / V | same | same |
| S30 | Arm B primary shapes (S/V/L) | rise / rise / rise-then-fall | rise-then-fall / rise / rise-then-fall | different | same |
| S30 | Arm B primary cumulative-R/mass rank | V / L / S | V / L / S | same | same |
| VL25 | Arm A Sand–Vertisol distinctness | True | True | same | same |
| VL25 | Arm A peak-height order | S / L / V | S / L / V | same | same |
| VL25 | Arm A peak-time order | S / L / V | S / L / V | same | same |
| VL25 | Arm B primary shapes (S/V/L) | rise-then-fall / rise / rise-then-fall | rise-then-fall / rise / rise-then-fall | same | same |
| VL25 | Arm B primary cumulative-R/mass rank | V / L / S | V / L / S | same | same |
| VL50 | Arm A Sand–Vertisol distinctness | True | True | same | same |
| VL50 | Arm A peak-height order | L / S / V | L / S / V | different | different |
| VL50 | Arm A peak-time order | S / L / V | S / L / V | same | same |
| VL50 | Arm B primary shapes (S/V/L) | rise-then-fall / rise / rise-then-fall | rise-then-fall / rise-then-fall / rise-then-fall | same | different |
| VL50 | Arm B primary cumulative-R/mass rank | V / L / S | V / L / S | same | same |
| CONV | Arm A Sand–Vertisol distinctness | False | False | different | different |
| CONV | Arm A peak-height order | L / S / V | L / S / V | different | different |
| CONV | Arm A peak-time order | S / L / V | S / L / V | same | same |
| CONV | Arm B primary shapes (S/V/L) | rise / rise / rise-then-fall | rise-then-fall / rise-then-fall / rise-then-fall | different | different |
| CONV | Arm B primary cumulative-R/mass rank | V / L / S | V / L / S | same | same |

## What flipped

- S04: Arm B primary shapes flat / rise / rise-then-fall → rise-then-fall / rise / rise-then-fall.
- S04: base comparison changed for Arm B primary shapes (S/V/L).
- S16: Sand–Vertisol distinctness False → True.
- S16: base comparison changed for Arm A Sand–Vertisol distinctness.
- S30: Arm B primary shapes rise / rise / rise-then-fall → rise-then-fall / rise / rise-then-fall.
- S30: base comparison changed for Arm B primary shapes (S/V/L).
- VL50: Arm B primary shapes rise-then-fall / rise / rise-then-fall → rise-then-fall / rise-then-fall / rise-then-fall.
- VL50: base comparison changed for Arm B primary shapes (S/V/L).
- CONV: Arm B primary shapes rise / rise / rise-then-fall → rise-then-fall / rise-then-fall / rise-then-fall.

## Notable criterion change

S16 Arm A Sand–Vertisol distinctness changed False → True. The 0–1500 soil-mean distance moved from 0.23075 to 0.337067 (criterion `> 0.3`); maximum within distance was 0.0748676 → 0.0720477, and minimum cross distance was 0.162959 → 0.226751.

The full corrected variants table and figures are in [`results/stage24_rfix/`](../../results/stage24_rfix/).
