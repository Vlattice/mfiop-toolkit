# Vibronic effective dipole of the tryptophan 1La transition (Herzberg–Teller)

**Result:** the vibronic transition dipole that governs the ~1020 cm⁻¹ infrared
couplings is **μ_vib ≈ 0.15 D** (range 0.11–0.16 D), computed by TD-DFT
Herzberg–Teller finite difference. This is far below the critical dipole
μ_crit ≈ 1.35 D, giving a detectability D ≈ 0.01: **the predicted infrared
collective signal is not resolvable for this chromophore model.** The geometric
eigenstructure is unaffected (it is independent of μ_eff).

## Why this calculation
The main MFIOP reference numbers use the *electronic* ¹La dipole (|μ_ge| ≈ 2.1 D).
But the quantity that actually drives the vibrational (~1020 cm⁻¹) couplings is the
*vibronic* dipole, set by Herzberg–Teller intensity borrowing from ¹La. This package
computes it directly, closing with a number what earlier versions left as an assumption.

## Method
- Chromophore: 3-methylindole (skatole), neutral singlet.
- Level: CAM-B3LYP / def2-SVP / D3BJ (same as the electronic reference).
- Step 1 — ground-state optimization + frequencies (`inputs/skatole_optfreq.inp`),
  giving normal modes (`skatole_optfreq.hess`). The two **in-plane** ring modes nearest
  1020 cm⁻¹ are at 999.5 and 1048.8 cm⁻¹. (The mode at 1026 cm⁻¹ is out-of-plane and
  does not couple to the in-plane ¹La.)
- Step 2 — each in-plane mode displaced ±0.1 Å along its mass-weighted normal coordinate
  (`inputs/skatole_td_m26_±.inp`, `skatole_td_m28_±.inp`), plus the equilibrium
  (`skatole_td_eq.inp`); ¹La transition dipole recomputed by TD-DFT at each geometry.
- Step 3 — μ_vib = |∂μ_ge/∂Q| · √(ħ/2ω), tracking ¹La by transition-dipole direction
  (phase-aligned). Robust across both in-plane modes and across the ¹La/¹Lb assignment.

## Files
- `inputs/` — the 6 ORCA input files (run in SERIAL, nprocs 1).
- `skatole_optfreq.hess` — Hessian / normal modes from step 1.
- `manifest_round2.json` — displacement amplitudes (ΔQ) and mode frequencies.
- `gen_displacements_step1_parse_hess.py` — parse the .hess, classify in/out-of-plane modes.
- `gen_displacements_step2_make_inputs.py` — generate the ±δ displaced TD-DFT inputs.
- `analyze_vibronic.py` — parse the TD-DFT outputs, compute dμ/dQ and μ_vib.
- `vibronic_results.json` — extracted ¹La dipoles per geometry + final μ_vib.

## Reproduce
1. Run the 6 inputs in `inputs/` with ORCA (serial).
2. `python analyze_vibronic.py` (point it at the 5 TD-DFT .out files) → μ_vib ≈ 0.15 D.

Electronic ¹La reference from the equilibrium calculation: |μ_ge| ≈ 1.85 D
(¹La/¹Lb pair 1.76–1.85 D), consistent with the literature 2.1 D.
