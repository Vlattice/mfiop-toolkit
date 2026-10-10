# ORCA calculation of the vibronic/electronic dipole (μ_eff)

First-principles check of the reference transition dipole μ = 2.1 D used
throughout MFIOP, on the tryptophan chromophore (3-methylindole).

- Method: B3LYP/def2-SVP D3BJ (geometry, frequencies); CAM-B3LYP/def2-SVP (TD-DFT).
- Result: the ¹Lₐ band appears near 245–249 nm with oscillator strength
  f = 0.065 (lowest state) to 0.128 (La/Lb band), i.e. |μ_ge| = 1.86–2.59 D.
  The literature value 2.1 D (Callis 1997) lies between these: the electronic
  reference dipole is confirmed from first principles.
- The ~1020 cm⁻¹ vibrational mode is intrinsically weak (~10⁻³ km/mol),
  consistent with the detectability-budget prediction that the IR window sits
  below the nano-FTIR resolvability threshold unless Herzberg–Teller borrowing
  is exceptional.

Files: 1_opt.out (geometry), 2_freq.out (frequencies), 3_tddft.out (excited states).
