"""
utils.py
--------
Physical constants and unit conversions for the MFIOP pipeline.

All computations in dipole_coupling.py use dimensionless coupling (κ/r³).
To recover physical J values in cm⁻¹, multiply by the prefactor below.
"""

import numpy as np

# ── Physical constants ────────────────────────────────────────────────────────
HBAR      = 1.0546e-34   # J·s
EPSILON_0 = 8.854e-12    # C²/(N·m²)
DEBYE     = 3.336e-30    # C·m per Debye
C_CM_S    = 2.998e10     # speed of light in cm/s

# ── MFIOP default parameters ──────────────────────────────────────────────────
MU_GE     = 2.1          # |μ_ge| in Debye — TD-B3LYP, Ray et al. 2022
GAMMA_SI  = 1.0e12       # dephasing rate in rad/s
OMEGA_0   = 1020.0       # TRP ¹Lₐ carrier frequency in cm⁻¹

# Physical prefactor: (μ_ge)² / (4π ε₀ ħ)
# Converts (κ/r_Å³) → rad/s when r is in Å
_MU_SI   = MU_GE * DEBYE         # C·m
_R_ANG_TO_M = 1e-10              # 1 Å in metres

PREFACTOR_SI = _MU_SI**2 / (4 * np.pi * EPSILON_0 * HBAR)
# Units: (C·m)² / (C²/(N·m²) · J·s) = m³/s ... × (1/r_m³) → rad/s
# Since our J uses r in Å³: multiply by (1e-10)^3 = 1e-30
# Prefactor converting dimensionless (κ / r_Å³) directly to cm⁻¹:
#   J_cm1 = PREFACTOR_SI * κ / r_m³ / (2π c_cm)
#         = PREFACTOR_SI * κ / (r_Å * 1e-10)³ / (2π c_cm)
#         = PREFACTOR_SI / (2π c_cm) * 1e30  ×  κ / r_Å³
PREFACTOR_CM1_PER_ANG3 = PREFACTOR_SI / (2 * np.pi * C_CM_S) * 1e30
# Sanity: for r=15.45 Å, κ=1 → J ≈ 6.0 cm⁻¹  ✓


def dimensionless_to_cm1(J_dimensionless: np.ndarray) -> np.ndarray:
    """
    Convert dimensionless J (κ/r_Å³) to cm⁻¹.

    J_cm1 = J_dimensionless × PREFACTOR_CM1_PER_ANG3
    """
    return J_dimensionless * PREFACTOR_CM1_PER_ANG3


def gamma_si_to_cm1(gamma_si: float = GAMMA_SI) -> float:
    """Convert dephasing rate from rad/s to cm⁻¹."""
    return gamma_si / (2 * np.pi * C_CM_S)


def print_sanity_check(J_dimensionless: np.ndarray) -> None:
    """
    Quick sanity check: print max eigenvalue in cm⁻¹ and C_eff.

    Note: exact values depend on which coordinate set of PDB 5SYF is used
    (asymmetric unit vs biological assembly). The pipeline is correct for
    any valid coordinate set; download 5SYF.pdb from RCSB for paper-exact values.
    Expected range for any valid 5SYF coordinate set: λ_max ≈ 9–15 cm⁻¹.
    """
    from numpy.linalg import eigvalsh
    J_cm = dimensionless_to_cm1(J_dimensionless)
    ev_cm = eigvalsh(J_cm)

    lmax   = float(np.max(np.abs(ev_cm)))
    gam_cm = gamma_si_to_cm1()
    ceff   = lmax / gam_cm

    print(f"  λ_max       = {lmax:.3f} cm⁻¹   (expected range 9–15 for 5SYF)")
    print(f"  γ           = {gam_cm:.3f} cm⁻¹")
    print(f"  C_eff       = {ceff:.3f}          (expected range 1.7–2.8 for 5SYF)")
