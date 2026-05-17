"""
hamiltonian.py
--------------
Diagonalize the dipolar coupling matrix J and return collective eigenmodes.

The system is described by a Frenkel-exciton (long-range XY) Hamiltonian:

    H = Σ_{i<j} J_ij (σ_i⁺ σ_j⁻ + h.c.) + Σ_i (ħΩ_i/2) σ_iz + H_diss

In the single-excitation subspace and with identical site energies (Ω_i = ω₀),
the eigenvalue problem reduces to diagonalizing J directly.

Eigenvalues λ_k give the collective mode splittings relative to ω₀:
    ω_k = ω₀ + λ_k
"""

import numpy as np


def diagonalize(J: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Diagonalize the coupling matrix J.

    Uses numpy.linalg.eigh (symmetric/Hermitian eigensolver) for
    numerical stability and guaranteed real eigenvalues.

    Parameters
    ----------
    J : np.ndarray, shape (N, N)
        Symmetric coupling matrix.

    Returns
    -------
    eigenvalues : np.ndarray, shape (N,)
        Collective mode splittings λ_k (sorted ascending).
    eigenvectors : np.ndarray, shape (N, N)
        Columns are eigenvectors v_k; eigenvectors[:, k] is mode k.
    """
    eigenvalues, eigenvectors = np.linalg.eigh(J)
    return eigenvalues, eigenvectors


def mode_summary(
    eigenvalues:  np.ndarray,
    eigenvectors: np.ndarray,
    dipoles:      np.ndarray,
    omega_0:      float = 1020.0,
    fk2_threshold: float = 0.5,
) -> list[dict]:
    """
    Summarise collective modes: frequency, oscillator strength projection, classification.

    Parameters
    ----------
    eigenvalues : np.ndarray
        Eigenvalues from diagonalize().
    eigenvectors : np.ndarray
        Eigenvectors from diagonalize().
    dipoles : np.ndarray, shape (N, 3)
        Unit dipole vectors; column 2 is the z (laboratory) axis.
    omega_0 : float
        Carrier frequency in cm⁻¹ (default: 1020 cm⁻¹ for TRP ¹Lₐ).
    fk2_threshold : float
        |f_k|² threshold to classify a mode as dominant.

    Returns
    -------
    modes : list of dict
        One entry per mode with keys: k, lambda_k, omega_k, fk2_z, dominant.
    """
    # z-axis projection: f_k = Σ_i (ê_i · ẑ) * v_ki
    z_proj = dipoles[:, 2]                        # ê_i · ẑ for each residue
    fk2    = (eigenvectors.T @ z_proj) ** 2       # |f_k|² for each mode

    modes = []
    for k in range(len(eigenvalues)):
        modes.append({
            "k":         k,
            "lambda_k":  float(eigenvalues[k]),
            "omega_k":   float(omega_0 + eigenvalues[k]),
            "fk2_z":     float(fk2[k]),
            "dominant":  bool(fk2[k] > fk2_threshold),
        })
    return modes
