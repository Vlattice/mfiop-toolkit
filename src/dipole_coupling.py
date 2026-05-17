"""
dipole_coupling.py
------------------
Build the dipolar coupling matrix J from TRP positions and dipole unit vectors.

The coupling between residues i and j is:

    J_ij = κ_ij / r_ij³

where the orientation factor is:

    κ_ij = ê_i · ê_j − 3 (ê_i · r̂_ij)(ê_j · r̂_ij)

Note: the physical prefactor μ²/(4πε₀ħ) is omitted here; eigenvalues are
therefore in units of [Å³]⁻¹. Scale by the prefactor when comparing to
experiment (see utils.py for the conversion).

Pairs beyond r_cutoff Å are ignored (contribute < 2% at typical geometries).
"""

import numpy as np
from itertools import combinations


def build_J(
    positions: np.ndarray,
    dipoles:   np.ndarray,
    r_cutoff:  float = 80.0,
) -> np.ndarray:
    """
    Construct the symmetric dipolar coupling matrix J.

    Parameters
    ----------
    positions : np.ndarray, shape (N, 3)
        Residue positions in Å.
    dipoles : np.ndarray, shape (N, 3)
        Unit dipole vectors (CG → CE2).
    r_cutoff : float
        Maximum inter-residue distance to include (Å).

    Returns
    -------
    J : np.ndarray, shape (N, N)
        Symmetric coupling matrix with zero diagonal.
    """
    N = len(positions)
    J = np.zeros((N, N), dtype=float)

    for i, j in combinations(range(N), 2):
        r_vec = positions[j] - positions[i]
        r     = np.linalg.norm(r_vec)

        if r > r_cutoff or r < 1e-6:
            continue

        r_hat = r_vec / r
        ei    = dipoles[i]
        ej    = dipoles[j]

        kappa = np.dot(ei, ej) - 3.0 * np.dot(ei, r_hat) * np.dot(ej, r_hat)
        Jij   = kappa / r**3

        J[i, j] = Jij
        J[j, i] = Jij  # enforce symmetry

    return J


def coupling_summary(positions: np.ndarray, r_cutoff: float = 30.0) -> list[dict]:
    """
    Return a sorted list of inter-residue distances below r_cutoff.

    Useful for identifying the minimal coupling core.
    """
    pairs = []
    for i, j in combinations(range(len(positions)), 2):
        r = float(np.linalg.norm(positions[j] - positions[i]))
        if r <= r_cutoff:
            pairs.append({"i": i, "j": j, "r_Ang": round(r, 3)})
    return sorted(pairs, key=lambda x: x["r_Ang"])
