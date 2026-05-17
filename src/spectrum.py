"""
spectrum.py
-----------
Compute the anisotropic infrared spectral response A(ω) from collective eigenmodes.

    A(ω) = Σ_k |f_k|² / [(ω − ω_k)² + γ²]

where:
    f_k = Σ_i (ê_i · ê_axis) * v_ki   (oscillator strength projection)
    γ   = dephasing linewidth in cm⁻¹   (default: 5.31 cm⁻¹ ≈ 10¹² rad/s)

|f_k|² values are unnormalized projections and may exceed 1 — they are not
probabilities but squared amplitudes along the chosen laboratory axis.

The anisotropy ratio A_z/(A_x+A_y) is evaluated at the peak of A_z
(canonical definition used in the MFIOP paper).
"""

import numpy as np


_DEFAULT_GAMMA   = 5.31    # cm⁻¹  (γ = 10¹² rad/s)
_DEFAULT_OMEGA_0 = 1020.0  # cm⁻¹  (TRP ¹Lₐ carrier frequency)


def compute_spectrum(
    eigenvalues:  np.ndarray,
    eigenvectors: np.ndarray,
    dipoles:      np.ndarray,
    omega_grid:   np.ndarray,
    gamma:        float = _DEFAULT_GAMMA,
    omega_0:      float = _DEFAULT_OMEGA_0,
) -> dict:
    """
    Compute A(ω) along x, y, z laboratory axes.

    Parameters
    ----------
    eigenvalues : np.ndarray, shape (N,)
    eigenvectors : np.ndarray, shape (N, N)
        Columns are eigenvectors.
    dipoles : np.ndarray, shape (N, 3)
        Unit dipole vectors in Cartesian [x, y, z] order.
    omega_grid : np.ndarray, shape (M,)
        Frequency grid in cm⁻¹.
    gamma : float
        Lorentzian half-width (cm⁻¹).
    omega_0 : float
        Carrier frequency (cm⁻¹).

    Returns
    -------
    result : dict with keys:
        'omega'      : frequency grid
        'A_z'        : A(ω) along z
        'A_x'        : A(ω) along x
        'A_y'        : A(ω) along y
        'omega_k'    : collective mode frequencies
        'fk2'        : dict of |f_k|² per axis ('x', 'y', 'z')
        'anisotropy' : A_z / (A_x + A_y) evaluated at peak of A_z
                       (canonical definition)
    """
    omega_k = omega_0 + eigenvalues                 # absolute mode frequencies

    # Lorentzian denominator — broadcast: (M, N)
    denom = (omega_grid[:, None] - omega_k[None, :]) ** 2 + gamma ** 2

    A     = {}
    fk2   = {}

    for axis_idx, axis_name in enumerate(["x", "y", "z"]):
        proj          = dipoles[:, axis_idx]        # (N,)
        f_k           = eigenvectors.T @ proj       # (N,) dot product per mode
        fk2[axis_name] = f_k ** 2                  # unnormalized |f_k|²

        # A(ω) = Σ_k |f_k|² / denom_k(ω)
        A[axis_name] = (fk2[axis_name][None, :] / denom).sum(axis=1)

    # Canonical anisotropy: evaluate A_x+A_y at the index of A_z peak
    peak_idx  = int(np.argmax(A["z"]))
    A_z_peak  = float(A["z"][peak_idx])
    A_xy_peak = float(A["x"][peak_idx] + A["y"][peak_idx])
    anisotropy = A_z_peak / (A_xy_peak + 1e-12)

    return {
        "omega":      omega_grid,
        "A_z":        A["z"],
        "A_x":        A["x"],
        "A_y":        A["y"],
        "omega_k":    omega_k,
        "fk2":        fk2,
        "anisotropy": round(anisotropy, 4),
    }


def make_omega_grid(
    center: float = _DEFAULT_OMEGA_0,
    half_width: float = 35.0,
    n_points: int = 4000,
) -> np.ndarray:
    """Return a frequency grid centred on `center` ± half_width cm⁻¹."""
    return np.linspace(center - half_width, center + half_width, n_points)
