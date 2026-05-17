# MFIOP Toolkit

**Minimal reproducible implementation of the MFIOP geometry → spectrum pipeline.**

This repository contains the core code behind the MFIOP framework (Mesoscopic Field-Infrared Oscillatory Processing). It maps a protein structure (PDB) to a dipolar coupling matrix, diagonalises it, and computes the resulting anisotropic infrared spectral response A(ω) — all directly from atomic coordinates, with no fitted parameters.

This is a teaching and verification toolkit. It reproduces the central geometry → spectrum mapping described in the MFIOP v2.0 paper. It does **not** include Monte Carlo robustness analysis, the ZB composite score, the proteome-scale scan, or the inhomogeneous broadening module. Those will be added in future releases.

For the full theoretical framework, see:

> Zumpano Blumenfeld, V. (2026). *MFIOP — Mesoscopic Field-Infrared Oscillatory Processing: A Geometry-Determined Bioinformatic Framework for Spectral Organization in Protein Dipolar Lattices.* Zenodo. https://doi.org/10.5281/zenodo.18166770

---

## Pipeline overview

```
PDB file
   │
   ▼
[pdb_parser.py]      Extract TRP positions + Cγ→CE2 dipole vectors
   │
   ▼
[dipole_coupling.py] Build J_ij = κ_ij / r_ij³
   │
   ▼
[hamiltonian.py]     Diagonalise J → eigenvalues λ_k, eigenvectors v_k
   │
   ▼
[spectrum.py]        A(ω) = Σ_k |f_k|² / [(ω−ω_k)² + γ²]
   │
   ▼
[run_5SYF.py]        Print mode summary + plot
```

---

## Quickstart

```bash
# 1. Clone
git clone https://github.com/Vlattice/mfiop-toolkit.git
cd mfiop-toolkit

# 2. Install dependencies
pip install -r requirements.txt

# 3. Download the reference PDB
#    Place 5SYF.pdb (or 5SYF.cif) in the data/ folder.
#    Direct download: https://files.rcsb.org/download/5SYF.pdb

# 4. Run the example
python examples/run_5SYF.py
```

Expected output for PDB 5SYF (μ = 2.1 D, γ = 5.31 cm⁻¹):

```
λ_max       = 9.49 cm⁻¹
C_eff       = 1.79
Dominant modes at λ_k: −9.49, −1.44, −0.34 cm⁻¹
Anisotropy A_z/(A_x+A_y) ≈ 0.652
n_dominant (in window) = 3
s_max / γ = 1.52   (> 1 → resolvable by nano-FTIR)
```

---

## Module reference

| File | Responsibility |
|------|----------------|
| `src/pdb_parser.py` | Load TRP positions and Cγ→CE2 unit vectors from PDB/CIF |
| `src/dipole_coupling.py` | Build symmetric J matrix (κ/r³, dimensionless) |
| `src/hamiltonian.py` | Diagonalise J; summarise modes with \|f_k\|² |
| `src/spectrum.py` | Compute A(ω) along x, y, z; canonical anisotropy |
| `src/utils.py` | Physical constants; convert dimensionless J to cm⁻¹ |
| `examples/run_5SYF.py` | End-to-end example with plot |

---

## Physical notes

**Dipole axis.** Cγ→CE2 is the standard approximation for the ¹Lₐ transition polarisation of tryptophan (Callis, 1997).

**Coupling matrix.** `dipole_coupling.py` returns dimensionless J (κ/r³ in Å⁻³). Use `utils.dimensionless_to_cm1(J)` to convert to cm⁻¹ using μ = 2.1 D as the reference electronic transition dipole.

**|f_k|².** Unnormalised squared projection of mode k onto the laboratory axis. Values greater than 1 are physically expected — these are oscillator strength projections, not probabilities.

**γ = 5.31 cm⁻¹.** Lower bound for vibrational dephasing in proteins (γ = 10¹² rad/s, Mukherjee et al. 2006). The near-degenerate pair (Δ ≈ 1.1 cm⁻¹) is unresolvable at room temperature; the dominant shoulder (Δ ≈ 8.05 cm⁻¹) is resolvable.

---

## Reproducibility

All analyses run on Python 3.10+. Pinned versions in `requirements.txt`. The example script uses no random number generation — output is deterministic for a given PDB file.

---

## Citation

If you use this code, please cite:

```bibtex
@software{ZumpanoBlumenfeld_MFIOP_2026,
  author       = {Zumpano Blumenfeld, Verónica},
  title        = {MFIOP — Mesoscopic Field-Infrared Oscillatory Processing},
  year         = {2026},
  publisher    = {Zenodo},
  version      = {v2.0},
  doi          = {10.5281/zenodo.18166770},
  url          = {https://doi.org/10.5281/zenodo.18166770}
}
```

---

## License

MIT — see [LICENSE](LICENSE) for details.

---

## Contact

Verónica Zumpano Blumenfeld
Montevideo, Uruguay
ORCID: [0009-0006-2030-1849](https://orcid.org/0009-0006-2030-1849)

Issues and pull requests welcome.
