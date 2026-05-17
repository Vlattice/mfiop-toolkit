"""
run_5SYF.py
-----------
Minimal example: geometry → spectrum for αβ-tubulin (PDB 5SYF).

Expected output (PDB 5SYF, μ = 2.1 D):
    λ_max  ≈  9.49 cm⁻¹
    C_eff  ≈  1.79
    Dominant modes at splittings: −9.49, −1.44, −0.34 cm⁻¹
    Anisotropy A_z/(A_x+A_y) ≈ 0.652

Usage:
    python examples/run_5SYF.py

Place 5SYF.pdb (or 5SYF.cif) in the data/ folder before running.
Download: https://www.rcsb.org/structure/5SYF
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
import matplotlib.pyplot as plt

from pdb_parser    import load_trp
from dipole_coupling import build_J, coupling_summary
from hamiltonian   import diagonalize, mode_summary
from spectrum      import compute_spectrum, make_omega_grid
from utils         import (dimensionless_to_cm1, gamma_si_to_cm1,
                            print_sanity_check, OMEGA_0)

# ── 1. Parse PDB ──────────────────────────────────────────────────────────────
PDB_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "5SYF.pdb")
if not os.path.exists(PDB_PATH):
    PDB_PATH = PDB_PATH.replace(".pdb", ".cif")

print("=" * 55)
print("MFIOP minimal — PDB 5SYF pipeline")
print("=" * 55)

positions, dipoles = load_trp(PDB_PATH)
N = len(positions)
print(f"\n[1] Loaded {N} TRP residues from {os.path.basename(PDB_PATH)}")

# ── 2. Coupling summary ───────────────────────────────────────────────────────
print("\n[2] Shortest TRP–TRP distances:")
pairs = coupling_summary(positions, r_cutoff=25.0)
for p in pairs[:5]:
    print(f"    residue {p['i']} – {p['j']}: {p['r_Ang']:.2f} Å")

# ── 3. Build J (dimensionless) and convert to cm⁻¹ ──────────────────────────
J_dim = build_J(positions, dipoles)
J     = dimensionless_to_cm1(J_dim)

print("\n[3] Coupling matrix sanity check:")
print_sanity_check(J_dim)

# ── 4. Diagonalise ────────────────────────────────────────────────────────────
eigenvalues, eigenvectors = diagonalize(J)

print("\n[4] Collective modes:")
print(f"    {'k':>4}  {'λ_k (cm⁻¹)':>12}  {'ω_k (cm⁻¹)':>12}  {'|f_k|²_z':>10}  dominant")
print("    " + "-" * 55)
modes = mode_summary(eigenvalues, eigenvectors, dipoles, omega_0=OMEGA_0)
for m in modes:
    flag = "  ✓" if m["dominant"] else ""
    print(f"    {m['k']:>4}  {m['lambda_k']:>12.4f}  {m['omega_k']:>12.3f}  "
          f"{m['fk2_z']:>10.4f}{flag}")

# ── 5. Compute A(ω) ───────────────────────────────────────────────────────────
gamma_cm  = gamma_si_to_cm1()
omega_grid = make_omega_grid(center=OMEGA_0, half_width=35.0, n_points=4000)
result     = compute_spectrum(eigenvalues, eigenvectors, dipoles,
                               omega_grid, gamma=gamma_cm, omega_0=OMEGA_0)

print(f"\n[5] Spectral response:")
print(f"    γ           = {gamma_cm:.3f} cm⁻¹")
print(f"    Anisotropy  = {result['anisotropy']:.4f}  (expected ≈ 0.652)")

# Dominant modes inside physical window [1005, 1035] cm⁻¹
fk2_z   = result["fk2"]["z"]
omega_k = result["omega_k"]
in_win  = (omega_k >= 1005) & (omega_k <= 1035) & (fk2_z > 0.5)
n_dom   = int(in_win.sum())
print(f"    n_dominant (in window) = {n_dom}  (expected = 3)")

if n_dom >= 2:
    gaps = np.diff(np.sort(omega_k[in_win]))
    s_max = float(np.max(gaps)) if len(gaps) > 0 else 0.0
    print(f"    s_max       = {s_max:.3f} cm⁻¹")
    print(f"    s_max / γ   = {s_max/gamma_cm:.3f}  (> 1 → resolvable by nano-FTIR)")

# ── 6. Plot ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.patch.set_facecolor("#F8F9FA")

# Panel A: A(ω)
ax = axes[0]
ax.set_facecolor("#F8F9FA")
norm = result["A_z"].max()
ax.plot(omega_grid, result["A_z"] / norm,
        color="#E63946", lw=2.0, label="$A_z$ (z-axis)")
ax.plot(omega_grid, result["A_x"] / norm,
        color="#2E75B6", lw=1.0, alpha=0.65, label="$A_x$")
ax.plot(omega_grid, result["A_y"] / norm,
        color="#457B9D", lw=1.0, alpha=0.65, label="$A_y$")

# Mark dominant modes
for k, (wk, f2) in enumerate(zip(omega_k, fk2_z)):
    if f2 > 0.5:
        ax.axvline(wk, color="#E63946", ls=":", lw=0.9, alpha=0.7)

ax.axvspan(1005, 1035, alpha=0.06, color="#E63946", label="TRP ¹Lₐ window")
ax.set_xlim(omega_grid[0], omega_grid[-1])
ax.set_xlabel("ω  (cm⁻¹)", fontsize=11)
ax.set_ylabel("A(ω)  (normalised)", fontsize=11)
ax.set_title(f"A.  Anisotropic spectral response — PDB 5SYF\n"
             f"Anisotropy $A_z/(A_x+A_y)$ = {result['anisotropy']:.4f}",
             fontsize=10, fontweight="bold")
ax.legend(fontsize=9)
ax.grid(True, alpha=0.2, lw=0.5)

# Panel B: eigenvalue spectrum (|f_k|² vs λ_k)
ax2 = axes[1]
ax2.set_facecolor("#F8F9FA")
colors = ["#E63946" if m["dominant"] else "#AAAAAA" for m in modes]
ax2.bar([m["lambda_k"] for m in modes],
        [m["fk2_z"]    for m in modes],
        width=0.3, color=colors, alpha=0.85, edgecolor="white")
ax2.axhline(0.5, color="#E63946", ls="--", lw=1.0, label="|f_k|² threshold = 0.5")
ax2.set_xlabel("λ_k  (cm⁻¹)", fontsize=11)
ax2.set_ylabel("|f_k|²  (unnormalised)", fontsize=11)
ax2.set_title("B.  Oscillator strength projections\nRed = dominant modes (|f_k|² > 0.5)",
              fontsize=10, fontweight="bold")
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.2, lw=0.5)

fig.suptitle(
    "MFIOP minimal — Geometry → Spectrum mapping  |  PDB 5SYF  |  "
    "Veronica Zumpano Blumenfeld  |  DOI: 10.5281/zenodo.18166770",
    fontsize=9, y=1.01, color="#1A1A2E"
)
plt.tight_layout()

out_path = os.path.join(os.path.dirname(__file__), "..", "mfiop_5SYF_spectrum.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor="#F8F9FA")
print(f"\n[6] Figure saved: {out_path}")
plt.show()

print("\nDone.")
