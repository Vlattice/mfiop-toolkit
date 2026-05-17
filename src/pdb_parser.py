"""
pdb_parser.py
-------------
Extract tryptophan (TRP) residue positions and dipole unit vectors from a PDB file.

Dipole axis: CG → CE2 (standard approximation for the ¹Lₐ transition of Trp).
Returns positions in Ångström and unit vectors (both as numpy arrays).
"""

import numpy as np
from Bio.PDB import PDBParser, MMCIFParser


def load_trp(pdb_path: str, model_index: int = 0):
    """
    Parse a PDB or mmCIF file and extract TRP residue data.

    Parameters
    ----------
    pdb_path : str
        Path to .pdb or .cif file.
    model_index : int
        Which MODEL to use (default: first model).

    Returns
    -------
    positions : np.ndarray, shape (N, 3)
        CG atom coordinates in Å.
    dipoles : np.ndarray, shape (N, 3)
        Unit vectors along CG → CE2 axis.

    Raises
    ------
    ValueError
        If fewer than 2 TRP residues with complete CG/CE2 atoms are found.
    """
    if pdb_path.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure("X", pdb_path)
    model = list(structure.get_models())[model_index]

    positions = []
    dipoles   = []

    for chain in model:
        for residue in chain:
            if residue.get_resname() != "TRP":
                continue
            try:
                cg  = residue["CG"].get_vector().get_array()
                ce2 = residue["CE2"].get_vector().get_array()
            except KeyError:
                continue  # skip residues missing atoms

            axis = ce2 - cg
            norm = np.linalg.norm(axis)
            if norm < 1e-6:
                continue

            positions.append(cg)
            dipoles.append(axis / norm)

    if len(positions) < 2:
        raise ValueError(
            f"Found only {len(positions)} complete TRP residue(s) in {pdb_path}. "
            "At least 2 are required."
        )

    return np.array(positions, dtype=float), np.array(dipoles, dtype=float)
