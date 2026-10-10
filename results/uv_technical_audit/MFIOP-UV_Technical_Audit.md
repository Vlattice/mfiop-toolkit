# MFIOP-UV — Technical Audit
**Date:** 2026-10-10
**Scope:** UV radiative decay matrix, geometric descriptor, null model, and IR–UV relationship.

## 1. Matrix implementation
The UV radiative decay matrix was checked against the analytical small-sample limit:

  kr → 0,   Γ_ij/γ₀ → ê_i · ê_j.

The reported error for PDB 5SYF is 0.003. This supports the implementation and confirms
that, in the regime considered, the normalized radiative matrix closely approximates the
Gram matrix of dipole orientations.

## 2. Interpretation of S
For 5SYF, the reported descriptor is S = 5.06. Within the model, this is the largest
eigenvalue of the normalized collective radiative decay matrix. In the kr ≪ 1 regime, it
primarily captures collective alignment of the assumed dipole orientations.
The result is consistent with orientation-driven collective decay in the Dicke
small-sample limit. It should not be attributed to near-field radiative coupling.

## 3. Random-orientation null
The random-orientation null yields approximately z = 1.42, corresponding to roughly the
90th–92nd percentile. The observed geometry therefore produces a modest elevation relative
to the specified null. This is a descriptive result, not strong statistical evidence of an
exceptional effect.

## 4. Relationship between IR and UV operators
The dominant-mode overlap is 0.59, while the alignment of the three-dimensional dominant
subspaces is 0.35. These values do not support an identity of eigenstructure.
The defensible interpretation is that the IR and UV operators depend on the same underlying
molecular geometry but produce distinct modal structures. The phrase
"one eigenstructure, two spectroscopic manifestations" should be replaced with language
describing a shared geometric basis and distinct operator responses.

## 5. Dipole-orientation limitation
The use of the Cγ→CE2 axis as the transition-dipole direction for each protein-bound
tryptophan remains an assumption requiring validation. The reported sensitivity to angular
perturbations reinforces the importance of this limitation. Residue-specific TD-DFT or an
equivalent validated approach would be needed to establish the actual ¹La transition-dipole
directions.

## 6. Audit conclusion
The geometric descriptor and its implementation remain viable within the stated model. The
necessary revisions concern physical interpretation, statistical strength, and the
distinction between shared geometry and shared eigenstructure. These findings motivate
targeted claim corrections, not wholesale rejection of the manuscript.

**Next step:** revise the relevant claims and limitations in the UV manuscript after
confirming the reported calculations and documenting the null-model procedure.

---
*Reproducible: audit_uv.py, audit_full.py (this folder). Numbers verified on PDB 5SYF.*
