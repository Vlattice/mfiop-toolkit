# Novelty control — 5SYF vs 3J6F (PRELIMINARY, not a validation)

A small-scale check of whether the tryptophan descriptor carries information
distinct from an amide-I backbone exciton model, on two tubulin structures.

## Honest status
- The raw range ratio (3J6F/5SYF ≈ 1.64) is largely a NETWORK-SIZE artifact
  (72 vs 8 Trp); once N is controlled, it does not indicate structural superiority.
- No statistically valid σ-level separation is claimed.
- The participation-ratio (localization) difference is interesting but not yet
  cleanly quantified.
- OPEN QUESTION: whether MFIOP adds information INDEPENDENT of amide-I cannot be
  answered with only 2 structures. It requires a multi-structure panel with
  descriptor-correlation / incremental-predictive-power analysis.

This folder documents the method and the negative/inconclusive result honestly;
it is NOT a validation of independence. Scripts: descriptor.py (arms 1+4),
topo.py (arm 2), amide3.py (arm 3), verify.py / verify2.py (error bars, N-control).
