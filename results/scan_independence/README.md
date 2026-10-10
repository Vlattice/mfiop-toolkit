# Internal descriptor independence — 643-protein scan

Correlations across the MFIOP v2 scan (MFIOP_scan_v23.csv, n=643 human proteins).

## Result: MFIOP descriptors are internally independent (not redundant with each other)
| pair | Pearson r | reading |
|------|-----------|---------|
| anisotropy vs lambda_range | +0.046 | independent axes |
| anisotropy vs N_trp        | -0.061 | NOT a size artifact |
| anisotropy vs n_tight_norm | +0.183 | weakly related |
| lambda_range vs N_trp      | +0.231 | range has mild size dependence (use normalized) |
| n_tight_norm vs N_trp      | -0.398 | normalization slightly over-corrects |

## Where 5SYF sits in its own scan (percentiles, n=643)
- n_tight_norm: 22.9 | lambda_range: 20.8 | anisotropy: 41.8 | **ZB_v2: 77.3**
- 5SYF is mid/low on single marginals and stands out ONLY in the composite ZB.
  This matches the v3.0 text: the distinction is the joint geometry, not any marginal.

## Still open (needs PDBs + network)
Direct MFIOP-vs-amide-I correlation across the panel requires computing the amide-I
descriptor on each structure; not done here (no structures downloaded).
A 5-structure pilot gave range/N Pearson r=+0.91 (redundant on range) and PR/N
r=+0.64 (n too small). The panel test on many structures is the real pending step.
