# Spin2CGPP

## Goal
Reproduce and check the mode equations in Siyang Ling et al., "Cosmological gravitational particle production of massive spin-2 particles" (`arXiv:2302.04390`), with end-to-end Mathematica scripts and detailed derivation notes.

## Structure
- `References/`: arXiv source for the target paper.
- `Code/script/`: Mathematica code that derives and verifies the published mode equations from the sector Lagrangians.
- `Code/test/`: regression checks for the derivations.
- `Notes/`: detailed derivation notes prepared for Overleaf upload.

## Next Steps
1. Simplify the final coupled scalar equations from `Code/script/derive_scalar_mode_eoms.wls` into the paper's published `{K_\varphi,K_B,L_2,L_1,L_0,\dots}` coefficient basis.
2. Use `Code/script/check_ds_scalar_limit.wls` and `Code/script/direct_metric_scalar_real.wls` to understand why naive Hessian-based elimination looks misleading even though the sequential `F -> A -> E` chain produces the expected two-field system.
3. Update the notes with the successful sequential elimination chain and add a comparison against the closed-form equations in the paper.
