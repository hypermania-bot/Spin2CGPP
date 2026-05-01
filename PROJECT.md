# Spin2CGPP

## Goal
Reproduce and check the mode equations in Siyang Ling et al., "Cosmological gravitational particle production of massive spin-2 particles" (`arXiv:2302.04390`), with end-to-end Mathematica scripts and detailed derivation notes.

## Structure
- `References/`: arXiv source for the target paper.
- `Code/script/`: Mathematica code that derives and verifies the published mode equations from the sector Lagrangians.
- `Code/test/`: regression checks for the derivations.
- `Notes/`: detailed derivation notes prepared for Overleaf upload.

## Next Steps
1. Reconstruct the minimally coupled scalar reduced-action basis itself, not just the final equation basis, so the full `{K_\varphi,M_\varphi,K_B,M_B,L_2,L_1,L_0}` match can be checked directly.
2. Extend `Code/script/paper_scalar_coeffs.wl` from the verified kinetic subset `{K_\varphi,K_B,L_2,L_1,\kappa,K_\mathcal{B}}` to the remaining mass-sector coefficients `{M_\varphi,M_B,L_0}`.
3. Use `Code/script/check_ds_scalar_limit.wls` and `Code/script/direct_metric_scalar_real.wls` to understand why naive Hessian-based elimination looks misleading even though the sequential `F -> A -> E` chain produces the expected two-field system.
