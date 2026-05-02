# Spin2CGPP

## Goal
Reproduce and check the mode equations in Siyang Ling et al., "Cosmological gravitational particle production of massive spin-2 particles" (`arXiv:2302.04390`), with end-to-end Mathematica scripts and detailed derivation notes.

## Structure
- `References/`: arXiv source for the target paper.
- `Code/script/`: Mathematica code that derives and verifies the published mode equations from the sector Lagrangians.
- `Code/test/`: regression checks for the derivations.
- `Notes/`: detailed derivation notes prepared for Overleaf upload.

## Next Steps
1. Treat `Code/script/derive_scalar_reduced_action.wls` as the canonical minimal-scalar reduced-action path; the full `{K_\varphi,M_\varphi,K_B,M_B,L_2,L_1,L_0}` match is checked symbolically.
2. Treat `Code/script/check_covariant_tensor_sector.wls` as the covariant minimal-tensor check; it removes the explicit boundary term and verifies the tensor mode equation.
3. Use `Code/script/check_ds_scalar_limit.wls` and `Code/script/direct_metric_scalar_real.wls` only as historical diagnostics for the misleading partially reduced Hessian route.
4. Optional next extension: add covariant-action vector-sector reductions. The tensor sectors and both scalar sectors now have covariant-to-paper regression checks.
