# Progress Log

## 2026-04-29
- Problem encountered: direct Lagrangian comparisons after canonical field redefinitions produced nonzero residuals because the transformed Lagrangians differ by total derivatives.
- Solution: changed the Mathematica checks to compare Euler-Lagrange equations and frequency identities instead of raw Lagrangian expressions.
- Avoid in future: when verifying canonical normalization steps, compare equations of motion or remove total derivatives explicitly before asserting equality.
- Commit ID: `535d51e`

## 2026-04-29
- Problem encountered: the minimally coupled scalar two-field system carries a derivative mixing term, and a naive real-field variation leaves a `\sigma_1'/2` convention shift relative to the paper's complex-mode presentation.
- Solution: documented the convention clearly in the notes and encoded the check in that same convention so the symbolic verification remains precise and reproducible.
- Avoid in future: for complex Fourier modes with off-diagonal derivative couplings, decide the conjugation and integration-by-parts convention before writing regression tests.
- Commit ID: `535d51e`

## 2026-04-30
- Problem encountered: a naive Hessian analysis of the covariant minimal scalar reduction incorrectly suggested that the scalar field `A` remained dynamical after eliminating `E,F`.
- Solution: switched to the correct sequential constraint chain `F_\pm -> A_\pm -> E_\pm`, implemented it in `Code/script/derive_scalar_mode_eoms.wls`, and added a regression test that confirms the final coupled `(B,\hat{\varphi}_v)` system has a rank-2 acceleration matrix.
- Avoid in future: for constrained spin-2 scalar systems, do not infer the propagating content from a partially reduced Hessian alone; first eliminate the multiplier equations in the order implied by the constraints.
- Commit ID: `2743131`
