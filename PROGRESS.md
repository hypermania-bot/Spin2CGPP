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
