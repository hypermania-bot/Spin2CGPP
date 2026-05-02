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

## 2026-05-01
- Problem encountered: the sequentially eliminated scalar equations were too large to compare directly to the paper's coefficient basis, and an initial helper put the published formulas into a private Mathematica package context, which made the numeric substitutions silently unreliable.
- Solution: added `Code/script/paper_scalar_coeffs.wl` with the published kinetic-sector formulas in the project `Global`` context, added `Code/script/check_scalar_paper_kinetics.wls` and `Code/test/test_scalar_paper_kinetics.wls`, and verified that the diagonal acceleration coefficients from the covariant scalar chain match the paper's `K_\varphi` and `K_B` after normalization. The same script also shows that lifting the acceleration matrix into the paper basis does not yet remove the first-derivative mismatch, which identifies the remaining gap as a reduced-action basis issue rather than a missing simplification.
- Avoid in future: when moving published formulas into a Mathematica package, bind explicitly to the intended symbol context before trusting any numeric regression output, and separate "equation basis" checks from "reduced-action basis" checks for constrained systems.
- Commit ID: `ac18f35`

## 2026-05-01
- Problem encountered: the earlier scalar reduction solved the `E,F` constraints from Euler equations and substituted them into a Lagrangian that still contained `E'` and `F'`, which obscured the reduced-action basis and left false coefficient mismatches.
- Solution: added `Code/script/derive_scalar_reduced_action.wls`, which first integrates derivatives off `E_\pm,F_\pm`, solves them algebraically, verifies the `A_+'A_-'` coefficient cancels, integrates derivatives off `A_\pm`, and then extracts the paper's coefficient basis. Extended `paper_scalar_coeffs.wl` to include `M_\varphi`, `M_B`, and `L_0`. The new full check proves zero symbolic residuals for `K_\varphi,M_\varphi,K_B,M_B,L_2,L_1,L_0`.
- Avoid in future: for auxiliary fields that appear through first derivatives, integrate by parts to make the auxiliary variables algebraic before substituting their constraints back into the action.
- Commit ID: `6fa5489`

## 2026-05-01
- Problem encountered: `Code/script/check_covariant_tensor_sector.wls` compared the covariant tensor reduction to the sector Lagrangian as a raw expression, so it printed a nonzero difference caused by an explicit boundary term and did not fail on the mismatch.
- Solution: rewrote the diagnostic to subtract the boundary term `d(a^3 H D^2)/d eta`, verify the reduced sector Lagrangian, verify the tensor mode Euler equation, and add the script to the regression suite.
- Avoid in future: covariant-to-sector checks should compare actions only after accounting for integrations by parts, and every diagnostic script that claims a check should return a nonzero exit code on failure.
- Commit ID: `859a5a5`

## 2026-05-01
- Problem encountered: the minimal scalar EOM note only summarized the covariant reduction and did not preserve the intermediate symbolic record needed to audit the derivation without re-running Mathematica.
- Solution: expanded `Notes/eom_notes.tex` with the scalar Fourier ansatz, background replacement rules, common denominator, reduced action basis, full coefficient polynomials, derived `\kappa` and `K_\mathcal{B}` combinations, symbolic residual identities, and a numerical spot check; synced and pushed the same note to Overleaf as commit `e6a9e62`.
- Avoid in future: when a long symbolic derivation is verified by scripts, record the compact definitions and polynomial coefficient basis in the note at the same time as the automated check.
- Commit ID: `70ad14b`

## 2026-05-02
- Problem encountered: the first nonminimal scalar coefficient check showed false `K_B` and `M_B` mismatches because top-level Mathematica assignments used leading `+` continuation lines; WolframScript evaluated those continuation lines as separate expressions, so the paper-side `P` and coefficient polynomials were silently truncated. A separate tensor diagnostic also initially used the wrong sign in the Friedmann substitution for `V(\bar\phi)`.
- Solution: rewrote the nonminimal paper coefficient transcription with explicit parenthesized right-hand sides, corrected the Friedmann substitution to `V = Mpl^2 (3 H^2 - Lambda) - phib^2/(2 a^2)`, and added covariant nonminimal tensor and scalar regression scripts. The scalar script now reduces `\mathcal L_{vv}^{(2)}` through the `E,F,A` constraints and verifies exact zero residuals for the paper's `K_B` and `M_B`.
- Avoid in future: in Wolfram scripts, wrap multi-line assigned expressions in parentheses or put binary operators at the end of the previous line; for background substitutions, derive `V` directly from the Friedmann equation before coding the rule.
- Commit ID: `18cf5af`
