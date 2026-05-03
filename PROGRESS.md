# Progress Log

## 2026-05-03
- Problem: Existing numerical code uses analytic eps→1.5 post-inflation continuation, which fails to capture oscillatory inflaton dynamics. This leads to: (a) high-k particle production that is too large (n_k rises instead of falling ∝ k^{-3/2}) and (b) absence of oscillatory spectral features from parametric resonance. Vector_nonmin scalar_nonmin omega² functions also had bugs in derivative computation.
- Solution: Fixed O2_vector_nonmin to use correct μ₁²=m²+H²(3+eps), μ₂²=m²+H²(3-2eps) and proper conformal-time derivatives (aH factor). Fixed _KB_MB to use correct ddotH = H³(2eps²-eps_p). Extended k-range to 0.03-80 a_e H_e with 100 modes for better peak resolution.
- Current status: tensor_min/vector_min/tensor_nonmin/vector_nonmin all give consistent results. vector_nonmin shows expected enhancement near gradient instability (m=2: ∫≈5.0, m=3: ∫≈1.3, m=4: ∫≈0.39, m=5: ∫≈0.28). Scalar_nonmin excluded due to ghost instability at accessible masses (paper only shows m≥7√2 H_inf).
- Remaining gap: High-k spectrum still rises with k (should fall ∝k^{-3/2}) due to missing inflaton oscillations in post-inflation background. Full oscillatory background is needed for paper-precision match. Peak n_k values are ~5-20× larger than paper's Fig. 3 due to flat high-k spectrum.
- Avoid in future: For nonminimal omega² formulas, always verify K'/K conversion between N-time and conformal time (factor of aH). For ddotH, the correct expression is ddotH = H³(2eps² - eps_p), not H³*2eps².

- Commit ID: `9fc6060`

## 2026-05-02 (afternoon)
- Problem: RK45 ODE integrator requires ~500k function evaluations per mode at rtol=1e-9, taking ~30s/mode. DOP853 and LSODA offer similar performance. Relaxing tolerance to rtol=1e-8 gives inaccurate |β|² (negative or wrong by factors of 10-100x at high k).
- Solution: Implemented exact-step (symplectic) midpoint integrator that evaluates Q² = ω²/(aH)² at each step midpoint and applies the exact harmonic oscillator solution. This gives perfect Wronskian conservation (Wr=1.000000) and runs ~0.5-3s per mode (10-60x faster than RK45). For k ≤ 10 a_e H_e, accuracy is within ~4% of the RK45 reference.
- Limitation: The fast integrator diverges from RK45 results at high k (k ≫ 10 a_e H_e) where Q varies rapidly. For these modes, the post-inflation analytic continuation may also miss oscillatory particle production features that depend on the full inflaton dynamics.
- Extended cgpp_numerics.py with new omega² functions: tensor_nonmin, vector_nonmin (via direct finite-difference), scalar_nonmin (KB/MB polynomial + finite-difference derivatives).
- Implemented new standalone code cgpp_full.py with all sectors and the fast symplectic integrator.
- Full spectra computed for 4 sectors × 4 masses (m/(√2 H_inf)=2,3,4,5) × 35 k-modes, saved to Code/output/cgpp_full_result.pkl.
- Key numerical results:
  - tensor_minimal ≈ vector_minimal ≈ tensor_nonmin: peak n_k ~ 6-8e-3 at k ≈ 25-50 a_e H_e
  - vector_nonmin: strongly enhanced for m/(√2 H_inf)=2 (peak n_k=1.1, 165x tensor) due to proximity to gradient instability. Enhancement decreases with mass: m=3→0.23, m=4→0.068, m=5→0.065.
  - Low-k spectra show oscillations in ln(k) from imaginary ν (m > 1.5 H_inf), consistent with paper's Appendix.
- Avoid in future: For oscillatory ODEs, use symplectic/exact-step methods rather than adaptive RK. Evaluate the frequency at the step midpoint for 2nd-order accuracy. For accurate |β|², the solution phase must be correct (simple Wronskian check is insufficient).
- Commit ID: `10182b8`

## 2026-05-02 (morning)
- Problem encountered: Numerical integration of CGPP mode equations is extremely challenging due to (a) enormous range of time scales (inflaton oscillation period ~10^5 Planck times vs Hubble time ~10^7), (b) rapidly oscillating mode functions at late times (m/H ≫ 1 in N-time), (c) conservation of the Wronskian is fragile with standard ODE integrators, and (d) the background ODE in N-time breaks down at the potential minimum (V→0 gives H²=0/0).
- Solutions: (1) Used dlnH/dN = -eps for background evolution, avoiding the zero-mass singularity. (2) Used Hamiltonian variable u = sqrt(aH) * chi which eliminates the friction term in the mode equation, resulting in perfect Wronskian conservation (W=1.000000). (3) Used an adiabatic stopping condition (|dΩ/dN|/Ω² < 0.01, N > N_end + 0.7) to stop integration when n_k has stabilized. (4) Started integration at N_end - 3 for all modes, balancing BD validity with computational efficiency.
- Avoid in future: Never use the V/(3-eps) formula for H² post-inflation; always evolve H independently. Never integrate mode equations in cosmic time with the friction term; use Hamiltonian variables or N-time + adiabatic stopping. Always check Wronskian conservation as a code-health metric.
- Commit ID: `0dd61a0` (main pipeline), `302ea72` (initial attempt)

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
