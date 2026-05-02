#!/usr/bin/env wolframscript

ClearAll["Global`*"];

eta = \[Eta];
u = uq;
spatial = Range[1, 3];
modeP = u;
modeM = 1/u;

ClearAll[PDz, PD, ModeField, zeroMode, normalizeWeirdDerivatives];
PDz[expr_] := Expand[I k u D[expr, u]];
PD[expr_, 0] := D[expr, eta];
PD[expr_, 1] := 0;
PD[expr_, 2] := 0;
PD[expr_, 3] := PDz[expr];
ModeField[qP_Symbol, qM_Symbol] := qP[eta] modeP + qM[eta] modeM;
zeroMode[expr_] := Coefficient[Expand[u^2 expr], u, 2] // Expand;
normalizeWeirdDerivatives[expr_] := expr //. {
  HoldPattern[Derivative[n_][inner_]] /; ! MatchQ[inner, _Symbol] :> D[inner, {eta, n}],
  Derivative[1][f_ /; ! MatchQ[f, _Symbol]][arg_] :> D[f /. eta -> x, x] /. x -> arg,
  Derivative[2][f_ /; ! MatchQ[f, _Symbol]][arg_] :> D[f /. eta -> x, {x, 2}] /. x -> arg
};

aeta = a[eta];
heta = H[eta];
hDot = H'[eta]/aeta;
hDDot = H''[eta]/aeta^2 - heta H'[eta]/aeta;

metricDD = DiagonalMatrix[aeta^2 {-1, 1, 1, 1}];
metricUU = DiagonalMatrix[{-1, 1, 1, 1}/aeta^2];

ClearAll[GammaUDD];
GammaUDD[0, 0, 0] := aeta'/aeta;
GammaUDD[0, i_, j_] /; MemberQ[spatial, i] && MemberQ[spatial, j] := If[i == j, aeta'/aeta, 0];
GammaUDD[i_, 0, j_] /; MemberQ[spatial, i] && MemberQ[spatial, j] := If[i == j, aeta'/aeta, 0];
GammaUDD[i_, j_, 0] /; MemberQ[spatial, i] && MemberQ[spatial, j] := If[i == j, aeta'/aeta, 0];
GammaUDD[_, _, _] := 0;

vDD = ConstantArray[0, {4, 4}];
vDD[[1, 1]] = aeta^2 ModeField[Ep, Em];
vDD[[1, 4]] = aeta^2 PDz[ModeField[Fp, Fm]];
vDD[[4, 1]] = vDD[[1, 4]];
vDD[[2, 2]] = aeta^2 ModeField[Ap, Am];
vDD[[3, 3]] = aeta^2 ModeField[Ap, Am];
vDD[[4, 4]] = aeta^2 (ModeField[Ap, Am] + PDz[PDz[ModeField[Bp, Bm]]]);

vUU = Table[
  Sum[metricUU[[mu + 1, alpha + 1]] metricUU[[nu + 1, beta + 1]] vDD[[alpha + 1, beta + 1]], {alpha, 0, 3}, {beta, 0, 3}],
  {mu, 0, 3}, {nu, 0, 3}
];
vUD = Table[
  Sum[metricUU[[mu + 1, alpha + 1]] vDD[[alpha + 1, nu + 1]], {alpha, 0, 3}],
  {mu, 0, 3}, {nu, 0, 3}
];
vDU = Table[
  Sum[vDD[[mu + 1, alpha + 1]] metricUU[[alpha + 1, nu + 1]], {alpha, 0, 3}],
  {mu, 0, 3}, {nu, 0, 3}
];
vTrace = Sum[metricUU[[mu + 1, nu + 1]] vDD[[mu + 1, nu + 1]], {mu, 0, 3}, {nu, 0, 3}] // Expand;

ClearAll[CovDDown, CovDUp, CovDMixed];
CovDDown[lam_, mu_, nu_] := CovDDown[lam, mu, nu] = Expand[
  PD[vDD[[mu + 1, nu + 1]], lam]
  - Sum[GammaUDD[rho, lam, mu] vDD[[rho + 1, nu + 1]], {rho, 0, 3}]
  - Sum[GammaUDD[rho, lam, nu] vDD[[mu + 1, rho + 1]], {rho, 0, 3}]
];
CovDUp[lam_, mu_, nu_] := CovDUp[lam, mu, nu] = Expand[
  PD[vUU[[mu + 1, nu + 1]], lam]
  + Sum[GammaUDD[mu, lam, rho] vUU[[rho + 1, nu + 1]], {rho, 0, 3}]
  + Sum[GammaUDD[nu, lam, rho] vUU[[mu + 1, rho + 1]], {rho, 0, 3}]
];
CovDMixed[lam_, mu_, nu_] := CovDMixed[lam, mu, nu] = Expand[
  PD[vUD[[mu + 1, nu + 1]], lam]
  + Sum[GammaUDD[mu, lam, rho] vUD[[rho + 1, nu + 1]], {rho, 0, 3}]
  - Sum[GammaUDD[rho, lam, nu] vUD[[mu + 1, rho + 1]], {rho, 0, 3}]
];

gradTrace = Table[PD[vTrace, mu], {mu, 0, 3}];
divUp = Table[Sum[CovDUp[mu, mu, nu], {mu, 0, 3}], {nu, 0, 3}];

term1 = -1/2 Sum[metricUU[[lam + 1, sig + 1]] CovDDown[lam, mu, nu] CovDUp[sig, mu, nu], {lam, 0, 3}, {sig, 0, 3}, {mu, 0, 3}, {nu, 0, 3}];
term2 = Sum[CovDUp[mu, nu, lam] CovDMixed[nu, mu, lam], {mu, 0, 3}, {nu, 0, 3}, {lam, 0, 3}];
term3 = -Sum[divUp[[nu + 1]] gradTrace[[nu + 1]], {nu, 0, 3}];
term4 = 1/2 Sum[metricUU[[mu + 1, nu + 1]] gradTrace[[mu + 1]] gradTrace[[nu + 1]], {mu, 0, 3}, {nu, 0, 3}];

(* The nonminimal scalar coefficients quoted in the paper set Lambda = 0. *)
curvA = DiagonalMatrix[
  {
    -3 aeta^2 heta^2/2 - 3 aeta H'[eta] - aeta^2 Lambda/2,
    aeta^2 (3 heta^2 + Lambda)/2,
    aeta^2 (3 heta^2 + Lambda)/2,
    aeta^2 (3 heta^2 + Lambda)/2
  }
];
curvB = DiagonalMatrix[
  {
    -3 aeta H'[eta] - aeta^2 Lambda,
    aeta^2 Lambda - aeta H'[eta],
    aeta^2 Lambda - aeta H'[eta],
    aeta^2 Lambda - aeta H'[eta]
  }
];
vSquaredMixed = Table[
  Sum[vUU[[mu + 1, lam + 1]] vDU[[lam + 1, nu + 1]], {lam, 0, 3}],
  {mu, 0, 3}, {nu, 0, 3}
];
term5 = Sum[curvA[[mu + 1, nu + 1]] vSquaredMixed[[mu + 1, nu + 1]], {mu, 0, 3}, {nu, 0, 3}];
term6 = -1/2 Sum[curvB[[mu + 1, nu + 1]] vUU[[mu + 1, nu + 1]] vTrace, {mu, 0, 3}, {nu, 0, 3}];
term7 = -1/2 m^2 (Sum[vUU[[mu + 1, nu + 1]] vDD[[mu + 1, nu + 1]], {mu, 0, 3}, {nu, 0, 3}] - vTrace^2);

lagK = zeroMode[Expand[aeta^4 (term1 + term2 + term3 + term4 + term5 + term6 + term7)]];

backgroundRules = {
  Lambda -> 0,
  Derivative[1][a][s_] :> a[s]^2 H[s],
  Derivative[2][a][s_] :> 2 a[s]^3 H[s]^2 + a[s]^2 Derivative[1][H][s],
  Derivative[3][a][s_] :> 6 a[s]^4 H[s]^3 + 6 a[s]^3 H[s] Derivative[1][H][s] + a[s]^2 Derivative[2][H][s]
};

ClearAll[simplifyOnBackground, removeSecondDerivativeTerm, cleanSecondDerivatives, removeFirstDerivativeTerm, removeDerivatives];
simplifyOnBackground[expr_] := Cancel[Together[normalizeWeirdDerivatives[expr //. backgroundRules]]];

removeSecondDerivativeTerm[expr_, q_] := Module[
  {coeff = simplifyOnBackground[Coefficient[expr, q''[eta]]]},
  simplifyOnBackground[expr - coeff q''[eta] - D[coeff, eta] q'[eta]]
];

cleanSecondDerivatives[expr_, fields_] := FixedPoint[
  Function[e, Fold[removeSecondDerivativeTerm, e, fields]],
  simplifyOnBackground[expr],
  6
];

removeFirstDerivativeTerm[expr_, q_] := Module[
  {coeff = simplifyOnBackground[Coefficient[expr, q'[eta]]]},
  simplifyOnBackground[expr - coeff q'[eta] - D[coeff, eta] q[eta]]
];

removeDerivatives[expr_, fields_] := FixedPoint[
  Function[e, simplifyOnBackground[Fold[removeFirstDerivativeTerm, e, fields]]],
  expr,
  4
];

lagKclean = cleanSecondDerivatives[lagK, {Ap, Am, Bp, Bm, Ep, Em, Fp, Fm}];

lagAuxiliaryAlgebraic = removeDerivatives[lagKclean, {Ep, Em, Fp, Fm}];
auxiliaryEquations = Table[
  D[lagAuxiliaryAlgebraic, q[eta]],
  {q, {Em, Fm, Ep, Fp}}
] // Expand;
auxiliarySolution = Solve[
  Thread[auxiliaryEquations == 0],
  {Ep[eta], Fp[eta], Em[eta], Fm[eta]}
][[1]];

lagAfterAuxiliaries = simplifyOnBackground[Expand[lagAuxiliaryAlgebraic /. auxiliarySolution]];
aKineticResidual = simplifyOnBackground[Coefficient[lagAfterAuxiliaries, Ap'[eta] Am'[eta]]];

lagAAlgebraic = removeDerivatives[lagAfterAuxiliaries, {Ap, Am}];
aEquations = Table[D[lagAAlgebraic, q[eta]], {q, {Am, Ap}}] // Expand;
aSolution = Solve[Thread[aEquations == 0], {Ap[eta], Am[eta]}][[1]];

lagScalarReduced = simplifyOnBackground[Expand[lagAAlgebraic /. aSolution]];
bPrimeB = simplifyOnBackground[Coefficient[lagScalarReduced, Bp'[eta] Bm[eta]]];

nonminimalScalarReducedCoefficients = <|
  "KB" -> simplifyOnBackground[Coefficient[lagScalarReduced, Bp'[eta] Bm'[eta]]/2],
  "MB" -> simplifyOnBackground[-(Coefficient[lagScalarReduced, Bp[eta] Bm[eta]] - D[bPrimeB, eta])/2]
|>;

paperP = (
  4 (m^2 + 3 heta^2 + 3 hDot) k^4
  + 12 aeta^2 ((m^2 + heta^2) (m^2 + 3 heta^2) + 2 (m^2 + heta^2) hDot - hDot^2) k^2
  + 9 aeta^4 (m^2 + heta^2) (m^2 + 3 heta^2 - hDot) (m^2 + 3 heta^2 + 2 hDot)
);

paperKB = aeta^4 (
  -4 hDot^2 k^6
  + 3 aeta^2 (m^2 + heta^2) (m^2 + 3 heta^2 - hDot) (m^2 + 3 heta^2 + 2 hDot) k^4
)/paperP;

paperC10 = (
  12 (m^2 + heta^2) (m^2 + 3 heta^2)^3
  + 16 (m^2 + 3 heta^2)^2 (6 m^2 + 7 heta^2) hDot
  + 4 (m^2 + 3 heta^2) (63 m^2 + 71 heta^2) hDot^2
  + 8 (25 m^2 + 27 heta^2) hDot^3
  - 48 hDot^4
  - 32 heta (m^2 + 3 heta^2) hDot hDDot
  - 48 heta hDot^2 hDDot
);

paperC8 = (
  12 aeta^2 (m^2 + 3 heta^2 + 2 hDot) (
    2 (m^2 + heta^2) (m^2 + 3 heta^2)^2 (2 m^2 + 5 heta^2)
    + (m^2 + 3 heta^2) (19 m^4 + 64 m^2 heta^2 + 49 heta^4) hDot
    + 2 (7 m^4 + 20 m^2 heta^2 + 17 heta^4) hDot^2
    - (23 m^2 + 25 heta^2) hDot^3
    + 2 hDot^4
    - 2 heta (m^2 + heta^2) (m^2 + 3 heta^2) hDDot
    - 4 heta (m^2 + heta^2) hDot hDDot
  )
);

paperC6 = (
  9 aeta^4 (m^2 + heta^2) (m^2 + 3 heta^2 - hDot) (m^2 + 3 heta^2 + 2 hDot)^2
  * (7 (m^2 + heta^2) (m^2 + 3 heta^2) + 17 (m^2 + heta^2) hDot - 8 hDot^2)
);

paperC4 = (
  27 aeta^6 (m^2 + heta^2)^2 (m^2 + 3 heta^2 - hDot)^2 (m^2 + 3 heta^2 + 2 hDot)^3
);

paperMB = aeta^6 (paperC10 k^10 + paperC8 k^8 + paperC6 k^6 + paperC4 k^4)/paperP^2;

paperNonminimalScalarCoefficients = <|
  "P" -> simplifyOnBackground[paperP],
  "KB" -> simplifyOnBackground[paperKB],
  "MB" -> simplifyOnBackground[paperMB]
|>;

nonminimalScalarCoefficientResiduals = <|
  "KB" -> simplifyOnBackground[nonminimalScalarReducedCoefficients["KB"] - paperNonminimalScalarCoefficients["KB"]],
  "MB" -> simplifyOnBackground[nonminimalScalarReducedCoefficients["MB"] - paperNonminimalScalarCoefficients["MB"]]
|>;

If[$InputFileName =!= "" && Environment["SPIN2CGPP_LIBRARY_MODE"] =!= "1",
  Print["auxiliary derivative flags: ", Table[! FreeQ[lagAuxiliaryAlgebraic, q'[eta]], {q, {Ep, Em, Fp, Fm}}]];
  Print["A kinetic residual zero: ", TrueQ[aKineticResidual === 0]];
  Print["reduced action leaf count: ", LeafCount[lagScalarReduced]];
  Print["coefficient residuals zero: ", Map[TrueQ[# === 0] &, nonminimalScalarCoefficientResiduals]];
  If[And @@ (TrueQ[# === 0] & /@ Values[nonminimalScalarCoefficientResiduals]),
    Exit[0],
    Exit[1]
  ];
];
