#!/usr/bin/env wolframscript

ClearAll["Global`*"];

eta = \[Eta];
u = uq;
spatial = Range[1, 3];

modeP = u;
modeM = 1/u;

ClearAll[PDz];
PDz[expr_] := Expand[I k u D[expr, u]];

ClearAll[PD];
PD[expr_, 0] := D[expr, eta];
PD[expr_, 1] := 0;
PD[expr_, 2] := 0;
PD[expr_, 3] := PDz[expr];

ClearAll[ModeField];
ModeField[qP_Symbol, qM_Symbol] := qP[eta] modeP + qM[eta] modeM;

aeta = a[eta];
metricDD = DiagonalMatrix[aeta^2 {-1, 1, 1, 1}];
metricUU = DiagonalMatrix[{-1, 1, 1, 1}/aeta^2];

ClearAll[GammaUDD];
GammaUDD[0, 0, 0] := aeta'/aeta;
GammaUDD[0, i_, j_] /; MemberQ[spatial, i] && MemberQ[spatial, j] := If[i == j, aeta'/aeta, 0];
GammaUDD[i_, 0, j_] /; MemberQ[spatial, i] && MemberQ[spatial, j] := If[i == j, aeta'/aeta, 0];
GammaUDD[i_, j_, 0] /; MemberQ[spatial, i] && MemberQ[spatial, j] := If[i == j, aeta'/aeta, 0];
GammaUDD[_, _, _] := 0;

ClearAll[zeroMode];
zeroMode[expr_] := Coefficient[Expand[u^2 expr], u, 2] // Expand;

normalizeWeirdDerivatives[expr_] := expr //. {
  HoldPattern[Derivative[n_][inner_]] /; ! MatchQ[inner, _Symbol] :> D[inner, {eta, n}],
  Derivative[1][f_ /; ! MatchQ[f, _Symbol]][arg_] :> D[f /. eta -> x, x] /. x -> arg,
  Derivative[2][f_ /; ! MatchQ[f, _Symbol]][arg_] :> D[f /. eta -> x, {x, 2}] /. x -> arg
};

vDD = ConstantArray[0, {4, 4}];
vDD[[1, 1]] = aeta^2 ModeField[Ep, Em];
vDD[[1, 4]] = aeta^2 PDz[ModeField[Fp, Fm]];
vDD[[4, 1]] = vDD[[1, 4]];
vDD[[2, 2]] = aeta^2 ModeField[Ap, Am];
vDD[[3, 3]] = aeta^2 ModeField[Ap, Am];
vDD[[4, 4]] = aeta^2 (ModeField[Ap, Am] + PDz[PDz[ModeField[Bp, Bm]]]);

varphi = ModeField[vp, vm];

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

gradPhiBar = {phib[eta], 0, 0, 0};
gradVar = Table[PD[varphi, mu], {mu, 0, 3}];
gradTrace = Table[PD[vTrace, mu], {mu, 0, 3}];
divUp = Table[Sum[CovDUp[mu, mu, nu], {mu, 0, 3}], {nu, 0, 3}];

term1 = -1/2 Sum[metricUU[[lam + 1, sig + 1]] CovDDown[lam, mu, nu] CovDUp[sig, mu, nu], {lam, 0, 3}, {sig, 0, 3}, {mu, 0, 3}, {nu, 0, 3}];
term2 = Sum[CovDUp[mu, nu, lam] CovDMixed[nu, mu, lam], {mu, 0, 3}, {nu, 0, 3}, {lam, 0, 3}];
term3 = -Sum[divUp[[nu + 1]] gradTrace[[nu + 1]], {nu, 0, 3}];
term4 = 1/2 Sum[metricUU[[mu + 1, nu + 1]] gradTrace[[mu + 1]] gradTrace[[nu + 1]], {mu, 0, 3}, {nu, 0, 3}];

(* On the Einstein-scalar FRW background, R_{mu nu} - Mpl^-2 d_mu phi d_nu phi is proportional to g_{mu nu}. *)
curvaturePrefactor = 3 H[eta]^2 - phib[eta]^2/(2 aeta^2 Mpl^2);
curvatureTensor = Table[curvaturePrefactor metricDD[[mu + 1, nu + 1]], {mu, 0, 3}, {nu, 0, 3}];
massiveQuadratic = Table[Sum[vUU[[mu + 1, lam + 1]] vDU[[lam + 1, nu + 1]], {lam, 0, 3}] - 1/2 vUU[[mu + 1, nu + 1]] vTrace, {mu, 0, 3}, {nu, 0, 3}];
term5 = Sum[curvatureTensor[[mu + 1, nu + 1]] massiveQuadratic[[mu + 1, nu + 1]], {mu, 0, 3}, {nu, 0, 3}];
term6 = -1/2 m^2 (Sum[vUU[[mu + 1, nu + 1]] vDD[[mu + 1, nu + 1]], {mu, 0, 3}, {nu, 0, 3}] - vTrace^2);
stressLike = Table[vUU[[mu + 1, nu + 1]] - 1/2 metricUU[[mu + 1, nu + 1]] vTrace, {mu, 0, 3}, {nu, 0, 3}];
term7 = (1/Mpl) Sum[(gradPhiBar[[mu + 1]] gradVar[[nu + 1]] + gradPhiBar[[nu + 1]] gradVar[[mu + 1]]) stressLike[[mu + 1, nu + 1]], {mu, 0, 3}, {nu, 0, 3}] - (Vp[eta] varphi vTrace)/Mpl;
term8 = -1/2 Sum[metricUU[[mu + 1, nu + 1]] gradVar[[mu + 1]] gradVar[[nu + 1]], {mu, 0, 3}, {nu, 0, 3}] - 1/2 Vpp[eta] varphi^2;

lagK = zeroMode[Expand[aeta^4 (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8)]];

backgroundRules = {
  Derivative[1][a][s_] :> a[s]^2 H[s],
  Derivative[2][a][s_] :> 2 a[s]^3 H[s]^2 + a[s]^2 Derivative[1][H][s],
  Derivative[1][phib][s_] :> -2 a[s] H[s] phib[s] - a[s]^2 Vp[s],
  Derivative[1][H][s_] :> -phib[s]^2/(2 a[s] Mpl^2),
  Derivative[1][Vp][s_] :> Vpp[s] phib[s]
};

lagKhat = Expand[
  lagK /. {
    vp -> Function[{s}, php[s] + phib[s] Ap[s]/(Mpl a[s] H[s])],
    vm -> Function[{s}, phm[s] + phib[s] Am[s]/(Mpl a[s] H[s])]
  }
];
lagKhat = normalizeWeirdDerivatives[lagKhat];
lagKhat = lagKhat //. backgroundRules;
lagKhat = normalizeWeirdDerivatives[lagKhat];

removeSecondDerivativeTerm[expr_, q_] := Module[
  {coeff = Coefficient[expr, q''[eta]]},
  Expand[expr - coeff q''[eta] - D[coeff, eta] q'[eta]]
];

cleanSecondDerivatives[expr_, fields_] := FixedPoint[
  Function[e, Fold[removeSecondDerivativeTerm, e, fields]],
  expr,
  6
];

lagKclean = cleanSecondDerivatives[lagKhat, {Ap, Am, Bp, Bm, Ep, Em, Fp, Fm, php, phm}];
lagKclean = normalizeWeirdDerivatives[lagKclean];
lagKclean = lagKclean //. backgroundRules;
lagKclean = normalizeWeirdDerivatives[lagKclean];

ClearAll[EulerMode];
EulerMode[L_, qminus_Symbol] := Numerator[Together[Expand[D[D[L, qminus'[eta]], eta] - D[L, qminus[eta]]]]];

eqE = EulerMode[lagKclean, Em];
eqF = EulerMode[lagKclean, Fm];
eqEp = Numerator[Together[Expand[D[D[lagKclean, Ep'[eta]], eta] - D[lagKclean, Ep[eta]]]]];
eqFp = Numerator[Together[Expand[D[D[lagKclean, Fp'[eta]], eta] - D[lagKclean, Fp[eta]]]]];

solEFPlus = Solve[{eqE == 0, eqF == 0}, {Ep[eta], Fp[eta]}];
solEFPlus = normalizeWeirdDerivatives[solEFPlus];
solEFPlus = solEFPlus //. backgroundRules;
solEFPlus = normalizeWeirdDerivatives[solEFPlus];

solEFMinus = Solve[{eqEp == 0, eqFp == 0}, {Em[eta], Fm[eta]}];
solEFMinus = normalizeWeirdDerivatives[solEFMinus];
solEFMinus = solEFMinus //. backgroundRules;
solEFMinus = normalizeWeirdDerivatives[solEFMinus];

lagAfterEF = Expand[lagKclean /. solEFPlus[[1]] /. solEFMinus[[1]]];
lagAfterEF = normalizeWeirdDerivatives[lagAfterEF];
lagAfterEF = lagAfterEF //. backgroundRules;
lagAfterEF = normalizeWeirdDerivatives[lagAfterEF];

hessian = Table[D[D[lagAfterEF, v1], v2], {v1, {Ap'[eta], Bp'[eta], php'[eta]}}, {v2, {Am'[eta], Bm'[eta], phm'[eta]}}] // Expand;
hessian = normalizeWeirdDerivatives[hessian];
hessian = hessian //. backgroundRules;
hessian = normalizeWeirdDerivatives[hessian];

eqA = EulerMode[lagAfterEF, Am];
eqA = normalizeWeirdDerivatives[eqA];
eqA = eqA //. backgroundRules;
eqA = normalizeWeirdDerivatives[eqA];

Print["lagK leaf count: ", LeafCount[lagKhat]];
Print["lagKclean leaf count: ", LeafCount[lagKclean]];
Print["lagKclean contains Ap'': ", ! FreeQ[lagKclean, Ap''[eta]]];
Print["solEFPlus count: ", Length[solEFPlus]];
Print["solEFMinus count: ", Length[solEFMinus]];
Print["Hessian(1,1): ", hessian[[1, 1]]];
Print["Hessian row1: ", Short[hessian[[1]], 10]];
Print["eqA contains Ap'': ", ! FreeQ[eqA, Ap''[eta]]];
Print["eqA contains Ap': ", ! FreeQ[eqA, Ap'[eta]]];
