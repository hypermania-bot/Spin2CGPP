#!/usr/bin/env wolframscript

ClearAll["Global`*"];

eta = \[Eta];
z = zc;
spatial = Range[1, 3];

aeta = a[eta];
metricDD = DiagonalMatrix[aeta^2 {-1, 1, 1, 1}];
metricUU = DiagonalMatrix[{-1, 1, 1, 1}/aeta^2];

ClearAll[GammaUDD];
GammaUDD[0, 0, 0] := aeta'/aeta;
GammaUDD[0, i_, j_] /; MemberQ[spatial, i] && MemberQ[spatial, j] := If[i == j, aeta'/aeta, 0];
GammaUDD[i_, 0, j_] /; MemberQ[spatial, i] && MemberQ[spatial, j] := If[i == j, aeta'/aeta, 0];
GammaUDD[i_, j_, 0] /; MemberQ[spatial, i] && MemberQ[spatial, j] := If[i == j, aeta'/aeta, 0];
GammaUDD[_, _, _] := 0;

ClearAll[PD];
PD[expr_, 0] := D[expr, eta];
PD[expr_, 1] := 0;
PD[expr_, 2] := 0;
PD[expr_, 3] := D[expr, z];

mode = Cos[k z];

vDD = ConstantArray[0, {4, 4}];
vDD[[2, 3]] = aeta^2 TT[eta] mode;
vDD[[3, 2]] = vDD[[2, 3]];

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

ClearAll[RicciDD];
RicciDD[mu_, nu_] := RicciDD[mu, nu] = Expand[
  Sum[PD[GammaUDD[rho, mu, nu], rho], {rho, 0, 3}]
  - Sum[PD[GammaUDD[rho, mu, rho], nu], {rho, 0, 3}]
  + Sum[GammaUDD[rho, rho, sigma] GammaUDD[sigma, mu, nu], {rho, 0, 3}, {sigma, 0, 3}]
  - Sum[GammaUDD[rho, nu, sigma] GammaUDD[sigma, mu, rho], {rho, 0, 3}, {sigma, 0, 3}]
];

gradPhiBar = {phib[eta], 0, 0, 0};
gradTrace = Table[PD[vTrace, mu], {mu, 0, 3}];
divUp = Table[Sum[CovDUp[mu, mu, nu], {mu, 0, 3}], {nu, 0, 3}];

term1 = -1/2 Sum[metricUU[[lam + 1, sig + 1]] CovDDown[lam, mu, nu] CovDUp[sig, mu, nu], {lam, 0, 3}, {sig, 0, 3}, {mu, 0, 3}, {nu, 0, 3}];
term2 = Sum[CovDUp[mu, nu, lam] CovDMixed[nu, mu, lam], {mu, 0, 3}, {nu, 0, 3}, {lam, 0, 3}];
term3 = -Sum[divUp[[nu + 1]] gradTrace[[nu + 1]], {nu, 0, 3}];
term4 = 1/2 Sum[metricUU[[mu + 1, nu + 1]] gradTrace[[mu + 1]] gradTrace[[nu + 1]], {mu, 0, 3}, {nu, 0, 3}];

curvatureTensor = Table[RicciDD[mu, nu] - gradPhiBar[[mu + 1]] gradPhiBar[[nu + 1]]/Mpl^2, {mu, 0, 3}, {nu, 0, 3}];
massiveQuadratic = Table[Sum[vUU[[mu + 1, lam + 1]] vDU[[lam + 1, nu + 1]], {lam, 0, 3}] - 1/2 vUU[[mu + 1, nu + 1]] vTrace, {mu, 0, 3}, {nu, 0, 3}];
term5 = Sum[curvatureTensor[[mu + 1, nu + 1]] massiveQuadratic[[mu + 1, nu + 1]], {mu, 0, 3}, {nu, 0, 3}];
term6 = -1/2 m^2 (Sum[vUU[[mu + 1, nu + 1]] vDD[[mu + 1, nu + 1]], {mu, 0, 3}, {nu, 0, 3}] - vTrace^2);

lagDensity = Expand[aeta^4 (term1 + term2 + term3 + term4 + term5 + term6)];
lagReduced = FullSimplify[(k/(2 Pi)) Integrate[lagDensity, {z, 0, 2 Pi/k}], Assumptions -> {k > 0, aeta > 0, Mpl > 0, m > 0}];

backgroundRules = {
  Derivative[1][a][s_] :> a[s]^2 H[s],
  Derivative[2][a][s_] :> 2 a[s]^3 H[s]^2 + a[s]^2 Derivative[1][H][s]
};

lagReduced = Expand[lagReduced /. {
    Derivative[1][a[eta]] :> Derivative[1][a][eta],
    Derivative[1][Derivative[1][a][eta]] :> Derivative[2][a][eta]
  }] //. backgroundRules;

ClearAll[eulerEquation, simplifyBackground];
simplifyBackground[expr_] := FullSimplify[Expand[expr] //. backgroundRules];
eulerEquation[lag_, q_] := simplifyBackground[D[D[lag, q'[eta]], eta] - D[lag, q[eta]]];

boundaryTerm = simplifyBackground[D[a[eta]^3 H[eta] TT[eta]^2, eta]];
sectorLagrangian = 1/2 a[eta]^2 TT'[eta]^2 - 1/2 a[eta]^2 (k^2 + a[eta]^2 m^2) TT[eta]^2;
sectorDifference = simplifyBackground[lagReduced - boundaryTerm - sectorLagrangian];

modeEquation = eulerEquation[sectorLagrangian, TT];
expectedModeEquation = simplifyBackground[
  a[eta]^2 (TT''[eta] + 2 a[eta] H[eta] TT'[eta] + (k^2 + a[eta]^2 m^2) TT[eta])
];
modeEquationDifference = simplifyBackground[modeEquation - expectedModeEquation];

Print["sector difference after boundary term: ", sectorDifference];
Print["mode equation difference: ", modeEquationDifference];

If[! TrueQ[sectorDifference === 0] || ! TrueQ[modeEquationDifference === 0],
  Exit[1],
  Exit[0]
];
