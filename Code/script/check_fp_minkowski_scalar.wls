#!/usr/bin/env wolframscript

ClearAll["Global`*"];

eta = \[Eta];
u = uq;
modeP = u;
modeM = 1/u;

ClearAll[PDz, PD, ModeField, zeroMode];
PDz[expr_] := Expand[I k u D[expr, u]];
PD[expr_, 0] := D[expr, eta];
PD[expr_, 1] := 0;
PD[expr_, 2] := 0;
PD[expr_, 3] := PDz[expr];
ModeField[qP_Symbol, qM_Symbol] := qP[eta] modeP + qM[eta] modeM;
zeroMode[expr_] := Coefficient[Expand[u^2 expr], u, 2] // Expand;

metricDD = DiagonalMatrix[{-1, 1, 1, 1}];
metricUU = metricDD;

vDD = ConstantArray[0, {4, 4}];
vDD[[1, 1]] = ModeField[Ep, Em];
vDD[[1, 4]] = PDz[ModeField[Fp, Fm]];
vDD[[4, 1]] = vDD[[1, 4]];
vDD[[2, 2]] = ModeField[Ap, Am];
vDD[[3, 3]] = ModeField[Ap, Am];
vDD[[4, 4]] = ModeField[Ap, Am] + PDz[PDz[ModeField[Bp, Bm]]];

vUU = Table[
  Sum[metricUU[[mu + 1, al + 1]] metricUU[[nu + 1, be + 1]] vDD[[al + 1, be + 1]], {al, 0, 3}, {be, 0, 3}],
  {mu, 0, 3}, {nu, 0, 3}
];
vUD = Table[
  Sum[metricUU[[mu + 1, al + 1]] vDD[[al + 1, nu + 1]], {al, 0, 3}],
  {mu, 0, 3}, {nu, 0, 3}
];
vTrace = Sum[metricUU[[mu + 1, nu + 1]] vDD[[mu + 1, nu + 1]], {mu, 0, 3}, {nu, 0, 3}] // Expand;

CovDDown[lam_, mu_, nu_] := PD[vDD[[mu + 1, nu + 1]], lam];
CovDUp[lam_, mu_, nu_] := PD[vUU[[mu + 1, nu + 1]], lam];
CovDMixed[lam_, mu_, nu_] := PD[vUD[[mu + 1, nu + 1]], lam];

divUp = Table[Sum[CovDUp[mu, mu, nu], {mu, 0, 3}], {nu, 0, 3}];
gradTrace = Table[PD[vTrace, mu], {mu, 0, 3}];

term1 = -1/2 Sum[metricUU[[lam + 1, sig + 1]] CovDDown[lam, mu, nu] CovDUp[sig, mu, nu], {lam, 0, 3}, {sig, 0, 3}, {mu, 0, 3}, {nu, 0, 3}];
term2 = Sum[CovDUp[mu, nu, lam] CovDMixed[nu, mu, lam], {mu, 0, 3}, {nu, 0, 3}, {lam, 0, 3}];
term3 = -Sum[divUp[[nu + 1]] gradTrace[[nu + 1]], {nu, 0, 3}];
term4 = 1/2 Sum[metricUU[[mu + 1, nu + 1]] gradTrace[[mu + 1]] gradTrace[[nu + 1]], {mu, 0, 3}, {nu, 0, 3}];
term6 = -1/2 m^2 (Sum[vUU[[mu + 1, nu + 1]] vDD[[mu + 1, nu + 1]], {mu, 0, 3}, {nu, 0, 3}] - vTrace^2);

lag = zeroMode[Expand[term1 + term2 + term3 + term4 + term6]];

ClearAll[EulerMode];
EulerMode[L_, qminus_Symbol] := Numerator[Together[Expand[D[D[L, qminus'[eta]], eta] - D[L, qminus[eta]]]]];

eqA = EulerMode[lag, Am];
hess = Table[Coefficient[Coefficient[lag, q1], q2], {q1, {Ap'[eta], Bp'[eta]}}, {q2, {Am'[eta], Bm'[eta]}}] // Expand;

Print["Minkowski scalar Hessian = ", hess];
Print["Minkowski scalar determinant = ", Det[hess] // Expand];
Print["Minkowski eqA contains Ap'': ", ! FreeQ[eqA, Ap''[eta]]];
