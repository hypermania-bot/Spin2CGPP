#!/usr/bin/env wolframscript

ClearAll["Global`*"];

eta = \[Eta];
z = zc;

coords = {eta, x, y, z};
spatial = Range[1, 3];

(* Background quantities *)
aeta = a[eta];
phiBarP = phib[eta];
VpEta = Vp[eta];
VppEta = Vpp[eta];

metricDD = DiagonalMatrix[aeta^2 {-1, 1, 1, 1}];
metricUU = DiagonalMatrix[{-1, 1, 1, 1}/aeta^2];

ClearAll[GammaUDD];
GammaUDD[0, 0, 0] := aeta' / aeta;
GammaUDD[0, i_, j_] /; MemberQ[spatial, i] && MemberQ[spatial, j] := If[i == j, aeta' / aeta, 0];
GammaUDD[i_, 0, j_] /; MemberQ[spatial, i] && MemberQ[spatial, j] := If[i == j, aeta' / aeta, 0];
GammaUDD[i_, j_, 0] /; MemberQ[spatial, i] && MemberQ[spatial, j] := If[i == j, aeta' / aeta, 0];
GammaUDD[_, _, _] := 0;

ClearAll[PD];
PD[expr_, 0] := D[expr, eta];
PD[expr_, 1] := 0;
PD[expr_, 2] := 0;
PD[expr_, 3] := D[expr, z];

mode = Cos[k z];
varphi = pv[eta] mode;

ClearAll[vDD];
vDD = ConstantArray[0, {4, 4}];
vDD[[1, 1]] = aeta^2 EE[eta] mode;
vDD[[1, 4]] = aeta^2 D[FF[eta] mode, z];
vDD[[4, 1]] = vDD[[1, 4]];
vDD[[2, 2]] = aeta^2 AA[eta] mode;
vDD[[3, 3]] = aeta^2 AA[eta] mode;
vDD[[4, 4]] = aeta^2 (AA[eta] mode + D[D[BB[eta] mode, z], z]);

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

gradPhiBar = {phiBarP, 0, 0, 0};
gradVar = Table[PD[varphi, mu], {mu, 0, 3}];
gradTrace = Table[PD[vTrace, mu], {mu, 0, 3}];
divUp = Table[Sum[CovDUp[mu, mu, nu], {mu, 0, 3}], {nu, 0, 3}];

term1 = -1/2 Sum[
  metricUU[[lam + 1, sig + 1]] CovDDown[lam, mu, nu] CovDUp[sig, mu, nu],
  {lam, 0, 3}, {sig, 0, 3}, {mu, 0, 3}, {nu, 0, 3}
];

term2 = Sum[
  CovDUp[mu, nu, lam] CovDMixed[nu, mu, lam],
  {mu, 0, 3}, {nu, 0, 3}, {lam, 0, 3}
];

term3 = -Sum[divUp[[nu + 1]] gradTrace[[nu + 1]], {nu, 0, 3}];

term4 = 1/2 Sum[
  metricUU[[mu + 1, nu + 1]] gradTrace[[mu + 1]] gradTrace[[nu + 1]],
  {mu, 0, 3}, {nu, 0, 3}
];

(* On the Einstein-scalar FRW background, R_{mu nu} - Mpl^-2 d_mu phi d_nu phi is proportional to g_{mu nu}. *)
curvaturePrefactor = 3 H[eta]^2 - phib[eta]^2/(2 aeta^2 Mpl^2);
curvatureTensor = Table[curvaturePrefactor metricDD[[mu + 1, nu + 1]], {mu, 0, 3}, {nu, 0, 3}];

massiveQuadratic = Table[
  Sum[vUU[[mu + 1, lam + 1]] vDU[[lam + 1, nu + 1]], {lam, 0, 3}] - 1/2 vUU[[mu + 1, nu + 1]] vTrace,
  {mu, 0, 3}, {nu, 0, 3}
];

term5 = Sum[curvatureTensor[[mu + 1, nu + 1]] massiveQuadratic[[mu + 1, nu + 1]], {mu, 0, 3}, {nu, 0, 3}];

term6 = -1/2 m^2 (Sum[vUU[[mu + 1, nu + 1]] vDD[[mu + 1, nu + 1]], {mu, 0, 3}, {nu, 0, 3}] - vTrace^2);

stressLike = Table[
  vUU[[mu + 1, nu + 1]] - 1/2 metricUU[[mu + 1, nu + 1]] vTrace,
  {mu, 0, 3}, {nu, 0, 3}
];

term7 = (1/Mpl) Sum[
  (gradPhiBar[[mu + 1]] gradVar[[nu + 1]] + gradPhiBar[[nu + 1]] gradVar[[mu + 1]]) stressLike[[mu + 1, nu + 1]],
  {mu, 0, 3}, {nu, 0, 3}
] - (VpEta varphi vTrace)/Mpl;

term8 = -1/2 Sum[metricUU[[mu + 1, nu + 1]] gradVar[[mu + 1]] gradVar[[nu + 1]], {mu, 0, 3}, {nu, 0, 3}] - 1/2 VppEta varphi^2;

lagrangianDensity = Expand[aeta^4 (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8)];

averageRules = {
  Integrate[Cos[k z]^2, {z, 0, 2 Pi/k}] -> Pi/k,
  Integrate[Sin[k z]^2, {z, 0, 2 Pi/k}] -> Pi/k,
  Integrate[Cos[k z] Sin[k z], {z, 0, 2 Pi/k}] -> 0,
  Integrate[Sin[k z] Cos[k z], {z, 0, 2 Pi/k}] -> 0
};

reducedLagrangianRaw = FullSimplify[
  (k/(2 Pi)) Integrate[lagrangianDensity, {z, 0, 2 Pi/k}],
  Assumptions -> {k > 0, aeta > 0, Mpl > 0, m > 0}
];

cleanupRules = {
  Derivative[1][a[eta]] :> Derivative[1][a][eta],
  Derivative[1][Derivative[1][a][eta]] :> Derivative[2][a][eta]
};

reducedLagrangian = Expand[reducedLagrangianRaw /. cleanupRules];

prime[expr_] := D[expr, eta];
EulerEquation[L_, q_Symbol] := FullSimplify[D[D[L, q'[eta]], eta] - D[L, q[eta]]];

backgroundRules = {
  Derivative[1][a][s_] :> a[s]^2 H[s],
  Derivative[2][a][s_] :> 2 a[s]^3 H[s]^2 + a[s]^2 Derivative[1][H][s],
  Derivative[1][phib][s_] :> -2 a[s] H[s] phib[s] - a[s]^2 Vp[s],
  Derivative[1][H][s_] :> -phib[s]^2/(2 a[s] Mpl^2),
  Derivative[1][Vp][s_] :> Vpp[s] phib[s]
};

beta = phib[eta]/(Mpl a[eta] H[eta]);
hatRule = pv -> Function[{s}, ph[s] + phib[s] AA[s]/(Mpl a[s] H[s])];
lagHat = FullSimplify[
  Expand[reducedLagrangian /. hatRule] //. backgroundRules,
  Assumptions -> {k > 0, a[eta] > 0, Mpl > 0, m > 0}
];

eqE = Numerator[Together[FullSimplify[EulerEquation[lagHat, EE]]]];
eqF = Numerator[Together[FullSimplify[EulerEquation[lagHat, FF]]]];
eqA0 = Numerator[Together[Expand[D[D[lagHat, AA'[eta]], eta] - D[lagHat, AA[eta]]]]];
eqB0 = Numerator[Together[Expand[D[D[lagHat, BB'[eta]], eta] - D[lagHat, BB[eta]]]]];
eqPh0 = Numerator[Together[Expand[D[D[lagHat, ph'[eta]], eta] - D[lagHat, ph[eta]]]]];

solEF = FullSimplify[
  Solve[{eqE == 0, eqF == 0}, {EE[eta], FF[eta]}] //. backgroundRules,
  Assumptions -> {k > 0, a[eta] > 0, Mpl > 0, m > 0}
];

Print["Lag size: ", LeafCount[lagHat]];
Print["eqE leaf count: ", LeafCount[eqE]];
Print["eqF leaf count: ", LeafCount[eqF]];
Print["Solutions for E,F found: ", Length[solEF]];
eqAeff = Numerator[Together[Expand[(eqA0 /. solEF[[1]]) //. backgroundRules]]];
Print["eqAeff leaf count: ", LeafCount[eqAeff]];
Print["Coeff[AA'']: ", Coefficient[eqAeff, AA''[eta]]];
Print["Coeff[AA']: ", Short[Coefficient[eqAeff, AA'[eta]], 12]];
solA = Solve[eqAeff == 0, AA[eta]] //. backgroundRules;
Print["Solutions for A found: ", Length[solA]];
If[Length[solA] > 0,
  Print["solA contains AA': ", ! FreeQ[solA[[1, 1, 2]], AA'[eta]]];
  Print["solA contains AA'': ", ! FreeQ[solA[[1, 1, 2]], AA''[eta]]];
];

removeSecondDerivativeTerm[expr_, q_] := Module[
  {coeff = Coefficient[expr, q''[eta]]},
  Expand[expr - coeff q''[eta] - D[coeff, eta] q'[eta]]
];

lagAfterEFRaw = Expand[(lagHat /. solEF[[1]]) //. backgroundRules];
lagAfterEFClean = Nest[removeSecondDerivativeTerm[#, AA] &, lagAfterEFRaw, 3];
Print["lagAfterEFRaw size: ", LeafCount[lagAfterEFRaw]];
Print["lagAfterEFClean size: ", LeafCount[lagAfterEFClean]];
Print["Coeff lag clean [AA'']: ", Coefficient[lagAfterEFClean, AA''[eta]]];
Print["Coeff lag clean [AA'^2]: ", Short[FullSimplify[Coefficient[lagAfterEFClean, AA'[eta]^2]], 20]];
eqAfromClean = Numerator[Together[Expand[D[D[lagAfterEFClean, AA'[eta]], eta] - D[lagAfterEFClean, AA[eta]]]]];
Print["Coeff eqAclean [AA'']: ", Coefficient[eqAfromClean, AA''[eta]]];
solAclean = Solve[eqAfromClean == 0, AA[eta]] //. backgroundRules;
Print["Solutions for A from clean lagrangian: ", Length[solAclean]];
If[Length[solAclean] > 0,
  Print["solAclean contains AA': ", ! FreeQ[solAclean[[1, 1, 2]], AA'[eta]]];
  Print["solAclean contains AA'': ", ! FreeQ[solAclean[[1, 1, 2]], AA''[eta]]];
];

eqBeff = Numerator[Together[Expand[(((eqB0 /. solEF[[1]]) /. solA[[1]])) //. backgroundRules]]];
eqPheff = Numerator[Together[Expand[(((eqPh0 /. solEF[[1]]) /. solA[[1]])) //. backgroundRules]]];
Print["eqBeff leaf count: ", LeafCount[eqBeff]];
Print["eqPheff leaf count: ", LeafCount[eqPheff]];
