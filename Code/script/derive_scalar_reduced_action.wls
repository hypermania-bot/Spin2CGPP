#!/usr/bin/env wolframscript

ClearAll["Global`*"];

Block[{Print = (Null &)},
  Get["/root/Agents/Spin2CGPP/Code/script/scalar_equation_rank_check.wls"]
];
Get["/root/Agents/Spin2CGPP/Code/script/paper_scalar_coeffs.wl"];

ClearAll[simplifyOnBackground, removeFirstDerivativeTerm, removeDerivatives];
simplifyOnBackground[expr_] :=
  Cancel[Together[normalizeWeirdDerivatives[expr //. backgroundRules]]];

removeFirstDerivativeTerm[expr_, q_] := Module[
  {coeff = simplifyOnBackground[Coefficient[expr, q'[eta]]]},
  normalizeWeirdDerivatives[Expand[expr - coeff q'[eta] - D[coeff, eta] q[eta]]]
];

removeDerivatives[expr_, fields_] := FixedPoint[
  Function[e, simplifyOnBackground[Fold[removeFirstDerivativeTerm, e, fields]]],
  expr,
  4
];

(* E and F have no quadratic kinetic terms, but they do appear through first
   derivatives. Move those derivatives by parts before solving the constraints. *)
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

(* After E,F elimination the A kinetic term cancels. Move the remaining A
   derivatives by parts, solve A algebraically, and leave the physical
   derivatives on B and phihat so the paper coefficient basis can be read off. *)
lagAAlgebraic = removeDerivatives[lagAfterAuxiliaries, {Ap, Am}];
aEquations = Table[D[lagAAlgebraic, q[eta]], {q, {Am, Ap}}] // Expand;
aSolution = Solve[Thread[aEquations == 0], {Ap[eta], Am[eta]}][[1]];

lagScalarReduced = simplifyOnBackground[Expand[lagAAlgebraic /. aSolution]];

phiPrimeB = simplifyOnBackground[Coefficient[lagScalarReduced, php'[eta] Bm[eta]]];
phiPrimeBMinus = simplifyOnBackground[Coefficient[lagScalarReduced, phm'[eta] Bp[eta]]];
phiPrimePhi = simplifyOnBackground[Coefficient[lagScalarReduced, php'[eta] phm[eta]]];
bPrimeB = simplifyOnBackground[Coefficient[lagScalarReduced, Bp'[eta] Bm[eta]]];

scalarReducedCoefficients = <|
  "Kphi" -> simplifyOnBackground[Coefficient[lagScalarReduced, php'[eta] phm'[eta]]/2],
  "Mphi" -> simplifyOnBackground[
    -(Coefficient[lagScalarReduced, php[eta] phm[eta]] - D[phiPrimePhi, eta])/2
  ],
  "KB" -> simplifyOnBackground[Coefficient[lagScalarReduced, Bp'[eta] Bm'[eta]]/2],
  "MB" -> simplifyOnBackground[
    -(Coefficient[lagScalarReduced, Bp[eta] Bm[eta]] - D[bPrimeB, eta])/2
  ],
  "L2" -> simplifyOnBackground[Coefficient[lagScalarReduced, php'[eta] Bm'[eta]]],
  "L2Minus" -> simplifyOnBackground[Coefficient[lagScalarReduced, phm'[eta] Bp'[eta]]],
  "L1" -> simplifyOnBackground[Coefficient[lagScalarReduced, php[eta] Bm'[eta]] - phiPrimeB],
  "L1Minus" -> simplifyOnBackground[Coefficient[lagScalarReduced, phm[eta] Bp'[eta]] - phiPrimeBMinus],
  "L0" -> simplifyOnBackground[
    -(Coefficient[lagScalarReduced, php[eta] Bm[eta]] - D[phiPrimeB, eta])
  ],
  "L0Minus" -> simplifyOnBackground[
    -(Coefficient[lagScalarReduced, phm[eta] Bp[eta]] - D[phiPrimeBMinus, eta])
  ]
|>;

paperCoefficients = PaperScalarCoeffs`MinimalScalarKineticCoefficients[eta];
scalarCoefficientResiduals = <|
  "Kphi" -> simplifyOnBackground[scalarReducedCoefficients["Kphi"] - paperCoefficients["Kphi"]],
  "Mphi" -> simplifyOnBackground[scalarReducedCoefficients["Mphi"] - paperCoefficients["Mphi"]],
  "KB" -> simplifyOnBackground[scalarReducedCoefficients["KB"] - paperCoefficients["KB"]],
  "MB" -> simplifyOnBackground[scalarReducedCoefficients["MB"] - paperCoefficients["MB"]],
  "L2" -> simplifyOnBackground[scalarReducedCoefficients["L2"] - paperCoefficients["L2"]],
  "L2Minus" -> simplifyOnBackground[scalarReducedCoefficients["L2Minus"] - paperCoefficients["L2"]],
  "L1" -> simplifyOnBackground[scalarReducedCoefficients["L1"] - paperCoefficients["L1"]],
  "L1Minus" -> simplifyOnBackground[scalarReducedCoefficients["L1Minus"] - paperCoefficients["L1"]],
  "L0" -> simplifyOnBackground[scalarReducedCoefficients["L0"] - paperCoefficients["L0"]],
  "L0Minus" -> simplifyOnBackground[scalarReducedCoefficients["L0Minus"] - paperCoefficients["L0"]]
|>;

If[$InputFileName =!= "",
  Print["auxiliary derivative flags: ", Table[! FreeQ[lagAuxiliaryAlgebraic, q'[eta]], {q, {Ep, Em, Fp, Fm}}]];
  Print["A kinetic residual zero: ", TrueQ[aKineticResidual === 0]];
  Print["reduced action leaf count: ", LeafCount[lagScalarReduced]];
  Print["coefficient residuals zero: ", Map[TrueQ[# === 0] &, scalarCoefficientResiduals]];
];
