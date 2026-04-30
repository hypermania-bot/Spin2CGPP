#!/usr/bin/env wolframscript

ClearAll["Global`*"];

scriptPath = "/root/Agents/Spin2CGPP/Code/script/scalar_equation_rank_check.wls";

Block[{Print = (Null &)}, Get[scriptPath]];

ClearAll[swapPM];
swapPM[expr_] := expr /. {
  Ap -> Am, Am -> Ap,
  Bp -> Bm, Bm -> Bp,
  php -> phm, phm -> php,
  Ep -> Em, Em -> Ep,
  Fp -> Fm, Fm -> Fp
};

subF = Join[
  Solve[eqF == 0, Fp[eta]][[1]],
  Solve[eqFp == 0, Fm[eta]][[1]]
];

eE2 = normalizeWeirdDerivatives[Expand[(eqE /. subF) //. backgroundRules]];
eEp2 = normalizeWeirdDerivatives[Expand[(eqEp /. subF) //. backgroundRules]];
eA2 = normalizeWeirdDerivatives[Expand[(eqA /. subF) //. backgroundRules]];
eB2 = normalizeWeirdDerivatives[Expand[(eqB /. subF) //. backgroundRules]];
eP2 = normalizeWeirdDerivatives[Expand[(eqP /. subF) //. backgroundRules]];

subA = Solve[{eE2 == 0, eEp2 == 0}, {Ap[eta], Am[eta]}][[1]];
eA3 = normalizeWeirdDerivatives[Expand[(eA2 /. subA) //. backgroundRules]];
eAm3 = normalizeWeirdDerivatives[Expand[(swapPM[eqA] /. subF /. subA) //. backgroundRules]];
subE = Solve[{eA3 == 0, eAm3 == 0}, {Ep[eta], Em[eta]}][[1]];

finalB = normalizeWeirdDerivatives[Expand[(eB2 /. subA /. subE) //. backgroundRules]];
finalP = normalizeWeirdDerivatives[Expand[(eP2 /. subA /. subE) //. backgroundRules]];

accelMat = Table[
  Coefficient[{finalB, finalP}[[i]], {Bp''[eta], php''[eta]}[[j]]],
  {i, 1, 2}, {j, 1, 2}
] // Expand;

sampleRules = {
  a[eta] -> 2, H[eta] -> 3/10, phib[eta] -> 1/5, Mpl -> 5, m -> 7/10, k -> 11/10,
  Vp[eta] -> 13/10, Vpp[eta] -> 17/10,
  Bp[eta] -> 19/10, Bm[eta] -> 23/10, php[eta] -> 29/10, phm[eta] -> 31/10,
  Bp'[eta] -> 37/10, Bm'[eta] -> 41/10, php'[eta] -> 43/10, phm'[eta] -> 47/10,
  Bp''[eta] -> 53/10, php''[eta] -> 61/10
};

If[FreeQ[finalB, Bp''[eta]] || FreeQ[finalB, php''[eta]],
  Print["FAIL: finalB is not a second-order equation in both fields."];
  Exit[1];
];

If[FreeQ[finalP, Bp''[eta]] || FreeQ[finalP, php''[eta]],
  Print["FAIL: finalP is not a second-order equation in both fields."];
  Exit[1];
];

If[MatrixRank[N[accelMat /. sampleRules]] =!= 2,
  Print["FAIL: final two-field acceleration matrix does not have rank 2."];
  Exit[1];
];

Print["PASS: sequential constraint elimination yields a coupled two-field scalar system."];
