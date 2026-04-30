#!/usr/bin/env wolframscript

ClearAll["Global`*"];

scriptPath = "/root/Agents/Spin2CGPP/Code/script/scalar_equation_rank_check.wls";

(* Load the covariant scalar reduction quietly; this defines eqE, eqF, eqEp, eqFp,
   eqA, eqB, eqP, backgroundRules, and normalizeWeirdDerivatives. *)
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

(* The E equations are constraints that solve A_\pm without introducing A_\pm''. *)
subA = Solve[{eE2 == 0, eEp2 == 0}, {Ap[eta], Am[eta]}][[1]];

eA3 = normalizeWeirdDerivatives[Expand[(eA2 /. subA) //. backgroundRules]];
eAm3 = normalizeWeirdDerivatives[Expand[(swapPM[eqA] /. subF /. subA) //. backgroundRules]];

(* The A equations then solve E_\pm. *)
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

Print["finalB contains Bp'': ", ! FreeQ[finalB, Bp''[eta]]];
Print["finalB contains php'': ", ! FreeQ[finalB, php''[eta]]];
Print["finalP contains Bp'': ", ! FreeQ[finalP, Bp''[eta]]];
Print["finalP contains php'': ", ! FreeQ[finalP, php''[eta]]];
Print["accelMat sample rank: ", MatrixRank[N[accelMat /. sampleRules]]];
Print["accelMat sample det: ", Det[N[accelMat /. sampleRules]]];
Print["finalB leaf count: ", LeafCount[finalB]];
Print["finalP leaf count: ", LeafCount[finalP]];
