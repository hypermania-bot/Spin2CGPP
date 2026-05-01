#!/usr/bin/env wolframscript

ClearAll["Global`*"];

root = DirectoryName[DirectoryName[$InputFileName]];
Get[FileNameJoin[{root, "script", "paper_scalar_coeffs.wl"}]];

Block[{Print = (Null &)},
  Get[FileNameJoin[{root, "script", "derive_scalar_mode_eoms.wls"}]]
];

coeffs = PaperScalarCoeffs`MinimalScalarKineticCoefficients[eta];

accelCurrent = {
  {
    Coefficient[finalP, php''[eta]],
    Coefficient[finalP, Bp''[eta]]
  },
  {
    Coefficient[finalB, php''[eta]],
    Coefficient[finalB, Bp''[eta]]
  }
};

scaleP = FullSimplify[accelCurrent[[1, 1]]/coeffs["Kphi"]];
scaleB = FullSimplify[accelCurrent[[2, 2]]/coeffs["KB"]];

sampleRules = {
  a[eta] -> 2,
  H[eta] -> 3/10,
  phib[eta] -> 1/5,
  Mpl -> 5,
  m -> 7/10,
  k -> 11/10,
  Vp[eta] -> 13/10,
  Vpp[eta] -> 17/10,
  Derivative[1][a][eta] -> 12/10,
  Derivative[1][H][eta] -> -1/750,
  Derivative[1][phib][eta] -> -16/25,
  Derivative[1][Vp][eta] -> 17/250
};

tol = 10^-8;

residuals = {
  N[(accelCurrent[[1, 1]] - scaleP coeffs["Kphi"]) /. sampleRules],
  N[(accelCurrent[[2, 2]] - scaleB coeffs["KB"]) /. sampleRules]
};

If[! VectorQ[residuals, NumericQ],
  Print["FAIL: scalar paper-kinetic residuals did not evaluate numerically."];
  Exit[1];
];

If[Max[Abs[residuals]] > tol,
  Print["FAIL: diagonal acceleration coefficients do not match the published Kphi/KB formulas."];
  Exit[1];
];

Print["PASS: diagonal acceleration coefficients match the published Kphi and KB formulas after normalization."];
