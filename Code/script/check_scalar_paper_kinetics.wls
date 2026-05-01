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

accelPaper = {
  {coeffs["Kphi"], coeffs["L2"]},
  {coeffs["L2"], coeffs["KB"]}
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

accelCurrentNum = N[accelCurrent /. sampleRules];
accelPaperNum = N[accelPaper /. sampleRules];
e1CurrentNum = N[
  {
    {
      Coefficient[finalP, php'[eta]],
      Coefficient[finalP, Bp'[eta]]
    },
    {
      Coefficient[finalB, php'[eta]],
      Coefficient[finalB, Bp'[eta]]
    }
  } /. sampleRules
];
e1PaperNum = N[
  {
    {
      D[coeffs["Kphi"], eta],
      D[coeffs["L2"], eta] - coeffs["L1"]
    },
    {
      D[coeffs["L2"], eta] + coeffs["L1"],
      D[coeffs["KB"], eta]
    }
  } /. sampleRules
];

lift = accelPaperNum . Inverse[accelCurrentNum];
accelResidual = Chop[lift . accelCurrentNum - accelPaperNum];
e1Residual = Chop[lift . e1CurrentNum - e1PaperNum];

Print["scaleP = ", N[scaleP /. sampleRules]];
Print["scaleB = ", N[scaleB /. sampleRules]];
Print["row-P Kphi residual = ", N[(accelCurrent[[1, 1]] - scaleP coeffs["Kphi"]) /. sampleRules]];
Print["row-B KB residual = ", N[(accelCurrent[[2, 2]] - scaleB coeffs["KB"]) /. sampleRules]];
Print["lifted acceleration residual = ", MatrixForm[accelResidual]];
Print["lifted first-derivative residual = ", MatrixForm[e1Residual]];
Print["kappa sample = ", N[coeffs["Kappa"] /. sampleRules]];
Print["KBcal sample = ", N[coeffs["KBcal"] /. sampleRules]];
