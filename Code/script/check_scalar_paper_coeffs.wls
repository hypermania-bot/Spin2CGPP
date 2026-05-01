#!/usr/bin/env wolframscript

ClearAll["Global`*"];

Block[{Print = (Null &)},
  Get["/root/Agents/Spin2CGPP/Code/script/derive_scalar_reduced_action.wls"]
];

zeroFlags = Map[TrueQ[# === 0] &, scalarCoefficientResiduals];

Print["A kinetic residual zero: ", TrueQ[aKineticResidual === 0]];
Print["coefficient residuals zero: ", zeroFlags];

sampleRules = {
  a[eta] -> 2,
  H[eta] -> 3/10,
  phib[eta] -> 1/5,
  Mpl -> 5,
  m -> 7/10,
  k -> 11/10,
  Vp[eta] -> 13/10,
  Vpp[eta] -> 17/10
};

Print[
  "sample coefficients {ours, paper}: ",
  AssociationMap[
    N[{scalarReducedCoefficients[#], paperCoefficients[#]} /. sampleRules] &,
    {"Kphi", "Mphi", "KB", "MB", "L2", "L1", "L0"}
  ]
];

If[! TrueQ[aKineticResidual === 0] || ! And @@ Values[zeroFlags],
  Exit[1],
  Exit[0]
];
