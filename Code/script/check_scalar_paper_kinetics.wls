#!/usr/bin/env wolframscript

ClearAll["Global`*"];

Block[{Print = (Null &)},
  Get["/root/Agents/Spin2CGPP/Code/script/derive_scalar_reduced_action.wls"]
];

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

Print["A kinetic residual zero: ", TrueQ[aKineticResidual === 0]];
Print[
  "kinetic residuals zero: ",
  KeyTake[Map[TrueQ[# === 0] &, scalarCoefficientResiduals], {"Kphi", "KB", "L2", "L2Minus", "L1", "L1Minus"}]
];
Print[
  "sample kinetic coefficients {ours, paper}: ",
  AssociationMap[
    N[{scalarReducedCoefficients[#], paperCoefficients[#]} /. sampleRules] &,
    {"Kphi", "KB", "L2", "L1"}
  ]
];
Print["kappa sample = ", N[paperCoefficients["Kappa"] /. sampleRules]];
Print["KBcal sample = ", N[paperCoefficients["KBcal"] /. sampleRules]];
