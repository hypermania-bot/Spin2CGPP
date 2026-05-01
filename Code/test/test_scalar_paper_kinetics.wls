#!/usr/bin/env wolframscript

ClearAll["Global`*"];

Block[{Print = (Null &)},
  Get["/root/Agents/Spin2CGPP/Code/script/derive_scalar_reduced_action.wls"]
];

kineticResiduals = KeyTake[
  Map[TrueQ[# === 0] &, scalarCoefficientResiduals],
  {"Kphi", "KB", "L2", "L2Minus", "L1", "L1Minus"}
];

If[! TrueQ[aKineticResidual === 0] || ! And @@ Values[kineticResiduals],
  Print["FAIL: scalar kinetic coefficients do not match the published basis."];
  Print[kineticResiduals];
  Exit[1];
];

Print["PASS: scalar kinetic coefficients match the published basis."];
