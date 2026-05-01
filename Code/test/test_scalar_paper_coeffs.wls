#!/usr/bin/env wolframscript

ClearAll["Global`*"];

Block[{Print = (Null &)},
  Get["/root/Agents/Spin2CGPP/Code/script/derive_scalar_reduced_action.wls"]
];

If[! TrueQ[aKineticResidual === 0],
  Print["FAIL: A kinetic term did not cancel after action-level E/F elimination."];
  Exit[1];
];

If[! And @@ Values[Map[TrueQ[# === 0] &, scalarCoefficientResiduals]],
  Print["FAIL: reduced scalar action coefficients do not match the paper."];
  Print[Map[TrueQ[# === 0] &, scalarCoefficientResiduals]];
  Exit[1];
];

Print["PASS: reduced scalar action coefficients match the published scalar-sector basis."];
