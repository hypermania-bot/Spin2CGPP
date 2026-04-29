#!/usr/bin/env wolframscript

Get[FileNameJoin[{DirectoryName[$InputFileName], "Spin2CGPP.wl"}]];

ClearAll[reportAssociation];
reportAssociation[name_, assoc_] := Module[{},
  Print["== ", name, " =="];
  KeyValueMap[Print[#1, ": ", #2] &, assoc];
  Print[""];
];

reportAssociation["Minimal tensor", Spin2CGPP`CheckMinimalTensor[\[Eta]]];
reportAssociation["Minimal vector", Spin2CGPP`CheckMinimalVector[\[Eta]]];
reportAssociation["Minimal scalar canonical system", Spin2CGPP`CheckMinimalScalarCanonical[\[Eta]]];
reportAssociation["Nonminimal tensor", Spin2CGPP`CheckNonminimalTensor[\[Eta]]];
reportAssociation["Nonminimal vector", Spin2CGPP`CheckNonminimalVector[\[Eta]]];
reportAssociation["Nonminimal scalar", Spin2CGPP`CheckNonminimalScalar[\[Eta]]];
