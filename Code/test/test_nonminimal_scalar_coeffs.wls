#!/usr/bin/env wolframscript

scriptPath = FileNameJoin[{DirectoryName[DirectoryName[$InputFileName]], "script", "derive_nonminimal_scalar_reduced_action.wls"}];

exitCode = RunProcess[{"wolframscript", "-file", scriptPath}, "ExitCode"];
If[exitCode =!= 0,
  Print["FAIL: nonminimal scalar covariant coefficient check failed."];
  Exit[1];
];

Exit[0];
