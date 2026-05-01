#!/usr/bin/env wolframscript

Get[FileNameJoin[{DirectoryName[DirectoryName[$InputFileName]], "script", "Spin2CGPP.wl"}]];

ClearAll[zeroQ];
zeroQ[expr_] := TrueQ[FullSimplify[expr] === 0];

results = {
  VerificationTest[
    And @@ (zeroQ /@ Values[Spin2CGPP`CheckMinimalTensor[\[Eta]]]),
    True,
    TestID -> "minimal-tensor"
  ],
  VerificationTest[
    And @@ (zeroQ /@ Values[Spin2CGPP`CheckMinimalVector[\[Eta]]]),
    True,
    TestID -> "minimal-vector"
  ],
  VerificationTest[
    And @@ (zeroQ /@ Values[Spin2CGPP`CheckMinimalScalarCanonical[\[Eta]]]),
    True,
    TestID -> "minimal-scalar-canonical"
  ],
  VerificationTest[
    And @@ (zeroQ /@ Values[Spin2CGPP`CheckNonminimalTensor[\[Eta]]]),
    True,
    TestID -> "nonminimal-tensor"
  ],
  VerificationTest[
    And @@ (zeroQ /@ Values[Spin2CGPP`CheckNonminimalVector[\[Eta]]]),
    True,
    TestID -> "nonminimal-vector"
  ],
  VerificationTest[
    And @@ (zeroQ /@ Values[Spin2CGPP`CheckNonminimalScalar[\[Eta]]]),
    True,
    TestID -> "nonminimal-scalar"
  ],
  VerificationTest[
    RunProcess[
      {
        "wolframscript",
        "-file",
        FileNameJoin[{DirectoryName[$InputFileName], "test_scalar_mode_chain.wls"}]
      },
      "ExitCode"
    ],
    0,
    TestID -> "minimal-scalar-covariant-chain"
  ],
  VerificationTest[
    RunProcess[
      {
        "wolframscript",
        "-file",
        FileNameJoin[{DirectoryName[$InputFileName], "test_scalar_paper_kinetics.wls"}]
      },
      "ExitCode"
    ],
    0,
    TestID -> "minimal-scalar-paper-kinetics"
  ]
};

report = TestReport[results];
Print[report];
If[report["TestsFailed"] > 0, Exit[1], Exit[0]];
