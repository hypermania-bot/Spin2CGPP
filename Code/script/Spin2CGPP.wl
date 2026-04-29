BeginPackage["Spin2CGPP`"];

CanonicalOmegaSquared::usage = "CanonicalOmegaSquared[K, M, eta] returns the canonical frequency squared for L = K q'^2 - M q^2.";
CanonicalLagrangian::usage = "CanonicalLagrangian[K, M, q, chi, eta] returns the canonically normalized single-field Lagrangian after q -> chi/Sqrt[2 K].";
EulerEquation::usage = "EulerEquation[L, q, eta] returns the Euler-Lagrange equation for q[eta].";
MinimalTensorFrequency::usage = "MinimalTensorFrequency[a, H, m, k, eta] returns the tensor-sector omega_k^2.";
MinimalVectorKinetic::usage = "MinimalVectorKinetic[a, m, k] returns the vector-sector kinetic coefficient.";
MinimalVectorFrequency::usage = "MinimalVectorFrequency[a, m, k, eta] returns the minimal vector omega_k^2.";
NonminimalTensorFrequency::usage = "NonminimalTensorFrequency[a, H, Hp, m, lambda, k] returns the nonminimal tensor omega_k^2.";
NonminimalVectorCoefficients::usage = "NonminimalVectorCoefficients[a, mu1, mu2, k] returns <|\"K\" -> ..., \"M\" -> ...|>.";
NonminimalVectorSoundSpeed::usage = "NonminimalVectorSoundSpeed[mu1, mu2] returns mu2^2/mu1^2.";
CoupledScalarEOMs::usage = "CoupledScalarEOMs[omegaPi, omegaB, sigma1, sigma0, eta] returns the two coupled scalar equations from eq. (3.52).";
CheckMinimalTensor::usage = "CheckMinimalTensor[eta] verifies the minimal tensor derivation.";
CheckMinimalVector::usage = "CheckMinimalVector[eta] verifies the minimal vector derivation.";
CheckMinimalScalarCanonical::usage = "CheckMinimalScalarCanonical[eta] verifies the coupled scalar Euler-Lagrange equations.";
CheckNonminimalTensor::usage = "CheckNonminimalTensor[eta] verifies the nonminimal tensor derivation.";
CheckNonminimalVector::usage = "CheckNonminimalVector[eta] verifies the nonminimal vector derivation.";
CheckNonminimalScalar::usage = "CheckNonminimalScalar[eta] verifies the nonminimal scalar derivation.";

Begin["`Private`"];

ClearAll[prime, realExpand];
prime[expr_, eta_] := D[expr, eta];
realExpand[expr_] := ComplexExpand[expr, TargetFunctions -> {Re, Im}];

ClearAll[CanonicalOmegaSquared];
CanonicalOmegaSquared[K_, M_, eta_] :=
  FullSimplify[(4 K M + prime[K, eta]^2 - 2 K prime[K, {eta, 2}])/(4 K^2)];

ClearAll[CanonicalLagrangian];
CanonicalLagrangian[K_, M_, q_Symbol, chi_Symbol, eta_] := Module[
  {lagrangian, transformed},
  lagrangian = K prime[q[eta], eta]^2 - M q[eta]^2;
  transformed = Expand[lagrangian /. q[eta] -> chi[eta]/Sqrt[2 K] /. Derivative[1][q][eta] -> prime[chi[eta]/Sqrt[2 K], eta]];
  transformed
];

ClearAll[EulerEquation];
EulerEquation[L_, q_Symbol, eta_] := FullSimplify[D[D[L, q'[eta]], eta] - D[L, q[eta]]];

ClearAll[MinimalTensorFrequency];
MinimalTensorFrequency[a_, H_, m_, k_, eta_] := k^2 + a[eta]^2 m^2 - 2 a[eta]^2 H[eta]^2 - a[eta] H'[eta];

ClearAll[MinimalVectorKinetic, MinimalVectorFrequency];
MinimalVectorKinetic[a_, m_, k_] := a^4 k^2 m^2/(k^2 + a^2 m^2);
MinimalVectorFrequency[a_, m_, k_, eta_] := Module[
  {f},
  f = a[eta]^2/Sqrt[k^2 + a[eta]^2 m^2];
  FullSimplify[k^2 + a[eta]^2 m^2 - prime[f, {eta, 2}]/f]
];

ClearAll[NonminimalTensorFrequency];
NonminimalTensorFrequency[a_, H_, Hp_, m_, lambda_, k_] := k^2 + a^2 m^2 - a^2 lambda + a^2 H^2 + a Hp;

ClearAll[NonminimalVectorCoefficients, NonminimalVectorSoundSpeed];
NonminimalVectorCoefficients[a_, mu1_, mu2_, k_] := <|
  "K" -> a^4 k^2 mu1^2/(k^2 + a^2 mu1^2),
  "M" -> a^4 k^2 mu2^2
|>;
NonminimalVectorSoundSpeed[mu1_, mu2_] := FullSimplify[mu2^2/mu1^2];

ClearAll[CoupledScalarEOMs];
CoupledScalarEOMs[omegaPi_, omegaB_, sigma1_, sigma0_, eta_] := Module[
  {chiPi, chiB, lagrangian, eq1, eq2},
  lagrangian =
    1/2 prime[chiPi[eta], eta]^2 - 1/2 omegaPi chiPi[eta]^2 +
    1/2 prime[chiB[eta], eta]^2 - 1/2 omegaB chiB[eta]^2 +
    sigma1 chiPi[eta] prime[chiB[eta], eta] -
    (sigma0 - prime[sigma1, eta]/2) chiPi[eta] chiB[eta];
  eq1 = EulerEquation[lagrangian, chiPi, eta];
  eq2 = EulerEquation[lagrangian, chiB, eta];
  <|"chiPi" -> eq1, "chiB" -> eq2|>
];

ClearAll[CheckMinimalTensor];
CheckMinimalTensor[eta_] := Module[
  {a, m, k, Df, chif, lagD, lagChi, eqChi, expected},
  lagD = 1/2 a[eta]^2 (prime[Df[eta], eta]^2 - (k^2 + a[eta]^2 m^2) Df[eta]^2);
  lagChi = FullSimplify[
    Expand[lagD /. Df[eta] -> chif[eta]/a[eta] /. Derivative[1][Df][eta] -> prime[chif[eta]/a[eta], eta]],
    Assumptions -> a[eta] > 0
  ];
  expected = prime[chif[eta], {eta, 2}] + (k^2 + a[eta]^2 m^2 - a''[eta]/a[eta]) chif[eta];
  eqChi = FullSimplify[EulerEquation[lagChi, chif, eta]];
  <|
    "EquationDifference" -> FullSimplify[eqChi - expected],
    "FRWFrequencyDifference" -> FullSimplify[(k^2 + a[eta]^2 m^2 - a''[eta]/a[eta]) - MinimalTensorFrequency[a, Function[x, a'[x]/a[x]^2], m, k, eta]]
  |>
];

ClearAll[CheckMinimalVector];
CheckMinimalVector[eta_] := Module[
  {a, m, k, C, chi, K, M, lagC, lagChi, omega, eqChi},
  K = MinimalVectorKinetic[a[eta], m, k];
  M = a[eta]^4 k^2 m^2;
  lagC = K prime[C[eta], eta]^2 - M C[eta]^2;
  omega = CanonicalOmegaSquared[K, M, eta];
  lagChi = FullSimplify[
    Expand[lagC /. C[eta] -> chi[eta]/Sqrt[2 K] /. Derivative[1][C][eta] -> prime[chi[eta]/Sqrt[2 K], eta]]
  ];
  eqChi = FullSimplify[EulerEquation[lagChi, chi, eta]];
  <|
    "EquationDifference" -> FullSimplify[eqChi - (chi''[eta] + omega chi[eta])],
    "FrequencyDifference" -> FullSimplify[omega - MinimalVectorFrequency[a, m, k, eta]]
  |>
];

ClearAll[CheckMinimalScalarCanonical];
CheckMinimalScalarCanonical[eta_] := Module[
  {omegaPi, omegaB, sigma1, sigma0, chiPi, chiB, lagrangian, eq1, eq2},
  lagrangian =
    1/2 prime[chiPi[eta], eta]^2 - 1/2 omegaPi[eta] chiPi[eta]^2 +
    1/2 prime[chiB[eta], eta]^2 - 1/2 omegaB[eta] chiB[eta]^2 +
    sigma1[eta] chiPi[eta] prime[chiB[eta], eta] -
    (sigma0[eta] - sigma1'[eta]/2) chiPi[eta] chiB[eta];
  eq1 = FullSimplify[EulerEquation[lagrangian, chiPi, eta]];
  eq2 = FullSimplify[EulerEquation[lagrangian, chiB, eta]];
  <|
    "chiPiDifference" -> FullSimplify[eq1 - (chiPi''[eta] + omegaPi[eta] chiPi[eta] - sigma1[eta] chiB'[eta] + (sigma0[eta] - sigma1'[eta]/2) chiB[eta])],
    "chiBDifference" -> FullSimplify[eq2 - (chiB''[eta] + omegaB[eta] chiB[eta] + sigma1[eta] chiPi'[eta] + (sigma0[eta] + sigma1'[eta]/2) chiPi[eta])]
  |>
];

ClearAll[CheckNonminimalTensor];
CheckNonminimalTensor[eta_] := Module[
  {a, h, hp, m, lambda, k, Df, chif, lagD, lagChi, omega, eqChi, relationRules},
  lagD = 1/2 a[eta]^2 (prime[Df[eta], eta]^2 - (k^2 + a[eta]^2 (m^2 - lambda + 3 h[eta]^2 + 2 hp[eta]/a[eta])) Df[eta]^2);
  lagChi = FullSimplify[
    Expand[lagD /. Df[eta] -> chif[eta]/a[eta] /. Derivative[1][Df][eta] -> prime[chif[eta]/a[eta], eta]]
  ];
  omega = k^2 + a[eta]^2 (m^2 - lambda + 3 h[eta]^2 + 2 hp[eta]/a[eta]) - a''[eta]/a[eta];
  eqChi = FullSimplify[EulerEquation[lagChi, chif, eta]];
  relationRules = {
    Derivative[1][a][eta] -> a[eta]^2 h[eta],
    Derivative[2][a][eta] -> 2 a[eta]^3 h[eta]^2 + a[eta]^2 hp[eta]
  };
  <|
    "EquationDifference" -> FullSimplify[eqChi - (chif''[eta] + omega chif[eta])],
    "FrequencyDifference" -> FullSimplify[(omega - NonminimalTensorFrequency[a[eta], h[eta], hp[eta], m, lambda, k]) /. relationRules]
  |>
];

ClearAll[CheckNonminimalVector];
CheckNonminimalVector[eta_] := Module[
  {a, mu1, mu2, k, C, chi, coeffs, K, M, lagC, lagChi, omega, eqChi},
  coeffs = NonminimalVectorCoefficients[a[eta], mu1[eta], mu2[eta], k];
  K = coeffs["K"];
  M = coeffs["M"];
  lagC = K prime[C[eta], eta]^2 - M C[eta]^2;
  omega = CanonicalOmegaSquared[K, M, eta];
  lagChi = FullSimplify[
    Expand[lagC /. C[eta] -> chi[eta]/Sqrt[2 K] /. Derivative[1][C][eta] -> prime[chi[eta]/Sqrt[2 K], eta]]
  ];
  eqChi = FullSimplify[EulerEquation[lagChi, chi, eta]];
  <|
    "EquationDifference" -> FullSimplify[eqChi - (chi''[eta] + omega chi[eta])],
    "SoundSpeedDifference" -> FullSimplify[NonminimalVectorSoundSpeed[mu1[eta], mu2[eta]] - mu2[eta]^2/mu1[eta]^2]
  |>
];

ClearAll[CheckNonminimalScalar];
CheckNonminimalScalar[eta_] := Module[
  {K, M, B, chi, lagB, lagChi, omega, eqChi},
  K = Kb[eta];
  M = Mb[eta];
  lagB = K prime[B[eta], eta]^2 - M B[eta]^2;
  omega = CanonicalOmegaSquared[K, M, eta];
  lagChi = FullSimplify[
    Expand[lagB /. B[eta] -> chi[eta]/Sqrt[2 K] /. Derivative[1][B][eta] -> prime[chi[eta]/Sqrt[2 K], eta]]
  ];
  eqChi = FullSimplify[EulerEquation[lagChi, chi, eta]];
  <|
    "EquationDifference" -> FullSimplify[eqChi - (chi''[eta] + omega chi[eta])]
  |>
];

End[];
EndPackage[];
