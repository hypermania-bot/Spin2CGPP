#!/usr/bin/env wolframscript

ClearAll["Global`*"];

$Assumptions = k > 0 && m > 0 && Mpl > 0;

eta = \[Eta];
z = zc;
eps = \[Epsilon];
aeta = a[eta];
mode = Cos[k z];

A0 = AA[eta] mode;
B0 = BB[eta] mode;
E0 = EE[eta] mode;
F0 = FF[eta] mode;
ph0 = pv[eta] mode;

bg = aeta^2 DiagonalMatrix[{-1, 1, 1, 1}];
h = ConstantArray[0, {4, 4}];
h[[1, 1]] = aeta^2 E0;
h[[1, 4]] = aeta^2 D[F0, z];
h[[4, 1]] = h[[1, 4]];
h[[2, 2]] = aeta^2 A0;
h[[3, 3]] = aeta^2 A0;
h[[4, 4]] = aeta^2 (A0 + D[D[B0, z], z]);

(* The paper's quadratic action is written in terms of the canonically normalized
   spin-2 field, so the metric perturbation enters as delta g = 2 h / Mpl. *)
g = bg + (2 eps/Mpl) h;
ginv = Normal[Series[Inverse[g], {eps, 0, 2}]] // Expand;
detg = Normal[Series[Det[g], {eps, 0, 2}]] // Expand;
sqrtg = Normal[Series[Sqrt[-detg], {eps, 0, 2}]] // Expand;

coords = {eta, x, y, z};

ClearAll[Conn];
Conn[mu_, nu_, rho_] := Conn[mu, nu, rho] = Expand[
  1/2 Sum[
    ginv[[mu, lam]] (D[g[[lam, nu]], coords[[rho]]] + D[g[[lam, rho]], coords[[nu]]] - D[g[[nu, rho]], coords[[lam]]]),
    {lam, 1, 4}
  ]
];

ClearAll[Ricci];
Ricci[mu_, nu_] := Ricci[mu, nu] = Expand[
  Sum[D[Conn[lam, mu, nu], coords[[lam]]], {lam, 1, 4}]
  - Sum[D[Conn[lam, mu, lam], coords[[nu]]], {lam, 1, 4}]
  + Sum[Conn[lam, lam, sig] Conn[sig, mu, nu], {lam, 1, 4}, {sig, 1, 4}]
  - Sum[Conn[lam, nu, sig] Conn[sig, mu, lam], {lam, 1, 4}, {sig, 1, 4}]
];

Rscalar = Expand[Sum[ginv[[mu, nu]] Ricci[mu, nu], {mu, 1, 4}, {nu, 1, 4}]];

phi = phib[eta] + eps ph0;
kinphi = Expand[-1/2 Sum[ginv[[mu, nu]] D[phi, coords[[mu]]] D[phi, coords[[nu]]], {mu, 1, 4}, {nu, 1, 4}]];
potphi = -(V0[eta] + eps Vp[eta] ph0 + eps^2/2 Vpp[eta] ph0^2);

Lfull = Expand[sqrtg (Mpl^2/2 Rscalar + kinphi + potphi)];
L2grav = Coefficient[Normal[Series[Lfull, {eps, 0, 2}]], eps, 2] // Expand;

bgInv = Inverse[bg];
traceh = Expand[Sum[bgInv[[mu, nu]] h[[mu, nu]], {mu, 1, 4}, {nu, 1, 4}]];
hUU = Table[
  Sum[bgInv[[mu, al]] bgInv[[nu, be]] h[[al, be]], {al, 1, 4}, {be, 1, 4}],
  {mu, 1, 4}, {nu, 1, 4}
] // Expand;

massTerm = Expand[(sqrtg /. eps -> 0) (-1/2 m^2) (Sum[hUU[[mu, nu]] h[[mu, nu]], {mu, 1, 4}, {nu, 1, 4}] - traceh^2)];
L2 = Expand[L2grav + massTerm];

reduced = FullSimplify[(k/(2 Pi)) Integrate[L2, {z, 0, 2 Pi/k}], Assumptions -> $Assumptions && aeta > 0] // Expand;

cleanupRules = {
  Derivative[1][a[eta]] :> Derivative[1][a][eta],
  Derivative[1][Derivative[1][a][eta]] :> Derivative[2][a][eta],
  Derivative[1][phib[eta]] :> Derivative[1][phib][eta],
  Derivative[1][Vp[eta]] :> Derivative[1][Vp][eta]
};
reduced = reduced /. cleanupRules;

backgroundRules = {
  V0[s_] :> 3 Mpl^2 H[s]^2 - phib[s]^2/(2 a[s]^2),
  Derivative[1][a][s_] :> a[s]^2 H[s],
  Derivative[2][a][s_] :> 2 a[s]^3 H[s]^2 + a[s]^2 Derivative[1][H][s],
  Derivative[1][phib][s_] :> -2 a[s] H[s] phib[s] - a[s]^2 Vp[s],
  Derivative[1][H][s_] :> -phib[s]^2/(2 a[s] Mpl^2),
  Derivative[1][Vp][s_] :> Vpp[s] phib[s]
};
reduced = Expand[reduced //. backgroundRules];

lagHat = Expand[
  reduced /. {
    pv -> Function[{s}, ph[s] + phib[s] AA[s]/(Mpl a[s] H[s])]
  }
];
lagHat = Expand[lagHat //. backgroundRules];

ClearAll[EulerEq];
EulerEq[L_, q_Symbol] := Numerator[Together[Expand[D[D[L, q'[eta]], eta] - D[L, q[eta]]]]];

eqE = EulerEq[lagHat, EE];
eqF = EulerEq[lagHat, FF];
solEF = Solve[{eqE == 0, eqF == 0}, {EE[eta], FF[eta]}];

lagAfterEF = Expand[lagHat /. solEF[[1]] //. backgroundRules];
eqA = EulerEq[lagAfterEF, AA] //. backgroundRules;

Print["reduced leaf count: ", LeafCount[reduced]];
Print["solEF count: ", Length[solEF]];
Print["lagAfterEF contains AA'^2: ", ! FreeQ[lagAfterEF, AA'[eta]^2]];
Print["eqA contains AA'': ", ! FreeQ[eqA, AA''[eta]]];
