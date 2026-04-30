#!/usr/bin/env wolframscript

Get["/root/Agents/Spin2CGPP/Code/script/minimal_scalar_covariant_complex.wls"];

mat = hessian /. {
  phib[eta] -> 0,
  Vp[eta] -> 0,
  Vpp[eta] -> 0
};

mat = FullSimplify[
  mat,
  Assumptions -> {a[eta] > 0, H[eta] > 0, Mpl > 0, m > 0, k > 0}
];

Print["de Sitter scalar Hessian = ", mat];
Print["de Sitter scalar determinant = ", FullSimplify[Det[mat], Assumptions -> {a[eta] > 0, H[eta] > 0, Mpl > 0, m > 0, k > 0}]];
Print["de Sitter sample rank = ", MatrixRank[N[mat /. {a[eta] -> 2, H[eta] -> 3/10, Mpl -> 5, m -> 7/10, k -> 11/10}]]];
