BeginPackage["PaperScalarCoeffs`"];

MinimalScalarKineticCoefficients::usage =
  "MinimalScalarKineticCoefficients[eta] returns an association with the \
published minimally coupled scalar-sector coefficients \
<|\"mHsq\",\"Den\",\"Kphi\",\"KB\",\"L2\",\"L1\",\"Kappa\",\"KBcal\"|>.";

Begin["`Private`"];

ClearAll[MinimalScalarKineticCoefficients];
MinimalScalarKineticCoefficients[eta_] := Module[
  {
    aa = Global`a,
    hh = Global`H,
    pb = Global`phib,
    vpl = Global`Vp,
    mpl = Global`Mpl,
    mm = Global`m,
    kk = Global`k,
    mHsq, den, kphi, kB, l2, l1, kappa, kBcal
  },
  mHsq = 2 hh[eta]^2 - pb[eta]^2/(aa[eta]^2 mpl^2);
  den =
    hh[eta]^2 kk^4
    + 3 aa[eta]^2 (mm^2 - mHsq) hh[eta]^2 kk^2
    + (3/8) aa[eta]^4 mm^2 (6 mm^2 hh[eta]^2 - 4 hh[eta]^2 mHsq - mHsq^2);

  kphi =
    (aa[eta]^2/2)
    (
      hh[eta]^2 kk^4
      + 3 aa[eta]^2 (mm^2 - mHsq) hh[eta]^2 kk^2
      + (9/4) aa[eta]^4 mm^2 (mm^2 - mHsq) hh[eta]^2
    )/den;

  kB =
    (aa[eta]^6 mm^2/8)
    ((8 mm^2 hh[eta]^2 - 6 hh[eta]^2 mHsq - mm^2 mHsq) kk^4)/den;

  l2 =
    (aa[eta]^3 mm^2 pb[eta]/(2 mpl hh[eta]))
    (
      hh[eta]^2 kk^4
      + (3/2) aa[eta]^2 (mm^2 - mHsq) hh[eta]^2 kk^2
    )/den;

  l1 =
    -(aa[eta]^4 mm^2 pb[eta]/mpl)
    (
      (
        hh[eta]^2
        - mHsq/4
        - (aa[eta] hh[eta] vpl[eta])/(2 pb[eta])
      ) kk^4
      - (3/2) aa[eta]^2 (mm^2 - mHsq)
        (
          hh[eta]^2
          + mHsq/4
          + (aa[eta] hh[eta] vpl[eta])/(2 pb[eta])
        ) kk^2
    )/den;

  kappa = -l2/(2 kk^2 kphi);
  kBcal = (4 kphi kB - l2^2)/(4 kk^4 kphi);

  <|
    "mHsq" -> mHsq,
    "Den" -> den,
    "Kphi" -> kphi,
    "KB" -> kB,
    "L2" -> l2,
    "L1" -> l1,
    "Kappa" -> kappa,
    "KBcal" -> kBcal
  |>
];

End[];
EndPackage[];
