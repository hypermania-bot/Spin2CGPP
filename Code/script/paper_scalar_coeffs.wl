BeginPackage["PaperScalarCoeffs`"];

MinimalScalarKineticCoefficients::usage =
  "MinimalScalarKineticCoefficients[eta] returns an association with the \
published minimally coupled scalar-sector coefficients \
<|\"mHsq\",\"Den\",\"Kphi\",\"Mphi\",\"KB\",\"MB\",\"L2\",\"L1\",\"L0\",\
\"Kappa\",\"KBcal\"|>.";

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
    vpp = Global`Vpp,
    mHsq, den, kphi, mphi, kB, mB, l2, l1, l0, kappa, kBcal,
    cp10, cp8, cp6, cp4, cp2, cp0, cb10, cb8, cb6, cb4,
    cl10, cl8, cl6, cl4, cl2
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

  cp10 = hh[eta]^4;
  cp8 =
    (1/2) aa[eta]^2 hh[eta]^2 (
      12 mm^2 hh[eta]^2 + 8 hh[eta]^4 - 14 hh[eta]^2 mHsq - mHsq^2
      + 4 hh[eta] vpl[eta] pb[eta]/(aa[eta] mpl^2)
      + 2 hh[eta]^2 vpp[eta]
    );
  cp6 =
    (3/8) aa[eta]^4 hh[eta]^2 (
      36 mm^4 hh[eta]^2 + 72 mm^2 hh[eta]^4 - 82 mm^2 hh[eta]^2 mHsq
      - 64 hh[eta]^4 mHsq - 7 mm^2 mHsq^2 + 40 hh[eta]^2 mHsq^2
      + 8 mHsq^3
      + 8 (3 mm^2 - 4 mHsq) hh[eta] vpl[eta] pb[eta]/(aa[eta] mpl^2)
      + 16 (mm^2 - mHsq) hh[eta]^2 vpp[eta]
    );
  cp4 =
    (3/8) aa[eta]^6 (
      4 hh[eta]^2 (
        9 mm^6 hh[eta]^2 + 36 mm^4 hh[eta]^4 + 16 mm^2 hh[eta]^6
        - 30 mm^4 hh[eta]^2 mHsq - 76 mm^2 hh[eta]^4 mHsq
        - 3 mm^4 mHsq^2 + 31 mm^2 hh[eta]^2 mHsq^2
        + 24 hh[eta]^4 mHsq^2 + 6 mm^2 mHsq^3
        - 6 hh[eta]^2 mHsq^3 - 3 mHsq^4
      )
      - 4 mm^2 hh[eta]^2 (hh[eta]^2 - mHsq) vpl[eta]^2/mpl^2
      + (
        36 mm^4 hh[eta]^2 + 8 mm^2 hh[eta]^4
        - 94 mm^2 hh[eta]^2 mHsq + mm^2 mHsq^2
        + 48 hh[eta]^2 mHsq^2
      ) hh[eta] vpl[eta] pb[eta]/(aa[eta] mpl^2)
      + (
        36 mm^4 hh[eta]^2 - 58 mm^2 hh[eta]^2 mHsq
        - mm^2 mHsq^2 + 24 hh[eta]^2 mHsq^2
      ) hh[eta]^2 vpp[eta]
    );
  cp2 =
    (9/32) aa[eta]^8 mm^2 (
      hh[eta]^2 (
        18 mm^6 hh[eta]^2 + 120 mm^4 hh[eta]^4 + 128 mm^2 hh[eta]^6
        - 78 mm^4 hh[eta]^2 mHsq - 384 mm^2 hh[eta]^4 mHsq
        - 9 mm^4 mHsq^2 + 132 mm^2 hh[eta]^2 mHsq^2
        + 128 hh[eta]^4 mHsq^2 + 23 mm^2 mHsq^3
        - 32 hh[eta]^2 mHsq^3 - 16 mHsq^4
      )
      - 8 hh[eta]^2 (2 mm^2 hh[eta]^2 - 2 mm^2 mHsq + mHsq^2) vpl[eta]^2/mpl^2
      + 4 (
        6 mm^4 hh[eta]^2 - 22 mm^2 hh[eta]^2 mHsq
        + mm^2 mHsq^2 + 14 hh[eta]^2 mHsq^2
      ) hh[eta] vpl[eta] pb[eta]/(aa[eta] mpl^2)
      + 4 (mm^2 - mHsq) (12 mm^2 hh[eta]^2 - 10 hh[eta]^2 mHsq - mHsq^2)
        hh[eta]^2 vpp[eta]
    );
  cp0 =
    (27/32) aa[eta]^10 mm^4 (
      -2 hh[eta]^2 (2 mm^2 hh[eta]^2 - 2 mm^2 mHsq + mHsq^2) vpl[eta]^2/mpl^2
      - mm^2 (2 hh[eta]^2 - mHsq) (4 hh[eta]^2 + mHsq)
        hh[eta] vpl[eta] pb[eta]/(aa[eta] mpl^2)
      + (mm^2 - mHsq) (6 mm^2 hh[eta]^2 - 4 hh[eta]^2 mHsq - mHsq^2)
        hh[eta]^2 vpp[eta]
    );
  mphi =
    (aa[eta]^2/2) (cp10 kk^10 + cp8 kk^8 + cp6 kk^6 + cp4 kk^4 + cp2 kk^2 + cp0)/den^2;

  cb10 = hh[eta]^2 (8 mm^2 hh[eta]^2 - 8 hh[eta]^4 - 2 hh[eta]^2 mHsq - mm^2 mHsq);
  cb8 =
    aa[eta]^2 hh[eta]^2 (
      30 mm^4 hh[eta]^2 + 32 mm^2 hh[eta]^4 - 96 hh[eta]^6
      - 3 mm^4 mHsq - 56 mm^2 hh[eta]^2 mHsq + 48 hh[eta]^4 mHsq
      + 5 mm^2 mHsq^2 + 6 hh[eta]^2 mHsq^2
      + (4 mm^2 - 24 hh[eta]^2) hh[eta] vpl[eta] pb[eta]/(aa[eta] mpl^2)
    );
  cb6 =
    (3/8) aa[eta]^4 mm^2 (
      96 mm^4 hh[eta]^4 + 144 mm^2 hh[eta]^6 - 6 mm^4 hh[eta]^2 mHsq
      - 252 mm^2 hh[eta]^4 mHsq - 192 hh[eta]^6 mHsq
      + 8 mm^2 hh[eta]^2 mHsq^2 + 200 hh[eta]^4 mHsq^2
      - 10 hh[eta]^2 mHsq^3 - mm^2 mHsq^3
      + (8 mm^2 mHsq - 16 hh[eta]^2 mHsq)
        hh[eta] vpl[eta] pb[eta]/(aa[eta] mpl^2)
    );
  cb4 =
    (3/8) aa[eta]^6 mm^4 (
      36 mm^4 hh[eta]^4 - 48 mm^2 hh[eta]^6 + 64 hh[eta]^8
      - 12 mm^2 hh[eta]^4 mHsq - 32 hh[eta]^6 mHsq
      - 12 mm^2 hh[eta]^2 mHsq^2 + 4 hh[eta]^4 mHsq^2
      + 12 hh[eta]^2 mHsq^3 - 3 mm^2 mHsq^3 + 2 mHsq^4
      - (24 mm^2 hh[eta]^2 - 16 hh[eta]^4 - 12 mm^2 mHsq
        - 8 hh[eta]^2 mHsq + 8 mHsq^2)
        hh[eta] vpl[eta] pb[eta]/(aa[eta] mpl^2)
    );
  mB =
    (aa[eta]^6 mm^2/8) (cb10 kk^10 + cb8 kk^8 + cb6 kk^6 + cb4 kk^4)/den^2;

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

  cl10 = hh[eta]^4;
  cl8 =
    (1/2) aa[eta]^2 hh[eta]^4 (
      9 mm^2 + 12 hh[eta]^2 - 13 mHsq
      - 4 aa[eta] hh[eta] vpl[eta]/pb[eta]
    );
  cl6 =
    (3/8) aa[eta]^4 hh[eta]^2 (
      18 mm^4 hh[eta]^2 + 32 mm^2 hh[eta]^4 + 64 hh[eta]^6
      - 48 mm^2 hh[eta]^2 mHsq - 64 hh[eta]^4 mHsq
      + mm^2 mHsq^2 + 28 hh[eta]^2 mHsq^2
      + 8 (-4 mm^2 hh[eta]^2 + 4 hh[eta]^4 + mm^2 mHsq)
        aa[eta] hh[eta] vpl[eta]/pb[eta]
    );
  cl4 =
    (3/16) aa[eta]^6 mm^2 hh[eta]^2 (
      18 mm^4 hh[eta]^2 - 24 mm^2 hh[eta]^4 + 256 hh[eta]^6
      - 54 mm^2 hh[eta]^2 mHsq - 160 hh[eta]^4 mHsq
      + 9 mm^2 mHsq^2 + 60 hh[eta]^2 mHsq^2 - 7 mHsq^3
      + 4 (-30 mm^2 hh[eta]^2 + 32 hh[eta]^4 + 12 mm^2 mHsq
        + 4 hh[eta]^2 mHsq - 7 mHsq^2)
        aa[eta] hh[eta] vpl[eta]/pb[eta]
    );
  cl2 =
    (9/16) aa[eta]^8 mm^4 hh[eta]^2 (2 hh[eta]^2 - mHsq) (
      -(4 hh[eta]^2 + mHsq) (3 mm^2 - 4 hh[eta]^2 - mHsq)
      + 4 (-3 mm^2 + 2 hh[eta]^2 + 2 mHsq)
        aa[eta] hh[eta] vpl[eta]/pb[eta]
    );
  l0 =
    (aa[eta]^3 mm^2 pb[eta]/(2 mpl hh[eta]))
    (cl10 kk^10 + cl8 kk^8 + cl6 kk^6 + cl4 kk^4 + cl2 kk^2)/den^2;

  kappa = -l2/(2 kk^2 kphi);
  kBcal = (4 kphi kB - l2^2)/(4 kk^4 kphi);

  <|
    "mHsq" -> mHsq,
    "Den" -> den,
    "Kphi" -> kphi,
    "Mphi" -> mphi,
    "KB" -> kB,
    "MB" -> mB,
    "L2" -> l2,
    "L1" -> l1,
    "L0" -> l0,
    "Kappa" -> kappa,
    "KBcal" -> kBcal
  |>
];

End[];
EndPackage[];
