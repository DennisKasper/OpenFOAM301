dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties
);

dimensionedScalar rho
(
    "rho",
    dimDensity,
    transportProperties
);

dimensionedScalar sigma
(
    "sigma",
    dimless,
    0.66666
);

dimensionedScalar kappa
(
    "kappa",
    dimless,
    0.41
);

dimensionedScalar cb1
(
    "cb1",
    dimless,
    0.1355
);

dimensionedScalar cb2
(
    "cb2",
    dimless,
    0.622
);

dimensionedScalar cs
(
    "cs",
    dimless,
    0.3
);

dimensionedScalar cv1
(
    "cv1",
    dimless,
    7.1
);

dimensionedScalar cw1 = SpalartAllmaras::cw1(cb1, cb2, kappa, sigma);

dimensionedScalar cw2
(
    "cw2",
    dimless,
    0.3
);

dimensionedScalar cw3
(
    "cw3",
    dimless,
    2.0
);