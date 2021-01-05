namespace NavierStokes
{

// Defining the derivative of the Navier-Stokes coupling term
volScalarField calcDNuTdNuTilda
(
    const volScalarField& nuTilda,
    const dimensionedScalar nu,
    const dimensionedScalar cv1
)
{
    volScalarField chi = SpalartAllmaras::chi(nuTilda, nu);
    volScalarField fv1 = SpalartAllmaras::fv1(chi, cv1);

    return 3*pow3(cv1)*pow3(chi)/sqr(pow3(cv1) + pow3(chi)) + fv1;
}

} // End namespace NavierStokes