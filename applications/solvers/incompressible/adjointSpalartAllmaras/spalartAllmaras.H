namespace SpalartAllmaras
{

// Defining Spalart-Allmaras functions
dimensionedScalar cw1
(
    const dimensionedScalar cb1,
    const dimensionedScalar cb2,
    const dimensionedScalar kappa,
    const dimensionedScalar sigma
)
{
    return (cb1/sqr(kappa) + (1 + cb2)/sigma);
};

tmp<volScalarField> chi
(
    const volScalarField& nuTilda,
    const dimensionedScalar nu
)
{
    return nuTilda/nu;
}


tmp<volScalarField> fv1
(
    const volScalarField& chi,
    const dimensionedScalar cv1
)
{
    const volScalarField chi3(pow3(chi));
 
    return chi3/(chi3 + pow3(cv1));
}


tmp<volScalarField> fv2
(
    const volScalarField& chi,
    const volScalarField& fv1
)
{
    return 1.0 - chi/(1.0 + chi*fv1);
}


tmp<volScalarField> Stilda
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volScalarField& nuTilda,
    const volVectorField& U,
    const volScalarField& d,
    const dimensionedScalar kappa,
    const dimensionedScalar cs
)
{
    volScalarField Omega = Foam::sqrt(2.0)*mag(skew(fvc::grad(U)));

    return
    (
        max
        (
            Omega + SpalartAllmaras::fv2(chi, fv1)*nuTilda/sqr(kappa*d),
            cs*Omega
        )
    );
}


tmp<volScalarField> nuTildaEff
(
    const volScalarField& nuTilda,
    const dimensionedScalar sigma,
    const dimensionedScalar nu
)
{
    return (nu + nuTilda)/sigma;
}


tmp<volScalarField> fw
(
    const volScalarField& Stilda,
    const volScalarField& nuTilda,
    const volScalarField& d,
    const dimensionedScalar kappa,
    const dimensionedScalar cw2,
    const dimensionedScalar cw3
)
{
    volScalarField r
    (
        min // min function due to floating point exception, see below pow6(r)
        (
            nuTilda
           /(
               max // max function to prevent dividing by zero
               (
                   Stilda,
                   dimensionedScalar("SMALL", Stilda.dimensions(), SMALL)
               )
              *sqr(kappa*d)
            ),
            scalar(10.0)
        )
    );
    r.boundaryField() == 0.0;

    const volScalarField g = r + cw2*(pow6(r) - r);

    return g*pow((1.0 + pow6(cw3))/(pow6(g) + pow6(cw3)), 1.0/6.0);
}


const volScalarField& d
(
    const fvMesh& mesh
)
{
    return wallDist::New(mesh).y();
}



// Defining partial derivatives

dimensionedScalar pChipNuTilda
(
    const dimensionedScalar nu
)
{
    return 1.0/nu;
}


tmp<volScalarField> pFv1pChi
(
    const volScalarField& chi,
    const dimensionedScalar cv1
)
{
    const volScalarField chi2 = Foam::pow(chi, 2);
    const volScalarField chi3 = Foam::pow3(chi);
    
    return 3.0*pow3(cv1)*chi2/pow(pow3(cv1) + chi3, 2);
}


tmp<volScalarField> pFv2pFv1
(
    const volScalarField& chi,
    const volScalarField& fv1
)
{   
    const volScalarField chi2 = pow(chi, 2);
    
    return chi2/pow(1.0 + chi*fv1, 2);
}


tmp<volScalarField> pStildapFv2
(
    const volScalarField& nuTilda,
    const volScalarField& d,
    const dimensionedScalar kappa
)
{
    return nuTilda/sqr(kappa*d);
}


tmp<volScalarField> pRpStilda
(
    const volScalarField& nuTilda,
    const volScalarField& d,
    const volScalarField& Stilda,
    const dimensionedScalar kappa
)
{
    return -nuTilda/sqr(kappa*d*Stilda);
}


tmp<volScalarField> pGpR
(
    const volScalarField& r,
    const dimensionedScalar cw2
)
{
    const volScalarField r5 = pow5(r);
    
    return 1.0 + cw2*(-1.0 + 6.0*r5);
}


tmp<volScalarField> pFwpG
(
    const volScalarField& g,
    const dimensionedScalar cw3
)
{   
    const volScalarField g6 = pow6(g);
    
    return pow6(cw3)*pow((1.0 + pow6(cw3))/(g6 + pow6(cw3)), 7/6)
          /(1.0 + pow6(cw3));
}


tmp<volScalarField> pFv2pChi
(
    const volScalarField& chi,
    const volScalarField& fv1
)
{
    return -1.0/sqr(1.0 + chi*fv1);
}


tmp<volScalarField> pStildapNuTilda
(
    const volScalarField& fv2,
    const volScalarField& d,
    const dimensionedScalar kappa
)
{
    return fv2/sqr(kappa*d);
}


tmp<volScalarField> pRpNuTilda
(
    const volScalarField& Stilda,
    const volScalarField& d,
    const dimensionedScalar kappa
)
{
    return 1.0
          /(
                max // avoid division by zero with SMALL
                (
                    Stilda,
                    dimensionedScalar("SMALL", Stilda.dimensions(), SMALL)
                )
               *sqr(kappa*d)
            );
}
           

// Defining derivatives for Spalart-Allmaras

// Defining the total derivative of term T2
tmp<volScalarField> calcDT2dNuTilda
(
    const volScalarField& d,
    const volScalarField& nuTilda,
    const volVectorField& U,
    const dimensionedScalar nu,
    const dimensionedScalar kappa,
    const dimensionedScalar cb1,
    const dimensionedScalar cv1
)
{
    const volScalarField Omega = Foam::sqrt(2.0)*mag(skew(fvc::grad(U)));
    
    return cb1/(nu*sqr(kappa*d))
          *(
               nu*Omega*sqr(kappa*d) 
             + 2.0*nu*nuTilda
             - 3.0*sqr(nuTilda)
             + 7.0*pow6(nuTilda)/(pow4(nu)*pow3(cv1) + pow3(nuTilda)*(nu + nuTilda))
             - (4.0*nuTilda + 3.0*nu)*pow(nuTilda, 9.0)/pow(pow3(nuTilda)*(nuTilda + nu) + pow4(nu)*pow3(cv1), 2.0)
            );
}


// Defining the total derivative of term T5
tmp<volScalarField> calcDT5dNuTilda
(
    const volScalarField& d,
    const volScalarField& nuTilda,
    const volVectorField& U,
    const dimensionedScalar nu,
    const dimensionedScalar cv1,
    const dimensionedScalar kappa,
    const dimensionedScalar cw1,
    const dimensionedScalar cw2,
    const dimensionedScalar cw3,
    const dimensionedScalar cs
)
{
    const volScalarField chi = SpalartAllmaras::chi(nuTilda, nu);
    const volScalarField fv1 = SpalartAllmaras::fv1(chi, cv1);
    const volScalarField fv2 = SpalartAllmaras::fv2(chi, fv1);
    const volScalarField Stilda = SpalartAllmaras::Stilda(chi, fv1, nuTilda, U, d, kappa, cs);
    const volScalarField fw = SpalartAllmaras::fw(Stilda, nuTilda, d, kappa, cw2, cw3);

    volScalarField r
    (
        min
        (
            nuTilda
           /(
               max
               (
                   Stilda,
                   dimensionedScalar("SMALL", Stilda.dimensions(), SMALL)
               )
              *sqr(kappa*d)
            ),
            scalar(10.0)
        )
    );
    r.boundaryField() == 0.0;

    const volScalarField g = r + cw2*(pow6(r) - r);

    return cw1*((pFwpG(g, cw3)*pGpR(r, cw2)*pRpNuTilda(Stilda, d, kappa)
          + pFwpG(g, cw3)*pGpR(r, cw2)*pRpStilda(nuTilda, d, Stilda, kappa)*pStildapNuTilda(fv2, d, kappa)
          + pFwpG(g, cw3)*pGpR(r, cw2)*pRpStilda(nuTilda, d, Stilda, kappa)*pStildapFv2(nuTilda, d, kappa)*pFv2pChi(chi, fv1)*pChipNuTilda(nu)
          + pFwpG(g, cw3)*pGpR(r, cw2)*pRpStilda(nuTilda, d, Stilda, kappa)*pStildapFv2(nuTilda, d, kappa)*pFv2pFv1(chi, fv1)*pFv1pChi(chi, cv1)*pChipNuTilda(nu)
                   )*sqr(nuTilda/d) + 2.0*fw/sqr(d)*nuTilda);
}

} // End namespace SpalartAllmaras

