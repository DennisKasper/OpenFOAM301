scalar calcAdjointDotResidual
(
    const fvMesh& mesh, 
    const Time& runTime, 
    const volScalarField& phi,
    const volScalarField& adjPhi
)
{
    // Compute residual
    volScalarField R = fvc::laplacian(phi);
    
    forAll(R,cI)
    {
        R[cI] *= mesh.V()[cI]; // Multiplication with volume due to fv dicretisation
    }

    scalar toReturn = 0.0;
    
    forAll(R, cI)
    {
        toReturn += (R[cI] * adjPhi[cI]);
    }

    return toReturn;
}


scalar calcPhiFD
(
    const fvMesh& mesh,
    const Time& runTime,
    const volScalarField& phi
)
{
    scalar phiLastCell = 0.0;
    
    fvScalarMatrix phiEqnFD(fvm::laplacian(phi));
    phiEqnFD.solve();
    
    phiLastCell = phi[phiEqnFD.source().size() - 1];

    return phiLastCell;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //