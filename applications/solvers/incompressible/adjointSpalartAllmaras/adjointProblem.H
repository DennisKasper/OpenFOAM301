// Computing the derivative of T2
const volScalarField dT2dNuTilda = SpalartAllmaras::calcDT2dNuTilda(d, nuTilda, U, nu, kappa, cb1, cv1);

// Computing the derivative of T5
const volScalarField dT5dNuTilda = SpalartAllmaras::calcDT5dNuTilda(d, nuTilda, U, nu, cv1, kappa, cw1, cw2, cw3, cs);

// Compute adjoint matrix
fvScalarMatrix adjnuTildaEqn
(   
    // T1
    fvm::div(phi, adjnuTilda)  // OK
    
    // T2
    - fvm::Sp(dT2dNuTilda, adjnuTilda)  // OK
    
    // T3a
    - (1.0 + cb2)*fvm::laplacian(SpalartAllmaras::nuTildaEff(nuTilda, sigma, nu), adjnuTilda)  // OK

    // T3b
    - (1.0 + cb2)*fvm::div(fvc::snGrad(nuTilda)*mesh.magSf()/sigma, adjnuTilda)  // OK

    // T4a
    + cb2*SpalartAllmaras::nuTildaEff(nuTilda, sigma, nu)*fvm::laplacian(adjnuTilda)  // OK
    
    // T4b
    + cb2*fvm::Sp(fvc::laplacian(nuTilda)/sigma, adjnuTilda)  // OK
    
    // T5
    + fvm::Sp(dT5dNuTilda, adjnuTilda)  // OK
);

// Get the off-diagonal terms
scalarField& upper = adjnuTildaEqn.upper();
scalarField& lower = adjnuTildaEqn.lower();

// Transpose the matrix
forAll(lower, i)
{
    scalar uI = upper[i];
    scalar lI = lower[i];
    lower[i] = uI;
    upper[i] = lI;
}

// Get the source term (RHS(i) = boundaryCoefficientFromDiscretisation(i) + sourceTerm(i))
scalarField& sourceCoeff = adjnuTildaEqn.source();

// Set the source term to zero to get rid of the boundary coeff from discrtization
forAll(sourceCoeff, sourceCoeffI)
    sourceCoeff[sourceCoeffI] = scalar(0.0);

if(objectiveFunction.match("artificial"))
{
    // Set four cells on top to -1
    if(caseExample.match("channel"))
    {
        sourceCoeff[157288] = scalar(-1.0);
        sourceCoeff[157367] = scalar(-1.0);
        sourceCoeff[157446] = scalar(-1.0);
        sourceCoeff[157525] = scalar(-1.0);
    }
    else if(caseExample.match("airFoil2D"))
    {
        forAll(optimizationPatchIDs, optimizationPatchIDsI)
        {
            const labelUList& cellId = mesh.boundaryMesh()[optimizationPatchIDs[optimizationPatchIDsI]].faceCells();
        
            forAll(cellId, cellIdI)
            {
                sourceCoeff[cellId[cellIdI]] = scalar(-1.0);
            }
        }
    }
    else if(caseExample.match("wing3D"))
    {
        // Set the first layer of cells around the wing to -1
        forAll(optimizationPatchIDs, optimizationPatchIDsI)
        {
            const labelUList& cellId = mesh.boundaryMesh()[optimizationPatchIDs[optimizationPatchIDsI]].faceCells();
        
            forAll(cellId, cellIdI)
            {
                sourceCoeff[cellId[cellIdI]] = scalar(-1.0);
            }
        }
    }
    else if(caseExample.match("unitCubeSpalartAllmaras"))
    {
        // Set the source term value of the last cell to -1
        sourceCoeff[sourceCoeff.size() - 1] = scalar(-1.0);   
    }
}
else if(objectiveFunction.match("coupling"))
{
    // Compute derivative of residual vector of Navier-Stokes
        const volVectorField pRNSpNut = - fvc::laplacian(U) - fvc::div(T(fvc::grad(U)));
//         pRNSpNut = -(fvc::laplacian(U));// + fvc::div(dev2(T(fvc::grad(U)))));
//         fvScalarMatrix whatever2(-fvm::laplacian(nut));
    
//         fvVectorMatrix test( -fvm::div(U,nut));
    Info<<"pRNSpNut[91126] = " << pRNSpNut[91126] << endl;

//         forAll(pRNSpNut, cI)
//         {
//             pRNSpNut[cI] *= mesh.V()[cI];
//         }
    
    // Computing the derivative for the RANS term
    pNutpNuTilda = NavierStokes::calcDNuTdNuTilda(nuTilda, nu, cv1);
    Info<<"pNutpNuTilda[91126] = " << pNutpNuTilda[91126] << endl;
    
    // Computing partial RNS / partial nuTilda
    const volVectorField pRNSpNuTilda  = pRNSpNut*pNutpNuTilda;
    
    
//         forAll(sourceCoeff, sourceCoeffI)
//         {
//             pRNSpNuTilde[sourceCoeffI] *= mesh.V()[sourceCoeffI];
//         }
    
    // Assign value to source term
    sourceCoeff = -(dStar & pRNSpNuTilda);
}
else if(objectiveFunction.match("force"))
{
        pNutpNuTilda = NavierStokes::calcDNuTdNuTilda(nuTilda, nu, cv1);
        
        volSymmTensorField devRhoReff = -1.0*rho*dev(twoSymm(fvc::grad(U))); 
        volSymmTensorField::GeometricBoundaryField& devRhoReffb = devRhoReff.boundaryField();
        
    forAll(optimizationPatchIDs, optimizationPatchIDsI)  // Loop over all optimization patches
    {
        label patchID = optimizationPatchIDs[optimizationPatchIDsI];
        vectorField  d_fT_d_nut = mesh.Sf().boundaryField()[patchID] & devRhoReffb[patchID];
        
        forAll(d_fT_d_nut,fI)
        {
            label cellID = mesh.boundaryMesh()[patchID].faceCells()[fI];
            sourceCoeff[cellID] = -1.0*(d_fT_d_nut[fI] & vectorForce)*pNutpNuTilda[cellID];
        }
    }
}
else
{
    FatalErrorIn("adjointProblem.H or scalar SensitivityAnalysis::calcnuTildaFD(...)")
        << "Wrong objectiveFunction in optimzationProperties"
        << abort(FatalError);
}

// Set the boundary coefficients to zero
forAll(adjnuTilda.boundaryField(), patchI)
{
    forAll(adjnuTildaEqn.boundaryCoeffs()[patchI], faceI)
    {
        adjnuTildaEqn.boundaryCoeffs()[patchI][faceI] = scalar(0.0);
    }
}

// Assign the source term to the field adjSource
forAll(sourceCoeff, sourceCoeffI)
{
    adjSource[sourceCoeffI] = sourceCoeff[sourceCoeffI];
}

// Solve for the adjoint variables
adjnuTildaEqn.solve();

// Write to disk
runTime.write();
