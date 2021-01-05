namespace SensitivityAnalysis
{
    
scalar calcAdjointDotResidual
(
    const fvMesh& mesh, 
    const Time& runTime, 
    const volScalarField& nuTilda,
    const volScalarField& nut,
    const volScalarField& adjnuTilda,
    const volVectorField& U,
    const volScalarField& p,
    const volVectorField& dStar,
    const dimensionedScalar rho,
    const dimensionedScalar cb1,
    const dimensionedScalar cb2,
    const dimensionedScalar kappa,
    const dimensionedScalar sigma,
    const dimensionedScalar cs,
    const dimensionedScalar nu,
    const dimensionedScalar cw1,
    const dimensionedScalar cw2,
    const dimensionedScalar cw3,
    const dimensionedScalar cv1,
    const word objectiveFunction,
    const vector vectorForce,
    const label optimizationPatchIDsI
)
{
    const volScalarField d = SpalartAllmaras::d(mesh);
    const volScalarField chi = SpalartAllmaras::chi(nuTilda, nu);
    const volScalarField fv1 = SpalartAllmaras::fv1(chi, cv1);
    const volScalarField fv2 = SpalartAllmaras::fv2(chi, fv1);
    const volScalarField Stilda = SpalartAllmaras::Stilda(chi, fv1, nuTilda, U, d, kappa, cs);
    
     const surfaceScalarField phi = linearInterpolate(U) & mesh.Sf();
    
    // Compute Spalart-Allmaras residual   
    volScalarField  RSA = 
        // T1
        fvc::div(phi, nuTilda)
        
        // T2
        - cb1*Stilda*nuTilda  // original       

        // T3
        - (1.0 + cb2)*fvc::laplacian(SpalartAllmaras::nuTildaEff(nuTilda, sigma, nu), nuTilda)
        
        // T4
        +  cb2*SpalartAllmaras::nuTildaEff(nuTilda, sigma, nu)*fvc::laplacian(nuTilda)
        
        // T5
        +  cw1*SpalartAllmaras::fw(Stilda, nuTilda, d, kappa, cw2, cw3)*sqr(nuTilda/d);

    forAll(RSA, cI)
    {
        RSA[cI] *= mesh.V()[cI]; // multiplication with volume due to fv dicretisation
    }    
        
    scalar toReturn = 0.0;
    
    if(objectiveFunction.match("artificial"))
    {
        forAll(RSA, cI)
        {
            toReturn += adjnuTilda[cI]*RSA[cI];
        }
    }
    else if(objectiveFunction.match("coupling"))
    {
        // Compute Navier-Stokes residual
        volVectorField RNS = fvc::div(phi, U) + fvc::grad(p) - fvc::laplacian((nu + nut), U) - fvc::div((nu + nut)*(T(fvc::grad(U))));
        
        forAll(RNS, cI)
        {
            RNS[cI] *= mesh.V()[cI];
        }
        
        // Compute objective function
        const volScalarField innerProduct = dStar & RNS;
        
        // Compute sensitivity
        forAll(RSA, cI)
        {
            toReturn += innerProduct[cI] + adjnuTilda[cI]*RSA[cI];
        }
    }
    else if(objectiveFunction.match("force"))
    {
        //calculation of fT
        volSymmTensorField devRhoReff = -1.0*rho*(nu + nut)*dev(twoSymm(fvc::grad(U)));
        volSymmTensorField::GeometricBoundaryField& devRhoReffb = devRhoReff.boundaryField();
        
        scalar fT_sum = gSum(mesh.Sf().boundaryField()[optimizationPatchIDsI] & devRhoReffb[optimizationPatchIDsI] & vectorForce);
     
        forAll(RSA, cI)
        {
            toReturn += adjnuTilda[cI]*RSA[cI];
        }
        
        toReturn += fT_sum;
    }

    return toReturn;
}



scalar calcnuTildaFD
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    const volScalarField& p,
    const volVectorField& dStar,
    volScalarField& nut,
    volScalarField& nuTilda,
    const dimensionedScalar rho,
    const dimensionedScalar cb1,
    const dimensionedScalar cb2,
    const dimensionedScalar kappa,
    const dimensionedScalar sigma,
    const dimensionedScalar cs,
    const dimensionedScalar nu,
    const dimensionedScalar cw1,
    const dimensionedScalar cw2,
    const dimensionedScalar cw3,
    const dimensionedScalar cv1,
    const word objectiveFunction,
    const label maxItr,
    const scalar eps,
    const scalar relaxFac,
    const word caseExample,
    const vector vectorForce,
    const label optimizationPatchIDsI,
    const labelList optimizationPatchIDs
)
{
    const volScalarField d = SpalartAllmaras::d(mesh);
    
    const surfaceScalarField phi = linearInterpolate(U) & mesh.Sf();
    
    label iter = 0;
    scalar maxResidual = 1.0;
       
    while (iter < maxItr && maxResidual > eps)
    {
        iter += 1.0;
        Info<< "\nIteration = " << iter << endl;
        
        nuTilda.storePrevIter();
        
        const volScalarField chi = SpalartAllmaras::chi(nuTilda, nu);
        const volScalarField fv1 = SpalartAllmaras::fv1(chi, cv1);        
        const volScalarField Stilda = SpalartAllmaras::Stilda(chi, fv1, nuTilda, U, d, kappa, cs);
        
        fvScalarMatrix nuTildaEqnFD
        (
            // T1
            fvm::div(phi, nuTilda) 
                   
            // T2
            - fvm::SuSp(cb1*Stilda, nuTilda)
            
            // T3
            - (1.0 + cb2)*fvm::laplacian(SpalartAllmaras::nuTildaEff(nuTilda, sigma, nu), nuTilda)
            
            // T4
            + cb2*SpalartAllmaras::nuTildaEff(nuTilda, sigma, nu)*fvm::laplacian(nuTilda)
            
            // T5
            + fvm::SuSp(cw1*SpalartAllmaras::fw(Stilda, nuTilda, d, kappa, cw2, cw3)*nuTilda/sqr(d), nuTilda)
        );
        
        nuTildaEqnFD.relax(relaxFac);
        maxResidual = solve(nuTildaEqnFD, mesh.solver("nuTilda")).initialResidual();
        bound(nuTilda, dimensionedScalar("0", nuTilda.dimensions(), 0.0));
        nuTilda.correctBoundaryConditions();
        
        nut = nuTilda*fv1;
        nut.correctBoundaryConditions();
    }
    
    scalar toReturn = 0.0;
    
    if(objectiveFunction.match("artificial"))
    {
        if(caseExample.match("channel"))
        {
            toReturn = nuTilda[157288] + nuTilda[157367] + nuTilda[157446] + nuTilda[157525];
        }
        else if(caseExample.match("airFoil2D"))
        {
            forAll(optimizationPatchIDs, optimizationPatchIDsI)
            {
                const labelUList& cellId = mesh.boundaryMesh()[optimizationPatchIDs[optimizationPatchIDsI]].faceCells();
            
                forAll(cellId, cellIdI)
                {
                    toReturn += nuTilda[cellId[cellIdI]];
                }
            }
        }
        else if(caseExample.match("wing3D"))
        {
            forAll(optimizationPatchIDs, optimizationPatchIDsI)
            {
                const labelUList& cellId = mesh.boundaryMesh()[optimizationPatchIDs[optimizationPatchIDsI]].faceCells();
            
                forAll(cellId, cellIdI)
                {
                    toReturn += nuTilda[cellId[cellIdI]];
                }
            }
        }
        else if(caseExample.match("unitCubeSpalartAllmaras"))
        {
            toReturn = nuTilda[nuTilda.size() - 1];    
        }        
    }
    else if(objectiveFunction.match("coupling"))
    {
        // Compute residual Navier-Stokes vector
        const volVectorField RNS = fvc::div(phi, U) + fvc::grad(p) - fvc::laplacian((nu + nut), U) - fvc::div((nu + nut)*(T(fvc::grad(U))));
        
//             forAll(RNS, cI)
//             {
//                 RNS[cI] *= mesh.V()[cI];
//             }
        
        // Compute inner product
        const volScalarField innerProduct = dStar & RNS;
    
        // Compute the sum
        forAll(innerProduct, i)
        {
            toReturn += innerProduct[i];
        }
    }
    else if(objectiveFunction.match("force"))
    {
        //calculation of fT
        volSymmTensorField devRhoReff = -1.0*rho*(nu + nut)*dev(twoSymm(fvc::grad(U))); 
        volSymmTensorField::GeometricBoundaryField& devRhoReffb = devRhoReff.boundaryField();
        
//         scalar fT_sum = gSum(mesh.Sf().boundaryField()[optimizationPatchIDsI] & devRhoReffb[optimizationPatchIDsI] & vectorForce);
        toReturn = gSum(mesh.Sf().boundaryField()[optimizationPatchIDsI] & devRhoReffb[optimizationPatchIDsI] & vectorForce);
    }
    
    return toReturn;
}
    
} // End namespace SensitivityAnalysis
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //