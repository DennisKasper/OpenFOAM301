// Initialize parameters
label iter = 0;
scalar maxResidual = 1.0;

// Solve primal problem
while (iter < maxItr && maxResidual > eps)
{
    iter += 1.0;
    Info<<"Iteration " << iter << endl;
    
    nuTilda.storePrevIter();
        
    const volScalarField chi = SpalartAllmaras::chi(nuTilda, nu);
    const volScalarField fv1 = SpalartAllmaras::fv1(chi, cv1);
    const volScalarField Stilda = SpalartAllmaras::Stilda(chi, fv1, nuTilda, U, d, kappa, cs);
    
    fvScalarMatrix nuTildaEqn
    (
        // T1
        fvm::div(phi, nuTilda)
                    
        // T2
        - fvm::SuSp(cb1*Stilda, nuTilda)
//             - cb1*Stilda*nuTilda
        // T3
        - (1.0 + cb2)*fvm::laplacian(SpalartAllmaras::nuTildaEff(nuTilda, sigma, nu), nuTilda)
        
        // T4
        + cb2*SpalartAllmaras::nuTildaEff(nuTilda, sigma, nu)*fvm::laplacian(nuTilda)

        // T5
        + fvm::SuSp(cw1*SpalartAllmaras::fw(Stilda, nuTilda, d, kappa, cw2, cw3)*nuTilda/sqr(d), nuTilda)
//             + cw1*SpalartAllmaras::fw(Stilda, nuTilda, d, kappa, cw2, cw3)*sqr(nuTilda/d)
    );

    nuTildaEqn.relax(relaxFac);      
    maxResidual = solve(nuTildaEqn, mesh.solver("nuTilda")).initialResidual();
    bound(nuTilda, dimensionedScalar("0", nuTilda.dimensions(), 0.0));
    nuTilda.correctBoundaryConditions();
    
    nut = nuTilda*fv1;
    nut.correctBoundaryConditions();
}

// Check for residuals
//     checkRNS = fvc::div(phi, U) + fvc::grad(p) - fvc::laplacian((nu + nut), U) - fvc::div((nu + nut)*dev2(T(fvc::grad(U))));

// Write to disk
runTime.write(); 
