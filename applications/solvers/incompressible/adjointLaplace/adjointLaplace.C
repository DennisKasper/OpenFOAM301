/*---------------------------------------------------------------------------*\

Application
    adjointLaplace.C

Description
    Solves the primal and adjoint problem for the Laplace equation and computes
    the sensitivities.

Author
    Dennis Kasper (dennis.m.kasper@gmail.com)

\*---------------------------------------------------------------------------*/

// Include OpenFOAM directive
#include "fvCFD.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "valuePointPatchFields.H"

// Include  C++ directive
#include <time.h>
#include <fstream>

// Include User directive
#include "sensitivityAnalysis.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"    // Set the correct path
    #include "createTime.H"     // Create the time
    #include "createMesh.H"     // Create the mesh
    #include "createFields.H"   // Create the fields

    Info << "\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << endl;
    Info << "// * * * * * * * * * *       Solving primal problem      * * * * * * * * * * //\n" << endl;
    
    // Advance in pseudo time
    runTime++;
    Info << "runTime = " << runTime.timeName() << nl << endl;

    // Setup equation system for phi and solve
    fvScalarMatrix phiEqn(fvm::laplacian(phi));
    phiEqn.solve();
    
    // Write to disk
    runTime.write();
    

    Info << "\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << endl;
    Info << "// * * * * * * * * * *      Solving adjoint problem      * * * * * * * * * * //\n" << endl;

    // Setup equation system for adjPhi
    fvScalarMatrix adjPhiEqn(fvm::laplacian(adjPhi));
    
    // Transpose the matrix
    scalarField& upper = adjPhiEqn.upper();
    scalarField& lower = adjPhiEqn.lower();
    
    forAll(lower, i)
    {
        scalar uI = upper[i];
        scalar lI = lower[i];
        lower[i] = uI;
        upper[i] = lI;
    }

    // Get the source term (RHS(i) = boundaryCoefficientFromDiscretisation(i) + sourceTerm(i))
    scalarField& sourceCoeff = adjPhiEqn.source();

    // Set the source term to zero
    forAll(sourceCoeff, i)
        sourceCoeff[i] = scalar(0.0);

    // Set the source term value of the last cell to -1
    sourceCoeff[sourceCoeff.size() - 1] = scalar(-1.0);

    // Assign the source term to the field adjSource
    forAll(sourceCoeff, i)
    {
        adjSource[i] = sourceCoeff[i];
    }

    // Set the boundary coefficients to zero
    forAll(adjPhi.boundaryField(), patchI)
    {
        forAll(adjPhiEqn.boundaryCoeffs()[patchI], faceI)
        {
            adjPhiEqn.boundaryCoeffs()[patchI][faceI] = scalar(0.0); //temp;
        }
    }

    // Solve for the adjoint variables
    adjPhiEqn.solve();


    Info << "\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << endl;
    Info << "// * * * * * * * * * *       Sensitivity analysis        * * * * * * * * * * //\n" << endl;

    // Get optimization patch ID
    label patchID = mesh.boundaryMesh().findPatchID("optPatch");

    // Get global point IDs of the optimization patch
    const labelList& myPoints = mesh.boundaryMesh()[patchID].meshPoints();  // mapping  to global Ids (mypoints are global)
    
    // Get the point normal vectors of the optPatch
    const vectorField normalPointVector = mesh.boundaryMesh()[patchID].pointNormals();
    
    // Save the initial mesh point positions (i.e. the unperturbed mesh point positions)
    const pointField initMeshPoints = mesh.points();

    // Parameter for the step sizes study
    scalar FDStepStart = 1e-1;
    scalar FDStepStop = 1e-10;
    scalar FDStepSize = 1e-1;

    // Account for computing time
    clock_t clockAdj = std::clock();
            
    // Adjoint sensitivity analysis
    if (true)
    {
        Info << "\n// * * * * * * * * * *          Adjoint method           * * * * * * * * * * //\n" << endl;

        // Create file stream object to write out the sensitivity data
        std::ofstream ofsAdj("sensitivityDataAdj.txt", std::ofstream::out);

        // Write header
        ofsAdj 
            << "stepSize \t xCoord \t yCoord \t zCoord \t adjSens \n";

        // Point counter
        label pointCounter = myPoints.size();
        
        forAll(myPoints, myPointI)    // Loop over all points on the optimization patch
        {
            // Report remaining points to compute
            Info<< "\nPoints remain = " << pointCounter-- << endl;
            
            // Get global point ID i on patch
            label globalID = myPoints[myPointI];
            Info<< "\nglobalPointID = " << globalID << endl;

            {
                ofsAdj 
                    << "PointID \t" << globalID << "\n";

                for(scalar FDStepLoop = FDStepStart; FDStepLoop >= FDStepStop; FDStepLoop *= FDStepSize)	// loop over the step sizes
                {
                    scalar FDStep = FDStepLoop;
                    Info<< "FDStep = " << FDStep << endl;

                    // Write to file
                    ofsAdj
                        << FDStep << "\t" << initMeshPoints[globalID].x() << "\t" << initMeshPoints[globalID].y() << "\t" << initMeshPoints[globalID].z() << "\t";

                    // Forward perturbation
                    Info<<"Forward"<<endl;
                        
                    // Store the initial mesh points coordinates
                    pointField pertMeshPoints = initMeshPoints;

                    // Save initial position of point
                    vector oldPos = initMeshPoints[globalID];

                    // Define vector of perturbation direction
                    vector pertVec = normalPointVector[myPointI];

                    // Assign new position of point myPointI
                    vector newPos = oldPos + FDStep*pertVec;

                    // Update position in the mesh
                    pertMeshPoints[globalID] = newPos;

                    // Move mesh position
                    mesh.movePoints(pertMeshPoints);

                    // Calculate adjoint dot residual
                    scalar lambdaDotResForw = calcAdjointDotResidual(mesh, runTime, phi, adjPhi);


                    // Backward perturbation
                    Info<<"Backward"<<endl;

                    // Store the initial mesh points coordinates
                    pertMeshPoints = initMeshPoints;

                    // Save initial position of point
                    oldPos = initMeshPoints[globalID];

                    // Assign new position of point i
                    newPos = oldPos - FDStep*pertVec;

                    // Update position in the mesh
                    pertMeshPoints[globalID] = newPos;

                    // Move mesh position
                    mesh.movePoints(pertMeshPoints);
                    
                    // Compute adjoint dot residual
                    scalar lambdaDotResBack = calcAdjointDotResidual(mesh, runTime, phi, adjPhi);

                    
                    // Compute sensitivity
                    scalar delta = 2.0*FDStep;
                    scalar adjSens = (lambdaDotResForw - lambdaDotResBack) / delta;

                    // Assign values to adjSensitivities field
                    adjSensField[globalID] =  adjSens;

                    // Write to file
                    ofsAdj
                        << adjSens << "\n";

                } // Close the loop over the step size

            clockAdj = std::clock() - clockAdj;

            // set mesh to initial state
            mesh.movePoints(initMeshPoints);

            } // close the if statemnent loop
                        
        } // Close the loop over the points on the opt patch

        // Close the ofstream object
        ofsAdj.close();

        // Write to disk
        runTime.write();
        
    } // Close adjoint sensitivity analysis



    // Sensitivity total finite difference: * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
            
    // Acoount for computation time
    clock_t clockFD = std::clock();

    // FD sensitivity analysis
    if (true)
    {
        Info << "\n// * * * * * * * * * *       Total finite difference method      * * * * * * * * * * //\n" << endl;

        // Create file stream object to write out the sensitivity data
        std::ofstream ofsFD("sensitivityDataFD.txt", std::ofstream::out);

        // Write file header
        ofsFD 
            << "stepSize \t xCoord \t y Coord \t zCoord \t fdSens \n";
        
        // Point counter
        label pointCounter = myPoints.size();
        
        forAll(myPoints, myPointI)  // loop over all points on the optPatch
        {
            // Report remaining points to compute
            Info<< "\nPoints remain = " << pointCounter-- << endl;
            
            // get global point ID i on patch
            label globalID = myPoints[myPointI];
            Info<< "\nglobalPointID = " << globalID << endl;

            {
                // Write to file
                ofsFD
                    << "PointID \t" << globalID << "\n";

                for(scalar FDStepLoop = FDStepStart; FDStepLoop >= FDStepStop; FDStepLoop *= FDStepSize)	// loop over the step sizes
                {
                    scalar FDStep = FDStepLoop;
                    Info<< "FDStep = " << FDStep << endl;

                    // Write to file
                    ofsFD 
                        << FDStep << "\t" << initMeshPoints[globalID].x() << "\t" << initMeshPoints[globalID].y() << "\t" << initMeshPoints[globalID].z() << "\t";

                        
                    // Forward perturbation
                    Info<<"Forward"<<endl;

                    // Save the initial mesh points coordinates
                    pointField pertMeshPoints = initMeshPoints;

                    // Save initial position of point
                    vector oldPos = initMeshPoints[globalID];

                    // Define vector of perturbation direction
                    vector pertVec = normalPointVector[myPointI];

                    // Assign new position of point i
                    vector newPos = oldPos + FDStep*pertVec;

                    // Update position in the mesh
                    pertMeshPoints[globalID] = newPos;

                    // Move mesh position
                    mesh.movePoints(pertMeshPoints);  // assign new position to mesh

                    // Compute phi
                    scalar phiForw = calcPhiFD(mesh, runTime, phi);

                    
                    // Backward perturbation
                    Info<<"Backward"<<endl;

                    // Save the initial mesh points coordinates
                    pertMeshPoints = initMeshPoints;

                    // Save initial position of point
                    oldPos = initMeshPoints[globalID];

                    // Assign new position of point i
                    newPos = oldPos - FDStep*pertVec;  // calculate pertubation in negative x-direction

                    // Update position in the mesh
                    pertMeshPoints[globalID] = newPos;

                    // Move mesh position
                    mesh.movePoints(pertMeshPoints);  // assign new position to mesh

                    // Compute phi
                    scalar phiBack = calcPhiFD(mesh, runTime, phi);

                    
                    // Compute sensitivity
                    scalar delta = 2.0*FDStep;
                    scalar fdSens = (phiForw - phiBack) / delta;
                   // dfdxFD[myPointI].z() = sensitivityZFD;

                    // Assign values to FDsensitivity field
                    fdSensField[globalID] =  fdSens;

                    // Write to file
                    ofsFD 
                        << fdSens << "\n";

                    // set mesh to initial state
                    mesh.movePoints(initMeshPoints);

                } // step size loop
            
            } // if statement loop

        } // points on opti patch loop

        // Compute the cpu time
        clockFD = std::clock() - clockFD;

        // Close the ofstream object
        ofsFD.close();

        // Write to disk
        runTime.write();

    } // Close FD sensitivity analysis

    // Write to disk
    runTime.write();

    // Report computation time
    Info<< nl << "Execution time Adjoint solver = " << clockAdj / CLOCKS_PER_SEC << " s" << nl << endl;
    Info<< nl << "Execution time FD solver = " << clockFD / CLOCKS_PER_SEC << " s" << nl << endl;
    Info<< nl << "Total executionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //? Info<< "nu = " << nu << endl;
    Info<< nl << "End" << nl << endl;
 
    return 0;
}

// ************************************************************************* //


