if (adjAnalysis.match("all") || adjAnalysis.match("selective"))
{
    Info << "\n// * * * * * * * * * *          Adjoint method           * * * * * * * * * * //\n" << endl;
    
    // Create file stream object to write out the sensitivity data
    std::ofstream ofsAdj("sensitivityDataAdj.dat", std::ofstream::out);

    // Write header
    ofsAdj 
        << "stepSize \t xCoord \t yCoord \t zCoord \t adjSens \n";
    
    forAll(optimizationPatchIDs, optimizationPatchIDsI)  // Loop over all optimization patches
    {
        // Get global point IDs of the optimization patch
        labelList myPoints;
        if(adjAnalysis.match("all"))
        {
            myPoints = mesh.boundaryMesh()[optimizationPatchIDs[optimizationPatchIDsI]].meshPoints();  // Global Ids on the Patch
        }
        else if(adjAnalysis.match("selective"))
        {
            myPoints = globalIDsAdj;
        }

        // Get the point normal vectors of the optPatch
        const vectorField normalPointVector = mesh.boundaryMesh()[optimizationPatchIDs[optimizationPatchIDsI]].pointNormals();

        // Save the initial mesh point positions (i.e. the unperturbed mesh point positions)
        const pointField initMeshPoints = mesh.points();
        
        // Point counter
        label pointCounter = myPoints.size();
        
        forAll(myPoints, myPointI)  // loop over all points on the optimization patch
        {
            // Report remaining points to compute
            Info<< "\nPoints to compute " << pointCounter-- << " on patch " << optimizationPatchIDsI + 1 << " from " << optimizationPatchIDs.size() << " patch(es)" <<  endl;
            
            // Get global point ID
            label globalID = myPoints[myPointI];
            Info<< "globalPointID = " << globalID << endl;

            ofsAdj 
                << "PointID \t" << globalID << "\n";

            for(scalar FDStepLoop = FDStepStart; FDStepLoop <= FDStepStop; FDStepLoop *= FDStepSize)    // loop over the step sizes
            {
                scalar FDStep = FDStepLoop;
                Info<< "FDStep = " << FDStep << endl;

                // Write to file
                ofsAdj
                    << FDStep << "\t" << initMeshPoints[globalID].x() << "\t" << initMeshPoints[globalID].y() << "\t" << initMeshPoints[globalID].z() << "\t";
        

                // Forward perturbation
                
                // store the initial mesh points coordinates
                pointField pertMeshPoints = initMeshPoints;

                // save initial position of point
                vector oldPos = initMeshPoints[globalID];

                // define vector of perturbation direction
                vector pertVec = normalPointVector[myPointI];

                // assign new position of point myPointI
                vector newPos = oldPos + FDStep*pertVec;

                // update position in the mesh
                pertMeshPoints[globalID] = newPos;

                // move mesh position
                mesh.movePoints(pertMeshPoints);

                // calculate adjoint dot residual
                scalar lambdaDotResForw = SensitivityAnalysis::calcAdjointDotResidual(mesh, runTime, nuTilda, nut, adjnuTilda, U, p, dStar, rho, cb1, cb2, kappa, sigma, cs, nu, cw1, cw2, cw3, cv1, objectiveFunction, vectorForce, optimizationPatchIDsI);
                
                
                // Backward perturbation

                // store the initial mesh points coordinates
                pertMeshPoints = initMeshPoints;

                // save initial position of point
                oldPos = initMeshPoints[globalID];

                // assign new position of point i
                newPos = oldPos - FDStep*pertVec;

                // Update position in the mesh
                pertMeshPoints[globalID] = newPos;

                // Move mesh position
                mesh.movePoints(pertMeshPoints);

                // Calculate adjoint dot residual
                scalar lambdaDotResBack = SensitivityAnalysis::calcAdjointDotResidual(mesh, runTime, nuTilda, nut, adjnuTilda, U, p, dStar, rho, cb1, cb2, kappa, sigma, cs, nu, cw1, cw2, cw3, cv1, objectiveFunction, vectorForce, optimizationPatchIDsI);

                // Calculate sensitivity
                scalar adjSens = (lambdaDotResForw - lambdaDotResBack)/(2.0*FDStep);
        
                // Assign values to adjSensitivities field
                adjSensField[globalID] = adjSens;
                
                // Write to file
                ofsAdj
                    << adjSens << "\n";
                    
                // set mesh to initial state
                mesh.movePoints(initMeshPoints);
                
            } // Close the loop over the step size 
        } // Close the loop over the points on the opt patch
    }  // Close the loop over the patches

    // Close the ofstream object
    ofsAdj.close();

    // Write to disk
    runTime.write();
    
    Info << "\n// * * * * * * * * * *         End Adjoint method            * * * * * * * * //\n" << endl;
    
} // Close adjoint sensitivity analysis

