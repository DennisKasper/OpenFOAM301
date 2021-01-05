if(fdAnalysis.match("all") || fdAnalysis.match("selective"))
{
    Info << "\n// * * * * * * * * * *       Total finite difference method      * * * * * * * * * * //\n" << endl;

    // Create file stream object to write out the sensitivity data
    std::ofstream ofsFD("sensitivityDataFD.dat", std::ofstream::out);

    // Write file header
    ofsFD 
        << "stepSize \t xCoord \t yCoord \t zCoord \t fdSens \n";

    forAll(optimizationPatchIDs, optimizationPatchIDsI)  // Loop over all optimization patches
    {
        // Get global point IDs of the optimization patch
        labelList myPoints;
        if(fdAnalysis.match("all"))
        {
            myPoints = mesh.boundaryMesh()[optimizationPatchIDs[optimizationPatchIDsI]].meshPoints();  // Global Ids on the Patch
        }
        else if(fdAnalysis.match("selective"))
        {
            myPoints = globalIDsFD;
        }

        // Get the point normal vectors of the optPatch
        const vectorField normalPointVector = mesh.boundaryMesh()[optimizationPatchIDs[optimizationPatchIDsI]].pointNormals();

        // Save the initial mesh point positions (i.e. the unperturbed mesh point positions)
        const pointField initMeshPoints = mesh.points();
        
        // Point counter
        label pointCounter = myPoints.size();
        
        forAll(myPoints, myPointI)  // loop over all points on the optimization patch 
        {
            //Report remaining points to compute
            Info<< "\nPoints to compute " << pointCounter-- << " on patch " << optimizationPatchIDsI + 1 << " from " << optimizationPatchIDs.size() << " patch(es)" <<  endl;

            // Get global point ID
            label globalID = myPoints[myPointI];
            Info<< "globalPointID = " << globalID << endl;
            
            ofsFD
                << "PointID \t" << globalID << "\n";

            for(scalar FDStepLoop = FDStepStart; FDStepLoop <= FDStepStop; FDStepLoop *= FDStepSize)    // loop over the step sizes
            {
                scalar FDStep = FDStepLoop;
                Info<< "\nFDStep = " << FDStep << endl;

                // Write to file
                ofsFD 
                    << FDStep << "\t" << initMeshPoints[globalID].x() << "\t" << initMeshPoints[globalID].y() << "\t" << initMeshPoints[globalID].z() << "\t";


                // Forward perturbation
                
                // Save the initial mesh points coordinates
                pointField pertMeshPoints = initMeshPoints;

                // Save initial position of point
                vector oldPos = initMeshPoints[globalID];

                // Define vector of perturbation direction
                vector pertVec = normalPointVector[myPointI];

                // Assign new position of point
                vector newPos = oldPos + FDStep*pertVec;  // calculate pertubation in y direction

                // update position in the mesh
                pertMeshPoints[globalID] = newPos;

                // move mesh position
                mesh.movePoints(pertMeshPoints);  // assign new position to mesh
                
                // compute nuTilda
                scalar nuTildaForw = SensitivityAnalysis::calcnuTildaFD(mesh, runTime, U, p, dStar, nut, nuTilda, rho, cb1, cb2, kappa, sigma, cs, nu, cw1, cw2, cw3, cv1, objectiveFunction, maxItr, eps, relaxFac, caseExample, vectorForce, optimizationPatchIDsI, optimizationPatchIDs);
                
                // Backward perturbation

                // store the initial mesh points coordinates
                pertMeshPoints = initMeshPoints;

                // save initial position of point
                oldPos = initMeshPoints[globalID];

                // assign new position of point i
                newPos = oldPos - FDStep*pertVec;  // calculate pertubation in negative x-direction

                // update position in the mesh
                pertMeshPoints[globalID] = newPos;

                // move mesh position
                mesh.movePoints(pertMeshPoints);  // assign new position to mesh

                // Compute nuTilda
                scalar nuTildaBack = SensitivityAnalysis::calcnuTildaFD(mesh, runTime, U, p, dStar, nut, nuTilda, rho, cb1, cb2, kappa, sigma, cs, nu, cw1, cw2, cw3, cv1, objectiveFunction, maxItr, eps, relaxFac, caseExample, vectorForce, optimizationPatchIDsI, optimizationPatchIDs);

                // calculate sensitivity
                scalar fdSens = (nuTildaForw - nuTildaBack)/(2.0*FDStep);


                // Assign values to adjSensitivities field
                fdSensField[globalID] = fdSens;
                
                // Write to file
                ofsFD 
                    << fdSens << "\n";
                
                // set mesh to initial state
                mesh.movePoints(initMeshPoints);
        
            } // step size loop   
        } // points on opti patch loop
    }  // patches loop
        
    // Close the ofstream object
    ofsFD.close();

    // Write to disk
    runTime.write();

    Info << "\n// * * * * * * * *   End Total finite difference method    * * * * * * * * * //\n" << endl;

} // Close FD sensitivity analysis
