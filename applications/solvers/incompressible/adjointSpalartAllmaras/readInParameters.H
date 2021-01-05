// Get the analysis types
const word adjAnalysis = optimizationProperties.lookup("adjoint");
const word fdAnalysis = optimizationProperties.lookup("finiteDifference");

// Get the optimization patch IDs
const wordList optimizationPatchNames = optimizationProperties.lookup("optimizationPatchNames");
const label numOptimizationPatches = optimizationPatchNames.size();
labelList optimizationPatchIDs(numOptimizationPatches);
forAll(optimizationPatchNames, patchI)
{
    optimizationPatchIDs[patchI] = mesh.boundaryMesh().findPatchID(optimizationPatchNames[patchI]);
}

// Get the global IDs of the nodes
const labelList globalIDsAdj = optimizationProperties.lookup("globalIDsAdj");
const labelList globalIDsFD = optimizationProperties.lookup("globalIDsFD");

// Get the parameters for the step sizes study
const scalar FDStepStart = readScalar(optimizationProperties.lookup("stepStart"));
const scalar FDStepStop = readScalar(optimizationProperties.lookup("stepStop"));    
const scalar FDStepSize = readScalar(optimizationProperties.lookup("stepSize"));

// Read in the parameters for nonlinear equation solving
const label maxItr = readScalar(optimizationProperties.lookup("maxItr"));
const scalar eps = readScalar(optimizationProperties.lookup("eps"));
const scalar relaxFac = readScalar(optimizationProperties.lookup("relaxFac"));

// Read in the vector for the force direction
const vector vectorForce = optimizationProperties.lookup("vector");

// Get the objective function
const word objectiveFunction = optimizationProperties.lookup("objectiveFunction");

// Get the case example
const word caseExample = optimizationProperties.lookup("caseExample");

