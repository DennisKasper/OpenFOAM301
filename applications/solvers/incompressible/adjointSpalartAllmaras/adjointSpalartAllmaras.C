/*---------------------------------------------------------------------------*\

Application
    adjointSpalartAllmaras.C

Description
    Adjoint solver for the Spalart-Allmaras turbulence model

Author
    Dennis Kasper (dennis.m.kasper@gmail.com)

\*---------------------------------------------------------------------------*/

// Include OpenFOAM directive
#include "fvCFD.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "wallDist.H"
#include "bound.H"
// Added due to the nutUSpalartAllmaras wall function error
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
// to create the scalarPointField
#include "valuePointPatchFields.H"

// Include  C++ directive
#include <fstream>

// Include User directive
#include "spalartAllmaras.H"
#include "navierStokes.H"
#include "sensitivityAnalysis.H"


int main(int argc, char *argv[])
{
    #include "setRootCase.H"      // Set the correct path
    #include "createTime.H"       // Create the time
    #include "createMesh.H"       // Create the mesh
    #include "createFields.H"     // Create the fields
    #include "createConstants.H"  // Create the Spalart-Allmaras constants
    #include "readInParameters.H" // Read in the paramters
  
    // Advance in pseudo time
    runTime++;
    Info << "runTime = " << runTime.timeName() << nl << endl;


    Info << "\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << endl;
    Info << "// * * * * * * * * * *       Solving primal problem      * * * * * * * * * * //\n" << endl;
  
    #include "primalProblem.H"
  
    Info << "\n// * * * * * * * * *     End solving primal problem    * * * * * * * * * * //" << endl;
    Info << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n" << endl;


    Info << "\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << endl;
    Info << "// * * * * * * * * * *      Solving adjoint problem      * * * * * * * * * * //\n" << endl;

    #include "adjointProblem.H"
    
    Info << "\n// * * * * * * * * *    End solving adjoint problem    * * * * * * * * * * //" << endl;
    Info << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n" << endl;
    
    
    Info << "\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << endl;
    Info << "// * * * * * * * * * *       Sensitivity analysis        * * * * * * * * * * //\n" << endl;
    
    #include "adjointMethod.H"
    
    #include "FDMethod.H"

    Info << "// * * * * * * * * * *     End sensitivity analysis      * * * * * * * * * * //" << endl;
    Info << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n" << endl;

    
    Info<< nl << "Total executionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;
    Info<< nl << "End" << nl << endl;
       
    return 0;
}


