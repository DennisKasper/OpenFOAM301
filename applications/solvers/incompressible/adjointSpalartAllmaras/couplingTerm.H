
// Double check the extra terms total finite difference
/*    
// Advance in pseudo time
runTime++;
Info << "runTime = " << runTime.timeName() << nl << endl;

// Check partial RNS / partial nut --------------------------------------- Not possible!?

// Save the initial nut field
const volScalarField nutTemp = nut;

// Create Delta nut
const scalar deltaNut = 0.000000001;

// perturb Forward
nut[91126] += deltaNut;

// Solve the momentum residual forward
volVectorField RmForw = - fvc::laplacian(nut, U); //- fvc::div((nu + nut)*dev2(T(fvc::grad(U))));
    
// Reset the nut field
nut = nutTemp;

// perturb backward
nut[91126] -= deltaNut;

// Solve the momentum residual backward
volVectorField RmBack = - fvc::laplacian(nut, U); // - fvc::div((nu + nut)*dev2(T(fvc::grad(U))));

// Compute total finite difference
pRNSpNutFD = (RmForw - RmBack)/(2.0*deltaNut);
volScalarField dStarDotpRNSpNutFD = dStar & pRNSpNutFD;
scalar whatever = gSum(dStarDotpRNSpNutFD);

Info<<"whatever = " << whatever << endl; 

// Reset the nut field
nut = nutTemp;


// Check partial nut / partial nuTilda ----------------------------------- OK
// Cell with high nuTilda value 91126

// Save the initial nuTilda field
const volScalarField nuTildaTemp = nuTilda;

// Create Delta nutilda
const scalar deltaNuTilda = 0.000000001;

// Perturb Forward
nuTilda[91126] += deltaNuTilda;

// Compute chi and fv1
volScalarField chi = SpalartAllmaras::chi(nuTilda, nu);
volScalarField fv1 = SpalartAllmaras::fv1(chi, cv1);

// Solve nut forward
const volScalarField nutForw = nuTilda*fv1;

// Reset the fields
nut = nutTemp;
nuTilda = nuTildaTemp;

// Perturb backward
nuTilda[91126] -= deltaNuTilda;

// Compute chi and fv1
chi = SpalartAllmaras::chi(nuTilda, nu);
fv1 = SpalartAllmaras::fv1(chi, cv1);

// Solve nut backward
const volScalarField nutBack = nuTilda*fv1;

// Compute sensitivity
pNutpNuTildaFD = (nutForw - nutBack)/(2.0*deltaNuTilda);
Info<<"pNutpNuTildaFD[91126] = " << pNutpNuTildaFD[91126] << endl;

// Reset the fields
nut = nutTemp;
nuTilda = nuTildaTemp;   

runTime.write();

*/
