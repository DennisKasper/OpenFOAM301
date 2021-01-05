Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info << "Reading field nut\n" << endl;

volScalarField nut
(
    IOobject
    (
        "nut",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field nuTilda\n" << endl;

volScalarField nuTilda
(
    IOobject
    (
        "nuTilda",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field adjnuTilda\n" << endl;

volScalarField adjnuTilda
(
    IOobject
    (
        "adjnuTilda",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading optimizationProperties\n" << endl;

IOdictionary optimizationProperties
(
    IOobject
    (
        "optimizationProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading transportProperties\n" << endl;
 
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info << "Create field adjSource\n" << endl;

volScalarField adjSource
(
    IOobject
    (
        "adjSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("adjSource", dimless, scalar(0.0))
);

// create field for adjoint sensitivities
Info << "Create field adjSensField\n" << endl;

// pointVectorField adjSensitivities
// (
//     IOobject
//     (
//         "adjSensitivities",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     pointMesh::New(mesh),
//     dimensionedVector("vector", dimless, vector::zero)
// );

pointScalarField adjSensField
(
    IOobject
    (
        "adjSensField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pointMesh::New(mesh),
    dimensionedScalar("adjSensField", dimless, scalar(0.0)),
    valuePointPatchField<scalar>::typeName
);

// create field for finite difference sensitivities
Info << "Create field fdSensField\n" << endl;

// pointVectorField FDSensitivities
// (
//     IOobject
//     (
//         "FDSensitivities",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     pointMesh::New(mesh),
//     dimensionedVector("vector", dimless, vector::zero)
// );

pointScalarField fdSensField
(
    IOobject
    (
        "fdSensField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pointMesh::New(mesh),
    dimensionedScalar("fdSensField", dimless, scalar(0.0)),
    valuePointPatchField<scalar>::typeName
);

// Needed to call nut.correctBoundaryConditions() in the primal problem? not sure why its needed?
singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

volVectorField dStar
(
    IOobject
    (
        "dStar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    vector(1, 1, 0)
);

// Computing field d
const volScalarField d = SpalartAllmaras::d(mesh);







// fields for checking
volVectorField checkRNS
(
    IOobject
    (
        "checkRNS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("checkRNS", dimAcceleration, vector(10, 10, 10))
);

volVectorField pRNSpNut
(
    IOobject
    (
        "pRNSpNut",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("pRNSpNut", dimensionSet(0, -1, -1, 0, 0, 0, 0), vector(10, 10, 10))
//     dimensionedVector("pRNSpNut", dimAcceleration, vector(10, 10, 10))    
);

volVectorField pRNSpNutFD
(
    IOobject
    (
        "pRNSpNutFD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("pRNSpNutFD", dimAcceleration, vector(10, 10, 10))
);

volScalarField pNutpNuTilda
(
    IOobject
    (
        "pNutpNuTilda",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("pNutpNuTilda", dimless, scalar(10.0))
);

volScalarField pNutpNuTildaFD
(
    IOobject
    (
        "pNutpNuTildaFD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("pNutpNuTildaFD", dimViscosity, scalar(10.0))
);