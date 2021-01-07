Info << "Reading field phi" << endl;

volScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field adjPhi" << endl;

volScalarField adjPhi
(
    IOobject
    (
        "adjPhi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Create field adjSource" << endl;

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
    dimensionedScalar("adjSource", dimless, 0.0)
);

Info << "Create field adjSens" << endl;

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
    dimensionedScalar("adjSensField", dimless, 0.0),
    valuePointPatchField<scalar>::typeName
);

Info << "Create field FDSensitivities" << endl;

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
    dimensionedScalar("fdSensField", dimless, 0.0),
    valuePointPatchField<scalar>::typeName
);
