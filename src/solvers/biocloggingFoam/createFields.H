
const float ppi = 3.14159264; 

Info<< "Reading field: Hydraulic head (h) \n" << endl;
volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

fvScalarMatrix InitFlowEquation
(
    fvm::laplacian(h)
);
    InitFlowEquation.solve();

Info<< "Reading field: Porosity (theta) \n" << endl;
volScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Store a copy of the initial porosity distribution
volScalarField theta0(theta);

Info<< "Reading field: Suspended clay (Cw) \n" << endl;
volScalarField Cw
(
    IOobject
    (
        "Cw",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field: Stopped clay (Cs) \n" << endl;
volScalarField Cs
(
    IOobject
    (
        "Cs",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Read the transport properties dictionary
#include "readParameters.H"

// Calculate diffusion coefficient from Einstein-Stokes equation
//   >> use this in the Péclet number calculation
dimensionedScalar D_diffusion
(
    (kBoltz*Temp)/(3.0*ppi*mu*dp)
);

// Calculate constant nondimensional filtration numbers
#include "filtrationNumbers.H"

volScalarField K 
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (sqr(ds) * pow(theta,3) * rho * g) / (180.0 * pow(1-theta,2) * mu)
);


volScalarField percentClay
(
    IOobject
    (
        "percentClay",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Calculate site maximum deposited mass
volScalarField Smax
(
    IOobject
    (   
        "Smax",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho_clay * (theta0 - thetaMin)
);

volScalarField PhiB
(
    IOobject
    (
        "PhiB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    1.0 - Cs/Smax
);

#include "createFvOptions.H"
