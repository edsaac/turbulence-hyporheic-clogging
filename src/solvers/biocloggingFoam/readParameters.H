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

// Diameters
Info<< "Reading bed grain size (ds)\n" << endl;
dimensionedScalar ds(transportProperties.lookup("ds"));

Info<< "Reading clay size (dp)\n" << endl;
dimensionedScalar dp(transportProperties.lookup("dp"));

// Densities
Info<< "Reading water density (rho)\n" << endl;
dimensionedScalar rho(transportProperties.lookup("rho"));

Info<< "Reading clay density (rho_clay)\n" << endl;
dimensionedScalar rho_clay(transportProperties.lookup("rho_clay"));

Info<< "Reading sand density (rho_sand)\n" << endl;
dimensionedScalar rho_sand(transportProperties.lookup("rho_sand"));

//Other physical stuff
Info<< "Reading gravity (g)\n" << endl;
dimensionedScalar g(transportProperties.lookup("g"));

Info<< "Reading dynamic viscosity (mu)\n" << endl;
dimensionedScalar mu(transportProperties.lookup("mu"));

Info<< "Reading Slope\n" << endl;
dimensionedVector Slope(transportProperties.lookup("Slope"));

//Dispersion
Info<< "Dispersion tensor\n" << endl;
dimensionedTensor DispTensor(transportProperties.lookup("DispTensor"));

//Filtration & clogging
Info<< "Reading attachement efficiency (alpha)\n" << endl;
dimensionedScalar alpha(transportProperties.lookup("alpha"));

Info<< "Reading min porosity (thetaMin)\n" << endl;
dimensionedScalar thetaMin(transportProperties.lookup("thetaMin"));

Info<< "Reading settling velocity (uSettling)\n" << endl;
dimensionedVector uSettling(transportProperties.lookup("uSettling"));

//Chemical parameters
Info<< "Reading Hamaker constant (Hamaker)\n" << endl;
dimensionedScalar Hamaker(transportProperties.lookup("Hamaker"));

Info<< "Reading Boltzmann constan (kBoltz)\n" << endl;
dimensionedScalar kBoltz(transportProperties.lookup("kBoltz"));

Info<< "Reading Temperature (Temp)\n" << endl;
dimensionedScalar Temp(transportProperties.lookup("Temp"));

Info<< "Detachment first order rate\n" << endl;
dimensionedScalar kdet(transportProperties.lookup("kdet"));

