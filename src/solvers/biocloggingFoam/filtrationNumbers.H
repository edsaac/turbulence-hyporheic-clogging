// Happel model for filtration

//   Check HappelPolyFit.ipynb on how this replaces the
//   original Happel expression. 
float As_Coeffs[3] = { 8.99992117, -7.49203318, 0.4119361 };

//Dimensionless numbers for filtration model
// Aspect Ratio
dimensionedScalar log_N_Ratio
(
    log(dp/ds)
);
Info<< "Aspect Ratio:"<< log_N_Ratio.value() <<"\n" << endl;

// van der Waals 
dimensionedScalar log_N_vdW
( 
    log(Hamaker/(kBoltz*Temp))
);
Info<< "van der Waals:"<< log_N_vdW.value() <<"\n" << endl;

// Gravitational 
dimensionedScalar log_N_gravit 
( 
    log(
    ppi * pow(dp,4) * (rho_clay - rho) * g / ( 12.0 * kBoltz * Temp )
    )

);
Info<< "Gravitational:"<< log_N_gravit.value() <<"\n" << endl;