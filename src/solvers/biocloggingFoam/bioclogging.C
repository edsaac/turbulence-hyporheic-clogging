/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    laplacianFoam

Description
    Solves a simple Laplace equation, e.g. for thermal diffusion in a solid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

	simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	Info<< "\nCalculating...\n" << endl;

	while (simple.loop(runTime))
	{

		// Field of conductivity [L/T] from Kozeny-Carman eq.
    	K = (sqr(ds) * pow(theta,3) * rho * g) / (180.0 * pow(1-theta,2) * mu);

		Info<< "Time = " << runTime.timeName() << nl << endl;

		//Calculate hydraulic head using mass balance + Darcy's equation
		while (simple.correctNonOrthogonal())
		{
			fvScalarMatrix FlowEquation
			(
				fvm::laplacian(K, h)
				==
				fvc::ddt(theta)
			);
			fvOptions.constrain(FlowEquation);
			FlowEquation.solve();
			fvOptions.correct(h);
		}

		// U in this code refers to the Darcy velocity (a.k.a q)
		volVectorField U
		(
			IOobject
			(
				"U",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			-K * fvc::grad(h) + (Slope*K)
		);

		#include "createPhi.H"

		surfaceScalarField phi_set
		(
		     IOobject
		     (
		         "phi_set",
		         runTime.timeName(),
		         mesh,
		         IOobject::READ_IF_PRESENT,
		         IOobject::AUTO_WRITE
		     ),
		     fvc::flux(U + uSettling)
		 );

		#include "CourantNo.H"

		// Field of Happel parameter As
		volScalarField log_As
		(
			log(As_Coeffs[0]*sqr(1/theta) + As_Coeffs[1]*(1/theta) + As_Coeffs[2])
		);

		// Peclét number field
		volScalarField log_N_Peclet
		(
			log((mag(U)*ds)/D_diffusion)
		);

		// Field of filtration efficiency
		volScalarField eta
		(
			// Diffusion mechanism
			exp(
		      0.875468737
		    + (0.333 * log_As)
		    - (0.081 * log_N_Ratio)
		    - (0.715 * log_N_Peclet)
		    + (0.052 * log_N_vdW)
		    )

		    + 

		    // Interception mechanism
			exp(
		    - 0.597837001
		    + log_As
		    + (1.550 * log_N_Ratio)
		    - (0.125 * log_N_Peclet)
		    + (0.125 * log_N_vdW)
		    )

			    + 

		    // Gravitational deposition mechanism
			exp(
		    - 0.744440475
		    - (1.350 * log_N_Ratio)
		    - (1.110 * log_N_Peclet)
		    + (0.053 * log_N_vdW)
		    + (1.110 * log_N_gravit)
		    )

		);

		// Field of filtration coefficient
		volScalarField Lambda
		(
			IOobject
			(
				"Lambda",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
		    (3.0 * (1-theta) * alpha * eta)/(2.0*ds)
		);

		// Field of filtration coefficient
		Info<< "Calculate site blocking \n" << endl;
		PhiB = 1.0 - Cs/Smax;

		// Transport equations
		while (simple.correctNonOrthogonal())
		{
			fvScalarMatrix MovingStuffEquation
			(
				fvm::ddt(theta,Cw)
				+ fvm::div(phi_set, Cw)
				- fvm::laplacian(mag(U)*DispTensor, Cw)
				==
				- fvm::SuSp(Lambda*mag(U)*PhiB,Cw)
				+ (kdet*Cs)
			);
			MovingStuffEquation.relax();
			fvOptions.constrain(MovingStuffEquation);
			MovingStuffEquation.solve();
			fvOptions.correct(Cw);


			fvScalarMatrix StoppedStuffEquation
			(
				fvm::ddt(Cs)
				==
				- fvm::Sp(kdet,Cs)
				+ (Lambda*mag(U)*PhiB*Cw)
				
			);
			StoppedStuffEquation.relax();
			fvOptions.constrain(StoppedStuffEquation);
			StoppedStuffEquation.solve();
			fvOptions.correct(Cs);	
		}

		// Update porosity field
		theta = theta0 - Cs/rho_clay;

		// Calculate clay fraction
		percentClay = (Cs + (Cw*theta)) / (Cs + rho_sand*(1-theta0));

		//End bits
		runTime.write();

		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
		<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
		<< nl << endl;

	}

	Info<< "End\n" << endl;

	return 0;
}


// ************************************************************************* //
