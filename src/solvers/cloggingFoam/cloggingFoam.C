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
    cloggingFoam

Description
    Fine particle filtration

    Clogging is introduced with Kozeny-Carman relationship affecting permeability
    k = max(k0*(n/no)^3*((1-n0)/(1-n))^2, kmin*k0)

    Solves two species:
    - Suspended particles (C) {C}
    - Deposited particles (S) {S}

    Eqs:
    d(nC)/dt = - n*katt*C + kdet*C
    d(S)/dt  =   n*katt*C - kdet*C

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

#include "cloggingModel.H"
#include "attachmentModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Foam::Info << "\nCalculating...\n" << endl;

    while (simple.loop(runTime))
    {
        Foam::Info<< "Time = " << runTime.timeName() << nl << endl;

        //Info << "\nUpdate clog space limitation" << endl;
        clogLimiter = 1.0 - depositedClay/XMAX;
        
        Foam::Info << "\n Update porosity field (n):" << endl;
        n  = clogging->nRef() - depositedClay/rho_clog;

        //Calculate hydraulic head using mass balance + Darcy's equation
        if (cloggingSwitch) 
        { 
            clogging->calcPerm();
        }

        hydraulicCond = perm * WATERDENSITY * GRAVITY / WATERDYNAMICVISCOSITY;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix FlowEquation
            (
                fvm::laplacian(hydraulicCond, h)
            );
            fvOptions.constrain(FlowEquation);
            FlowEquation.solve();
            fvOptions.correct(h);
        }

        // Update flow field
        U   = - hydraulicCond * fvc::grad(h);
        phi = fvc::flux(U);

        #include "CourantNo.H"
        
        // Transport equations
        while (simple.correctNonOrthogonal())
        {
            // Calculate attachment rate
            attachment->calcAttachment();

            //Info << "\nDeposited fine particles (depositedClay)" << endl;
            fvScalarMatrix FiltratedEq
            (
                fvm::ddt(depositedClay)
                ==
                  n * katt * clogLimiter * suspendedClay
                - kdet * depositedClay
            );
            fvOptions.constrain(FiltratedEq);
            FiltratedEq.solve();
            fvOptions.correct(depositedClay);

            //Info << "\nSuspended fine particles" << endl;
            fvScalarMatrix SuspendedEq
            (
                fvm::ddt(n,suspendedClay)
                + fvm::div(phi, suspendedClay)
                - fvm::laplacian(n*(mag(U)*DispTensor + molDiff), suspendedClay)
                ==
                - n * katt * clogLimiter * suspendedClay 
                + kdet * depositedClay
            );
            fvOptions.constrain(SuspendedEq);
            SuspendedEq.solve();
            fvOptions.correct(suspendedClay);
        }

        //End bits
        runTime.write();

        Foam::Info << "ExecutionTime = " << runTime.elapsedCpuTime()   << " s"
                   << "    ClockTime = " << runTime.elapsedClockTime() << " s"
                   << nl << endl;

        //break;
    }

    Foam::Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
