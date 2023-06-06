/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude
#line 83 "/home/edsaac/Repos/turbulence-hyporheic-clogging/subsurface/simulations/turbulent-pressure/0/h.boundaryField.top"
#include "fvCFD.H"
            #include <iostream>
            #include <math.h>
            #include "calc_signal.H"
//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = aaf11e018367024fdb3ee7c28c1d65792717caab
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void sineWaves_aaf11e018367024fdb3ee7c28c1d65792717caab(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    sineWavesFixedValueFvPatchScalarField
);


const char* const sineWavesFixedValueFvPatchScalarField::SHA1sum =
    "aaf11e018367024fdb3ee7c28c1d65792717caab";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sineWavesFixedValueFvPatchScalarField::
sineWavesFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct sineWaves sha1: aaf11e018367024fdb3ee7c28c1d65792717caab"
            " from patch/DimensionedField\n";
    }
}


sineWavesFixedValueFvPatchScalarField::
sineWavesFixedValueFvPatchScalarField
(
    const sineWavesFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct sineWaves sha1: aaf11e018367024fdb3ee7c28c1d65792717caab"
            " from patch/DimensionedField/mapper\n";
    }
}


sineWavesFixedValueFvPatchScalarField::
sineWavesFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct sineWaves sha1: aaf11e018367024fdb3ee7c28c1d65792717caab"
            " from patch/dictionary\n";
    }
}


sineWavesFixedValueFvPatchScalarField::
sineWavesFixedValueFvPatchScalarField
(
    const sineWavesFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info<<"construct sineWaves sha1: aaf11e018367024fdb3ee7c28c1d65792717caab"
            " as copy\n";
    }
}


sineWavesFixedValueFvPatchScalarField::
sineWavesFixedValueFvPatchScalarField
(
    const sineWavesFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct sineWaves sha1: aaf11e018367024fdb3ee7c28c1d65792717caab "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

sineWavesFixedValueFvPatchScalarField::
~sineWavesFixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy sineWaves sha1: aaf11e018367024fdb3ee7c28c1d65792717caab\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void sineWavesFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs sineWaves sha1: aaf11e018367024fdb3ee7c28c1d65792717caab\n";
    }

//{{{ begin code
    #line 51 "/home/edsaac/Repos/turbulence-hyporheic-clogging/subsurface/simulations/turbulent-pressure/0/h.boundaryField.top"
// Call geometry
            const fvPatch& boundaryPatch = patch(); 
            // Initialize field
            const vectorField& Cf = boundaryPatch.Cf(); 
            scalarField& field = *this; 

            // Call time
            scalar t = this->db().time().value();
            scalar X = 0;

            forAll(Cf, faceI)
            {
                X = Cf[faceI].x();
                // head = p_rho*toCentimeters/g
                field[faceI] = signalBuild(X,t);
                field[faceI] /= 9.81;

            }
//}}} end code

    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

