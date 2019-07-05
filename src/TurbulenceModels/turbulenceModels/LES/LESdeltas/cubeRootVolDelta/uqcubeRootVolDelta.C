/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "uqcubeRootVolDelta.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(uqcubeRootVolDelta, 0);
    addToRunTimeSelectionTable(uqLESdelta, uqcubeRootVolDelta, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::uqcubeRootVolDelta::calcDelta()
{
    const fvMesh& mesh = turbulenceModel_.mesh();

    label nD = mesh.nGeometricD();

    if (nD == 3)
    {
        delta_.primitiveFieldRef() = deltaCoeff_*pow(mesh.V(), 1.0/3.0);
    }
    else if (nD == 2)
    {
        WarningInFunction
            << "Case is 2D, LES is not strictly applicable\n"
            << endl;

        const Vector<label>& directions = mesh.geometricD();

        scalar thickness = 0.0;
        for (direction dir=0; dir<directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = mesh.bounds().span()[dir];
                break;
            }
        }

        delta_.primitiveFieldRef() = deltaCoeff_*sqrt(mesh.V()/thickness);
    }
    else
    {
        if (debug)
        {
            FatalErrorInFunction
                << "Case is not 3D or 2D, LES is not applicable"
                << exit(FatalError);
        }
    }

    // Handle coupled boundaries
    delta_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::uqcubeRootVolDelta::uqcubeRootVolDelta
(
    const word& name,
    const uqturbulenceModel& turbulence,
    const dictionary& dict
)
:
    uqLESdelta(name, turbulence),
    deltaCoeff_
    (
        dict.optionalSubDict(type() + "Coeffs").lookupOrDefault<scalar>
        (
            "deltaCoeff",
            1
        )
    )
{
    //Info<< "@uqcubeRootVolDelta" << endl;
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::uqcubeRootVolDelta::read(const dictionary& dict)
{
    dict.optionalSubDict(type() + "Coeffs").readIfPresent<scalar>
    (
        "deltaCoeff",
        deltaCoeff_
    );

    calcDelta();
}


void Foam::LESModels::uqcubeRootVolDelta::correct()
{
    if (turbulenceModel_.mesh().changing())
    {
        calcDelta();
    }
}


// ************************************************************************* //
