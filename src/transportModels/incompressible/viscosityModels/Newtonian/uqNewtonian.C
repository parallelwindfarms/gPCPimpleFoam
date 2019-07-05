/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd
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

#include "uqNewtonian.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace uqviscosityModels
{
    defineTypeNameAndDebug(uqNewtonian, 0);
    addToRunTimeSelectionTable(uqviscosityModel, uqNewtonian, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uqviscosityModels::uqNewtonian::uqNewtonian
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    uqviscosityModel(name, viscosityProperties, U, phi),
    nu0_(  viscosityProperties_.found(name)                            ?
           dimensionedScalar(name, dimViscosity, viscosityProperties_) :
           dimensionedScalar(name, dimViscosity, 0.0)//VSMALL)
        ),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh(),
        nu0_
    )
{
    //Info<< "@uqNewtonian" << endl;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::uqviscosityModels::uqNewtonian::read
(
    const dictionary& viscosityProperties
)
{
    //Info<< "@uqNewtonian::read" << endl;
    uqviscosityModel::read(viscosityProperties);

    viscosityProperties_.lookup(nu_.name()) >> nu0_;
    nu_ = nu0_;

    return true;
}


// ************************************************************************* //
