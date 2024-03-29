/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "uqsinglePhaseTransportModel.H"
#include "uqviscosityModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(uqsinglePhaseTransportModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
/*
Foam::uqsinglePhaseTransportModel::uqsinglePhaseTransportModel
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    viscosityModelPtr_(uqviscosityModel::New("nu", *this, U, phi))
{}
*/

Foam::uqsinglePhaseTransportModel::uqsinglePhaseTransportModel
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word& nuNode
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    viscosityModelPtr_(uqviscosityModel::New(nuNode, *this, U, phi))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::uqsinglePhaseTransportModel::~uqsinglePhaseTransportModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::uqsinglePhaseTransportModel::nu() const
{
    return viscosityModelPtr_->nu();
}


Foam::tmp<Foam::scalarField>
Foam::uqsinglePhaseTransportModel::nu(const label patchi) const
{
    return viscosityModelPtr_->nu(patchi);
}


void Foam::uqsinglePhaseTransportModel::correct()
{
    viscosityModelPtr_->correct();
}


bool Foam::uqsinglePhaseTransportModel::read()
{
    if (regIOobject::read())
    {
        return viscosityModelPtr_->read(*this);
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
