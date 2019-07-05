/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "uqincompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(uqincompressibleTurbulenceModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
/*
Foam::uqincompressibleTurbulenceModel::uqincompressibleTurbulenceModel
(
    const geometricOneField&,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const word& propertiesName
)
:
    uqturbulenceModel
    (
        U,
        alphaRhoPhi,
        phi,
        propertiesName
    )
{Info<< "????????????????????????" << endl;}
*/
Foam::uqincompressibleTurbulenceModel::uqincompressibleTurbulenceModel
(
    const geometricOneField&,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const word& propertiesName,
    const label& uqNode
)
:
    uqturbulenceModel
    (
        U,
        alphaRhoPhi,
        phi,
        propertiesName,
        uqNode
    )
{
  //Info<< "@uqincompressibleTurbulenceModel" << endl;
}

Foam::tmp<Foam::volScalarField>
Foam::uqincompressibleTurbulenceModel::mu() const
{
    return nu();
}


Foam::tmp<Foam::scalarField>
Foam::uqincompressibleTurbulenceModel::mu(const label patchi) const
{
    return nu(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::uqincompressibleTurbulenceModel::mut() const
{
    return nut();
}


Foam::tmp<Foam::scalarField>
Foam::uqincompressibleTurbulenceModel::mut(const label patchi) const
{
    return nut(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::uqincompressibleTurbulenceModel::muEff() const
{
    return nuEff();
}


Foam::tmp<Foam::scalarField>
Foam::uqincompressibleTurbulenceModel::muEff(const label patchi) const
{
    return nuEff(patchi);
}


// ************************************************************************* //
