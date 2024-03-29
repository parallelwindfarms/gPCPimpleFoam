/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "uqLESdelta.H"
#include "calculatedFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(uqLESdelta, 0);
    defineRunTimeSelectionTable(uqLESdelta, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uqLESdelta::uqLESdelta
(
    const word& name,
    const uqturbulenceModel& turbulence
)
:
    turbulenceModel_(turbulence),
    delta_
    (
        IOobject
        (
            name,
            turbulence.mesh().time().timeName(),
            turbulence.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulence.mesh(),
        dimensionedScalar(name, dimLength, SMALL),
        calculatedFvPatchScalarField::typeName
    )
{
  //Info<< "@uqLESdelta" << endl;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::uqLESdelta> Foam::uqLESdelta::New
(
    const word& name,
    const uqturbulenceModel& turbulence,
    const dictionary& dict,
    const word& lookupName
)
{
    const word deltaType(dict.lookup(lookupName));

    Info<< "Selecting LES " << lookupName << " type " << deltaType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(deltaType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown LESdelta type "
            << deltaType << nl << nl
            << "Valid LESdelta types :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<uqLESdelta>(cstrIter()(name, turbulence, dict));
}


Foam::autoPtr<Foam::uqLESdelta> Foam::uqLESdelta::New
(
    const word& name,
    const uqturbulenceModel& turbulence,
    const dictionary& dict,
    const dictionaryConstructorTable& additionalConstructors,
    const word& lookupName
)
{
    const word deltaType(dict.lookup(lookupName));

    Info<< "Selecting LES " << lookupName << " type " << deltaType << endl;

    // First any additional ones
    {
        auto cstrIter = additionalConstructors.cfind(deltaType);

        if (cstrIter.found())
        {
            return autoPtr<uqLESdelta>(cstrIter()(name, turbulence, dict));
        }
    }

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(deltaType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown LESdelta type "
            << deltaType << nl << nl
            << "Valid LESdelta types :" << endl
            << additionalConstructors.sortedToc()
            << " and "
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<uqLESdelta>(cstrIter()(name, turbulence, dict));
}


// ************************************************************************* //
