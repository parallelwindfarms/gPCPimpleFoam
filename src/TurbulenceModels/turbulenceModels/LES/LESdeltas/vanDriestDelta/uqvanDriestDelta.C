/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "uqvanDriestDelta.H"
#include "wallFvPatch.H"
#include "wallDistData.H"
#include "wallPointYPlus.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(uqvanDriestDelta, 0);
    addToRunTimeSelectionTable(uqLESdelta, uqvanDriestDelta, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::uqvanDriestDelta::calcDelta()
{
    const fvMesh& mesh = turbulenceModel_.mesh();

    const volVectorField& U = turbulenceModel_.U();
    const tmp<volScalarField> tnu = turbulenceModel_.nu();
    const volScalarField& nu = tnu();
    tmp<volScalarField> nuSgs = turbulenceModel_.nut();

    volScalarField ystar
    (
        IOobject
        (
            "ystar",
            mesh.time().constant(),
            mesh
        ),
        mesh,
        dimensionedScalar("ystar", dimLength, GREAT)
    );

    const fvPatchList& patches = mesh.boundary();
    volScalarField::Boundary& ystarBf = ystar.boundaryFieldRef();

    forAll(patches, patchi)
    {
        if (isA<wallFvPatch>(patches[patchi]))
        {
            const fvPatchVectorField& Uw = U.boundaryField()[patchi];
            const scalarField& nuw = nu.boundaryField()[patchi];
            const scalarField& nuSgsw = nuSgs().boundaryField()[patchi];

            ystarBf[patchi] =
                    (nuw + VSMALL)/sqrt((nuw + nuSgsw + VSMALL)*mag(Uw.snGrad()) + VSMALL);
                    //nuw/sqrt((nuw + nuSgsw)*mag(Uw.snGrad()) + VSMALL);
        }
    }

    scalar cutOff = wallPointYPlus::yPlusCutOff;
    wallPointYPlus::yPlusCutOff = 500;
    wallDistData<wallPointYPlus> y(mesh, ystar);
    wallPointYPlus::yPlusCutOff = cutOff;

    delta_ = min
    (
        static_cast<const volScalarField&>(geometricDelta_()),
        (kappa_/Cdelta_)*((scalar(1) + SMALL) - exp(-y/ystar/Aplus_))*y
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::uqvanDriestDelta::uqvanDriestDelta
(
    const word& name,
    const uqturbulenceModel &turbulence,
    const dictionary& dict
)
:
    uqLESdelta(name, turbulence),
    geometricDelta_
    (
        uqLESdelta::New
        (
            IOobject::groupName("geometricDelta", turbulence.U().group()),
            turbulence,
            // Note: cannot use optionalSubDict - if no *Coeffs dict the
            // code will get stuck in a loop attempting to read the delta entry
            // - consider looking up "geometricDelta" instead of "delta"?
            dict.subDict(type() + "Coeffs")
        )
    ),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    Aplus_
    (
        dict.optionalSubDict(type() + "Coeffs").lookupOrDefault<scalar>
        (
            "Aplus",
            26.0
        )
    ),
    Cdelta_
    (
        dict.optionalSubDict(type() + "Coeffs").lookupOrDefault<scalar>
        (
            "Cdelta",
            0.158
        )
    ),
    calcInterval_
    (
        dict.optionalSubDict(type() + "Coeffs").lookupOrDefault<label>
        (
            "calcInterval",
            1
        )
    )
{
    delta_ = geometricDelta_();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::uqvanDriestDelta::read(const dictionary& dict)
{
    const dictionary& coeffsDict(dict.optionalSubDict(type() + "Coeffs"));

    geometricDelta_().read(coeffsDict);
    dict.readIfPresent<scalar>("kappa", kappa_);
    coeffsDict.readIfPresent<scalar>("Aplus", Aplus_);
    coeffsDict.readIfPresent<scalar>("Cdelta", Cdelta_);
    coeffsDict.readIfPresent<label>("calcInterval", calcInterval_);

    calcDelta();
}


void Foam::LESModels::uqvanDriestDelta::correct()
{
    if (turbulenceModel_.mesh().time().timeIndex() % calcInterval_ == 0)
    {
        geometricDelta_().correct();
        calcDelta();
    }
}


// ************************************************************************* //
