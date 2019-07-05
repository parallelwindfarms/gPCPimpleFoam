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

Class
    Foam::uqvanDriestDelta

Description
    Simple cube-root of cell volume delta used in incompressible LES models.

SourceFiles
    uqvanDriestDelta.C

\*---------------------------------------------------------------------------*/

#ifndef uqvanDriestDelta_H
#define uqvanDriestDelta_H

#include "uqLESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class uqvanDriestDelta Declaration
\*---------------------------------------------------------------------------*/

class uqvanDriestDelta
:
    public uqLESdelta
{
    // Private data

        autoPtr<uqLESdelta> geometricDelta_;
        scalar kappa_;
        scalar Aplus_;
        scalar Cdelta_;
        label calcInterval_;


    // Private Member Functions

        //- No copy construct
        uqvanDriestDelta(const uqvanDriestDelta&) = delete;

        //- No copy assignment
        void operator=(const uqvanDriestDelta&) = delete;

        // Calculate the delta values
        void calcDelta();


public:

    //- Runtime type information
    TypeName("uqvanDriest");


    // Constructors

        //- Construct from name, turbulenceModel and dictionary
        uqvanDriestDelta
        (
            const word& name,
            const uqturbulenceModel& turbulence,
            const dictionary&
        );


    //- Destructor
    virtual ~uqvanDriestDelta()
    {}


    // Member Functions

        //- Read the uqLESdelta dictionary
        virtual void read(const dictionary&);

        // Correct values
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
