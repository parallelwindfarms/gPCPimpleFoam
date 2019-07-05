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
    Foam::LESdelta

Description
    Abstract base class for LES deltas

SourceFiles
    LESdelta.C

\*---------------------------------------------------------------------------*/

#ifndef uqLESdelta_H
#define uqLESdelta_H

#include "uqturbulenceModel.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class uqLESdelta Declaration
\*---------------------------------------------------------------------------*/

class uqLESdelta
{

protected:

    // Protected data

        const uqturbulenceModel& turbulenceModel_;

        volScalarField delta_;


    // Protected Member Functions

        //- No copy construct
        uqLESdelta(const uqLESdelta&) = delete;

        //- No copy assignment
        void operator=(const uqLESdelta&) = delete;


public:

    //- Runtime type information
    TypeName("uqLESdelta");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            uqLESdelta,
            dictionary,
            (
                const word& name,
                const uqturbulenceModel& turbulence,
                const dictionary& dict
            ),
            (name, turbulence, dict)
        );


    // Constructors

        //- Construct from name, uqturbulenceModel and dictionary
        uqLESdelta
        (
            const word& name,
            const uqturbulenceModel& turbulence
        );


    // Selectors

        //- Return a reference to the selected LES delta
        static autoPtr<uqLESdelta> New
        (
            const word& name,
            const uqturbulenceModel& turbulence,
            const dictionary& dict,
            const word& lookupName = "delta"
        );

        //- Return a reference to the selected LES delta
        static autoPtr<uqLESdelta> New
        (
            const word& name,
            const uqturbulenceModel& turbulence,
            const dictionary& dict,
            const dictionaryConstructorTable& additionalConstructors,
            const word& lookupName = "delta"
        );


    //- Destructor
    virtual ~uqLESdelta()
    {}


    // Member Functions

        //- Return uqturbulenceModel reference
        const uqturbulenceModel& turbulence() const
        {
            return turbulenceModel_;
        }

        //- Read the uqLESdelta dictionary
        virtual void read(const dictionary&) = 0;

        // Correct values
        virtual void correct() = 0;


    // Member Operators

        virtual operator const volScalarField&() const
        {
            return delta_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
