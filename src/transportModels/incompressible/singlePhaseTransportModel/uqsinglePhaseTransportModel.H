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

Class
    Foam::uqsinglePhaseTransportModel

Description
    A simple single-phase transport model based on viscosityModel.

    Used by the incompressible single-phase solvers like simpleFoam,
    pimpleFoam etc.

SourceFiles
    uqsinglePhaseTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef uqsinglePhaseTransportModel_H
#define uqsinglePhaseTransportModel_H

//#include "incompressible/transportModel/transportModel.H"
#include "uqtransportModel.H"
#include "IOdictionary.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class uqviscosityModel;

/*---------------------------------------------------------------------------*\
                Class uqsinglePhaseTransportModel Declaration
\*---------------------------------------------------------------------------*/

class uqsinglePhaseTransportModel
:
    public IOdictionary,
    public uqtransportModel
{
    // Private Data

        autoPtr<uqviscosityModel> viscosityModelPtr_;


    // Private Member Functions

        //- No copy construct
        //uqsinglePhaseTransportModel(const uqsinglePhaseTransportModel&) = delete;

        //- No copy assignment
        //void operator=(const uqsinglePhaseTransportModel&) = delete;

public:

    //- Runtime type information
    TypeName("uqsinglePhaseTransportModel");

    // Constructors

        //- Construct from components
        uqsinglePhaseTransportModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );
        
        //- Construct from uq components
        uqsinglePhaseTransportModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word& nuNode
        );

    //- Destructor
    virtual ~uqsinglePhaseTransportModel();


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const;

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const;

        //- Correct the laminar viscosity
        virtual void correct();

        //- Read transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
