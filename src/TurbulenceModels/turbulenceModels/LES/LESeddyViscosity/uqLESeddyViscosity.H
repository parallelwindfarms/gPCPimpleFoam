/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::LESModels::uqLESeddyViscosity

Group
    grpLESTurbulence

Description
    Eddy viscosity LES SGS model base class

SourceFiles
    uqLESeddyViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef uqLESeddyViscosity_H
#define uqLESeddyViscosity_H

#include "uqLESModel.H"
#include "uqeddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                      Class uqLESeddyViscosity Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class uqLESeddyViscosity
:
    public uqeddyViscosity<uqLESModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        uqLESeddyViscosity(const uqLESeddyViscosity&) = delete;

        //- No copy assignment
        void operator=(const uqLESeddyViscosity&) = delete;


protected:

    // Protected data

        dimensionedScalar Ce_;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    // Constructors

        //- Construct from components
        uqLESeddyViscosity
        (
            const word& type,
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = uqturbulenceModel::propertiesName,
            const label& uqNode = 0
        );


    //- Destructor
    virtual ~uqLESeddyViscosity()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return sub-grid dissipation rate
        virtual tmp<volScalarField> epsilon() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "uqLESeddyViscosity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
