/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
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
    Foam::uqTurbulenceModel

Description
    Templated abstract base class for turbulence models

SourceFiles
    uqTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef uqTurbulenceModel_H
#define uqTurbulenceModel_H

#include "uqturbulenceModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class uqTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

template
<
    class Alpha,
    class Rho,
    class BasicTurbulenceModel,
    class TransportModel
>
class uqTurbulenceModel
:
    public BasicTurbulenceModel
{

public:

    typedef Alpha alphaField;
    typedef Rho rhoField;
    typedef TransportModel transportModel;


protected:

    // Protected data

        const alphaField& alpha_;
        const transportModel& transport_;


private:

    // Private Member Functions

        //- No copy construct
        uqTurbulenceModel(const uqTurbulenceModel&) = delete;

        //- No copy assignment
        void operator=(const uqTurbulenceModel&) = delete;


public:

    // Declare run-time constructor selection table
/*
        declareRunTimeNewSelectionTable
        (
            autoPtr,
            uqTurbulenceModel,
            dictionary,
            (
                const alphaField& alpha,
                const rhoField& rho,
                const volVectorField& U,
                const surfaceScalarField& alphaRhoPhi,
                const surfaceScalarField& phi,
                const transportModel& transport,
                const word& propertiesName
            ),
            (alpha, rho, U, alphaRhoPhi, phi, transport, propertiesName)
        );

*/
        declareRunTimeNewSelectionTable
        (
            autoPtr,
            uqTurbulenceModel,
            dictionary,
            (
                const alphaField& alpha,
                const rhoField& rho,
                const volVectorField& U,
                const surfaceScalarField& alphaRhoPhi,
                const surfaceScalarField& phi,
                const transportModel& transport,
                const word& propertiesName,
                const label& uqNode
            ),
            (alpha, rho, U, alphaRhoPhi, phi, transport, propertiesName, uqNode)
        );

    // Constructors

        //- Construct
/*        uqTurbulenceModel
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName
        );
*/
        uqTurbulenceModel
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName,
            const label& uqNode
        );

    // Selectors

        //- Return a reference to the selected turbulence model
/*        static autoPtr<uqTurbulenceModel> New
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName
        );
*/
        static autoPtr<uqTurbulenceModel> New
        (
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
    virtual ~uqTurbulenceModel()
    {}


    // Member Functions

        //- Access function to phase fraction
        const alphaField& alpha() const
        {
            return alpha_;
        }

        //- Access function to incompressible transport model
        const transportModel& transport() const
        {
            return transport_;
        }

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return transport_.nu();
        }

        //- Return the laminar viscosity on patchi
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return transport_.nu(patchi);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "uqTurbulenceModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
