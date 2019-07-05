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
    Foam::uqIncompressibleTurbulenceModel

Description
    Templated abstract base class for single-phase incompressible
    turbulence models.

SourceFiles
    uqIncompressibleTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef uqIncompressibleTurbulenceModel_H
#define uqIncompressibleTurbulenceModel_H

#include "uqTurbulenceModel.H"
#include "uqincompressibleTurbulenceModel.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class uqIncompressibleTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

template<class TransportModel>
class uqIncompressibleTurbulenceModel
:
    public uqTurbulenceModel
    <
        geometricOneField,
        geometricOneField,
        uqincompressibleTurbulenceModel,
        TransportModel
    >
{

public:

    typedef geometricOneField alphaField;
    typedef geometricOneField rhoField;
    typedef TransportModel transportModel;


    // Constructors

        //- Construct

/*        uqIncompressibleTurbulenceModel
        (
            const word& type,
            const geometricOneField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const TransportModel& transport,
            const word& propertiesName
        );
*/
        uqIncompressibleTurbulenceModel
        (
            const word& type,
            const geometricOneField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const TransportModel& transport,
            const word& propertiesName,
            const label& uqNode
        );

    // Selectors

        //- Return a reference to the selected turbulence model
/*        static autoPtr<uqIncompressibleTurbulenceModel> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const TransportModel& transportModel,
            const word& propertiesName = turbulenceModel::propertiesName
        );
*/
        static autoPtr<uqIncompressibleTurbulenceModel> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const TransportModel& transportModel,
            const word& propertiesName = uqturbulenceModel::propertiesName,
            const label& uqNode = 0
        );

    //- Destructor
    virtual ~uqIncompressibleTurbulenceModel()
    {}


    // Member Functions

        //- Return the effective stress tensor
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the effective stress tensor
        virtual tmp<volSymmTensorField> devRhoReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "uqIncompressibleTurbulenceModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
