/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
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

Typedef
    Foam::incompressible::turbulenceModel

Typedef
    Foam::incompressible::RASModel

Typedef
    Foam::incompressible::LESModel

Description
    Typedefs for turbulence, RAS and LES models for incompressible flow
    based on the standard laminar transport package.

SourceFiles
    uqturbulentTransportModel.C
    turbulentTransportModels.C

\*---------------------------------------------------------------------------*/

#ifndef uqturbulentTransportModel_H
#define uqturbulentTransportModel_H

#include "uqIncompressibleTurbulenceModel.H"
//#include "laminarModel.H"
//#include "RASModel.H"
//#include "uqLESModel.H"
//#include "incompressible/transportModel/transportModel.H"
#include "uqtransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace uqincompressible
    {
        typedef uqIncompressibleTurbulenceModel<uqtransportModel> uqturbulenceModel;

        //typedef laminarModel<turbulenceModel> laminarModel;
        //typedef RASModel<turbulenceModel> RASModel;
        //typedef uqLESModel<uqturbulenceModel> LESModel;

/*        template<class BasicCompressibleTurbulenceModel>
        autoPtr<BasicCompressibleTurbulenceModel> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const typename BasicCompressibleTurbulenceModel::transportModel&
                transport,
            const word& propertiesName = uqturbulenceModel::propertiesName
        );
*/
        template<class BasicCompressibleTurbulenceModel>
        autoPtr<BasicCompressibleTurbulenceModel> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const typename BasicCompressibleTurbulenceModel::uqtransportModel&
                transport,
            const word& propertiesName = uqturbulenceModel::propertiesName,
            const label& uqNode = 0
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "uqturbulentTransportModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
