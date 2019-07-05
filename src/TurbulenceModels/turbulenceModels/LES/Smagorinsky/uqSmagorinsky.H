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
    Foam::LESModels::uqSmagorinsky

Group
    grpLESTurbulence

Description
    The uqSmagorinsky SGS model.

    Reference:
    \verbatim
        uqSmagorinsky, J. (1963).
        General circulation experiments with the primitive equations: I.
        The basic experiment*.
        Monthly weather review, 91(3), 99-164.
    \endverbatim

    The form of the uqSmagorinsky model implemented is obtained from the
    k-equation model assuming local equilibrium which provides estimates of both
    k and epsilon separate from the sub-grid scale viscosity:

    \verbatim
        B = 2/3*k*I - 2*nuSgs*dev(D)

    where

        D = symm(grad(U));
        k from D:B + Ce*k^3/2/delta = 0
        nuSgs = Ck*sqrt(k)*delta
    \endverbatim

    The default model coefficients are
    \verbatim
        uqSmagorinskyCoeffs
        {
            Ck                  0.094;
            Ce                  1.048;
        }
    \endverbatim

See also
    Foam::LESModels::kEqn

SourceFiles
    uqSmagorinsky.C

\*---------------------------------------------------------------------------*/

#ifndef uqSmagorinsky_H
#define uqSmagorinsky_H

#include "uqLESModel.H"
#include "uqLESeddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class uqSmagorinsky Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class uqSmagorinsky
:
    public uqLESeddyViscosity<BasicTurbulenceModel>
{
    // Private Member Functions

        //- No copy construct
        uqSmagorinsky(const uqSmagorinsky&) = delete;

        //- No copy assignment
        void operator=(const uqSmagorinsky&) = delete;


protected:

    // Protected data

        dimensionedScalar Ck_;

        dimensionedScalar Cs_;

    // Protected Member Functions

        //- Return SGS kinetic energy
        //  calculated from the given velocity gradient
        tmp<volScalarField> k(const tmp<volTensorField>& gradU) const;

        //- Update the SGS eddy viscosity
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("uqSmagorinsky");


    // Constructors

        uqSmagorinsky
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = uqturbulenceModel::propertiesName,
            const label& uqNode = 0,
            const word& type = typeName
        );

    //- Destructor
    virtual ~uqSmagorinsky()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k(fvc::grad(this->U_));
        }

        //- Return sub-grid dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Correct Eddy-Viscosity and related properties
        virtual void correct();

        //- Return Cs_ (needed for UQ)
        virtual dimensionedScalar getCs(){ return Cs_; }

        //- Return delta (needed for UQ)
        virtual volScalarField getDelta(){ return this->delta(); }

        //- Return model name
        virtual word getModelName(){ return this->typeName_(); }

        //- Set uqNut
        virtual void setNut(volScalarField uqNut) { this->nut_ = uqNut; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace uqLESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "uqSmagorinsky.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
