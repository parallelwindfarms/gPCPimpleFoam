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

Application
    pimpleFoam.C

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids
    on a moving mesh.

    \heading Solver details
    The solver uses the PIMPLE (merged PISO-SIMPLE) algorithm to solve the
    continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}} + \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p + \vec{S}_U
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
        \vec{S}_U | Momentum source
    \endvartable

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "uqsinglePhaseTransportModel.H"
#include "uqturbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //#include "postProcess.H"                  // Need some fixing

    #include "addCheckCaseOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "uqCreateUfIfPresent.H"
    #include "U0CourantNo.H"
    #include "setInitialDeltaT.H"

    forAll(U, k)
        turbulence[k]->validate();

    Info<< "\nStarting time loop\n" << endl;

    // --- Time loop (solving P+1 systems evey time iteration)
    while (runTime.run())
    {
        runTime++;

        Info<< "Time = " << runTime.timeName() << endl;

        #include "readDyMControls.H"
        #include "U0CourantNo.H"
        #include "setDeltaT.H"

        // --- Explicit iteration cycles b/w systems of Uk-pk
        while (expItr < expItrMax)
        {
            expItr++;

            Info<< nl << "expItr: cycle " << expItr << endl;

            // --- Loop over velocity modes (reverse)
            forAllReverse(U, k)
            {
                // --- Efficiently updating UQ modes
                if(expItrMax-expItr >= (P-k) || expItr <= k+1)
                {
                    Info<< "gPC: mode " << k << endl;

                    // --- Pressure-velocity PIMPLE corrector loop
                    while (pimple.loop())
                    {
                          #include "UkCourantNo.H"

                          // --- MRF calculations
                          #include "MRFcalculations.H"

                          // --- Momentum predictor
                          #include "UEqn.H"

                          // --- Pressure corrector loop
                          while (pimple.correct())
                          {
                              #include "pEqn.H"
                          }

                          // --- Update/Solve for turbulence
                          if (pimple.turbCorr())
                          {
                              laminarTransport[k]->correct();
                              //turbulence[k]->correct();
                          }

                    } // --- End of PIMPLE loop

                    Info<< "------------------------------------" << endl;

                }

            } // --- End of velocity modes loop

            #include "uqS_k.H"                          // Need some fixing
            #include "uqNut_k.H"

        } // --- End of explicit cycles


        // --- Calculating the mean and st.dev. for UQ
        #include "uqPostProcess.H"

        runTime.write();
        runTime.printExecutionTime(Info);

    } // --- End of time loop

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
