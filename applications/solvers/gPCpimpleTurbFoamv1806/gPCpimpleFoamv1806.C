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

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    bool restarted1(false);
    bool perfectRestart1(false);
    bool restarted2(false);
    bool perfectRestart2(false);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    // --- Time loop (solving P+1 systems evey time iteration)
    while (runTime.run())
    {
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readDyMControls.H"
        #include "setDeltaT.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            // --- Loop over velocity modes
            forAll(U, k)
            {
                  Info<< "gPC: node " << k << endl;

                  #include "uqCourantNo.H"

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
                      turbulence[k]->correct();
                  }

                  if(k<P)
                      Info<< "------------------------------------" << endl;
                  if(k==P)
                      Info<< endl;

            } // --- End of velocity modes loop

        } // --- End of PIMPLE loop

        // --- Calculating the mean and st.dev. for UQ
        #include "uqPostProcess.H"
        #include "calcRk.H"

        runTime.write();
        runTime.printExecutionTime(Info);

    } // --- End of time loop

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
