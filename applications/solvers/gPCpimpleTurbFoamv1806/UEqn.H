// --- Solve the Momentum equation

Ckk = M1[k];

MRF.correctBoundaryVelocity(U[k]);

fvVectorMatrix UEqn (Ckk * fvm::ddt(U[k]));

if ( FS )
{
    forAll(U, j)
    {
        forAll(U, i)
        {
            if ( (truncated && (i+j+k <= 2*P)) || (truncated != 1) )
            {
                Cijk = M3[k][j][i];

                if ( Cijk != 0.0 )
                {
                    if ( j==k )
                        { UEqn += Cijk * fvm::div(phi[i], U[j]); }
                    //else if ( j!=k && i==k )
                    //    { UEqn += Cijk * fvm::div(phi[j], U[i]); } // fvm -> Sp. Convg.
                    else
                        { UEqn += Cijk * fvc::div(phi[i], U[j]); }
                }
            }
            else { Cijk = 0.0; }

            #include "uqNut_FS.H"

            /*/volScalarField uqNut (uqNut_FS+turbulence[i]->nu());
            volScalarField uqNut (uqNut_FS);
            Info<< "kji " << k << j << i << "\tuqNut_i"
                << "\tavg "   << uqNut.internalField().weightedAverage(mesh.V()).value()
                << "\tmax " << max(uqNut.internalField()).value()
                << "\tmin " << min(uqNut.internalField()).value()
                << endl;//*/

            if ( j==k ) { UEqn -= fvm::laplacian(uqNut_FS+Cijk*turbulence[i]->nu(), U[j]); }
            else        { UEqn -= fvc::laplacian(uqNut_FS+Cijk*turbulence[i]->nu(), U[j]); }

            UEqn -= fvc::div((uqNut_FS+Cijk*turbulence[i]->nu())*dev2(T(fvc::grad(U[j]))));
        }
    }
}

else if ( PS1 || PS2 )
{
    forAll(U, j)
    {
        forAll(U, i)
        {
            if ( (truncated && (i+j+k <= 2*P)) || (truncated != 1) )
            {
                Cijk = M3[k][j][i];

                if ( Cijk != 0.0 )
                {
                    if ( j==k )
                    {
                        UEqn += Cijk * fvm::div(phi[i], U[j]);
                        UEqn -= Cijk * fvm::laplacian(turbulence[i]->nuEff(), U[j]);
                    }
                    /*else if ( j!=k && i==k )
                    {
                        UEqn += Cijk * fvm::div(phi[j], U[i]); // fvm -> Sp. Convg.
                        UEqn -= Cijk * fvc::laplacian(turbulence[i]->nuEff(), U[j]);
                    }*/
                    else
                    {
                        UEqn += Cijk * fvc::div(phi[i], U[j]);
                        UEqn -= Cijk * fvc::laplacian(turbulence[i]->nuEff(), U[j]);
                    }

                    UEqn -= Cijk * fvc::div((turbulence[i]->nuEff())*dev2(T(fvc::grad(U[j]))));

                    /*/volScalarField uqNut (turbulence[i]->nuEff());
                    volScalarField uqNut (turbulence[i]->nut());
                    Info<< "kji " << k << j << i << "\tuqNut_i"
                        << "\tavg "<< uqNut.weightedAverage(mesh.V()).value()
                        << "\tmax " << max(uqNut).value()
                        << "\tmin " << min(uqNut).value()
                        << endl;//*/
                }
            }
        }
    }
}

UEqn += MRF.DDt(U[k]);

// Explicit pressure forcing for periodicity
if ( k==0 ){ UEqn += Ckk*GradP; }

UEqn == fvOptions(U[k]);

UEqn *= 1.0/Ckk;

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == - fvc::grad(p[k]));

    fvOptions.correct(U[k]);
}
