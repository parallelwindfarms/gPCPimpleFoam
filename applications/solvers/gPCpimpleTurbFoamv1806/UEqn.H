// --- Solve the Momentum equation

Ckk = M1[P-k];

MRF.correctBoundaryVelocity(U[P-k]);

fvVectorMatrix UEqn (Ckk * fvm::ddt(U[P-k]));

if ( FS )
{
    forAll(U, j)
    {
        forAll(U, i)
        {
            if ( (truncated && (i+j+P-k <= 2*P)) || (truncated != 1) )
            {
                Cijk = M3[i][j][P-k];

                if ( Cijk != 0 )
                {
                    if ( j==(P-k) ) { UEqn += Cijk * fvm::div(phi[i], U[j]); }
                    else            { UEqn += Cijk * fvc::div(phi[i], U[j]); }
                }
            }
            else { Cijk = 0.0; }

            #include "uqNutCalc.H"

            if ( j==(P-k) ) { UEqn -= fvm::laplacian(uqNut_i+Cijk*turbulence[i]->nu(), U[j]); }
            else            { UEqn -= fvc::laplacian(uqNut_i+Cijk*turbulence[i]->nu(), U[j]); }

            UEqn -= fvc::div((uqNut_i+Cijk*turbulence[i]->nu())*dev2(T(fvc::grad(U[j]))));
        }
    }
}

else if ( PS1 || PS2 )
{
    forAll(U, j)
    {
        forAll(U, i)
        {
            if ( (truncated && (i+j+P-k <= 2*P)) || (truncated != 1) )
            {
                Cijk = M3[i][j][P-k];

                if ( Cijk != 0 )
                {
                    if ( j==(P-k) )
                    {
                        UEqn += Cijk * fvm::div(phi[i], U[j]);
                        UEqn -= Cijk * fvm::laplacian(turbulence[i]->nuEff(), U[j]);
                    }
                    else
                    {
                        UEqn += Cijk * fvc::div(phi[i], U[j]);
                        UEqn -= Cijk * fvc::laplacian(turbulence[i]->nuEff(), U[j]);
                    }

                    UEqn -= Cijk * fvc::div((turbulence[i]->nuEff())*dev2(T(fvc::grad(U[j]))));
                }
            }
        }
    }
}

UEqn += MRF.DDt(U[P-k]);
UEqn == fvOptions(U[P-k]);

UEqn *= 1.0/Ckk;

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == - fvc::grad(p[P-k]));

    fvOptions.correct(U[P-k]);
}
