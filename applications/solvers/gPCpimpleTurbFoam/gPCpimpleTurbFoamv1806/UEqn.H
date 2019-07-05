// --- Solve the Momentum equation

Ckk = M1[k];

MRF.correctBoundaryVelocity(U[k]);

fvVectorMatrix UEqn (Ckk * fvm::ddt(U[k]));

forAll(U, j)
{
    forAll(U, i)
    {
        Cijk = M3[i][j][k];

        if ( j==k ) { UEqn += Cijk * fvm::div(phi[i], U[j]); }
        else        { UEqn += Cijk * fvc::div(phi[i], U[j]); }

        #include "uqNutCalc.H"

        if ( j==k ) { UEqn -= fvm::laplacian(uqNut_i+Cijk*turbulence[i]->nu(), U[j]); }
        else        { UEqn -= fvc::laplacian(uqNut_i+Cijk*turbulence[i]->nu(), U[j]); }

        UEqn -= fvc::div((uqNut_i+Cijk*turbulence[i]->nu())*dev2(T(fvc::grad(U[j]))));
    }
}

UEqn += MRF.DDt(U[k]);
UEqn == fvOptions(U[k]);

UEqn *= 1.0/Ckk;

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == - fvc::grad(p[k]));

    fvOptions.correct(U[k]);
}
