// --- Calculate nuEff for UQ

volScalarField      uqNut_i     (turbulence[0]->nut() * 0.0);

volSymmTensorField  Di          (symm(fvc::grad(U[i])));
volScalarField      Si          (sqrt(2*(Di && Di)));
volScalarField      delta2      (sqr(turbulence[0]->getDelta()));
volScalarField      delta2S     (delta2 * Si);

if ( FS )
{
    forAll(U, l)
    {
        Cijkl    = M4[i][j][k][l];
        uqNut_i += Cijkl
                 * turbulence[i]->getCs()*turbulence[l]->getCs()
                 * delta2S;
    }
}














if ( PS )
{
    uqNut_i = Cijk * turbulence[i]->nut();
}













































/*
    volSymmTensorField  Dp(0.0*Di);
    volSymmTensorField  Dq(0.0*Di);
    volScalarField      Sipq(0.0*Si);

    forAll(U, p)
    {
        Dp = symm(fvc::grad(U[p]));

        forAll(U, q)
        {
            Dq      = symm(fvc::grad(U[q]));
            Sipq    += M3[i][p][q] * (Dp && Dq);
        }
    }

    Sipq    = sqrt(2*Sipq/M1[i]);
    delta2S = sqr(turbulence[0]->getDelta()) * Sipq;
//*/



/*
if ( uqSmagorinsky )
{
    volSymmTensorField  D(symm(fvc::grad(U[i])));
    volScalarField      S(sqrt(2*(D && D)));
    volScalarField      delta2S(sqr(turbulence[0]->getDelta()) * S);

    if ( FS )
    {
        forAll(U, l)
        {
          Cijkl = M4[i][j][k][l];

          uqNut_i += Cijkl * turbulence[i]->getCs()*turbulence[l]->getCs()*delta2S;
        }
        Info << max(uqNut_i) << endl;
    }

    if ( PS )
    {
        uqNut_i = sqr(turbulence[i]->getCs())*delta2S;
    }
}


if ( uqWALE )
{
    if ( FS || PS )
    {
        uqNut_i = turbulence[i]->nut();
    }
}
*/
