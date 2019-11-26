// --- Calculate nuEff for UQ

scalar              tmpCijklm   (0.0);
volScalarField      uqNut_FS    (turbulence[0]->nut() * 0.0);
volScalarField      delta2      (sqr(turbulence[0]->getDelta()));

if ( FS )
{
    forAll(U, l)
    {
        forAll(U, m)
        {
            if ( (truncated && (i+j+k+l+m <= 4*P)) || (truncated != 1) )
            {
                tmpCijklm = M5[k][j][i][l][m];
                if ( tmpCijklm != 0 )
                {
                    uqNut_FS += tmpCijklm
                              * turbulence[l]->getCs()*turbulence[m]->getCs()
                              * delta2*S[i];
                }
            }

        }

    }
}
/*
else if ( PS1 )
{
    forAll(U, l)
    {
        forAll(U, m)
        {
            forAll(U, n)
            {
                tmpCijk = M4[l][m][n][i];
                if ( tmpCijk != 0 )
                {
                    uqNut_i += tmpCijk
                             * turbulence[l]->getCs()*turbulence[m]->getCs()
                             * delta2*S[n];
                }
            }
        }
    }
    uqNut_i *= 1/M1[i];
}

else if ( PS2 )
{
    forAll(U, l)
    {
        forAll(U, m)
        {
            tmpCijk = M3[l][m][i];
            if ( tmpCijk != 0 )
            {
                uqNut_i += tmpCijk
                         * turbulence[l]->getCs()*turbulence[l]->getCs()
                         * delta2*S[m];
            }
        }
    }
    uqNut_i *= 1/M1[i];
}
*/



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
