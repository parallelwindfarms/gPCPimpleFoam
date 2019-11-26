scalar              tmpCijk     (0.0);
scalar              tmpCijkl    (0.0);
volScalarField      uqNut_k     (turbulence[0]->nut() * 0.0);
volScalarField      delta2      (sqr(turbulence[0]->getDelta()));

// LES Smag - Pseado Spectral Method 1
if ( PS1 )
{
    forAll(U, k)
    {
        uqNut_k *= 0.0;

        forAll(U, l)
        {
            forAll(U, m)
            {
                forAll(U, n)
                {
                    if ( (truncated && (l+m+n+k <= 3*P)) || (truncated == 0) )
                    {
                        tmpCijkl    = M4[l][m][n][k];
                        if ( tmpCijkl != 0 )
                        {
                            uqNut_k += tmpCijkl
                                     * turbulence[l]->getCs()*turbulence[m]->getCs()
                                     * delta2*S[n];
                        }
                    }
                }
            }
        }
        uqNut_k *= 1/M1[k];
        turbulence[k]->setNut(uqNut_k);
    }
}

// LES Smag - Pseado Spectral Method 2
else if ( PS2 )
{
    dimensionedScalar  tmpCs2_l("tmpCs2_l", dimless, 0);

    forAll(U, k)
    {
        uqNut_k *= 0.0;

        forAll(U, l)
        {
            tmpCs2_l *= 0.0;
            forAll(U, p)
            {
                forAll(U,q)
                {
                    if ( (truncated && (p+q+l <= 2*P)) || (truncated == 0) )
                    {
                        tmpCijk = M3[p][q][l];
                        if ( tmpCijk != 0 )
                        {
                            tmpCs2_l += tmpCijk
                                      * turbulence[p]->getCs()*turbulence[q]->getCs();
                        }
                    }
                }
            }
            tmpCs2_l *= 1/M1[l];

            forAll(U, m)
            {
                if ( (truncated && (l+m+k <= 2*P)) || (truncated == 0) )
                {
                    tmpCijk = M3[l][m][k];
                    if ( tmpCijk != 0 )
                    {
                        uqNut_k += tmpCijk
                                 * tmpCs2_l
                                 * delta2*S[m];
                    }
                }
            }
        }
        uqNut_k *= 1/M1[k];
        turbulence[k]->setNut(uqNut_k);
    }
}
