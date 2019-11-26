scalar              tmpCijk     (0.0);
scalar              tmpCijkl    (0.0);
volScalarField      uqNut_k     (turbulence[0]->nut() * 0.0);
volScalarField      delta2      (sqr(turbulence[0]->getDelta()));

// LES Smag - Full Spectral or Pseado Spectral Method 1
if ( FS || PS1 )
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
                        tmpCijkl = M4[k][l][m][n];
                        if ( tmpCijkl != 0.0 )
                        {
                            uqNut_k += tmpCijkl
                                     * turbulence[l]->getCs() * turbulence[m]->getCs()
                                     * delta2*S[n];
                        }
                    }
                }
            }
        }
        uqNut_k *= 1.0/M1[k];
        turbulence[k]->uqCorrect(uqNut_k);

        volScalarField uqNut (turbulence[k]->nut());
        Info<< "uqNut_" << k
            << "\tavg " << uqNut.weightedAverage(mesh.V()).value()
            << "\tmax " << max(uqNut).value()
            << "\tmin " << min(uqNut).value()
            << endl;//*/
    }
    Info<< endl;
}

// LES Smag - Psuado Spectral Method 2
//else if ( PS2 )
else if ( 0 )
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
                        tmpCijk = M3[l][p][q];
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
                    tmpCijk = M3[k][l][m];
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
        turbulence[k]->uqCorrect(uqNut_k);
    }
}


// LES Smag - Psuado Spectral Method 2 with (Cs^2)_i = (Cs_i)^2
else if ( PS2 )
{
    forAll(U, k)
    {
        uqNut_k *= 0.0;
        forAll(U,i)
        {
            forAll(U,j)
            {
                tmpCijk = M3[k][j][i];
                if ( tmpCijk!=0 )
                {
                uqNut_k += tmpCijk
                        * turbulence[i]->getCs() * turbulence[i]->getCs()
                        * delta2*S[j];
                }
            }

        }
        uqNut_k *= 1.0/M1[k];
        turbulence[k]->uqCorrect(uqNut_k);
    }

}
