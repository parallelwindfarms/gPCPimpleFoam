// Setting number of Picard Iteration
label  picItr(0);
label  picItrMax(runTime.controlDict().lookupOrDefault<label>("picItrMax",100));
scalar relax_Sk(runTime.controlDict().lookupOrDefault<scalar>("relax_Sk",1.0));
scalar scale_Sk(runTime.controlDict().lookupOrDefault<scalar>("scale_Sk",1.0));

// Initialsing Lk and Sk using Uk
//if (runTime.value() == runTime.deltaT().value())
//{
        forAll(S, k)
        {
            L[k] = 2*(symm(fvc::grad(U[k])) && symm(fvc::grad(U[k])));
            S[k] = sqrt(L[k]);
            if ( k>0 ) S[k] = S[k] * pow(scale_Sk,k);
            Info<< "|S|"    << k
                << "\tavg " << S[k].internalField().weightedAverage(mesh.V()).value()
                << "\tmax " << max(S[k].internalField()).value()
                << "\tmin " << min(S[k].internalField()).value()
                << endl;
        }
//}


// Need some fixing or a new method
if(picItrMax != 0)
{
    // Estimating Lk using Ui and Uj
    forAll(L, k)
    {
        Ckk  = M1[k];
        L[k] *= 0.0;

        forAll(L, i)
        {
            forAll(L, j)
            {
                Cijk = M3[i][j][k];
                L[k] += (symm(fvc::grad(U[i])) && symm(fvc::grad(U[j]))) * Cijk;
            }
        }
        L[k] = 2*L[k]/Ckk;
        //L[k] = mag(L[k]);
        //Info<< max(L[k]) << endl;
    }

    // Using picard iterations to estimate Sk
    volScalarField unitDenom
    (
        IOobject
        (
            "unitDenom",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("unitDenom",pow(dimTime,-1), 1.0)
    );
    volScalarField numerator(sqr(unitDenom)*0.0);
    volScalarField denominator(unitDenom*0.0);

    picItr = 1;
    dimensionedScalar deltaPicTol("deltaPicTol", dimless, 1e-03);
    dimensionedScalar Sk_prev("Sk_prev", pow(dimTime,-1), 0.0);
    dimensionedScalar deltaPicSkRel("deltaPicSk", dimless, 1.0);
    scalar r = relax_Sk;

    while( (picItr <= picItrMax) && (deltaPicSkRel > deltaPicTol) )
    {
        Sk_prev = S[0].weightedAverage(mesh.V());
        forAll(S, k)
        {
            Ckk  = M1[k];
            numerator   *= 0.0;
            denominator *= 0.0;

            forAll(S, i)
            {
                if(i != k)
                {
                    forAll(S, j)
                    {
                        if(j != k)
                        {
                            Cijk = M3[i][j][k];
                            if (Cijk != 0.0)
                                numerator += S[i]*S[j]*Cijk;
                        }
                    }
                    denominator += 2.0*S[i]*M3[i][k][k];
                }
            }
            S[k] = (1-r)*S[k] + r*(L[k]*Ckk - numerator)/(S[k]*M3[k][k][k] + denominator + SMALL*unitDenom);
            //S[k] = (-denominator+sqrt(sqr(denominator)+4*M3[k][k][k]*(L[k]*Ckk-numerator)))/(2*M3[k][k][k]);
            //S[k] = sqrt((L[k]*Ckk - numerator - denominator*S[k])/(M3[k][k][k]+SMALL));
            //S[k] = (S[k] + mag(S[k]))/2.0;
        }

        deltaPicSkRel = mag(S[0].weightedAverage(mesh.V()) - Sk_prev) / mag(Sk_prev);
        if( (picItr == picItrMax) || (deltaPicSkRel < deltaPicTol))
        {
            Info<< "\nMax Picard iterations used : " << picItr << endl;

            forAll(S, k)
            {
                Info<< "mean " << S[k].weightedAverage(mesh.V())
                    << "\t" << max(S[k])
                    << "\t" << min(S[k])
                    << endl;
            }
        }
        picItr++;
    }
}

Info << endl;


