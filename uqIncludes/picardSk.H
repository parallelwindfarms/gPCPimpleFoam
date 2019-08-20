// Setting number of Picard Iteration
label picItr(0);
label picItrMax(runTime.controlDict().lookupOrDefault<label>("picItrMax",10));

// Initialsing Sk using Uk
forAll(S, k)
{
    L[k] = 2*(symm(fvc::grad(U[k])) && symm(fvc::grad(U[k])));
    S[k] = sqrt(L[k]);
    Info<< "mean " << S[k].weightedAverage(mesh.V()) << "\t" << max(S[k]) << endl;
}


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
    dimensionedScalar deltaPicSk("deltaPicSk", pow(dimTime,-1), 1.0);
    dimensionedScalar deltaPicSkRel("deltaPicSk", dimless, 1.0);


    while( (picItr <= picItrMax) && (deltaPicSkRel > deltaPicTol) )
    {
        deltaPicSk = S[0].weightedAverage(mesh.V());
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
            S[k] = (L[k]*Ckk - numerator)/(S[k]*M3[k][k][k] + denominator + SMALL*unitDenom);
        }

        deltaPicSkRel = mag(S[0].weightedAverage(mesh.V()) - deltaPicSk) / mag(deltaPicSk);
        if( (picItr == picItrMax) || (deltaPicSkRel < deltaPicTol))
        {
            forAll(S, k)
            {
                S[k] = (S[k] + mag(S[k]))/2.0;
                Info<< "mean " << S[k].weightedAverage(mesh.V()) << "\t" << max(S[k]) << endl;
            }
            Info<< "\nMax Picard iterations used : " << picItr << endl;
        }
        picItr++;
    }
}

Info << endl;
