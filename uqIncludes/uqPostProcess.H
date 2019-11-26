/****************************** Restart settings  **************************************/

scalar dt   = runTime.deltaT().value();
scalar time = runTime.value();
scalar beta = dt/time;

scalar restartAvgTime
(
    runTime.controlDict().lookupOrDefault<scalar>
        ("restartAvgTime",runTime.endTime().value()+1)
);

// No restart for averaging, i.e. start averaging from t=0
if (restartAvgTime == 0.0)
{
    restartAvgTime = runTime.endTime().value()+1;
}

// If case is restarted at a time after restartAvgTime
if (time > restartAvgTime)
{
   restarted1 = true;
   restarted2 = true;
}

/********************** nut, p, U - Mean, Sigma, MeanSigma **********************/

if ((time-restartAvgTime >= 0) && (restarted1 == true))
{
    if (perfectRestart1)
        beta = dt/(time-restartAvgTime+dt);
    else
        beta = dt/(time-restartAvgTime);
}
if ((time-restartAvgTime >= 0) && (restarted1 == false))
{
    beta = 1;
    restarted1 = true;
    if (time-restartAvgTime < SMALL)
        perfectRestart1 = true;
}

forAll(UMean, k)
{
    nutMean[k]  = (1-beta)*nutMean[k] + beta*turbulence[k]->nut();
    pMean[k]    = (1-beta)*pMean[k]   + beta*p[k];
    UMean[k]    = (1-beta)*UMean[k]   + beta*U[k];
}

volScalarField      tmpnutVar     (0.0*sqr(turbulence[0]->nut()));
volScalarField      tmppVar       (0.0*sqr(p[0]));
volVectorField      tmpUVar       (0.0*cmptMultiply(U[0], U[0]));
volScalarField      tmpnutMeanVar (0.0*sqr(nutMean[0]));
volScalarField      tmppMeanVar   (0.0*sqr(pMean[0]));
volVectorField      tmpUMeanVar   (0.0*cmptMultiply(UMean[0], UMean[0]));

forAll(U, k)
{
    if ( k > 0 )
    {
        tmpnutVar       += M1[k] * sqr(turbulence[k]->nut());
        tmppVar         += M1[k] * sqr(p[k]);
        tmpUVar         += M1[k] * cmptMultiply(U[k], U[k]);
        tmpnutMeanVar   += M1[k] * sqr(nutMean[k]);
        tmppMeanVar     += M1[k] * sqr(pMean[k]);
        tmpUMeanVar     += M1[k] * cmptMultiply(UMean[k], UMean[k]);
    }
}

nutSigma     = sqrt(tmpnutVar);
nutMeanSigma = sqrt(tmpnutMeanVar);

pSigma       = sqrt(tmppVar);
pMeanSigma   = sqrt(tmppMeanVar);

for (label dir=0; dir<3; dir++)
{
    USigma.replace(dir, sqrt(tmpUVar.component(dir)));
    UMeanSigma.replace(dir, sqrt(tmpUMeanVar.component(dir)));
}

/*************************** R - Mean, Sigma, MeanSigma *************************/

if(time > dt)
{

    beta = dt/(time-dt);

    if ((time-restartAvgTime >= 0) && (restarted2 == true))
    {
        if (perfectRestart2)
            beta = dt/(time-restartAvgTime + dt*1e-3);
        else
            beta = dt/(time-dt-restartAvgTime + dt*1e-3);
    }
    if ((time-restartAvgTime >= 0) && (restarted2 == false))
    {
        beta = 1;
        restarted2 = true;
        if (time-restartAvgTime < SMALL)
            perfectRestart2 = true;
    }

    if (abs(beta) >= 1) { beta = 1; }

    volSymmTensorField  tmpRk   (0.0*RMean[0]);
    volSymmTensorField  tmpRVar (0.0*cmptMultiply(RMean[0], RMean[0]));

    forAll(RMean, k)
    {
        tmpRk *= 0.0;

        forAll(RMean, i)
        {
            forAll(RMean, j)
            {
                tmpRk +=  M3[i][j][k] * symm((U[i]-UMean[i])*(U[j]-UMean[j]));
            }
        }

        RMean[k] = (1-beta)*RMean[k] + beta*tmpRk/M1[k];
        if ( k > 0)
            tmpRVar += M1[k] * cmptMultiply(tmpRk, tmpRk);
    }

    for (label dir=0; dir<6; dir++)
        RSigma.replace(dir, sqrt(tmpRVar.component(dir)));
}

volSymmTensorField  tmpRMeanVar (0.0*cmptMultiply(RMean[0], RMean[0]));
forAll(U, k)
{
    if ( k > 0 )
        tmpRMeanVar += M1[k] * cmptMultiply(RMean[k], RMean[k]);
}

for (label dir=0; dir<6; dir++)
    RMeanSigma.replace(dir, sqrt(tmpRMeanVar.component(dir)));
