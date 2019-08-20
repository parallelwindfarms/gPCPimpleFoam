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

/******************* pMean, UMean, UMeanSigma, UMeanSigmaMean  *************************/

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
    pMean[k] = (1-beta)*pMean[k] + beta*p[k];
    UMean[k] = (1-beta)*UMean[k] + beta*U[k];
}

volScalarField      tmppMeanSigma(0.0*sqr(pMean[0]));
volVectorField      tmpUMeanSigma(0.0*cmptMultiply(UMean[0], UMean[0]));
forAll(U, k)
{
    if ( k > 0 )
    {
        tmppMeanSigma += M1[k] * sqr(pMean[k]);
        tmpUMeanSigma += M1[k] * cmptMultiply(UMean[k], UMean[k]);
    }
}

pMeanSigma = sqrt(tmppMeanSigma);
for (label dir=0; dir<3; dir++)
    UMeanSigma.replace(dir, sqrt(tmpUMeanSigma.component(dir)));

//UMeanSigmaMean = (1-beta)*UMeanSigmaMean + beta*UMeanSigma;

/********************************** RMean ***************************************/

if(time > dt)
{
    beta = dt/(time-dt);

    if ((time-restartAvgTime >= 0) && (restarted2 == true))
    {
        if (perfectRestart2)
            beta = dt/(time-dt-restartAvgTime+dt);
        else
            beta = dt/(time-dt-restartAvgTime);
    }
    if ((time-restartAvgTime >= 0) && (restarted2 == false))
    {
        beta = 1;
        restarted2 = true;
        if (time-restartAvgTime < SMALL)
            perfectRestart2 = true;
    }

    if (abs(beta) > 1) { beta = 1; }

    volSymmTensorField  tmpRMeank(0.0*RMean[0]);
    forAll(RMean, k)
    {
        tmpRMeank *= 0.0;

        forAll(RMean, i)
        {
            forAll(RMean, j)
            {
                tmpRMeank +=  M3[i][j][k] * symm((U[i]-UMean[i])*(U[j]-UMean[j]));
            }
        }

        RMean[k] = (1-beta)*RMean[k] + beta*tmpRMeank/M1[k];
    }
}

/*************** RMeanSigma, UMeanSigma2, UMeanSigmaPrim2Mean *******************/

volSymmTensorField  tmpRMeanSigma(0.0*cmptMultiply(RMean[0], RMean[0]));
forAll(U, k)
{
    if ( k > 0 )
        tmpRMeanSigma += M1[k] * cmptMultiply(RMean[k], RMean[k]);
}
for (label dir=0; dir<6; dir++)
    RMeanSigma.replace(dir, sqrt(tmpRMeanSigma.component(dir)));

/*
UMeanSigma2 = symm(UMeanSigma*UMeanSigma);

if (time > dt)
{
    UMeanSigmaPrime2Mean = (1-beta)*UMeanSigmaPrime2Mean
                         +  beta*symm((UMeanSigma-UMeanSigmaMean)
                                     *(UMeanSigma-UMeanSigmaMean));
}
*/


/*
volScalarField      tmppVar       (0.0*sqr(p[0]));
volScalarField      tmpnutVar     (0.0*sqr(turbulence[0]->nut()));
volVectorField      tmpUVar       (0.0*cmptMultiply(U[0], U[0]));

forAll(U, k)
{
    Ckk = M1[k];

    if ( k > 0 )
    {
        tmppVar    += Ckk * sqr(p[k]);
        tmpnutVar  += Ckk * sqr(turbulence[k]->nut());
        tmpUVar    += Ckk * cmptMultiply(U[k], U[k]);
    }
}

pSigma   = sqrt(tmppVar);
nutSigma = sqrt(tmpnutVar);

for (label dir=0; dir<3; dir++)
    USigma.replace(dir, sqrt(tmpUVar.component(dir)));
*/
