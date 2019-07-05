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






































/*

//ModEffuqMean    = - 2*(turbulence[0]->nuEff())*symm(fvc::grad(UuqMean));


// -- UQ Std.dev. Calculation

volScalarField      tmppVar       (0.0*sqr(p[0]));
volScalarField      tmpnutVar     (0.0*sqr(turbulence[0]->nut()));
volVectorField      tmpUVar       (0.0*cmptMultiply(U[0], U[0]));
//volSymmTensorField  tmpModEffVar    (0.0*cmptMultiply(ModEffuqSigma, ModEffuqSigma));
//volSymmTensorField  tmpModEff       (ModEffuqSigma);

forAll(U, k)
{
    Ckk = M1[k];

    if ( k > 0 )
    {
        tmppVar    += Ckk * sqr(p[k]);
        tmpnutVar  += Ckk * sqr(turbulence[k]->nut());
        tmpUVar    += Ckk * cmptMultiply(U[k], U[k]);
        //tmpModEff     = - 2*(turbulence[k]->nuEff())*symm(fvc::grad(U[k]));
        //tmpModEffVar += Ckk * cmptMultiply(tmpModEff, tmpModEff);
    }
}

pSigma   = sqrt(tmppVar);
nutSigma = sqrt(tmpnutVar);

USigma.replace(0, sqrt(tmpUVar.component(0)));
USigma.replace(1, sqrt(tmpUVar.component(1)));
USigma.replace(2, sqrt(tmpUVar.component(2)));

//ModEffuqSigma.replace(0, sqrt(tmpModEffVar.component(0)));
//ModEffuqSigma.replace(1, sqrt(tmpModEffVar.component(1)));
//ModEffuqSigma.replace(2, sqrt(tmpModEffVar.component(2)));

*/
