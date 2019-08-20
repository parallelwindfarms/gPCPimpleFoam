// UQ pre-processing
Info<< "--Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "--Reading polynomial order (n)\n" << endl;
label order (readLabel(transportProperties.lookup("order")));

Info<< "--Reading parametric dimension (dim)\n" << endl;
label dim (readLabel(transportProperties.lookup("dim")));

/*
Info<< "--Reading random variable(s) distribution type(s)\n" << endl;
word distType (transportProperties.lookup("distType"));
*/

Info<< "--Calculating P" << endl;
label P = factorial(order + dim) / (factorial(order)*factorial(dim)) - 1;
Info<< "\n  P = " << P << " => k = [0," << P << "]\n" << endl;

Info<< "--Reading Cijk, Cijkl (UQ coeffs) calculated by chaospy\n" << endl;
#include "createUQfields.H"

/*
Info<< "--Reading viscosity nodes nu[k]\n" << endl;

dimensionedScalar zeroViscDimSclr
("zeroViscDimSclr",dimensionSet(0,2,-1,0,0,0,0),0.);
PtrList<dimensionedScalar> nu(P+1);
forAll(nu, k)
{
    Info<< "Reading nu[" << k << "]" << endl;
    nu.set
    (
        k,
        new dimensionedScalar
        (
            transportProperties.lookupOrDefault<dimensionedScalar>
            (("nu"+std::to_string(k)),zeroViscDimSclr)
        )
    );
}
*/
/*********************************************************************/

// Reading fields
PtrList<volScalarField> p(P+1);

Info<< "--Reading field p[k]\n" << endl;
forAll(p, k)
{
    Info << "Reading field p" << k << endl;
    p.set(k,new volScalarField
    (
        IOobject
        (
            "p" + std::to_string(k),
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ));
}

PtrList<volVectorField> U(P+1);

Info<< "\n--Reading field U[k]\n" << endl;
forAll(U, k)
{
    Info << "Reading field U" << k << endl;
    U.set(k,new volVectorField
    (
        IOobject
        (
            "U" + std::to_string(k),
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ));
}

PtrList<surfaceScalarField> phi(P+1);

Info<< "\n--Reading field phi[k]\n" << endl;
forAll(phi, k)
{
    Info << "Reading field phi" << k << endl;
    phi.set(k,new surfaceScalarField
    (
        IOobject
        (
            "phi" + std::to_string(k),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(U[k])
    ));
}

label pRefCell = 0;
scalar pRefValue = 0.0;

forAll(p, k)
{
    setRefCell(p[k], pimple.dict(), pRefCell, pRefValue);
    mesh.setFluxRequired(p[k].name());
}


#include "createMRF.H"
#include "createFvOptions.H"

/*********************************************************************/

// Transport and turbulence model(s) initialization
Info<< "\n--Reading transport and turbulence models\n" << endl;

uqsinglePhaseTransportModel* laminarTransport[P+1];
List<autoPtr<uqincompressible::uqturbulenceModel>> turbulence(P+1);

forAll(turbulence, k)
{
    Info<<"\n-Models for gPC node : " << k << endl;
    laminarTransport[k] = new uqsinglePhaseTransportModel
        (U[k], phi[k], "nu"+std::to_string(k));

    turbulence[k] = uqincompressible::uqturbulenceModel::New
        (U[k], phi[k], *laminarTransport[k], "turbulenceProperties", k);
}
#include "uqTurbModelName.H"
volScalarField uqNu_0(turbulence[0]->nu());

/*
forAll(turbulence, k)
  {
    //Info<< "nu" << k << "\n" << laminarTransport[k]->nu() << endl;
    //Info<< "nut" << k << "\n" << turbulence[k]->nut() << endl;
    //Info<< "nuEff" << k << "\n" << turbulence[k]->nuEff() << endl;
    //Info<< turbulence[k]->getCs() << endl;
  }
Info<< endl;
//*/

// Initializing UQ solution modes with 10% sd based on deterministic solution
/*
scalar sd(0.1);
forAll(p, k)
{
    if ( k>0 )
    {
        p[k]    = pow(sd, k+1) * p[0];
        U[k]    = pow(sd, k+1) * U[0];
        phi[k]  = pow(sd, k+1) * phi[0];
    }
}
*/
/*********************************************************************/

// UQ mean, st.dev. and their ModEff initialization (for postProcessing)

// Sigma -> Mean -> Prime2Mean
/*
Info<< "\n--Creating pSigma field pSigma\n" << endl;

volScalarField pSigma
(
    IOobject
    (
        "pSigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p[1]*0.0
);

Info<< "\n--Creating nutSigma field nutSigma\n" << endl;

volScalarField nutSigma
(
    IOobject
    (
        "nutSigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    turbulence[1]->nut()*0.0
);


Info<< "\n--Creating USigma field USigma\n" << endl;
volVectorField USigma
(
    IOobject
    (
        "USigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U[1]*0.0
);
*/
// Mean -> Sigma, Prime2Mean -> Sigma
PtrList<volScalarField> pMean(P+1);
Info<< "\n--Reading field p[k]Mean\n" << endl;
forAll(pMean, k)
{
    pMean.set(k,new volScalarField
    (
        IOobject
        (
            "pMean" + std::to_string(k),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        p[k]*0.0
    ));
}

Info<< "\n--Creating field pMeanSigma\n" << endl;
volScalarField pMeanSigma
(
    IOobject
    (
        "pMeanSigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    p[1]*0.0
);

PtrList<volVectorField> UMean(P+1);
Info<< "\n--Reading field U[k]Mean\n" << endl;
forAll(UMean, k)
{
    UMean.set(k,new volVectorField
    (
        IOobject
        (
            "UMean" + std::to_string(k),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U[k]*0.0
    ));
}

Info<< "\n--Creating field UMeanSigma\n" << endl;
volVectorField UMeanSigma
(
    IOobject
    (
        "UMeanSigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U[1]*0.0
);

PtrList<volSymmTensorField> RMean(P+1);
Info<< "\n--Creating field RMean[k]\n" << endl;
forAll(RMean, k)
{
    RMean.set(k,new volSymmTensorField
    (
        IOobject
        (
            "RMean" + std::to_string(k),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        symm(U[0]*U[0])*0.0
    ));
}

Info<< "\n--Creating field RMeanSigma\n" << endl;
volSymmTensorField RMeanSigma
(
    IOobject
    (
        "RMeanSigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    symm(U[0]*U[0])*0.0
);


/*
Info<< "\n--Creating UMeanSigmaMean field UMeanSigmaMean\n" << endl;
volVectorField UMeanSigmaMean
(
    IOobject
    (
        "UMeanSigmaMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U[1]*0.0
);

Info<< "\n--Creating UMeanSigma2 field UMeanSigma2\n" << endl;
volSymmTensorField UMeanSigma2
(
    IOobject
    (
        "UMeanSigma2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    symm(UMeanSigma*UMeanSigma)
);

Info<< "\n--Creating UMeanSigmaPrime2Mean field UMeanSigmaPrime2Mean\n" << endl;
volSymmTensorField UMeanSigmaPrime2Mean
(
    IOobject
    (
        "UMeanSigmaPrime2Mean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    symm(UMeanSigma*UMeanSigma)
);
*/


/***************************************************************/

// On the fly calculations of non-polynomial terms

PtrList<volScalarField> L(P+1);
Info<< "\n--Creating field L[k]\n" << endl;
forAll(L, k)
{
    L.set(k,new volScalarField
    (
        IOobject
        (
            "L" + std::to_string(k),
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        2*(symm(fvc::grad(U[k])) && symm(fvc::grad(U[k])))
    ));
}

PtrList<volScalarField> S(P+1);
Info<< "\n--Creating field S[k]\n" << endl;
forAll(S, k)
{
    S.set(k,new volScalarField
    (
        IOobject
        (
            "S" + std::to_string(k),
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sqrt(L[k])
    ));
    Info<< max(S[k]) << endl;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Restart booleans for UQ

bool restarted1(false);
bool restarted2(false);
bool perfectRestart1(false);
bool perfectRestart2(false);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //





































/*
Info<< "\n--Creating ModEffuqMean field ModEff\n" << endl;
volSymmTensorField ModEffuqMean
(
    IOobject
    (
        "ModEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    turbulence[0]->R()*0.0
);

Info<< "\n--Creating ModEffuqSigma field ModEffSigma\n" << endl;
volSymmTensorField ModEffuqSigma
(
    IOobject
    (
        "ModEffSigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    turbulence[1]->R()*0.0
);
*/
