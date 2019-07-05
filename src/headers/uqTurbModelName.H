IOdictionary turbulenceProperties
(
    IOobject
    (
        "turbulenceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

/**************************** uqTurbulence model ******************************/

//bool uqSmagorinsky  ( turbulence[0]->getModelName() == "uqSmagorinsky" ? 1:0 );
//bool uqWALE         ( turbulence[0]->getModelName() == "uqWALE" ? 1:0 );

/****************************** Spectral Mode  *******************************/

bool FS(0), PS(0);
word mode (turbulenceProperties.subDict("LES").lookup("mode"));

// Full Spectral (FS)
if ( mode == "FS" )
    { FS = 1; Info << "--Using Full Spectral (FS) Mode" << endl; }

// Pseudo Spectral (PS)
if ( mode == "PS" )
    { PS = 1; Info << "--Using Pseudo Spectral (PS) Mode" << endl; }
