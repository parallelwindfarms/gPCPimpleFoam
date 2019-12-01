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

bool FS(0), PS1(0), PS2(0), truncated(0), FS_modelForm(0);
word mode (turbulenceProperties.subDict("LES").lookup("mode"));
word Truncated (turbulenceProperties.subDict("LES").lookup("truncated"));

// Full Spectral (FS)
if      ( mode == "FS"  ){ FS = 1; Info << "\n--Using Full Spectral (FS) Mode\n" << endl; }
// Pseudo Spectral (PS1)
else if ( mode == "PS1" ){ PS1 = 1; Info << "\n--Using Pseudo Spectral (PS1) Mode\n" << endl; }
// Pseudo Spectral (PS2)
else if ( mode == "PS2" ){ PS2 = 1; Info << "\n--Using Pseudo Spectral (PS2) Mode\n" << endl; }
// Full Spectral (FS) Model Form
else if ( mode == "FS_modelForm" ){ FS_modelForm = 1; Info << "\n--Using Full Spectral (FS) Model Form Mode\n" << endl; }
else
{
    Info<< "\n--> FOAM FATAL IO ERROR:\n\n"
        << "Invalid entry for spectral mode in constant/turbulenceProperties\n\n"
        << "FOAM exiting\n"
        << endl;
    return 0;
}


// Truncation of higher order terms
if ( Truncated == "on" || Truncated == "yes" || Truncated == "true" || Truncated == "1" )
{ truncated = 1; Info << "\n--Truncation of higher order terms is on\n" << endl; }
else if ( Truncated == "off" || Truncated == "no" || Truncated == "false" || Truncated == "0" )
{ truncated = 0; Info << "\n--Truncation of higher order terms is off\n" << endl; }
else
{
    Info<< "\n--> FOAM FATAL IO ERROR:\n\n"
        << "Invalid entry for truncated mode in constant/turbulenceProperties\n\n"
        << "FOAM exiting\n"
        << endl;
    return 0;
}
