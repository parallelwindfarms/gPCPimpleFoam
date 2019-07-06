
IOdictionary UQPdt
(
    IOobject
    (
        "gPCcoeffs",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

label Pplus1 = P+1;
scalar Ckk(0);
scalar M1[Pplus1];

for (label k=0; k< Pplus1; k++)
  M1[k] = UQPdt.lookupOrDefault<scalar>
               ("M["+std::to_string(k)+"]",0.);


scalar Cijk(0);
scalar M3[Pplus1][Pplus1][Pplus1];

for (label i=0; i < Pplus1; i++)
  for (label j=0; j< Pplus1; j++)
    for (label k=0; k< Pplus1; k++)
      M3[i][j][k] = UQPdt.lookupOrDefault<scalar>
                   ("M["+std::to_string(i)+
                    "]["+std::to_string(j)+
                    "]["+std::to_string(k)+"]",0.);


scalar Cijkl(0);
scalar M4[Pplus1][Pplus1][Pplus1][Pplus1];

for (label i=0; i < Pplus1; i++)
  for (label j=0; j< Pplus1; j++)
    for (label k=0; k< Pplus1; k++)
      for (label l=0; l< Pplus1; l++)
        M4[i][j][k][l] = UQPdt.lookupOrDefault<scalar>
                        ("M["+std::to_string(i)+
                         "]["+std::to_string(j)+
                         "]["+std::to_string(k)+
                         "]["+std::to_string(l)+"]",0.);


/*
for (label i=0; i < Pplus1; i++)
    for (label j=0; j<Pplus1; j++)
        for (label k=0; k<Pplus1; k++)
            if (M3[i][j][k] != 0)
            Info<< "M3["<<i<<"]["<<j<<"]["<<k<<"] = "
                << M3[i][j][k] << endl;

for (label i=0; i < Pplus1; i++)
    for (label j=0; j< Pplus1; j++)
        for (label k=0; k< Pplus1; k++)
            for (label l=0; l< Pplus1; l++)
                if (M4[i][j][k][l] != 0)
                Info<< "M4["<<i<<"]["<<j<<"]["<<k<<"]["<<l<<"] = "
                    <<  M4[i][j][k][l] << endl;

//*/