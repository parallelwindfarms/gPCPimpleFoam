volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U[P-k], p[P-k]));
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
  + MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U[P-k], phi[P-k], Uf[P-k]))
);

MRF.makeRelative(phiHbyA);

if (p[P-k].needReference())
{
    fvc::makeRelative(phiHbyA, U[P-k]);
    adjustPhi(phiHbyA, U[P-k], p[P-k]);
    fvc::makeAbsolute(phiHbyA, U[P-k]);
}

tmp<volScalarField> rAtU(rAU);

if (pimple.consistent())
{
    rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p[P-k])*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p[P-k]);
}

if (pimple.nCorrPISO() <= 1)
{
    //tUEqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p[P-k], U[P-k], phiHbyA, rAtU(), MRF);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p[P-k]) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p[P-k].select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi[P-k] = phiHbyA - pEqn.flux();
    }
}

#include "uqContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
p[P-k].relax();

U[P-k] = HbyA - rAtU*fvc::grad(p[P-k]);
U[P-k].correctBoundaryConditions();
fvOptions.correct(U[P-k]);

// Correct Uf if the mesh is moving
fvc::correctUf(Uf[P-k], U[P-k], phi[P-k]);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi[P-k], U[P-k]);
