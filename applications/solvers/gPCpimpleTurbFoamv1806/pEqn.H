volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U[k], p[k]));
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
  + MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U[k], phi[k], Uf[k]))
);

MRF.makeRelative(phiHbyA);

if (p[k].needReference())
{
    fvc::makeRelative(phiHbyA, U[k]);
    adjustPhi(phiHbyA, U[k], p[k]);
    fvc::makeAbsolute(phiHbyA, U[k]);
}

tmp<volScalarField> rAtU(rAU);

if (pimple.consistent())
{
    rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p[k])*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p[k]);
}

if (pimple.nCorrPISO() <= 1)
{
    //tUEqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p[k], U[k], phiHbyA, rAtU(), MRF);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p[k]) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p[k].select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi[k] = phiHbyA - pEqn.flux();
    }
}

#include "uqContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
p[k].relax();

U[k] = HbyA - rAtU*fvc::grad(p[k]);
U[k].correctBoundaryConditions();
fvOptions.correct(U[k]);

// Correct Uf if the mesh is moving
fvc::correctUf(Uf[k], U[k], phi[k]);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi[k], U[k]);
