if (pimple.firstIter() || moveMeshOuterCorrectors)
{
    mesh.update();

    if (mesh.changing())
    {
        MRF.update();

        if (correctPhi)
        {
            // Calculate absolute flux
            // from the mapped surface velocity
            phi[P-k] = mesh.Sf() & Uf[P-k]();

            #include "correctPhi.H"

            // Make the flux relative to the mesh motion
            fvc::makeRelative(phi[P-k], U[P-k]);
        }

        if (checkMeshCourantNo)
        {
            #include "uqMeshCourantNo.H"
        }
    }
}
