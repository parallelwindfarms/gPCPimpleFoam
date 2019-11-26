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
            phi[k] = mesh.Sf() & Uf[k]();

            #include "correctPhi.H"

            // Make the flux relative to the mesh motion
            fvc::makeRelative(phi[k], U[k]);
        }

        if (checkMeshCourantNo)
        {
            #include "uqMeshCourantNo.H"
        }
    }
}
