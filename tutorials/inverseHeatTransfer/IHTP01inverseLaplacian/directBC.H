fvMesh& mesh = _mesh();
volScalarField& T(_T());
// Setting BC at the cold side
coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
label coldSideSize = T.boundaryField()[coldSide_ind].size();
Tf.resize(coldSideSize);
refGrad.resize(coldSideSize);
valueFraction.resize(coldSideSize);
forAll(Tf, faceI)
{
    scalar faceZ =
        mesh.boundaryMesh()[coldSide_ind].faceCentres()[faceI].z();
    scalar faceY =
        mesh.boundaryMesh()[coldSide_ind].faceCentres()[faceI].y();
    scalar faceX =
        mesh.boundaryMesh()[coldSide_ind].faceCentres()[faceI].x();
    Tf[faceI] = k / H * (b * faceX + c) + a * faceX * faceX +
                b * faceX * faceY + c * faceY - a * faceZ * faceZ + c;
}
// Setting BC at gammaEx1
gammaEx1_ind = mesh.boundaryMesh().findPatchID("gammaEx1");
label gammaEx1Size = T.boundaryField()[gammaEx1_ind].size();
heatFlux_gammaEx1.resize(gammaEx1Size);
forAll(heatFlux_gammaEx1, faceI)
{
    scalar faceZ =
        mesh.boundaryMesh()[gammaEx1_ind].faceCentres()[faceI].z();
    heatFlux_gammaEx1[faceI] = k * 2 * a * faceZ;
}
// Setting BC at gammaEx2
gammaEx2_ind = mesh.boundaryMesh().findPatchID("gammaEx2");
label gammaEx2Size = T.boundaryField()[gammaEx2_ind].size();
heatFlux_gammaEx2.resize(gammaEx2Size);
forAll(heatFlux_gammaEx2, faceI)
{
    scalar faceX =
        mesh.boundaryMesh()[gammaEx2_ind].faceCentres()[faceI].x();
    scalar faceY =
        mesh.boundaryMesh()[gammaEx2_ind].faceCentres()[faceI].y();
    heatFlux_gammaEx2[faceI] = - k * (2 * a * faceX + b * faceY);
}
// Setting BC at gammaEx3
gammaEx3_ind = mesh.boundaryMesh().findPatchID("gammaEx3");
label gammaEx3Size = T.boundaryField()[gammaEx3_ind].size();
heatFlux_gammaEx3.resize(gammaEx3Size);
forAll(heatFlux_gammaEx3, faceI)
{
    heatFlux_gammaEx3[faceI] = 0;
}
// Setting BC at gammaEx4
gammaEx4_ind = mesh.boundaryMesh().findPatchID("gammaEx4");
label gammaEx4Size = T.boundaryField()[gammaEx4_ind].size();
heatFlux_gammaEx4.resize(gammaEx4Size);
forAll(heatFlux_gammaEx4, faceI)
{
    scalar faceY =
        mesh.boundaryMesh()[gammaEx4_ind].faceCentres()[faceI].y();
    heatFlux_gammaEx4[faceI] = k * (b * faceY);
}
set_valueFraction();
forAll(mesh.boundaryMesh(), patchI)
{
    if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
    {
        ITHACAutilities::assignMixedBC(T, patchI, Tf, refGrad, valueFraction);
    }

    if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
    {
        ITHACAutilities::assignBC(T, patchI, - g / k);
    }

    if (patchI == mesh.boundaryMesh().findPatchID("gammaEx1"))
    {
        ITHACAutilities::assignBC(T, patchI, -heatFlux_gammaEx1 / k);
    }

    if (patchI == mesh.boundaryMesh().findPatchID("gammaEx2"))
    {
        ITHACAutilities::assignBC(T, patchI, -heatFlux_gammaEx2 / k);
    }

    if (patchI == mesh.boundaryMesh().findPatchID("gammaEx3"))
    {
        ITHACAutilities::assignBC(T, patchI, -heatFlux_gammaEx3 / k);
    }

    if (patchI == mesh.boundaryMesh().findPatchID("gammaEx4"))
    {
        ITHACAutilities::assignBC(T, patchI, -heatFlux_gammaEx4 / k);
    }
}

