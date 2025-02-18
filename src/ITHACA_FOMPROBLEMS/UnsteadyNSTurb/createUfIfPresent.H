
if (mesh.dynamic())
{
    Info << "Constructing face velocity Uf\n" << endl;
    _Uf.reset
    (
        new surfaceVectorField
        (
            IOobject
            (
                "Uf",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(U)
        )
    );
}
