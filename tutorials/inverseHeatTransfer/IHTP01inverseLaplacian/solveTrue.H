volScalarField& T = _T();
fvMesh& mesh = _mesh();
Foam::Time& runTime = _runTime();
simpleControl simple(mesh);
Info << "DT = " << DT << endl;
assignDirectBC();

#if defined(OFVER) && (OFVER == 6)
while (simple.loop(runTime))
#else
while (simple.loop())
#endif
{
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix TEqn
        (
            fvm::laplacian(DT, T)
        );
        TEqn.solve();
    }
}

Info << "Exporting" << endl;
ITHACAstream::exportSolution(T, "1", outputFolder, T.name());
//Reinitializing runTime
instantList Times = runTime.times();
runTime.setTime(Times[1], 1);
