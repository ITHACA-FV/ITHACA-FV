while (runTime.loop())
{
    scalar time = runTime.value();
    //Info<< "Time = " << runTime.timeName() << nl << endl;
    updateBC();

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix TEqn
        (
            fvm::ddt(T) - fvm::laplacian(diffusivity, T)
        );

        fvOptions.constrain(TEqn);
        TEqn.solve();
        fvOptions.correct(T);
    }
    ITHACAstream::exportSolution(T, runTime.timeName(),
                     folderOffline,
                     "T");


    runTime.printExecutionTime(Info);
    runTime.write();
}
