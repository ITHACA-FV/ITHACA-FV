Info << endl;
Info << "*********************************************************" << endl;
Info << "Testing CG inverse solver with NOISY data" << endl;
Info << "Performing " << Ntests << " tests." << endl;
Info << endl;
word outputFolder = "./ITHACAoutput/CGnoiseTest/";
volScalarField gTrueField = example_CG.list2Field(example_CG.gTrue);
ITHACAstream::exportSolution(gTrueField,
                             "1", outputFolder,
                             "gTrue");
example_CG.saveSolInLists = 1;
Eigen::VectorXd TmeasOrig = example_CG.Tmeas;

for (label i = 0; i < Ntests; i++)
{
    Info << "Test " << i << endl;
    example_CG.addNoise(noiseLevel);
    Info << "Stopping for J < " << example_CG.Jtol << endl;

    if (example_CG.conjugateGradient())
    {
        Info << "CG converged" << endl;
        volScalarField heatFluxField = example_CG.list2Field(
                                           example_CG.gList[example_CG.gList.size() - 1]);
        ITHACAstream::exportSolution(heatFluxField,
                                     std::to_string(i + 1), outputFolder,
                                     "g_CG");
        Info << "************************************" << endl;
        Info << endl << endl;
    }
    else
    {
        Info << "CG did not converged" << endl;
        Info << "************************************" << endl;
        Info << endl << endl;
    }

    example_CG.Tmeas = TmeasOrig;
}

example_CG.postProcess(outputFolder, "g_CG");
Info << "*********************************************************" << endl;
Info << endl;
