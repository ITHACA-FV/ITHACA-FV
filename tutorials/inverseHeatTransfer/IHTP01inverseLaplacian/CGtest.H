Info << endl;
Info << "*********************************************************" << endl;
Info << "Performing test for the CG inverse solver" << endl;
Info << endl;
word outputFolder = "./ITHACAoutput/CGtest/";
volScalarField gTrueField = example_CG.list2Field(example_CG.gTrue);
ITHACAstream::exportSolution(gTrueField,
                             "1", outputFolder,
                             "gTrue");
example_CG.saveSolInLists = 1;
auto t1 = std::chrono::high_resolution_clock::now();

if (example_CG.conjugateGradient())
{
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>
                    ( t2 - t1 ).count() / 1e6;
    std::cout << "Duration = " << duration << " seconds" << std::endl;
    Info << "CG converged" << endl;
    PtrList<volScalarField> heatFluxField;
    forAll(example_CG.gList, solutionI)
    {
        heatFluxField.append(example_CG.list2Field(
                                 example_CG.gList[solutionI]).clone());
        ITHACAstream::exportSolution(heatFluxField[solutionI],
                                     std::to_string(solutionI + 1), outputFolder,
                                     "g_CG");
    }
    example_CG.postProcess(outputFolder, "g_CG");
}
else
{
    Info << "CG did not converged" << endl;
}

Info << "*********************************************************" << endl;
Info << endl;
