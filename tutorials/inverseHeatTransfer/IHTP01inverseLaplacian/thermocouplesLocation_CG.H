Info << endl;
Info << "*********************************************************" << endl;
Info << "Performing test for the CG inverse solver" << endl;
Info << "Changint the position of the thermocouples on y" << endl;
Info << endl;
word outputFolder = "./ITHACAoutput/thermocouplesLocationTest_CG/";
volScalarField gTrueField = example_CG.list2Field(example_paramBC.gTrue);
ITHACAstream::exportSolution(gTrueField,
                             "1", outputFolder,
                             "gTrue");
example_CG.saveSolInLists = 1;


int testSize = 20;
Eigen::VectorXd TCplane_Y = Eigen::VectorXd::LinSpaced(testSize, 0.1, 0.9);
ITHACAstream::exportMatrix(TCplane_Y, "TCplane_Y", "eigen",
                           outputFolder);

for (int i = 0; i < testSize; i++)
{
    Info << "\n**************************************************** " << endl;
    Info << "\nTest " << i << endl;
    Info << "\n**************************************************** " << endl;
    forAll(example_CG.thermocouplesPos, TCi)
    {
        example_CG.thermocouplesPos[TCi][1] = TCplane_Y(i);
    }

    /// Rereading thermocouples
    if ( remove( "./thermocouplesCellsID_mat.txt" ) != 0 )
    {
        perror( "\nError deleting thermocouples file\n" );
    }
    else
    {
        puts( "\nThermocouples file successfully deleted\n" );
    }

    example_CG.thermocouplesRead = 0;
    example_CG.readThermocouples();
    example_CG.Tmeas = example_CG.fieldValueAtThermocouples(T_true);

    if (example_CG.conjugateGradient())
    {
        Info << "\nCG converged" << endl;
        volScalarField heatFluxField = example_CG.list2Field(
                                           example_CG.gList[example_CG.gList.size() - 1]);
        ITHACAstream::exportSolution(heatFluxField,
                                     std::to_string(i + 1), outputFolder,
                                     "g");
        Info << "************************************" << endl;
        Info << endl << endl;
    }
    else
    {
        Info << "CG did not converged" << endl;
        Info << "************************************" << endl;
        Info << endl << endl;
    }
}

example_CG.postProcess(outputFolder, "g");

/// Rereading thermocouples
if ( remove( "./thermocouplesCellsID_mat.txt" ) != 0 )
{
    perror( "\nError deleting thermocouples file\n" );
}
else
{
    puts( "\nThermocouples file successfully deleted\n" );
}
