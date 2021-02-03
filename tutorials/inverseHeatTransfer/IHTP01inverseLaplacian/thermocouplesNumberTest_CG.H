Info << endl;
Info << "*********************************************************" << endl;
Info << "Performing test for the CG inverse solver" << endl;
Info << "Changint the position of the thermocouples on y" << endl;
Info << endl;
word outputFolder = "./ITHACAoutput/thermocouplesNumberTest_CG/";
volScalarField gTrueField = example_CG.list2Field(example_paramBC.gTrue);
ITHACAstream::exportSolution(gTrueField,
                             "1", outputFolder,
                             "gTrue");
example_CG.saveSolInLists = 1;

int testSize = 10;

Eigen::VectorXi numberTCperAxis = Eigen::VectorXi::LinSpaced(testSize, 1,
                                  testSize);
ITHACAstream::exportMatrix(numberTCperAxis, "numberTCperAxis", "eigen",
                           outputFolder);

for (int i = 0; i < testSize; i++)
{
    Info << "\n**************************************************** " << endl;
    Info << "\nTest " << i << endl;
    Info << "\n**************************************************** " << endl;
    Eigen::VectorXd TCcoord_temp = Eigen::VectorXd::LinSpaced(numberTCperAxis(
                                       i) + 2, 0.0, 1.0);
    Eigen::VectorXd TCcoord = TCcoord_temp.segment(1, numberTCperAxis(i));
    example_CG.thermocouplesPos.resize(numberTCperAxis(i) * numberTCperAxis(i));
    double TCyCoord = 0.2;
    int TCi = 0;

    for (int Xi = 0; Xi < numberTCperAxis(i); Xi++)
    {
        for (int Zi = 0; Zi < numberTCperAxis(i); Zi++)
        {
            example_CG.thermocouplesPos[TCi][0] = TCcoord(Xi);
            example_CG.thermocouplesPos[TCi][1] = TCyCoord;
            example_CG.thermocouplesPos[TCi][2] = TCcoord(Zi);
            Info << "thermocouplesPos[" << TCi << "] = " << example_CG.thermocouplesPos[TCi]
                 << endl;
            TCi++;
        }
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
