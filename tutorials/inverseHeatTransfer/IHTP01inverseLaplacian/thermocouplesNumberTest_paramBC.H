Info << endl;
Info << "*********************************************************" << endl;
Info << "Performing test for the paramBC inverse solver" << endl;
Info << "Changing the position of the thermocouples on y" << endl;
Info << endl;
word outputFolder = "./ITHACAoutput/thermocouplesNumberTest_paramBC/";
word outputFolder_bestFit =
    "./ITHACAoutput/thermocouplesNumberTest_paramBC_bestFit/";
word outputFolder_bestInt =
    "./ITHACAoutput/thermocouplesNumberTest_paramBC_bestInt/";
volScalarField gTrueField = example_paramBC.list2Field(example_paramBC.gTrue);
ITHACAstream::exportSolution(gTrueField,
                             "1", outputFolder,
                             "gTrue");

int testSize = 40;

Eigen::VectorXi numberTCperAxis = Eigen::VectorXi::LinSpaced(testSize, 1,
                                  testSize + 1);
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
    example_paramBC.thermocouplesPos.resize(numberTCperAxis(i) * numberTCperAxis(
            i));
    double TCyCoord = 0.2;
    int TCi = 0;

    for (int Zi = 0; Zi < numberTCperAxis(i); Zi++)
    {
        for (int Xi = 0; Xi < numberTCperAxis(i); Xi++)
        {
            example_paramBC.thermocouplesPos[TCi][0] = TCcoord(Xi);
            example_paramBC.thermocouplesPos[TCi][1] = TCyCoord;
            example_paramBC.thermocouplesPos[TCi][2] = TCcoord(Zi);
            Info << "thermocouplesPos[" << TCi << "] = " <<
                 example_paramBC.thermocouplesPos[TCi] << endl;
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

    volScalarField Ttrue(example_paramBC._T());

    for (label i = 0; i < Ttrue.internalField().size(); i++)
    {
        auto cx = Ttrue.mesh().C()[i].component(vector::X);
        auto cy = Ttrue.mesh().C()[i].component(vector::Y);
        auto cz = Ttrue.mesh().C()[i].component(vector::Z);
        Ttrue.ref()[i] = a * cx * cx + b * cx * cy + c * cy - a * cz * cz + c;
    }

    example_paramBC.thermocouplesRead = 0;
    example_paramBC.readThermocouples();
    example_paramBC.Tmeas = Eigen::VectorXd::Zero(numberTCperAxis(
                                i) * numberTCperAxis(i));
    example_paramBC.Tmeas = example_paramBC.fieldValueAtThermocouples(Ttrue);
    example_paramBC.set_gParametrized("rbf", rbfShapePar);
    example_paramBC.parameterizedBCoffline(1);
    example_paramBC.parameterizedBC("fullPivLU");
    volScalarField gParametrizedField = example_paramBC.list2Field(
                                            example_paramBC.g);
    ITHACAstream::exportSolution(gParametrizedField,
                                 std::to_string(i + 1),
                                 outputFolder,
                                 "g");
    volScalarField& T(example_paramBC._T());
    ITHACAstream::exportSolution(T,
                                 std::to_string(i + 1),
                                 outputFolder,
                                 "T");
    // Best fit
    example_paramBC.parameterizedBCpostProcess(example_paramBC.bestApproximator());
    gParametrizedField = example_paramBC.list2Field(
                             example_paramBC.g);
    ITHACAstream::exportSolution(gParametrizedField,
                                 std::to_string(i + 1),
                                 outputFolder_bestFit,
                                 "g");
    volScalarField& T2(example_paramBC._T());
    ITHACAstream::exportSolution(T2,
                                 std::to_string(i + 1),
                                 outputFolder_bestFit,
                                 "T");
    // Interpolation
    example_paramBC.parameterizedBCpostProcess(example_paramBC.bestInterpolation());
    gParametrizedField = example_paramBC.list2Field(
                             example_paramBC.g);
    ITHACAstream::exportSolution(gParametrizedField,
                                 std::to_string(i + 1),
                                 outputFolder_bestInt,
                                 "g");
    volScalarField& T3(example_paramBC._T());
    ITHACAstream::exportSolution(T3,
                                 std::to_string(i + 1),
                                 outputFolder_bestInt,
                                 "T");
}

example_paramBC.postProcess(outputFolder, "g");
example_paramBC.postProcess(outputFolder_bestFit, "g");
example_paramBC.postProcess(outputFolder_bestInt, "g");

/// Rereading thermocouples
if ( remove( "./thermocouplesCellsID_mat.txt" ) != 0 )
{
    perror( "\nError deleting thermocouples file\n" );
}
else
{
    puts( "\nThermocouples file successfully deleted\n" );
}

