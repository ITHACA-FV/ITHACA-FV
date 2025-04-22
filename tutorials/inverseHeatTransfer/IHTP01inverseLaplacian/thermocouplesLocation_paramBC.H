Info << endl;
Info << "*********************************************************" << endl;
Info << "Performing test for the parameterized BC inverse solver" << endl;
Info << "Changing the position of the thermocouples on y" << endl;
Info << endl;
word outputFolder = "./ITHACAoutput/thermocouplesLocationTest_paramBC/";
volScalarField gTrueField = example_paramBC.list2Field(example_paramBC.gTrue);
ITHACAstream::exportSolution(gTrueField,
                             "1", outputFolder,
                             "gTrue");


int testSize = 20;
Eigen::VectorXd TCplane_Y = Eigen::VectorXd::LinSpaced(testSize, 0.1, 0.9);
ITHACAstream::exportMatrix(TCplane_Y, "TCplane_Y", "eigen",
                           outputFolder);

for (int i = 0; i < testSize; i++)
{
    Info << "\n**************************************************** " << endl;
    Info << "\nTest " << i << endl;
    Info << "\n**************************************************** " << endl;
    forAll(example_paramBC.thermocouplesPos, TCi)
    {
        example_paramBC.thermocouplesPos[TCi][1] = TCplane_Y(i);
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
}

example_paramBC.postProcess(outputFolder, "g");

/// Rereading thermocouples
if ( remove( "./thermocouplesCellsID_mat.txt" ) != 0 )
{
    perror( "\nError deleting thermocouples file\n" );
}
else
{
    puts( "\nThermocouples file successfully deleted\n" );
}
