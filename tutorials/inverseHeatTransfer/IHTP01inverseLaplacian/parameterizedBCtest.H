Info << endl;
Info << "*********************************************************" << endl;
Info << "Performing test for the parameterized BC inverse solver" << endl;
Info << endl;
word outputFolder = "./ITHACAoutput/parameterizedBCtest/";
volScalarField gTrueField = example_paramBC.list2Field(example_paramBC.gTrue);
ITHACAstream::exportSolution(gTrueField,
                             "1", outputFolder,
                             "gTrue");
List<word> linSys_solvers;
linSys_solvers.resize(1);
linSys_solvers[0] = "fullPivLU";
Eigen::VectorXd residualNorms;
residualNorms.resize(linSys_solvers.size());
auto tO1 = std::chrono::high_resolution_clock::now();
example_paramBC.set_gParametrized("rbf", rbfShapePar);
example_paramBC.parameterizedBCoffline();
auto tO2 = std::chrono::high_resolution_clock::now();
auto duration = std::chrono::duration_cast<std::chrono::microseconds>
                ( tO2 - tO1 ).count() / 1e6;
std::cout << "Duration offline part = " << duration << " seconds" << std::endl;
forAll(linSys_solvers, solverI)
{
    Info << "Solver " << linSys_solvers[solverI] << endl;
    Info << endl;
    auto t1 = std::chrono::high_resolution_clock::now();
    example_paramBC.parameterizedBC(linSys_solvers[solverI], 6);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>
                    ( t2 - t1 ).count() / 1e6;
    std::cout << "Duration online part = " << duration << " seconds" << std::endl;
    volScalarField gParametrizedField = example_paramBC.list2Field(
                                            example_paramBC.g);
    ITHACAstream::exportSolution(gParametrizedField,
                                 std::to_string(solverI + 1),
                                 outputFolder,
                                 "gParametrized");
    volScalarField& T(example_paramBC._T());
    ITHACAstream::exportSolution(T,
                                 std::to_string(solverI + 1),
                                 outputFolder,
                                 "T");
    residualNorms(solverI) = Foam::sqrt(
                                 example_paramBC.residual.squaredNorm());
}
Eigen::MatrixXd A = example_paramBC.Theta.transpose() * example_paramBC.Theta;
ITHACAstream::exportMatrix(residualNorms, "residuals2norm", "eigen",
                           outputFolder);
example_paramBC.postProcess(outputFolder, "gParametrized");
Info << "*********************************************************" << endl;
Info << endl;
