Info << endl;
Info << "*********************************************************" << endl;
Info << "Performing test for the parameterized BC inverse solver" << endl;
Info << endl;
word outputFolder = "./ITHACAoutput/parameterizedBCtest_RBFparameter/";
volScalarField gTrueField = example_paramBC.list2Field(example_paramBC.gTrue);
ITHACAstream::exportSolution(gTrueField,
                             "1", outputFolder,
                             "gTrue");

Eigen::VectorXd rbfWidth = EigenFunctions::ExpSpaced(0.001, 100,
                           rbfWidthTest_size);
ITHACAstream::exportMatrix(rbfWidth, "rbfShapeParameters", "eigen",
                           outputFolder);

Eigen::VectorXd residualNorms;
residualNorms.resize(rbfWidthTest_size);
Eigen::VectorXd heatFluxL2norm(rbfWidthTest_size);
Eigen::VectorXd heatFluxLinfNorm = heatFluxL2norm;
Eigen::VectorXd condNumber = heatFluxL2norm;
Eigen::MatrixXd singVal;

for (int i = 0; i < rbfWidthTest_size; i++)
{
    Info << "*********************************************************" << endl;
    Info << "RBF parameter " << rbfWidth(i) << endl;
    Info << "Test " << i + 1 << endl;
    Info << endl;
    example_paramBC.set_gParametrized("rbf", rbfWidth(i));
    example_paramBC.parameterizedBCoffline(1);
    example_paramBC.parameterizedBC("fullPivLU");
    volScalarField gParametrizedField = example_paramBC.list2Field(
                                            example_paramBC.g);
    ITHACAstream::exportSolution(gParametrizedField,
                                 std::to_string(i + 1),
                                 outputFolder,
                                 "gParametrized");
    volScalarField& T(example_paramBC._T());
    ITHACAstream::exportSolution(T,
                                 std::to_string(i + 1),
                                 outputFolder,
                                 "T");
    residualNorms(i) = Foam::sqrt(
                           example_paramBC.residual.squaredNorm());
    Eigen::MatrixXd A = example_paramBC.Theta.transpose() * example_paramBC.Theta;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A,
                                          Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd singularValues = svd.singularValues();
    singVal.conservativeResize(singularValues.rows(), singVal.cols() + 1);
    singVal.col(i) = singularValues;
    double conditionNumber;

    if (singularValues.minCoeff() > 0)
    {
        conditionNumber = singularValues.maxCoeff() / singularValues.minCoeff();
    }
    else
    {
        conditionNumber = DBL_MAX;
    }

    Info << "Condition number = " << conditionNumber << endl;
    condNumber(i) = conditionNumber;
    volScalarField gDiffField = gParametrizedField - gTrueField;
    scalar EPS = 1e-6;
    volScalarField relativeErrorField(gTrueField);

    for (label i = 0; i < relativeErrorField.internalField().size(); i++)
    {
        if (std::abs(gTrueField.ref()[i]) < EPS)
        {
            relativeErrorField.ref()[i] = (std::abs(gDiffField.ref()[i])) / EPS;
        }
        else
        {
            relativeErrorField.ref()[i] = (std::abs(gDiffField.ref()[i])) /
                                          gTrueField.ref()[i];
        }
    }

    ITHACAstream::exportSolution(relativeErrorField,
                                 std::to_string(i + 1), outputFolder,
                                 "relativeErrorField");
    heatFluxL2norm(i) = ITHACAutilities::L2normOnPatch(mesh, relativeErrorField,
                        "hotSide");
    heatFluxLinfNorm(i) = ITHACAutilities::LinfNormOnPatch(mesh, relativeErrorField,
                          "hotSide");
}

ITHACAstream::exportMatrix(condNumber, "condNumber", "eigen",
                           outputFolder);
ITHACAstream::exportMatrix(heatFluxL2norm, "relError_L2norm", "eigen",
                           outputFolder);
ITHACAstream::exportMatrix(heatFluxLinfNorm, "relError_LinfNorm", "eigen",
                           outputFolder);
ITHACAstream::exportMatrix(singVal, "singularValues", "eigen",
                           outputFolder);
ITHACAstream::exportMatrix(residualNorms, "residuals2norm", "eigen",
                           outputFolder);
example_paramBC.postProcess(outputFolder, "gParametrized");
Info << "*********************************************************" << endl;
Info << "*********************************************************" << endl;
Info << endl;
