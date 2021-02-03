Info << "Computing errors" << endl;
gField.resize(0);
ITHACAstream::read_fields(gField, heatFluxFieldName,
                          folder);
volScalarField gTrueField = list2Field(gTrue, innerField);
label Nsolutions = gField.size();
Eigen::MatrixXd heatFluxL2norm;
heatFluxL2norm.resize(Nsolutions, 1);
Eigen::MatrixXd heatFluxLinfNorm = heatFluxL2norm;
Eigen::VectorXd relativeErrorNorm_L2(Nsolutions);
Eigen::VectorXd relativeErrorNorm_Linf(Nsolutions);

forAll(gField, solutionI)
{
    volScalarField gDiffField = gField[solutionI] - gTrueField;
    ITHACAstream::exportSolution(gDiffField,
                                 std::to_string(solutionI + 1), folder,
                                 "gDiffField");
    ITHACAstream::exportSolution(gTrueField,
                                 std::to_string(solutionI + 1), folder,
                                 "gTrueField");
    relativeErrorNorm_L2(solutionI) = ITHACAutilities::L2normOnPatch(mesh,
                                      gDiffField,
                                      "hotSide") / ITHACAutilities::L2normOnPatch(mesh, gTrueField, "hotSide");
    relativeErrorNorm_Linf(solutionI) = ITHACAutilities::LinfNormOnPatch(mesh,
                                        gDiffField,
                                        "hotSide") / ITHACAutilities::LinfNormOnPatch(mesh, gTrueField, "hotSide");
    scalar EPS = 1e-6;
    volScalarField relativeErrorField(gDiffField);

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
                                 std::to_string(solutionI + 1), folder,
                                 "relativeErrorField");
    heatFluxL2norm(solutionI) = ITHACAutilities::L2normOnPatch(mesh,
                                relativeErrorField,
                                "hotSide");
    heatFluxLinfNorm(solutionI) = ITHACAutilities::LinfNormOnPatch(mesh,
                                  relativeErrorField,
                                  "hotSide");
}
ITHACAstream::exportMatrix(heatFluxL2norm, "relError_L2norm", "eigen",
                           folder);
ITHACAstream::exportMatrix(heatFluxLinfNorm, "relError_LinfNorm", "eigen",
                           folder);
ITHACAstream::exportMatrix(relativeErrorNorm_L2, "relativeErrorNorm_L2",
                           "eigen",
                           folder);
ITHACAstream::exportMatrix(relativeErrorNorm_Linf, "relativeErrorNorm_Linf",
                           "eigen",
                           folder);
