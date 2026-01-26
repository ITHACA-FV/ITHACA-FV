IOdictionary regularizationDict
(
    IOobject
    (
        "regularizationDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

linSys_solver =
    regularizationDict.lookupOrDefault<word>("regularizationTechnique",
        "fullPivLU");
TSVD_filter = regularizationDict.lookupOrDefault<int>("TSVD_filter", 0);
Tikhonov_filter = regularizationDict.lookupOrDefault<double>("Tikhonov_filter",
                  0.0);
CG_Nsteps = regularizationDict.lookupOrDefault<int>("CG_Nsteps", 0);
