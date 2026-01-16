/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------

  License
  This file is part of ITHACA-FV

  ITHACA-FV is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ITHACA-FV is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "UnsteadyNSTurb.H"

/// \file
/// Source file of the unsteadyNS class.

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Construct Null
UnsteadyNSTurb::UnsteadyNSTurb() {}

// Construct from zero
UnsteadyNSTurb::UnsteadyNSTurb(int argc, char* argv[])
{
    _args = autoPtr<argList>
            (
                new argList(argc, argv)
            );

    if (!_args->checkRootCase())
    {
        Foam::FatalError.exit();
    }

    argList& args = _args();
#include "createTime.H"
#include "createMesh.H"
    _pimple = autoPtr<pimpleControl>
              (
                  new pimpleControl
                  (
                      mesh
                  )
              );
#include "createFields.H"
#include "createUfIfPresent.H"
#include "createFvOptions.H"
    ITHACAdict = new IOdictionary
    (
        IOobject
        (
            "ITHACAdict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    tolerance = ITHACAdict->lookupOrDefault<scalar>("tolerance", 1e-5);
    maxIter = ITHACAdict->lookupOrDefault<scalar>("maxIter", 1000);
    bcMethod = ITHACAdict->lookupOrDefault<word>("bcMethod", "lift");
    timeDerivativeSchemeOrder =
        ITHACAdict->lookupOrDefault<word>("timeDerivativeSchemeOrder", "second");
    M_Assert(bcMethod == "lift" || bcMethod == "penalty",
             "The BC method must be set to lift or penalty in ITHACAdict");
    M_Assert(timeDerivativeSchemeOrder == "first"
             || timeDerivativeSchemeOrder == "second",
                                          "The time derivative approximation must be set to either first or second order scheme in ITHACAdict");
    para = ITHACAparameters::getInstance(mesh, runTime);
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    supex = ITHACAutilities::check_sup();
    /// Number of velocity modes to be calculated
    NUmodesOut = para->ITHACAdict->lookupOrDefault<label>("NmodesUout", 15);
    /// Number of pressure modes to be calculated
    NPmodesOut = para->ITHACAdict->lookupOrDefault<label>("NmodesPout", 15);
    /// Number of nut modes to be calculated
    NNutModesOut = para->ITHACAdict->lookupOrDefault<label>("NmodesNutOut", 15);
    /// Number of velocity modes used for the projection
    NUmodes = para->ITHACAdict->lookupOrDefault<label>("NmodesUproj", 10);
    /// Number of supremizers modes used for the projection
    NSUPmodes = para->ITHACAdict->lookupOrDefault<label>("NmodesSUPproj", 10);
    /// Number of pressure modes used for the projection
    NPmodes = para->ITHACAdict->lookupOrDefault<label>("NmodesPproj", 10);
    /// Number of nut modes used for the projection
    NNutModes = para->ITHACAdict->lookupOrDefault<label>("NmodesNutProj", 0);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void UnsteadyNSTurb::truthSolve(List<scalar> mu_now, std::string& offlinepath)
{
    Time& runTime = _runTime();
    surfaceScalarField& phi = _phi();
    fvMesh& mesh = _mesh();
#include "initContinuityErrs.H"
    fv::options& fvOptions = _fvOptions();
    pimpleControl& pimple = _pimple();
    volScalarField& p = _p();
    volVectorField& U = _U();
    volScalarField& nut = _nut();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    instantList Times = runTime.times();
    label& pRefCell = _pRefCell;
    scalar& pRefValue = _pRefValue;
    mesh.setFluxRequired(p.name());
    runTime.setEndTime(finalTime);
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime + writeEvery;
    label nsnapshots = 0;

    // Start the time loop
    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        ++runTime;
        Info << "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
#include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
#include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;

        if (checkWrite(runTime))
        {
            nsnapshots += 1;
            // Produces error when uncommented
            // volScalarField nut = turbulence->nut().ref();
            nut = turbulence->nut();
            ITHACAstream::exportSolution(U, name(counter), offlinepath);
            ITHACAstream::exportSolution(p, name(counter), offlinepath);
            ITHACAstream::exportSolution(nut, name(counter), offlinepath);
            std::ofstream of(offlinepath + name(counter) + "/" +
                             runTime.timeName());
            Ufield.append(tmp<volVectorField>(U));
            Pfield.append(tmp<volScalarField>(p));
            nutFields.append(tmp<volScalarField>(nut));
            counter++;
            nextWrite += writeEvery;
            writeMu(mu_now);
            // --- Fill in the mu_samples with parameters (time, mu) to be used for the PODI sample points
            mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size() + 1);
            mu_samples(mu_samples.rows() - 1, 0) = atof(runTime.timeName().c_str());

            for (label i = 0; i < mu_now.size(); i++)
            {
                mu_samples(mu_samples.rows() - 1, i + 1) = mu_now[i];
            }
        }
    }

    // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
    if (mu.cols() == 0)
    {
        mu.resize(1, 1);
    }

    if (mu_samples.rows() == nsnapshots * mu.cols())
    {
        ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
                                   offlinepath);
    }
}


// ====== SUP Full Tensor 1 ======
Eigen::Tensor<double, 3>
UnsteadyNSTurb::turbulenceTensor1(label NUmodes, label NSUPmodes, label nNutModes)
{
    const label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct1Tensor(cSize, nNutModes, cSize);

    for (label i = 0; i < cSize; ++i)
        for (label j = 0; j < nNutModes; ++j)
            for (label k = 0; k < cSize; ++k)
                ct1Tensor(i, j, k) =
                    fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::laplacian(nutModes[j], L_U_SUPmodes[k])).value();

    ITHACAstream::SaveDenseTensor(ct1Tensor, "./ITHACAoutput/Matrices/",
        "ct1_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(nNutModes) + "_t");
    return ct1Tensor;
}

// ====== SUP Average Tensor 1 ======
Eigen::Tensor<double, 3>
UnsteadyNSTurb::turbulenceAveTensor1(label NUmodes, label NSUPmodes)
{
    const label cSize = NUmodes + NSUPmodes + liftfield.size();
    const label nAvg  = nutAve.size();
    Eigen::Tensor<double, 3> ct1AveTensor(cSize, nAvg, cSize);

    for (label i = 0; i < cSize; ++i)
        for (label j = 0; j < nAvg;  ++j)
            for (label k = 0; k < cSize; ++k)
                ct1AveTensor(i, j, k) =
                    fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::laplacian(nutAve[j], L_U_SUPmodes[k])).value();

    ITHACAstream::SaveDenseTensor(ct1AveTensor, "./ITHACAoutput/Matrices/",
        "ct1Ave_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_t");
    return ct1AveTensor;
}

// ====== SUP Fluctuation Tensor 1 ======
Eigen::Tensor<double, 3>
UnsteadyNSTurb::turbulenceFluctTensor1(label NUmodes, label NSUPmodes)
{
    const label cSize  = NUmodes + NSUPmodes + liftfield.size();
    const label nFluct = nutFluctModes.size();
    Eigen::Tensor<double, 3> ct1FluctTensor(cSize, nFluct, cSize);

    for (label i = 0; i < cSize; ++i)
        for (label j = 0; j < nFluct; ++j)
            for (label k = 0; k < cSize; ++k)
                ct1FluctTensor(i, j, k) =
                    fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::laplacian(nutFluctModes[j], L_U_SUPmodes[k])).value();

    ITHACAstream::SaveDenseTensor(ct1FluctTensor, "./ITHACAoutput/Matrices/",
        "ct1Fluct_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_t");
    return ct1FluctTensor;
}

// ====== PPE Full Tensor 1 ======
Eigen::Tensor<double, 3>
UnsteadyNSTurb::turbulencePPETensor1(label NUmodes, label NSUPmodes, label NPmodes, label nNutModes)
{
    const label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct1PPETensor(NPmodes, nNutModes, cSize);

    for (label i = 0; i < NPmodes; ++i)
        for (label j = 0; j < nNutModes; ++j)
            for (label k = 0; k < cSize; ++k)
                ct1PPETensor(i, j, k) =
                    fvc::domainIntegrate(fvc::grad(Pmodes[i]) & (fvc::laplacian(nutModes[j], L_U_SUPmodes[k]))).value();

    ITHACAstream::SaveDenseTensor(ct1PPETensor, "./ITHACAoutput/Matrices/",
        "ct1PPE_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPmodes) + "_" + name(nNutModes) + "_t");
    return ct1PPETensor;
}

// ====== PPE Average Tensor 1 ======
Eigen::Tensor<double, 3>
UnsteadyNSTurb::turbulencePPEAveTensor1(label NUmodes, label NSUPmodes, label NPmodes)
{
    const label cSize = NUmodes + NSUPmodes + liftfield.size();
    const label nAvg  = nutAve.size();
    Eigen::Tensor<double, 3> ct1PPEAveTensor(NPmodes, nAvg, cSize);

    for (label i = 0; i < NPmodes; ++i)
        for (label j = 0; j < nAvg;  ++j)
            for (label k = 0; k < cSize; ++k)
                ct1PPEAveTensor(i, j, k) =
                    fvc::domainIntegrate(fvc::grad(Pmodes[i]) & (fvc::laplacian(nutAve[j], L_U_SUPmodes[k]))).value();

    ITHACAstream::SaveDenseTensor(ct1PPEAveTensor, "./ITHACAoutput/Matrices/",
        "ct1PPEAve_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPmodes) + "_t");
    return ct1PPEAveTensor;
}

// ====== PPE Fluctuation Tensor 1 ======
Eigen::Tensor<double, 3>
UnsteadyNSTurb::turbulencePPEFluctTensor1(label NUmodes, label NSUPmodes, label NPmodes)
{
    const label cSize  = NUmodes + NSUPmodes + liftfield.size();
    const label nFluct = nutFluctModes.size();
    Eigen::Tensor<double, 3> ct1PPEFluctTensor(NPmodes, nFluct, cSize);

    for (label i = 0; i < NPmodes; ++i)
        for (label j = 0; j < nFluct; ++j)
            for (label k = 0; k < cSize; ++k)
                ct1PPEFluctTensor(i, j, k) =
                    fvc::domainIntegrate(fvc::grad(Pmodes[i]) & (fvc::laplacian(nutFluctModes[j], L_U_SUPmodes[k]))).value();

    ITHACAstream::SaveDenseTensor(ct1PPEFluctTensor, "./ITHACAoutput/Matrices/",
        "ct1PPEFluct_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPmodes) + "_t");
    return ct1PPEFluctTensor;
}

// ====== SUP Full Tensor 2 ======
Eigen::Tensor<double, 3>
UnsteadyNSTurb::turbulenceTensor2(label NUmodes, label NSUPmodes, label nNutModes)
{
    const label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct2Tensor(cSize, nNutModes, cSize);

    for (label i = 0; i < cSize; ++i)
        for (label j = 0; j < nNutModes; ++j)
            for (label k = 0; k < cSize; ++k)
                ct2Tensor(i, j, k) = fvc::domainIntegrate(
                    L_U_SUPmodes[i] &
                    (fvc::div(nutModes[j] * dev2((fvc::grad(L_U_SUPmodes[k]))().T())))
                ).value();

    ITHACAstream::SaveDenseTensor(ct2Tensor, "./ITHACAoutput/Matrices/",
        "ct2_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(nNutModes) + "_t");
    return ct2Tensor;
}

// ====== SUP Average Tensor 2 ======
Eigen::Tensor<double, 3>
UnsteadyNSTurb::turbulenceAveTensor2(label NUmodes, label NSUPmodes)
{
    const label cSize = NUmodes + NSUPmodes + liftfield.size();
    const label nAvg  = nutAve.size();
    Eigen::Tensor<double, 3> ct2AveTensor(cSize, nAvg, cSize);

    for (label i = 0; i < cSize; ++i)
        for (label j = 0; j < nAvg;  ++j)
            for (label k = 0; k < cSize; ++k)
                ct2AveTensor(i, j, k) = fvc::domainIntegrate(
                    L_U_SUPmodes[i] &
                    (fvc::div(nutAve[j] * dev2((fvc::grad(L_U_SUPmodes[k]))().T())))
                ).value();

    ITHACAstream::SaveDenseTensor(ct2AveTensor, "./ITHACAoutput/Matrices/",
        "ct2Ave_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_t");
    return ct2AveTensor;
}

// ====== SUP Fluctuation Tensor 2 ======
Eigen::Tensor<double, 3>
UnsteadyNSTurb::turbulenceFluctTensor2(label NUmodes, label NSUPmodes)
{
    const label cSize  = NUmodes + NSUPmodes + liftfield.size();
    const label nFluct = nutFluctModes.size();
    Eigen::Tensor<double, 3> ct2FluctTensor(cSize, nFluct, cSize);

    for (label i = 0; i < cSize; ++i)
        for (label j = 0; j < nFluct; ++j)
            for (label k = 0; k < cSize; ++k)
                ct2FluctTensor(i, j, k) = fvc::domainIntegrate(
                    L_U_SUPmodes[i] &
                    (fvc::div(nutFluctModes[j] * dev2((fvc::grad(L_U_SUPmodes[k]))().T())))
                ).value();

    ITHACAstream::SaveDenseTensor(ct2FluctTensor, "./ITHACAoutput/Matrices/",
        "ct2Fluct_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_t");
    return ct2FluctTensor;
}

// ====== PPE Full Tensor 2 ======
Eigen::Tensor<double, 3>
UnsteadyNSTurb::turbulencePPETensor2(label NUmodes, label NSUPmodes, label NPmodes, label nNutModes)
{
    const label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct2PPETensor(NPmodes, nNutModes, cSize);

    for (label i = 0; i < NPmodes; ++i)
        for (label j = 0; j < nNutModes; ++j)
            for (label k = 0; k < cSize; ++k)
                ct2PPETensor(i, j, k) = fvc::domainIntegrate(
                    fvc::grad(Pmodes[i]) &
                    (fvc::div(nutModes[j] * dev2((fvc::grad(L_U_SUPmodes[k]))().T())))
                ).value();

    ITHACAstream::SaveDenseTensor(ct2PPETensor, "./ITHACAoutput/Matrices/",
        "ct2PPE_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPmodes) + "_" + name(nNutModes) + "_t");
    return ct2PPETensor;
}

// ====== PPE Average Tensor 2 ======
Eigen::Tensor<double, 3>
UnsteadyNSTurb::turbulencePPEAveTensor2(label NUmodes, label NSUPmodes, label NPmodes)
{
    const label cSize = NUmodes + NSUPmodes + liftfield.size();
    const label nAvg  = nutAve.size();
    Eigen::Tensor<double, 3> ct2PPEAveTensor(NPmodes, nAvg, cSize);

    for (label i = 0; i < NPmodes; ++i)
        for (label j = 0; j < nAvg;  ++j)
            for (label k = 0; k < cSize; ++k)
                ct2PPEAveTensor(i, j, k) = fvc::domainIntegrate(
                    fvc::grad(Pmodes[i]) &
                    (fvc::div(nutAve[j] * dev2((fvc::grad(L_U_SUPmodes[k]))().T())))
                ).value();

    ITHACAstream::SaveDenseTensor(ct2PPEAveTensor, "./ITHACAoutput/Matrices/",
        "ct2PPEAve_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPmodes) + "_t");
    return ct2PPEAveTensor;
}

// ====== PPE Fluctuation Tensor 2 ======
Eigen::Tensor<double, 3>
UnsteadyNSTurb::turbulencePPEFluctTensor2(label NUmodes, label NSUPmodes, label NPmodes)
{
    const label cSize  = NUmodes + NSUPmodes + liftfield.size();
    const label nFluct = nutFluctModes.size();
    Eigen::Tensor<double, 3> ct2PPEFluctTensor(NPmodes, nFluct, cSize);

    for (label i = 0; i < NPmodes; ++i)
        for (label j = 0; j < nFluct; ++j)
            for (label k = 0; k < cSize; ++k)
                ct2PPEFluctTensor(i, j, k) = fvc::domainIntegrate(
                    fvc::grad(Pmodes[i]) &
                    (fvc::div(nutFluctModes[j] * dev2((fvc::grad(L_U_SUPmodes[k]))().T())))
                ).value();

    ITHACAstream::SaveDenseTensor(ct2PPEFluctTensor, "./ITHACAoutput/Matrices/",
        "ct2PPEFluct_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPmodes) + "_t");
    return ct2PPEFluctTensor;
}

// ====== (unchanged) BT Turbulence Matrix (SUP) ======
Eigen::MatrixXd
UnsteadyNSTurb::btTurbulence(label NUmodes, label NSUPmodes)
{
    const label btSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd btMatrix(btSize, btSize);
    btMatrix = btMatrix * 0;

    for (label i = 0; i < btSize; ++i)
        for (label j = 0; j < btSize; ++j)
            btMatrix(i, j) =
                fvc::domainIntegrate(L_U_SUPmodes[i] &
                    (fvc::div(dev2((T(fvc::grad(L_U_SUPmodes[j])))))))
                .value();

    ITHACAstream::SaveDenseMatrix(btMatrix, "./ITHACAoutput/Matrices/",
        "bt_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));
    return btMatrix;
}

// ====== NEW: SUP Continuity Matrix P: P_{ij} = (χ_i, div(φ_j)) ======
Eigen::MatrixXd
UnsteadyNSTurb::continuity_matrix(label NUmodes, label NSUPmodes, label NPmodes)
{
    const label cSize = NUmodes + NSUPmodes + liftfield.size();  // φ_j basis: [lift, U, SUP]
    M_Assert(cSize > 0,   "continuity_matrix: cSize=0");
    M_Assert(NPmodes > 0, "continuity_matrix: NPmodes=0");

    Eigen::MatrixXd Pmat(NPmodes, cSize);
    Pmat.setZero();

    for (label i = 0; i < NPmodes; ++i)       // χ_i
        for (label j = 0; j < cSize;   ++j)   // φ_j
            Pmat(i, j) =
                fvc::domainIntegrate(Pmodes[i] * fvc::div(L_U_SUPmodes[j])).value();

    ITHACAstream::SaveDenseMatrix(Pmat, "./ITHACAoutput/Matrices/",
        "P_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPmodes));
    return Pmat;
}

// ====== PPE L vector (boundary RHS): L_i = ∫_Γ χ_i (n · R_t) dΓ ======
Eigen::MatrixXd
UnsteadyNSTurb::pressurePPE_L(label NPmodes)
{
    // Result is a column vector (NPmodes x 1)
    Eigen::MatrixXd Lvec(NPmodes, 1);
    Lvec.setZero();

    // --- get dt at the beginning of the first μ-block (same log used offline)
    const word muPath = "./ITHACAoutput/Offline/mu_samples_mat.txt";
    M_Assert(ITHACAutilities::check_file(muPath),
             "[PPE-L] mu_samples_mat.txt not found; cannot compute R_t.");
    Eigen::MatrixXd muMat = ITHACAstream::readMatrix(muPath);
    M_Assert(muMat.rows() >= 2, "[PPE-L] mu/time log has <2 rows; need two snapshots for R_t.");
    const double dt0 = muMat(1, 0) - muMat(0, 0);
    M_Assert(std::abs(dt0) > SMALL, "[PPE-L] dt0 ~ 0; check mu_samples_mat.txt.");

    // --- R_t ≈ (U(t1) - U(t0)) / dt0 from the first two FOM snapshots
    M_Assert(Uomfield.size() >= 2, "[PPE-L] Need at least 2 velocity snapshots to build R_t.");
    const volVectorField& U0 = Uomfield[0];
    const volVectorField& U1 = Uomfield[1];
    tmp<volVectorField> tUdot = (U1 - U0) / dt0;
    const volVectorField& Udot = tUdot();

    // --- boundary integration: ∫_Γ χ_i (Udot · n) dΓ = Σ_faces χ_i (Udot · S_f)
    const fvMesh& mesh = Pmodes[0]().mesh();

    for (label i = 0; i < NPmodes; ++i)
    {
        scalar Li = 0.0;
        forAll(mesh.boundary(), patchi)
        {
            const scalarField& chiF  = Pmodes[i].boundaryField()[patchi];
            const vectorField& UdotF = Udot.boundaryField()[patchi];
            const vectorField& Sf    = mesh.Sf().boundaryField()[patchi];
            Li += gSum( chiF * (UdotF & Sf) );
        }
        Lvec(i, 0) = Li;
    }

    ITHACAstream::SaveDenseMatrix(Lvec, "./ITHACAoutput/Matrices/",
        "L_" + name(NPmodes));
    return Lvec;
}


void UnsteadyNSTurb::projectSUP(fileName folder,
                                label NU, label NP, label NSUP,
                                label Nnut, bool rbfInterp)
{
    std::cout << "[DEBUG] Entered projectSUP()" << std::endl;

    // --- set sizes ---
    NUmodes    = NU;
    NPmodes    = NP;
    NSUPmodes  = NSUP;
    nNutModes  = Nnut;

    // --- build concatenated basis φ_j = [lift, U, SUP] used in all SUP integrals ---
    L_U_SUPmodes.resize(0);
    if (liftfield.size() != 0)
        for (label k = 0; k < liftfield.size(); ++k)
            L_U_SUPmodes.append(tmp<volVectorField>(liftfield[k]));
    if (NUmodes != 0)
        for (label k = 0; k < NUmodes; ++k)
            L_U_SUPmodes.append(tmp<volVectorField>(Umodes[k]));
    if (NSUPmodes != 0)
        for (label k = 0; k < NSUPmodes; ++k)
            L_U_SUPmodes.append(tmp<volVectorField>(supmodes[k]));

    // Convenience
    const label cSize = liftfield.size() + NUmodes + NSUPmodes;

    // ---------- Read-if-exists, otherwise assemble (SUP-only objects) ----------
    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {
        // --- Mass (SUP M) ---
        {
            word M_str = "M_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes);
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + M_str))
                ITHACAstream::ReadDenseMatrix(M_matrix, "./ITHACAoutput/Matrices/", M_str);
            else
                M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes); // uses SUP basis
        }

        // --- Diffusion matrix (SUP B) ---
        {
            word B_str = "B_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes);
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + B_str))
                ITHACAstream::ReadDenseMatrix(B_matrix, "./ITHACAoutput/Matrices/", B_str);
            else
                B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        }

        // --- Turbulent correction to diffusion (bt) ---
        {
            word btStr = "bt_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes);
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + btStr))
                ITHACAstream::ReadDenseMatrix(btMatrix, "./ITHACAoutput/Matrices/", btStr);
            else
                btMatrix = btTurbulence(NUmodes, NSUPmodes);
        }

        // --- Pressure gradient coupling (K) ---
        {
            word K_str = "K_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPmodes);
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + K_str))
                ITHACAstream::ReadDenseMatrix(K_matrix, "./ITHACAoutput/Matrices/", K_str);
            else
                K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        }

        // --- NEW: Continuity matrix (P) ---
        {
            word P_str = "P_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPmodes);
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + P_str))
                ITHACAstream::ReadDenseMatrix(P_matrix, "./ITHACAoutput/Matrices/", P_str);
            else
                P_matrix = continuity_matrix(NUmodes, NSUPmodes, NPmodes); // P: (NP x cSize)
        }

        // --- Convective tensor (C) ---
        {
            word C_str = "C_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_t";
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + C_str))
                ITHACAstream::ReadDenseTensor(C_tensor, "./ITHACAoutput/Matrices/", C_str);
            else
                C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        }

        // --- Eddy-viscosity tensors (ct1, ct2) ---
        {
            word ct1Str = "ct1_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(nNutModes) + "_t";
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct1Str))
                ITHACAstream::ReadDenseTensor(ct1Tensor, "./ITHACAoutput/Matrices/", ct1Str);
            else
                ct1Tensor = turbulenceTensor1(NUmodes, NSUPmodes, nNutModes);

            word ct2Str = "ct2_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(nNutModes) + "_t";
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct2Str))
                ITHACAstream::ReadDenseTensor(ct2Tensor, "./ITHACAoutput/Matrices/", ct2Str);
            else
                ct2Tensor = turbulenceTensor2(NUmodes, NSUPmodes, nNutModes);
        }

        // --- Optional: Average (μ-linear) and Fluctuation (RBF) parts for ν_t ---
        if (nutAve.size() != 0)
        {
            word ct1AveStr = "ct1Ave_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_t";
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct1AveStr))
                ITHACAstream::ReadDenseTensor(ct1AveTensor, "./ITHACAoutput/Matrices/", ct1AveStr);
            else
                ct1AveTensor = turbulenceAveTensor1(NUmodes, NSUPmodes);

            word ct2AveStr = "ct2Ave_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_t";
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct2AveStr))
                ITHACAstream::ReadDenseTensor(ct2AveTensor, "./ITHACAoutput/Matrices/", ct2AveStr);
            else
                ct2AveTensor = turbulenceAveTensor2(NUmodes, NSUPmodes);
        }

        if (nutFluctModes.size() != 0)
        {
            word ct1FluctStr = "ct1Fluct_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_t";
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct1FluctStr))
                ITHACAstream::ReadDenseTensor(ct1FluctTensor, "./ITHACAoutput/Matrices/", ct1FluctStr);
            else
                ct1FluctTensor = turbulenceFluctTensor1(NUmodes, NSUPmodes);

            word ct2FluctStr = "ct2Fluct_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_t";
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct2FluctStr))
                ITHACAstream::ReadDenseTensor(ct2FluctTensor, "./ITHACAoutput/Matrices/", ct2FluctStr);
            else
                ct2FluctTensor = turbulenceFluctTensor2(NUmodes, NSUPmodes);
        }

        // --- Optional velocity BC penalty assembly ---
        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }
    }
    else
    {
        // No cache folder: assemble everything needed for SUP
        M_matrix  = mass_term(NUmodes, NPmodes, NSUPmodes);
        B_matrix  = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        btMatrix  = btTurbulence(NUmodes, NSUPmodes);
        K_matrix  = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        P_matrix  = continuity_matrix(NUmodes, NSUPmodes, NPmodes);
        C_tensor  = convective_term_tens(NUmodes, NPmodes, NSUPmodes);

        ct1Tensor = turbulenceTensor1(NUmodes, NSUPmodes, nNutModes);
        ct2Tensor = turbulenceTensor2(NUmodes, NSUPmodes, nNutModes);

        if (nutAve.size() != 0)
        {
            ct1AveTensor = turbulenceAveTensor1(NUmodes, NSUPmodes);
            ct2AveTensor = turbulenceAveTensor2(NUmodes, NSUPmodes);
        }
        if (nutFluctModes.size() != 0)
        {
            ct1FluctTensor = turbulenceFluctTensor1(NUmodes, NSUPmodes);
            ct2FluctTensor = turbulenceFluctTensor2(NUmodes, NSUPmodes);
        }

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }
    }

    // ---------- Combine linear diffusion (physical + turbulent) ----------
    bTotalMatrix = B_matrix + btMatrix;  // size: cSize x cSize

    // ---------- Combine ν_t tensors ----------
    cTotalTensor.resize(cSize, nNutModes, cSize);
    cTotalTensor = ct1Tensor + ct2Tensor;

    if (nutAve.size() != 0)
    {
        cTotalAveTensor.resize(cSize, nutAve.size(), cSize);
        cTotalAveTensor = ct1AveTensor + ct2AveTensor;
    }
    if (nutFluctModes.size() != 0)
    {
        cTotalFluctTensor.resize(cSize, nutFluctModes.size(), cSize);
        cTotalFluctTensor = ct1FluctTensor + ct2FluctTensor;
    }

    // ------------ Export block (Python, if enabled) ------------
    if (para->exportPython)
    {
        ITHACAstream::exportMatrix(M_matrix,   "M",   "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(B_matrix,   "B",   "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(btMatrix,   "bt",  "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(bTotalMatrix,"bTot","python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix,   "K",   "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(P_matrix,   "P",   "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor,   "C",   "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct1Tensor,  "ct1", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct2Tensor,  "ct2", "python", "./ITHACAoutput/Matrices/");
        if (nutAve.size() != 0)
        {
            ITHACAstream::exportTensor(ct1AveTensor, "ct1Ave", "python", "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2AveTensor, "ct2Ave", "python", "./ITHACAoutput/Matrices/");
        }
        if (nutFluctModes.size() != 0)
        {
            ITHACAstream::exportTensor(ct1FluctTensor, "ct1Fluct", "python", "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2FluctTensor, "ct2Fluct", "python", "./ITHACAoutput/Matrices/");
        }
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportMatrix(M_matrix,   "M",   "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(B_matrix,   "B",   "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(btMatrix,   "bt",  "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(bTotalMatrix,"bTot","matlab","./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix,   "K",   "matlab","./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(P_matrix,   "P",   "matlab","./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor,   "C",   "matlab","./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct1Tensor,  "ct1", "matlab","./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct2Tensor,  "ct2", "matlab","./ITHACAoutput/Matrices/");
        if (nutAve.size() != 0)
        {
            ITHACAstream::exportTensor(ct1AveTensor, "ct1Ave", "matlab", "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2AveTensor, "ct2Ave", "matlab", "./ITHACAoutput/Matrices/");
        }
        if (nutFluctModes.size() != 0)
        {
            ITHACAstream::exportTensor(ct1FluctTensor, "ct1Fluct", "matlab", "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2FluctTensor, "ct2Fluct", "matlab", "./ITHACAoutput/Matrices/");
        }
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportMatrix(M_matrix,   "M",   "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(B_matrix,   "B",   "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(btMatrix,   "bt",  "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(bTotalMatrix,"bTot","eigen","./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix,   "K",   "eigen","./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(P_matrix,   "P",   "eigen","./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor,   "C",   "eigen","./ITHACAoutput/Matrices/C");
        ITHACAstream::exportTensor(ct1Tensor,  "ct1_","eigen","./ITHACAoutput/Matrices/ct1");
        ITHACAstream::exportTensor(ct2Tensor,  "ct2_","eigen","./ITHACAoutput/Matrices/ct2");
        if (nutAve.size() != 0)
        {
            ITHACAstream::exportTensor(ct1AveTensor, "ct1Ave_", "eigen", "./ITHACAoutput/Matrices/ct1Ave");
            ITHACAstream::exportTensor(ct2AveTensor, "ct2Ave_", "eigen", "./ITHACAoutput/Matrices/ct2Ave");
        }
        if (nutFluctModes.size() != 0)
        {
            ITHACAstream::exportTensor(ct1FluctTensor, "ct1Fluct", "eigen", "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2FluctTensor, "ct2Fluct", "eigen", "./ITHACAoutput/Matrices/");
        }
    }

    // --- quick sanity logs ---
    std::cout << "[DEBUG] SUP sizes: cSize=" << cSize
              << ", NP=" << NPmodes << ", nNut=" << nNutModes << '\n';
    std::cout << "[DEBUG] Shapes: M(" << M_matrix.rows() << "x" << M_matrix.cols()
              << ") B("   << B_matrix.rows() << "x" << B_matrix.cols()
              << ") bt("  << btMatrix.rows() << "x" << btMatrix.cols()
              << ") K("   << K_matrix.rows() << "x" << K_matrix.cols()
              << ") P("   << P_matrix.rows() << "x" << P_matrix.cols() << ")\n";
    std::cout << "[DEBUG] projectSUP() done. SUP matrices/tensors exported.\n";



// ---- OFFLINE RBF (SUP variant): AVG linear-μ + FLUCT RBF, weights -> weightsSUP ----
if (rbfInterp && (!Pstream::parRun()))
{
    std::cout << "\n==== [RBF DEBUG SUP] ENTERING OFFLINE RBF CONSTRUCTION (AVG linear-μ + FLUCT RBF) ====\n";

    const word coeffDir = "./ITHACAoutput/Coefficients/";
    const word muPath   = "./ITHACAoutput/Offline/mu_samples_mat.txt";

    // --- 1. LOAD mu and time vectors ---
    Eigen::MatrixXd muMat = ITHACAstream::readMatrix(muPath);
    std::cout << "[RBF DEBUG SUP] muMat shape: " << muMat.rows() << " x " << muMat.cols() << std::endl;
    const int nPar = 12, nSnapshotsPerPar = 200;
    Eigen::VectorXd timeVec = muMat.col(0), muVec = muMat.col(1);

    // --- 2. LOAD avg/fluct νt coeffs ---
    Eigen::MatrixXd coeffNutAvg   = ITHACAstream::readMatrix(coeffDir + "Nut_avg_coeffs_mat.txt");
    Eigen::MatrixXd coeffNutFluct = ITHACAstream::readMatrix(coeffDir + "Nut_fluct_coeffs_mat.txt");
    std::cout << "[RBF DEBUG SUP] Loaded coeffNutAvg shape: "   << coeffNutAvg.rows()   << " x " << coeffNutAvg.cols()   << std::endl;
    std::cout << "[RBF DEBUG SUP] Loaded coeffNutFluct shape: " << coeffNutFluct.rows() << " x " << coeffNutFluct.cols() << std::endl;
    int nNutAvgModes   = coeffNutAvg.rows();
    int nNutFluctModes = coeffNutFluct.rows();
    int nUniqueMu      = nPar;

    // --- 3. Velocity projection coefficients (PHYSICAL modes only) ---
    Eigen::MatrixXd a = ITHACAutilities::getCoeffs(Uomfield, Umodes);
    a.transposeInPlace(); // [nSnapshots, nVelModes]

    // --- 4. Snapshot grouping per parameter ---
    Eigen::VectorXd initSnapInd(nPar), timeSnap(nPar);
    for (int i = 0; i < nPar; ++i) {
        int start = i * nSnapshotsPerPar;
        initSnapInd(i) = start;
        timeSnap(i) = timeVec(start + 1) - timeVec(start);
        std::cout << "[RBF DEBUG SUP] i=" << i << ", start=" << start << ", dt=" << timeSnap(i) << std::endl;
    }

    // --- 5. Unique μ for AVG part ---
    Eigen::VectorXd muVecUnique(nUniqueMu);
    for (int i = 0; i < nUniqueMu; ++i)
        muVecUnique(i) = muVec(static_cast<int>(initSnapInd(i)));
    std::cout << "[RBF DEBUG SUP] muVecUnique: [" << muVecUnique(0) << " ... " << muVecUnique(muVecUnique.size()-1)
              << "] (M=" << muVecUnique.size() << ")\n";

    // --- 6. Build fluctuation interpolation input: [a, aDot] ---
    std::cout << "[RBF DEBUG SUP] Calling velDerivativeCoeff() for fluctuation part...\n";
    Eigen::MatrixXd Gfluct = coeffNutFluct.transpose();
    List<Eigen::MatrixXd> interpDataFluct = velDerivativeCoeff(a, Gfluct, initSnapInd, timeSnap);
    Eigen::MatrixXd velRBF_fluct = interpDataFluct[0];  // [a, aDot]
    Eigen::MatrixXd coeffs_fluct = interpDataFluct[1];
    std::cout << "[RBF DEBUG SUP] velRBF_fluct shape: " << velRBF_fluct.rows() << " x " << velRBF_fluct.cols() << std::endl;
    std::cout << "[RBF DEBUG SUP] coeffs_fluct shape: " << coeffs_fluct.rows() << " x " << coeffs_fluct.cols() << std::endl;

    // --- 7. Optional debug exports ---
    int nSnapshots = velRBF_fluct.rows();
    int nModes = velRBF_fluct.cols() / 2;
    Eigen::MatrixXd adot_only(nSnapshots, nModes);
    for (int i = 0; i < nSnapshots; ++i)
        adot_only.row(i) = velRBF_fluct.row(i).tail(nModes);
    ITHACAstream::exportMatrix(adot_only,    "adot_coeffs",   "python", coeffDir);
    ITHACAstream::exportMatrix(velRBF_fluct, "a_adot_concat", "python", coeffDir);
    ITHACAstream::exportMatrix(velRBF_fluct, "a_adot_concat", "eigen",  coeffDir);

    // --- 8. Shape parameters (kept as before) ---
    const double eAvg   = 3;       // unused for AVG after linear-μ switch
    const double eFluct = 0.05;
    Eigen::MatrixXd radiiAvg, radiiFluct;

    if (ITHACAutilities::check_file("./radii_avg.txt")) {
        radiiAvg = ITHACAstream::readMatrix("./radii_avg.txt");
        M_Assert(radiiAvg.size() == nNutAvgModes, "radiiAvg size mismatch");
    } else {
        radiiAvg = Eigen::MatrixXd::Ones(nNutAvgModes, 1) * eAvg;
    }
    if (ITHACAutilities::check_file("./radii_fluct.txt")) {
        radiiFluct = ITHACAstream::readMatrix("./radii_fluct.txt");
        M_Assert(radiiFluct.size() == nNutFluctModes, "radiiFluct size mismatch");
    } else {
        radiiFluct = Eigen::MatrixXd::Ones(nNutFluctModes, 1) * eFluct;
    }

    // --- 9. AVG νt uses linear μ: persist μ-grid & coeff table (no RBF) ---
    List<SPLINTER::DataTable*> samplesNutAvg;  samplesNutAvg.resize(0);
    List<SPLINTER::RBFSpline*> rbfSplinesNutAvg; rbfSplinesNutAvg.resize(0);
    ITHACAstream::exportMatrix(muVecUnique, "NutAvg_mu_unique",    "eigen", coeffDir);
    ITHACAstream::exportMatrix(coeffNutAvg, "NutAvg_coeffs_by_mu", "eigen", coeffDir);

    // --- 10. FLUCT νt: build/load RBFs; SAVE/READ from ./ITHACAoutput/weightsSUP/ ---
    List<SPLINTER::DataTable*> samplesNutFluct;
    List<SPLINTER::RBFSpline*> rbfSplinesNutFluct;
    samplesNutFluct.resize(nNutFluctModes);
    rbfSplinesNutFluct.resize(nNutFluctModes);

    std::cout << ">>> [SUP] Building nut_fluct RBF splines...\n";
    for (label i = 0; i < nNutFluctModes; ++i) {
        word weightName = "wRBF_NUTFLUCT_" + name(i + 1);

        samplesNutFluct[i] = new SPLINTER::DataTable(velRBF_fluct.cols(), 1);
        for (label j = 0; j < velRBF_fluct.rows(); ++j)
            samplesNutFluct[i]->addSample(velRBF_fluct.row(j), coeffs_fluct(j, i));

        Eigen::MatrixXd weights;
        if (ITHACAutilities::check_file("./ITHACAoutput/weightsSUP/" + weightName)) {
            ITHACAstream::ReadDenseMatrix(weights, "./ITHACAoutput/weightsSUP/", weightName);
            rbfSplinesNutFluct[i] = new SPLINTER::RBFSpline(*samplesNutFluct[i],
                                        SPLINTER::RadialBasisFunctionType::GAUSSIAN,
                                        weights, radiiFluct(i));
            std::cout << "   [SUP] nut_fluct RBF " << i + 1 << "/" << nNutFluctModes
                      << " loaded from ./ITHACAoutput/weightsSUP/" << weightName << "\n";
        } else {
            rbfSplinesNutFluct[i] = new SPLINTER::RBFSpline(*samplesNutFluct[i],
                                        SPLINTER::RadialBasisFunctionType::GAUSSIAN,
                                        false, radiiFluct(i));
            ITHACAstream::SaveDenseMatrix(rbfSplinesNutFluct[i]->weights,
                                          "./ITHACAoutput/weightsSUP/", weightName);
            std::cout << "   [SUP] nut_fluct RBF " << i + 1 << "/" << nNutFluctModes
                      << " fitted & saved to ./ITHACAoutput/weightsSUP/" << weightName << "\n";
        }
    }

    this->rbfSplinesNutAvg   = rbfSplinesNutAvg;    // empty (avg via linear μ)
    this->rbfSplinesNutFluct = rbfSplinesNutFluct;
    this->samplesNutAvg      = samplesNutAvg;       // empty
    this->samplesNutFluct    = samplesNutFluct;

    std::cout << ">>> [SUP] Finished AVG linear-μ table export + FLUCT RBF build.\n";
}

}
void UnsteadyNSTurb::projectPPE(fileName folder, label NU, label NP, label NSUP,
                                label Nnut, bool rbfInterp)
{
    NUmodes   = NU;
    NPmodes   = NP;
    NSUPmodes = 0;
    nNutModes = Nnut;
    L_U_SUPmodes.resize(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            L_U_SUPmodes.append(tmp<volVectorField>(liftfield[k]));
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            L_U_SUPmodes.append(tmp<volVectorField>(Umodes[k]));
        }
    }

    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {
        word B_str = "B_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes);
        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + B_str))
            ITHACAstream::ReadDenseMatrix(B_matrix, "./ITHACAoutput/Matrices/", B_str);
        else
            B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);

        word btStr = "bt_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes);
        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + btStr))
            ITHACAstream::ReadDenseMatrix(btMatrix, "./ITHACAoutput/Matrices/", btStr);
        else
            btMatrix = btTurbulence(NUmodes, NSUPmodes);

        word K_str = "K_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPmodes);
        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + K_str))
            ITHACAstream::ReadDenseMatrix(K_matrix, "./ITHACAoutput/Matrices/", K_str);
        else
            K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);

        word M_str = "M_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes);
        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + M_str))
            ITHACAstream::ReadDenseMatrix(M_matrix, "./ITHACAoutput/Matrices/", M_str);
        else
            M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);

        word D_str = "D_" + name(NPmodes);
        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + D_str))
            ITHACAstream::ReadDenseMatrix(D_matrix, "./ITHACAoutput/Matrices/", D_str);
        else
            D_matrix = laplacian_pressure(NPmodes);

        word bc1_str = "BC1_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NPmodes);
        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + bc1_str))
            ITHACAstream::ReadDenseMatrix(BC1_matrix, "./ITHACAoutput/Matrices/", bc1_str);
        else
            BC1_matrix = pressure_BC1(NUmodes, NPmodes);

        word bc2_str = "BC2_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPmodes) + "_t";
        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + bc2_str))
            ITHACAstream::ReadDenseTensor(bc2Tensor, "./ITHACAoutput/Matrices/", bc2_str);
        else
            bc2Tensor = pressureBC2(NUmodes, NPmodes);

        word bc3_str = "BC3_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NPmodes);
        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + bc3_str))
            ITHACAstream::ReadDenseMatrix(BC3_matrix, "./ITHACAoutput/Matrices/", bc3_str);
        else
            BC3_matrix = pressure_BC3(NUmodes, NPmodes);

        word C_str = "C_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_t";
        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + C_str))
            ITHACAstream::ReadDenseTensor(C_tensor, "./ITHACAoutput/Matrices/", C_str);
        else
            C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);

        word ct1Str = "ct1_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(nNutModes) + "_t";
        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct1Str))
            ITHACAstream::ReadDenseTensor(ct1Tensor, "./ITHACAoutput/Matrices/", ct1Str);
        else
            ct1Tensor = turbulenceTensor1(NUmodes, NSUPmodes, nNutModes);

        word ct2Str = "ct2_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(nNutModes) + "_t";
        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct2Str))
            ITHACAstream::ReadDenseTensor(ct2Tensor, "./ITHACAoutput/Matrices/", ct2Str);
        else
            ct2Tensor = turbulenceTensor2(NUmodes, NSUPmodes, nNutModes);

        word G_str = "G_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPmodes) + "_t";
        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + G_str))
            ITHACAstream::ReadDenseTensor(gTensor, "./ITHACAoutput/Matrices/", G_str);
        else
            gTensor = divMomentum(NUmodes, NPmodes);

        // ==== Average tensors ====
        if (nutAve.size() != 0)
        {
            word ct1AveStr = "ct1Ave_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_t";
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct1AveStr))
                ITHACAstream::ReadDenseTensor(ct1AveTensor, "./ITHACAoutput/Matrices/", ct1AveStr);
            else
                ct1AveTensor = turbulenceAveTensor1(NUmodes, NSUPmodes);

            word ct2AveStr = "ct2Ave_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_t";
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct2AveStr))
                ITHACAstream::ReadDenseTensor(ct2AveTensor, "./ITHACAoutput/Matrices/", ct2AveStr);
            else
                ct2AveTensor = turbulenceAveTensor2(NUmodes, NSUPmodes);
        }

        // ==== Fluctuation tensors ====
        if (nutFluctModes.size() != 0)
        {
            word ct1FluctStr = "ct1Fluct_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_t";
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct1FluctStr))
                ITHACAstream::ReadDenseTensor(ct1FluctTensor, "./ITHACAoutput/Matrices/", ct1FluctStr);
            else
                ct1FluctTensor = turbulenceFluctTensor1(NUmodes, NSUPmodes);

            word ct2FluctStr = "ct2Fluct_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_t";
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct2FluctStr))
                ITHACAstream::ReadDenseTensor(ct2FluctTensor, "./ITHACAoutput/Matrices/", ct2FluctStr);
            else
                ct2FluctTensor = turbulenceFluctTensor2(NUmodes, NSUPmodes);
        }

        // ===== Full PPE tensors =====
        word ct1PPEStr = "ct1PPE_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPmodes) + "_" + name(nNutModes) + "_t";
        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct1PPEStr))
            ITHACAstream::ReadDenseTensor(ct1PPETensor, "./ITHACAoutput/Matrices/", ct1PPEStr);
        else
            ct1PPETensor = turbulencePPETensor1(NUmodes, NSUPmodes, NPmodes, nNutModes);

        word ct2PPEStr = "ct2PPE_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPmodes) + "_" + name(nNutModes) + "_t";
        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct2PPEStr))
            ITHACAstream::ReadDenseTensor(ct2PPETensor, "./ITHACAoutput/Matrices/", ct2PPEStr);
        else
            ct2PPETensor = turbulencePPETensor2(NUmodes, NSUPmodes, NPmodes, nNutModes);

        // ==== PPE Average tensors ====
        if (nutAve.size() != 0)
        {
            word ct1PPEAveStr = "ct1PPEAve_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPmodes) + "_t";
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct1PPEAveStr))
                ITHACAstream::ReadDenseTensor(ct1PPEAveTensor, "./ITHACAoutput/Matrices/", ct1PPEAveStr);
            else
                ct1PPEAveTensor = turbulencePPEAveTensor1(NUmodes, NSUPmodes, NPmodes);

            word ct2PPEAveStr = "ct2PPEAve_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPmodes) + "_t";
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct2PPEAveStr))
                ITHACAstream::ReadDenseTensor(ct2PPEAveTensor, "./ITHACAoutput/Matrices/", ct2PPEAveStr);
            else
                ct2PPEAveTensor = turbulencePPEAveTensor2(NUmodes, NSUPmodes, NPmodes);
        }

        // ==== PPE Fluctuation tensors ====
        if (nutFluctModes.size() != 0)
        {
            word ct1PPEFluctStr = "ct1PPEFluct_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPmodes) + "_t";
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct1PPEFluctStr))
                ITHACAstream::ReadDenseTensor(ct1PPEFluctTensor, "./ITHACAoutput/Matrices/", ct1PPEFluctStr);
            else
                ct1PPEFluctTensor = turbulencePPEFluctTensor1(NUmodes, NSUPmodes, NPmodes);

            word ct2PPEFluctStr = "ct2PPEFluct_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPmodes) + "_t";
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct2PPEFluctStr))
                ITHACAstream::ReadDenseTensor(ct2PPEFluctTensor, "./ITHACAoutput/Matrices/", ct2PPEFluctStr);
            else
                ct2PPEFluctTensor = turbulencePPEFluctTensor2(NUmodes, NSUPmodes, NPmodes);
        }

        // ===== PPE L vector (load or build via pressurePPE_L) =====
        if (NPmodes > 0)
        {
            word L_str = "L_" + name(NPmodes);
            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + L_str))
            {
                ITHACAstream::ReadDenseMatrix(L_vector, "./ITHACAoutput/Matrices/", L_str);
            }
            else
            {
                L_vector = pressurePPE_L(NPmodes); // builds & saves L internally
            }
        }

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }
    }
    else
    {
        // If no matrices folder, just build everything
        B_matrix   = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        C_tensor   = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        M_matrix   = mass_term(NUmodes, NPmodes, NSUPmodes);
        K_matrix   = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        D_matrix   = laplacian_pressure(NPmodes);
        gTensor    = divMomentum(NUmodes, NPmodes);
        BC1_matrix = pressure_BC1(NUmodes, NPmodes);
        bc2Tensor  = pressureBC2(NUmodes, NPmodes);
        BC3_matrix = pressure_BC3(NUmodes, NPmodes);
        btMatrix   = btTurbulence(NUmodes, NSUPmodes);
        ct1Tensor  = turbulenceTensor1(NUmodes, NSUPmodes, nNutModes);
        ct2Tensor  = turbulenceTensor2(NUmodes, NSUPmodes, nNutModes);
        ct1PPETensor = turbulencePPETensor1(NUmodes, NSUPmodes, NPmodes, nNutModes);
        ct2PPETensor = turbulencePPETensor2(NUmodes, NSUPmodes, NPmodes, nNutModes);

        if (nutAve.size() != 0)
        {
            ct1AveTensor    = turbulenceAveTensor1(NUmodes, NSUPmodes);
            ct2AveTensor    = turbulenceAveTensor2(NUmodes, NSUPmodes);
            ct1PPEAveTensor = turbulencePPEAveTensor1(NUmodes, NSUPmodes, NPmodes);
            ct2PPEAveTensor = turbulencePPEAveTensor2(NUmodes, NSUPmodes, NPmodes);
        }
        if (nutFluctModes.size() != 0)
        {
            ct1FluctTensor    = turbulenceFluctTensor1(NUmodes, NSUPmodes);
            ct2FluctTensor    = turbulenceFluctTensor2(NUmodes, NSUPmodes);
            ct1PPEFluctTensor = turbulencePPEFluctTensor1(NUmodes, NSUPmodes, NPmodes);
            ct2PPEFluctTensor = turbulencePPEFluctTensor2(NUmodes, NSUPmodes, NPmodes);
        }

        if (NPmodes > 0)
        {
            L_vector = pressurePPE_L(NPmodes); // builds and saves
        }

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }
    }

    // ------------ Export block (Python, if enabled) ------------
    if (para->exportPython)
    {
        ITHACAstream::exportMatrix(B_matrix,   "B",   "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix,   "K",   "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(D_matrix,   "D",   "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix,   "M",   "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC1_matrix, "BC1", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC3_matrix, "BC3", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor,   "C",   "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(gTensor,    "G",   "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(bc2Tensor,  "BC2", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct1Tensor,  "ct1", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct2Tensor,  "ct2", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct1PPETensor, "ct1PPE", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct2PPETensor, "ct2PPE", "python", "./ITHACAoutput/Matrices/");
        if (nutAve.size() != 0)
        {
            ITHACAstream::exportTensor(ct1AveTensor,    "ct1Ave",    "python", "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2AveTensor,    "ct2Ave",    "python", "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct1PPEAveTensor, "ct1PPEAve", "python", "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2PPEAveTensor, "ct2PPEAve", "python", "./ITHACAoutput/Matrices/");
        }
        if (nutFluctModes.size() != 0)
        {
            ITHACAstream::exportTensor(ct1FluctTensor,    "ct1Fluct",    "python", "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2FluctTensor,    "ct2Fluct",    "python", "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct1PPEFluctTensor, "ct1PPEFluct", "python", "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2PPEFluctTensor, "ct2PPEFluct", "python", "./ITHACAoutput/Matrices/");
        }
        if (NPmodes > 0)
        {
            ITHACAstream::exportMatrix(L_vector, "L", "python", "./ITHACAoutput/Matrices/");
        }
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportMatrix(B_matrix,   "B",   "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix,   "K",   "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(D_matrix,   "D",   "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix,   "M",   "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC1_matrix, "BC1", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC3_matrix, "BC3", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor,   "C",   "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(gTensor,    "G",   "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(bc2Tensor,  "BC2", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct1Tensor,  "ct1", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct2Tensor,  "ct2", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct1PPETensor, "ct1PPE", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct2PPETensor, "ct2PPE", "matlab", "./ITHACAoutput/Matrices/");
        if (nutAve.size() != 0)
        {
            ITHACAstream::exportTensor(ct1AveTensor,    "ct1Ave",    "matlab", "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2AveTensor,    "ct2Ave",    "matlab", "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct1PPEAveTensor, "ct1PPEAve", "matlab", "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2PPEAveTensor, "ct2PPEAve", "matlab", "./ITHACAoutput/Matrices/");
        }
        if (nutFluctModes.size() != 0)
        {
            ITHACAstream::exportTensor(ct1FluctTensor,    "ct1Fluct",    "matlab", "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2FluctTensor,    "ct2Fluct",    "matlab", "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct1PPEFluctTensor, "ct1PPEFluct", "matlab", "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2PPEFluctTensor, "ct2PPEFluct", "matlab", "./ITHACAoutput/Matrices/");
        }
        if (NPmodes > 0)
        {
            ITHACAstream::exportMatrix(L_vector, "L", "matlab", "./ITHACAoutput/Matrices/");
        }
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportMatrix(B_matrix,   "B",    "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix,   "K",    "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(D_matrix,   "D",    "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix,   "M",    "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC1_matrix, "BC1",  "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC3_matrix, "BC3",  "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor,   "C",    "eigen", "./ITHACAoutput/Matrices/C");
        ITHACAstream::exportTensor(gTensor,    "G",    "eigen", "./ITHACAoutput/Matrices/G");
        ITHACAstream::exportTensor(bc2Tensor,  "BC2_", "eigen", "./ITHACAoutput/Matrices/BC2");
        ITHACAstream::exportMatrix(btMatrix,   "bt",   "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct1Tensor,  "ct1_", "eigen", "./ITHACAoutput/Matrices/ct1");
        ITHACAstream::exportTensor(ct2Tensor,  "ct2_", "eigen", "./ITHACAoutput/Matrices/ct2");
        ITHACAstream::exportTensor(ct1PPETensor, "ct1PPE_", "eigen", "./ITHACAoutput/Matrices/ct1PPE");
        ITHACAstream::exportTensor(ct2PPETensor, "ct2PPE_", "eigen", "./ITHACAoutput/Matrices/ct2PPE");
        if (nutAve.size() != 0)
        {
            ITHACAstream::exportTensor(ct1AveTensor,    "ct1Ave_",    "eigen", "./ITHACAoutput/Matrices/ct1Ave");
            ITHACAstream::exportTensor(ct2AveTensor,    "ct2Ave_",    "eigen", "./ITHACAoutput/Matrices/ct2Ave");
            ITHACAstream::exportTensor(ct1PPEAveTensor, "ct1PPEAve_", "eigen", "./ITHACAoutput/Matrices/ct1PPEAve");
            ITHACAstream::exportTensor(ct2PPEAveTensor, "ct2PPEAve_", "eigen", "./ITHACAoutput/Matrices/ct2PPEAve");
        }
        if (nutFluctModes.size() != 0)
        {
            ITHACAstream::exportTensor(ct1FluctTensor,    "ct1Fluct",    "eigen", "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2FluctTensor,    "ct2Fluct",    "eigen", "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct1PPEFluctTensor, "ct1PPEFluct", "eigen", "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2PPEFluctTensor, "ct2PPEFluct", "eigen", "./ITHACAoutput/Matrices/");
        }
        if (NPmodes > 0)
        {
            ITHACAstream::exportMatrix(L_vector, "L", "eigen", "./ITHACAoutput/Matrices/");
        }
    }

    // --------- Compose total tensors ----------
    bTotalMatrix = B_matrix + btMatrix;
    label cSize = NUmodes + NSUPmodes + liftfield.size();

    cTotalTensor.resize(cSize, nNutModes, cSize);
    cTotalTensor = ct1Tensor + ct2Tensor;

    cTotalPPETensor.resize(NPmodes, nNutModes, cSize);
    cTotalPPETensor = ct1PPETensor + ct2PPETensor;

    // Average (if available)
    if (nutAve.size() != 0)
    {
        cTotalAveTensor.resize(cSize, nutAve.size(), cSize);
        cTotalAveTensor = ct1AveTensor + ct2AveTensor;

        cTotalPPEAveTensor.resize(NPmodes, nutAve.size(), cSize);
        cTotalPPEAveTensor = ct1PPEAveTensor + ct2PPEAveTensor;
    }

    // Fluctuation (if available)
    if (nutFluctModes.size() != 0)
    {
        cTotalFluctTensor.resize(cSize, nutFluctModes.size(), cSize);
        cTotalFluctTensor = ct1FluctTensor + ct2FluctTensor;

        cTotalPPEFluctTensor.resize(NPmodes, nutFluctModes.size(), cSize);
        cTotalPPEFluctTensor = ct1PPEFluctTensor + ct2PPEFluctTensor;
    }



if (rbfInterp && (!Pstream::parRun()))
{
    std::cout << "\n==== [RBF DEBUG] ENTERING OFFLINE RBF CONSTRUCTION (AVG linear-μ + FLUCT RBF) ====\n";

    const word coeffDir = "./ITHACAoutput/Coefficients/";
    const word muPath   = "./ITHACAoutput/Offline/mu_samples_mat.txt";

    // --- 1. LOAD mu and time vectors ---
    Eigen::MatrixXd muMat = ITHACAstream::readMatrix(muPath);
    std::cout << "[RBF DEBUG] muMat shape: " << muMat.rows() << " x " << muMat.cols() << std::endl;
    const int nPar = 12, nSnapshotsPerPar = 200;
    Eigen::VectorXd timeVec = muMat.col(0), muVec = muMat.col(1);

    // --- 2. LOAD avg/fluct nut coeffs ---
    Eigen::MatrixXd coeffNutAvg   = ITHACAstream::readMatrix(coeffDir + "Nut_avg_coeffs_mat.txt");
    Eigen::MatrixXd coeffNutFluct = ITHACAstream::readMatrix(coeffDir + "Nut_fluct_coeffs_mat.txt");
    std::cout << "[RBF DEBUG] Loaded coeffNutAvg shape: "   << coeffNutAvg.rows()   << " x " << coeffNutAvg.cols()   << std::endl;
    std::cout << "[RBF DEBUG] Loaded coeffNutFluct shape: " << coeffNutFluct.rows() << " x " << coeffNutFluct.cols() << std::endl;
    int nNutAvgModes   = coeffNutAvg.rows();
    int nNutFluctModes = coeffNutFluct.rows();
    int nUniqueMu      = nPar;

    // --- 3. Compute velocity projection coefficients ---
    Eigen::MatrixXd a = ITHACAutilities::getCoeffs(Uomfield, Umodes);
    a.transposeInPlace(); // [nSnapshots, nVelModes]

    // --- 4. Prepare snapshot grouping for each parameter ---
    Eigen::VectorXd initSnapInd(nPar), timeSnap(nPar);
    for (int i = 0; i < nPar; ++i) {
        int start = i * nSnapshotsPerPar;
        initSnapInd(i) = start;
        timeSnap(i) = timeVec(start + 1) - timeVec(start);
        std::cout << "[RBF DEBUG] i=" << i << ", start=" << start << ", dt=" << timeSnap(i) << std::endl;
    }

    // --- 5. Prepare unique muVec for average part ---
    Eigen::VectorXd muVecUnique(nUniqueMu);
    for (int i = 0; i < nUniqueMu; ++i)
        muVecUnique(i) = muVec(static_cast<int>(initSnapInd(i)));
    std::cout << "[RBF DEBUG] muVecUnique: [" << muVecUnique(0) << " ... " << muVecUnique(muVecUnique.size()-1) << "] (M=" << muVecUnique.size() << ")\n";

    // --- 6. Prepare fluctuation interpolation input: [a, aDot] ---
    std::cout << "[RBF DEBUG] Calling velDerivativeCoeff() for fluctuation part...\n";
    Eigen::MatrixXd Gfluct = coeffNutFluct.transpose();
    List<Eigen::MatrixXd> interpDataFluct = velDerivativeCoeff(a, Gfluct, initSnapInd, timeSnap);
    Eigen::MatrixXd velRBF_fluct = interpDataFluct[0];  // [a, aDot]
    Eigen::MatrixXd coeffs_fluct = interpDataFluct[1];
    std::cout << "[RBF DEBUG] velRBF_fluct shape: " << velRBF_fluct.rows() << " x " << velRBF_fluct.cols() << std::endl;
    std::cout << "[RBF DEBUG] coeffs_fluct shape: " << coeffs_fluct.rows() << " x " << coeffs_fluct.cols() << std::endl;

    // --- 7. Extract adot_only and export adot/a_adot_concat ---
    int nSnapshots = velRBF_fluct.rows();
    int nModes = velRBF_fluct.cols() / 2;
    Eigen::MatrixXd adot_only(nSnapshots, nModes);
    for (int i = 0; i < nSnapshots; ++i)
        adot_only.row(i) = velRBF_fluct.row(i).tail(nModes);
    ITHACAstream::exportMatrix(adot_only,    "adot_coeffs",   "python", coeffDir);
    ITHACAstream::exportMatrix(velRBF_fluct, "a_adot_concat", "python", coeffDir);
    ITHACAstream::exportMatrix(velRBF_fluct, "a_adot_concat", "eigen",  coeffDir);

    // --- 8. DEBUG PRINTS FOR FIRST ROWS ---
    std::cout << "\n[RBF DEBUG] First row of a (velocity coeffs): ";
    for (int k = 0; k < a.cols(); ++k) std::cout << a(0, k) << " ";
    std::cout << std::endl;
    std::cout << "[RBF DEBUG] First row of velRBF_fluct ([a, aDot]): ";
    for (int k = 0; k < velRBF_fluct.cols(); ++k) std::cout << velRBF_fluct(0, k) << " ";
    std::cout << std::endl;
    std::cout << "[RBF DEBUG] First row of adot_only (aDot): ";
    for (int k = 0; k < adot_only.cols(); ++k) std::cout << adot_only(0, k) << " ";
    std::cout << std::endl;

    // --- 9. SHAPE PARAMETERS (only used for FLUCT below; AVG will be linear μ) ---
    const double eAvg   = 3;      // kept for completeness; not used for AVG now
    const double eFluct = 0.05;

    Eigen::MatrixXd radiiAvg, radiiFluct;

    // AVG nut (not used after switch to linear μ; keep for compatibility)
    if (ITHACAutilities::check_file("./radii_avg.txt")) {
        radiiAvg = ITHACAstream::readMatrix("./radii_avg.txt");
        std::cout << "[RBF DEBUG] Loaded radii_avg.txt, shape: " << radiiAvg.rows() << " x " << radiiAvg.cols() << std::endl;
        M_Assert(radiiAvg.size() == nNutAvgModes, "radiiAvg size mismatch");
    } else {
        radiiAvg = Eigen::MatrixXd::Ones(nNutAvgModes, 1) * eAvg;
        std::cout << "[RBF DEBUG] Set default radii for nutAvg (unused with linear μ), e=" << eAvg << std::endl;
    }

    // FLUCT nut
    if (ITHACAutilities::check_file("./radii_fluct.txt")) {
        radiiFluct = ITHACAstream::readMatrix("./radii_fluct.txt");
        std::cout << "[RBF DEBUG] Loaded radii_fluct.txt, shape: " << radiiFluct.rows() << " x " << radiiFluct.cols() << std::endl;
        M_Assert(radiiFluct.size() == nNutFluctModes, "radiiFluct size mismatch");
    } else {
        radiiFluct = Eigen::MatrixXd::Ones(nNutFluctModes, 1) * eFluct;
        std::cout << "[RBF DEBUG] Set default radii for nutFluct, e=" << eFluct << std::endl;
    }

    // --- 10. AVG νt: persist μ-grid & coeff table for linear interpolation (NO RBF) ---
    List<SPLINTER::DataTable*> samplesNutAvg;  samplesNutAvg.resize(0);
    List<SPLINTER::RBFSpline*> rbfSplinesNutAvg; rbfSplinesNutAvg.resize(0);

    // Write exactly the filenames you’ll read online later
    ITHACAstream::exportMatrix(muVecUnique, "NutAvg_mu_unique",        "eigen", coeffDir);
    ITHACAstream::exportMatrix(coeffNutAvg, "NutAvg_coeffs_by_mu",     "eigen", coeffDir);

    std::cout << ">>> Persisted νt_avg tables for linear μ-interpolation: "
              << "mu size=" << muVecUnique.size()
              << ", coeff table=" << coeffNutAvg.rows() << "x" << coeffNutAvg.cols() << "\n";

    // --- 11. BUILD FLUCT RBFs ([a, aDot] input) --- (UNCHANGED)
    List<SPLINTER::DataTable*> samplesNutFluct;
    List<SPLINTER::RBFSpline*> rbfSplinesNutFluct;
    samplesNutFluct.resize(nNutFluctModes);
    rbfSplinesNutFluct.resize(nNutFluctModes);

    std::cout << ">>> Building nut_fluct RBF splines...\n";
    for (label i = 0; i < nNutFluctModes; ++i) {
        word weightName = "wRBF_NUTFLUCT_" + name(i + 1);
        samplesNutFluct[i] = new SPLINTER::DataTable(velRBF_fluct.cols(), 1);
        for (label j = 0; j < velRBF_fluct.rows(); ++j)
            samplesNutFluct[i]->addSample(velRBF_fluct.row(j), coeffs_fluct(j, i));
        Eigen::MatrixXd weights;
        if (ITHACAutilities::check_file("./ITHACAoutput/weightsPPE/" + weightName)) {
            ITHACAstream::ReadDenseMatrix(weights, "./ITHACAoutput/weightsPPE/", weightName);
            rbfSplinesNutFluct[i] = new SPLINTER::RBFSpline(*samplesNutFluct[i],
                                            SPLINTER::RadialBasisFunctionType::GAUSSIAN,
                                            weights, radiiFluct(i));
            std::cout << "   nut_fluct RBF " << i + 1 << "/" << nNutFluctModes << " loaded from weights.\n";
        } else {
            rbfSplinesNutFluct[i] = new SPLINTER::RBFSpline(*samplesNutFluct[i],
                                            SPLINTER::RadialBasisFunctionType::GAUSSIAN,
                                            false, radiiFluct(i));
            ITHACAstream::SaveDenseMatrix(rbfSplinesNutFluct[i]->weights,
                                          "./ITHACAoutput/weightsPPE/", weightName);
            std::cout << "   nut_fluct RBF " << i + 1 << "/" << nNutFluctModes << " fitted & saved.\n";
        }
    }

    std::cout << "[OFFLINE] Built nutAvgSplines:  "  << rbfSplinesNutAvg.size()   << " (expected 0 for linear μ)\n";
    std::cout << "[OFFLINE] Built nutFluctSplines:"  << rbfSplinesNutFluct.size() << std::endl;
    std::cout << ">>> Finished AVG linear-μ table export + FLUCT RBF build.\n";

    this->rbfSplinesNutAvg   = rbfSplinesNutAvg;    // empty (avg uses linear μ)
    this->rbfSplinesNutFluct = rbfSplinesNutFluct;
    this->samplesNutAvg      = samplesNutAvg;       // empty
    this->samplesNutFluct    = samplesNutFluct;

    std::cout << "[OFFLINE] Built nutAvgSplines: "   << rbfSplinesNutAvg.size()   << std::endl;
    std::cout << "[OFFLINE] Built nutFluctSplines: " << rbfSplinesNutFluct.size() << std::endl;
}
}

List<Eigen::MatrixXd> UnsteadyNSTurb::velDerivativeCoeff(
    Eigen::MatrixXd A,
    Eigen::MatrixXd G,
    Eigen::VectorXd initSnapInd, // [start index per parameter]
    Eigen::VectorXd timeSnap     // [dt per parameter]
)
{
    List<Eigen::MatrixXd> newCoeffs;
    newCoeffs.setSize(2);

    label velCoeffsNum = A.cols();
    label snapshotsNum = A.rows();
    label parsSamplesNum = initSnapInd.size();
    label timeSnapshotsPerSample = snapshotsNum / parsSamplesNum;
    label newColsNum = 2 * velCoeffsNum;
    label newRowsNum = snapshotsNum - parsSamplesNum; // Each param loses 1 sample for derivative

    newCoeffs[0].resize(newRowsNum, newColsNum); // [a, adot]
    newCoeffs[1].resize(newRowsNum, G.cols());   // target coeffs

    int rowCount = 0;

    for (label j = 0; j < parsSamplesNum; ++j)
    {
        int i0 = j * timeSnapshotsPerSample; // starting index for this μ
        int N  = timeSnapshotsPerSample;

        // Backward difference within the block: n = 1..N-1
        for (int n = 1; n < N; ++n, ++rowCount)
        {
            Eigen::RowVectorXd a_now  = A.row(i0 + n);      // a(n)
            Eigen::RowVectorXd a_prev = A.row(i0 + n - 1);  // a(n-1)
            Eigen::RowVectorXd adot   = (a_now - a_prev) / timeSnap(j);

            // Build [a(n), ȧ(n)] and target coeffs at n
            newCoeffs[0].row(rowCount) << a_now, adot;
            newCoeffs[1].row(rowCount)  = G.row(i0 + n);
        }
    }

    interChoice = 3; // your flag
    return newCoeffs;
}

List < Eigen::MatrixXd > UnsteadyNSTurb::velParCoeff(Eigen::MatrixXd A,
        Eigen::MatrixXd G)
{
    List < Eigen::MatrixXd > newCoeffs;
    newCoeffs.setSize(2);
    Eigen::MatrixXd pars;
    pars = z.leftCols(z.cols() - 1);
    newCoeffs[0].resize(A.rows(), A.cols() + z.cols() - 1);
    newCoeffs[1].resize(G.rows(), G.cols());
    newCoeffs[0] << pars, A;
    newCoeffs[1] = G;
    interChoice = 2;
    return newCoeffs;
}
List<Eigen::MatrixXd> UnsteadyNSTurb::velParDerivativeCoeff(
    Eigen::MatrixXd A, Eigen::MatrixXd G,
    Eigen::VectorXd initSnapInd, Eigen::VectorXd timeSnap)
{
    std::cout << "[velParDerivativeCoeff] ENTER" << std::endl;
    List<Eigen::MatrixXd> newCoeffs;
    newCoeffs.setSize(2);

    label velCoeffsNum = A.cols(); // number of velocity modes
    label snapshotsNum = A.rows(); // total number of snapshots
    label parsSamplesNum = initSnapInd.size(); // number of parameter samples
    label timeSnapshotsPerSample = snapshotsNum / parsSamplesNum; // =200 in your case

    std::cout << "[velParDerivativeCoeff] A shape: " << snapshotsNum << " x " << velCoeffsNum << std::endl;
    std::cout << "[velParDerivativeCoeff] G shape: " << G.rows() << " x " << G.cols() << std::endl;
    std::cout << "[velParDerivativeCoeff] nPars: " << parsSamplesNum << ", timeSnapshotsPerSample: " << timeSnapshotsPerSample << std::endl;

    // Construct parameter matrix (mu for each snapshot)
    Eigen::MatrixXd pars(snapshotsNum, 1);
    for (label j = 0; j < parsSamplesNum; ++j)
        pars.block(j * timeSnapshotsPerSample, 0, timeSnapshotsPerSample, 1) = 
            Eigen::VectorXd::Constant(timeSnapshotsPerSample, mu(j)); // mu(j) should be known

    label newColsNum = 2 * velCoeffsNum;
    label newRowsNum = snapshotsNum - parsSamplesNum; // each group loses 1 row due to diff

    std::cout << "[velParDerivativeCoeff] newCoeffs[0] size: " << newRowsNum << " x " << (newColsNum + pars.cols()) << std::endl;
    std::cout << "[velParDerivativeCoeff] newCoeffs[1] size: " << newRowsNum << " x " << G.cols() << std::endl;

    newCoeffs[0].resize(newRowsNum, newColsNum + pars.cols());
    newCoeffs[1].resize(newRowsNum, G.cols());

    int totalRowsWritten = 0;

    for (label j = 0; j < parsSamplesNum; j++)
    {
        int start = j * timeSnapshotsPerSample;
        int rowOffset = j * (timeSnapshotsPerSample - 1);

        std::cout << "[velParDerivativeCoeff] Group " << j << ": start=" << start << ", rowOffset=" << rowOffset << std::endl;

        // b0: a(t), size (timeSnapshotsPerSample-1, velCoeffsNum)
        std::cout << "[velParDerivativeCoeff] b0: rows " << start << " to " << (start + timeSnapshotsPerSample - 2) << std::endl;
        if (start + timeSnapshotsPerSample - 1 > snapshotsNum) {
            std::cerr << "[velParDerivativeCoeff][ERROR] b0 access out of bounds! (start=" << start << ", timeSnapshotsPerSample-1=" << (timeSnapshotsPerSample-1) << ", snapshotsNum=" << snapshotsNum << ")" << std::endl;
            abort();
        }
        Eigen::MatrixXd b0 = A.middleRows(start, timeSnapshotsPerSample - 1);

        // b1: a(t+dt)
        std::cout << "[velParDerivativeCoeff] b1: rows " << (start+1) << " to " << (start+timeSnapshotsPerSample-1) << std::endl;
        if (start+1 + timeSnapshotsPerSample - 2 >= snapshotsNum) {
            std::cerr << "[velParDerivativeCoeff][ERROR] b1 access out of bounds! (start+1=" << (start+1) << ", timeSnapshotsPerSample-1=" << (timeSnapshotsPerSample-1) << ", snapshotsNum=" << snapshotsNum << ")" << std::endl;
            abort();
        }
        Eigen::MatrixXd b1 = A.middleRows(start+1, timeSnapshotsPerSample - 1);

        Eigen::MatrixXd bNew(b0.rows(), b0.cols() + b1.cols());
        bNew << b1, (b1 - b0) / timeSnap(j);

        // [mu, a, adot]
        std::cout << "[velParDerivativeCoeff] pars block for input: rows " << (start+1) << " to " << (start+timeSnapshotsPerSample-1) << std::endl;
        if (start+1 + timeSnapshotsPerSample - 2 >= snapshotsNum) {
            std::cerr << "[velParDerivativeCoeff][ERROR] pars access out of bounds!" << std::endl;
            abort();
        }
        newCoeffs[0].block(rowOffset, 0, timeSnapshotsPerSample - 1, pars.cols()) =
            pars.middleRows(start+1, timeSnapshotsPerSample - 1);

        // [a, adot]
        newCoeffs[0].block(rowOffset, pars.cols(), timeSnapshotsPerSample - 1, newColsNum) = bNew;

        // nut coeffs (target)
        std::cout << "[velParDerivativeCoeff] G block for output: rows " << (start+1) << " to " << (start+timeSnapshotsPerSample-1) << std::endl;
        if (start+1 + timeSnapshotsPerSample - 2 >= G.rows()) {
            std::cerr << "[velParDerivativeCoeff][ERROR] G access out of bounds!" << std::endl;
            abort();
        }
        newCoeffs[1].middleRows(rowOffset, timeSnapshotsPerSample - 1) =
            G.middleRows(start+1, timeSnapshotsPerSample - 1);

        totalRowsWritten += timeSnapshotsPerSample - 1;
        std::cout << "[velParDerivativeCoeff] Finished group " << j << ", totalRowsWritten so far: " << totalRowsWritten << std::endl;
    }

    std::cout << "[velParDerivativeCoeff] FINISHED. Total rows written: " << totalRowsWritten << "/" << newRowsNum << std::endl;
    interChoice = 4;
    return newCoeffs;
}


Eigen::MatrixXd UnsteadyNSTurb::velParDerivativeCoeff(Eigen::MatrixXd A,
        Eigen::VectorXd par, double timeSnap)
{
    Eigen::MatrixXd newCoeffs;
    label velCoeffsNum = A.cols();
    label snapshotsNum = A.rows();
    label parsSamplesNum = par.size();
    label newColsNum = 2 * velCoeffsNum + parsSamplesNum;
    label newRowsNum = snapshotsNum - 1;
    newCoeffs.resize(newRowsNum, newColsNum);
    Eigen::MatrixXd b0 = A.topRows(A.rows() - 1);
    Eigen::MatrixXd b1 = A.bottomRows(A.rows() - 1);
    Eigen::MatrixXd bNew(b0.rows(), b0.cols() + b1.cols());
    bNew << b1, ((b1 - b0) / (timeSnap));
    newCoeffs.leftCols(parsSamplesNum) = Eigen::MatrixXd::Ones(newRowsNum,
                                         parsSamplesNum) * par;
    newCoeffs.rightCols(newColsNum - parsSamplesNum) = bNew;
    return newCoeffs;
}
