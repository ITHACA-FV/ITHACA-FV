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

#include "UnsteadyNSTurbIntrusive.H"

/// \file
/// Source file of the unsteadyNS class.

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Construct Null
UnsteadyNSTurbIntrusive::UnsteadyNSTurbIntrusive() {}

// Construct from zero
UnsteadyNSTurbIntrusive::UnsteadyNSTurbIntrusive(int argc, char* argv[])
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
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void UnsteadyNSTurbIntrusive::truthSolve(List<scalar> mu_now)
{
    Time& runTime = _runTime();
    surfaceScalarField& phi = _phi();
    fvMesh& mesh = _mesh();
#include "initContinuityErrs.H"
    fv::options& fvOptions = _fvOptions();
    pimpleControl& pimple = _pimple();
    volScalarField p = _p();
    volVectorField U = _U();
    volScalarField nut = _nut();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    // Perform a TruthSolve
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;
    // Initialize Nsnapshots
    int nsnapshots = 0;

    // Start the time loop
    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime.setEndTime(finalTime + timeStep);
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
            volScalarField nut = turbulence->nut();
            ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(nut, name(counter), "./ITHACAoutput/Offline/");
            std::ofstream of("./ITHACAoutput/Offline/" + name(counter) + "/" +
                             runTime.timeName());
            Ufield.append(U);
            Pfield.append(p);
            nutFields.append(nut);
            counter++;
            nextWrite += writeEvery;
            writeMu(mu_now);
            // --- Fill in the mu_samples with parameters (time, mu) to be used for the PODI sample points
            mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size() + 1);
            mu_samples(mu_samples.rows() - 1, 0) = atof(runTime.timeName().c_str());

            for (int i = 0; i < mu_now.size(); i++)
            {
                mu_samples(mu_samples.rows() - 1, i + 1) = mu_now[i];
            }
        }

        runTime++;
    }

    // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
    if (mu.cols() == 0)
    {
        mu.resize(1, 1);
    }

    if (mu_samples.rows() == nsnapshots * mu.cols())
    {
        ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
                                   "./ITHACAoutput/Offline");
    }
}

Eigen::Tensor<double, 3> UnsteadyNSTurbIntrusive::turbulenceTensor1(
    label nModes)
{
    Eigen::Tensor<double, 3> ct1Tensor;
    ct1Tensor.resize(nModes, nModes, nModes);

    for (label i = 0; i < nModes; i++)
    {
        for (label j = 0; j < nModes; j++)
        {
            for (label k = 0; k < nModes; k++)
            {
                ct1Tensor(i, j, k) = fvc::domainIntegrate(Umodes[i] & fvc::laplacian(
                                         nutModes[j], Umodes[k])).value();
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(ct1Tensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(ct1Tensor, "./ITHACAoutput/Matrices/",
                                  "ct1_" + name(nModes) + "_t");
    return ct1Tensor;
}

Eigen::Tensor<double, 3> UnsteadyNSTurbIntrusive::turbulenceTensor2(
    label nModes)
{
    Eigen::Tensor<double, 3> ct2Tensor;
    ct2Tensor.resize(nModes, nModes, nModes);

    for (label i = 0; i < nModes; i++)
    {
        for (label j = 0; j < nModes; j++)
        {
            for (label k = 0; k < nModes; k++)
            {
                ct2Tensor(i, j, k) = fvc::domainIntegrate(Umodes[i] & (fvc::div(
                                         nutModes[j] * dev((fvc::grad(Umodes[k]))().T())))).value();
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(ct2Tensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(ct2Tensor, "./ITHACAoutput/Matrices/",
                                  "ct2_" + name(nModes) + "_t");
    return ct2Tensor;
}

Eigen::MatrixXd UnsteadyNSTurbIntrusive::btTurbulence(label nModes)
{
    Eigen::MatrixXd btMatrix(nModes, nModes);
    btMatrix = btMatrix * 0;

    // Project everything
    for (label i = 0; i < nModes; i++)
    {
        for (label j = 0; j < nModes; j++)
        {
            btMatrix(i, j) = fvc::domainIntegrate(Umodes[i] & (fvc::div(dev((T(
                    fvc::grad(
                        Umodes[j]))))))).value();
        }
    }

    // Export the matrix
    ITHACAstream::SaveDenseMatrix(btMatrix, "./ITHACAoutput/Matrices/",
                                  "bt_" + name(nModes) + "_t");
    return btMatrix;
}

void UnsteadyNSTurbIntrusive::project(fileName folder, label nModes)
{
    nModesOnline = nModes;
    L_U_SUPmodes.resize(0);

    if (nModes != 0)
    {
        for (label k = 0; k < nModes; k++)
        {
            L_U_SUPmodes.append(Umodes[k]);
        }
    }

    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {
        word bStr = "b_" + name(nModes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + bStr))
        {
            ITHACAstream::ReadDenseMatrix(bMatrix, "./ITHACAoutput/Matrices/", bStr);
        }
        else
        {
            bMatrix = diffusiveTerm(nModes);
        }

        word btStr = "bt_" + name(nModes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + btStr))
        {
            ITHACAstream::ReadDenseMatrix(btMatrix, "./ITHACAoutput/Matrices/", btStr);
        }
        else
        {
            btMatrix = btTurbulence(nModes);
        }

        word kStr = "k_" + name(nModes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + kStr))
        {
            ITHACAstream::ReadDenseMatrix(kMatrix, "./ITHACAoutput/Matrices/", kStr);
        }
        else
        {
            kMatrix = pressureGradientTerm(nModes);
        }

        word cStr = "c_" + name(nModes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + cStr))
        {
            ITHACAstream::ReadDenseTensor(convTensor, "./ITHACAoutput/Matrices/", cStr);
        }
        else
        {
            convTensor = convectiveTerm(nModes);
        }

        word ct1Str = "ct1_" + name(nModes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct1Str))
        {
            ITHACAstream::ReadDenseTensor(ct1Tensor, "./ITHACAoutput/Matrices/", ct1Str);
        }
        else
        {
            ct1Tensor = turbulenceTensor1(nModes);
        }

        word ct2Str = "ct2_" + name(nModes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct2Str))
        {
            ITHACAstream::ReadDenseTensor(ct2Tensor, "./ITHACAoutput/Matrices/", ct2Str);
        }
        else
        {
            ct2Tensor = turbulenceTensor2(nModes);
        }

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(nModes);
            bcVelMat = bcVelocityMat(nModes);
        }
    }
    else
    {
        bMatrix = diffusiveTerm(nModes);
        convTensor = convectiveTerm(nModes);
        kMatrix = pressureGradientTerm(nModes);
        btMatrix = btTurbulence(nModes);
        ct1Tensor = turbulenceTensor1(nModes);
        ct2Tensor = turbulenceTensor2(nModes);

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(nModes);
            bcVelMat = bcVelocityMat(nModes);
        }
    }

    // Export the matrices
    if (para->exportPython)
    {
        ITHACAstream::exportMatrix(bMatrix, "b", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(kMatrix, "k", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(btMatrix, "bt", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(convTensor, "c", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct1Tensor, "ct1", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct2Tensor, "ct2", "python",
                                   "./ITHACAoutput/Matrices/");
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportMatrix(bMatrix, "b", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(kMatrix, "k", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(btMatrix, "bt", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(convTensor, "c", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct1Tensor, "ct1", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct2Tensor, "ct2", "matlab",
                                   "./ITHACAoutput/Matrices/");
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportMatrix(bMatrix, "b", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(kMatrix, "k", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(btMatrix, "bt", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(convTensor, "c", "eigen",
                                   "./ITHACAoutput/Matrices/c");
        ITHACAstream::exportTensor(ct1Tensor, "ct1_", "eigen",
                                   "./ITHACAoutput/Matrices/ct1");
        ITHACAstream::exportTensor(ct2Tensor, "ct2_", "eigen",
                                   "./ITHACAoutput/Matrices/ct2");
    }

    bTotalMatrix = bMatrix + btMatrix;
    cTotalTensor.resize(nModes, nModes, nModes);
    cTotalTensor = convTensor - ct1Tensor - ct2Tensor;
}

void UnsteadyNSTurbIntrusive::projectPPE(fileName folder, label nUModes,
        label nPModes)
{
    NUmodes = nUModes;
    NPmodes = nPModes;
    L_U_SUPmodes.resize(0);

    if (nUModes != 0)
    {
        for (label k = 0; k < nUModes; k++)
        {
            L_U_SUPmodes.append(Umodes[k]);
        }
    }

    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {
        word bStr = "b_" + name(nUModes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + bStr))
        {
            ITHACAstream::ReadDenseMatrix(bMatrix, "./ITHACAoutput/Matrices/", bStr);
        }
        else
        {
            bMatrix = diffusiveTerm(nUModes);
        }

        word btStr = "bt_" + name(nUModes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + btStr))
        {
            ITHACAstream::ReadDenseMatrix(btMatrix, "./ITHACAoutput/Matrices/", btStr);
        }
        else
        {
            btMatrix = btTurbulence(nUModes);
        }

        word kStr = "k_" + name(nUModes) + "_" + name(nPModes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + kStr))
        {
            ITHACAstream::ReadDenseMatrix(kMatrix, "./ITHACAoutput/Matrices/", kStr);
        }
        else
        {
            kMatrix = pressureGradientTerm(nUModes, nPModes);
        }

        word D_str = "D_" + name(nPModes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + D_str))
        {
            ITHACAstream::ReadDenseMatrix(D_matrix, "./ITHACAoutput/Matrices/", D_str);
        }
        else
        {
            D_matrix = laplacianPressure(nPModes);
        }

        word bc1_str = "BC1_" + name(nUModes) + "_" + name(nPModes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + bc1_str))
        {
            ITHACAstream::ReadDenseMatrix(BC1_matrix, "./ITHACAoutput/Matrices/", bc1_str);
        }
        else
        {
            BC1_matrix = pressureBC1(nUModes, nPModes);
        }

        word bc2_str = "BC2_" + name(nUModes) + "_" + name(nPModes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + bc2_str))
        {
            ITHACAstream::ReadDenseTensor(bc2Tensor, "./ITHACAoutput/Matrices/", bc2_str);
        }
        else
        {
            bc2Tensor = pressureBC2(nUModes, nPModes);
        }

        word bc3_str = "BC3_" + name(nUModes) + "_" + name(nPModes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + bc3_str))
        {
            ITHACAstream::ReadDenseMatrix(BC3_matrix, "./ITHACAoutput/Matrices/", bc3_str);
        }
        else
        {
            BC3_matrix = pressureBC3(nUModes, nPModes);
        }

        word cStr = "c_" + name(nUModes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + cStr))
        {
            ITHACAstream::ReadDenseTensor(convTensor, "./ITHACAoutput/Matrices/", cStr);
        }
        else
        {
            convTensor = convectiveTerm(nUModes);
        }

        word ct1Str = "ct1_" + name(nUModes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct1Str))
        {
            ITHACAstream::ReadDenseTensor(ct1Tensor, "./ITHACAoutput/Matrices/", ct1Str);
        }
        else
        {
            ct1Tensor = turbulenceTensor1(nUModes);
        }

        word ct2Str = "ct2_" + name(nUModes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct2Str))
        {
            ITHACAstream::ReadDenseTensor(ct2Tensor, "./ITHACAoutput/Matrices/", ct2Str);
        }
        else
        {
            ct2Tensor = turbulenceTensor2(nUModes);
        }

        word ct1PPEStr = "ct1PPE_" + name(nUModes) + "_" + name(nPModes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct1PPEStr))
        {
            ITHACAstream::ReadDenseTensor(ct1PPETensor, "./ITHACAoutput/Matrices/",
                                          ct1PPEStr);
        }
        else
        {
            ct1PPETensor = turbulencePPETensor1(nUModes, nPModes);
        }

        word ct2PPEStr = "ct2PPE_" + name(nUModes) + "_" + name(nPModes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct2PPEStr))
        {
            ITHACAstream::ReadDenseTensor(ct2PPETensor, "./ITHACAoutput/Matrices/",
                                          ct2PPEStr);
        }
        else
        {
            ct2PPETensor = turbulencePPETensor2(nUModes, nPModes);
        }

        word G_str = "G_" + name(nUModes) + "_" + name(nPModes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + G_str))
        {
            ITHACAstream::ReadDenseTensor(gTensor, "./ITHACAoutput/Matrices/", G_str);
        }
        else
        {
            gTensor = divMomentum(nUModes, nPModes);
        }

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(nUModes);
            bcVelMat = bcVelocityMat(nUModes);
        }
    }
    else
    {
        bMatrix = diffusiveTerm(nUModes);
        convTensor = convectiveTerm(nUModes);
        kMatrix = pressureGradientTerm(nUModes, nPModes);
        btMatrix = btTurbulence(nUModes);
        ct1Tensor = turbulenceTensor1(nUModes);
        ct2Tensor = turbulenceTensor2(nUModes);
        D_matrix = laplacianPressure(nPModes);
        ct1PPETensor = turbulencePPETensor1(nUModes, nPModes);
        ct2PPETensor = turbulencePPETensor2(nUModes, nPModes);
        gTensor = divMomentum(nUModes, nPModes);
        BC1_matrix = pressureBC1(nUModes, nPModes);
        bc2Tensor = pressureBC2(nUModes, nPModes);
        BC3_matrix = pressureBC3(nUModes, nPModes);

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(nUModes);
            bcVelMat = bcVelocityMat(nUModes);
        }
    }

    // Export the matrices
    if (para->exportPython)
    {
        ITHACAstream::exportMatrix(bMatrix, "b", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(kMatrix, "k", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(btMatrix, "bt", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(D_matrix, "D", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC1_matrix, "BC1", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC3_matrix, "BC3", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct1PPETensor, "ct1PPE", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct2PPETensor, "ct2PPE", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(gTensor, "G", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(bc2Tensor, "BC2", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(convTensor, "c", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct1Tensor, "ct1", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct2Tensor, "ct2", "python",
                                   "./ITHACAoutput/Matrices/");
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportMatrix(bMatrix, "b", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(kMatrix, "k", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(btMatrix, "bt", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(D_matrix, "D", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC1_matrix, "BC1", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC3_matrix, "BC3", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct1PPETensor, "ct1PPE", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct2PPETensor, "ct2PPE", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(gTensor, "G", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(bc2Tensor, "BC2", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(convTensor, "c", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct1Tensor, "ct1", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct2Tensor, "ct2", "matlab",
                                   "./ITHACAoutput/Matrices/");
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportMatrix(bMatrix, "b", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(kMatrix, "k", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(btMatrix, "bt", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(D_matrix, "D", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC1_matrix, "BC1", "eigen",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC3_matrix, "BC3", "eigen",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct1PPETensor, "ct1PPE_", "eigen",
                                   "./ITHACAoutput/Matrices/ct1PPE");
        ITHACAstream::exportTensor(ct2PPETensor, "ct2PPE_", "eigen",
                                   "./ITHACAoutput/Matrices/ct2PPE");
        ITHACAstream::exportTensor(gTensor, "G", "eigen",
                                   "./ITHACAoutput/Matrices/G");
        ITHACAstream::exportTensor(bc2Tensor, "BC2_", "eigen",
                                   "./ITHACAoutput/Matrices/BC2");
        ITHACAstream::exportTensor(convTensor, "c", "eigen",
                                   "./ITHACAoutput/Matrices/c");
        ITHACAstream::exportTensor(ct1Tensor, "ct1_", "eigen",
                                   "./ITHACAoutput/Matrices/ct1");
        ITHACAstream::exportTensor(ct2Tensor, "ct2_", "eigen",
                                   "./ITHACAoutput/Matrices/ct2");
    }

    bTotalMatrix = bMatrix + btMatrix;
    cTotalTensor.resize(nUModes, nUModes, nUModes);
    cTotalTensor = convTensor - ct1Tensor - ct2Tensor;
    cTotalPPETensor.resize(nPModes, nUModes, nUModes);
    cTotalPPETensor = ct1PPETensor + ct2PPETensor;
}

Eigen::MatrixXd UnsteadyNSTurbIntrusive::diffusiveTerm(label nModes)
{
    Eigen::MatrixXd bMatrix;
    bMatrix.resize(nModes, nModes);

    // Project everything
    for (label i = 0; i < nModes; i++)
    {
        for (label j = 0; j < nModes; j++)
        {
            bMatrix(i, j) = fvc::domainIntegrate(Umodes[i] & fvc::laplacian(
                    dimensionedScalar("1", dimless, 1), Umodes[j])).value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(bMatrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(bMatrix, "./ITHACAoutput/Matrices/",
                                  "b_" + name(nModes));
    return bMatrix;
}

Eigen::Tensor<double, 3> UnsteadyNSTurbIntrusive::convectiveTerm(label nModes)
{
    Eigen::Tensor<double, 3> convTensor;
    convTensor.resize(nModes, nModes, nModes);

    for (label i = 0; i < nModes; i++)
    {
        for (label j = 0; j < nModes; j++)
        {
            for (label k = 0; k < nModes; k++)
            {
                convTensor(i, j, k) = fvc::domainIntegrate(Umodes[i] & fvc::div(
                                          linearInterpolate(Umodes[j]) & Umodes[j].mesh().Sf(),
                                          Umodes[k])).value();
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(convTensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(convTensor, "./ITHACAoutput/Matrices/",
                                  "c_" + name(nModes) + "_t");
    return convTensor;
}

Eigen::MatrixXd UnsteadyNSTurbIntrusive::pressureGradientTerm(label nModes)
{
    Eigen::MatrixXd kMatrix(nModes, nModes);

    // Project everything
    for (label i = 0; i < nModes; i++)
    {
        for (label j = 0; j < nModes; j++)
        {
            kMatrix(i, j) = fvc::domainIntegrate(Umodes[i] & fvc::grad(
                    Pmodes[j])).value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(kMatrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(kMatrix, "./ITHACAoutput/Matrices/",
                                  "k_" + name(nModes));
    return kMatrix;
}

Eigen::MatrixXd UnsteadyNSTurbIntrusive::pressureGradientTerm(label nUModes,
        label nPModes)
{
    Eigen::MatrixXd kMatrix(nUModes, nPModes);

    // Project everything
    for (label i = 0; i < nUModes; i++)
    {
        for (label j = 0; j < nPModes; j++)
        {
            kMatrix(i, j) = fvc::domainIntegrate(Umodes[i] & fvc::grad(
                    Pmodes[j])).value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(kMatrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(kMatrix, "./ITHACAoutput/Matrices/",
                                  "k_" + name(nUModes) + "_" + name(nPModes));
    return kMatrix;
}

Eigen::MatrixXd UnsteadyNSTurbIntrusive::laplacianPressure(label nPModes)
{
    Eigen::MatrixXd dMatrix(nPModes, nPModes);

    // Project everything
    for (label i = 0; i < nPModes; i++)
    {
        for (label j = 0; j < nPModes; j++)
        {
            dMatrix(i, j) = fvc::domainIntegrate(fvc::grad(Pmodes[i])&fvc::grad(
                    Pmodes[j])).value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(dMatrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(dMatrix, "./ITHACAoutput/Matrices/",
                                  "D_" + name(nPModes));
    return dMatrix;
}

Eigen::Tensor<double, 3> UnsteadyNSTurbIntrusive::divMomentum(label nUModes,
        label nPModes)
{
    label g1Size = nPModes;
    label g2Size = nUModes;
    Eigen::Tensor<double, 3> gTensor;
    gTensor.resize(g1Size, g2Size, g2Size);

    for (label i = 0; i < g1Size; i++)
    {
        for (label j = 0; j < g2Size; j++)
        {
            for (label k = 0; k < g2Size; k++)
            {
                gTensor(i, j, k) = fvc::domainIntegrate(fvc::grad(Pmodes[i]) & (fvc::div(
                        fvc::interpolate(Umodes[j]) & Umodes[j].mesh().Sf(),
                        Umodes[k]))).value();
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(gTensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(gTensor, "./ITHACAoutput/Matrices/",
                                  "G_" + name(nUModes) + "_" + name(
                                      nPModes) + "_t");
    return gTensor;
}

Eigen::MatrixXd UnsteadyNSTurbIntrusive::pressureBC1(label nUModes,
        label nPModes)
{
    label P_BC1size = nPModes;
    label P_BC2size = nUModes;
    Eigen::MatrixXd BC1_matrix(P_BC1size, P_BC2size);
    fvMesh& mesh = _mesh();

    for (label i = 0; i < P_BC1size; i++)
    {
        for (label j = 0; j < P_BC2size; j++)
        {
            surfaceScalarField lpl((fvc::interpolate(fvc::laplacian(
                                        Umodes[j]))&mesh.Sf())*fvc::interpolate(Pmodes[i]));
            double s = 0;

            for (label k = 0; k < lpl.boundaryField().size(); k++)
            {
                s += gSum(lpl.boundaryField()[k]);
            }

            BC1_matrix(i, j) = s;
        }
    }

    if (Pstream::parRun())
    {
        reduce(BC1_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(BC1_matrix, "./ITHACAoutput/Matrices/",
                                  "BC1_" + name(nUModes) + "_" + name(nPModes));
    return BC1_matrix;
}

Eigen::Tensor<double, 3 > UnsteadyNSTurbIntrusive::pressureBC2(label nUModes,
        label nPModes)
{
    label pressureBC1Size = nPModes;
    label pressureBC2Size = nUModes;
    Eigen::Tensor<double, 3 > bc2Tensor;
    fvMesh& mesh = _mesh();
    bc2Tensor.resize(pressureBC1Size, pressureBC2Size, pressureBC2Size);

    for (label i = 0; i < pressureBC1Size; i++)
    {
        for (label j = 0; j < pressureBC2Size; j++)
        {
            for (label k = 0; k < pressureBC2Size; k++)
            {
                surfaceScalarField div_m(fvc::interpolate(fvc::div(fvc::interpolate(
                                             Umodes[j]) & mesh.Sf(),
                                         Umodes[k]))&mesh.Sf()*fvc::interpolate(Pmodes[i]));
                double s = 0;

                for (label k = 0; k < div_m.boundaryField().size(); k++)
                {
                    s += gSum(div_m.boundaryField()[k]);
                }

                bc2Tensor(i, j, k) = s;
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(bc2Tensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(bc2Tensor, "./ITHACAoutput/Matrices/",
                                  "BC2_" + name(nUModes) + "_" + name(nPModes) + "_t");
    return bc2Tensor;
}

Eigen::MatrixXd UnsteadyNSTurbIntrusive::pressureBC3(label nUModes,
        label nPModes)
{
    label P3_BC1size = nPModes;
    label P3_BC2size = nUModes;
    Eigen::MatrixXd BC3_matrix(P3_BC1size, P3_BC2size);
    fvMesh& mesh = _mesh();
    surfaceVectorField n(mesh.Sf() / mesh.magSf());

    for (label i = 0; i < P3_BC1size; i++)
    {
        for (label j = 0; j < P3_BC2size; j++)
        {
            surfaceVectorField BC3 = fvc::interpolate(fvc::curl(Umodes[j]));
            surfaceVectorField BC4 = n ^ fvc::interpolate(fvc::grad(Pmodes[i]));
            surfaceScalarField BC5 = (BC3 & BC4) * mesh.magSf();
            double s = 0;

            for (label k = 0; k < BC5.boundaryField().size(); k++)
            {
                s += gSum(BC5.boundaryField()[k]);
            }

            BC3_matrix(i, j) = s;
        }
    }

    if (Pstream::parRun())
    {
        reduce(BC3_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(BC3_matrix, "./ITHACAoutput/Matrices/",
                                  "BC3_" + name(nUModes) + "_" + name(nPModes));
    return BC3_matrix;
}

Eigen::MatrixXd UnsteadyNSTurbIntrusive::pressureBC4(label nUModes,
        label nPModes)
{
    label P4_BC1size = nPModes;
    label P4_BC2size = nUModes;
    Eigen::MatrixXd BC4_matrix(P4_BC1size, P4_BC2size);
    fvMesh& mesh = _mesh();
    surfaceVectorField n(mesh.Sf() / mesh.magSf());

    for (label i = 0; i < P4_BC1size; i++)
    {
        for (label j = 0; j < P4_BC2size; j++)
        {
            surfaceScalarField BC3 = fvc::interpolate(Pmodes[i]);
            surfaceScalarField BC4 = n & fvc::interpolate(Umodes[j]);
            surfaceScalarField BC5 = (BC3 * BC4) * mesh.magSf();
            double s = 0;

            for (label k = 0; k < BC5.boundaryField().size(); k++)
            {
                s += gSum(BC5.boundaryField()[k]);
            }

            BC4_matrix(i, j) = s;
        }
    }

    if (Pstream::parRun())
    {
        reduce(BC4_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(BC4_matrix, "./ITHACAoutput/Matrices/",
                                  "BC4_" + name(nUModes) + "_" + name(nPModes));
    return BC4_matrix;
}

List< Eigen::MatrixXd > UnsteadyNSTurbIntrusive::bcVelocityVec(label nModes)
{
    List < Eigen::MatrixXd > bcVelVec(inletIndex.rows());

    for (label j = 0; j < inletIndex.rows(); j++)
    {
        bcVelVec[j].resize(nModes, 1);
    }

    for (label k = 0; k < inletIndex.rows(); k++)
    {
        label BCind = inletIndex(k, 0);
        label BCcomp = inletIndex(k, 1);

        for (label i = 0; i < nModes; i++)
        {
            bcVelVec[k](i, 0) = gSum(Umodes[i].boundaryField()[BCind].component(
                                         BCcomp));
        }
    }

    ITHACAstream::exportMatrix(bcVelVec, "bcVelVec", "eigen",
                               "./ITHACAoutput/Matrices/bcVelVec");
    return bcVelVec;
}

List< Eigen::MatrixXd > UnsteadyNSTurbIntrusive::bcVelocityMat(label nModes)
{
    label BCUsize = inletIndex.rows();
    List < Eigen::MatrixXd > bcVelMat(BCUsize);

    for (label j = 0; j < inletIndex.rows(); j++)
    {
        bcVelMat[j].resize(nModes, nModes);
    }

    for (label k = 0; k < inletIndex.rows(); k++)
    {
        label BCind = inletIndex(k, 0);
        label BCcomp = inletIndex(k, 1);

        for (label i = 0; i < nModes; i++)
        {
            for (label j = 0; j < nModes; j++)
            {
                bcVelMat[k](i, j) = gSum(Umodes[i].boundaryField()[BCind].component(
                                             BCcomp) *
                                         Umodes[j].boundaryField()[BCind].component(BCcomp));
            }
        }
    }

    ITHACAstream::exportMatrix(bcVelMat, "bcVelMat", "eigen",
                               "./ITHACAoutput/Matrices/bcVelMat");
    return bcVelMat;
}

Eigen::Tensor<double, 3> UnsteadyNSTurbIntrusive::turbulencePPETensor1(
    label nUModes, label nPModes)
{
    Eigen::Tensor<double, 3> ct1PPETensor;
    ct1PPETensor.resize(nPModes, nUModes, nUModes);

    for (label i = 0; i < nPModes; i++)
    {
        for (label j = 0; j < nUModes; j++)
        {
            for (label k = 0; k < nUModes; k++)
            {
                ct1PPETensor(i, j, k) = fvc::domainIntegrate(fvc::grad(Pmodes[i]) & (
                                            fvc::laplacian(
                                                nutModes[j], Umodes[k]))).value();
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(ct1PPETensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(ct1PPETensor, "./ITHACAoutput/Matrices/",
                                  "ct1PPE_" + name(nUModes) + "_" + name(nPModes) + "_t");
    return ct1PPETensor;
}

Eigen::Tensor<double, 3> UnsteadyNSTurbIntrusive::turbulencePPETensor2(
    label nUModes, label nPModes)
{
    Eigen::Tensor<double, 3> ct2PPETensor;
    ct2PPETensor.resize(nPModes, nUModes, nUModes);

    for (label i = 0; i < nPModes; i++)
    {
        for (label j = 0; j < nUModes; j++)
        {
            for (label k = 0; k < nUModes; k++)
            {
                ct2PPETensor(i, j, k) = fvc::domainIntegrate(fvc::grad(Pmodes[i]) & ((fvc::div(
                                            nutModes[j] * dev2((fvc::grad(Umodes[k]))().T()))))).value();
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(ct2PPETensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(ct2PPETensor, "./ITHACAoutput/Matrices/",
                                  "ct2PPE_" + name(nUModes) + "_" + name(nPModes) + "_t");
    return ct2PPETensor;
}
