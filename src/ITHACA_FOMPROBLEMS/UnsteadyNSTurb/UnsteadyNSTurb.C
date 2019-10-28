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
    para = new ITHACAparameters;
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    supex = ITHACAutilities::check_sup();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void UnsteadyNSTurb::truthSolve(List<scalar> mu_now)
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

List < Eigen::MatrixXd > UnsteadyNSTurb::turbulenceTerm1(label NUmodes,
        label NSUPmodes, label nNutModes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    List < Eigen::MatrixXd > ct1Matrix;
    ct1Matrix.setSize(cSize);

    for (label j = 0; j < cSize; j++)
    {
        ct1Matrix[j].resize(nNutModes, cSize);
        ct1Matrix[j] = ct1Matrix[j] * 0;
    }

    for (label i = 0; i < cSize; i++)
    {
        Info << "Filling layer number " << i + 1 << " in the matrix ct1Matrix" << endl;

        for (label j = 0; j < nNutModes; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                ct1Matrix[i](j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::laplacian(
                                         nutModes[j], L_U_SUPmodes[k])).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(ct1Matrix, "ct1Matrix", "python",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(ct1Matrix, "ct1Matrix", "matlab",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(ct1Matrix, "ct1Matrix", "eigen",
                               "./ITHACAoutput/Matrices/ct1");
    return ct1Matrix;
}

Eigen::Tensor<double, 3> UnsteadyNSTurb::turbulenceTensor1(label NUmodes,
        label NSUPmodes, label nNutModes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct1Tensor;
    ct1Tensor.resize(cSize, nNutModes, cSize);

    for (label i = 0; i < cSize; i++)
    {
        for (label j = 0; j < nNutModes; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                ct1Tensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::laplacian(
                                         nutModes[j], L_U_SUPmodes[k])).value();
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(ct1Tensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(ct1Tensor, "./ITHACAoutput/Matrices/",
                                  "ct1_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(nNutModes) + "_t");
    return ct1Tensor;
}





List < Eigen::MatrixXd > UnsteadyNSTurb::turbulenceTerm2(label NUmodes,
        label NSUPmodes, label nNutModes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    List < Eigen::MatrixXd > ct2Matrix;
    ct2Matrix.setSize(cSize);

    for (label j = 0; j < cSize; j++)
    {
        ct2Matrix[j].resize(nNutModes, cSize);
        ct2Matrix[j] = ct2Matrix[j] * 0;
    }

    for (label i = 0; i < cSize; i++)
    {
        Info << "Filling layer number " << i + 1 << " in the matrix ct2Matrix" << endl;

        for (label j = 0; j < nNutModes; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                ct2Matrix[i](j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & (fvc::div(
                                         nutModes[j] * dev((fvc::grad(L_U_SUPmodes[k]))().T())))).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(ct2Matrix, "ct2Matrix", "python",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(ct2Matrix, "ct2Matrix", "matlab",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(ct2Matrix, "ct2Matrix", "eigen",
                               "./ITHACAoutput/Matrices/ct2");
    return ct2Matrix;
}

Eigen::Tensor<double, 3> UnsteadyNSTurb::turbulenceTensor2(label NUmodes,
        label NSUPmodes, label nNutModes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct2Tensor;
    ct2Tensor.resize(cSize, nNutModes, cSize);

    for (label i = 0; i < cSize; i++)
    {
        for (label j = 0; j < nNutModes; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                ct2Tensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & (fvc::div(
                                         nutModes[j] * dev((fvc::grad(L_U_SUPmodes[k]))().T())))).value();
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(ct2Tensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(ct2Tensor, "./ITHACAoutput/Matrices/",
                                  "ct2_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(nNutModes) + "_t");
    return ct2Tensor;
}

Eigen::MatrixXd UnsteadyNSTurb::btTurbulence(label NUmodes, label NSUPmodes)
{
    label btSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd btMatrix(btSize, btSize);
    btMatrix = btMatrix * 0;

    // Project everything
    for (label i = 0; i < btSize; i++)
    {
        for (label j = 0; j < btSize; j++)
        {
            btMatrix(i, j) = fvc::domainIntegrate(L_U_SUPmodes[i] & (fvc::div(dev((T(
                    fvc::grad(
                        L_U_SUPmodes[j]))))))).value();
        }
    }

    // Export the matrix
    ITHACAstream::SaveDenseMatrix(btMatrix, "./ITHACAoutput/Matrices/",
                                  "bt_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));
    return btMatrix;
}

void UnsteadyNSTurb::projectSUP(fileName folder, label NU, label NP, label NSUP,
                                label Nnut)
{
    NUmodes = NU;
    NPmodes = NP;
    NSUPmodes = NSUP;
    nNutModes = Nnut;
    L_U_SUPmodes.resize(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            L_U_SUPmodes.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            L_U_SUPmodes.append(Umodes[k]);
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            L_U_SUPmodes.append(supmodes[k]);
        }
    }

    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {
        word bStr = "B_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                        NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + bStr))
        {
            ITHACAstream::ReadDenseMatrix(B_matrix, "./ITHACAoutput/Matrices/", bStr);
        }
        else
        {
            B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        }

        word btStr = "bt_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + btStr))
        {
            ITHACAstream::ReadDenseMatrix(btMatrix, "./ITHACAoutput/Matrices/", btStr);
        }
        else
        {
            btMatrix = btTurbulence(NUmodes, NSUPmodes);
        }

        word kStr = "K_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                        NSUPmodes) + "_" + name(NPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + kStr))
        {
            ITHACAstream::ReadDenseMatrix(K_matrix, "./ITHACAoutput/Matrices/", kStr);
        }
        else
        {
            K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        }

        word pStr = "P_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                        NSUPmodes) + "_" + name(NPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + pStr))
        {
            ITHACAstream::ReadDenseMatrix(P_matrix, "./ITHACAoutput/Matrices/", pStr);
        }
        else
        {
            P_matrix = divergence_term(NUmodes, NPmodes, NSUPmodes);
        }

        word mStr = "M_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                        NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + mStr))
        {
            ITHACAstream::ReadDenseMatrix(M_matrix, "./ITHACAoutput/Matrices/", mStr);
        }
        else
        {
            M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
        }

        word C_str = "C_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + C_str))
        {
            ITHACAstream::ReadDenseTensor(C_tensor, "./ITHACAoutput/Matrices/", C_str);
        }
        else
        {
            C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        }

        word ct1Str = "ct1_" + name(liftfield.size()) + "_" + name(
                          NUmodes) + "_" + name(
                          NSUPmodes) + "_" + name(nNutModes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct1Str))
        {
            ITHACAstream::ReadDenseTensor(ct1Tensor, "./ITHACAoutput/Matrices/", ct1Str);
        }
        else
        {
            ct1Tensor = turbulenceTensor1(NUmodes, NSUPmodes, nNutModes);
        }

        word ct2Str = "ct2_" + name(liftfield.size()) + "_" + name(
                          NUmodes) + "_" + name(
                          NSUPmodes) + "_" + name(nNutModes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct2Str))
        {
            ITHACAstream::ReadDenseTensor(ct2Tensor, "./ITHACAoutput/Matrices/", ct2Str);
        }
        else
        {
            ct2Tensor = turbulenceTensor2(NUmodes, NSUPmodes, nNutModes);
        }

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }
    }
    else
    {
        B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        P_matrix = divergence_term(NUmodes, NPmodes, NSUPmodes);
        M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
        btMatrix = btTurbulence(NUmodes, NSUPmodes);
        ct1Tensor = turbulenceTensor1(NUmodes, NSUPmodes, nNutModes);
        ct2Tensor = turbulenceTensor2(NUmodes, NSUPmodes, nNutModes);

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }
    }

    // Export the matrices
    if (para->exportPython)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(P_matrix, "P", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(btMatrix, "bt", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor, "C", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct1Tensor, "ct1", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct2Tensor, "ct2", "python",
                                   "./ITHACAoutput/Matrices/");
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(P_matrix, "P", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(btMatrix, "bt", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor, "C", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct1Tensor, "ct1", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct2Tensor, "ct2", "matlab",
                                   "./ITHACAoutput/Matrices/");
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(P_matrix, "P", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(btMatrix, "bt", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor, "C", "eigen", "./ITHACAoutput/Matrices/C");
        ITHACAstream::exportTensor(ct1Tensor, "ct1_", "eigen",
                                   "./ITHACAoutput/Matrices/ct1");
        ITHACAstream::exportTensor(ct2Tensor, "ct2_", "eigen",
                                   "./ITHACAoutput/Matrices/ct2");
    }

    bTotalMatrix = B_matrix + btMatrix;
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    cTotalTensor.resize(cSize, nNutModes, cSize);
    cTotalTensor = ct1Tensor + ct2Tensor;
    // Get the coeffs for interpolation (the orthonormal one is used because basis are orthogonal)
    coeffL2 = ITHACAutilities::get_coeffs_ortho(nutFields,
              nutModes, nNutModes);
    ITHACAstream::exportMatrix(coeffL2, "coeffL2", "python",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(coeffL2, "coeffL2", "matlab",
                               "./ITHACAoutput/Matrices/");
    // Export the matrix
    ITHACAstream::SaveDenseMatrix(coeffL2, "./ITHACAoutput/Matrices/",
                                  "coeffL2_nut_" + name(nNutModes));
    samples.resize(nNutModes);
    rbfSplines.resize(nNutModes);
    Eigen::MatrixXd weights;

    for (int i = 0; i < nNutModes; i++)
    {
        word weightName = "wRBF_N" + name(i + 1) + "_" + name(liftfield.size()) + "_"
                          + name(NUmodes) + "_" + name(NSUPmodes) ;

        if (ITHACAutilities::check_file("./ITHACAoutput/weightsSUP/" + weightName))
        {
            samples[i] = new SPLINTER::DataTable(1, 1);

            for (int j = 0; j < coeffL2.cols(); j++)
            {
                samples[i]->addSample(mu.row(j), coeffL2(i, j));
            }

            ITHACAstream::ReadDenseMatrix(weights, "./ITHACAoutput/weightsSUP/",
                                          weightName);
            rbfSplines[i] = new SPLINTER::RBFSpline(*samples[i],
                                                    SPLINTER::RadialBasisFunctionType::GAUSSIAN, weights);
            std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
        }
        else
        {
            samples[i] = new SPLINTER::DataTable(1, 1);

            for (int j = 0; j < coeffL2.cols(); j++)
            {
                samples[i]->addSample(mu.row(j), coeffL2(i, j));
            }

            rbfSplines[i] = new SPLINTER::RBFSpline(*samples[i],
                                                    SPLINTER::RadialBasisFunctionType::GAUSSIAN);
            ITHACAstream::SaveDenseMatrix(rbfSplines[i]->weights,
                                          "./ITHACAoutput/weightsSUP/", weightName);
            std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
        }
    }
}


void UnsteadyNSTurb::projectPPE(fileName folder, label NU, label NP, label NSUP,
                                label Nnut)
{
    NUmodes = NU;
    NPmodes = NP;
    NSUPmodes = 0;
    nNutModes = Nnut;
    L_U_SUPmodes.resize(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            L_U_SUPmodes.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            L_U_SUPmodes.append(Umodes[k]);
        }
    }

    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {
        word B_str = "B_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + B_str))
        {
            ITHACAstream::ReadDenseMatrix(B_matrix, "./ITHACAoutput/Matrices/", B_str);
        }
        else
        {
            B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        }

        word btStr = "bt_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + btStr))
        {
            ITHACAstream::ReadDenseMatrix(btMatrix, "./ITHACAoutput/Matrices/", btStr);
        }
        else
        {
            btMatrix = btTurbulence(NUmodes, NSUPmodes);
        }

        word K_str = "K_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_" + name(NPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + K_str))
        {
            ITHACAstream::ReadDenseMatrix(K_matrix, "./ITHACAoutput/Matrices/", K_str);
        }
        else
        {
            K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        }

        word M_str = "M_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + M_str))
        {
            ITHACAstream::ReadDenseMatrix(M_matrix, "./ITHACAoutput/Matrices/", M_str);
        }
        else
        {
            M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
        }

        word D_str = "D_" + name(NPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + D_str))
        {
            ITHACAstream::ReadDenseMatrix(D_matrix, "./ITHACAoutput/Matrices/", D_str);
        }
        else
        {
            D_matrix = laplacian_pressure(NPmodes);
        }

        word bc1_str = "BC1_" + name(liftfield.size()) + "_" + name(
                           NUmodes) + "_" + name(NPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + bc1_str))
        {
            ITHACAstream::ReadDenseMatrix(BC1_matrix, "./ITHACAoutput/Matrices/", bc1_str);
        }
        else
        {
            BC1_matrix = pressure_BC1(NUmodes, NPmodes);
        }

        word bc2_str = "BC2_" + name(liftfield.size()) + "_" + name(
                           NUmodes) + "_" + name(
                           NSUPmodes) + "_" + name(NPmodes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + bc2_str))
        {
            ITHACAstream::ReadDenseTensor(bc2Tensor, "./ITHACAoutput/Matrices/", bc2_str);
        }
        else
        {
            bc2Tensor = pressureBC2(NUmodes, NPmodes);
        }

        word bc3_str = "BC3_" + name(liftfield.size()) + "_" + name(
                           NUmodes) + "_" + name(NPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + bc3_str))
        {
            ITHACAstream::ReadDenseMatrix(BC3_matrix, "./ITHACAoutput/Matrices/", bc3_str);
        }
        else
        {
            BC3_matrix = pressure_BC3(NUmodes, NPmodes);
        }

        word C_str = "C_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + C_str))
        {
            ITHACAstream::ReadDenseTensor(C_tensor, "./ITHACAoutput/Matrices/", C_str);
        }
        else
        {
            C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        }

        word ct1Str = "ct1_" + name(liftfield.size()) + "_" + name(
                          NUmodes) + "_" + name(
                          NSUPmodes) + "_" + name(nNutModes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct1Str))
        {
            ITHACAstream::ReadDenseTensor(ct1Tensor, "./ITHACAoutput/Matrices/", ct1Str);
        }
        else
        {
            ct1Tensor = turbulenceTensor1(NUmodes, NSUPmodes, nNutModes);
        }

        word ct2Str = "ct2_" + name(liftfield.size()) + "_" + name(
                          NUmodes) + "_" + name(
                          NSUPmodes) + "_" + name(nNutModes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct2Str))
        {
            ITHACAstream::ReadDenseTensor(ct2Tensor, "./ITHACAoutput/Matrices/", ct2Str);
        }
        else
        {
            ct2Tensor = turbulenceTensor2(NUmodes, NSUPmodes, nNutModes);
        }

        word G_str = "G_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_" + name(NPmodes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + G_str))
        {
            ITHACAstream::ReadDenseTensor(gTensor, "./ITHACAoutput/Matrices/", G_str);
        }
        else
        {
            gTensor = divMomentum(NUmodes, NPmodes);
        }

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }
    }
    else
    {
        B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        D_matrix = laplacian_pressure(NPmodes);
        gTensor = divMomentum(NUmodes, NPmodes);
        BC1_matrix = pressure_BC1(NUmodes, NPmodes);
        bc2Tensor = pressureBC2(NUmodes, NPmodes);
        BC3_matrix = pressure_BC3(NUmodes, NPmodes);
        btMatrix = btTurbulence(NUmodes, NSUPmodes);
        ct1Tensor = turbulenceTensor1(NUmodes, NSUPmodes, nNutModes);
        ct2Tensor = turbulenceTensor2(NUmodes, NSUPmodes, nNutModes);

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }
    }

    // Export the matrices
    if (para->exportPython)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(D_matrix, "D", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC1_matrix, "BC1", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC3_matrix, "BC3", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor, "C", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(gTensor, "G", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(bc2Tensor, "BC2", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct1Tensor, "ct1", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct2Tensor, "ct2", "python",
                                   "./ITHACAoutput/Matrices/");
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(D_matrix, "D", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC1_matrix, "BC1", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC3_matrix, "BC3", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor, "C", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(gTensor, "G", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(bc2Tensor, "BC2", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct1Tensor, "ct1", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct2Tensor, "ct2", "matlab",
                                   "./ITHACAoutput/Matrices/");
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(D_matrix, "D", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC1_matrix, "BC1", "eigen",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC3_matrix, "BC3", "eigen",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor, "C", "eigen",
                                   "./ITHACAoutput/Matrices/C");
        ITHACAstream::exportTensor(gTensor, "G", "eigen",
                                   "./ITHACAoutput/Matrices/G");
        ITHACAstream::exportTensor(bc2Tensor, "BC2_", "eigen",
                                   "./ITHACAoutput/Matrices/BC2");
        ITHACAstream::exportMatrix(btMatrix, "bt", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct1Tensor, "ct1_", "eigen",
                                   "./ITHACAoutput/Matrices/ct1");
        ITHACAstream::exportTensor(ct2Tensor, "ct2_", "eigen",
                                   "./ITHACAoutput/Matrices/ct2");
    }

    bTotalMatrix = B_matrix + btMatrix;
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    cTotalTensor.resize(cSize, nNutModes, cSize);
    cTotalTensor = ct1Tensor + ct2Tensor;
    // Get the coeffs for interpolation (the orthonormal one is used because basis are orthogonal)
    coeffL2 = ITHACAutilities::get_coeffs_ortho(nutFields,
              nutModes, nNutModes);
    ITHACAstream::exportMatrix(coeffL2, "coeffL2", "python",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(coeffL2, "coeffL2", "matlab",
                               "./ITHACAoutput/Matrices/");
    // Export the matrix
    ITHACAstream::SaveDenseMatrix(coeffL2, "./ITHACAoutput/Matrices/",
                                  "coeffL2_nut_" + name(nNutModes));
    samples.resize(nNutModes);
    rbfSplines.resize(nNutModes);
    Eigen::MatrixXd weights;

    for (int i = 0; i < nNutModes; i++)
    {
        word weightName = "wRBF_N" + name(i + 1) + "_" + name(liftfield.size()) + "_"
                          + name(NUmodes) + "_" + name(NSUPmodes) ;

        if (ITHACAutilities::check_file("./ITHACAoutput/weightsPPE/" + weightName))
        {
            samples[i] = new SPLINTER::DataTable(1, 1);

            for (int j = 0; j < coeffL2.cols(); j++)
            {
                samples[i]->addSample(mu.row(j), coeffL2(i, j));
            }

            ITHACAstream::ReadDenseMatrix(weights, "./ITHACAoutput/weightsPPE/",
                                          weightName);
            rbfSplines[i] = new SPLINTER::RBFSpline(*samples[i],
                                                    SPLINTER::RadialBasisFunctionType::GAUSSIAN, weights);
            std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
        }
        else
        {
            samples[i] = new SPLINTER::DataTable(1, 1);

            for (int j = 0; j < coeffL2.cols(); j++)
            {
                samples[i]->addSample(mu.row(j), coeffL2(i, j));
            }

            rbfSplines[i] = new SPLINTER::RBFSpline(*samples[i],
                                                    SPLINTER::RadialBasisFunctionType::GAUSSIAN);
            ITHACAstream::SaveDenseMatrix(rbfSplines[i]->weights,
                                          "./ITHACAoutput/weightsPPE/", weightName);
            std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
        }
    }
}
