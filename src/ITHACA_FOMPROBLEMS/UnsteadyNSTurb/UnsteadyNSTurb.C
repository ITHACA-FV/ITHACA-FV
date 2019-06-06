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
UnsteadyNSTurb::UnsteadyNSTurb() {};

// Construct from zero
UnsteadyNSTurb::UnsteadyNSTurb(int argc, char* argv[])
{
#include "setRootCase.H"
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
    supex = ITHACAutilities::check_sup();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void UnsteadyNSTurb::truthSolve(List<scalar> mu_now)
{
#include "initContinuityErrs.H"
    Time& runTime = _runTime();
    surfaceScalarField& phi = _phi();
    fvMesh& mesh = _mesh();
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

    PtrList<volVectorField> together(0);

    // Create PTRLIST with lift, velocities and supremizers

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            together.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            together.append(Umodes[k]);
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            together.append(supmodes[k]);
        }
    }

    for (label i = 0; i < cSize; i++)
    {
        Info << "Filling layer number " << i + 1 << " in the matrix ct1Matrix" << endl;

        for (label j = 0; j < nNutModes; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                ct1Matrix[i](j, k) = fvc::domainIntegrate(together[i] & fvc::laplacian(
                                         nutModes[j], together[k])).value();
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

    PtrList<volVectorField> together(0);

    // Create PTRLIST with lift, velocities and supremizers

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            together.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            together.append(Umodes[k]);
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            together.append(supmodes[k]);
        }
    }

    for (label i = 0; i < cSize; i++)
    {
        Info << "Filling layer number " << i + 1 << " in the matrix ct2Matrix" << endl;

        for (label j = 0; j < nNutModes; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                ct2Matrix[i](j, k) = fvc::domainIntegrate(together[i] & (fvc::div(
                                         nutModes[j] * dev((fvc::grad(together[k]))().T())))).value();
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

Eigen::MatrixXd UnsteadyNSTurb::btTurbulence(label NUmodes, label NSUPmodes)
{
    label btsize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd btMatrix(btsize, btsize);
    btMatrix = btMatrix * 0;
    // Create PTRLIST with lift, velocities and supremizers
    PtrList<volVectorField> together(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            together.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            together.append(Umodes[k]);
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            together.append(supmodes[k]);
        }
    }

    // Project everything
    for (label i = 0; i < btsize; i++)
    {
        for (label j = 0; j < btsize; j++)
        {
            btMatrix(i, j) = fvc::domainIntegrate(together[i] & (fvc::div(dev((T(fvc::grad(
                    together[j]))))))).value();
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(btMatrix, "btMatrix", "python",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(btMatrix, "btMatrix", "matlab",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(btMatrix, "btMatrix", "eigen",
                               "./ITHACAoutput/Matrices/");
    return btMatrix;
}

void UnsteadyNSTurb::projectSUP(fileName folder, label NU, label NP, label NSUP,
                                label Nnut)
{
    NUmodes = NU;
    NPmodes = NP;
    NSUPmodes = NSUP;
    nNutModes = Nnut;

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

        C_matrix = convective_term(NUmodes, NPmodes, NSUPmodes);
        ct1Matrix = turbulenceTerm1(NUmodes, NSUPmodes, nNutModes);
        ct2Matrix = turbulenceTerm2(NUmodes, NSUPmodes, nNutModes);

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }
    }
    else
    {
        B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        C_matrix = convective_term(NUmodes, NPmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        P_matrix = divergence_term(NUmodes, NPmodes, NSUPmodes);
        M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
        btMatrix = btTurbulence(NUmodes, NSUPmodes);
        ct1Matrix = turbulenceTerm1(NUmodes, NSUPmodes, nNutModes);
        ct2Matrix = turbulenceTerm2(NUmodes, NSUPmodes, nNutModes);

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
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(P_matrix, "P", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "matlab", "./ITHACAoutput/Matrices/");
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(P_matrix, "P", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "eigen", "./ITHACAoutput/Matrices/");
    }

    bTotalMatrix = B_matrix + btMatrix;
    cTotalMatrix.setSize(C_matrix.size());

    for (label i = 0; i < C_matrix.size(); i++)
    {
        cTotalMatrix[i] =  ct2Matrix[i] + ct1Matrix[i];
    }

    // Get the coeffs for interpolation (the orthonormal one is used because basis are orthogonal)
    Eigen::MatrixXd coeffL2 = ITHACAutilities::get_coeffs_ortho(nutFields,
                              nutModes);
    ITHACAstream::exportMatrix(coeffL2, "coeffL2", "python",
                               "./ITHACAoutput/Matrices/");
    NUmodes = NU;
    NPmodes = NP;
    NSUPmodes = NSUP;
    nNutModes = Nnut;
    samples.resize(nNutModes);
    rbfSplines.resize(nNutModes);

    for (int i = 0; i < nNutModes; i++)
    {
        samples[i] = new SPLINTER::DataTable(1, 1);

        for (int j = 0; j < coeffL2.cols(); j++)
        {
            samples[i]->addSample(mu.row(j), coeffL2(i, j));
        }

        rbfSplines[i] = new SPLINTER::RBFSpline(*samples[i],
                                                SPLINTER::RadialBasisFunctionType::GAUSSIAN);
        std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
    }
}


void UnsteadyNSTurb::projectPPE(fileName folder, label NU, label NP, label NSUP,
                                label Nnut)
{
    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {
        B_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/B_mat.txt");
        C_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/C", "C");
        K_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/K_mat.txt");
        M_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/M_mat.txt");
        D_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/D_mat.txt");
        G_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/G", "G");
        BC1_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/BC1_mat.txt");
        BC2_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/BC2/", "BC2");
        BC3_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/BC3_mat.txt");
        btMatrix =
            ITHACAstream::readMatrix("./ITHACAoutput/Matrices/btMatrix_mat.txt");
        ct1Matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/ct1",
                                             "ct1Matrix");
        ct2Matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/ct2",
                                             "ct2Matrix");
    }
    else
    {
        NUmodes = NU;
        NPmodes = NP;
        NSUPmodes = 0;
        nNutModes = Nnut;
        B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        C_matrix = convective_term(NUmodes, NPmodes, NSUPmodes);
        M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        D_matrix = laplacian_pressure(NPmodes);
        G_matrix = div_momentum(NUmodes, NPmodes);
        BC1_matrix = pressure_BC1(NUmodes, NPmodes);
        BC2_matrix = pressure_BC2(NUmodes, NPmodes);
        BC3_matrix = pressure_BC3(NUmodes, NPmodes);
        btMatrix = btTurbulence(NUmodes, NSUPmodes);
        ct1Matrix = turbulenceTerm1(NUmodes, NSUPmodes, nNutModes);
        ct2Matrix = turbulenceTerm2(NUmodes, NSUPmodes, nNutModes);
    }

    bTotalMatrix = B_matrix + btMatrix;
    cTotalMatrix.setSize(C_matrix.size());

    for (label i = 0; i < C_matrix.size(); i++)
    {
        cTotalMatrix[i] =  ct2Matrix[i] + ct1Matrix[i];
    }

    // Get the coeffs for interpolation (the orthonormal one is used because basis are orthogonal)
    Eigen::MatrixXd coeffL2 = ITHACAutilities::get_coeffs_ortho(nutFields,
                              nutModes);
    ITHACAstream::exportMatrix(coeffL2, "coeffL2", "python",
                               "./ITHACAoutput/Matrices/");
    NUmodes = NU;
    NPmodes = NP;
    NSUPmodes = 0;
    nNutModes = Nnut;
    samples.resize(nNutModes);
    rbfSplines.resize(nNutModes);

    for (int i = 0; i < nNutModes; i++)
    {
        samples[i] = new SPLINTER::DataTable(1, 1);

        for (int j = 0; j < coeffL2.cols(); j++)
        {
            samples[i]->addSample(mu.row(j), coeffL2(i, j));
        }

        rbfSplines[i] = new SPLINTER::RBFSpline(*samples[i],
                                                SPLINTER::RadialBasisFunctionType::GAUSSIAN);
        std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
    }

    for (int i = 0; i < nNutModes; i++)
    {
        for (int j = 0; j < coeffL2.cols(); j++)
        {
            Info << rbfSplines[i]->eval(mu.row(j)) - coeffL2(i, j) << endl;
        }
    }
}
