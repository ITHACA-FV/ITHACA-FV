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

// Small construct to access member functions linked to Smagorinsky diffusion   
UnsteadyNSTurb::UnsteadyNSTurb(const ITHACAPOD::Parameters* myParameters):
  ithacaFVParameters(static_cast<const ITHACAPOD::PODParameters*>(myParameters))
{
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

    // Export the tensor
    ITHACAstream::SaveDenseTensor(ct1Tensor, "./ITHACAoutput/Matrices/",
                                  "ct1_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(nNutModes) + "_t");
    return ct1Tensor;
}


Eigen::Tensor<double, 3> UnsteadyNSTurb::turbulenceAveTensor1(label NUmodes,
        label NSUPmodes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct1AveTensor;
    label samplesNumber = nutAve.size();
    ct1AveTensor.resize(cSize, samplesNumber, cSize);

    for (label i = 0; i < cSize; i++)
    {
        for (label j = 0; j < samplesNumber; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                ct1AveTensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::laplacian(
                                            nutAve[j], L_U_SUPmodes[k])).value();
            }
        }
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(ct1AveTensor, "./ITHACAoutput/Matrices/",
                                  "ct1Ave_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_t");
    return ct1AveTensor;
}

Eigen::Tensor<double, 3> UnsteadyNSTurb::turbulencePPETensor1(label NUmodes,
        label NSUPmodes, label NPmodes, label nNutModes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct1PPETensor;
    ct1PPETensor.resize(NPmodes, nNutModes, cSize);

    for (label i = 0; i < NPmodes; i++)
    {
        for (label j = 0; j < nNutModes; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                // ct1PPETensor(i, j, k) = fvc::domainIntegrate(2 * Pmodes[i] * (fvc::laplacian(
                //                             L_U_SUPmodes[k]) & fvc::grad(nutModes[j]))).value();
                // ct1PPETensor(i, j, k) = fvc::domainIntegrate(Pmodes[i] * (fvc::div(
                //     fvc::laplacian(
                //         nutModes[j], L_U_SUPmodes[k])))).value();
                ct1PPETensor(i, j, k) = fvc::domainIntegrate(fvc::grad(Pmodes[i]) & (
                                            fvc::laplacian(
                                                nutModes[j], L_U_SUPmodes[k]))).value();
            }
        }
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(ct1PPETensor, "./ITHACAoutput/Matrices/",
                                  "ct1PPE_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(NPmodes) + "_" + name(nNutModes) + "_t");
    return ct1PPETensor;
}

Eigen::Tensor<double, 3> UnsteadyNSTurb::turbulencePPEAveTensor1(label NUmodes,
        label NSUPmodes, label NPmodes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct1PPEAveTensor;
    label samplesNumber = nutAve.size();
    ct1PPEAveTensor.resize(NPmodes, samplesNumber, cSize);

    for (label i = 0; i < NPmodes; i++)
    {
        for (label j = 0; j < samplesNumber; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                // ct1PPEAveTensor(i, j, k) = fvc::domainIntegrate(2 * Pmodes[i] * (fvc::laplacian(
                //                                L_U_SUPmodes[k]) & fvc::grad(nutAve[j]))).value();
                // ct1PPEAveTensor(i, j, k) = fvc::domainIntegrate(Pmodes[i] * (fvc::div(
                //  fvc::laplacian(
                //      nutAve[j], L_U_SUPmodes[k])))).value();
                ct1PPEAveTensor(i, j, k) = fvc::domainIntegrate(fvc::grad(Pmodes[i]) & (
                                               fvc::laplacian(
                                                   nutAve[j], L_U_SUPmodes[k]))).value();
            }
        }
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(ct1PPEAveTensor, "./ITHACAoutput/Matrices/",
                                  "ct1PPEAve_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(NPmodes) + "_t");
    return ct1PPEAveTensor;
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
                                         nutModes[j] * dev2((fvc::grad(L_U_SUPmodes[k]))().T())))).value();
            }
        }
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(ct2Tensor, "./ITHACAoutput/Matrices/",
                                  "ct2_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(nNutModes) + "_t");
    return ct2Tensor;
}

Eigen::Tensor<double, 3> UnsteadyNSTurb::turbulenceAveTensor2(label NUmodes,
        label NSUPmodes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct2AveTensor;
    label samplesNumber = nutAve.size();
    ct2AveTensor.resize(cSize, samplesNumber, cSize);

    for (label i = 0; i < cSize; i++)
    {
        for (label j = 0; j < samplesNumber; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                ct2AveTensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & (fvc::div(
                                            nutAve[j] * dev2((fvc::grad(L_U_SUPmodes[k]))().T())))).value();
            }
        }
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(ct2AveTensor, "./ITHACAoutput/Matrices/",
                                  "ct2Ave_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_t");
    return ct2AveTensor;
}

Eigen::Tensor<double, 3> UnsteadyNSTurb::turbulencePPETensor2(label NUmodes,
        label NSUPmodes, label NPmodes, label nNutModes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct2PPETensor;
    ct2PPETensor.resize(NPmodes, nNutModes, cSize);

    for (label i = 0; i < NPmodes; i++)
    {
        for (label j = 0; j < nNutModes; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                // ct2PPETensor(i, j, k) = fvc::domainIntegrate(Pmodes[i] * (fvc::grad(fvc::grad(
                //                             nutModes[j])) && (dev2((fvc::grad(L_U_SUPmodes[k]))() + fvc::grad(
                //                                         L_U_SUPmodes[k]))().T()))).value();
                // ct2PPETensor(i, j, k) = fvc::domainIntegrate(Pmodes[i] * ((fvc::div(fvc::div(
                //     nutModes[j] * dev2((fvc::grad(L_U_SUPmodes[k]))().T())))))).value();
                ct2PPETensor(i, j, k) = fvc::domainIntegrate(fvc::grad(Pmodes[i]) & ((fvc::div(
                                            nutModes[j] * dev2((fvc::grad(L_U_SUPmodes[k]))().T()))))).value();
            }
        }
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(ct2PPETensor, "./ITHACAoutput/Matrices/",
                                  "ct2PPE_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(NPmodes) + "_" + name(nNutModes) + "_t");
    return ct2PPETensor;
}

Eigen::Tensor<double, 3> UnsteadyNSTurb::turbulencePPEAveTensor2(label NUmodes,
        label NSUPmodes, label NPmodes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct2PPEAveTensor;
    label samplesNumber = nutAve.size();
    ct2PPEAveTensor.resize(NPmodes, samplesNumber, cSize);

    for (label i = 0; i < NPmodes; i++)
    {
        for (label j = 0; j < samplesNumber; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                // ct2PPEAveTensor(i, j, k) = fvc::domainIntegrate(Pmodes[i] * (fvc::grad(
                //                                fvc::grad(
                //                                    nutAve[j])) && (dev2((fvc::grad(L_U_SUPmodes[k]))() + fvc::grad(
                //                                            L_U_SUPmodes[k]))().T()))).value();
                // ct2PPEAveTensor(i, j, k) = fvc::domainIntegrate(Pmodes[i] * ((fvc::div(fvc::div(
                //  nutAve[j] * dev2((fvc::grad(L_U_SUPmodes[k]))().T())))))).value();
                ct2PPEAveTensor(i, j, k) = fvc::domainIntegrate(fvc::grad(Pmodes[i]) & ((
                                               fvc::div(
                                                   nutAve[j] * dev2((fvc::grad(L_U_SUPmodes[k]))().T()))))).value();
            }
        }
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(ct2PPEAveTensor, "./ITHACAoutput/Matrices/",
                                  "ct2PPEAve_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(NPmodes) + "_t");
    return ct2PPEAveTensor;
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
            btMatrix(i, j) = fvc::domainIntegrate(L_U_SUPmodes[i] & (fvc::div(dev2((T(
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
                                label Nnut, bool rbfInterp)
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

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            L_U_SUPmodes.append(tmp<volVectorField>(supmodes[k]));
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

        if (nutAve.size() != 0)
        {
            word ct1AveStr = "ct1Ave_" + name(liftfield.size()) + "_" + name(
                                 NUmodes) + "_" + name(
                                 NSUPmodes) + "_t";

            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct1AveStr))
            {
                ITHACAstream::ReadDenseTensor(ct1AveTensor, "./ITHACAoutput/Matrices/",
                                              ct1AveStr);
            }
            else
            {
                ct1AveTensor = turbulenceAveTensor1(NUmodes, NSUPmodes);
            }

            word ct2AveStr = "ct2Ave_" + name(liftfield.size()) + "_" + name(
                                 NUmodes) + "_" + name(
                                 NSUPmodes) + "_t";

            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct2AveStr))
            {
                ITHACAstream::ReadDenseTensor(ct2AveTensor, "./ITHACAoutput/Matrices/",
                                              ct2AveStr);
            }
            else
            {
                ct2AveTensor = turbulenceAveTensor2(NUmodes, NSUPmodes);
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
        B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        P_matrix = divergence_term(NUmodes, NPmodes, NSUPmodes);
        M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
        btMatrix = btTurbulence(NUmodes, NSUPmodes);
        ct1Tensor = turbulenceTensor1(NUmodes, NSUPmodes, nNutModes);
        ct2Tensor = turbulenceTensor2(NUmodes, NSUPmodes, nNutModes);

        if (nutAve.size() != 0)
        {
            ct1AveTensor = turbulenceAveTensor1(NUmodes, NSUPmodes);
            ct2AveTensor = turbulenceAveTensor2(NUmodes, NSUPmodes);
        }

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

        if (nutAve.size() != 0)
        {
            ITHACAstream::exportTensor(ct1AveTensor, "ct1Ave", "python",
                                       "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2AveTensor, "ct2Ave", "python",
                                       "./ITHACAoutput/Matrices/");
        }
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

        if (nutAve.size() != 0)
        {
            ITHACAstream::exportTensor(ct1AveTensor, "ct1Ave", "matlab",
                                       "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2AveTensor, "ct2Ave", "matlab",
                                       "./ITHACAoutput/Matrices/");
        }
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

        if (nutAve.size() != 0)
        {
            ITHACAstream::exportTensor(ct1AveTensor, "ct1Ave_", "eigen",
                                       "./ITHACAoutput/Matrices/ct1Ave");
            ITHACAstream::exportTensor(ct2AveTensor, "ct2Ave_", "eigen",
                                       "./ITHACAoutput/Matrices/ct2Ave");
        }
    }

    bTotalMatrix = B_matrix + btMatrix;
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    cTotalTensor.resize(cSize, nNutModes, cSize);
    cTotalTensor = ct1Tensor + ct2Tensor;
    // Define coeffL2
    coeffL2 = ITHACAutilities::getCoeffs(nutFields,
                                         nutModes, nNutModes);
    ITHACAstream::exportMatrix(coeffL2, "coeffL2", "python",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(coeffL2, "coeffL2", "matlab",
                               "./ITHACAoutput/Matrices/");

    if (nutAve.size() != 0)
    {
        cTotalAveTensor.resize(cSize, nutAve.size(), cSize);
        cTotalAveTensor = ct1AveTensor + ct2AveTensor;
    }

    if (rbfInterp == true && (!Pstream::parRun()))
    {
        if (ITHACAutilities::check_file("./radii.txt"))
        {
            radii = ITHACAstream::readMatrix("./radii.txt");
            M_Assert(radii.size() == nNutModes,
                     "Thes size of the shape parameters vector must be equal to the number of eddy viscosity modes nNutModes");
        }
        else
        {
            radii = Eigen::MatrixXd::Ones(nNutModes,
                                          1) * e;
        }

        samples.resize(nNutModes);
        rbfSplines.resize(nNutModes);
        Eigen::MatrixXd weights;

        for (label i = 0; i < nNutModes; i++)
        {
            word weightName = "wRBF_N" + name(i + 1) + "_" + name(liftfield.size()) + "_"
                              + name(NUmodes) + "_" + name(NSUPmodes) ;

            if (ITHACAutilities::check_file("./ITHACAoutput/weightsSUP/" + weightName))
            {
                samples[i] = new SPLINTER::DataTable(1, 1);

                for (label j = 0; j < coeffL2.cols(); j++)
                {
                    samples[i]->addSample(velRBF.row(j), coeffL2(i, j));
                }

                ITHACAstream::ReadDenseMatrix(weights, "./ITHACAoutput/weightsSUP/",
                                              weightName);
                rbfSplines[i] = new SPLINTER::RBFSpline(* samples[i],
                                                        SPLINTER::RadialBasisFunctionType::GAUSSIAN, weights, radii(i));
                std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
            }
            else
            {
                samples[i] = new SPLINTER::DataTable(1, 1);

                for (label j = 0; j < coeffL2.cols(); j++)
                {
                    samples[i]->addSample(velRBF.row(j), coeffL2(i, j));
                }

                rbfSplines[i] = new SPLINTER::RBFSpline(* samples[i],
                                                        SPLINTER::RadialBasisFunctionType::GAUSSIAN, false, radii(i));
                ITHACAstream::SaveDenseMatrix(rbfSplines[i]->weights,
                                              "./ITHACAoutput/weightsSUP/", weightName);
                std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
            }
        }
    }
}



void UnsteadyNSTurb::projectPPE(fileName folder, label NU, label NP, label NSUP,
                                label Nnut, bool rbfInterp)
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

        if (nutAve.size() != 0)
        {
            word ct1AveStr = "ct1Ave_" + name(liftfield.size()) + "_" + name(
                                 NUmodes) + "_" + name(
                                 NSUPmodes) + "_t";

            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct1AveStr))
            {
                ITHACAstream::ReadDenseTensor(ct1AveTensor, "./ITHACAoutput/Matrices/",
                                              ct1AveStr);
            }
            else
            {
                ct1AveTensor = turbulenceAveTensor1(NUmodes, NSUPmodes);
            }

            word ct2AveStr = "ct2Ave_" + name(liftfield.size()) + "_" + name(
                                 NUmodes) + "_" + name(
                                 NSUPmodes) + "_t";

            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct2AveStr))
            {
                ITHACAstream::ReadDenseTensor(ct2AveTensor, "./ITHACAoutput/Matrices/",
                                              ct2AveStr);
            }
            else
            {
                ct2AveTensor = turbulenceAveTensor2(NUmodes, NSUPmodes);
            }
        }

        word ct1PPEStr = "ct1PPE_" + name(liftfield.size()) + "_" + name(
                             NUmodes) + "_" + name(
                             NSUPmodes) + "_" + name(NPmodes) + "_" + name(nNutModes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct1PPEStr))
        {
            ITHACAstream::ReadDenseTensor(ct1PPETensor, "./ITHACAoutput/Matrices/",
                                          ct1PPEStr);
        }
        else
        {
            ct1PPETensor = turbulencePPETensor1(NUmodes, NSUPmodes, NPmodes, nNutModes);
        }

        word ct2PPEStr = "ct2PPE_" + name(liftfield.size()) + "_" + name(
                             NUmodes) + "_" + name(
                             NSUPmodes) + "_" + name(NPmodes) + "_" + name(nNutModes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct2PPEStr))
        {
            ITHACAstream::ReadDenseTensor(ct2PPETensor, "./ITHACAoutput/Matrices/",
                                          ct2PPEStr);
        }
        else
        {
            ct2PPETensor = turbulencePPETensor2(NUmodes, NSUPmodes, NPmodes, nNutModes);
        }

        if (nutAve.size() != 0)
        {
            word ct1PPEAveStr = "ct1PPEAve_" + name(liftfield.size()) + "_" + name(
                                    NUmodes) + "_" + name(
                                    NSUPmodes) + "_" + name(NPmodes) + "_t";

            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct1PPEAveStr))
            {
                ITHACAstream::ReadDenseTensor(ct1PPEAveTensor, "./ITHACAoutput/Matrices/",
                                              ct1PPEAveStr);
            }
            else
            {
                ct1PPEAveTensor = turbulencePPEAveTensor1(NUmodes, NSUPmodes, NPmodes);
            }

            word ct2PPEAveStr = "ct2PPEAve_" + name(liftfield.size()) + "_" + name(
                                    NUmodes) + "_" + name(
                                    NSUPmodes) + "_" + name(NPmodes) + "_t";

            if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct2PPEAveStr))
            {
                ITHACAstream::ReadDenseTensor(ct2PPEAveTensor, "./ITHACAoutput/Matrices/",
                                              ct2PPEAveStr);
            }
            else
            {
                ct2PPEAveTensor = turbulencePPEAveTensor2(NUmodes, NSUPmodes, NPmodes);
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
        ct1PPETensor = turbulencePPETensor1(NUmodes, NSUPmodes, NPmodes, nNutModes);
        ct2PPETensor = turbulencePPETensor2(NUmodes, NSUPmodes, NPmodes, nNutModes);

        if (nutAve.size() != 0)
        {
            ct1AveTensor = turbulenceAveTensor1(NUmodes, NSUPmodes);
            ct2AveTensor = turbulenceAveTensor2(NUmodes, NSUPmodes);
            ct1PPEAveTensor = turbulencePPEAveTensor1(NUmodes, NSUPmodes, NPmodes);
            ct2PPEAveTensor = turbulencePPEAveTensor2(NUmodes, NSUPmodes, NPmodes);
        }

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
        ITHACAstream::exportTensor(ct1PPETensor, "ct1PPE", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct2PPETensor, "ct2PPE", "python",
                                   "./ITHACAoutput/Matrices/");

        if (nutAve.size() != 0)
        {
            ITHACAstream::exportTensor(ct1AveTensor, "ct1Ave", "python",
                                       "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2AveTensor, "ct2Ave", "python",
                                       "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct1PPEAveTensor, "ct1PPEAve", "python",
                                       "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2PPEAveTensor, "ct2PPEAve", "python",
                                       "./ITHACAoutput/Matrices/");
        }
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
        ITHACAstream::exportTensor(ct1PPETensor, "ct1PPE", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ct2PPETensor, "ct2PPE", "matlab",
                                   "./ITHACAoutput/Matrices/");

        if (nutAve.size() != 0)
        {
            ITHACAstream::exportTensor(ct1AveTensor, "ct1Ave", "matlab",
                                       "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2AveTensor, "ct2Ave", "matlab",
                                       "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct1PPEAveTensor, "ct1PPEAve", "matlab",
                                       "./ITHACAoutput/Matrices/");
            ITHACAstream::exportTensor(ct2PPEAveTensor, "ct2PPEAve", "matlab",
                                       "./ITHACAoutput/Matrices/");
        }
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
        ITHACAstream::exportTensor(ct1PPETensor, "ct1PPE_", "eigen",
                                   "./ITHACAoutput/Matrices/ct1PPE");
        ITHACAstream::exportTensor(ct2PPETensor, "ct2PPE_", "eigen",
                                   "./ITHACAoutput/Matrices/ct2PPE");

        if (nutAve.size() != 0)
        {
            ITHACAstream::exportTensor(ct1AveTensor, "ct1Ave_", "eigen",
                                       "./ITHACAoutput/Matrices/ct1Ave");
            ITHACAstream::exportTensor(ct2AveTensor, "ct2Ave_", "eigen",
                                       "./ITHACAoutput/Matrices/ct2Ave");
            ITHACAstream::exportTensor(ct1PPEAveTensor, "ct1PPEAve_", "eigen",
                                       "./ITHACAoutput/Matrices/ct1PPEAve");
            ITHACAstream::exportTensor(ct2PPEAveTensor, "ct2PPEAve_", "eigen",
                                       "./ITHACAoutput/Matrices/ct2PPEAve");
        }
    }

    bTotalMatrix = B_matrix + btMatrix;
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    cTotalTensor.resize(cSize, nNutModes, cSize);
    cTotalTensor = ct1Tensor + ct2Tensor;
    cTotalPPETensor.resize(NPmodes, nNutModes, cSize);
    cTotalPPETensor = ct1PPETensor + ct2PPETensor;

    if (nutAve.size() != 0)
    {
        cTotalAveTensor.resize(cSize, nutAve.size(), cSize);
        cTotalAveTensor = ct1AveTensor + ct2AveTensor;
        cTotalPPEAveTensor.resize(NPmodes, nutAve.size(), cSize);
        cTotalPPEAveTensor = ct1PPEAveTensor + ct2PPEAveTensor;
    }

    if (rbfInterp == true && (!Pstream::parRun()))
    {
        if (ITHACAutilities::check_file("./radii.txt"))
        {
            radii = ITHACAstream::readMatrix("./radii.txt");
            M_Assert(radii.size() ==  nNutModes,
                     "Thes size of the shape parameters vector must be equal to the number of eddy viscosity modes nNutModes");
        }
        else
        {
            radii = Eigen::MatrixXd::Ones(nNutModes,
                                          1) * e;
        }

        samples.resize(nNutModes);
        rbfSplines.resize(nNutModes);
        Eigen::MatrixXd weights;

        for (label i = 0; i < nNutModes; i++)
        {
            word weightName = "wRBF_N" + name(i + 1) + "_" + name(liftfield.size()) + "_"
                              + name(NUmodes) + "_" + name(NSUPmodes) ;

            if (ITHACAutilities::check_file("./ITHACAoutput/weightsPPE/" + weightName))
            {
                samples[i] = new SPLINTER::DataTable(1, 1);

                for (label j = 0; j < coeffL2.cols(); j++)
                {
                    samples[i]->addSample(velRBF.row(j), coeffL2(i, j));
                }

                ITHACAstream::ReadDenseMatrix(weights, "./ITHACAoutput/weightsPPE/",
                                              weightName);
                rbfSplines[i] = new SPLINTER::RBFSpline(* samples[i],
                                                        SPLINTER::RadialBasisFunctionType::GAUSSIAN, weights, radii(i));
                std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
            }
            else
            {
                samples[i] = new SPLINTER::DataTable(1, 1);

                for (label j = 0; j < coeffL2.cols(); j++)
                {
                    samples[i]->addSample(velRBF.row(j), coeffL2(i, j));
                }

                rbfSplines[i] = new SPLINTER::RBFSpline(* samples[i],
                                                        SPLINTER::RadialBasisFunctionType::GAUSSIAN, false, radii(i));
                ITHACAstream::SaveDenseMatrix(rbfSplines[i]->weights,
                                              "./ITHACAoutput/weightsPPE/", weightName);
                std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
            }
        }
    }
}

List < Eigen::MatrixXd > UnsteadyNSTurb::velDerivativeCoeff(Eigen::MatrixXd A,
        Eigen::MatrixXd G,
        Eigen::VectorXd initSnapInd, Eigen::VectorXd timeSnap)
{
    List < Eigen::MatrixXd > newCoeffs;
    newCoeffs.setSize(2);
    label velCoeffsNum = A.cols();
    label snapshotsNum = A.rows();
    Eigen::MatrixXd pars;
    label parsSamplesNum = initSnapInd.size();
    label timeSnapshotsPerSample = snapshotsNum / parsSamplesNum;
    label newColsNum = 2 * velCoeffsNum;
    label newRowsNum = snapshotsNum - parsSamplesNum;
    newCoeffs[0].resize(newRowsNum, newColsNum);
    newCoeffs[1].resize(newRowsNum, G.cols());

    for (label j = 0; j < parsSamplesNum; j++)
    {
        Eigen::MatrixXd b0 = A.middleRows(j * timeSnapshotsPerSample,
                                          timeSnapshotsPerSample - 1);
        Eigen::MatrixXd b2 = A.middleRows(j * timeSnapshotsPerSample + 1,
                                          timeSnapshotsPerSample - 1);
        Eigen::MatrixXd bNew(b0.rows(), b0.cols() + b2.cols());
        bNew << b2, (b2 - b0) / (timeSnap(j, 0));
        newCoeffs[0].block(j * timeSnapshotsPerSample - j, 0,
                           timeSnapshotsPerSample - 1, newColsNum) = bNew;
        newCoeffs[1].middleRows(j * timeSnapshotsPerSample - j,
                                timeSnapshotsPerSample - 1) = G.middleRows(j * timeSnapshotsPerSample + 1,
                                    timeSnapshotsPerSample - 1);
    }

    interChoice = 3;
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

List < Eigen::MatrixXd > UnsteadyNSTurb::velParDerivativeCoeff(
    Eigen::MatrixXd A, Eigen::MatrixXd G,
    Eigen::VectorXd initSnapInd, Eigen::VectorXd timeSnap)
{
    List < Eigen::MatrixXd > newCoeffs;
    newCoeffs.setSize(2);
    label velCoeffsNum = A.cols();
    label snapshotsNum = A.rows();
    Eigen::MatrixXd pars;
    pars = z.leftCols(z.cols() - 1);
    label parsSamplesNum = initSnapInd.size();
    label timeSnapshotsPerSample = snapshotsNum / parsSamplesNum;
    label newColsNum = 2 * velCoeffsNum;
    label newRowsNum = snapshotsNum - parsSamplesNum;
    newCoeffs[0].resize(newRowsNum, newColsNum + z.cols() - 1);
    newCoeffs[1].resize(newRowsNum, G.cols());

    for (label j = 0; j < parsSamplesNum; j++)
    {
        Eigen::MatrixXd b0 = A.middleRows(j * timeSnapshotsPerSample,
                                          timeSnapshotsPerSample - 1);
        Eigen::MatrixXd b2 = A.middleRows(j * timeSnapshotsPerSample + 1,
                                          timeSnapshotsPerSample - 1);
        Eigen::MatrixXd bNew(b0.rows(), b0.cols() + b2.cols());
        bNew << b2, (b2 - b0) / (timeSnap(j, 0));
        newCoeffs[0].block(j * timeSnapshotsPerSample - j, 0,
                           timeSnapshotsPerSample - 1, z.cols() - 1) = pars.middleRows(
                                   j * timeSnapshotsPerSample + 1,
                                   timeSnapshotsPerSample - 1);
        newCoeffs[0].block(j * timeSnapshotsPerSample - j, z.cols() - 1,
                           timeSnapshotsPerSample - 1, newColsNum) = bNew;
        newCoeffs[1].middleRows(j * timeSnapshotsPerSample - j,
                                timeSnapshotsPerSample - 1) = G.middleRows(j * timeSnapshotsPerSample + 1,
                                    timeSnapshotsPerSample - 1);
    }

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



// =========================================================================
// ============================= Smagorinsky ===============================
// =========================================================================

void UnsteadyNSTurb::computeNonLinearSnapshot_at_time(const word& snap_time, 
                 volVectorField& Smag, std::optional<PtrList<volVectorField>> modesU)
{
    if (ithacaFVParameters->get_DEIMInterpolatedField() == "fullStressFunction"
      || ITHACAutilities::containsSubstring(ithacaFVParameters->get_DEIMInterpolatedField(), "reducedFullStressFunction"))
    {
        Smag = computeSmagTerm_at_time(snap_time, modesU) ;
    }
    else
    {
      Info << "Error : DEIMInterpolatedField not valid : "
           << ithacaFVParameters->get_DEIMInterpolatedField() << endl;
      Info << "DEIM is available for fullStressFunction and nut only." << endl;
      abort();
    }
}

// Init Smagorinsky term
volVectorField UnsteadyNSTurb::initSmagFunction()
{
    return volVectorField("fullStressFunction", computeSmagTerm_at_time("0") );
}

// Init Smagorinsky term
volScalarField UnsteadyNSTurb::initSmagPhiFunction(const volScalarField template_field_phi)
{
    return volScalarField(ithacaFVParameters->get_DEIMInterpolatedField(), computeSmagTermPhi_at_time("0",template_field_phi));
}

volVectorField UnsteadyNSTurb::computeSmagTerm_at_time(const word& snap_time, 
                                        std::optional<PtrList<volVectorField>> modesU)
{
    // Read the j-th field
    volVectorField snapshotj = ithacaFVParameters->get_template_field_U();
    ITHACAstream::read_snapshot(snapshotj, snap_time, ithacaFVParameters->get_casenameData());

    if (modesU)
    {
        Eigen::MatrixXd coeffs_field_U = ITHACAutilities::getCoeffs(snapshotj, modesU.value(), modesU.value().size(), 1);
        PtrList<volVectorField> reduced_field_U = ITHACAutilities::reconstructFromCoeff(modesU.value(), 
                                                                coeffs_field_U, modesU.value().size());
        snapshotj = ithacaFVParameters->get_meanU() + reduced_field_U[0];
    }

    return computeSmagTerm_fromU(snapshotj);
}

volScalarField UnsteadyNSTurb::computeSmagTermPhi_at_time(const word& snap_time, const volScalarField template_field_phi)
{
    // Read the j-th field
    volScalarField phij = template_field_phi;
    ITHACAstream::read_snapshot(phij, snap_time, ithacaFVParameters->get_casenameData());
    volVectorField snapshotj = ithacaFVParameters->get_template_field_U();
    ITHACAstream::read_snapshot(snapshotj, snap_time, ithacaFVParameters->get_casenameData());
    return computeSmagTermPhi_fromUPhi(snapshotj,phij);
}

volVectorField UnsteadyNSTurb::computeSmagTerm_fromU(const volVectorField& snapshotj)
{
    volTensorField S=computeS_fromU(snapshotj);
    volScalarField nut = computeNut_fromS(S);
    return (fvc::div(2*nut*dev(S)));
}

volScalarField UnsteadyNSTurb::computeSmagTermPhi_fromUPhi(const volVectorField& snapshotj, const volScalarField& phij)
{
    volTensorField S=computeS_fromU(snapshotj);
    volScalarField nut = computeNut_fromS(S);
    return (fvc::div(nut*fvc::grad(phij)));
}

template<typename T>
volScalarField UnsteadyNSTurb::diffusion(const T& coefDiff, const volScalarField& phi)
{
    return ( fvc::laplacian(coefDiff , phi) );
}
template volScalarField UnsteadyNSTurb::diffusion(const volScalarField& coefDiff, const volScalarField& u);
template volScalarField UnsteadyNSTurb::diffusion(const volTensorField& coefDiff, const volScalarField& u);

template<typename T>
volVectorField UnsteadyNSTurb::diffusion(const T& coefDiff, const volVectorField& u)
{
    volTensorField S_u=computeS_fromU(u);
    return ( fvc::div( 2 * ITHACAutilities::tensorFieldProduct( coefDiff , dev(S_u) ) ) );
}
template volVectorField UnsteadyNSTurb::diffusion(const volScalarField& coefDiff, const volVectorField& u);
template volVectorField UnsteadyNSTurb::diffusion(const volTensorField& coefDiff, const volVectorField& u);

volVectorField UnsteadyNSTurb::computeSmagTermPhi_fromUPhi(const volVectorField& snapshotj, const volVectorField& phij)
{
    return computeSmagTerm_fromU(snapshotj);
}


// =========================================================================
// ======================= projFullStressFunction ==========================
// =========================================================================

void UnsteadyNSTurb::computeNonLinearSnapshot_at_time(const word& snap_time, volScalarField& phi, volVectorField& modeU)
{
    if (ITHACAutilities::containsSubstring(ithacaFVParameters->get_HRSnapshotsField(), "projFullStressFunction") 
        || ITHACAutilities::containsSubstring(ithacaFVParameters->get_HRSnapshotsField(), "projReducedFullStressFunction"))
    {
        phi = computeProjSmagTerm_at_time_fromMode(snap_time, modeU);
    }
}

// Init projected Smagorinsky term
volScalarField UnsteadyNSTurb::initProjSmagFunction()
{
    return volScalarField("projFullStressFunction",
        ithacaFVParameters->get_template_field_fullStressFunction() & ithacaFVParameters->get_template_field_U());
}

volScalarField UnsteadyNSTurb::computeProjSmagTerm_at_time_fromMode(const word& snap_time, const volVectorField& mode)
{
    // Read the j-th field
    volVectorField Smagj = ithacaFVParameters->get_template_field_fullStressFunction();  
    if (!ITHACAutilities::containsSubstring(ithacaFVParameters->get_HRSnapshotsField(), "projReducedFullStressFunction"))
    {
        ITHACAstream::read_snapshot(Smagj, snap_time, ithacaFVParameters->get_casenameData());
    }
    else
    {
        word path = "./ITHACAoutput/Hyperreduction/reducedFullStressFunction/";
        ITHACAstream::read_snapshot(Smagj, snap_time, path);
    }
    
    return Smagj & mode;
}

void UnsteadyNSTurb::computeProjSmagTerm_fromSmag_fromMode(volScalarField& phi, const volVectorField& Smag, 
                                                                                 const volVectorField& mode)
{
    phi = Smag & mode;
}


// =========================================================================
// ========================== projSmagFromNut ==============================
// =========================================================================

void UnsteadyNSTurb::computeNonLinearSnapshot_at_time(const word& snap_time, volScalarField& phi, 
                                                                 volVectorField& modeU_proj , volVectorField& modeU_grad)
{
    if (ITHACAutilities::containsSubstring(ithacaFVParameters->get_HRSnapshotsField(), "projSmagFromNut") 
        || ITHACAutilities::containsSubstring(ithacaFVParameters->get_HRSnapshotsField(), "projSmagFromReducedNut"))
    {
        phi = computeProjSmagFromNut_at_time_fromModes(snap_time, modeU_proj, modeU_grad);
    }
}

// Init projected Smagorinsky term from nut
volScalarField UnsteadyNSTurb::initProjSmagFromNutFunction()
{
    return volScalarField("projSmagFromNut", initSmagFunction() & ithacaFVParameters->get_template_field_U());
}

volScalarField UnsteadyNSTurb::computeProjSmagFromNut_at_time_fromModes
                  (const word& snap_time, const volVectorField& modeU_proj, const volVectorField& modeU_grad)
{
    // Read the j-th field
    volScalarField Nutj(ithacaFVParameters->get_DEIMInterpolatedField(), ithacaFVParameters->get_template_field_nut());
    if (ITHACAutilities::containsSubstring(ithacaFVParameters->get_HRSnapshotsField(), "projSmagFromNut") )
    {
        ITHACAstream::read_snapshot(Nutj, snap_time, ithacaFVParameters->get_casenameData());
    }
    else if (ITHACAutilities::containsSubstring(ithacaFVParameters->get_HRSnapshotsField(), "projSmagFromReducedNut"))
    {
        word path = "./ITHACAoutput/Hyperreduction/reducedNut/";
        ITHACAstream::read_snapshot(Nutj, snap_time, path);
    }
    return projDiffusionIBP(Nutj, modeU_grad, modeU_proj);
}

void UnsteadyNSTurb::computeProjSmagFromNut_fromNut_fromModes(volScalarField& phi, 
        const volScalarField& Nut, const volVectorField& modeU_proj, const volVectorField& modeU_grad)
{
    phi = projDiffusionIBP(Nut, modeU_grad, modeU_proj);
}

volScalarField UnsteadyNSTurb::projDiffusionIBP(const volScalarField& coefDiff, 
                               const volVectorField& u, const volVectorField& v)
{
    volTensorField symGradU = computeS_fromU(u);
    volTensorField symGradV = computeS_fromU(v);
    return 2*coefDiff*(dev(symGradU) && dev(symGradV));
}


// =========================================================================
// ============================= NUT =======================================
// =========================================================================

void UnsteadyNSTurb::computeNonLinearSnapshot_at_time(const word& snap_time, volScalarField& phi, 
                                                 std::optional<PtrList<volVectorField>> modesU)
{
    if (ithacaFVParameters->get_DEIMInterpolatedField() == "nut"
        || ITHACAutilities::containsSubstring(ithacaFVParameters->get_DEIMInterpolatedField(), "reducedNut"))
    {
        phi = computeNut_at_time(snap_time, modesU) ;
    }
    else if (ithacaFVParameters->get_DEIMInterpolatedField() == "fullStressFunction_K")
    {
        phi = computeSmagTermPhi_at_time(snap_time, ithacaFVParameters->get_template_field_k());
    }
    else if (ithacaFVParameters->get_DEIMInterpolatedField() == "fullStressFunction_Omega")
    {
        phi = computeSmagTermPhi_at_time(snap_time, ithacaFVParameters->get_template_field_omega());
    }
    else
    {
        Info << "Error : DEIMInterpolatedField not valid : "
        << ithacaFVParameters->get_DEIMInterpolatedField() << endl;
        Info << "DEIM is available for fullStressFunction and nut only." << endl;
        abort();
    }
}

volScalarField UnsteadyNSTurb::initNutFunction()
{
    return volScalarField(ithacaFVParameters->get_template_field_nut());
}

volScalarField UnsteadyNSTurb::computeNut_at_time(const word& snap_time, std::optional<PtrList<volVectorField>> modesU)
{
    // Read the j-th field
    volVectorField snapshotj = ithacaFVParameters->get_template_field_U();
    ITHACAstream::read_snapshot(snapshotj, snap_time, ithacaFVParameters->get_casenameData());

    if (modesU)
    {
    Eigen::MatrixXd coeffs_field_U = ITHACAutilities::getCoeffs(snapshotj, modesU.value(), modesU.value().size(), 1);
    PtrList<volVectorField> reduced_field_U = ITHACAutilities::reconstructFromCoeff(modesU.value(), 
                                                                coeffs_field_U, modesU.value().size());
    snapshotj = ithacaFVParameters->get_meanU() + reduced_field_U[0];
    }

    return computeNut_fromU(snapshotj);
}

volScalarField UnsteadyNSTurb::computeNut_fromU(const volVectorField& snapshotj)
{
    volTensorField S=computeS_fromU(snapshotj);
    return (computeNut_fromS(S));
}

volScalarField UnsteadyNSTurb::computeNut_fromS(const volTensorField& S)
{

    volScalarField delta = ithacaFVParameters->get_delta();
    float Ck = ithacaFVParameters->get_Ck();
    float Ce = ithacaFVParameters->get_Ce();

    // // Incompressible flows:
    // float Cs = std::pow(Ck*std::pow(Ck/Ce,0.5),0.5);
    // return (pow(Cs*delta,2)*sqrt(2*S&&S));

    // OpenFOAM like version
    // Piece of code strongly inspired by Smagorinsky OpenFOAM Methods
    // see https://develop.openfoam.com/Development/openfoam/blob/OpenFOAM-v2012/src/TurbulenceModels/turbulenceModels/LES/Smagorinsky/Smagorinsky.C

    volScalarField a(Ce/delta);
    volScalarField b((2.0/3.0)*tr(S));
    volScalarField c(2*Ck*delta*(dev(S) && S));

    volScalarField k(sqr((-b + sqrt(sqr(b) + 4*a*c))/(2*a)));

    volScalarField nut = ithacaFVParameters->get_template_field_nut();

    // Loop over inner values to conserve BC informations
    for (int i=0; i<nut.size(); i++)
        {
        nut[i] = Ck*delta[i]*std::pow(k[i],0.5);
        }

    // TO DO : make it works
    nut.correctBoundaryConditions();
    return nut;
 }


 // void UnsteadyNSTurb::initTurbModel()
 // {
 //   const volVectorField& U(ithacaFVParameters->get_template_field_U());
 //
 //   Foam::Time runTime(Foam::Time::controlDictName, ".", ithacaFVParameters->get_casenameData());
 //   wordList boundaryTypeNut = ithacaFVParameters->get_template_field_nut().boundaryField().types();
 //
 //   const surfaceScalarField phi
 //     (
 //     IOobject
 //       (
 //       "phi",
 //       runTime.timeName(),
 //       ithacaFVParameters->get_mesh(),
 //       IOobject::READ_IF_PRESENT,
 //       IOobject::AUTO_WRITE
 //       ),
 //     ithacaFVParameters->get_mesh(),
 //     dimensionedScalar("phi", dimensionSet(0,0,0,0,0,0,0), 1.0),
 //     boundaryTypeNut
 //     );
 //
 //    const Foam::singlePhaseTransportModel transportModel_(U, phi);
 //    turbModel = autoPtr<incompressible::turbulenceModel>
 //      (
 //          incompressible::turbulenceModel::New(U, phi, transportModel_)
 //      );
 //
 //    // turbModel->correct();
 //    // volScalarField nut_test( turbModel->nut());
 //
 //
 //  }



 volTensorField UnsteadyNSTurb::computeS_fromU(const volVectorField& snapshotj)
 {
    volTensorField gradV=fvc::grad(snapshotj);
    return 0.5*(gradV+gradV.T());
 }

 volVectorField UnsteadyNSTurb::computeS_fromU(const volScalarField& phij)
 {
    volVectorField gradV=fvc::grad(phij);
    return gradV;
 }
