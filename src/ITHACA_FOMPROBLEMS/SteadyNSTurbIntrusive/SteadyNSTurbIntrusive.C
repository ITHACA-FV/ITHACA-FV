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

#include "steadyNS.H"
#include "SteadyNSTurbIntrusive.H"
#include "viscosityModel.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
SteadyNSTurbIntrusive::SteadyNSTurbIntrusive() {}
SteadyNSTurbIntrusive::SteadyNSTurbIntrusive(int argc, char* argv[])
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
    _simple = autoPtr<simpleControl>
              (
                  new simpleControl
                  (
                      mesh
                  )
              );
    simpleControl& simple = _simple();
#include "createFields.H"
#include "createFvOptions.H"
    //
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include "initContinuityErrs.H"
#pragma GCC diagnostic pop
    //
    turbulence->validate();
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
    M_Assert(bcMethod == "lift" || bcMethod == "penalty",
             "The BC method must be set to lift or penalty in ITHACAdict");
    para = new ITHACAparameters;
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    supex = ITHACAutilities::check_sup();
}


// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //

// Method to performa a truthSolve
void SteadyNSTurbIntrusive::truthSolve(List<scalar> mu_now)
{
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
    volScalarField& p = _p();
    volVectorField& U = _U();
    surfaceScalarField& phi = _phi();
    fv::options& fvOptions = _fvOptions();
    simpleControl& simple = _simple();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
#include "NLsolve.H"
    ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
    ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
    volScalarField _nut(turbulence->nut());
    ITHACAstream::exportSolution(_nut, name(counter), "./ITHACAoutput/Offline/");
    Ufield.append(U);
    Pfield.append(p);
    nutFields.append(_nut);
    counter++;
    writeMu(mu_now);
    // --- Fill in the mu_samples with parameters (mu) to be used for the PODI sample points
    mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size());

    for (int i = 0; i < mu_now.size(); i++)
    {
        mu_samples(mu_samples.rows() - 1, i) = mu_now[i];
    }

    // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
    if (mu.cols() == 0)
    {
        mu.resize(1, 1);
    }

    if (mu_samples.rows() == mu.cols())
    {
        ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
                                   "./ITHACAoutput/Offline");
    }
}

Eigen::Tensor<double, 3> SteadyNSTurbIntrusive::turbulenceTensor1(label nModes)
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

Eigen::Tensor<double, 3> SteadyNSTurbIntrusive::turbulenceTensor2(label nModes)
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

Eigen::MatrixXd SteadyNSTurbIntrusive::btTurbulence(label nModes)
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

void SteadyNSTurbIntrusive::project(fileName folder, label nModes)
{
    nModesOnline = nModes;

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

Eigen::MatrixXd SteadyNSTurbIntrusive::diffusiveTerm(label nModes)
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

Eigen::Tensor<double, 3> SteadyNSTurbIntrusive::convectiveTerm(label nModes)
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

Eigen::MatrixXd SteadyNSTurbIntrusive::pressureGradientTerm(label nModes)
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

List< Eigen::MatrixXd > SteadyNSTurbIntrusive::bcVelocityVec(label nModes)
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

List< Eigen::MatrixXd > SteadyNSTurbIntrusive::bcVelocityMat(label nModes)
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