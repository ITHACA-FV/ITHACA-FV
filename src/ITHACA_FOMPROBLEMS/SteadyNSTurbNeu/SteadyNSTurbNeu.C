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
#include "SteadyNSTurbNeu.H"
#include "viscosityModel.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
SteadyNSTurbNeu::SteadyNSTurbNeu() {}

SteadyNSTurbNeu::SteadyNSTurbNeu(int argc, char* argv[])
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
    viscCoeff = ITHACAdict->lookupOrDefault<word>("viscCoeff", "RBF");
    rbfParams = ITHACAdict->lookupOrDefault<word>("rbfParams", "vel");
    M_Assert(rbfParams == "vel" || rbfParams == "params",
             "The rbfParams must be set to vel or params in ITHACAdict");
    rbfKernel = ITHACAdict->lookupOrDefault<word>("rbfKernel", "linear");
    para = ITHACAparameters::getInstance(mesh, runTime);
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    supex = ITHACAutilities::check_sup();
    neumannMethod = ITHACAdict->lookupOrDefault<word>("neumannMethod", "none");
    M_Assert(neumannMethod == "penalty" || neumannMethod == "NeuTerm" || neumannMethod == "none",
             "The neumann BC method must be set to penalty or NeuTerm or none in ITHACAdict");
    nonUniformbc = ITHACAdict->lookupOrDefault<bool>("nonUniformbc", false);
}


// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //

// Method to performa a truthSolve
void SteadyNSTurbNeu::truthSolve(List<scalar> mu_now)
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
#include "NLsolveSteadyNSTurbNeu.H"
    ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
    ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
    volScalarField _nut(turbulence->nut());
    ITHACAstream::exportSolution(_nut, name(counter), "./ITHACAoutput/Offline/");
    Ufield.append(U.clone());
    Pfield.append(p.clone());
    nutFields.append(_nut.clone());
    counter++;
    writeMu(mu_now);
    // --- Fill in the mu_samples with parameters (mu) to be used for the PODI sample points
    mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size());

    for (label i = 0; i < mu_now.size(); i++)
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

Eigen::Tensor<double, 3> SteadyNSTurbNeu::turbulenceTensor1(label NUmodes,
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

Eigen::Tensor<double, 3> SteadyNSTurbNeu::turbulenceTensor1_cache_mem(label NUmodes,
        label NSUPmodes, label nNutModes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct1Tensor(cSize, nNutModes, cSize);

    for (label j = 0; j < nNutModes; ++j)
    {
        // Cache only one row (j fixed)
        List<tmp<volVectorField>> lapRow(cSize);
        for (label k = 0; k < cSize; ++k)
        {
            lapRow[k] = fvc::laplacian(nutModes[j], L_U_SUPmodes[k]);
        }

        for (label i = 0; i < cSize; ++i)
        {
            for (label k = 0; k < cSize; ++k)
            {
                const volVectorField& lapField = lapRow[k]();
                ct1Tensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & lapField).value();
            }
        }
    }

    // Export the tensor
    if (Pstream::master())
    {
        ITHACAstream::SaveDenseTensor(ct1Tensor, "./ITHACAoutput/Matrices/",
                                      "ct1_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                          NSUPmodes) + "_" + name(nNutModes) + "_t");
    }
    return ct1Tensor;
}

Eigen::Tensor<double, 3> SteadyNSTurbNeu::turbulenceTensor2(label NUmodes,
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

    // Export the tensor
    ITHACAstream::SaveDenseTensor(ct2Tensor, "./ITHACAoutput/Matrices/",
                                  "ct2_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(nNutModes) + "_t");
    return ct2Tensor;
}

Eigen::Tensor<double, 3> SteadyNSTurbNeu::turbulenceTensor2_cache_mem(label NUmodes,
        label NSUPmodes, label nNutModes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct2Tensor(cSize, nNutModes, cSize);

    for (label j = 0; j < nNutModes; ++j)
    {
        // Cache only one j-row
        List<tmp<volVectorField>> divRow(cSize);
        for (label k = 0; k < cSize; ++k)
        {
            divRow[k] = fvc::div(nutModes[j] * dev((fvc::grad(L_U_SUPmodes[k]))().T()));
        }

        for (label i = 0; i < cSize; ++i)
        {
            for (label k = 0; k < cSize; ++k)
            {
                const volVectorField& divRowField = divRow[k]();
                ct2Tensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & divRowField).value();
            }
        }
    }

    // Export the tensor
    if (Pstream::master())
    {
        ITHACAstream::SaveDenseTensor(ct2Tensor, "./ITHACAoutput/Matrices/",
                                      "ct2_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                          NSUPmodes) + "_" + name(nNutModes) + "_t");
    }
    return ct2Tensor;
}

Eigen::MatrixXd SteadyNSTurbNeu::btTurbulence(label NUmodes, label NSUPmodes)
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

Eigen::MatrixXd SteadyNSTurbNeu::diffusive_term_sym(label NUmodes, label NSUPmodes)
{
    label Bsize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd B_matrix_sym;
    B_matrix_sym.resize(Bsize, Bsize);

    // Use PtrList for storing volTensorField pointers
    PtrList<volTensorField> grads(Bsize);

    for (label i = 0; i < Bsize; ++i)
    {
        grads.set(i, new volTensorField(fvc::grad(L_U_SUPmodes[i]))); // Explicit evaluation
    }

    // Compute only upper triangle and mirror
    for (label i = 0; i < Bsize; ++i)
    {
        for (label j = i; j < Bsize; ++j)
        {
            scalar val = fvc::domainIntegrate(grads[i] && grads[j]).value();

            B_matrix_sym(i, j) = val;
            B_matrix_sym(j, i) = val;
        }
    }

    if (Pstream::master())
    {
        ITHACAstream::SaveDenseMatrix(B_matrix_sym, "./ITHACAoutput/Matrices/",
                                      "B_sym_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));
    }

    return B_matrix_sym;
}

Eigen::MatrixXd SteadyNSTurbNeu::bc1_diffusive_term_sym(label NUmodes, label NSUPmodes)
{
    label Bsize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd bc1_B_matrix_sym;
    bc1_B_matrix_sym.resize(Bsize, NeumannFields.size());
    
    label NeumannIndex = outletIndex(0, 0);

    const scalarField& magSf = _mesh().boundary()[NeumannIndex].magSf();

    for (label i = 0; i < Bsize; ++i)
    {
        const vectorField& U_mode = L_U_SUPmodes[i].boundaryField()[NeumannIndex];
        for (label j = 0; j < NeumannFields.size(); ++j)
        {
            scalar val = gSum(U_mode & NeumannFields[j] * magSf);
            bc1_B_matrix_sym(i, j) = val;
        }
    }

    if (Pstream::master())
    {
        ITHACAstream::SaveDenseMatrix(bc1_B_matrix_sym, "./ITHACAoutput/Matrices/",
                                      "bc1_B_sym_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));
    }

    return bc1_B_matrix_sym;
}

Eigen::MatrixXd SteadyNSTurbNeu::bc2_diffusive_term_sym(label NUmodes, label NSUPmodes)
{
    label Bsize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd bc2_B_matrix_sym;
    bc2_B_matrix_sym.resize(Bsize, Bsize);

    label NeumannIndex = outletIndex(0, 0);

    // Use PtrList for storing surfaceVectorField pointers
    PtrList<surfaceVectorField> grads(Bsize);
    for (label i = 0; i < Bsize; ++i)
    {
        grads.set(i, new surfaceVectorField(_mesh().Sf() & fvc::interpolate(fvc::grad(L_U_SUPmodes[i])))); // Explicit evaluation
    }

    // Compute the matrix
    for (label i = 0; i < Bsize; ++i)
    {
        const vectorField& U_mode = fvc::interpolate(L_U_SUPmodes[i]);
        for (label j = 0; j < Bsize; ++j)
        {
            scalar val = gSum(U_mode & grads[j]);
            bc2_B_matrix_sym(i, j) = val;
        }
    }

    const scalarField& magSf = _mesh().boundary()[NeumannIndex].magSf();

    List<vectorField> gradsNeu(Bsize);

    for (label i = 0; i < Bsize; ++i)
    {
        gradsNeu[i] = L_U_SUPmodes[i].boundaryField()[NeumannIndex].snGrad();
    }

    for (label i = 0; i < Bsize; ++i)
    {
        const vectorField& U_mode = L_U_SUPmodes[i].boundaryField()[NeumannIndex];
        for (label j = 0; j < Bsize; ++j)
        {
            scalar val = gSum(U_mode & gradsNeu[j] * magSf);
            bc2_B_matrix_sym(i, j) = bc2_B_matrix_sym(i, j) - val;
        }
    }

    if (Pstream::master())
    {
        ITHACAstream::SaveDenseMatrix(bc2_B_matrix_sym, "./ITHACAoutput/Matrices/",
                                      "bc2_B_sym_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));
    }

    return bc2_B_matrix_sym;
}

Eigen::Tensor<double, 3> SteadyNSTurbNeu::turbulenceTensor1_sym(label NUmodes,
        label NSUPmodes, label nNutModes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct1Tensor_sym(cSize, nNutModes, cSize);

    // Use PtrList for storing volTensorField pointers
    PtrList<volTensorField> grads(cSize);

    for (label i = 0; i < cSize; ++i)
    {
        grads.set(i, new volTensorField(fvc::grad(L_U_SUPmodes[i])())); // Explicit evaluation
    }

    for (label i = 0; i < cSize; ++i)
    {
        for (label j = 0; j < nNutModes; ++j)
        {
            for (label k = i; k < cSize; ++k)
            {
                // Compute only upper triangle and mirror
                scalar val = fvc::domainIntegrate(grads[i] && (nutModes[j] * grads[k])).value();
                ct1Tensor_sym(i, j, k) = val;
                ct1Tensor_sym(k, j, i) = val; // Symmetric part
            }
        }
    }

    // Export the tensor
    if (Pstream::master())
    {
        ITHACAstream::SaveDenseTensor(ct1Tensor_sym, "./ITHACAoutput/Matrices/",
                                      "ct1_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                          NSUPmodes) + "_" + name(nNutModes) + "_t_sym");
    }
    return ct1Tensor_sym;
}

Eigen::Tensor<double, 3> SteadyNSTurbNeu::turbulenceTensor2_sym(label NUmodes,
        label NSUPmodes, label nNutModes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct2Tensor_sym(cSize, nNutModes, cSize);

    // Use PtrList for storing volTensorField pointers
    PtrList<volTensorField> grads(cSize);

    for (label i = 0; i < cSize; ++i)
    {
        grads.set(i, new volTensorField(fvc::grad(L_U_SUPmodes[i])())); // Explicit evaluation
    }

    for (label i = 0; i < cSize; ++i)
    {
        for (label j = 0; j < nNutModes; ++j)
        {
            for (label k = 0; k < cSize; ++k)
            {
                // Compute only upper triangle and mirror
                scalar val = fvc::domainIntegrate(grads[i] && (nutModes[j] * grads[k].T())).value();
                ct2Tensor_sym(i, j, k) = val;
            }
        }
    }

    // Export the tensor
    if (Pstream::master())
    {
        ITHACAstream::SaveDenseTensor(ct2Tensor_sym, "./ITHACAoutput/Matrices/",
                                      "ct2_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                          NSUPmodes) + "_" + name(nNutModes) + "_t_sym");
    }
    return ct2Tensor_sym;
}

Eigen::Tensor<double, 3> SteadyNSTurbNeu::turbulenceTensor_BC(label NUmodes,
        label NSUPmodes, label nNutModes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ctTensor_BC(cSize, nNutModes, cSize);

    label NeumannIndex = outletIndex(0,0);

    // Use PtrList for storing volTensorField pointers
    PtrList<volTensorField> grads(cSize);

    for (label i = 0; i < cSize; ++i)
    {
        grads.set(i, new volTensorField(fvc::grad(L_U_SUPmodes[i])())); // Explicit evaluation
    }

    const vectorField& Sf = _mesh().boundary()[NeumannIndex].Sf();

    for (label i = 0; i < cSize; ++i)
    {
        const vectorField& U_mode = L_U_SUPmodes[i].boundaryField()[NeumannIndex];
        
        for (label j = 0; j < nNutModes; ++j)
        {
            const scalarField& nut_mode = nutModes[j].boundaryField()[NeumannIndex];

            for (label k = 0; k < cSize; ++k)
            {                
                scalarField faceVals(U_mode & (nut_mode * grads[k].boundaryField()[NeumannIndex].T() & Sf));                
                ctTensor_BC(i, j, k) = gSum(faceVals);
            }
        }
    }

    // Export the tensor
    if (Pstream::master())
    {
        ITHACAstream::SaveDenseTensor(ctTensor_BC, "./ITHACAoutput/Matrices/",
                                      "ct_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                          NSUPmodes) + "_" + name(nNutModes) + "_t_BC");
    }
    return ctTensor_BC;
}

Eigen::Tensor<double, 3> SteadyNSTurbNeu::bc_turbulenceTensor(label NUmodes, label NSUPmodes, 
                                                     label nNutModes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> bc_ctTensor;
    bc_ctTensor.resize(cSize, nNutModes, NeumannFields.size());

    label NeumannIndex = outletIndex(0, 0);

    const scalarField& magSf = _mesh().boundary()[NeumannIndex].magSf();

    for (label i = 0; i < cSize; ++i)
    {
        const vectorField& U_mode = L_U_SUPmodes[i].boundaryField()[NeumannIndex];
        
        for (label j = 0; j < nNutModes; ++j)
        {
            const scalarField& nut_mode = nutModes[j].boundaryField()[NeumannIndex];

            for (label k = 0; k < NeumannFields.size(); ++k)
            {
                scalar var = gSum(U_mode & (nut_mode * NeumannFields[k]) * magSf);
                bc_ctTensor(i, j, k) = var;
            }
        }
    }

    if (Pstream::master())
    {
        ITHACAstream::SaveDenseTensor(bc_ctTensor, "./ITHACAoutput/Matrices/",
                                      "bc_ct_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));
    }

    return bc_ctTensor;
}

List<Eigen::MatrixXd> SteadyNSTurbNeu::bcGradVelocityVec(label NUmodes,
        label NSUPmodes)
{    
    label BCsize = NUmodes + NSUPmodes + liftfield.size();
    List <Eigen::MatrixXd> bcGradVelVec(outletIndex.rows());

    M_Assert(outletIndex.rows() != 0,
             "The outletIndex should be assigned.");

    for (label j = 0; j < outletIndex.rows(); j++)
    {
        bcGradVelVec[j].resize(BCsize, NeumannFields.size());
    }

    for (label k = 0; k < outletIndex.rows(); k++)
    {
        label BCind = outletIndex(k, 0);
        const scalarField& magSf = _mesh().boundary()[BCind].magSf();

        for (label i = 0; i < BCsize; i++)
        {
            for (label j = 0; j < NeumannFields.size(); j++)
            {
                // Compute the gradient of the velocity modes at the outlet
                bcGradVelVec[k](i, j) = gSum(L_U_SUPmodes[i].boundaryField()[BCind] &
                    NeumannFields[j] * magSf);
            }
        }
    }

    if (Pstream::master())
    {
        ITHACAstream::exportMatrix(bcGradVelVec, "bcGradVelVec", "eigen",
                                   "./ITHACAoutput/Matrices");
    }
    
    return bcGradVelVec;
}

List<Eigen::MatrixXd> SteadyNSTurbNeu::bcGradVelocityMat(label NUmodes,
        label NSUPmodes)
{        
    label BCsize = NUmodes + NSUPmodes + liftfield.size();
    List <Eigen::MatrixXd> bcGradVelMat(outletIndex.rows());

    M_Assert(outletIndex.rows() != 0,
             "The outletIndex should be assigned.");

    for (label j = 0; j < outletIndex.rows(); j++)
    {
        bcGradVelMat[j].resize(BCsize, BCsize);
    }

    for (label k = 0; k < outletIndex.rows(); k++)
    {
        label BCind = outletIndex(k, 0);
        const scalarField& magSf = _mesh().boundary()[BCind].magSf();
        
        for (label i = 0; i < BCsize; i++)
        {
            for (label j = 0; j < BCsize; j++)
            {
                bcGradVelMat[k](i, j) = gSum(L_U_SUPmodes[i].boundaryField()[BCind] & 
                            L_U_SUPmodes[j].boundaryField()[BCind].snGrad() * magSf);
            }
        }
    }

    if (Pstream::master())
    {
        ITHACAstream::exportMatrix(bcGradVelMat, "bcGradVelMat", "eigen",
                                   "./ITHACAoutput/Matrices");
    }

    return bcGradVelMat;
}

void SteadyNSTurbNeu::getRBFType(const word& viscCoeff, const word& rbfKernel)
{
    if (viscCoeff == "RBF")
    {
        if (rbfKernel == "linear")
        {
            rbfType = SPLINTER::RadialBasisFunctionType::LINEAR;
        }
        else if (rbfKernel == "thinPlate")
        {
            rbfType = SPLINTER::RadialBasisFunctionType::THIN_PLATE_SPLINE;
        }
        else if (rbfKernel == "multiQuadric")
        {
            rbfType = SPLINTER::RadialBasisFunctionType::MULTIQUADRIC;
        }
        else if (rbfKernel == "inverseQuadric")
        {
            rbfType = SPLINTER::RadialBasisFunctionType::INVERSE_MULTIQUADRIC;
        }
        else if (rbfKernel == "inverseMultiQuadric")
        {
            rbfType = SPLINTER::RadialBasisFunctionType::INVERSE_MULTIQUADRIC;
        }
        else if (rbfKernel == "gaussian")
        {
            rbfType = SPLINTER::RadialBasisFunctionType::GAUSSIAN;
        }
        else
        {
            Info<< "Available RBF kernels are: linear, thinPlate, multiQuadric, "
                << "inverseQuadric, inverseMultiQuadric, gaussian." << endl;
            Info<< "Current rbfKernel is: " << rbfKernel << endl;
            FatalError << "Unknown RBF kernel type: " << rbfKernel << endl;
            FatalError.exit();
        }
    }
    else
    {
        FatalError << "Unknown viscCoeff type: " << viscCoeff << endl;
        FatalError.exit();
    }
}

void SteadyNSTurbNeu::projectSUP(fileName folder, label NU, label NP, label NSUP,
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
            L_U_SUPmodes.append(liftfield[k].clone());
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            L_U_SUPmodes.append(Umodes[k].clone());
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            L_U_SUPmodes.append(supmodes[k].clone());
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
            C_tensor = convective_term_tens_cache_mem(NUmodes, NPmodes, NSUPmodes);
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
            ct1Tensor = turbulenceTensor1_cache_mem(NUmodes, NSUPmodes, nNutModes);
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
            ct2Tensor = turbulenceTensor2_cache_mem(NUmodes, NSUPmodes, nNutModes);
        }

        word bSymStr = "B_sym_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                           NSUPmodes);
        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + bSymStr))
        {
            ITHACAstream::ReadDenseMatrix(B_matrix_sym, "./ITHACAoutput/Matrices/", bSymStr);
        }
        else
        {
            B_matrix_sym = diffusive_term_sym(NUmodes, NSUPmodes);            
        }

        word bc1BSymStr = "bc1_B_sym_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                             NSUPmodes);
        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + bc1BSymStr))
        {
            ITHACAstream::ReadDenseMatrix(bc1_B_matrix_sym, "./ITHACAoutput/Matrices/", bc1BSymStr);
        }
        else
        {
            bc1_B_matrix_sym = bc1_diffusive_term_sym(NUmodes, NSUPmodes);
        }

        word bc2BSymStr = "bc2_B_sym_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                             NSUPmodes);
        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + bc2BSymStr))
        {
            ITHACAstream::ReadDenseMatrix(bc2_B_matrix_sym, "./ITHACAoutput/Matrices/", bc2BSymStr);
        }
        else
        {
            bc2_B_matrix_sym = bc2_diffusive_term_sym(NUmodes, NSUPmodes);
        }

        word ct1SymStr = "ct1_" + name(liftfield.size()) + "_" + name(
                          NUmodes) + "_" + name(
                          NSUPmodes) + "_" + name(nNutModes) + "_t_sym";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct1SymStr))
        {
            ITHACAstream::ReadDenseTensor(ct1Tensor_sym, "./ITHACAoutput/Matrices/", ct1SymStr);
        }
        else
        {
            ct1Tensor_sym = turbulenceTensor1_sym(NUmodes, NSUPmodes, nNutModes);
        }

        word ct2SymStr = "ct2_" + name(liftfield.size()) + "_" + name(
                          NUmodes) + "_" + name(
                          NSUPmodes) + "_" + name(nNutModes) + "_t_sym";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ct2SymStr))
        {
            ITHACAstream::ReadDenseTensor(ct2Tensor_sym, "./ITHACAoutput/Matrices/", ct2SymStr);
        }
        else
        {
            ct2Tensor_sym = turbulenceTensor2_sym(NUmodes, NSUPmodes, nNutModes);
        }

        word ctBCStr = "ct_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                          NSUPmodes) + "_" + name(nNutModes) + "_t_BC";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + ctBCStr))
        {
            ITHACAstream::ReadDenseTensor(ctTensor_BC, "./ITHACAoutput/Matrices/", ctBCStr);
        }
        else
        {
            ctTensor_BC = turbulenceTensor_BC(NUmodes, NSUPmodes, nNutModes);
        }

        word bcctStr = "bc_ct_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes)  + "_" + name(nNutModes) + "_t";
        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + bcctStr))
        {
            ITHACAstream::ReadDenseTensor(bc_ctTensor, "./ITHACAoutput/Matrices/", bcctStr);
        }
        else
        {
            bc_ctTensor = bc_turbulenceTensor(NUmodes, NSUPmodes, nNutModes);
        }

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }

        if (neumannMethod == "penalty")
        {
            bcGradVelVec = bcGradVelocityVec(NUmodes, NSUPmodes);
            bcGradVelMat = bcGradVelocityMat(NUmodes, NSUPmodes);
        }
    }
    else
    {
        L_U_SUPmodes.resize(0);

        if (liftfield.size() != 0)
        {
            for (label k = 0; k < liftfield.size(); k++)
            {
                L_U_SUPmodes.append(liftfield[k].clone());
            }
        }

        if (NUmodes != 0)
        {
            for (label k = 0; k < NUmodes; k++)
            {
                L_U_SUPmodes.append(Umodes[k].clone());
            }
        }

        if (NSUPmodes != 0)
        {
            for (label k = 0; k < NSUPmodes; k++)
            {
                L_U_SUPmodes.append(supmodes[k].clone());
            }
        }

        B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        C_tensor = convective_term_tens_cache_mem(NUmodes, NPmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        P_matrix = divergence_term(NUmodes, NPmodes, NSUPmodes);
        M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
        btMatrix = btTurbulence(NUmodes, NSUPmodes);
        ct1Tensor = turbulenceTensor1_cache_mem(NUmodes, NSUPmodes, nNutModes);
        ct2Tensor = turbulenceTensor2_cache_mem(NUmodes, NSUPmodes, nNutModes);

        B_matrix_sym = diffusive_term_sym(NUmodes, NSUPmodes);
        ct1Tensor_sym = turbulenceTensor1_sym(NUmodes, NSUPmodes, nNutModes);
        ct2Tensor_sym = turbulenceTensor2_sym(NUmodes, NSUPmodes, nNutModes);
        ctTensor_BC = turbulenceTensor_BC(NUmodes, NSUPmodes, nNutModes);
        bc1_B_matrix_sym = bc1_diffusive_term_sym(NUmodes, NSUPmodes);
        bc2_B_matrix_sym = bc2_diffusive_term_sym(NUmodes, NSUPmodes);
        bc_ctTensor = bc_turbulenceTensor(NUmodes, NSUPmodes, nNutModes);

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }

        if (neumannMethod == "penalty")
        {
            bcGradVelVec = bcGradVelocityVec(NUmodes, NSUPmodes);
            bcGradVelMat = bcGradVelocityMat(NUmodes, NSUPmodes);
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

        ITHACAstream::exportMatrix(B_matrix_sym, "B_sym", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ctTensor_BC, "ct", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(bc1_B_matrix_sym, "bc1_B_sym", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(bc2_B_matrix_sym, "bc2_B_sym", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(bc_ctTensor, "bc_ct", "python",
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

        ITHACAstream::exportMatrix(B_matrix_sym, "B_sym", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ctTensor_BC, "ct", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(bc1_B_matrix_sym, "bc1_B_sym", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(bc2_B_matrix_sym, "bc2_B_sym", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(bc_ctTensor, "bc_ct", "matlab",
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
        
        ITHACAstream::exportMatrix(B_matrix_sym, "B_sym", "eigen",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(ctTensor_BC, "ct_", "eigen",
                                   "./ITHACAoutput/Matrices/ct");
        ITHACAstream::exportMatrix(bc1_B_matrix_sym, "bc1_B_sym", "eigen",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(bc2_B_matrix_sym, "bc2_B_sym", "eigen",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(bc_ctTensor, "bc_ct_", "eigen",
                                   "./ITHACAoutput/Matrices/bc_ct");
    }

    bTotalMatrix = B_matrix + btMatrix;
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    cTotalTensor.resize(cSize, nNutModes, cSize);
    cTotalTensor = ct1Tensor + ct2Tensor;
    cTotalTensor_sym.resize(cSize, nNutModes, cSize);
    cTotalTensor_sym = - ct1Tensor_sym - ct2Tensor_sym + ctTensor_BC;
    // Get the coeffs for interpolation (the orthonormal one is used because basis are orthogonal)
    coeffL2 = ITHACAutilities::getCoeffs(nutFields,
                                         nutModes, nNutModes);
    if (rbfParams == "vel")
    {
        coeffL2_vel = ITHACAutilities::getCoeffs(Uomfield, Umodes, NUmodes);
        if (bcMethod == "lift")
        {
            Eigen::MatrixXd coeffL2_tmp (NUmodes + liftfield.size(), coeffL2_vel.cols());
            coeffL2_tmp.topRows(liftfield.size()) = coeffL2_lift;
            coeffL2_tmp.bottomRows(NUmodes) = coeffL2_vel;
            coeffL2_vel = coeffL2_tmp;
        }
    }
    if (Pstream::master())
    {
        if (para->exportPython)
        {
            ITHACAstream::exportMatrix(coeffL2, "coeffL2", "python",
                                       "./ITHACAoutput/Matrices/");
        }
        if (para->exportMatlab)
        {
            ITHACAstream::exportMatrix(coeffL2, "coeffL2", "matlab",
                                       "./ITHACAoutput/Matrices/");
        }
        if (para->exportTxt)
        {
            ITHACAstream::exportMatrix(coeffL2, "coeffL2", "eigen",
                                       "./ITHACAoutput/Matrices/");
        }
        // Save the coeffs for interpolation
        ITHACAstream::SaveDenseMatrix(coeffL2, "./ITHACAoutput/Matrices/",
                                      "coeffL2_nut_" + name(nNutModes));

        if (rbfParams == "vel" && para->exportPython)
        {
            ITHACAstream::exportMatrix(coeffL2_vel, "coeffL2_vel", "python",
                                       "./ITHACAoutput/Matrices/");
        }
        if (rbfParams == "vel" && para->exportMatlab)   
        {
            ITHACAstream::exportMatrix(coeffL2_vel, "coeffL2_vel", "matlab",
                                       "./ITHACAoutput/Matrices/");
        }
        if (rbfParams == "vel" && para->exportTxt)
        {
            ITHACAstream::exportMatrix(coeffL2_vel, "coeffL2_vel", "eigen",
                                       "./ITHACAoutput/Matrices/");
        }
        // Save the coeffs for interpolation of velocity
        ITHACAstream::SaveDenseMatrix(coeffL2_vel, "./ITHACAoutput/Matrices/",
                                      "coeffL2_vel_" + name(NUmodes));
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
                if (rbfParams == "vel")
                {
                    samples[i]->addSample(coeffL2_vel.col(j), coeffL2(i, j));
                }
                else if (rbfParams == "params")
                {
                    samples[i]->addSample(mu.row(j), coeffL2(i, j));
                }
                else
                {
                    FatalError << "Unknown rbfParams type: " << rbfParams << endl;
                    FatalError.exit();
                }
            }

            if (Pstream::master())
            {
                ITHACAstream::ReadDenseMatrix(weights, "./ITHACAoutput/weightsSUP/",
                                              weightName);
            }
            rbfSplines[i] = new SPLINTER::RBFSpline( * samples[i], rbfType, weights);
            Info << "dim of rbfSplines[" << i << "] = " << rbfSplines[i]->getNumVariables() << endl;
            Info << "Constructing RadialBasisFunction for mode " << i + 1 << endl;
        }
        else
        {
            samples[i] = new SPLINTER::DataTable(1, 1);

            for (label j = 0; j < coeffL2.cols(); j++)
            {
                if (rbfParams == "vel")
                {
                    samples[i]->addSample(coeffL2_vel.col(j), coeffL2(i, j));
                }
                else if (rbfParams == "params")
                {
                    samples[i]->addSample(mu.row(j), coeffL2(i, j));
                }
                else
                {
                    FatalError << "Unknown rbfParams type: " << rbfParams << endl;
                    FatalError.exit();
                }
            }

            rbfSplines[i] = new SPLINTER::RBFSpline( * samples[i], rbfType);
            std::cout << "dim of rbfSplines[" << i << "] = " << rbfSplines[i]->getNumVariables() << std::endl;
            if (Pstream::master())
            {
                ITHACAstream::SaveDenseMatrix(rbfSplines[i]->weights,
                    "./ITHACAoutput/weightsSUP/", weightName);
            }
            Info << "Constructing RadialBasisFunction for mode " << i + 1 << endl;
        }
    }
}
