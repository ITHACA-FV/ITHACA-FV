/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2020 by the ITHACA-FV authors
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

#include "UnsteadyNSExplicitTurb.H"
#include "fvCFD.H"

/// \file
/// Source file of the UnsteadyNSExplicitTurb class.

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Construct Null
UnsteadyNSExplicitTurb::UnsteadyNSExplicitTurb() {}

// Construct from zero
UnsteadyNSExplicitTurb::UnsteadyNSExplicitTurb(int argc, char* argv[])
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
#include "createFields.H"
    para = ITHACAparameters::getInstance(mesh, runTime);
    bcMethod = ITHACAdict->lookupOrDefault<word>("bcMethod", "none");
    M_Assert(bcMethod == "lift" || bcMethod == "penalty" || bcMethod == "none",
             "The BC method must be set to lift or penalty or none in ITHACAdict");
    fluxMethod = ITHACAdict->lookupOrDefault<word>("fluxMethod", "inconsistent");
    M_Assert(fluxMethod == "inconsistent" || fluxMethod == "consistent",
             "The flux method must be set to inconsistent or consistent in ITHACAdict");
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    supex = ITHACAutilities::check_sup();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void UnsteadyNSExplicitTurb::truthSolve(List<scalar> mu_now, fileName folder)
{
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
    volScalarField& p = _p();
    volVectorField& U = _U();
    surfaceScalarField& phi = _phi();

    if (fluxMethod == "inconsistent")
    {
        phi = fvc::flux(U);
    }

#include "initContinuityErrs.H"
    dimensionedScalar nu = _nu() * mu_now[0];
    dimensionedScalar dt = timeStep * _dt();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    // Perform a TruthSolve
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;
    // Export and store the initial conditions for velocity, pressure and flux
    ITHACAstream::exportSolution(U, name(counter), folder);
    ITHACAstream::exportSolution(p, name(counter), folder);
    ITHACAstream::exportSolution(phi, name(counter), folder);
    volScalarField _nut(turbulence->nut());
    ITHACAstream::exportSolution(_nut, name(counter), folder);
    std::ofstream of(folder + name(counter) + "/" +
                     runTime.timeName());
    Ufield.append(U.clone());
    Pfield.append(p.clone());
    Phifield.append(phi.clone());
    nutFields.append(_nut.clone());
    counter++;
    nextWrite += writeEvery;

    // Start the time loop
    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime.setEndTime(finalTime);
        runTime++;
        Info << "Time = " << runTime.timeName() << nl << endl;

        if (fluxMethod == "inconsistent")
        {
#include "IFM.H"
        }
        else if (fluxMethod == "consistent")
        {
#include "CFM.H"
        }

        turbulence->correct();
#include "initContinuityErrs.H"
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;

        if (checkWrite(runTime))
        {
            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(p, name(counter), folder);
            ITHACAstream::exportSolution(phi, name(counter), folder);
            volScalarField _nut(turbulence->nut());
            ITHACAstream::exportSolution(_nut, name(counter), folder);
            std::ofstream of(folder + name(counter) + "/" +
                             runTime.timeName());
            Ufield.append(U.clone());
            Pfield.append(p.clone());
            Phifield.append(phi.clone());
            nutFields.append(_nut.clone());
            counter++;
            nextWrite += writeEvery;
        }
    }
}

//////////////   AGGIUNTA DELLA TURBOLENZA (NUT)   //////////////


// Eigen::Tensor<double, 3> UnsteadyNSExplicitTurb::diffusive_term_consistent_turb(label NUmodes,
//     label NPmodes,
//     label Nnutmodes,
//     label NSUPmodes)
// {
//     label DFsize = NUmodes + NSUPmodes + Nnutmodes + liftfield.size();
//     Eigen::Tensor<double, 3> DF_tensor;
//     DF_tensor.resize(Nnutmodes, DFsize, NUmodes);
//     surfaceScalarField phi_tmp("Phi_tmp", _phi());

//     for (label i = 0; i < Nnutmodes; i++)
//     {
//         for (label j = 0; j < DFsize; j++)
//         {
//             for (label k = 0; k < NUmodes; k++)
//             {
//                 phi_tmp = dt_dummy * nu_dummy() * fvc::flux(fvc::laplacian(
//                 nutModes[i] , L_U_SUPmodes[k]));

//                 volVectorField CoeffB = fvc::reconstruct(phi_tmp).ref();
//                 volVectorField CoeffA = fvc::reconstruct(L_PHImodes[j]).ref();

//                 DF_tensor(i, j, k) = fvc::domainIntegrate(CoeffA & CoeffB).value();
//             }
//         }
//     }

//     if (Pstream::parRun())
//     {
//         reduce(DF_tensor, sumOp<Eigen::Tensor<double, 3 >> ());
//     }

//     if (Pstream::master())
//     {
//         ITHACAstream::SaveDenseTensor(DF_tensor, "./ITHACAoutput/Matrices/",
//                                   "DF_tensor_" + name(NUmodes) + "_" + name(
//                                       NSUPmodes) + "_" + name(Nnutmodes));
//         // Salvataggio della matrice in formato pyhton
//         cnpy::save(DF_tensor, "./ITHACAoutput/Matrices/DF_tensor_" + name(
//         NUmodes) + "_" + name(NSUPmodes) + "_" + name(Nnutmodes) + ".npy");
//     }

//     return DF_tensor;

// }



// Eigen::Tensor<double, 3> UnsteadyNSExplicitTurb::diffusive_term_turb(label NUmodes,
//     label NPmodes,
//     label Nnutmodes,
//     label NSUPmodes)
// {
//     label Bsize = NUmodes + NSUPmodes + Nnutmodes + liftfield.size();
//     Eigen::Tensor<double, 3> B_tensor;
//     B_tensor.resize(Bsize, Bsize, Bsize);
//     // Eigen::MatrixXd B_matrix;
//     // B_matrix.resize(Bsize, Bsize);

//     // Project everything
//     for (label i = 0; i < Bsize; i++)
//     {
//         for (label j = 0; j < Bsize; j++)
//         {
//             for (label k = 0; k < Bsize; k++)
//             {
//                 B_tensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::laplacian(
//                     nutModes[j], L_U_SUPmodes[k])).value();
//             }
//         }
//     }

//     if (Pstream::parRun())
//     {
//         reduce(B_tensor, sumOp<Eigen::Tensor<double, 3 >> ());
//     }

//     if (Pstream::master())
//     {
//         ITHACAstream::SaveDenseTensor(B_tensor, "./ITHACAoutput/Matrices/",
//                                   "B_tensor_" + name(liftfield.size()) + "_" + name(
//                                     NUmodes) + "_" + name(NSUPmodes) + "_" + name(
//                                     Nnutmodes));
//         // Salvataggio della matrice in formato python 
//         cnpy::save(B_tensor, "./ITHACAoutput/Matrices/B_tensor_" + name(
//         liftfield.size()) + "_" + name(NUmodes) + "_" + name(
//         NSUPmodes) + "_" + name(Nnutmodes) + ".npy");
//     }

//     return B_tensor;

// }



// Eigen::Tensor<double, 3> UnsteadyNSExplicitTurb::diffusive_term_sym_turb(label NUmodes,
//     label NPmodes,
//     label Nnutmodes,
//     label NSUPmodes)
// {
//     label Bsize = NUmodes + NSUPmodes + Nnutmodes + liftfield.size();
//     Eigen::Tensor<double, 3> B_tensor;
//     B_tensor.resize(Bsize, Bsize, Bsize);
//     // Eigen::MatrixXd B_matrix;
//     // B_matrix.resize(Bsize, Bsize);

//     // Project everything
//     for (label i = 0; i < Bsize; i++)
//     {
//         for (label j = 0; j < Bsize; j++)
//         {
//             for (label k = 0; k < Bsize; k++)
//             {
//                 B_tensor(i, j, k) = - fvc::domainIntegrate(fvc::grad(L_U_SUPmodes[i])
//                             && fvc::grad(L_U_SUPmodes[k]) && nutModes[j]).value();
//             }
//         }
//     }

//     if (Pstream::parRun())
//     {
//          reduce(B_tensor, sumOp<Eigen::Tensor<double, 3 >> ());
//     }

//     if (Pstream::master())
//     {
//          ITHACAstream::SaveDenseTensor(B_tensor, "./ITHACAoutput/Matrices/",
//                                   "B_tensor_" + name(liftfield.size()) + "_" + name(
//                                     NUmodes) + "_" + name(NSUPmodes) + "_" + name(
//                                     Nnutmodes));
//         // Salvataggio della matrice in formato python
//          cnpy::save(B_tensor, "./ITHACAoutput/Matrices/B_tensor_" + name(
//             liftfield.size()) + "_" + name(NUmodes) + "_" + name(
//              NSUPmodes)+ "_" + name(Nnutmodes) + ".npy");
//     }

//     return B_tensor;
// }



// Eigen::Tensor<double, 3> UnsteadyNSExplicitTurb::diffusive_term_flux_method_turb(label NUmodes,
//     label NPmodes,
//     label Nnutmodes,
//     label NSUPmodes)
// {
//     label BPsize1 = NPmodes;
//     label BPsize2 = NUmodes + NSUPmodes + Nnutmodes + liftfield.size();
//     Eigen::Tensor<double, 3> BP_tensor;
//     BP_tensor.resize(BPsize1, BPsize2, BPsize2);
//     // Eigen::MatrixXd BP_matrix(BPsize1, BPsize2);
//     volVectorField L_U_SUPmodesaux(L_U_SUPmodes[0]);

//     for (label i = 0; i < BPsize1; i++)
//     {
//         for (label j = 0; j < BPsize2; j++)
//         {
//             for (label k = 0; k < BPsize2; k++)
//             {
//                 L_U_SUPmodesaux = dt_dummy * fvc::laplacian(
//                                 nu_dummy() * nutModes[j], L_U_SUPmodes[k]);
//                 BP_tensor(i, j, k) = fvc::domainIntegrate(Pmodes[i] *
//                                 fvc::div(L_U_SUPmodesaux)).value();
//             }
//         }
//     }

//     if (Pstream::parRun())
//     {
//         reduce(BP_tensor, sumOp<Eigen::Tensor<double, 3 >> ());
//     }

//     if (Pstream::master())
//     {
//         ITHACAstream::SaveDenseTensor(BP_tensor, "./ITHACAoutput/Matrices/",
//                                     "BP_tensor_" + name(NUmodes) + "_" + name(
//                                     NPmodes) + "_" + name(Nnutmodes));
//         // Salvataggio della matrice in formato python
//         cnpy::save(BP_tensor, "./ITHACAoutput/Matrices/BP_tensor_" + name(NUmodes) + "_" + name(
//                                     NPmodes) + "_" + name(Nnutmodes) + ".npy");
//     }

//     return BP_tensor;
// }



// List<Eigen::Tensor<double, 3>> UnsteadyNSExplicitTurb::boundary_vector_diffusion_turb(
//     label NUmodes,
//     label NPmodes,
//     label Nnutmodes,
//     label NSUPmodes)
// {
//     Eigen::VectorXd ModeVector;
//     label BCsize = inletIndex.rows();
//     label RDsize = NUmodes + NSUPmodes + Nnutmodes;
//     List<Eigen::Tensor<double, 3>> RD_tensors(BCsize);
//     Eigen::SparseMatrix<double> A;
//     Eigen::VectorXd b;
//     surfaceScalarField phi_tmp("Phi_tmp", _phi());

//     for (label i = 0; i < BCsize; i++)
//     {
//         label BCind = inletIndex(i, 0);
//         Vector<double> v(0, 0, 0);
//         v[inletIndex(i, 1)] = 1;
//         volVectorField Upara(Uinl());
//         assignBC(Upara, BCind, v);

//         fvVectorMatrix UEqn
//         (
//             -fvm::laplacian(nu_dummy(), Upara)
//         );

//         Foam2Eigen::fvMatrix2Eigen(UEqn, A, b);
//         Eigen::Tensor<double, 3> RD_tensor(RDsize, RDsize, 1);
//         RD_tensors[i] = RD_tensor;

//         for (label j = 0; j < RDsize; j++)
//         {
//             phi_tmp = fvc::flux(fvc::laplacian(nutModes[j], Upara)); 

//             for (label k = 0; k < RDsize; k++)
//             {
//                 ModeVector = Foam2Eigen::field2Eigen(L_U_SUPmodes[k]);   
//                 RD_tensor(j,k,0) = ModeVector.dot(b.col(0));
//             }

//         }

//         if (Pstream::parRun())
//         {
//             reduce(RD_tensors[i], sumOp<Eigen::Tensor<double, 3>>());
//         }

//         if (Pstream::master())
//         {
//             ITHACAstream::SaveDenseTensor(RD_tensor, "./ITHACAoutput/Matrices/RD/",
//                                       "RD_tensor" + name(i) + "_" + name(
//                                         NUmodes) + "_" + name(NSUPmodes) + "_" + name(Nnutmodes));

//             // Salvataggio della matrice in formato python e Conversione del tensore in un array C-style per salvarlo con cnpyy
//             std::vector<double> tensor_data(RDsize * RDsize);
//             for (label j = 0; j < RDsize; j++) 
//             {
//                 for (label k = 0; k < RDsize; k++) 
//                 {
//                     tensor_data[j * RDsize + k] = RD_tensors[i](j, k, 0);
//                 }
//             }

//             cnpy::npy_save("./ITHACAoutput/Matrices/RD/RD_tensor" + name(i) + "_" + name(
//                 NUmodes) + "_" + name(NSUPmodes) + "_" + name(Nnutmodes) + ".npy", 
//             tensor_data.data(), {static_cast<size_t>(RDsize), static_cast<size_t>(RDsize)}, "w");
//             // cnpy::save(RD_tensor[i], "./ITHACAoutput/Matrices/RD/RD" + name(i) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + ".npy");
//         }
//     }

//     return RD_tensors;
// }



// List<Eigen::Tensor<double, 3>> UnsteadyNSExplicitTurb::boundary_vector_diffusion_consistent_turb(
//     label NUmodes,
//     label NPmodes,
//     label Nnutmodes,
//     label NSUPmodes)
// {
//     label BCsize = inletIndex.rows();
//     label SDsize = NUmodes + NSUPmodes + Nnutmodes;
//     List<Eigen::Tensor<double, 3>> SD_tensors(SDsize);
//     surfaceScalarField phi_tmp("Phi_tmp", _phi());

//     for (label i = 0; i < BCsize; i++)
//     {
//         label BCind = inletIndex(i, 0);
//         Vector<double> v(0, 0, 0);
//         v[inletIndex(i, 1)] = 1;
//         volVectorField Upara(Uinl());
//         assignBC(Upara, BCind, v);
//         Eigen::Tensor<double, 3> SD_tensor(SDsize, SDsize, 1);
//         SD_tensors[i] = SD_tensor;

//         // Project everything
//         for (label j = 0; j < SDsize; j++)
//         {
//             for (label k = 0; k < SDsize; k++)
//             {
//                 phi_tmp = dt_dummy * fvc::flux(fvc::laplacian(nu_dummy() * nutModes[j], Upara));

//                 volVectorField CoeffB = fvc::reconstruct(L_PHImodes[k]).ref();
//                 volVectorField CoeffA = fvc::reconstruct(phi_tmp).ref();

//                 SD_tensor(i,j,k) = fvc::domainIntegrate(CoeffA & CoeffB).value();
//             }
//         }

//         ITHACAstream::SaveDenseTensor(SD_tensor, "./ITHACAoutput/Matrices/SD/",
//                                       "SD_tensor" + name(i) + "_" + name(NUmodes) + "_" + name(
//                                         NSUPmodes) + "_" + name(Nnutmodes));
//         // salvataggio della matrice in formato pyhton 
//         std::vector<double> tensor_data(SDsize * SDsize);
//         for (label j = 0; j < SDsize; j++) 
//         {
//             for (label k = 0; k < SDsize; k++) 
//             {
//                 tensor_data[j * SDsize + k] = SD_tensors[i](j, k, 0);
//             }
//         }
//         cnpy::npy_save("./ITHACAoutput/Matrices/SD/SD_tensor" + name(i) + "_" + name(NUmodes) + "_" + name(
//                                       NSUPmodes) + "_" + name(Nnutmodes) + ".npy", 
//         tensor_data.data(), {static_cast<size_t>(SDsize), static_cast<size_t>(SDsize)}, "w");
//         // cnpy::save(SD_tensor[i], "./ITHACAoutput/Matrices/SD/SD" + name(i) + "_" + name(NUmodes) + "_" + name(NSUPmodes)+".npy");

//         if (Pstream::parRun())
//         {
//             reduce(SD_tensors[i], sumOp<Eigen::Tensor<double, 3>>());
//         }

//     }

//     return SD_tensors;
// }

////////////////////////////////////////////


