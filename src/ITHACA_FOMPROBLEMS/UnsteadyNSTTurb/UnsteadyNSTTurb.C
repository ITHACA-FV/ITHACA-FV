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

#include "UnsteadyNSTTurb.H"
#include "viscosityModel.H"
#include "alphatJayatillekeWallFunctionFvPatchScalarField.H"


/// \file
/// Source file of the UnsteadyNSTTurb class.

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Construct Null
UnsteadyNSTTurb::UnsteadyNSTTurb() {};

// Construct from zero
UnsteadyNSTTurb::UnsteadyNSTTurb(int argc, char* argv[])
{
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
    _piso = autoPtr<pisoControl>
            (
                new pisoControl
                (
                    mesh
                )
            );
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include "createFields.H"
#include "createFvOptions.H"
#pragma GCC diagnostic pop
    supex = ITHACAutilities::check_sup();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void UnsteadyNSTTurb::truthSolve(List<scalar> mu_now)
{
    Time& runTime = _runTime();
    surfaceScalarField& phi = _phi();
    fvMesh& mesh = _mesh();
#include "initContinuityErrs.H"
    fv::options& fvOptions = _fvOptions();
    pisoControl& piso = _piso();
    volScalarField p = _p();
    volVectorField U = _U();
    volScalarField T = _T();
    volScalarField _nut(turbulence->nut());
    dimensionedScalar nu = _nu();
    dimensionedScalar Pr = _Pr();
    dimensionedScalar Prt = _Prt();
    volScalarField alphat = _alphat();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    //dimensionedScalar alphaEff = _alphaEff();
    IOMRFZoneList& MRF = _MRF();
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    // Perform a TruthSolve
    bool WRITE;
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Start the time loop
    while (runTime.run())
    {
#include "CourantNo.H"
        runTime.setEndTime(finalTime);
        runTime++;
        Info << "Time = " << runTime.timeName() << nl << endl;
        // --- Pressure-velocity PIMPLE corrector loop
        {
#include "UEqn.H"

            // --- Pressure corrector loop
            while (piso.correct())
            {
#include "pEqn.H"
            }
        }
        {
#include "TEqn.H"
            //TEqn.solve();
        }
        volScalarField _nut(turbulence->nut());
        laminarTransport.correct();
        turbulence->correct();
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
        WRITE = checkWrite(runTime);

        if (WRITE)
        {
            ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(_nut, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(T, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(alphat, name(counter), "./ITHACAoutput/Offline/");
            std::ofstream of("./ITHACAoutput/Offline/" + name(counter) + "/" +
                             runTime.timeName());
            Ufield.append(U);
            Pfield.append(p);
            nutFields.append(_nut);
            Tfield.append(T);
            counter++;
            nextWrite += writeEvery;
            writeMu(mu_now);
        }
    }
}

bool UnsteadyNSTTurb::checkWrite(Time& timeObject)
{
    scalar diffnow = mag(nextWrite - atof(timeObject.timeName().c_str()));
    scalar diffnext = mag(nextWrite - atof(timeObject.timeName().c_str()) -
                          timeObject.deltaTValue());

    if ( diffnow < diffnext)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void UnsteadyNSTTurb::liftSolveT()
{
    for (label k = 0; k < inletIndexT.rows(); k++)
    {
        Time& runTime = _runTime();
        fvMesh& mesh = _mesh();
        volScalarField& T = _T();
        simpleControl simple(mesh);
        //IOMRFZoneList& MRF = _MRF();
        //singlePhaseTransportModel& laminarTransport = _laminarTransport();
        //volScalarField& nut = _nut();
        volScalarField& alphat = _alphat();
        //dimensionedScalar& nu = _nu();
        dimensionedScalar& Pr = _Pr();
        dimensionedScalar& Prt = _Prt();
        label BCind = inletIndexT(k, 0);
        volScalarField Tlift("Tlift" + name(k), T);
        instantList Times = runTime.times();
        runTime.setTime(Times[1], 1);
        Info << "Solving a lifting Problem" << endl;
        scalar t1 = 1;
        scalar t0 = 0;

        for (label j = 0; j < T.boundaryField().size(); j++)
        {
            if (j == BCind)
            {
                assignBC(Tlift, j, t1);
            }
            else if (T.boundaryField()[BCind].type() == "zeroGradient")
            {
                assignBC(Tlift, j, t0);
                assignIF(Tlift, t1);
            }
            else if (T.boundaryField()[BCind].type() == "fixedValue")
            {
                assignBC(Tlift, j, t0);
                assignIF(Tlift, t0);
            }
            else
            {
            }
        }

        while (simple.correctNonOrthogonal())
        {
            alphat = turbulence->nut() / Prt;
            alphat.correctBoundaryConditions();
            volScalarField alphaEff("alphaEff", turbulence->nu() / Pr + alphat);
            fvScalarMatrix TEqn
            (
                fvm::ddt(Tlift)
                - fvm::laplacian(alphaEff, Tlift)
            );
            TEqn.solve();
            Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                 << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                 << nl << endl;
        }

        scalar totalTime = mesh.time().value();
        scalar dt = mesh.time().deltaTValue();
        forAll(Tlift, i)
        {
            Tlift[i] = (totalTime * Tlift[i] + dt * Tlift[i] ) / (totalTime + dt);
        }
        Tlift.write();
        liftfieldT.append(Tlift);
    }
}

List <Eigen::MatrixXd> UnsteadyNSTTurb::turbulenceTerm1(label NUmodes,
        label NSUPmodes, label Nnutmodes)
{
    label Csize = NUmodes + NSUPmodes + liftfield.size();
    List <Eigen::MatrixXd> CT1_matrix;
    CT1_matrix.setSize(Csize);

    for (label j = 0; j < Csize; j++)
    {
        CT1_matrix[j].resize(Nnutmodes, Csize);
        CT1_matrix[j] = CT1_matrix[j] * 0;
    }

    PtrList<volVectorField> Together(0);

    // Create PTRLIST with lift, velocities and supremizers

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k]);
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            Together.append(supmodes[k]);
        }
    }

    for (label i = 0; i < Csize; i++)
    {
        Info << "Filling layer number " << i + 1 << " in the matrix CT1_matrix" << endl;

        for (label j = 0; j < Nnutmodes; j++)
        {
            for (label k = 0; k < Csize; k++)
            {
                CT1_matrix[i](j, k) = fvc::domainIntegrate(Together[i] & fvc::laplacian(
                                          nuTmodes[j], Together[k])).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(CT1_matrix, "CT1_matrix", "python",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(CT1_matrix, "CT1_matrix", "matlab",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(CT1_matrix, "CT1_matrix", "eigen",
                               "./ITHACAoutput/Matrices/CT1");
    return CT1_matrix;
}

List <Eigen::MatrixXd> UnsteadyNSTTurb::turbulenceTerm2(label NUmodes,
        label NSUPmodes, label Nnutmodes)
{
    label Csize = NUmodes + NSUPmodes + liftfield.size();
    List <Eigen::MatrixXd> CT2_matrix;
    CT2_matrix.setSize(Csize);

    for (label j = 0; j < Csize; j++)
    {
        CT2_matrix[j].resize(Nnutmodes, Csize);
        CT2_matrix[j] = CT2_matrix[j] * 0;
    }

    PtrList<volVectorField> Together(0);

    // Create PTRLIST with lift, velocities and supremizers

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k]);
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            Together.append(supmodes[k]);
        }
    }

    for (label i = 0; i < Csize; i++)
    {
        Info << "Filling layer number " << i + 1 << " in the matrix CT2_matrix" << endl;

        for (label j = 0; j < Nnutmodes; j++)
        {
            for (label k = 0; k < Csize; k++)
            {
                CT2_matrix[i](j, k) = fvc::domainIntegrate(Together[i] & (fvc::div(
                                          nuTmodes[j] * dev((fvc::grad(Together[k]))().T())))).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(CT2_matrix, "CT2_matrix", "python",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(CT2_matrix, "CT2_matrix", "matlab",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(CT2_matrix, "CT2_matrix", "eigen",
                               "./ITHACAoutput/Matrices/CT2");
    return CT2_matrix;
}

Eigen::MatrixXd UnsteadyNSTTurb::BTturbulence(label NUmodes, label NSUPmodes)
{
    label BTsize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd BT_matrix(BTsize, BTsize);
    BT_matrix = BT_matrix * 0;
    // Create PTRLIST with lift, velocities and supremizers
    PtrList<volVectorField> Together(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k]);
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            Together.append(supmodes[k]);
        }
    }

    // Project everything
    for (label i = 0; i < BTsize; i++)
    {
        for (label j = 0; j < BTsize; j++)
        {
            BT_matrix(i, j) = fvc::domainIntegrate(Together[i] & (fvc::div(dev((T(fvc::grad(
                    Together[j]))))))).value();
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(BT_matrix, "BT_matrix", "python",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(BT_matrix, "BT_matrix", "matlab",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(BT_matrix, "BT_matrix", "eigen",
                               "./ITHACAoutput/Matrices/");
    return BT_matrix;
}

List <Eigen::MatrixXd> UnsteadyNSTTurb::temperatureTurbulenceTerm(
    label NTmodes, label Nnutmodes)
{
    label Stsize = NTmodes + liftfieldT.size();
    List <Eigen::MatrixXd> S_matrix;
    S_matrix.setSize(Stsize);

    for (label j = 0; j < Stsize; j++)
    {
        S_matrix[j].resize(Nnutmodes, Stsize);
    }

    PtrList<volScalarField> Together(0);

    // Create PTRLIST with lift, velocities and supremizers

    if (liftfieldT.size() != 0)
    {
        for (label k = 0; k < liftfieldT.size(); k++)
        {
            Together.append(liftfieldT[k]);
        }
    }

    if (NTmodes != 0)
    {
        for (label k = 0; k < NTmodes; k++)
        {
            Together.append(Tmodes[k]);
        }
    }

    for (label i = 0; i < Stsize; i++)
    {
        Info << "Filling layer number " << i + 1 << " in the matrix S_matrix" << endl;

        for (label j = 0; j < Nnutmodes; j++)
        {
            for (label k = 0; k < Stsize; k++)
            {
                S_matrix[i](j, k) = fvc::domainIntegrate(Together[i] * fvc::laplacian(
                                        nuTmodes[j], Together[k])).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(S_matrix, "S_matrix", "python",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(S_matrix, "S_matrix", "matlab",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(S_matrix, "S_matrix", "eigen",
                               "./ITHACAoutput/Matrices/S");
    return S_matrix;
}

void UnsteadyNSTTurb::projectSUP(fileName folder, label NU, label NP,
                                 label NSUP, label Nnut, label NT)
{
    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {
        B_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/B_mat.txt");
        C_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/C", "C");
        K_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/K_mat.txt");
        P_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/P_mat.txt");
        M_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/M_mat.txt");
        BT_matrix =
            ITHACAstream::readMatrix("./ITHACAoutput/Matrices/BT_matrix_mat.txt");
        CT1_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/CT1",
                                              "CT1_matrix");
        CT2_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/CT2",
                                              "CT2_matrix");
        Q_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/Q", "Q");
        Y_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/Y_mat.txt");
        S_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/S_mat.txt");
        MT_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/MT_mat.txt");
    }
    else
    {
        NUmodes = NU;
        NPmodes = NP;
        NSUPmodes = NSUP;
        Nnutmodes = Nnut;
        NTmodes = NT;
        B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        C_matrix = convective_term(NUmodes, NPmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        P_matrix = divergence_term(NUmodes, NPmodes, NSUPmodes);
        M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
        BT_matrix = BTturbulence(NUmodes, NSUPmodes);
        CT1_matrix = turbulenceTerm1(NUmodes, NSUPmodes, Nnutmodes);
        CT2_matrix = turbulenceTerm2(NUmodes, NSUPmodes, Nnutmodes);
        S_matrix = temperatureTurbulenceTerm(NTmodes, Nnutmodes);
        Q_matrix = convective_term_temperature(NUmodes, NTmodes, NSUPmodes);
        Y_matrix = diffusive_term_temperature(NUmodes, NTmodes, NSUPmodes);
        MT_matrix = mass_term_temperature(NUmodes, NTmodes, NSUPmodes);
    }

    B_total_matrix = B_matrix + BT_matrix;
    C_total_matrix.setSize(C_matrix.size());

    for (label i = 0; i < C_matrix.size(); i++)
    {
        C_total_matrix[i] =  CT2_matrix[i] + CT1_matrix[i];
    }

    NUmodes = NU;
    NPmodes = NP;
    NSUPmodes = NSUP;
    Nnutmodes = Nnut;
    // Get the coeffs for interpolation (the orthonormal one is used because basis are orthogonal)
    Eigen::MatrixXd Ncoeff = ITHACAutilities::getCoeffs(nutFields, nuTmodes);
    ITHACAstream::exportMatrix(Ncoeff, "Ncoeff", "python",
                               "./ITHACAoutput/Matrices/");
    SAMPLES.resize(Nnutmodes);
    rbfsplines.resize(Nnutmodes);

    for (int i = 0; i < Nnutmodes; i++)
    {
        SAMPLES[i] = new SPLINTER::DataTable(1, 1);

        for (int j = 0; j < Ncoeff.cols(); j++)
        {
            SAMPLES[i]->addSample(mu.row(j), Ncoeff(i, j));
        }

        rbfsplines[i] = new SPLINTER::RBFSpline(*SAMPLES[i],
                                                SPLINTER::RadialBasisFunctionType::GAUSSIAN);
        std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
    }
}
