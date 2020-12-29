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


/// \file
/// Source file of the unsteadyNST class.


#include "unsteadyNST.H"


unsteadyNST::unsteadyNST() {}

// Construct from zero
unsteadyNST::unsteadyNST(int argc, char* argv[])
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
#include "fvOptions.H"
    _piso = autoPtr<pisoControl>
            (
                new pisoControl
                (
                    mesh
                )
            );
#include "createFields.H"
#include "createFvOptions.H"
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    supex = ITHACAutilities::check_sup();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
#include "fvCFD.H"

void unsteadyNST::truthSolve(List<scalar> mu_now)

{
    Time& runTime = _runTime();
    surfaceScalarField& phi = _phi();
    fvMesh& mesh = _mesh();
#include "initContinuityErrs.H"
    pisoControl& piso = _piso();
    volScalarField& p = _p();
    volVectorField& U = _U();
    volScalarField& T = _T();
    dimensionedScalar DT = _DT();
    dimensionedScalar nu = _nu();
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    // Perform a TruthSolve
    bool WRITE;
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;

    // Start the time loop
    while (runTime.run())
    {
        runTime.setEndTime(finalTime);
        runTime++;
        Info << "Time = " << runTime.timeName() << nl << endl;
#include "CourantNo.H"
        // Momentum predictor
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
            + fvm::div(phi, U)
            - fvm::laplacian(nu, U)
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- Pressure-velocity PISO corrector loop
        while (piso.correct())
        {
            {
#include "pEqn.H"
            }
        }

        {
#include "TEqn.H"
            TEqn.solve();
        }

        //runTime.write();
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
        //laminarTransport.correct();
        //turbulence->correct();
        WRITE = checkWrite(runTime);

        if (WRITE)
        {
            ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(T, name(counter), "./ITHACAoutput/Offline/");
            std::ofstream of("./ITHACAoutput/Offline/" + name(counter) + "/" +
                             runTime.timeName());
            Ufield.append(tmp<volVectorField>(U));
            Pfield.append(tmp<volScalarField>(p));
            Tfield.append(tmp<volScalarField>(T));
            counter++;
            nextWrite += writeEvery;
            writeMu(mu_now);
        }
    }
}

bool unsteadyNST::checkWrite(Time& timeObject)
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

// Method to compute the lifting function for temperature
void unsteadyNST::liftSolveT()
{
    for (label k = 0; k < inletIndexT.rows(); k++)
    {
        Time& runTime = _runTime();
        fvMesh& mesh = _mesh();
        volScalarField T = _T();
        volScalarField p = _p();
        volVectorField U = _U();
        dimensionedScalar DT = _DT();
        label BCind = inletIndexT(k, 0);
        volScalarField Tlift("Tlift" + name(k), T);
        volVectorField Ulift("Ulift" + name(k), U);
        pisoControl potentialFlow(mesh, "potentialFlow");
        instantList Times = runTime.times();
        runTime.setTime(Times[1], 1);
        simpleControl simple(mesh);
        dimensionedScalar nu = _nu();
        runTime.setEndTime(finalTime);
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
            fvScalarMatrix TEqn
            (
                fvm::ddt(Tlift)
                - fvm::laplacian(DT, Tlift)
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
        liftfieldT.append(tmp<volScalarField>(Tlift));
    }
}


// * * * * * * * * * * * * * * Projection Methods * * * * * * * * * * * * * * //

void unsteadyNST::projectSUP(fileName folder, label NU, label NP, label NT,
                             label NSUP)
{
    NUmodes = NU;
    NPmodes = NP;
    NTmodes = NT;
    NSUPmodes = NSUP;

    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {
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

        word MT_str = "MT_" + name(liftfieldT.size()) + "_" + name(NTmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + MT_str))
        {
            ITHACAstream::ReadDenseMatrix(MT_matrix, "./ITHACAoutput/Matrices/", MT_str);
        }
        else
        {
            MT_matrix = mass_term_temperature(NUmodes, NTmodes, NSUPmodes);
        }

        Q_matrix = convective_term_temperature(NUmodes, NTmodes, NSUPmodes);
        word Y_str = "Y_" + name(liftfieldT.size()) + "_" + name(NTmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + Y_str))
        {
            ITHACAstream::ReadDenseMatrix(Y_matrix, "./ITHACAoutput/Matrices/", Y_str);
        }
        else
        {
            Y_matrix = diffusive_term_temperature(NUmodes, NTmodes, NSUPmodes);
        }

        word P_str = "P_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_" + name(NPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + P_str))
        {
            ITHACAstream::ReadDenseMatrix(P_matrix, "./ITHACAoutput/Matrices/", P_str);
        }
        else
        {
            P_matrix = divergence_term(NUmodes, NPmodes, NSUPmodes);
        }
    }
    else
    {
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

        L_T_modes.resize(0);

        if (liftfieldT.size() != 0)
        {
            for (label k = 0; k < liftfieldT.size(); k++)
            {
                L_T_modes.append(tmp<volScalarField>(liftfieldT[k]));
            }
        }

        if (NTmodes != 0)
        {
            for (label k = 0; k < NTmodes; k++)
            {
                L_T_modes.append(tmp<volScalarField>(Tmodes[k]));
            }
        }

        B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        C_matrix = convective_term(NUmodes, NPmodes, NSUPmodes);
        Q_matrix = convective_term_temperature(NUmodes, NTmodes, NSUPmodes);
        Y_matrix = diffusive_term_temperature(NUmodes, NTmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        P_matrix = divergence_term(NUmodes, NPmodes, NSUPmodes);
        M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
        MT_matrix = mass_term_temperature(NUmodes, NTmodes, NSUPmodes);
    }
}


List< Eigen::MatrixXd > unsteadyNST::convective_term_temperature(label NUmodes,
        label NTmodes, label NSUPmodes)
{
    label Qsize = NUmodes + liftfield.size() + NSUPmodes;
    label Qsizet = NTmodes + liftfieldT.size() ;
    List < Eigen::MatrixXd > Q_matrix;
    Q_matrix.setSize(Qsizet);

    for (label j = 0; j < Qsizet; j++)
    {
        Q_matrix[j].resize(Qsize, Qsizet);
    }

    {
        for (label i = 0; i < Qsizet; i++)
        {
            for (label j = 0; j < Qsize; j++)
            {
                for (label k = 0; k < Qsizet; k++)
                {
                    Q_matrix[i](j, k) = fvc::domainIntegrate(L_T_modes[i] * fvc::div(
                                            fvc::interpolate(L_U_SUPmodes[j]) & L_U_SUPmodes[j].mesh().Sf(),
                                            L_T_modes[k])).value();
                }
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(Q_matrix, "Q", "python", "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(Q_matrix, "Q", "matlab", "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(Q_matrix, "Q", "eigen", "./ITHACAoutput/Matrices/Q");
    return Q_matrix;
}

Eigen::MatrixXd unsteadyNST::diffusive_term_temperature(label NUmodes,
        label NTmodes, label NSUPmodes)
{
    label Ysize = NTmodes  + liftfieldT.size();
    Eigen::MatrixXd Y_matrix(Ysize, Ysize);
    dimensionedScalar DT = _DT();
    {
        for (label i = 0; i < Ysize; i++)
        {
            for (label j = 0; j < Ysize; j++)
            {
                Y_matrix(i, j) = fvc::domainIntegrate(L_T_modes[i] * fvc::laplacian(
                        dimensionedScalar("1", dimless, 1), L_T_modes[j])).value();
            }
        }
    }
    // Export the matrix
    ITHACAstream::SaveDenseMatrix(Y_matrix, "./ITHACAoutput/Matrices/",
                                  "Y_" + name(liftfieldT.size()) + "_" + name(NTmodes));
    return Y_matrix;
}

Eigen::MatrixXd unsteadyNST::mass_term_temperature(label NUmodes, label NTmodes,
        label NSUPmodes)
{
    label Ysize = NTmodes  + liftfieldT.size();
    Eigen::MatrixXd MT_matrix(Ysize, Ysize);

    for (label i = 0; i < Ysize; i++)
    {
        for (label j = 0; j < Ysize; j++)
        {
            MT_matrix(i, j) = fvc::domainIntegrate(L_T_modes[i] * L_T_modes[j]).value();
        }
    }

    // Export the matrix
    ITHACAstream::SaveDenseMatrix(MT_matrix, "./ITHACAoutput/Matrices/",
                                  "MT_" + name(liftfieldT.size()) + "_" + name(NTmodes));
    return MT_matrix;
}
// Calculate lifting function for velocity
void unsteadyNST::liftSolve()
{
    for (label k = 0; k < inletIndex.rows(); k++)
    {
        Time& runTime = _runTime();
        surfaceScalarField& phi = _phi();
        fvMesh& mesh = _mesh();
        volScalarField p = _p();
        volVectorField U = _U();
        dimensionedScalar nu = _nu();
        IOMRFZoneList& MRF = _MRF();
        label BCind = inletIndex(k, 0);
        volVectorField Ulift("Ulift" + name(k), U);
        instantList Times = runTime.times();
        runTime.setTime(Times[1], 1);
        pisoControl potentialFlow(mesh, "potentialFlow");
        Info << "Solving a lifting Problem" << endl;
        Vector<double> v1(0, 0, 0);
        v1[inletIndex(k, 1)] = 1;
        Vector<double> v0(0, 0, 0);

        for (label j = 0; j < U.boundaryField().size(); j++)
        {
            if (j == BCind)
            {
                assignBC(Ulift, j, v1);
            }
            else if (U.boundaryField()[BCind].type() == "fixedValue")
            {
                assignBC(Ulift, j, v0);
            }
            else
            {
            }

            assignIF(Ulift, v0);
            phi = linearInterpolate(Ulift) & mesh.Sf();
        }

        Info << "Constructing velocity potential field Phi\n" << endl;
        volScalarField Phi
        (
            IOobject
            (
                "Phi",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Phi", dimLength * dimVelocity, 0),
            p.boundaryField().types()
        );
        label PhiRefCell = 0;
        scalar PhiRefValue = 0;
        setRefCell
        (
            Phi,
            potentialFlow.dict(),
            PhiRefCell,
            PhiRefValue
        );
        mesh.setFluxRequired(Phi.name());
        runTime.functionObjects().start();
        MRF.makeRelative(phi);
        adjustPhi(phi, Ulift, p);

        while (potentialFlow.correctNonOrthogonal())
        {
            fvScalarMatrix PhiEqn
            (
                fvm::laplacian(dimensionedScalar("1", dimless, 1), Phi)
                //fvm::laplacian(nu, Phi)
                ==
                fvc::div(phi)
            );
            PhiEqn.setReference(PhiRefCell, PhiRefValue);
            PhiEqn.solve();

            if (potentialFlow.finalNonOrthogonalIter())
            {
                phi -= PhiEqn.flux();
            }
        }

        MRF.makeAbsolute(phi);
        Info << "Continuity error = "
             << mag(fvc::div(phi))().weightedAverage(mesh.V()).value()
             << endl;
        Ulift = fvc::reconstruct(phi);
        Ulift.correctBoundaryConditions();
        Info << "Interpolated velocity error = "
             << (sqrt(sum(sqr((fvc::interpolate(U) & mesh.Sf()) - phi)))
                 / sum(mesh.magSf())).value()
             << endl;
        Ulift.write();
        liftfield.append(tmp<volVectorField>(Ulift));
    }
}



