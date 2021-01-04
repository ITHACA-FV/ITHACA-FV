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
/// Source file of the UnsteadyBB class.

#include "UnsteadyBB.H"
#include <cmath>

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
UnsteadyBB::UnsteadyBB() {}
UnsteadyBB::UnsteadyBB(int argc, char* argv[])
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
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include "createFields.H"
#pragma GCC diagnostic pop
#include "createFvOptions.H"
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    supex = ITHACAutilities::check_sup();
}


// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //
#include "fvCFD.H"

// Method to performa a truthSolve
void UnsteadyBB::truthSolve(List<scalar> mu_now)
{
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
#include "initContinuityErrs.H"
    volScalarField& p = _p();
    volVectorField& U = _U();
    volScalarField& p_rgh = _p_rgh();
    volScalarField& T = _T();
    volScalarField& alphat = _alphat();
    volScalarField& rhok = _rhok();
    volScalarField& gh = _gh();
    surfaceScalarField& ghf = _ghf();
    surfaceScalarField& phi = _phi();
    fv::options& fvOptions = _fvOptions();
    pimpleControl& pimple = _pimple();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    dimensionedScalar& beta = _beta();
    dimensionedScalar& TRef = _TRef();
    dimensionedScalar& Pr = _Pr();
    dimensionedScalar& Prt = _Prt();
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;
    // save initial condition in folder 0
    ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
    ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
    ITHACAstream::exportSolution(p_rgh, name(counter), "./ITHACAoutput/Offline/");
    ITHACAstream::exportSolution(T, name(counter), "./ITHACAoutput/Offline/");
    std::ofstream of("./ITHACAoutput/Offline/" + name(counter) + "/" +
                     runTime.timeName());
    Ufield.append(U.clone());
    Pfield.append(p.clone());
    Prghfield.append(p_rgh.clone());
    Tfield.append(T.clone());
    counter++;
    nextWrite += writeEvery;

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
#include "TEqn.H"

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
            ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(p_rgh, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(T, name(counter), "./ITHACAoutput/Offline/");
            std::ofstream of("./ITHACAoutput/Offline/" + name(counter) + "/" +
                             runTime.timeName());
            Ufield.append(U.clone());
            Pfield.append(p.clone());
            Prghfield.append(p_rgh.clone());
            Tfield.append(T.clone());
            counter++;
            nextWrite += writeEvery;
            writeMu(mu_now);
        }

        runTime++;
    }
}

void UnsteadyBB::truthSolve(fileName folder)
{
#include "initContinuityErrs.H"
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
    volScalarField& p = _p();
    volVectorField& U = _U();
    volScalarField& p_rgh = _p_rgh();
    volScalarField& T = _T();
    // volScalarField& nut = _nut();
    volScalarField& alphat = _alphat();
    volScalarField& rhok = _rhok();
    // dimensionedVector& g = _g();
    volScalarField& gh = _gh();
    surfaceScalarField& ghf = _ghf();
    surfaceScalarField& phi = _phi();
    fv::options& fvOptions = _fvOptions();
    pimpleControl& pimple = _pimple();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    // dimensionedScalar& nu = _nu();
    dimensionedScalar& beta = _beta();
    // dimensionedScalar& hRef = _hRef();
    // dimensionedScalar& ghRef = _ghRef();
    dimensionedScalar& TRef = _TRef();
    dimensionedScalar& Pr = _Pr();
    dimensionedScalar& Prt = _Prt();
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;
    // Save initial condition
    ITHACAstream::exportSolution(U, name(counter), folder);
    ITHACAstream::exportSolution(p, name(counter), folder);
    ITHACAstream::exportSolution(T, name(counter), folder);
    ITHACAstream::exportSolution(p_rgh, name(counter), folder);
    std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
    counter++;
    nextWrite += writeEvery;

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
#include "TEqn.H"

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
            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(p, name(counter), folder);
            ITHACAstream::exportSolution(T, name(counter), folder);
            ITHACAstream::exportSolution(p_rgh, name(counter), folder);
            std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
            counter++;
            nextWrite += writeEvery;
        }

        runTime++;
    }
}




bool UnsteadyBB::checkWrite(Time& timeObject)
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

// Method to solve the supremizer problem
void UnsteadyBB::solvesupremizer(word type)
{
    PtrList<volScalarField> P_sup;

    if (type == "modes")
    {
        P_sup = Prghmodes;
    }
    else if (type == "snapshots")
    {
        P_sup = Prghfield;
    }
    else
    {
        std::cout << "You must specify the variable type with either snapshots or modes"
                  << std::endl;
        exit(0);
    }

    if (supex == 1)
    {
        volVectorField U = _U();
        volVectorField Usup
        (
            IOobject
            (
                "Usup",
                U.time().timeName(),
                U.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            U.mesh(),
            dimensionedVector("zero", U.dimensions(), vector::zero)
        );

        if (type == "snapshots")
        {
            ITHACAstream::read_fields(supfield, Usup, "./ITHACAoutput/supfield/");
        }
        else if (type == "modes")
        {
            ITHACAstream::read_fields(supmodes, Usup, "./ITHACAoutput/supremizer/");
        }
        else
        {
        }
    }
    else
    {
        volVectorField U = _U();
        volVectorField Usup
        (
            IOobject
            (
                "Usup",
                U.time().timeName(),
                U.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            U.mesh(),
            dimensionedVector("zero", U.dimensions(), vector::zero)
        );
        dimensionedScalar nu_fake
        (
            "nu_fake",
            dimensionSet(0, 2, -1, 0, 0, 0, 0),
            scalar(1)
        );
        Vector<double> v(0, 0, 0);

        for (label i = 0; i < Usup.boundaryField().size(); i++)
        {
            ITHACAutilities::changeBCtype(Usup, "fixedValue", i);
            assignBC(Usup, i, v);
            assignIF(Usup, v);
        }

        if (type == "snapshots")
        {
            for (label i = 0; i < P_sup.size(); i++)
            {
                fvVectorMatrix u_sup_eqn
                (
                    - fvm::laplacian(nu_fake, Usup)
                );
                solve
                (
                    u_sup_eqn == fvc::grad(P_sup[i])
                );
                supfield.append(Usup.clone());
                ITHACAstream::exportSolution(Usup, name(i + 1), "./ITHACAoutput/supfield/");
            }

            label systemRet =
                system("ln -s ../../constant ./ITHACAoutput/supfield/constant");
            systemRet += system("ln -s ../../0 ./ITHACAoutput/supfield/0");
            systemRet += system("ln -s ../../system ./ITHACAoutput/supfield/system");

            if (systemRet < 0)
            {
                Info << "System Command Failed in steadyNS.C" << endl;
                exit(0);
            }
        }
        else
        {
            for (label i = 0; i < Prghmodes.size(); i++)
            {
                fvVectorMatrix u_sup_eqn
                (
                    - fvm::laplacian(nu_fake, Usup)
                );
                solve
                (
                    u_sup_eqn == fvc::grad(Prghmodes[i])
                );
                supmodes.append(Usup.clone());
                ITHACAstream::exportSolution(Usup, name(i + 1), "./ITHACAoutput/supremizer/");
            }

            label systemRet =
                system("ln -s ../../constant ./ITHACAoutput/supremizer/constant");
            systemRet += system("ln -s ../../0 ./ITHACAoutput/supremizer/0");
            systemRet += system("ln -s ../../system ./ITHACAoutput/supremizer/system");

            if (systemRet < 0)
            {
                Info << "System Command Failed in steadyNS.C" << endl;
                exit(0);
            }
        }
    }
}


// * * * * * * * * * * * * * * lifting function temperature * * * * * * * * * * * * * //


// Method to compute the lifting function for temperature
void UnsteadyBB::liftSolveT()
{
    for (label k = 0; k < inletIndexT.rows(); k++)
    {
        Time& runTime = _runTime();
        fvMesh& mesh = _mesh();
        volScalarField& T = _T();
        volVectorField& U = _U();
        surfaceScalarField& phi = _phi();
        phi = linearInterpolate(U) & mesh.Sf();
        simpleControl simple(mesh);
        // IOMRFZoneList& MRF = _MRF();
        // singlePhaseTransportModel& laminarTransport = _laminarTransport();
        // volScalarField& nut = _nut();
        volScalarField& alphat = _alphat();
        // dimensionedScalar& nu = _nu();
        dimensionedScalar& Pr = _Pr();
        dimensionedScalar& Prt = _Prt();
        label BCind = inletIndexT(k, 0);
        volScalarField Tlift("Tlift" + name(k), T);
        instantList Times = runTime.times();
        runTime.setTime(Times[1], 1);
        Info << "Solving a lifting Problem" << endl;
        scalar t1 = 1;
        scalar t0 = 0;
        alphat = turbulence->nut() / Prt;
        alphat.correctBoundaryConditions();
        volScalarField alphaEff("alphaEff", turbulence->nu() / Pr + alphat);

        for (label j = 0; j < T.boundaryField().size(); j++)
        {
            if (j == BCind)
            {
                assignBC(Tlift, j, t1);
                assignIF(Tlift, t0);
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
                fvm::div(phi, Tlift)
                - fvm::laplacian(alphaEff, Tlift)
            );
            TEqn.solve();
            Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                 << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                 << nl << endl;
        }

        Tlift.write();
        liftfieldT.append(Tlift.clone());
    }
}

// * * * * * * * * * * * * * * Projection Methods * * * * * * * * * * * * * * //

void UnsteadyBB::projectSUP(fileName folder, label NU, label NP, label NT,
                            label NSUP)
{
    NUmodes = NU;
    NTmodes = NT;
    NSUPmodes = NSUP;
    NPrghmodes = NP;

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
            M_matrix = mass_term(NUmodes, NPrghmodes, NSUPmodes);
        }

        word B_str = "B_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + B_str))
        {
            ITHACAstream::ReadDenseMatrix(B_matrix, "./ITHACAoutput/Matrices/", B_str);
        }
        else
        {
            B_matrix = diffusive_term(NUmodes, NPrghmodes, NSUPmodes);
        }

        word C_str = "C_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + C_str))
        {
            ITHACAstream::ReadDenseTensor(C_tensor, "./ITHACAoutput/Matrices/", C_str);
        }
        else
        {
            C_tensor = convective_term_tens(NUmodes, NPrghmodes, NSUPmodes);
        }

        word K_str = "K_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_" + name(NPrghmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + K_str))
        {
            ITHACAstream::ReadDenseMatrix(K_matrix, "./ITHACAoutput/Matrices/", K_str);
        }
        else
        {
            K_matrix = pressure_gradient_term(NUmodes, NPrghmodes, NSUPmodes);
        }

        word H_str = "H_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) +  "_" + name(liftfieldT.size()) + "_" + name(NTmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + H_str))
        {
            ITHACAstream::ReadDenseMatrix(H_matrix, "./ITHACAoutput/Matrices/", H_str);
        }
        else
        {
            H_matrix = buoyant_term(NUmodes, NPrghmodes, NSUPmodes);
        }

        word W_str = "W_" + name(liftfieldT.size()) + "_" + name(NTmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + W_str))
        {
            ITHACAstream::ReadDenseMatrix(W_matrix, "./ITHACAoutput/Matrices/", W_str);
        }
        else
        {
            W_matrix = mass_term_temperature(NUmodes, NTmodes, NSUPmodes);
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
                         NSUPmodes) + "_" + name(NPrghmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + P_str))
        {
            ITHACAstream::ReadDenseMatrix(P_matrix, "./ITHACAoutput/Matrices/", P_str);
        }
        else
        {
            P_matrix = divergence_term(NUmodes, NPrghmodes, NSUPmodes);
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

        L_T_modes.resize(0);

        if (liftfieldT.size() != 0)
        {
            for (label k = 0; k < liftfieldT.size(); k++)
            {
                L_T_modes.append(liftfieldT[k].clone());
            }
        }

        if (NTmodes != 0)
        {
            for (label k = 0; k < NTmodes; k++)
            {
                L_T_modes.append(Tmodes[k].clone());
            }
        }

        M_matrix = mass_term(NUmodes, NPrghmodes, NSUPmodes);
        B_matrix = diffusive_term(NUmodes, NPrghmodes, NSUPmodes);
        C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPrghmodes, NSUPmodes);
        H_matrix = buoyant_term(NUmodes, NTmodes, NSUPmodes);
        W_matrix = mass_term_temperature(NUmodes, NTmodes, NSUPmodes);
        Q_matrix = convective_term_temperature(NUmodes, NTmodes, NSUPmodes);
        Y_matrix = diffusive_term_temperature(NUmodes, NTmodes, NSUPmodes);
        P_matrix = divergence_term(NUmodes, NPrghmodes, NSUPmodes);
    }
}


void UnsteadyBB::projectPPE(fileName folder, label NU, label NP, label NT,
                            label NSUP)
{
    NUmodes = NU;
    NTmodes = NT;
    NSUPmodes = 0;
    NPrghmodes = NP;

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
            M_matrix = mass_term(NUmodes, NPrghmodes, NSUPmodes);
        }

        word B_str = "B_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + B_str))
        {
            ITHACAstream::ReadDenseMatrix(B_matrix, "./ITHACAoutput/Matrices/", B_str);
        }
        else
        {
            B_matrix = diffusive_term(NUmodes, NPrghmodes, NSUPmodes);;
        }

        ITHACAstream::exportMatrix(B_matrix, "B_matrix", "eigen",
                                   "./ITHACAoutput/Matrices/");
        word C_str = "C_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + C_str))
        {
            ITHACAstream::ReadDenseTensor(C_tensor, "./ITHACAoutput/Matrices/", C_str);
        }
        else
        {
            C_tensor = convective_term_tens(NUmodes, NPrghmodes, NSUPmodes);
        }

        word K_str = "K_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_" + name(NPrghmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + K_str))
        {
            ITHACAstream::ReadDenseMatrix(K_matrix, "./ITHACAoutput/Matrices/", K_str);
        }
        else
        {
            K_matrix = pressure_gradient_term(NUmodes, NPrghmodes, NSUPmodes);
        }

        word H_str = "H_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) +  "_" + name(liftfieldT.size()) + "_" + name(NTmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + H_str))
        {
            ITHACAstream::ReadDenseMatrix(H_matrix, "./ITHACAoutput/Matrices/", H_str);
        }
        else
        {
            H_matrix = buoyant_term(NUmodes, NPrghmodes, NSUPmodes);
        }

        word HP_str = "HP_" + name(NPrghmodes) + "_" + name(liftfieldT.size()) + "_" +
                      name(NTmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + HP_str))
        {
            ITHACAstream::ReadDenseMatrix(HP_matrix, "./ITHACAoutput/Matrices/", HP_str);
        }
        else
        {
            HP_matrix = buoyant_term_poisson(NPrghmodes, NTmodes);
        }

        word W_str = "W_" + name(liftfieldT.size()) + "_" + name(NTmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + W_str))
        {
            ITHACAstream::ReadDenseMatrix(W_matrix, "./ITHACAoutput/Matrices/", W_str);
        }
        else
        {
            W_matrix = mass_term_temperature(NUmodes, NTmodes, NSUPmodes);
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
                         NSUPmodes) + "_" + name(NPrghmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + P_str))
        {
            ITHACAstream::ReadDenseMatrix(P_matrix, "./ITHACAoutput/Matrices/", P_str);
        }
        else
        {
            P_matrix = divergence_term(NUmodes, NPrghmodes, NSUPmodes);
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

        L_T_modes.resize(0);

        if (liftfieldT.size() != 0)
        {
            for (label k = 0; k < liftfieldT.size(); k++)
            {
                L_T_modes.append(liftfieldT[k].clone());
            }
        }

        if (NTmodes != 0)
        {
            for (label k = 0; k < NTmodes; k++)
            {
                L_T_modes.append(Tmodes[k].clone());
            }
        }

        M_matrix = mass_term(NUmodes, NPrghmodes, NSUPmodes);
        B_matrix = diffusive_term(NUmodes, NPrghmodes, NSUPmodes);
        C_tensor = convective_term_tens(NUmodes, NPrghmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPrghmodes, NSUPmodes);
        H_matrix = buoyant_term(NUmodes, NTmodes, NSUPmodes);
        W_matrix = mass_term_temperature(NUmodes, NTmodes, NSUPmodes);
        Q_matrix = convective_term_temperature(NUmodes, NTmodes, NSUPmodes);
        Y_matrix = diffusive_term_temperature(NUmodes, NTmodes, NSUPmodes);
        P_matrix = divergence_term(NUmodes, NPrghmodes, NSUPmodes);
        D_matrix = laplacian_pressure(NPrghmodes);
        G_matrix = div_momentum(NUmodes, NPrghmodes);
        BC1_matrix = pressure_BC1(NUmodes, NPrghmodes);
        BC2_matrix = pressure_BC2(NUmodes, NPrghmodes);
        BC3_matrix = pressure_BC3(NUmodes, NPrghmodes);
        HP_matrix = buoyant_term_poisson(NPrghmodes, NTmodes);
    }
}
// * * * * * * * * * * * * * * Momentum Eq. Methods * * * * * * * * * * * * * //
Eigen::MatrixXd UnsteadyBB::pressure_gradient_term(label NUmodes,
        label NPrghmodes,
        label NSUPmodes)
{
    label K1size = NUmodes + NSUPmodes + liftfield.size();
    label K2size = NPrghmodes;
    Eigen::MatrixXd K_matrix(K1size, K2size);
    dimensionedVector g = _g();

    // Project everything
    for (label i = 0; i < K1size; i++)
    {
        for (label j = 0; j < K2size; j++)
        {
            K_matrix(i, j) = fvc::domainIntegrate(L_U_SUPmodes[i] &
                                                  fvc::reconstruct(fvc::snGrad(Prghmodes[j]) *
                                                          Prghmodes[j].mesh().magSf())).value();
        }
    }

    // Export the matrix
    ITHACAstream::SaveDenseMatrix(K_matrix, "./ITHACAoutput/Matrices/",
                                  "K_" + name(liftfield.size()) + "_" + name(NUmodes)
                                  + "_" + name(NSUPmodes) + "_" + name(NPrghmodes));
    return K_matrix;
}

// * * * * * * * * * * * * * * Continuity Eq. Methods * * * * * * * * * * * * * //

Eigen::MatrixXd UnsteadyBB::divergence_term(label NUmodes, label NPrghmodes,
        label NSUPmodes)
{
    label P1size = NPrghmodes;
    label P2size = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd P_matrix(P1size, P2size);

    // Project everything
    for (label i = 0; i < P1size; i++)
    {
        for (label j = 0; j < P2size; j++)
        {
            P_matrix(i, j) = fvc::domainIntegrate(Prghmodes[i] * fvc::div (
                    L_U_SUPmodes[j])).value();
        }
    }

    //Export the matrix
    ITHACAstream::SaveDenseMatrix(P_matrix, "./ITHACAoutput/Matrices/",
                                  "P_" + name(liftfield.size()) + "_" + name(NUmodes) + "_"
                                  + name(NSUPmodes) + "_" + name(NPrghmodes));
    return P_matrix;
}

Eigen::MatrixXd UnsteadyBB::buoyant_term(label NUmodes, label NTmodes,
        label NSUPmodes)
{
    label H1size = NUmodes + liftfield.size() + NSUPmodes;
    label H2size = NTmodes + liftfieldT.size() ;
    Eigen::MatrixXd H_matrix(H1size, H2size);
    dimensionedScalar beta = _beta();
    dimensionedScalar TRef = _TRef();
    dimensionedVector g = _g();
    // volScalarField& gh = _gh();
    surfaceScalarField& ghf = _ghf();

    // Project everything
    for (label i = 0; i < H1size; i++)
    {
        for (label j = 0; j < H2size; j++)
        {
            H_matrix(i, j) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::reconstruct(
                    ghf * fvc::snGrad(1.0 - (beta * (L_T_modes[j] - TRef)))
                    * L_T_modes[j].mesh().magSf())).value();
        }
    }

    //Export the matrix
    ITHACAstream::SaveDenseMatrix(H_matrix, "./ITHACAoutput/Matrices/",
                                  "H_" + name(liftfield.size()) + "_" + name(NUmodes) + "_"
                                  + name(NSUPmodes) +  "_" + name(liftfieldT.size()) + "_"
                                  + name(NTmodes));
    return H_matrix;
}

// * * * * * * * * * * * * * * Additional term for the PPE Method * * * * * * * * * * * * * //
Eigen::MatrixXd UnsteadyBB::buoyant_term_poisson(label NPrghmodes,
        label NTmodes)
{
    label H1size = NPrghmodes;
    label H2size = NTmodes + liftfieldT.size() ;
    Eigen::MatrixXd HP_matrix(H1size, H2size);
    // Create PTRLIST with lift, velocities and temperatures
    dimensionedScalar beta = _beta();
    dimensionedScalar TRef = _TRef();
    dimensionedVector g = _g();
    volScalarField& gh = _gh();
    surfaceScalarField& ghf = _ghf();

    // Project everything
    for (label i = 0; i < H1size; i++)
    {
        for (label j = 0; j < H2size; j++)
        {
            HP_matrix(i, j) = fvc::domainIntegrate(fvc::reconstruct(fvc::snGrad(
                    Prghmodes[i]) *
                                                   Prghmodes[i].mesh().magSf())  &fvc::reconstruct(
                                                           ghf * fvc::snGrad( -(beta * (L_T_modes[j])))
                                                           * L_T_modes[j].mesh().magSf())).value();
        }
    }

    ITHACAstream::SaveDenseMatrix(HP_matrix, "./ITHACAoutput/Matrices/",
                                  "HP_"
                                  + name(NPrghmodes) +  "_" + name(liftfieldT.size()) + "_"
                                  + name(NTmodes));
    return HP_matrix;
}


// * * * * * * * * * * * * * * Energy Eq. Methods * * * * * * * * * * * * * //
List<Eigen::MatrixXd> UnsteadyBB::convective_term_temperature(label NUmodes,
        label NTmodes, label NSUPmodes)
{
    label Qsize = NUmodes + liftfield.size() + NSUPmodes;
    label Qsizet = NTmodes + liftfieldT.size() ;
    List <Eigen::MatrixXd> Q_matrix;
    Q_matrix.setSize(Qsizet);

    for (label j = 0; j < Qsizet; j++)
    {
        Q_matrix[j].resize(Qsize, Qsizet);
    }

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

    // Export the matrix
    ITHACAstream::exportMatrix(Q_matrix, "Q", "python", "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(Q_matrix, "Q", "matlab", "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(Q_matrix, "Q", "eigen", "./ITHACAoutput/Matrices/Q");
    return Q_matrix;
}



Eigen::MatrixXd UnsteadyBB::diffusive_term_temperature(label NUmodes,
        label NTmodes, label NSUPmodes)
{
    label Ysize = NTmodes  + liftfieldT.size();
    Eigen::MatrixXd Y_matrix(Ysize, Ysize);

    for (label i = 0; i < Ysize; i++)
    {
        for (label j = 0; j < Ysize; j++)
        {
            Y_matrix(i, j) = fvc::domainIntegrate(L_T_modes[i] * fvc::laplacian(
                    dimensionedScalar("1", dimless, 1), L_T_modes[j])).value();
        }
    }

    // Export the matrix
    ITHACAstream::SaveDenseMatrix(Y_matrix, "./ITHACAoutput/Matrices/",
                                  "Y_" + name(liftfieldT.size()) + "_" + name(NTmodes));
    return Y_matrix;
}


Eigen::MatrixXd UnsteadyBB::mass_term_temperature(label NUmodes, label NTmodes,
        label NSUPmodes)
{
    label Wsize = NTmodes  + liftfieldT.size();
    Eigen::MatrixXd W_matrix(Wsize, Wsize);

    for (label i = 0; i < Wsize; i++)
    {
        for (label j = 0; j < Wsize; j++)
        {
            W_matrix(i, j) = fvc::domainIntegrate(L_T_modes[i] * L_T_modes[j]).value();
        }
    }

    // Export the matrix
    ITHACAstream::SaveDenseMatrix(W_matrix, "./ITHACAoutput/Matrices/",
                                  "W_" + name(liftfieldT.size()) + "_" + name(NTmodes));
    return W_matrix;
}


// Method to compute the lifting function for velocity
void UnsteadyBB::liftSolve()
{
    for (label k = 0; k < inletIndex.rows(); k++)
    {
        Time& runTime = _runTime();
        surfaceScalarField& phi = _phi();
        fvMesh& mesh = _mesh();
        // volScalarField& p = _p();
        volVectorField& U = _U();
        // volScalarField& p_rgh = _p_rgh();
        volScalarField& UliftBC = _UliftBC();
        // volScalarField& T = _T();
        // dimensionedScalar& nu = _nu();
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
            UliftBC.boundaryField().types()
        );
        label PhiRefCell = 0;
        scalar PhiRefValue = 0.0;
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
        adjustPhi(phi, Ulift, UliftBC);

        while (potentialFlow.correctNonOrthogonal())
        {
            fvScalarMatrix PhiEqn
            (
                fvm::laplacian(dimensionedScalar("1", dimless, 1), Phi)
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
        liftfield.append(Ulift.clone());
    }
}

void UnsteadyBB::change_viscosity(double mu)
{
    const volScalarField& nu =  _laminarTransport().nu();
    volScalarField& ciao = const_cast<volScalarField&>(nu);
    this->assignIF(ciao, mu);

    for (label i = 0; i < ciao.boundaryFieldRef().size(); i++)
    {
        this->assignBC(ciao, i, mu);
    }
}







