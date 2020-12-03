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
/// Source file of the steadyNS class.

#include "steadyNS.H"
#include "viscosityModel.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //
// Constructor
steadyNS::steadyNS() {}
steadyNS::steadyNS(int argc, char* argv[])
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
    M_Assert(bcMethod == "lift" || bcMethod == "penalty"|| bcMethod == "none",
             "The BC method must be set to lift or penalty or none in ITHACAdict");
    para = ITHACAparameters::getInstance(mesh, runTime);
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    supex = ITHACAutilities::check_sup();
}


// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //

// Method to perform a truthSolve
void steadyNS::truthSolve(List<scalar> mu_now)
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
    Ufield.append(U);
    Pfield.append(p);
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

// Method to solve the supremizer problem
void steadyNS::solvesupremizer(word type)
{
    M_Assert(type == "modes"
             || type == "snapshots",
             "You must specify the variable type with either snapshots or modes");
    PtrList<volScalarField> P_sup;

    if (type == "modes")
    {
        P_sup = Pmodes;
    }
    else if (type == "snapshots")
    {
        P_sup = Pfield;
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
        else
        {
            ITHACAstream::read_fields(supmodes, Usup, "./ITHACAoutput/supremizer/");
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
            if (Usup.boundaryField()[i].type() != "processor")
            {
                ITHACAutilities::changeBCtype(Usup, "fixedValue", i);
                assignBC(Usup, i, v);
                assignIF(Usup, v);
            }
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
                supfield.append(Usup);
                ITHACAstream::exportSolution(Usup, name(i + 1), "./ITHACAoutput/supfield/");
            }

            ITHACAutilities::createSymLink("./ITHACAoutput/supfield");
        }
        else
        {
            for (label i = 0; i < Pmodes.size(); i++)
            {
                fvVectorMatrix u_sup_eqn
                (
                    - fvm::laplacian(nu_fake, Usup)
                );
                solve
                (
                    u_sup_eqn == fvc::grad(Pmodes[i])
                );
                supmodes.append(Usup);
                ITHACAstream::exportSolution(Usup, name(i + 1), "./ITHACAoutput/supremizer/");
            }

            ITHACAutilities::createSymLink("./ITHACAoutput/supremizer");
        }
    }
}

// Method to compute the lifting function
void steadyNS::liftSolve()
{
    for (label k = 0; k < inletIndex.rows(); k++)
    {
        Time& runTime = _runTime();
        surfaceScalarField& phi = _phi();
        fvMesh& mesh = _mesh();
        volScalarField p = _p();
        volVectorField U = _U();
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
        liftfield.append(Ulift);
    }
}

// * * * * * * * * * * * * * * Projection Methods * * * * * * * * * * * * * * //

void steadyNS::projectPPE(fileName folder, label NU, label NP, label NSUP)
{
    NUmodes = NU;
    NPmodes = NP;
    NSUPmodes = 0;
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

        word bc4_str = "BC4_" + name(liftfield.size()) + "_" + name(
                           NUmodes) + "_" + name(NPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + bc4_str))
        {
            ITHACAstream::ReadDenseMatrix(BC4_matrix, "./ITHACAoutput/Matrices/", bc4_str);
        }
        else
        {
            BC4_matrix = pressure_BC4(NUmodes, NPmodes);
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
        BC3_matrix = pressure_BC3(NUmodes, NPmodes);
        BC4_matrix = pressure_BC4(NUmodes, NPmodes);

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
        ITHACAstream::exportMatrix(BC3_matrix, "BC3", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC4_matrix, "BC4", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor, "C", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(gTensor, "G", "python", "./ITHACAoutput/Matrices/");
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(D_matrix, "D", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC3_matrix, "BC3", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC4_matrix, "BC4", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor, "C", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(gTensor, "G", "matlab", "./ITHACAoutput/Matrices/");
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(D_matrix, "D", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC3_matrix, "BC3", "eigen",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC4_matrix, "BC4", "eigen",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor, "C", "eigen",
                                   "./ITHACAoutput/Matrices/C");
        ITHACAstream::exportTensor(gTensor, "G", "eigen",
                                   "./ITHACAoutput/Matrices/G");
    }
}

void steadyNS::projectSUP(fileName folder, label NU, label NP, label NSUP)
{
    NUmodes = NU;
    NPmodes = NP;
    NSUPmodes = NSUP;
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
        ITHACAstream::exportTensor(C_tensor, "C", "python", "./ITHACAoutput/Matrices/");
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(P_matrix, "P", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor, "C", "python", "./ITHACAoutput/Matrices/");
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(P_matrix, "P", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor, "C", "python",
                                   "./ITHACAoutput/Matrices/C");
    }
}

void steadyNS::discretizeThenProject(fileName folder, label NU, label NP, label NSUP)
{
    NUmodes = NU;
    NPmodes = NP;
    NSUPmodes = 0;
    L_U_SUPmodes.resize(0);

    Vector<double> inl(0, 0, 0);
    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            L_U_SUPmodes.append(Umodes[k]);
	    
	    // set homogenous boundary conditions
	    forAll(L_U_SUPmodes[k].mesh().boundary(), l)
	    {
		assignBC(L_U_SUPmodes[k],l,inl);
	    }
        }
    }
    
    // Dummy field with all Neumann Boundary conditions converted to homogeneous 
    // Dirichlet boundary conditions
    Uinl = new volVectorField("Uinl", _U());
    assignIF(Uinl(), inl);
    ITHACAutilities::changeNeumann2Dirichlet(Uinl(),inl);
    //changeNeumann2Dirichlet
    /*forAll( Uinl().mesh().boundary(), l)
    {
        if (Uinl().boundaryField()[l].type() == "zeroGradient" ||  
            Uinl().boundaryField()[l].type() ==  "fixedGradient")
        {
            ITHACAutilities::changeBCtype(Uinl(),"fixedValue",l);
            assignBC(Uinl(),l,inl);
        }
    }*/
    
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
        
        word Cf_str = "Cf_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + Cf_str))
        {
            ITHACAstream::ReadDenseTensor(Cf_tensor, "./ITHACAoutput/Matrices/", Cf_str);
        }
        else
        {
            Cf_tensor = convective_term_flux_tens(NUmodes, NPmodes, NSUPmodes);
        }
        
        word BP_str = "BP_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + BP_str))
        {
            ITHACAstream::ReadDenseMatrix(B_matrix, "./ITHACAoutput/Matrices/", BP_str);
        }
        else
        {
            BP_matrix = diffusive_term_flux_method(NUmodes, NPmodes, NSUPmodes);
        }
        
        word RD_str = "RD_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + RD_str))
        {
            ITHACAstream::ReadDenseMatrix(RD_matrix, "./ITHACAoutput/Matrices/", RD_str);
        }
        else
        {
            RD_matrix = boundary_vector_diffusion(NUmodes, NPmodes, NSUPmodes);
        }
        
        word RC_str = "RC_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + RC_str))
        {
            ITHACAstream::ReadDenseMatrix(RC_matrix, "./ITHACAoutput/Matrices/", RC_str);
        }
        else
        {
            RC_matrix = boundary_vector_convection(NUmodes, NPmodes, NSUPmodes);
        }
        
        LinSysDiv = pressure_gradient_term_linsys_div(NPmodes);
	LinSysDiff = pressure_gradient_term_linsys_diff(NPmodes);
	LinSysConv = pressure_gradient_term_linsys_conv(NPmodes);             
    }
    else
    {
        B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
	P_matrix = divergence_term(NUmodes, NPmodes, NSUPmodes);
	BP_matrix = diffusive_term_flux_method(NUmodes, NPmodes, NSUPmodes);
	RD_matrix = boundary_vector_diffusion(NUmodes, NPmodes, NSUPmodes);
	RC_matrix = boundary_vector_convection(NUmodes, NPmodes, NSUPmodes);
        Cf_tensor = convective_term_flux_tens(NUmodes, NPmodes, NSUPmodes);
	LinSysDiv = pressure_gradient_term_linsys_div(NPmodes);
	LinSysDiff = pressure_gradient_term_linsys_diff(NPmodes);
	LinSysConv = pressure_gradient_term_linsys_conv(NPmodes);
    }
}

// * * * * * * * * * * * * * * Momentum Eq. Methods * * * * * * * * * * * * * //

Eigen::MatrixXd steadyNS::diffusive_term(label NUmodes, label NPmodes,
        label NSUPmodes)
{
    label Bsize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd B_matrix;
    B_matrix.resize(Bsize, Bsize);

    // Project everything
    for (label i = 0; i < Bsize; i++)
    {
        for (label j = 0; j < Bsize; j++)
        {
            B_matrix(i, j) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::laplacian(
                    dimensionedScalar("1", dimless, 1), L_U_SUPmodes[j])).value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(B_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(B_matrix, "./ITHACAoutput/Matrices/",
                                  "B_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));
    return B_matrix;
}

Eigen::MatrixXd steadyNS::pressure_gradient_term(label NUmodes, label NPmodes,
        label NSUPmodes)
{
    label K1size = NUmodes + NSUPmodes + liftfield.size();
    label K2size = NPmodes;
    Eigen::MatrixXd K_matrix(K1size, K2size);

    // Project everything
    for (label i = 0; i < K1size; i++)
    {
        for (label j = 0; j < K2size; j++)
        {
            K_matrix(i, j) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::grad(
                    Pmodes[j])).value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(K_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(K_matrix, "./ITHACAoutput/Matrices/",
                                  "K_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(NPmodes));
    return K_matrix;
}

List <Eigen::MatrixXd> steadyNS::convective_term(label NUmodes, label NPmodes,
        label NSUPmodes)
{
    label Csize = NUmodes + NSUPmodes + liftfield.size();
    List <Eigen::MatrixXd> C_matrix;
    C_matrix.setSize(Csize);

    for (label j = 0; j < Csize; j++)
    {
        C_matrix[j].resize(Csize, Csize);
    }

    for (label i = 0; i < Csize; i++)
    {
        for (label j = 0; j < Csize; j++)
        {
            for (label k = 0; k < Csize; k++)
            {
                C_matrix[i](j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::div(
                                        linearInterpolate(L_U_SUPmodes[j]) & L_U_SUPmodes[j].mesh().Sf(),
                                        L_U_SUPmodes[k])).value();
            }
        }
    }

    if (Pstream::parRun())
    {
        for (label i = 0; i < Csize; i++)
        {
            List<double> vec(C_matrix[i].data(), C_matrix[i].data() + C_matrix[i].size());
            reduce(vec, sumOp<List<double>>());
            std::memcpy(C_matrix[i].data(), &vec[0], sizeof (double)*vec.size());
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(C_matrix, "C", "python", "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(C_matrix, "C", "matlab", "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(C_matrix, "C", "eigen", "./ITHACAoutput/Matrices/C");
    return C_matrix;
}

Eigen::Tensor<double, 3> steadyNS::convective_term_tens(label NUmodes,
        label NPmodes,
        label NSUPmodes)
{
    label Csize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> C_tensor;
    C_tensor.resize(Csize, Csize, Csize);

    for (label i = 0; i < Csize; i++)
    {
        for (label j = 0; j < Csize; j++)
        {
            for (label k = 0; k < Csize; k++)
            {
                C_tensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::div(
                                        linearInterpolate(L_U_SUPmodes[j]) & L_U_SUPmodes[j].mesh().Sf(),
                                        L_U_SUPmodes[k])).value();
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(C_tensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(C_tensor, "./ITHACAoutput/Matrices/",
                                  "C_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_t");
    return C_tensor;
}

Eigen::MatrixXd steadyNS::mass_term(label NUmodes, label NPmodes,
                                    label NSUPmodes)
{
    label Msize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd M_matrix(Msize, Msize);

    // Project everything
    for (label i = 0; i < Msize; i++)
    {
        for (label j = 0; j < Msize; j++)
        {
            M_matrix(i, j) = fvc::domainIntegrate(L_U_SUPmodes[i] &
                                                  L_U_SUPmodes[j]).value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(M_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(M_matrix, "./ITHACAoutput/Matrices/",
                                  "M_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));
    return M_matrix;
}

// * * * * * * * * * * * * * * Continuity Eq. Methods * * * * * * * * * * * * * //

Eigen::MatrixXd steadyNS::divergence_term(label NUmodes, label NPmodes,
        label NSUPmodes)
{
    label P1size = NPmodes;
    label P2size = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd P_matrix(P1size, P2size);

    // Project everything
    for (label i = 0; i < P1size; i++)
    {
        for (label j = 0; j < P2size; j++)
        {
            P_matrix(i, j) = fvc::domainIntegrate(Pmodes[i] * fvc::div (
                    L_U_SUPmodes[j])).value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(P_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(P_matrix, "./ITHACAoutput/Matrices/",
                                  "P_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(NPmodes));
    return P_matrix;
}


List <Eigen::MatrixXd> steadyNS::div_momentum(label NUmodes, label NPmodes)
{
    label G1size = NPmodes;
    label G2size = NUmodes + NSUPmodes + liftfield.size();
    List <Eigen::MatrixXd> G_matrix;
    G_matrix.setSize(G1size);

    for (label j = 0; j < G1size; j++)
    {
        G_matrix[j].resize(G2size, G2size);
    }

    for (label i = 0; i < G1size; i++)
    {
        for (label j = 0; j < G2size; j++)
        {
            for (label k = 0; k < G2size; k++)
            {
                G_matrix[i](j, k) = fvc::domainIntegrate(fvc::grad(Pmodes[i]) & (fvc::div(
                                        fvc::interpolate(L_U_SUPmodes[j]) & L_U_SUPmodes[j].mesh().Sf(),
                                        L_U_SUPmodes[k]))).value();
            }
        }
    }

    if (Pstream::parRun())
    {
        for (label i = 0; i < G2size; i++)
        {
            List<double> vec(G_matrix[i].data(), G_matrix[i].data() + G_matrix[i].size());
            reduce(vec, sumOp<List<double>>());
            std::memcpy(G_matrix[i].data(), &vec[0], sizeof (double)*vec.size());
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(G_matrix, "G", "python", "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(G_matrix, "G", "matlab", "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(G_matrix, "G", "eigen", "./ITHACAoutput/Matrices/G");
    return G_matrix;
}

Eigen::Tensor<double, 3> steadyNS::divMomentum(label NUmodes, label NPmodes)
{
    label g1Size = NPmodes;
    label g2Size = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> gTensor;
    gTensor.resize(g1Size, g2Size, g2Size);

    for (label i = 0; i < g1Size; i++)
    {
        for (label j = 0; j < g2Size; j++)
        {
            for (label k = 0; k < g2Size; k++)
            {
                gTensor(i, j, k) = fvc::domainIntegrate(fvc::grad(Pmodes[i]) & (fvc::div(
                        fvc::interpolate(L_U_SUPmodes[j]) & L_U_SUPmodes[j].mesh().Sf(),
                        L_U_SUPmodes[k]))).value();
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(gTensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(gTensor, "./ITHACAoutput/Matrices/",
                                  "G_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(NPmodes) + "_t");
    return gTensor;
}

Eigen::MatrixXd steadyNS::laplacian_pressure(label NPmodes)
{
    label Dsize = NPmodes;
    Eigen::MatrixXd D_matrix(Dsize, Dsize);

    // Project everything
    for (label i = 0; i < Dsize; i++)
    {
        for (label j = 0; j < Dsize; j++)
        {
            D_matrix(i, j) = fvc::domainIntegrate(fvc::grad(Pmodes[i])&fvc::grad(
                    Pmodes[j])).value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(D_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(D_matrix, "./ITHACAoutput/Matrices/",
                                  "D_" + name(NPmodes));
    return D_matrix;
}

Eigen::MatrixXd steadyNS::pressure_BC1(label NUmodes, label NPmodes)
{
    label P_BC1size = NPmodes;
    label P_BC2size = NUmodes + liftfield.size();
    Eigen::MatrixXd BC1_matrix(P_BC1size, P_BC2size);
    fvMesh& mesh = _mesh();

    for (label i = 0; i < P_BC1size; i++)
    {
        for (label j = 0; j < P_BC2size; j++)
        {
            surfaceScalarField lpl((fvc::interpolate(fvc::laplacian(
                                        L_U_SUPmodes[j]))&mesh.Sf())*fvc::interpolate(Pmodes[i]));
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
                                  "BC1_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NPmodes));
    return BC1_matrix;
}


List <Eigen::MatrixXd> steadyNS::pressure_BC2(label NUmodes, label NPmodes)
{
    label P2_BC1size = NPmodes;
    label P2_BC2size = NUmodes + NSUPmodes + liftfield.size();
    List <Eigen::MatrixXd> BC2_matrix;
    fvMesh& mesh = _mesh();
    BC2_matrix.setSize(P2_BC1size);

    for (label j = 0; j < P2_BC1size; j++)
    {
        BC2_matrix[j].resize(P2_BC2size, P2_BC2size);
    }

    for (label i = 0; i < P2_BC1size; i++)
    {
        for (label j = 0; j < P2_BC2size; j++)
        {
            for (label k = 0; k < P2_BC2size; k++)
            {
                surfaceScalarField div_m(fvc::interpolate(fvc::div(fvc::interpolate(
                                             L_U_SUPmodes[j]) & mesh.Sf(),
                                         L_U_SUPmodes[k]))&mesh.Sf()*fvc::interpolate(Pmodes[i]));
                double s = 0;

                for (label k = 0; k < div_m.boundaryField().size(); k++)
                {
                    s += gSum(div_m.boundaryField()[k]);
                }

                BC2_matrix[i](j, k) = s;
            }
        }
    }

    if (Pstream::parRun())
    {
        for (label i = 0; i < P2_BC1size; i++)
        {
            List<double> vec(BC2_matrix[i].data(),
                             BC2_matrix[i].data() + BC2_matrix[i].size());
            reduce(vec, sumOp<List<double>>());
            std::memcpy(BC2_matrix[i].data(), &vec[0], sizeof (double)*vec.size());
        }
    }

    return BC2_matrix;
}

Eigen::Tensor<double, 3 > steadyNS::pressureBC2(label NUmodes, label NPmodes)
{
    label pressureBC1Size = NPmodes;
    label pressureBC2Size = NUmodes + NSUPmodes + liftfield.size();
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
                                             L_U_SUPmodes[j]) & mesh.Sf(),
                                         L_U_SUPmodes[k]))&mesh.Sf()*fvc::interpolate(Pmodes[i]));
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
                                  "BC2_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(NPmodes) + "_t");
    return bc2Tensor;
}

Eigen::MatrixXd steadyNS::pressure_BC3(label NUmodes, label NPmodes)
{
    label P3_BC1size = NPmodes;
    label P3_BC2size = NUmodes + liftfield.size();
    Eigen::MatrixXd BC3_matrix(P3_BC1size, P3_BC2size);
    fvMesh& mesh = _mesh();
    surfaceVectorField n(mesh.Sf() / mesh.magSf());

    for (label i = 0; i < P3_BC1size; i++)
    {
        for (label j = 0; j < P3_BC2size; j++)
        {
            surfaceVectorField BC3 = fvc::interpolate(fvc::curl(L_U_SUPmodes[j]));
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
                                  "BC3_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NPmodes));
    return BC3_matrix;
}

Eigen::MatrixXd steadyNS::pressure_BC4(label NUmodes, label NPmodes)
{
    label P4_BC1size = NPmodes;
    label P4_BC2size = NUmodes + liftfield.size();
    Eigen::MatrixXd BC4_matrix(P4_BC1size, P4_BC2size);
    fvMesh& mesh = _mesh();
    surfaceVectorField n(mesh.Sf() / mesh.magSf());

    for (label i = 0; i < P4_BC1size; i++)
    {
        for (label j = 0; j < P4_BC2size; j++)
        {
            surfaceScalarField BC3 = fvc::interpolate(Pmodes[i]);
            surfaceScalarField BC4 = n & fvc::interpolate(L_U_SUPmodes[j]);
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
                                  "BC4_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NPmodes));
    return BC4_matrix;
}

List< Eigen::MatrixXd > steadyNS::bcVelocityVec(label NUmodes,
        label NSUPmodes)
{
    label BCsize = NUmodes + NSUPmodes;
    List < Eigen::MatrixXd > bcVelVec(inletIndex.rows());

    for (label j = 0; j < inletIndex.rows(); j++)
    {
        bcVelVec[j].resize(BCsize, 1);
    }

    for (label k = 0; k < inletIndex.rows(); k++)
    {
        label BCind = inletIndex(k, 0);
        label BCcomp = inletIndex(k, 1);

        for (label i = 0; i < BCsize; i++)
        {
            bcVelVec[k](i, 0) = gSum(L_U_SUPmodes[i].boundaryField()[BCind].component(
                                         BCcomp));
        }
    }

    ITHACAstream::exportMatrix(bcVelVec, "bcVelVec", "eigen",
                               "./ITHACAoutput/Matrices/bcVelVec");
    return bcVelVec;
}

List< Eigen::MatrixXd > steadyNS::bcVelocityMat(label NUmodes,
        label NSUPmodes)
{
    label BCsize = NUmodes + NSUPmodes;
    label BCUsize = inletIndex.rows();
    List < Eigen::MatrixXd > bcVelMat(BCUsize);

    for (label j = 0; j < inletIndex.rows(); j++)
    {
        bcVelMat[j].resize(BCsize, BCsize);
    }

    for (label k = 0; k < inletIndex.rows(); k++)
    {
        label BCind = inletIndex(k, 0);
        label BCcomp = inletIndex(k, 1);

        for (label i = 0; i < BCsize; i++)
        {
            for (label j = 0; j < BCsize; j++)
            {
                bcVelMat[k](i, j) = gSum(L_U_SUPmodes[i].boundaryField()[BCind].component(
                                             BCcomp) *
                                         L_U_SUPmodes[j].boundaryField()[BCind].component(BCcomp));
            }
        }
    }

    ITHACAstream::exportMatrix(bcVelMat, "bcVelMat", "eigen",
                               "./ITHACAoutput/Matrices/bcVelMat");
    return bcVelMat;
}

Eigen::MatrixXd steadyNS::diffusive_term_flux_method(label NUmodes, label NPmodes, label NSUPmodes)
{
    label BPsize1 = NPmodes;
    label BPsize2 = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd BP_matrix(BPsize1, BPsize2);

    dimensionedScalar dt_fake
    (
        "dt_fake",
        dimensionSet(0, 0, 1, 0, 0, 0, 0),
        scalar(1.0)
    );
    dimensionedScalar nu_fake
    (
        "nu_fake",
        dimensionSet(0, 2, -1, 0, 0, 0, 0),
        scalar(1.0)
    );

    volVectorField L_U_SUPmodesaux(L_U_SUPmodes[0]);
    for (label i = 0; i < BPsize1; i++)
    {
        for (label j = 0; j < BPsize2; j++)
        {
	    L_U_SUPmodesaux = dt_fake*fvc::laplacian(
                    nu_fake, L_U_SUPmodes[j]);

            BP_matrix(i, j) = fvc::domainIntegrate(Pmodes[i]*
		   fvc::div(L_U_SUPmodesaux)).value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(BP_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(BP_matrix, "./ITHACAoutput/Matrices/",
                                  "BP_" + name(NPmodes));
    return BP_matrix;
}

Eigen::MatrixXd steadyNS::boundary_vector_diffusion(label NUmodes, label NPmodes, label NSUPmodes)
{
    dimensionedScalar nu_fake
    (
        "nu_fake",
        dimensionSet(0, 2, -1, 0, 0, 0, 0),
        scalar(1.0)
    );
    // Determine boundary vector
    fvVectorMatrix UEqn
    (
	-fvm::laplacian(nu_fake,Uinl())
    );     
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Foam2Eigen::fvMatrix2Eigen(UEqn, A, b);

    label RDsize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd RD_matrix(RDsize,1);
    Eigen::VectorXd ModeVector;
 
    for (label i = 0; i < RDsize; i++)
    {
        ModeVector = Foam2Eigen::field2Eigen(L_U_SUPmodes[i]);
        RD_matrix(i,0) = ModeVector.dot(b.col(0));
    }

    if (Pstream::parRun())
    {
        reduce(RD_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(RD_matrix, "./ITHACAoutput/Matrices/",
                                  "RD_" + name(liftfield.size()) + "_" + 
                                  name(NUmodes) + "_" + name(NSUPmodes));
    return RD_matrix;
}

Eigen::MatrixXd steadyNS::boundary_vector_convection(label NUmodes, label NPmodes, label NSUPmodes)
{
    // Determine boundary vector
    fvVectorMatrix UEqn
    (
	fvm::div(fvc::flux(Uinl()),Uinl())
    );     
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Foam2Eigen::fvMatrix2Eigen(UEqn, A, b);

    label RCsize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd RC_matrix(RCsize,1);
    Eigen::VectorXd ModeVector;
 
    for (label i = 0; i < RCsize; i++)
    {
        ModeVector = Foam2Eigen::field2Eigen(L_U_SUPmodes[i]);
        RC_matrix(i,0) = ModeVector.dot(b.col(0));
    }

    if (Pstream::parRun())
    {
        reduce(RC_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(RC_matrix, "./ITHACAoutput/Matrices/",
                                  "RC_" + name(liftfield.size()) + "_" + 
                                  name(NUmodes) + "_" + name(NSUPmodes));
    return RC_matrix;
}

Eigen::Tensor<double, 3> steadyNS::convective_term_flux_tens(label NUmodes,
        label NPmodes, label NSUPmodes)
{
    dimensionedScalar dt_fake
    (
        "dt_fake",
        dimensionSet(0, 0, 1, 0, 0, 0, 0),
        scalar(1.0)
    );

    label Csize1 = NUmodes + NSUPmodes + liftfield.size();
    label Csize2 = NPmodes;
    Eigen::Tensor<double, 3> Cf_tensor;
    volVectorField L_U_SUPmodesaux(L_U_SUPmodes[0]);

    Cf_tensor.resize(Csize2, Csize1, Csize1);
    for (label i = 0; i < Csize2; i++)
    {
        for (label j = 0; j < Csize1; j++)
        {
            for (label k = 0; k < Csize1; k++)
            {
	        L_U_SUPmodesaux = dt_fake*(fvc::div(
                                      fvc::flux(L_U_SUPmodes[j]),
                                      L_U_SUPmodes[k]));

		Cf_tensor(i, j, k) = fvc::domainIntegrate(Pmodes[i] 
					* fvc::div(L_U_SUPmodesaux)).value();
	    }
	}
    }

    if (Pstream::parRun())
    {
        reduce(Cf_tensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    ITHACAstream::SaveDenseTensor(Cf_tensor, "./ITHACAoutput/Matrices/",
                                  "Cf_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_t");
    return Cf_tensor;
}

List<Eigen::MatrixXd> steadyNS::pressure_gradient_term_linsys_div(label NPmodes)
{
    dimensionedScalar dt_fake
    (
        "dt_fake",
        dimensionSet(0, 0, 1, 0, 0, 0, 0),
        scalar(1.0)
    );
    dimensionedScalar nu_fake
    (
        "nu_fake",
        dimensionSet(0, 2, -1, 0, 0, 0, 0),
        scalar(1.0)
    );
  
    volScalarField p(_p);
    fvScalarMatrix pEqn
    (
	fvm::laplacian(p) == (1/dt_fake)*fvc::div(Uinl())
    );
    pEqn.setReference(0, 0);
    LinSysDiv = Pmodes.project(pEqn,NPmodes);

    return LinSysDiv;
}

List<Eigen::MatrixXd> steadyNS::pressure_gradient_term_linsys_conv(label NPmodes)
{ 
    dimensionedScalar dt_fake
    (
        "dt_fake",
        dimensionSet(0, 0, 1, 0, 0, 0, 0),
        scalar(1.0)
    );
 
    volVectorField Caux(L_U_SUPmodes[0]);
    Caux = dt_fake * (-fvc::div(fvc::flux(Uinl()),Uinl()));
 
    volScalarField p(_p); 	
    fvScalarMatrix pEqn
    (
        fvm::laplacian(p) == (1/dt_fake)*fvc::div(Caux)
    );
    pEqn.setReference(0, 0);
    LinSysConv = Pmodes.project(pEqn,NPmodes);

    return LinSysConv;
}

List<Eigen::MatrixXd> steadyNS::pressure_gradient_term_linsys_diff(label NPmodes)
{
    dimensionedScalar nu_fake
    (
        "nu_fake",
        dimensionSet(0, 2, -1, 0, 0, 0, 0),
        scalar(1.0)
    );
    dimensionedScalar dt_fake
    (
        "dt_fake",
        dimensionSet(0, 0, 1, 0, 0, 0, 0),
        scalar(1.0)
    );

    volVectorField Daux(L_U_SUPmodes[0]);
    Daux = dt_fake * fvc::laplacian(nu_fake,Uinl());

    volScalarField p(_p);
    fvScalarMatrix pEqn
    (
         fvm::laplacian(p) == (1/dt_fake)*fvc::div(Daux)
    );
    pEqn.setReference(0, 0);
    LinSysDiff = Pmodes.project(pEqn,NPmodes);

    return LinSysDiff;
}

void steadyNS::change_viscosity(double mu)
{
    const volScalarField& nu =  _laminarTransport().nu();
    volScalarField& ciao = const_cast<volScalarField&>(nu);
    this->assignIF(ciao, mu);

    for (label i = 0; i < ciao.boundaryFieldRef().size(); i++)
    {
        this->assignBC(ciao, i, mu);
    }
}


void steadyNS::forcesMatrices(label NUmodes, label NPmodes, label NSUPmodes)
{
    tauMatrix.resize(L_U_SUPmodes.size(), 3);
    nMatrix.resize(NPmodes, 3);
    tauMatrix = tauMatrix * 0;
    nMatrix = nMatrix * 0;
    Time& runTime = _runTime();
    instantList Times = runTime.times();
    fvMesh& mesh = _mesh();
    volScalarField& p = _p();
    volVectorField& U = _U();
    //Read FORCESdict
    IOdictionary FORCESdict
    (
        IOobject
        (
            "FORCESdict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    word pName(FORCESdict.lookup("pName"));
    word UName(FORCESdict.lookup("UName"));
    functionObjects::ITHACAforces f("Forces", mesh, FORCESdict);

    for (label i = 0; i < L_U_SUPmodes.size(); i++)
    {
        U = L_U_SUPmodes[i];
        p = Pmodes[0];
        mesh.readUpdate();
        f.write();
        f.calcForcesMoment();

        for (label j = 0; j < 3; j++)
        {
            tauMatrix(i, j) = f.forceTau()[j];
        }
    }

    for (label i = 0; i < NPmodes; i++)
    {
        U = L_U_SUPmodes[0];
        p = Pmodes[i];
        mesh.readUpdate();
        f.write();
        f.calcForcesMoment();

        for (label j = 0; j < 3; j++)
        {
            nMatrix(i, j) = f.forcePressure()[j];
        }
    }

    if (Pstream::parRun())
    {
        reduce(tauMatrix, sumOp<Eigen::MatrixXd>());
    }

    if (Pstream::parRun())
    {
        reduce(nMatrix, sumOp<Eigen::MatrixXd>());
    }

    if (para->exportPython)
    {
        ITHACAstream::exportMatrix(tauMatrix, "tau", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(nMatrix, "n", "python", "./ITHACAoutput/Matrices/");
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportMatrix(tauMatrix, "tau", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(nMatrix, "n", "matlab", "./ITHACAoutput/Matrices/");
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportMatrix(tauMatrix, "tau", "eigen",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(nMatrix, "n", "eigen", "./ITHACAoutput/Matrices/");
    }
}

void steadyNS::forcesMatrices(label nModes)
{
    tauMatrix.resize(nModes, 3);
    nMatrix.resize(nModes, 3);
    tauMatrix = tauMatrix * 0;
    nMatrix = nMatrix * 0;
    Time& runTime = _runTime();
    instantList Times = runTime.times();
    fvMesh& mesh = _mesh();
    volScalarField& p = _p();
    volVectorField& U = _U();
    //Read FORCESdict
    IOdictionary FORCESdict
    (
        IOobject
        (
            "FORCESdict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    word pName(FORCESdict.lookup("pName"));
    word UName(FORCESdict.lookup("UName"));
    functionObjects::ITHACAforces f("Forces", mesh, FORCESdict);

    for (label i = 0; i < nModes; i++)
    {
        U = Umodes[i];
        p = Pmodes[0];
        mesh.readUpdate();
        f.write();
        f.calcForcesMoment();

        for (label j = 0; j < 3; j++)
        {
            tauMatrix(i, j) = f.forceTau()[j];
        }
    }

    for (label i = 0; i < nModes; i++)
    {
        U = Umodes[0];
        p = Pmodes[i];
        mesh.readUpdate();
        f.write();
        f.calcForcesMoment();

        for (label j = 0; j < 3; j++)
        {
            nMatrix(i, j) = f.forcePressure()[j];
        }
    }

    if (Pstream::parRun())
    {
        reduce(tauMatrix, sumOp<Eigen::MatrixXd>());
    }

    if (Pstream::parRun())
    {
        reduce(nMatrix, sumOp<Eigen::MatrixXd>());
    }

    if (para->exportPython)
    {
        ITHACAstream::exportMatrix(tauMatrix, "tau", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(nMatrix, "n", "python", "./ITHACAoutput/Matrices/");
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportMatrix(tauMatrix, "tau", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(nMatrix, "n", "matlab", "./ITHACAoutput/Matrices/");
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportMatrix(tauMatrix, "tau", "eigen",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(nMatrix, "n", "eigen", "./ITHACAoutput/Matrices/");
    }
}



void steadyNS::restart()
{
    _runTime().objectRegistry::clear();
    _mesh().objectRegistry::clear();
    // _mesh.clear();
    // _runTime.clear();
    _simple.clear();
    _p.clear();
    _U.clear();
    _phi.clear();
    turbulence.clear();
    _fvOptions.clear();
    argList& args = _args();
    Time& runTime = _runTime();
    runTime.setTime(0, 1);
    Foam::fvMesh& mesh = _mesh();
    _simple = autoPtr<simpleControl>
              (
                  new simpleControl
                  (
                      mesh
                  )
              );
    simpleControl& simple = _simple();
    Info << "ReReading field p\n" << endl;
    _p = autoPtr<volScalarField>
         (
             new volScalarField
             (
                 IOobject
                 (
                     "p",
                     runTime.timeName(),
                     mesh,
                     IOobject::MUST_READ,
                     IOobject::AUTO_WRITE
                 ),
                 mesh
             )
         );
    volScalarField& p = _p();
    Info << "ReReading field U\n" << endl;
    _U = autoPtr<volVectorField>
         (
             new volVectorField
             (
                 IOobject
                 (
                     "U",
                     runTime.timeName(),
                     mesh,
                     IOobject::MUST_READ,
                     IOobject::AUTO_WRITE
                 ),
                 mesh
             )
         );
    volVectorField& U = _U();
    Info << "ReReading/calculating face flux field phi\n" << endl;
    _phi = autoPtr<surfaceScalarField>
           (
               new surfaceScalarField
               (
                   IOobject
                   (
                       "phi",
                       runTime.timeName(),
                       mesh,
                       IOobject::READ_IF_PRESENT,
                       IOobject::AUTO_WRITE
                   ),
                   linearInterpolate(U) & mesh.Sf()
               )
           );
    surfaceScalarField& phi = _phi();
    pRefCell = 0;
    pRefValue = 0.0;
    setRefCell(p, simple.dict(), pRefCell, pRefValue);
    _laminarTransport = autoPtr<singlePhaseTransportModel>
                        (
                            new singlePhaseTransportModel( U, phi )
                        );
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    turbulence = autoPtr<incompressible::turbulenceModel>
                 (
                     incompressible::turbulenceModel::New(U, phi, laminarTransport)
                 );
    _MRF = autoPtr<IOMRFZoneList>
           (
               new IOMRFZoneList(mesh)
           );
    _fvOptions = autoPtr<fv::options>(new fv::options(mesh));
    turbulence->validate();
}
