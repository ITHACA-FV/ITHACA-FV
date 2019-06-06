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
#include "SteadyNSTurb.H"
#include "viscosityModel.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
SteadyNSTurb::SteadyNSTurb() {}
SteadyNSTurb::SteadyNSTurb(int argc, char* argv[])
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
    bcMethod = ITHACAdict->lookupOrDefault<word>("BCMethod", "lift");
    M_Assert(bcMethod == "lift" || bcMethod == "penalty",
             "The BC method must be set to lift or penalty in ITHACAdict");
    viscCoeff = ITHACAdict->lookupOrDefault<word>("viscCoeff", "RBF");
    para = new ITHACAparameters;
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    supex = ITHACAutilities::check_sup();
}


// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //

// Method to performa a truthSolve
void SteadyNSTurb::truthSolve(List<scalar> mu_now)
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

List < Eigen::MatrixXd > SteadyNSTurb::turbulenceTerm1(label NUmodes,
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







List < Eigen::MatrixXd > SteadyNSTurb::turbulenceTerm2(label NUmodes,
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

Eigen::MatrixXd SteadyNSTurb::btTurbulence(label NUmodes, label NSUPmodes)
{
    label btSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd btMatrix(btSize, btSize);
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
    for (label i = 0; i < btSize; i++)
    {
        for (label j = 0; j < btSize; j++)
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

void SteadyNSTurb::projectSUP(fileName folder, label NU, label NP, label NSUP,
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
    coeffL2 = ITHACAutilities::get_coeffs_ortho(nutFields, nutModes);
    ITHACAstream::exportMatrix(coeffL2, "coeffL2", "python",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(coeffL2, "coeffL2", "matlab",
                               "./ITHACAoutput/Matrices/");
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
