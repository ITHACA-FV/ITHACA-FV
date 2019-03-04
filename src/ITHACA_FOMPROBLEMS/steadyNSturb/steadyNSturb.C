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
#include "steadyNSturb.H"
#include "viscosityModel.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
steadyNSturb::steadyNSturb() {}
steadyNSturb::steadyNSturb(int argc, char* argv[])
{
    // #include "postProcess.H"
#include "setRootCase.H"
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
    supex = ITHACAutilities::check_sup();
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
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    supex = ITHACAutilities::check_sup();
}


// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //

// Method to performa a truthSolve
void steadyNSturb::truthSolve(List<scalar> mu_now)
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
    exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
    volScalarField _nut(turbulence->nut());
    exportSolution(_nut, name(counter), "./ITHACAoutput/Offline/");
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

List < Eigen::MatrixXd > steadyNSturb::turbulence_term1(label NUmodes,
        label NSUPmodes, label Nnutmodes)
{
    label Csize = NUmodes + NSUPmodes + liftfield.size();
    List < Eigen::MatrixXd > CT1_matrix;
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







List < Eigen::MatrixXd > steadyNSturb::turbulence_term2(label NUmodes,
        label NSUPmodes, label Nnutmodes)
{
    label Csize = NUmodes + NSUPmodes + liftfield.size();
    List < Eigen::MatrixXd > CT2_matrix;
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

Eigen::MatrixXd steadyNSturb::BT_turbulence(label NUmodes, label NSUPmodes)
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

void steadyNSturb::projectSUP(fileName folder, label NU, label NP, label NSUP,
                              label Nnut)
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
    }
    else
    {
        NUmodes = NU;
        NPmodes = NP;
        NSUPmodes = NSUP;
        Nnutmodes = Nnut;
        B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        C_matrix = convective_term(NUmodes, NPmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        P_matrix = divergence_term(NUmodes, NPmodes, NSUPmodes);
        M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
        BT_matrix = BT_turbulence(NUmodes, NSUPmodes);
        CT1_matrix = turbulence_term1(NUmodes, NSUPmodes, Nnutmodes);
        CT2_matrix = turbulence_term2(NUmodes, NSUPmodes, Nnutmodes);
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
    Eigen::MatrixXd Ncoeff = ITHACAutilities::get_coeffs_ortho(nutFields, nuTmodes);
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
