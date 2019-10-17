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
/// Source file of the reductionProblem class.


#include "reductionProblem.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
reductionProblem::reductionProblem()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Set the parameters
void reductionProblem::setParameters()
{
    mu.resize(Pnumber, Tnumber);
    mu_range.resize(Pnumber, 2);
}

// Generate Random Parameters
void reductionProblem::genRandPar()
{
    for (int k = 0; k < Pnumber; k++)
    {
        for (int n = 0; n < Tnumber; ++n)
        {
            mu(k, n) = mu_range(k, 0) + static_cast <float> (rand()) / static_cast <float>
                       (RAND_MAX / (mu_range(k, 1) - mu_range(k, 0))); // generate numbers
        }
    }
}

void reductionProblem::genRandPar(int Tsize)
{
    mu.resize(Tsize, Pnumber);
    std::srand(std::time(0));
    auto rand = Eigen::MatrixXd::Random(Tsize, Pnumber);
    auto dx = mu_range.col(1) - mu_range.col(0);
    Eigen::MatrixXd k;
    k = (rand.array() + double(1)) / 2;

    for (int i = 0; i < Pnumber; i ++)
    {
        k.col(i) = (k.col(i) * dx(i)).array() + mu_range(i, 0);
    }

    mu = k;
}

// Generate Equidistributed Parameters
void reductionProblem::genEquiPar()
{
    Eigen::VectorXd vec;

    for (int k = 0; k < Pnumber; k++)
    {
        mu.row(k) = vec.LinSpaced(Tnumber, mu_range(k, 0), mu_range(k, 1));
    }
}

// Perform a TruthSolve (To be overridden)
void reductionProblem::truthSolve()
{
    Info << "reductionProblem::truthSolve() is a Method to be overridden -> Exiting the code"
         << endl;
    exit(0);
}

// Assign a BC for a vector field
void reductionProblem::assignBC(volScalarField& s, label BC_ind, double& value)
{
    if (s.boundaryField()[BC_ind].type() == "fixedValue"
            || s.boundaryField()[BC_ind].type() == "processor")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "fixedGradient")
    {
        fixedGradientFvPatchScalarField& Tpatch =
            refCast<fixedGradientFvPatchScalarField>(s.boundaryFieldRef()[BC_ind]);
        scalarField& gradTpatch = Tpatch.gradient();
        forAll(gradTpatch, faceI)
        {
            gradTpatch[faceI] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "fixedFluxPressure")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "freestream")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }

        freestreamFvPatchField<scalar>& Tpatch =
            refCast<freestreamFvPatchField<scalar>>(s.boundaryFieldRef()[BC_ind]);
        scalarField& gradTpatch = Tpatch.freestreamValue();
        forAll(gradTpatch, faceI)
        {
            gradTpatch[faceI] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "calculated")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "empty")
    {
    }
}

// Assign a BC for a scalar field
void reductionProblem::assignBC(volVectorField& s, label BC_ind,
                                Vector<double>& value)
{
    if (s.boundaryField()[BC_ind].type() == "fixedValue")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "fixedGradient")
    {
        Info << "This Feature is not implemented for this boundary condition" << endl;
        exit(0);
    }
    else if (s.boundaryField()[BC_ind].type() == "freestream")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }

        freestreamFvPatchField<vector>& Tpatch =
            refCast<freestreamFvPatchField<vector>>(s.boundaryFieldRef()[BC_ind]);
        vectorField& gradTpatch = Tpatch.freestreamValue();
        forAll(gradTpatch, faceI)
        {
            gradTpatch[faceI] = value;
        }
    }
}

// Reconstruct using a Matrix of coefficients (vector field)
void reductionProblem::reconstructFromMatrix(PtrList<volVectorField>&
        rec_field2, PtrList<volVectorField>& modes, label Nmodes,
        Eigen::MatrixXd coeff_matrix)
{
    for (label k = 0; k < coeff_matrix.cols(); k++)
    {
        for (label i = 0; i < Nmodes; i++)
        {
            if ( i == 0)
            {
                rec_field2.append(coeff_matrix(i, k)*modes[i]);
            }
            else
            {
                rec_field2[k] += coeff_matrix(i, k) * modes[i];
            }
        }
    }
}


// Reconstruct using a Matrix of coefficients (vector field)
void reductionProblem::reconstructFromMatrix(PtrList<volScalarField>&
        rec_field2, PtrList<volScalarField>& modes, label Nmodes,
        Eigen::MatrixXd coeff_matrix)
{
    for (label k = 0; k < coeff_matrix.cols(); k++)
    {
        for (label i = 0; i < Nmodes; i++)
        {
            if ( i == 0)
            {
                rec_field2.append(coeff_matrix(i, k)*modes[i]);
            }
            else
            {
                rec_field2[k] += coeff_matrix(i, k) * modes[i];
            }
        }
    }
}

void reductionProblem::project()
{
    Info << "reductionProblem::project is a function to be overridden exiting the code"
         << endl;
    exit(0);
}

void reductionProblem::writeMu(List<scalar> mu_now)
{
    mkDir("./ITHACAoutput/Parameters");
    std::ofstream ofs;
    ofs.open ("./ITHACAoutput/Parameters/par",
              std::ofstream::out | std::ofstream::app);
    forAll(mu_now, i)
    {
        ofs << mu_now[i] << ' ';
    }
    ofs << "\n";
    ofs.close();
}

void reductionProblem::liftSolve()
{
    Info << "reductionProblem::liftSolve is a virtual function it must be overridden"
         << endl;
    exit(0);
}

void reductionProblem::liftSolveT()
{
    Info << "reductionProblem::liftSolveT is a virtual function it must be overridden"
         << endl;
    exit(0);
}

std::vector<SPLINTER::RBFSpline> reductionProblem::getCoeffManifoldRBF(
    PtrList<volVectorField> snapshots, PtrList<volVectorField>& modes,
    word rbfBasis)
{
    Eigen::MatrixXd coeff = ITHACAutilities::get_coeffs_ortho(snapshots, modes);
    M_Assert(mu_samples.rows() == coeff.cols() * mu.rows(),
             "The dimension of the coefficient matrix must correspond to the constructed parameter matrix 'mu_samples'");
    std::vector<SPLINTER::RBFSpline> rbfsplines;
    rbfsplines.reserve(coeff.rows());
    // Check that each column of matrix "mu_samples" does not satisfy the case where all elements hold same value (which might appear for e.g. unsteady non-parametric problems)
    // In which case the spline builder uses a refined matrix "muSampRefined" that discards those columns aforementioned, hence allows for a feasible manifold construction
    Eigen::MatrixXd muSampRefined(mu_samples.rows(), 0);

    for (label i = 0; i < mu_samples.cols(); i++)
    {
        double firstVal = mu_samples(0, i);
        bool isConst = (mu_samples.col(i).array() == firstVal).all();

        // Consider only columns with non-constant values
        if (isConst == 0)
        {
            muSampRefined.conservativeResize(muSampRefined.rows(),
                                             muSampRefined.cols() + 1);
            muSampRefined.col(muSampRefined.cols() - 1) = mu_samples.col(i);
        }
    }

    // --- Loop over Nmodes
    for (label i = 0; i < coeff.rows(); i++)
    {
        SPLINTER::DataTable samples;

        // --- Loop over Nsamples
        for (label j = 0; j < mu_samples.rows(); j++)
        {
            SPLINTER::DenseVector x(muSampRefined.cols());

            // --- Loop over Nparameters
            for (label k = 0; k < muSampRefined.cols(); k++)
            {
                x(k) = muSampRefined(j, k);
            }

            samples.addSample(x, coeff(i, j));
        }

        if (rbfBasis == "GAUSSIAN")
        {
            rbfsplines.push_back(SPLINTER::RBFSpline(samples,
                                 SPLINTER::RadialBasisFunctionType::GAUSSIAN));
        }
        else if (rbfBasis == "THIN_PLATE")
        {
            rbfsplines.push_back(SPLINTER::RBFSpline(samples,
                                 SPLINTER::RadialBasisFunctionType::THIN_PLATE_SPLINE));
        }
        else if (rbfBasis == "MULTI_QUADRIC")
        {
            rbfsplines.push_back(SPLINTER::RBFSpline(samples,
                                 SPLINTER::RadialBasisFunctionType::MULTIQUADRIC));
        }
        else if (rbfBasis == "INVERSE_QUADRIC")
        {
            rbfsplines.push_back(SPLINTER::RBFSpline(samples,
                                 SPLINTER::RadialBasisFunctionType::INVERSE_QUADRIC));
        }
        else if (rbfBasis == "INVERSE_MULTI_QUADRIC")
        {
            rbfsplines.push_back(SPLINTER::RBFSpline(samples,
                                 SPLINTER::RadialBasisFunctionType::INVERSE_MULTIQUADRIC));
        }
        else
        {
            std::cout <<
                      "Unknown string for rbfBasis. Valid types are 'GAUSSIAN', 'THIN_PLATE', 'MULTI_QUADRIC', 'INVERSE_QUADRIC', 'INVERSE_MULTI_QUADRIC'"
                      << std::endl;
            exit(0);
        }
    }

    return rbfsplines;
}


std::vector<SPLINTER::RBFSpline> reductionProblem::getCoeffManifoldRBF(
    PtrList<volScalarField> snapshots, PtrList<volScalarField>& modes,
    word rbfBasis)
{
    Eigen::MatrixXd coeff = ITHACAutilities::get_coeffs_ortho(snapshots, modes);
    M_Assert(mu_samples.rows() == coeff.cols() * mu.rows(),
             "The dimension of the coefficient matrix must correspond to the constructed parameter matrix 'mu_samples'");
    std::vector<SPLINTER::RBFSpline> rbfsplines;
    rbfsplines.reserve(coeff.rows());
    // Check that each column of matrix "mu_samples" does not satisfy the case where all elements hold same value (which might appear for e.g. unsteady non-parametric problems)
    // In which case the spline builder uses a refined matrix "muSampRefined" that discards those columns aforementioned, hence allows for a feasible manifold construction
    Eigen::MatrixXd muSampRefined(mu_samples.rows(), 0);

    for (label i = 0; i < mu_samples.cols(); i++)
    {
        double firstVal = mu_samples(0, i);
        bool isConst = (mu_samples.col(i).array() == firstVal).all();

        // Consider only columns with non-constant values
        if (isConst == 0)
        {
            muSampRefined.conservativeResize(muSampRefined.rows(),
                                             muSampRefined.cols() + 1);
            muSampRefined.col(muSampRefined.cols() - 1) = mu_samples.col(i);
        }
    }

    // --- Loop over Nmodes
    for (label i = 0; i < coeff.rows(); i++)
    {
        SPLINTER::DataTable samples;

        // --- Loop over Nsamples
        for (label j = 0; j < mu_samples.rows(); j++)
        {
            SPLINTER::DenseVector x(muSampRefined.cols());

            // --- Loop over Nparameters
            for (label k = 0; k < muSampRefined.cols(); k++)
            {
                x(k) = muSampRefined(j, k);
            }

            samples.addSample(x, coeff(i, j));
        }

        if (rbfBasis == "GAUSSIAN")
        {
            rbfsplines.push_back(SPLINTER::RBFSpline(samples,
                                 SPLINTER::RadialBasisFunctionType::GAUSSIAN));
        }
        else if (rbfBasis == "THIN_PLATE")
        {
            rbfsplines.push_back(SPLINTER::RBFSpline(samples,
                                 SPLINTER::RadialBasisFunctionType::THIN_PLATE_SPLINE));
        }
        else if (rbfBasis == "MULTI_QUADRIC")
        {
            rbfsplines.push_back(SPLINTER::RBFSpline(samples,
                                 SPLINTER::RadialBasisFunctionType::MULTIQUADRIC));
        }
        else if (rbfBasis == "INVERSE_QUADRIC")
        {
            rbfsplines.push_back(SPLINTER::RBFSpline(samples,
                                 SPLINTER::RadialBasisFunctionType::INVERSE_QUADRIC));
        }
        else if (rbfBasis == "INVERSE_MULTI_QUADRIC")
        {
            rbfsplines.push_back(SPLINTER::RBFSpline(samples,
                                 SPLINTER::RadialBasisFunctionType::INVERSE_MULTIQUADRIC));
        }
        else
        {
            std::cout <<
                      "Unknown string for rbfBasis. Valid types are 'GAUSSIAN', 'THIN_PLATE', 'MULTI_QUADRIC', 'INVERSE_QUADRIC', 'INVERSE_MULTI_QUADRIC'"
                      << std::endl;
            exit(0);
        }
    }

    return rbfsplines;
}


std::vector<SPLINTER::BSpline> reductionProblem::getCoeffManifoldSPL(
    PtrList<volVectorField> snapshots, PtrList<volVectorField>& modes, label splDeg)
{
    Eigen::MatrixXd coeff = ITHACAutilities::get_coeffs_ortho(snapshots, modes);
    M_Assert(mu_samples.rows() == coeff.cols() * mu.rows(),
             "The dimension of the coefficient matrix must correspond to the constructed parameter matrix 'mu_samples'");
    std::vector<SPLINTER::BSpline> bsplines;
    bsplines.reserve(coeff.rows());
    // Check that each column of matrix "mu_samples" does not satisfy the case where all elements hold same value (which might appear for e.g. unsteady non-parametric problems)
    // In which case the spline builder uses a refined matrix "muSampRefined" that discards those columns aforementioned, hence allows for a feasible manifold construction
    Eigen::MatrixXd muSampRefined(mu_samples.rows(), 0);

    for (label i = 0; i < mu_samples.cols(); i++)
    {
        double firstVal = mu_samples(0, i);
        bool isConst = (mu_samples.col(i).array() == firstVal).all();

        // Consider only columns with non-constant values
        if (isConst == 0)
        {
            muSampRefined.conservativeResize(muSampRefined.rows(),
                                             muSampRefined.cols() + 1);
            muSampRefined.col(muSampRefined.cols() - 1) = mu_samples.col(i);
        }
    }

    // --- Loop over Nmodes
    for (label i = 0; i < coeff.rows(); i++)
    {
        SPLINTER::DataTable samples;

        // --- Loop over Nsamples
        for (label j = 0; j < mu_samples.rows(); j++)
        {
            SPLINTER::DenseVector x(muSampRefined.cols());

            // --- Loop over Nparameters
            for (label k = 0; k < muSampRefined.cols(); k++)
            {
                x(k) = muSampRefined(j, k);
            }

            samples.addSample(x, coeff(i, j));
        }

        SPLINTER::BSpline bspline = SPLINTER::BSpline::Builder(samples).degree(
                                        splDeg).build();
        bsplines.push_back(bspline);
    }

    return bsplines;
}


std::vector<SPLINTER::BSpline> reductionProblem::getCoeffManifoldSPL(
    PtrList<volScalarField> snapshots, PtrList<volScalarField>& modes, label splDeg)
{
    Eigen::MatrixXd coeff = ITHACAutilities::get_coeffs_ortho(snapshots, modes);
    M_Assert(mu_samples.rows() == coeff.cols() * mu.rows(),
             "The dimension of the coefficient matrix must correspond to the constructed parameter matrix 'mu_samples'");
    std::vector<SPLINTER::BSpline> bsplines;
    bsplines.reserve(coeff.rows());
    // Check that each column of matrix "mu_samples" does not satisfy the case where all elements hold same value (which might appear for e.g. unsteady non-parametric problems)
    // In which case the spline builder uses a refined matrix "muSampRefined" that discards those columns aforementioned, hence allows for a feasible manifold construction
    Eigen::MatrixXd muSampRefined(mu_samples.rows(), 0);

    for (label i = 0; i < mu_samples.cols(); i++)
    {
        double firstVal = mu_samples(0, i);
        bool isConst = (mu_samples.col(i).array() == firstVal).all();

        // Consider only columns with non-constant values
        if (isConst == 0)
        {
            muSampRefined.conservativeResize(muSampRefined.rows(),
                                             muSampRefined.cols() + 1);
            muSampRefined.col(muSampRefined.cols() - 1) = mu_samples.col(i);
        }
    }

    // --- Loop over Nmodes
    for (label i = 0; i < coeff.rows(); i++)
    {
        SPLINTER::DataTable samples;

        // --- Loop over Nsamples
        for (label j = 0; j < mu_samples.rows(); j++)
        {
            SPLINTER::DenseVector x(muSampRefined.cols());

            // --- Loop over Nparameters
            for (label k = 0; k < muSampRefined.cols(); k++)
            {
                x(k) = muSampRefined(j, k);
            }

            samples.addSample(x, coeff(i, j));
        }

        SPLINTER::BSpline bspline = SPLINTER::BSpline::Builder(samples).degree(
                                        splDeg).build();
        bsplines.push_back(bspline);
    }

    return bsplines;
}

// ************************************************************************* //

