/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Author
    Frank Bos, TU Delft.  All rights reserved.
    Dubravko Matijasevic, FSB Zagreb.

\*---------------------------------------------------------------------------*/

#include "RBFInterpolation.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::scalarSquareMatrix& Foam::RBFInterpolation::B() const
{
    if (!BPtr_)
    {
        calcB();
    }

    return *BPtr_;
}

Foam::SquareMatrix<double> Foam::EigenInvert(Foam::SquareMatrix<double>& A)
{
    Foam::SquareMatrix<double> invMatrix = A;
    Eigen::MatrixXd Aeig(A.n(), A.n());
    Eigen::MatrixXd one = Eigen::MatrixXd::Identity(A.n(), A.n());

    for (int i = 0; i < A.n(); i++)
    {
        for (int k = 0; k < A.n(); k++)
        {
            Aeig(i, k) = A[i][k];
        }
    }

    //Eigen::MatrixXd invEig = Aeig.inverse();
    Eigen::MatrixXd invEig = Aeig.fullPivHouseholderQr().solve(one);

    for (int i = 0; i < invEig.rows(); i++)
    {
        for (int k = 0; k < invEig.cols(); k++)
        {
            invMatrix[i][k] = invEig(i, k);
        }
    }

    return invMatrix;
}

void Foam::RBFInterpolation::calcB() const
{
    // Determine inverse of boundary connectivity matrix
    label polySize(4);

    if (!polynomials_)
    {
        polySize = 0;
    }

    Eigen::MatrixXd Aeig = Eigen::MatrixXd::Zero(controlPoints_.size() + polySize,
                           controlPoints_.size() + polySize);
    const label nControlPoints = controlPoints_.size();

    for (label i = 0; i < nControlPoints; i++)
    {
        scalarField weights = RBF_->weights(controlPoints_, controlPoints_[i]);

        for (label col = 0; col < nControlPoints; col++)
        {
            Aeig(i, col) = weights[col];
        }
    }

    if (polynomials_)
    {
        for
        (
            label row = nControlPoints;
            row < nControlPoints + 1;
            row++
        )
        {
            for (label col = 0; col < nControlPoints; col++)
            {
                Aeig(col, row) = 1.0;
                Aeig(row, col) = 1.0;
            }
        }

        // Fill in X components of polynomial part of matrix
        for
        (
            label row = nControlPoints + 1;
            row < nControlPoints + 2;
            row++
        )
        {
            for (label col = 0; col < nControlPoints; col++)
            {
                Aeig(col, row) = controlPoints_[col].x();
                Aeig(row, col) = controlPoints_[col].x();
            }
        }

        // Fill in Y components of polynomial part of matrix
        for
        (
            label row = nControlPoints + 2;
            row < nControlPoints + 3;
            row++
        )
        {
            for (label col = 0; col < nControlPoints; col++)
            {
                Aeig(col, row) = controlPoints_[col].y();
                Aeig(row, col) = controlPoints_[col].y();
            }
        }

        // Fill in Z components of polynomial part of matrix
        for
        (
            label row = nControlPoints + 3;
            row < nControlPoints + 4;
            row++
        )
        {
            for (label col = 0; col < nControlPoints; col++)
            {
                Aeig(col, row) = controlPoints_[col].z();
                Aeig(row, col) = controlPoints_[col].z();
            }
        }

        // Fill 4x4 zero part of matrix
        for
        (
            label row = nControlPoints;
            row < nControlPoints + 4;
            row++
        )
        {
            for
            (
                label col = nControlPoints;
                col < nControlPoints + 4;
                col++
            )
            {
                Aeig(row, col) = 0.0;
            }
        }
    }

    Info << "Inverting RBF motion matrix" << endl;
    Eigen::MatrixXd InvAeig = Aeig.fullPivLu().inverse();
    simpleMatrix<scalar> InvA(controlPoints_.size() + polySize);

    for (int i = 0; i < InvAeig.rows(); i++)
    {
        for (int k = 0; k < InvAeig.cols(); k++)
        {
            InvA[i][k] = InvAeig(i, k);
        }
    }

    // HJ and FB (05 Jan 2009)
    // Collect ALL control points from ALL CPUs
    // Create an identical inverse for all CPUs
    BPtr_ = new scalarSquareMatrix(InvA);
}


void Foam::RBFInterpolation::clearOut()
{
    deleteDemandDrivenData(BPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBFInterpolation::RBFInterpolation
(
    const dictionary& dict,
    const vectorField& controlPoints,
    const vectorField& dataPoints
)
    :
    controlPoints_(controlPoints),
    dataPoints_(dataPoints),
    RBF_(RBFFunction::New(word(dict.lookup("RBF")), dict)),
    BPtr_(NULL),
    focalPoint_(dict.lookup("focalPoint")),
    innerRadius_(readScalar(dict.lookup("innerRadius"))),
    outerRadius_(readScalar(dict.lookup("outerRadius"))),
    polynomials_(dict.lookup("polynomials"))
{}


Foam::RBFInterpolation::RBFInterpolation
(
    const RBFInterpolation& rbf
)
    :
    controlPoints_(rbf.controlPoints_),
    dataPoints_(rbf.dataPoints_),
    RBF_(rbf.RBF_->clone()),
    BPtr_(NULL),
    focalPoint_(rbf.focalPoint_),
    innerRadius_(rbf.innerRadius_),
    outerRadius_(rbf.outerRadius_),
    polynomials_(rbf.polynomials_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBFInterpolation::~RBFInterpolation()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RBFInterpolation::movePoints()
{
    clearOut();
}


// ************************************************************************* //
