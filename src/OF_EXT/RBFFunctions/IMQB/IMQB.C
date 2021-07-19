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

\*---------------------------------------------------------------------------*/

#include "IMQB.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(IMQB, 0);
addToRunTimeSelectionTable(RBFFunction, IMQB, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IMQB::IMQB(const scalar radius)
    :
    RBFFunction(),
    radius_(radius)
{}


// Construct from dictionary
Foam::IMQB::IMQB(const dictionary& dict)
    :
    RBFFunction(),
    radius_(readScalar(dict.lookup("radius")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::IMQB::~IMQB()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField Foam::IMQB::weights
(
    const vectorField& controlPoints,
    const vector& dataPoint
) const
{
    // Algorithmic improvement, Matteo Lombardi.  21/Mar/2011
    scalarField sqrDist = magSqr(controlPoints - dataPoint).ref();
    return (1 / sqrt(sqrDist + sqr(radius_))).ref();
}


// ************************************************************************* //
