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
Class
    baseFilter
Description
    set of functions to convert Eigen object to Torch objects and viceversa
SourceFiles
    baseFilter.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the BaseFilter class.
///
///

#ifndef FILTER_H
#define FILTER_H
#include <torch/script.h>
#include <torch/torch.h>
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "tmp.H"
#include "autoPtr.H"
#include "primitiveFields.H"
#include "fvCFD.H"


namespace Foam
{

class Filter
{
        //- Disallow copy construct
        Filter(const Filter&);

        //- Disallow default bitwise assignment
        void operator=(const Filter&);

    public:

        //- Runtime type information
        TypeName("Filter");

        declareRunTimeSelectionTable
        (
            autoPtr,
            Filter,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );

        //- Construct null
        Filter()
        {};

        //- Create and return a clone
        virtual autoPtr<Filter> clone() const = 0;

        // Destructor

        virtual ~Filter()
        {};

        virtual List<scalarList> apply(const List<labelList>& cells, const List<point>& convPoints, const fvMesh& mesh) const = 0;

        static autoPtr<Filter> New
        (
            const word& type,
            const dictionary& dict
        );

};

}

#endif
