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
/// Header file of the Burgers class.
/// \dir
/// Directory containing the header and source files for the Burgers class.

#ifndef burgers_H
#define burgers_H
#include "fvCFD.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "reductionProblem.H"
#include "UnsteadyProblem.H"
#include "ITHACAstream.H"
#include "ITHACAparameters.H"
#include "Modes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class burgers Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a parametrized full order <b> Burgers </b> and preparation of the the reduced matrices for the online solve.
/** In this class are implemented the methods for the offline solve of a burgers problem
and the for the generation of the reduced matrices for subsequent online solve, this class is a son
of the reduction problem class */
class Burgers: public reductionProblem, public UnsteadyProblem
{

    public:
        // Constructors
        /// Null constructor
        Burgers();

        /// Construct with argc and argv
        Burgers(int argc, char* argv[]);

        ~Burgers() {};

        ITHACAparameters* para;

        /// List of pointers used to form the velocity snapshots matrix
        PtrList<volVectorField> Ufield;

        /// List of pointers used to form the velocity modes
        volVectorModes Umodes;

        /// Velocity field
        autoPtr<volVectorField> _U;

        /// Initial Velocity field (for restart purposes)
        autoPtr<volVectorField> _U0;

        /// Mesh
        mutable autoPtr<fvMesh> _mesh;

        /// Time
        autoPtr<Time> _runTime;

        /// Flux
        autoPtr<surfaceScalarField> _phi;

        /// Initial Flux (for restart purposes)
        autoPtr<surfaceScalarField> _phi0;

        autoPtr<IOdictionary> _transportProperties;

        autoPtr<dimensionedScalar> _nu;

        /// fvOptions
        autoPtr<fv::options> _fvOptions;

        /// simpleControl
        autoPtr<simpleControl> _simple;

        /// truthsolve
        void truthSolve(word folder = "./ITHACAoutput/Offline/");

        //----------------------------------------------------------------------------
        /// @brief      Function to restart the fields of the Burgers problem
        ///
        void restart();

};

#endif
