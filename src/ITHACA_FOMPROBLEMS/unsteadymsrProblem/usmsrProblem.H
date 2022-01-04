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

#ifndef usmsrProblem_H
#define usmsrProblem_H
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "IOporosityModelList.H"
#include "IOMRFZoneList.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "../msrProblem/msrProblem.H"
#include "UnsteadyProblem.H"
#include <iostream>

// us stands for unsteady
class usmsrProblem: public msrProblem, public UnsteadyProblem
{
    public:
        // Constructors
        /// Construct Null
        usmsrProblem();

        /// Construct with argc and argv
        usmsrProblem(int argc, char* argv[]);

        // Specific variable for the unstationary case
        autoPtr<dimensionedScalar> _tau;

        // Dummy variables to transform pimplefoam into a class
        /// pimpleControl
        autoPtr<pimpleControl> _pimple;
        autoPtr<pimpleControl> _npimple;

        /// Turbulence model
        autoPtr<incompressible::turbulenceModel> turbulence;

        /// adjustTimeStep
        bool adjustTimeStep;

        /// maxCourant
        scalar maxCo;

        /// maxDeltaT
        scalar maxDeltaT;



        //Methods:
        //--------------------------------------------------------------------------
        /// Function to check if the solution must be exported.
        ///
        /// @param      timeObject  The time object of OpenFOAM.
        ///
        /// @return     1 if we must write 0 elsewhere.
        ///
        bool checkWrite(Time& timeObject);

        //--------------------------------------------------------------------------
        /// @brief      Perform a truthsolve
        /// @param[in]  mu_now  The actual value of the parameter for this truthSolve. Used only
        /// to construct mu_interp matrix which is written out in the Offline folder, also for par file in
        /// the Parameters folder.
        ///
        void truthSolve(List<scalar> mu_now);
        //--------------------------------------------------------------------------
        /// @brief      Perform a truthsolve exporting the solution in std::string folder
        /// @param[in]  mu_now  The actual value of the parameter for this truthSolve. Used only
        /// to construct mu_interp matrix which is written out in the Offline folder, also for par file in
        /// the Parameters folder.
        void truthSolve(List<scalar> mu_now, std::string folder);

    protected:

        /// Method to change the precursors' "precinIndex" (i.e. the
        /// inlet) patch  type from any to fixedValue
        void changePrecsBC();
        //---------------------------------------------------------------

        /// Method to compute the precursors inlet bc as:
        /// prec_i(x=outlet,t)exp(-lam_i*tau)
        /// where prec_i(x=outlet,t-tau) is the averege over the outlet
        /// patch identified by "precoutIndex"
        /// values are saved in bc_prec(i,call)
        void computePrecsBC(label call);
        //---------------------------------------------------------------

        /// Method to assign the inlet BC to all precursors
        /// it assigns prec_i(x=inlet,t)=prec_i(x=outlet,t-tau)exp(-lam_i*tau)
        /// in practice bc_prec(i,call-Ntau) is assigned to prec_i
        /// patch identified by "precinIndex"
        void assignPrecsBC(label call, label Ntau);
        //---------------------------------------------------------------

        ///note: these methods are protected since they are very peculiar and problem
        ///dependent, in this way the user cannot use it in the main of the programm causing
        ///errors.

};

#endif



