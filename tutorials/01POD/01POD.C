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

/// \dir 01POD Folder of the turorial 1
/// \file 01POD/01POD.C fake C file for documentation of tutorial 01
/// 
/// \example 01POD.C
/// 
/// \section intro Introduction
/// @brief The perform_POD application is used to perform POD on a standard OpenFOAM case.
/// @details In this tutorial the POD is applied on a simple lid-driven cavity example. 
/// 
/// The case is solved using icoFoam and then the modes are extracted using the perform_POD utility. Run
/// the code using the standard icoFoam solver from OpenFOAM
/// 
/// \code{.sh}
/// icoFoam
/// \endcode
///  
/// Once the solution is performed you can the ITHACA-FV utility to extract the modes running:
/// 
/// \code{.sh}
/// perform_POD
/// \endcode
/// 
/// The modes, together with the eigenvalues and the cumulative eigenvalues are stored inside the ITHACAoutput/POD.
/// folder.
/// 
/// You can eventually use also the Allrun file to perform both operation together. For more details see also the perform_POD.C file 
/// and check the system/ITHACAPODdict file in the system folder to understand how this file should be prepared.
/// 
/// \section code The ITHACAPODdict file
/// 
/// The ITHACAPODdict file located in the system folder is used to define the characteristics of the fields on which 
/// you want to perform the POD and how many modes you want to extract. Let's have a detailed look to it.
/// 
/// \dontinclude ./01POD/system/ITHACAPODdict
/// 
/// With the following lines you decide on which fields you want to perform the POD, put the name
/// of the field followed by the _pod string:
/// 
/// \skip subdictionary
/// \skipline fields
/// \until )
/// 
/// Then for each field we have to specify the exact file name of the field, how many modes we want to extract
/// and the type of the field (vector or scalar).
/// In this case we have 10 modes for both velocity and pressure field:
/// 
/// \skipline U_pod
/// \until }
/// 
/// \skipline p_pod
/// \until }
/// 
/// Then you have to select the initialTime from which you want to start acquiring the snapshots
/// and the FinalTime
/// 
/// \skipline // Set
/// \until 0;
/// \skipline // Set
/// \until 50;
/// 
/// Eventually, instead of the finalTime you could define the number of snapshots
/// 
/// \skipline // Even
/// \skipline //N
/// 
/// \section plain The plain ITHACAPODdict dictionary
/// \include ./01POD/system/ITHACAPODdict

// a fake class just to print out the doxygen documentation
class fake