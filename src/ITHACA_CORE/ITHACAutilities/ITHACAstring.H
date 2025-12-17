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

#ifndef ITHACAstring_H
#define ITHACAstring_H

/// \file
/// Header file of the ITHACAstring file.

#include "ITHACAparameters.H"



namespace ITHACAutilities
{

//------------------------------------------------------------------------------
/// @brief      Trim a string cf. 
///             https://stackoverflow.com/questions/25829143/trim-whitespace-from-a-string
///
/// @param      String s
std::string str_trim(std::string const& s);

//------------------------------------------------------------------------------
/// @brief      Format a string
///
/// @param      s string
/// @param      Nmax max length of the string
void str_format_io( std::string const& s, unsigned int Nmax );

//------------------------------------------------------------------------------
/// @brief      Returns the double d in string format without keeping all the extra 0
///
/// @param      d double
std::string double2ConciseString(const double& d);

//------------------------------------------------------------------------------
/// @brief      Returns true if contained is a substring of contain, false otherwise
///
/// @param      contain string
/// @param      contained string
bool containsSubstring(std::string contain, std::string contained);

//------------------------------------------------------------------------------
/// @brief      Returns an array storing all the integers in the input string
///
/// @param      input string
std::vector<int> extractIntFromString(std::string input);

};

#endif
