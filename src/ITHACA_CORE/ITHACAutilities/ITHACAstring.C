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

#include "ITHACAstring.H"

/// \file
/// Source file of the ITHACAstring file.

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace ITHACAutilities
{

std::string str_trim(std::string const& s)
{
  auto const first{ s.find_first_not_of(' ') };
  if (first == std::string::npos) return {};
  auto const last{ s.find_last_not_of(' ') };
  return s.substr(first, (last - first + 1));
}

// io format function for files name
void str_format_io(std::string const& s, unsigned int nMax)
{
  if ( nMax > s.length() )
  {
    for (unsigned int n=0; n<s.length(); n++) std::cout << s[n];
    for (unsigned int n=s.length(); n<nMax; n++) std::cout << " ";
  }
  else
  {
    for (unsigned int n=0; n<nMax-3; n++) std::cout << s[n];
    for (unsigned int n=nMax-3; n<nMax; n++) std::cout << ".";
  }
}

std::string double2ConciseString(const double& d)
{
    std::string s = std::to_string(d);

    // Removing trailling zeros
    int dot_pos = s.find_first_of('.');
    if(dot_pos != std::string::npos)
    {
        int ipos = s.size()-1;
        while((s[ipos]=='0' || s[ipos]=='.') && ipos>dot_pos-1)
        {
            --ipos;
        }
        s.erase(ipos + 1, std::string::npos);
    }
  return s;
}

bool containsSubstring(std::string contain, std::string contained)
{
    std::transform(contain.begin(), contain.end(), contain.begin(), ::tolower);
    std::transform(contained.begin(), contained.end(), contained.begin(), ::tolower);
    return contain.find(contained) != std::string::npos;
}

std::vector<int> extractIntFromString(std::string input)
{
    std::stringstream ss;
    std::vector<int> numbers;
    int num;

    for (char c : input) 
    {
        if (isdigit(c)) ss << c;
        else ss << ' ';
    }
    while (ss >> num)
        numbers.push_back(num);
    return numbers;
}

}
