/*---------------------------------------------------------------------------*\
Copyright (C) 2017 by the ITHACA-FV authors

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
/// Template function file of the ITHACAstream class, it contains the implementation of
/// several methods for input output operations.

template<typename T>
void ITHACAstream::exportSolution(T& s, fileName subfolder, fileName folder, word fieldName)
{
    mkDir(folder+"/"+subfolder);
    T act(fieldName, s);
    fileName fieldname = folder+"/"+subfolder + "/" + fieldName;
    OFstream os(fieldname);
    act.writeHeader(os);
    os << act << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
