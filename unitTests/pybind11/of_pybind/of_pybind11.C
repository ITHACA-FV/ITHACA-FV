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
Description
    Example of a heat transfer Reduction Problem
SourceFiles
    02thermalBlock.C
\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "IOmanip.H"
#include "Time.H"
#include "laplacianProblem.H"
#include "ReducedLaplacian.H"
#include "ITHACAPOD.H"
#include "ITHACAutilities.H"
#include <cstddef>
#define _USE_MATH_DEFINES
#include <cmath>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>
#include "iostream"
#include "Foam2Eigen.H"

#if PY_VERSION_HEX < 0x03000000
#define MyPyText_AsString PyString_AsString
#else
#define MyPyText_AsString PyUnicode_AsUTF8
#endif

namespace py = pybind11;

/// \brief Class where the tutorial number 2 is implemented.
/// \details It is a child of the laplacianProblem class and some of its
/// functions are overridden to be adapted to the specific case.
class of_pybind11 : public laplacianProblem
{
public:
    of_pybind11(int argc, char* argv[])
    {
        _args = autoPtr<argList>(
            new argList(argc, argv, true, true, /*initialise=*/false));
        argList& args = _args();
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"
    }
    ~of_pybind11() {};
    Eigen::Map<Eigen::MatrixXd> getT()
    {
        volScalarField& T = _T();
        Eigen::Map<Eigen::MatrixXd> Teig(Foam2Eigen::field2EigenMap(T));
        return std::move(Teig);
    }

    void printT()
    {
        volScalarField& T = _T();
        Info << T << endl;
    }
};


PYBIND11_MODULE(of_pybind11, m)
{
    // bindings to Matrix class
    py::class_<of_pybind11>(m, "of_pybind11")
        .def(py::init([](
                          std::vector<std::string> args) {
            std::vector<char*> cstrs;
            cstrs.reserve(args.size());
            for (auto& s : args)
                cstrs.push_back(const_cast<char*>(s.c_str()));
            return new of_pybind11(cstrs.size(), cstrs.data());
        }),
        py::arg("args") = std::vector<std::string> { "." })
        .def("getT", &of_pybind11::getT, py::return_value_policy::reference_internal)
        .def("printT", &of_pybind11::printT);
}
