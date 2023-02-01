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

class of_pybind11_system : public laplacianProblem
{
public:
    of_pybind11_system(int argc, char* argv[])
    {
        _args = autoPtr<argList>(
            new argList(argc, argv, true, true, /*initialise=*/false));
        argList& args = _args();
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"
    }
    ~of_pybind11_system() {};
    Eigen::Map<Eigen::MatrixXd> getT()
    {
        Eigen::Map<Eigen::MatrixXd> Teig(Foam2Eigen::field2EigenMap(_T()));
        return std::move(Teig);
    }

    Eigen::Map<Eigen::MatrixXd> getS()
    {
        Eigen::Map<Eigen::MatrixXd> Teig(Foam2Eigen::field2EigenMap(_S()));
        return std::move(Teig);
    }

    Eigen::Map<Eigen::MatrixXd> get_residual(Eigen::VectorXd& T, Eigen::VectorXd& S)
    {
        Foam2Eigen::Eigen2field(_T(), T);
        Foam2Eigen::Eigen2field(_S(), S);
        _res() = fvc::laplacian(_nu(), _T()) - _S();
        Eigen::Map<Eigen::MatrixXd> res_eig(Foam2Eigen::field2EigenMap(_res()));
        return std::move(res_eig);
    }

    void printT()
    {
        Info << _T() << endl;
    }

    void printS()
    {
        Info << _S() << endl;
    }

    void printMatrix(Eigen::Ref<Eigen::MatrixXd>& a)
    {
        std::cerr << a << std::endl;
    }

    Eigen::SparseMatrix<double> get_system_matrix(Eigen::VectorXd& T, Eigen::VectorXd& S)
    {
        Foam2Eigen::Eigen2field(_T(), T);
        Foam2Eigen::Eigen2field(_S(), S);
        fvScalarMatrix TEqn(-fvm::laplacian(_nu(), _T()) == _S());
        Eigen::SparseMatrix<double> A;
        Eigen::VectorXd b;
        Foam2Eigen::fvMatrix2Eigen(TEqn, A, b);
        return A;
    }

    Eigen::VectorXd get_rhs(Eigen::VectorXd& T, Eigen::VectorXd& S)
    {
        _T() = Foam2Eigen::Eigen2field(_T(), T);
        _S() = Foam2Eigen::Eigen2field(_S(), S);
        fvScalarMatrix TEqn(-fvm::laplacian(_nu(), _T()) == _S());
        Eigen::SparseMatrix<double> A;
        Eigen::VectorXd b;
        Foam2Eigen::fvMatrix2Eigen(TEqn, A, b);
        return b;
    }

    void solve()
    {
        fvScalarMatrix TEqn(-fvm::laplacian(_nu(), _T()) == _S());
        TEqn.solve();
    }

    void exportT(std::string& subFolder, std::string& folder, std::string& fieldname)
    {
        ITHACAstream::exportSolution(_T(), subFolder, folder, fieldname);
    }

    void setT(Eigen::VectorXd T)
    {
        _T() = Foam2Eigen::Eigen2field(_T(), T);
    }


    autoPtr<volScalarField> _res;
};


PYBIND11_MODULE(of_pybind11_system, m)
{
    // bindings to Matrix class
    py::class_<of_pybind11_system>(m, "of_pybind11_system")
        .def(py::init([](
                          std::vector<std::string> args) {
            std::vector<char*> cstrs;
            cstrs.reserve(args.size());
            for (auto& s : args)
                cstrs.push_back(const_cast<char*>(s.c_str()));
            return new of_pybind11_system(cstrs.size(), cstrs.data());
        }),
            py::arg("args") = std::vector<std::string> { "." })
        .def("getT", &of_pybind11_system::getT, py::return_value_policy::reference_internal)
        .def("setT", &of_pybind11_system::setT, py::return_value_policy::reference_internal)
        .def("printT", &of_pybind11_system::printT)
        .def("getS", &of_pybind11_system::getS, py::return_value_policy::reference_internal)
        .def("printS", &of_pybind11_system::printS)
        .def("getResidual", &of_pybind11_system::get_residual)
        .def("printMatrix", &of_pybind11_system::printMatrix)
        .def("get_rhs", &of_pybind11_system::get_rhs)
        .def("get_system_matrix", &of_pybind11_system::get_system_matrix)
        .def("solve", &of_pybind11_system::solve)
        .def("exportT", &of_pybind11_system::exportT);
}
