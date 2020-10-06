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
    newton_argument
Description
    A general newton argument is implemented
SourceFiles
    newton_argument.H
\*---------------------------------------------------------------------------*/

/// \file
/// Header file for the implementation of the newton_argument class, a simple
/// general object for the resolution of nonlinear problems with the Eigen library.


#include <Eigen/Eigen>

#ifndef newton_argument_H
#define newton_argument_H

/// Template object created to solve non linear problems using the Eigen library
template<typename _Scalar, label NX = Eigen::Dynamic, label NY = Eigen::Dynamic>
class newton_argument
{
    public:
        typedef _Scalar Scalar;
        enum
        {
            InputsAtCompileTime = NX,
            ValuesAtCompileTime = NY
        };
        typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
        typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
        typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime>
        JacobianType;

        label m_inputs, m_values;


        /// @brief      Constructor
        ///
        newton_argument() : m_inputs(InputsAtCompileTime),
            m_values(ValuesAtCompileTime) {}

        /// @brief      Constructor
        ///
        /// @param[in]  inputs  The inputs
        /// @param[in]  values  The values
        ///
        newton_argument(label inputs, label values) : m_inputs(inputs),
            m_values(values) {}

        /// @brief Operator to compute the residual vector
        ///
        /// @param[in]  x     The x value used to compute the residual
        /// @param      fvec  The residual vector
        ///
        /// @return     0
        ///
        label operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
        {
            std::cout <<
                      "This is a virtual method, it needs to be overriden with a suitable residual function"
                      << std::endl;
            return 0;
        }

        /// @brief Operator to compute the jacobian matrix
        ///
        /// @param[in]  x     The x value used to compute the jacobian
        /// @param      fvec  The jacobian vector
        ///
        /// @return     0
        ///
        label df(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
        {
            std::cout <<
                      "This is a virtual method, it needs to be overriden with a suitable jacobian function"
                      << std::endl;
            return 0;
        }

        /// @brief      Return the inputs
        ///
        /// @return     inputs
        ///
        label inputs() const
        {
            return m_inputs;
        }


        /// @brief      Return the values
        ///
        /// @return     inputs
        ///
        label values() const
        {
            return m_values;
        }
};

#endif
