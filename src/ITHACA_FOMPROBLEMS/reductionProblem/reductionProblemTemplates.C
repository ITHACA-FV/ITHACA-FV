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
/// Template file of the reductionProblem class.

template<typename T, typename G>
void reductionProblem::assignIF(T& s, G& value)
{
    for (label i = 0; i < s.internalField().size(); i++)
    {
        s.ref()[i] = value;
    }
}

template<typename T>
void reductionProblem::computeLift(T& Lfield, T& liftfield, T& omfield)
{
    scalar u_bc;
    scalar u_lf;
    scalar area;

    for (label k = 0; k < inletIndex.rows(); k++)
    {
        label p = inletIndex(k, 0);
        label l = inletIndex(k, 1);
        area = gSum(Lfield[0].mesh().magSf().boundaryField()[p]);
        u_lf = gSum(liftfield[k].mesh().magSf().boundaryField()[p] *
                    liftfield[k].boundaryField()[p]).component(l) / area;
        M_Assert(std::abs(u_lf) > 1e-5,
                 "The lift cannot be computed. Please, check your inletIndex definition");

        for (label j = 0; j < Lfield.size(); j++)
        {
            if (k == 0)
            {
                u_bc = gSum(Lfield[j].mesh().magSf().boundaryField()[p] *
                            Lfield[j].boundaryField()[p]).component(l) / area;
                volVectorField C("U", Lfield[j] - liftfield[k]*u_bc / u_lf);
                omfield.append(C.clone());
            }
            else
            {
                u_bc = gSum(omfield[j].mesh().magSf().boundaryField()[p] *
                            omfield[j].boundaryField()[p]).component(l) / area;
                volVectorField C("U", omfield[j] - liftfield[k]*u_bc / u_lf);
                omfield.set(j, C.clone());
            }
        }
    }
}

template<typename T>
void reductionProblem::computeLiftT(T& Lfield, T& liftfield, T& omfield)
{
    scalar t_bc;
    scalar t_lf;
    scalar area;

    for (label k = 0; k < inletIndexT.rows(); k++)
    {
        label p = inletIndexT(k, 0);
        area = gSum(Lfield[0].mesh().magSf().boundaryField()[p]);
        t_lf = gSum(liftfield[k].mesh().magSf().boundaryField()[p] *
                    liftfield[k].boundaryField()[p]) / area;

        for (label j = 0; j < Lfield.size(); j++)
        {
            if (k == 0)
            {
                t_bc = gSum(Lfield[j].mesh().magSf().boundaryField()[p] *
                            Lfield[j].boundaryField()[p]) / area;
                volScalarField C(Lfield[0].name(), Lfield[j] - liftfield[k]*t_bc / t_lf);
                omfield.append(C.clone());
            }
            else
            {
                t_bc = gSum(omfield[j].mesh().magSf().boundaryField()[p] *
                            omfield[j].boundaryField()[p]) / area;
                volScalarField C(Lfield[0].name(), omfield[j] - liftfield[k]*t_bc / t_lf);
                omfield.set(j, C.clone());
            }
        }
    }
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
