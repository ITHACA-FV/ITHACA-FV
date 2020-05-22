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
#include "ITHACAassign.H"

namespace ITHACAutilities
{

template<class TypeField>
PtrList<TypeField> averageSubtract(PtrList<TypeField>
                                   fields, Eigen::MatrixXd ind, PtrList<TypeField>& ave)
{
    PtrList<TypeField> aveSubtracted;
    Eigen::VectorXd newInd;
    newInd.resize(ind.size() + 1);
    newInd.head(ind.size()) = ind;
    newInd(ind.size()) = fields.size();

    for (label i = 0; i < ind.size(); i++)
    {
        TypeField aveTemp("nut", fields[0] * 0);

        for (label j = newInd(i); j < newInd(i + 1); j++)
        {
            aveTemp += fields[j];
        }

        aveTemp /= newInd(i + 1) - newInd(i);
        ave.append(aveTemp);
    }

    for (label i = 0; i < ind.size(); i++)
    {
        for (label j = newInd(i); j < newInd(i + 1); j++)
        {
            TypeField newfield("nut", fields[0] * 0);
            newfield = fields[j] - ave[i];
            aveSubtracted.append(newfield);
        }
    }

    return aveSubtracted;
}

template PtrList<volScalarField> averageSubtract(
    PtrList<volScalarField>
    fields, Eigen::MatrixXd ind, PtrList<volScalarField>& ave);
template PtrList<volVectorField> averageSubtract(
    PtrList<volVectorField>
    fields, Eigen::MatrixXd ind, PtrList<volVectorField>& ave);

template<class TypeField>
TypeField computeAverage(PtrList<TypeField>& fields)
{
    TypeField av(fields[0]);

    for (int i = 1; i < fields.size(); i++)
    {
        av += fields[i];
    }

    av = av / fields.size();
    return av;
}

template volVectorField computeAverage(
    PtrList<volVectorField>& fields);
template volScalarField computeAverage(
    PtrList<volScalarField>& fields);


template<typename T>
void assignIF(GeometricField<T, fvPatchField, volMesh>& s,
              T& value)
{
    for (label i = 0; i < s.internalField().size(); i++)
    {
        s.ref()[i] = value;
    }
}

template void assignIF(
    GeometricField<scalar, fvPatchField, volMesh>& field, scalar& value);
template void assignIF(
    GeometricField<vector, fvPatchField, volMesh>& field, vector& value);

template<typename T>
void assignIF(GeometricField<T, fvPatchField, volMesh>& s,
              T& value, List<int>& indices)
{
    for (label i = 0; i < indices.size(); i++)
    {
        s.ref()[indices[i]] = value;
    }
}

template void assignIF(GeometricField<scalar, fvPatchField, volMesh>& s,
                       scalar& value, List<int>& indices);
template void assignIF(GeometricField<vector, fvPatchField, volMesh>& s,
                       vector& value, List<int>& indices);

template<typename T>
void assignIF(GeometricField<T, fvPatchField, volMesh>& s, T& value, int index)
{
    s.ref()[index] = value;
}

template void assignIF(GeometricField<scalar, fvPatchField, volMesh>& field,
                       scalar& value, int index);
template void assignIF(GeometricField<vector, fvPatchField, volMesh>& field,
                       vector& value, int index);

void assignONE(volScalarField& s, List<int>& L)
{
    for (label i = 0; i < L.size(); i++)
    {
        s.ref()[L[i]] = 1;
    }
}

// Assign a BC for a vector field
void assignBC(volScalarField& s, label BC_ind, double& value)
{
    if (s.boundaryField()[BC_ind].type() == "fixedValue"
            || s.boundaryField()[BC_ind].type() == "calculated")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "fixedGradient")
    {
        fixedGradientFvPatchScalarField& Tpatch =
            refCast<fixedGradientFvPatchScalarField>(s.boundaryFieldRef()[BC_ind]);
        scalarField& gradTpatch = Tpatch.gradient();
        forAll(gradTpatch, faceI)
        {
            gradTpatch[faceI] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "fixedFluxPressure")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "freestream")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }

        freestreamFvPatchField<scalar>& Tpatch =
            refCast<freestreamFvPatchField<scalar>>(s.boundaryFieldRef()[BC_ind]);
        scalarField& gradTpatch = Tpatch.freestreamValue();
        forAll(gradTpatch, faceI)
        {
            gradTpatch[faceI] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "calculated")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "processor")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "empty")
    {
    }
}

// Assign a BC for a scalar field
void assignBC(volScalarField& s, label BC_ind,
              List<double> value)
{
    if (s.boundaryField()[BC_ind].type() == "fixedValue")
    {
        s.boundaryFieldRef()[BC_ind] = value;
    }
    else if (s.boundaryField()[BC_ind].type() == "fixedGradient")
    {
        fixedGradientFvPatchScalarField& Tpatch =
            refCast<fixedGradientFvPatchScalarField>(s.boundaryFieldRef()[BC_ind]);
        scalarField& gradTpatch = Tpatch.gradient();
        gradTpatch = value;
    }
    else if (s.boundaryField()[BC_ind].type() == "empty")
    {
    }
    else
    {
        Info << "This type of boundary condition is not yet implemented, code will abort"
             << endl;
        exit(0);
    }
}

// Assign a BC for a scalar field
void assignBC(volVectorField& s, label BC_ind,
              Vector<double>& value)
{
    if (s.boundaryField()[BC_ind].type() == "fixedValue"
            || s.boundaryField()[BC_ind].type() == "processor"
            || s.boundaryField()[BC_ind].type() == "calculated")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "fixedGradient")
    {
        Info << "This Feature is not implemented for this boundary condition" << endl;
        exit(0);
    }
    else if (s.boundaryField()[BC_ind].type() == "freestream")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }

        freestreamFvPatchField<vector>& Tpatch =
            refCast<freestreamFvPatchField<vector>>(s.boundaryFieldRef()[BC_ind]);
        vectorField& gradTpatch = Tpatch.freestreamValue();
        forAll(gradTpatch, faceI)
        {
            gradTpatch[faceI] = value;
        }
    }
}

void assignBC(volScalarField& s, label BC_ind,
              Eigen::MatrixXd valueVec)
{
    word typeBC = s.boundaryField()[BC_ind].type();

    if (typeBC == "fixedValue" || typeBC == "calculated"
            || typeBC == "fixedFluxPressure" ||  typeBC == "processor")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            double value = valueVec(i);
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "fixedGradient")
    {
        fixedGradientFvPatchScalarField& Tpatch =
            refCast<fixedGradientFvPatchScalarField>(s.boundaryFieldRef()[BC_ind]);
        scalarField& gradTpatch = Tpatch.gradient();
        forAll(gradTpatch, faceI)
        {
            double value = valueVec(faceI);
            gradTpatch[faceI] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "freestream")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            double value = valueVec(i);
            s.boundaryFieldRef()[BC_ind][i] = value;
        }

        freestreamFvPatchField<scalar>& Tpatch =
            refCast<freestreamFvPatchField<scalar>>(s.boundaryFieldRef()[BC_ind]);
        scalarField& gradTpatch = Tpatch.freestreamValue();
        forAll(gradTpatch, faceI)
        {
            double value = valueVec(faceI);
            gradTpatch[faceI] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "empty")
    {
    }
}

// Assign a BC for a scalar field
void assignBC(volVectorField& s, label BC_ind,
              Eigen::MatrixXd valueVec)
{
    word typeBC = s.boundaryField()[BC_ind].type();
    int sizeBC = s.boundaryField()[BC_ind].size();

    if (typeBC == "fixedValue" || typeBC == "calculated" || typeBC == "processor")
    {
        for (label i = 0; i < sizeBC; i++)
        {
            Vector<double> value(valueVec(i), valueVec(i + sizeBC),
                                 valueVec(i + sizeBC * 2));
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "fixedGradient")
    {
        Info << "This Feature is not implemented for this boundary condition" << endl;
        exit(0);
    }
    else if (s.boundaryField()[BC_ind].type() == "freestream")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            Vector<double> value(valueVec(i), valueVec(i + sizeBC),
                                 valueVec(i + sizeBC * 2));
            s.boundaryFieldRef()[BC_ind][i] = value;
        }

        freestreamFvPatchField<vector>& Tpatch =
            refCast<freestreamFvPatchField<vector>>(s.boundaryFieldRef()[BC_ind]);
        vectorField& gradTpatch = Tpatch.freestreamValue();
        forAll(gradTpatch, faceI)
        {
            Vector<double> value(valueVec(faceI), valueVec(faceI + sizeBC),
                                 valueVec(faceI + sizeBC * 2));
            gradTpatch[faceI] = value;
        }
    }
}


template<typename T>
void setBoxToValue(GeometricField<T, fvPatchField, volMesh>& field,
                   Eigen::MatrixXd Box, T value)
{
    M_Assert(Box.rows() == 2
             && Box.cols() == 3,
             "The box must be a 2*3 matrix shaped in this way: \nBox = \t|x0, y0, z0|\n\t|x1, yi, z1|\n");

    for (label i = 0; i < field.internalField().size(); i++)
    {
        auto cx = field.mesh().C()[i].component(vector::X);
        auto cy = field.mesh().C()[i].component(vector::Y);
        auto cz = field.mesh().C()[i].component(vector::Z);

        if (cx >= Box(0, 0) && cy >= Box(0, 1) && cz >= Box(0, 2) && cx <= Box(1, 0)
                && cy <= Box(1, 1) && cz <= Box(1, 2) )
        {
            field.ref()[i] = value;
        }
    }

    for (label i = 0; i < field.boundaryField().size(); i++)
    {
        for (label j = 0; j < field.boundaryField()[i].size(); j++)
        {
            if (field.boundaryField()[i].type() == "fixedValue"
                    || field.boundaryField()[i].type() == "calculated")
            {
                auto cx = field.mesh().C().boundaryField()[i][j][0];
                auto cy = field.mesh().C().boundaryField()[i][j][1];
                auto cz = field.mesh().C().boundaryField()[i][j][2];

                if (cx >= Box(0, 0) && cy >= Box(0, 1) && cz >= Box(0, 2) && cx <= Box(1, 0)
                        && cy <= Box(1, 1) && cz <= Box(1, 2) )
                {
                    field.boundaryFieldRef()[i][j] = value;
                }
            }
        }
    }
}

template void setBoxToValue(GeometricField<scalar, fvPatchField, volMesh>&
                            field, Eigen::MatrixXd Box, scalar value);
template void setBoxToValue(GeometricField<vector, fvPatchField, volMesh>&
                            field, Eigen::MatrixXd Box, vector value);

template<typename T>
void setIndices2Value(labelList& ind2set, List<T>& value2set,
                      labelList& movingIDS, List<T>& originalList)
{
    M_Assert(ind2set.size() == value2set.size(),
             "The size of the indices must be equal to the size of the values list");
    M_Assert(originalList.size() >= value2set.size(),
             "The size of the original list of values must be bigger than the size of the list of values you want to set");
    labelList ind_ok(ind2set);

    for (int i = 0; i < ind2set.size(); i++)
    {
        for (int k = 0; k < movingIDS.size(); k++)
        {
            if (ind2set[i] == movingIDS[k])
            {
                ind_ok[i] = k;
                break;
            }
        }
    }

    for (int i = 0; i < ind2set.size(); i++)
    {
        originalList[ind_ok[i]] = value2set[i];
    }
}

template void setIndices2Value(labelList& ind2set, List<scalar>& value2set,
                               labelList& movingIDS, List<scalar>& originalList);
template void setIndices2Value(labelList& ind2set, List<vector>& value2set,
                               labelList& movingIDS, List<vector>& originalList);

template<class TypeField>
void changeBCtype(
    GeometricField<TypeField, fvPatchField, volMesh>& field, word BCtype,
    label BC_ind)
{
    field.boundaryFieldRef().set(BC_ind, fvPatchField<TypeField>::New(BCtype,
                                 field.mesh().boundary()[BC_ind], field));
}

template void changeBCtype<scalar>
(GeometricField<scalar, fvPatchField, volMesh>& field, word BCtype,
 label BC_ind);
template void changeBCtype<vector>
(GeometricField<vector, fvPatchField, volMesh>& field, word BCtype,
 label BC_ind);

void assignMixedBC(
    GeometricField<scalar, fvPatchField, volMesh>& field, label BC_ind,
    Eigen::MatrixXd& value, Eigen::MatrixXd& grad, Eigen::MatrixXd& valueFrac)
{
    std::string message = "Patch is NOT mixed. It is of type: " +
                          field.boundaryField()[BC_ind].type();
    M_Assert(field.boundaryField()[BC_ind].type() == "mixed", message.c_str());

    if (field.boundaryField()[BC_ind].type() == "mixed")
    {
        mixedFvPatchScalarField& Tpatch =
            refCast<mixedFvPatchScalarField>(field.boundaryFieldRef()[BC_ind]);
        scalarField& valueTpatch = Tpatch.refValue();
        scalarField& gradTpatch = Tpatch.refGrad();
        scalarField& valueFracTpatch = Tpatch.valueFraction();
        Foam2Eigen::Eigen2field(valueTpatch, value);
        Foam2Eigen::Eigen2field(gradTpatch, grad);
        Foam2Eigen::Eigen2field(valueFracTpatch, valueFrac);
    }
}

void assignMixedBC(
    GeometricField<vector, fvPatchField, volMesh>& field, label BC_ind,
    Eigen::MatrixXd& value, Eigen::MatrixXd& grad, Eigen::MatrixXd& valueFrac)
{
    std::string message = "Patch is NOT mixed. It is of type: " +
                          field.boundaryField()[BC_ind].type();
    M_Assert(field.boundaryField()[BC_ind].type() == "mixed", message.c_str());

    if (field.boundaryField()[BC_ind].type() == "mixed")
    {
        mixedFvPatchVectorField& Tpatch =
            refCast<mixedFvPatchVectorField>(field.boundaryFieldRef()[BC_ind]);
        vectorField& valueTpatch = Tpatch.refValue();
        vectorField& gradTpatch = Tpatch.refGrad();
        scalarField& valueFracTpatch = Tpatch.valueFraction();
        Foam2Eigen::Eigen2field(valueTpatch, value);
        Foam2Eigen::Eigen2field(gradTpatch, grad);
        Foam2Eigen::Eigen2field(valueFracTpatch, valueFrac);
    }
}

void assignMixedBC(
    GeometricField<scalar, fvPatchField, volMesh>& field, label BC_ind,
    List<scalar>& value, List<scalar>& grad, List<scalar>& valueFrac)
{
    std::string message = "Patch is NOT mixed. It is of type: " +
                          field.boundaryField()[BC_ind].type();
    M_Assert(field.boundaryField()[BC_ind].type() == "mixed", message.c_str());

    if (field.boundaryField()[BC_ind].type() == "mixed")
    {
        mixedFvPatchScalarField& Tpatch =
            refCast<mixedFvPatchScalarField>(field.boundaryFieldRef()[BC_ind]);
        scalarField& valueTpatch = Tpatch.refValue();
        scalarField& gradTpatch = Tpatch.refGrad();
        scalarField& valueFracTpatch = Tpatch.valueFraction();
        valueTpatch = value;
        gradTpatch = grad;
        valueFracTpatch = valueFrac;
    }
}

void assignMixedBC(
    GeometricField<vector, fvPatchField, volMesh>& field, label BC_ind,
    List<vector>& value, List<vector>& grad, List<scalar>& valueFrac)
{
    std::string message = "Patch is NOT mixed. It is of type: " +
                          field.boundaryField()[BC_ind].type();
    M_Assert(field.boundaryField()[BC_ind].type() == "mixed", message.c_str());

    if (field.boundaryField()[BC_ind].type() == "mixed")
    {
        mixedFvPatchVectorField& Tpatch =
            refCast<mixedFvPatchVectorField>(field.boundaryFieldRef()[BC_ind]);
        vectorField& valueTpatch = Tpatch.refValue();
        vectorField& gradTpatch = Tpatch.refGrad();
        scalarField& valueFracTpatch = Tpatch.valueFraction();
        valueTpatch = value;
        gradTpatch = grad;
        valueFracTpatch = valueFrac;
    }
}

}
