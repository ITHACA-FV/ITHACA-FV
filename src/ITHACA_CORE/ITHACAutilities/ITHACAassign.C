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
        ave.append(aveTemp.clone());
    }

    for (label i = 0; i < ind.size(); i++)
    {
        for (label j = newInd(i); j < newInd(i + 1); j++)
        {
            TypeField newfield("nut", fields[0] * 0);
            newfield = fields[j] - ave[i];
            aveSubtracted.append(newfield.clone());
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

    for (label i = 1; i < fields.size(); i++)
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


template<typename Type>
void assignIF(GeometricField<Type, fvPatchField, volMesh>& s,
              Type value)
{
    for (label i = 0; i < s.internalField().size(); i++)
    {
        s.ref()[i] = value;
    }
}

template void assignIF(
    GeometricField<scalar, fvPatchField, volMesh>& field, scalar value);
template void assignIF(
    GeometricField<vector, fvPatchField, volMesh>& field, vector value);

template<typename Type>
void assignIF(GeometricField<Type, fvPatchField, volMesh>& s,
              Type& value, List<label>& indices)
{
    for (label i = 0; i < indices.size(); i++)
    {
        s.ref()[indices[i]] = value;
    }
}

template void assignIF(GeometricField<scalar, fvPatchField, volMesh>& s,
                       scalar& value, List<label>& indices);
template void assignIF(GeometricField<vector, fvPatchField, volMesh>& s,
                       vector& value, List<label>& indices);

template<typename Type>
void assignIF(GeometricField<Type, fvPatchField, volMesh>& s,
              Type& value, label index)
{
    s.ref()[index] = value;
}

template void assignIF(GeometricField<scalar, fvPatchField, volMesh>& field,
                       scalar& value, label index);
template void assignIF(GeometricField<vector, fvPatchField, volMesh>& field,
                       vector& value, label index);

void assignONE(volScalarField& s, List<label>& L)
{
    for (label i = 0; i < L.size(); i++)
    {
        s.ref()[L[i]] = 1;
    }
}

void assignBC(GeometricField<scalar, fvPatchField, volMesh>& s, label BC_ind,
              double value)
{
    label sizeBC = s.boundaryField()[BC_ind].size();
    List<double> valueList(sizeBC);

    for (label i = 0; i < sizeBC; i++)
    {
        valueList[i] = value;
    }

    assignBC(s, BC_ind, valueList);
}

void assignBC(GeometricField<scalar, fvPatchField, volMesh>& s, label BC_ind,
              Eigen::MatrixXd valueVec)
{
    label sizeBC = s.boundaryField()[BC_ind].size();
    M_Assert(sizeBC == valueVec.size(),
             "The size of the given values matrix has to be equal to the dimension of the boundaryField");
    List<double> valueList(sizeBC);

    for (label i = 0; i < sizeBC; i++)
    {
        valueList[i] = valueVec(i);
    }

    assignBC(s, BC_ind, valueList);
}

// Assign a BC for a scalar field
void assignBC(GeometricField<scalar, fvPatchField, volMesh>& s, label BC_ind,
              List<double> valueList)
{
    word typeBC = s.boundaryField()[BC_ind].type();
    label sizeBC = s.boundaryField()[BC_ind].size();
    M_Assert(sizeBC == valueList.size(),
             "The size of the given values list has to be equal to the dimension of the boundaryField");

    if (typeBC == "fixedGradient")
    {
        fixedGradientFvPatchScalarField& Tpatch =
            refCast<fixedGradientFvPatchScalarField>(s.boundaryFieldRef()[BC_ind]);
        scalarField& gradTpatch = Tpatch.gradient();
        forAll(gradTpatch, faceI)
        {
            double value = valueList[faceI];
            gradTpatch[faceI] = value;
        }
    }
    else if (typeBC == "freestream")
    {
        for (label i = 0; i < sizeBC; i++)
        {
            double value = valueList[i];
            s.boundaryFieldRef()[BC_ind][i] = value;
        }

        freestreamFvPatchField<scalar>& Tpatch =
            refCast<freestreamFvPatchField<scalar>>(s.boundaryFieldRef()[BC_ind]);
        scalarField& gradTpatch = Tpatch.freestreamValue();
        forAll(gradTpatch, faceI)
        {
            double value = valueList[faceI];
            gradTpatch[faceI] = value;
        }
    }
    else if (typeBC == "empty" || typeBC == "zeroGradient")
    {}
    else
    {
        try
        {
            if (typeBC != "fixedGradient" && typeBC != "freestream" && typeBC != "empty"
                    && typeBC != "zeroGradient" && typeBC != "fixedValue" && typeBC != "calculated"
                    && typeBC != "fixedFluxPressure" && typeBC != "processor"
                    && typeBC != "nutkWallFunction" && typeBC != "mixedEnergy")
            {
                word message = "Pay attention, your typeBC " + typeBC + " for " + s.name() +
                               " is not included into the developed ones. Your BC will be treated as a classical fixedValue.";
                throw (message);
            }
        }
        catch (const word message)
        {
            std::cerr << "WARNING: " << message << endl;
        }

        for (label i = 0; i < sizeBC; i++)
        {
            double value = valueList[i];
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
}

void assignBC(GeometricField<vector, fvPatchField, volMesh>& s, label BC_ind,
              vector value)
{
    M_Assert(value.size() == 3,
             "The size of the given vector has to be equal to 3 for the 3 components");
    label sizeBC = s.boundaryField()[BC_ind].size();
    List<vector> valueList(sizeBC);

    for (label i = 0; i < sizeBC; i++)
    {
        valueList[i] = value;
    }

    assignBC(s, BC_ind, valueList);
}

void assignBC(GeometricField<vector, fvPatchField, volMesh>& s, label BC_ind,
              Eigen::MatrixXd valueVec)
{
    label sizeBC = s.boundaryField()[BC_ind].size();
    M_Assert(sizeBC * 3 == valueVec.size(),
             "The size of the given values matrix has to be equal to 3 times the dimension of the boundaryField");
    List<vector> valueList(sizeBC);

    for (label i = 0; i < sizeBC; i++)
    {
        valueList[i].component(0) = valueVec(i);
        valueList[i].component(1) = valueVec(i + sizeBC);
        valueList[i].component(2) = valueVec(i + sizeBC * 2);
    }

    assignBC(s, BC_ind, valueList);
}

void assignBC(GeometricField<scalar, fvsPatchField, surfaceMesh>& s,
              label BC_ind, Eigen::MatrixXd valueVec)
{
    label sizeBC = s.boundaryField()[BC_ind].size();
    M_Assert(sizeBC  == valueVec.rows() && valueVec.cols() == 1,
             "The given matrix must be a column one with the size equal to 3 times the dimension of the boundaryField");
    List<scalar> valueList(sizeBC);

    for (label i = 0; i < sizeBC; i++)
    {
        valueList[i] = valueVec(i);
    }

    assignBC(s, BC_ind, valueList);
}

void assignBC(GeometricField<vector, fvsPatchField, surfaceMesh>& s,
              label BC_ind, Eigen::MatrixXd valueVec)
{
    label sizeBC = s.boundaryField()[BC_ind].size();
    M_Assert(sizeBC * 3  == valueVec.rows() && valueVec.cols() == 1,
             "The given matrix must be a column one with the size equal to the dimension of the boundaryField");
    List<vector> valueList(sizeBC);

    for (label i = 0; i < sizeBC; i++)
    {
        valueList[i].component(0) = valueVec(i);
        valueList[i].component(1) = valueVec(i + sizeBC);
        valueList[i].component(2) = valueVec(i + sizeBC * 2);
    }

    assignBC(s, BC_ind, valueList);
}

// Assign a BC for a vector field
void assignBC(GeometricField<vector, fvPatchField, volMesh>& s, label BC_ind,
              List<vector> valueList)
{
    word typeBC = s.boundaryField()[BC_ind].type();
    label sizeBC = s.boundaryField()[BC_ind].size();
    M_Assert(sizeBC == valueList.size(),
             "The size of the given values list has to be equal to the dimension of the boundaryField");

    if (s.boundaryField()[BC_ind].type() == "fixedGradient")
    {
        Info << "This Feature is not implemented for this boundary condition" << endl;
        exit(0);
    }
    else if (s.boundaryField()[BC_ind].type() == "freestream")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = valueList[i];
        }

        freestreamFvPatchField<vector>& Tpatch =
            refCast<freestreamFvPatchField<vector>>(s.boundaryFieldRef()[BC_ind]);
        vectorField& gradTpatch = Tpatch.freestreamValue();
        forAll(gradTpatch, faceI)
        {
            gradTpatch[faceI] = valueList[faceI];
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "empty"
             || s.boundaryField()[BC_ind].type() == "zeroGradient")
    {}
    else
    {
        try
        {
            if (typeBC != "fixedGradient" && typeBC != "freestream" && typeBC != "empty"
                    && typeBC != "zeroGradient" && typeBC != "fixedValue" && typeBC != "calculated"
                    &&  typeBC != "processor")
            {
                word message = "Pay attention, your typeBC " + typeBC + " for " + s.name() +
                               " is not included into the developed ones. Your BC will be treated as a classical fixedValue.";
                throw (message);
            }
        }
        catch (const word message)
        {
            cerr << "WARNING: " << message << endl;
        }

        for (label i = 0; i < sizeBC; i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = valueList[i];
        }
    }
}

template<typename Type>
void assignBC(GeometricField<Type, fvsPatchField, surfaceMesh>& s, label BC_ind,
              List<Type>& valueList)
{
    word typeBC = s.boundaryField()[BC_ind].type();
    label sizeBC = s.boundaryField()[BC_ind].size();
    M_Assert(sizeBC == valueList.size(),
             "The size of the given values list has to be equal to the dimension of the boundaryField");

    if (s.boundaryField()[BC_ind].type() == "fixedGradient")
    {
        fixedGradientFvPatchField<Type>& Tpatch =
            refCast<fixedGradientFvPatchField<Type>>(s.boundaryFieldRef()[BC_ind]);
        Field<Type>& gradTpatch = Tpatch.gradient();
        forAll(gradTpatch, faceI)
        {
            gradTpatch[faceI] = valueList[faceI];
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "freestream")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = valueList[i];
        }

        freestreamFvPatchField<Type>& Tpatch =
            refCast<freestreamFvPatchField<Type>>(s.boundaryFieldRef()[BC_ind]);
        Field<Type>& gradTpatch = Tpatch.freestreamValue();
        forAll(gradTpatch, faceI)
        {
            gradTpatch[faceI] = valueList[faceI];
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "empty"
             || s.boundaryField()[BC_ind].type() == "zeroGradient")
    {}
    else
    {
        try
        {
            if (typeBC != "fixedGradient" && typeBC != "freestream" && typeBC != "empty"
                    || typeBC != "zeroGradient" && typeBC != "fixedValue" && typeBC != "calculated"
                    && typeBC != "fixedFluxPressure" &&  typeBC != "processor")
            {
                word message = "Pay attention, your typeBC " + typeBC + " for " + s.name() +
                               " is not included into the developed ones. Your BC will be treated as a classical fixedValue.";
                throw (message);
            }
        }
        catch (const word message)
        {
            cerr << "WARNING: " << message << endl;
        }

        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = valueList[i];
        }
    }
}

template void assignBC(
    GeometricField<scalar, fvsPatchField, surfaceMesh>& s, label BC_ind,
    List<scalar>& valueList);
template void assignBC(
    GeometricField<vector, fvsPatchField, surfaceMesh>& s, label BC_ind,
    List<vector>& valueList);

template<typename Type>
void assignBC(GeometricField<Type, fvsPatchField, surfaceMesh>& s, label BC_ind,
              Type& value)
{
    label sizeBC = s.boundaryField()[BC_ind].size();
    List<Type> valueList(sizeBC);

    for (label i = 0; i < sizeBC; i++)
    {
        valueList[i] = value;
    }

    assignBC(s, BC_ind, valueList);
}

template void assignBC(
    GeometricField<scalar, fvsPatchField, surfaceMesh>& s, label BC_ind,
    scalar& valueList);
template void assignBC(
    GeometricField<vector, fvsPatchField, surfaceMesh>& s, label BC_ind,
    vector& valueList);

template<typename Type>
void changeNeumann2Dirichlet(GeometricField<Type, fvPatchField, volMesh>& field,
                             Type& value)
{
    forAll(field.mesh().boundary(), i)
    {
        if (field.boundaryField()[i].type() == "zeroGradient" ||
                field.boundaryField()[i].type() ==  "fixedGradient")
        {
            ITHACAutilities::changeBCtype(field, "fixedValue", i);
            assignBC(field, i, value);
        }
    }
}

template void changeNeumann2Dirichlet(
    GeometricField<scalar, fvPatchField, volMesh>& field, scalar& value);
template void changeNeumann2Dirichlet(
    GeometricField<vector, fvPatchField, volMesh>& field, vector& value);

template<typename Type>
void assignZeroDirichlet(GeometricField<Type, fvPatchField, volMesh>& field)
{
    Type v;
    v = v * 0;
    assignIF(field, v);
    forAll(field.mesh().boundary(), i)
    {
        if (field.boundaryField()[i].type() == "fixedValue")
        {
            assignBC(field, i, v);
        }
    }
    changeNeumann2Dirichlet(field, v);
}
template void assignZeroDirichlet(
    GeometricField<vector, fvPatchField, volMesh>& field);
template void assignZeroDirichlet(
    GeometricField<scalar, fvPatchField, volMesh>& field);

template<typename Type>
void setBoxToValue(GeometricField<Type, fvPatchField, volMesh>& field,
                   Eigen::MatrixXd Box, Type value)
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

template<typename Type>
void setIndices2Value(labelList& ind2set, List<Type>& value2set,
                      labelList& movingIDS, List<Type>& originalList)
{
    M_Assert(ind2set.size() == value2set.size(),
             "The size of the indices must be equal to the size of the values list");
    M_Assert(originalList.size() >= value2set.size(),
             "The size of the original list of values must be bigger than the size of the list of values you want to set");
    labelList ind_ok(ind2set);

    for (label i = 0; i < ind2set.size(); i++)
    {
        for (label k = 0; k < movingIDS.size(); k++)
        {
            if (ind2set[i] == movingIDS[k])
            {
                ind_ok[i] = k;
                break;
            }
        }
    }

    for (label i = 0; i < ind2set.size(); i++)
    {
        originalList[ind_ok[i]] = value2set[i];
    }
}

template void setIndices2Value(labelList& ind2set, List<scalar>& value2set,
                               labelList& movingIDS, List<scalar>& originalList);
template void setIndices2Value(labelList& ind2set, List<vector>& value2set,
                               labelList& movingIDS, List<vector>& originalList);

template<class Type>
void changeBCtype(
    GeometricField<Type, fvPatchField, volMesh>& field, word BCtype,
    label BC_ind)
{
    field.boundaryFieldRef().set(BC_ind, fvPatchField<Type>::New(BCtype,
                                 field.mesh().boundary()[BC_ind], field));
}

template void changeBCtype<scalar>
(GeometricField<scalar, fvPatchField, volMesh>& field, word BCtype,
 label BC_ind);
template void changeBCtype<vector>
(GeometricField<vector, fvPatchField, volMesh>& field, word BCtype,
 label BC_ind);

template<typename Type>
void assignMixedBC(
    GeometricField<Type, fvPatchField, volMesh>& field, label BC_ind,
    List<Type>& value, List<Type>& grad, List<scalar>& valueFrac)
{
    std::string message = "Patch is NOT mixed. It is of type: " +
                          field.boundaryField()[BC_ind].type();
    M_Assert(field.boundaryField()[BC_ind].type() == "mixed", message.c_str());

    if (field.boundaryField()[BC_ind].type() == "mixed")
    {
        mixedFvPatchField<Type>& Tpatch =
            refCast<mixedFvPatchField<Type>>(field.boundaryFieldRef()[BC_ind]);
        Field<Type>& valueTpatch = Tpatch.refValue();
        Field<Type>& gradTpatch = Tpatch.refGrad();
        Field<scalar>& valueFracTpatch = Tpatch.valueFraction();
        valueTpatch = value;
        gradTpatch = grad;
        valueFracTpatch = valueFrac;
    }
}

template void assignMixedBC<scalar>(
    GeometricField<scalar, fvPatchField, volMesh>& field, label BC_ind,
    List<scalar>& value, List<scalar>& grad, List<scalar>& valueFrac);

template void assignMixedBC<vector>(
    GeometricField<vector, fvPatchField, volMesh>& field, label BC_ind,
    List<vector>& value, List<vector>& grad, List<scalar>& valueFrac);

template<typename Type>
void normalizeFields(
    PtrList<GeometricField<Type, fvPatchField, volMesh>>& fields)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());
    word normType = para->ITHACAdict->lookupOrDefault<word>("normalizationNorm",
                    "L2");
    M_Assert(normType == "L2" ||
             normType == "Frobenius", "The normalizationNorm can be only L2 or Frobenius" );
    Eigen::MatrixXd eigenFields = Foam2Eigen::PtrList2Eigen(fields);
    List<Eigen::MatrixXd> eigenFieldsBC = Foam2Eigen::PtrList2EigenBC(fields);

    for (label i = 0; i < fields.size(); i++)
    {
        double norm;

        if (normType == "L2")
        {
            norm = L2Norm(fields[i]);
        }
        else if (normType == "Frobenius")
        {
            norm = frobNorm(fields[i]);
        }

        GeometricField<Type, fvPatchField, volMesh> tmp2(fields[0].name(),
                fields[0] * 0);
        Eigen::VectorXd vec = eigenFields.col(i) / norm;
        tmp2 = Foam2Eigen::Eigen2field(tmp2, vec);

        // Adjusting boundary conditions
        for (label k = 0; k < tmp2.boundaryField().size(); k++)
        {
            Eigen::MatrixXd vec = eigenFieldsBC[k].col(i) / norm;
            assignBC(tmp2, k, vec);
        }

        fields.set(i, tmp2.clone());
    }
}

template void normalizeFields(
    PtrList<GeometricField<scalar, fvPatchField, volMesh>>& fields);
template void normalizeFields(
    PtrList<GeometricField<vector, fvPatchField, volMesh>>& fields);

}
