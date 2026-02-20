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

#include "ITHACAfieldsOperations.H"

/// \file
/// Source file of the ITHACAfieldsOperations file.

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace ITHACAutilities
{

template<typename T>
void multField(T& f1, double alpha)
{
  int NBC = f1.boundaryField().size();
  Eigen::VectorXd f1v = Foam2Eigen::field2Eigen(f1);
  List<Eigen::VectorXd> f1BC = Foam2Eigen::field2EigenBC(f1);
  for (label k = 0; k < f1v.size(); k++)
  {
    f1v(k) *= alpha;
  }
  for (label l = 0; l < NBC; l++)
  {
    for (label k = 0; k < f1BC[l].size(); k++)
    {
      f1BC[l](k) *= alpha;
    }
  }

  f1 = Foam2Eigen::Eigen2field(f1, f1v);
  for (int k = 0; k < f1BC.size(); k++)
  {
    assignBC(f1, k, f1BC[k]);
  }
}
template void multField(volScalarField& f1, double alpha);
template void multField(volVectorField& f1, double alpha);
template void multField(volTensorField& f1, double alpha);

template<typename T>
void multField(T &f1, const Eigen::VectorXd alphaVec)
{
    Eigen::VectorXd f1v = Foam2Eigen::field2Eigen(f1);
    List<Eigen::VectorXd> f1BC = Foam2Eigen::field2EigenBC(f1);
    for (label k = 0; k < f1v.size(); k++)
    {
        f1v(k) *= (alphaVec[k]);
    }
    f1 = Foam2Eigen::Eigen2field(f1, f1v);
    for (int k = 0; k < f1BC.size(); k++)
    {
        assignBC(f1, k, f1BC[k]);
    }
}
template void multField(volScalarField& f1, const Eigen::VectorXd alphaVec);
template void multField(volVectorField& f1, const Eigen::VectorXd alphaVec);
template void multField(volTensorField& f1, const Eigen::VectorXd alphaVec);

template<typename T>
void multField(PtrList<T> &f1, const Eigen::VectorXd alphaVec)
{
  for(label ith_field = 0 ; ith_field < f1.size() ; ith_field++){
    Eigen::VectorXd f1v = Foam2Eigen::field2Eigen(f1[ith_field]);
    List<Eigen::VectorXd> f1BC = Foam2Eigen::field2EigenBC(f1[ith_field]);
    for (label k = 0; k < f1v.size(); k++)
    {
        f1v(k) *= (alphaVec[k]);
    }
    f1[ith_field] = Foam2Eigen::Eigen2field(f1[ith_field], f1v);
    for (int k = 0; k < f1BC.size(); k++)
    {
        assignBC(f1[ith_field], k, f1BC[k]);
    }
  }
}
template void multField(PtrList<volScalarField>& f1, const Eigen::VectorXd alphaVec);
template void multField(PtrList<volVectorField>& f1, const Eigen::VectorXd alphaVec);
template void multField(PtrList<volTensorField>& f1, const Eigen::VectorXd alphaVec);


template<typename T>
void addFields(T& f1, const T& f2c, double alpha)
{
  T f2 = f2c;
  int NBC = f1.boundaryField().size();
  Eigen::VectorXd f1v = Foam2Eigen::field2Eigen(f1);
  List<Eigen::VectorXd> f1BC = Foam2Eigen::field2EigenBC(f1);
  Eigen::VectorXd f2v = Foam2Eigen::field2Eigen(f2);
  List<Eigen::VectorXd> f2BC = Foam2Eigen::field2EigenBC(f2);

  for (label k = 0; k < f1v.size(); k++)
  {
    f1v(k) += alpha * f2v(k);
  }
  for (label l = 0; l < NBC; l++)
  {
    for (label k = 0; k < f1BC[l].size(); k++)
    {
      f1BC[l](k) += alpha * f2BC[l](k);
    }
  }

  f1 = Foam2Eigen::Eigen2field(f1, f1v);
  for (int k = 0; k < f1BC.size(); k++)
  {
    assignBC(f1, k, f1BC[k]);
  }
}
template void addFields(volScalarField& f1, const volScalarField& f2c, double alpha);
template void addFields(volVectorField& f1, const volVectorField& f2c, double alpha);
template void addFields(volTensorField& f1, const volTensorField& f2c, double alpha);


template<typename T>
void subtractFields(T& f1, const T& f2)
{
  addFields(f1, f2, -1.0);
}
template void subtractFields(volScalarField& f1, const volScalarField& f2);
template void subtractFields(volVectorField& f1, const volVectorField& f2);
template void subtractFields(volTensorField& f1, const volTensorField& f2);


volTensorField tensorFieldProduct(const volScalarField& coef, const volTensorField& S)
{
  return (coef * S);
}

volTensorField tensorFieldProduct(const volTensorField& coef, const volTensorField& S)
{
  return (coef & S);
}


int dimensionField(const volTensorField& v)
{
  int d = 9;
  return d;
}

int dimensionField(const volVectorField& v)
{
  int d = 3;
  return d;
}

int dimensionField(const volScalarField& v)
{
  int d = 1;
  return d;
}


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
void normalizeFields(
    PtrList<GeometricField<Type, fvPatchField, volMesh >>& fields)
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
    PtrList<GeometricField<scalar, fvPatchField, volMesh >>& fields);
template void normalizeFields(
    PtrList<GeometricField<vector, fvPatchField, volMesh >>& fields);


template<typename Type>
Eigen::MatrixXd getValues(GeometricField<Type, fvPatchField,
                          volMesh>& field, labelList& indices)
{
    List<Type> list(indices.size());
    M_Assert(max(indices) < field.size(),
             "The list indices are too large respect to field dimension");

    for (label i = 0; i < indices.size(); i++)
    {
        list[i] = field[indices[i]];
    }

    return Foam2Eigen::field2Eigen(list);
}

template<>
Eigen::MatrixXd getValues(GeometricField<vector, fvPatchField,
                          volMesh>& field, labelList& indices, labelList* xyz)
{
    M_Assert(max(indices) < field.size(),
             "The list of indices is too large respect to field dimension. There is at least one value larger than the dimension of the list");

    if (xyz != NULL)
    {
        List<scalar> list;
        list.resize(indices.size());
        M_Assert(max(* xyz) <= 2,
                 "The list of xyz positions contains at list one value larger than 2");
        labelList l = * xyz;

        for (label i = 0; i < indices.size(); i++)
        {
            list[i] = field[indices[i]][l[i]];
        }

        return Foam2Eigen::field2Eigen(list);
    }
    else
    {
        List<vector> list;
        list.resize(indices.size());

        for (label i = 0; i < indices.size(); i++)
        {
            list[i] = field[indices[i]];
        }

        return Foam2Eigen::field2Eigen(list);
    }
}

template<>
Eigen::MatrixXd getValues(GeometricField<scalar, fvPatchField,
                          volMesh>& field, labelList& indices, labelList* xyz)
{
    M_Assert(max(indices) < field.size(),
             "The list of indices is too large respect to field dimension. There is at least one value larger than the dimension of the list");
    List<scalar> list;
    list.resize(indices.size());

    for (label i = 0; i < indices.size(); i++)
    {
        list[i] = field[indices[i]];
    }

    return Foam2Eigen::field2Eigen(list);
}

template<typename T>
Eigen::MatrixXd getValues(PtrList<GeometricField<T, fvPatchField,
                          volMesh >>& fields, labelList& indices, labelList* xyz)
{
    Eigen::MatrixXd out;
    Eigen::MatrixXd a = getValues(fields[0], indices, xyz);
    out.resize(a.rows(), fields.size());
    out.col(0) = a;
    for (label i = 1; i < fields.size(); i++)
    {
        out.col(i) = getValues(fields[i], indices, xyz);
    }

    return out;
}


template
Eigen::MatrixXd getValues(PtrList<GeometricField<scalar, fvPatchField,
                                  volMesh >>& fields, labelList& indices, labelList* xyz);
template
Eigen::MatrixXd getValues(PtrList<GeometricField<vector, fvPatchField,
                                  volMesh >>& fields, labelList& indices, labelList* xyz);

}
