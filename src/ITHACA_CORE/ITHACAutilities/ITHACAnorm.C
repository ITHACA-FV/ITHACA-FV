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

#include "ITHACAnorm.H"

/// \file
/// Source file of the ITHACAnorm file.

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace ITHACAutilities
{

template<typename T>
Eigen::MatrixXd dot_product_POD(PtrList<T>& v, PtrList<T>& w,
                                const word& hilbertSpacePOD,
                                const double& weightBC,
                                const word& patchBC,
                                const double& weightH1)
{
  Eigen::MatrixXd matrix_dot_product = Eigen::MatrixXd::Zero(v.size(),w.size());
  if (hilbertSpacePOD == "L2")
  {
    matrix_dot_product = getMassMatrix(v,w);
  }
  else if (hilbertSpacePOD == "dL2"){
    Eigen::VectorXd V = getMassMatrixFV(v[0]).array().pow(4.0/3.0);
    V = V.array().inverse();
    matrix_dot_product = getMassMatrix(v,w,V);
  }
  else
  {
    Info << "Error: hilbertSpacePOD " << hilbertSpacePOD << " is not valid." << endl;
    Info << "NOT CODED YET : dot_product_POD is available for L2 only." << endl;
    abort();
  }
  return matrix_dot_product;
}
// Specialization
template Eigen::MatrixXd dot_product_POD(PtrList<volTensorField>& v, PtrList<volTensorField>& w,
                                const word& hilbertSpacePOD,
                                const double& weightBC,
                                const word& patchBC,
                                const double& weightH1);
template Eigen::MatrixXd dot_product_POD(PtrList<volVectorField>& v, PtrList<volVectorField>& w,
                                const word& hilbertSpacePOD,
                                const double& weightBC,
                                const word& patchBC,
                                const double& weightH1);
template Eigen::MatrixXd dot_product_POD(PtrList<volScalarField>& v, PtrList<volScalarField>& w,
                                const word& hilbertSpacePOD,
                                const double& weightBC,
                                const word& patchBC,
                                const double& weightH1);


//Separation because dot product H1 can't be applied when using volTensorFields
double dot_product_POD(const volScalarField& v, const volScalarField& w,
                                const word& hilbertSpacePOD,
                                const double& weightBC,
                                const word& patchBC,
                                const double& weightH1)
{
  double integral(0);
  if (hilbertSpacePOD == "L2" || hilbertSpacePOD == "dL2")
  {
    integral = dot_product_L2(v,w);
  }
  else if (hilbertSpacePOD == "L2wBC")
  {
    integral = dot_product_L2wBC(v,w,weightBC, patchBC);
  }
  else if (hilbertSpacePOD == "H1")
  {
    integral = dot_product_H1(v,w);
  }
  else if (hilbertSpacePOD == "wH1")
  {
    integral = dot_product_H1(v,w,weightH1);
  }
  else
  {
    Info << "Error: hilbertSpacePOD " << hilbertSpacePOD << " is not valid." << endl;
    Info << "dot_product_POD is available for L2, L2wBC, H1 and wH1 only." << endl;
    abort();
  }
  return integral;
}

double dot_product_POD(const volVectorField& v, const volVectorField& w, const word& hilbertSpacePOD, const double& weightBC, const word& patchBC, const double& weightH1)
{
  double integral(0);
  if (hilbertSpacePOD == "L2")
  {
    integral = dot_product_L2(v,w);
  }
  else if (hilbertSpacePOD == "L2wBC")
  {
    integral = dot_product_L2wBC(v,w,weightBC, patchBC);
  }
  else if (hilbertSpacePOD == "H1")
  {
    integral = dot_product_H1(v,w);
  }
  else if (hilbertSpacePOD == "wH1")
  {
    integral = dot_product_H1(v,w,weightH1);
  }
  else
  {
    Info << "Error: hilbertSpacePOD " << hilbertSpacePOD << " is not valid." << endl;
    Info << "dot_product_POD is available for L2, L2wBC, H1 and wH1 only." << endl;
    abort();
  }
  return integral;
}

double dot_product_POD(const volTensorField& v, const volTensorField& w,
                                const word& hilbertSpacePOD,
                                const double& weightBC,
                                const word& patchBC,
                                const double& weightH1)
{
  double integral(0);
  if (hilbertSpacePOD == "L2")
  {
    integral = dot_product_L2(v,w);
  }
  else if (hilbertSpacePOD == "L2wBC")
  {
    integral = dot_product_L2wBC(v,w,weightBC, patchBC);
  }
  else if ((hilbertSpacePOD == "H1")||(hilbertSpacePOD == "wH1"))
  {
    Info << "Error : dot_product_POD cannot be computed between volTensorFields with the hilbertSpacePOD =" <<hilbertSpacePOD << " ." << endl;
    abort();
  }
  else
  {
    Info << "Error: hilbertSpacePOD " << hilbertSpacePOD << " is not valid." << endl;
    Info << "dot_product_POD is available for L2, L2wBC, H1 and wH1 only." << endl;
    abort();
  }
  return integral;
}

template<typename T>
double dot_product_H1(const T& v, const T& w, const double& weightH1)
{
  return ( dot_product_L2(v,w) + weightH1 * dot_product_L2(fvc::grad(v),fvc::grad(w)) );
}
// Specialization
template double dot_product_H1(const volVectorField& v, const volVectorField& w, const double& weightH1);
template double dot_product_H1(const volScalarField& v, const volScalarField& w, const double& weightH1);

double dot_product_L2(const volVectorField& v, const volVectorField& w)
{
  return fvc::domainIntegrate(v & w).value();

}

double dot_product_L2(const volScalarField& v, const volScalarField& w)
{
  return fvc::domainIntegrate(v * w).value();
}


double dot_product_L2(const volTensorField& v, const volTensorField& w)
{
  return fvc::domainIntegrate(v && w).value();
}

//Dot product at a boundary patch
double dot_product_patch(const Eigen::VectorXd& f1BC_i,const Eigen::VectorXd& f2BC_i, const scalarField& AreaFace, const int& d)
{
   double dot_prod_patch=0.;
   double dx = 1.;
   //Loop for the boundary cells
   for (label k = 0; k < f1BC_i.size()/d; k++)
   {
     double productBC=0.;
     //Loop for the coordinates separated in the Eigen::VectorXd
     for (label i=0;i<d;i++)
     {
       productBC+=f1BC_i(k+i*f1BC_i.size()/d) * f2BC_i(k+i*f1BC_i.size()/d);
     }
     dot_prod_patch+=productBC* AreaFace[k] * dx;
   }
  return dot_prod_patch;
}

template<typename T>
// Dot product at the boundary
double dot_product_boundary(const T& v,const T& w, const word& patchBC)
{
  T f1 = v;
  T f2 = w;
  double dot_prod_boundary=0.;

  int NBC = f1.boundaryField().size();
  List<Eigen::VectorXd> f1BC = Foam2Eigen::field2EigenBC(f1);
  List<Eigen::VectorXd> f2BC = Foam2Eigen::field2EigenBC(f2);

  label ind = f1.mesh().boundaryMesh().findPatchID(patchBC);
  //Dimension of the field to separate the coordinates after Eigen2Field
  int d = dimensionField(f1);
  //List of indexes of the considered boundary patches
  List<label> l;

  if (ind==-1)
  {
    for (label i=0;i<NBC;i++)
    {
      l.append(i);
    }
  }

  else
  {
    l.append(ind);
  }
  //Loop for the boundary patches
  for (label g = 0; g < l.size(); g++)
  {
    // Create a scalar field with area value
    scalarField AreaFace = f1.mesh().magSf().boundaryField()[l[g]];
    // Compute the dot product at each patch
    dot_prod_boundary+=dot_product_patch(f1BC[l[g]],f2BC[l[g]],AreaFace,d);
  }
  return dot_prod_boundary;
}
template double dot_product_boundary(const volTensorField& v, const volTensorField& w, const word& patchBC);
template double dot_product_boundary(const volVectorField& v, const volVectorField& w, const word& patchBC);
template double dot_product_boundary(const volScalarField& v, const volScalarField& w, const word& patchBC);

template<typename T>
// Dot product L2wBC : uses dot_product_L2 at the interior and dot_product_boundary at the boundary with a weight weightBC
double dot_product_L2wBC(const T& v, const T& w, const double& weightBC, const word& patchBC)
{
  double dot_prod;

  //Compute the dot product at the boundary
  double dp_boundary = dot_product_boundary(v,w,patchBC);

  //Add the dot product at the boundary to the L2 one inside of the domain
  dot_prod=dp_boundary*weightBC+dot_product_L2(v,w);
  return dot_prod;
}

template double dot_product_L2wBC(const volTensorField& v, const volTensorField& w, const double& weightBC, const word& patchBC);
template double dot_product_L2wBC(const volVectorField& v, const volVectorField& w, const double& weightBC, const word& patchBC);
template double dot_product_L2wBC(const volScalarField& v, const volScalarField& w, const double& weightBC, const word& patchBC);

template<typename T>
double norm_L2wBC(const T& v, const double& weightBC, const word& patchBC)
{
  return std::sqrt(dot_product_L2wBC(v,v, weightBC, patchBC));
}

template double norm_L2wBC(const volTensorField& v, const double& weightBC, const word& patchBC);
template double norm_L2wBC(const volVectorField& v, const double& weightBC, const word& patchBC);
template double norm_L2wBC(const volScalarField& v, const double& weightBC, const word& patchBC);

template<typename T>
double norm_POD(const T& v, const word& hilbertSpacePOD, const double& weightBC, const word& patchBC, const double& weightH1)
{
  return std::sqrt(dot_product_POD(v,v, hilbertSpacePOD, weightBC, patchBC, weightH1));
}
// Specialization
template double norm_POD(const volVectorField& v, const word& hilbertSpacePOD, const double& weightBC, const word& patchBC, const double& weightH1);
template double norm_POD(const volScalarField& v, const word& hilbertSpacePOD, const double& weightBC, const word& patchBC, const double& weightH1);

template<typename T>
double L2Norm(const T v)
{
  return Foam::sqrt(dot_product_L2(v,v));
}

// Specialization
template double L2Norm(const volTensorField& v);
template double L2Norm(const volVectorField& v);
template double L2Norm(const volScalarField& v);
template double L2Norm(const volTensorField v);
template double L2Norm(const volVectorField v);
template double L2Norm(const volScalarField v);
template double L2Norm(const tmp<volTensorField> v);
template double L2Norm(const tmp<volVectorField> v);
template double L2Norm(const tmp<volScalarField> v);


template<>
double LinfNorm(GeometricField<scalar, fvPatchField, volMesh>& field)
{
    double a;
    a = Foam::max(Foam::sqrt(field.internalField() *
                             field.internalField())).value();
    return a;
}

template<>
double LinfNorm(GeometricField<vector, fvPatchField, volMesh>& field)
{
    double a;
    Info << "LinfNorm(GeometricField<vector, fvPatchField, volMesh>& field) is still to be implemented"
         << endl;
    exit(12);
    return a;
}

template<>
double H1Seminorm(GeometricField<scalar, fvPatchField, volMesh>& field)
{
    double a;
    a = Foam::sqrt(fvc::domainIntegrate(fvc::grad(field) & fvc::grad(
                                            field)).value());
    return a;
}

template<>
double H1Seminorm(GeometricField<vector, fvPatchField, volMesh>& field)
{
    double a;
    a = Foam::sqrt(fvc::domainIntegrate(fvc::grad(field)
                                        && fvc::grad(field)).value());
    return a;
}

template<class Type, template<class> class PatchField, class GeoMesh>
double frobNorm(GeometricField<Type, PatchField, GeoMesh>& field)
{
    double norm(0);
    Eigen::VectorXd vF = Foam2Eigen::field2Eigen(field);
    norm = vF.norm();
    return norm;
}
template double frobNorm(GeometricField<scalar, fvPatchField, volMesh>& field);
template double frobNorm(GeometricField<scalar, fvsPatchField, surfaceMesh>& field);
template double frobNorm(GeometricField<vector, fvPatchField, volMesh>& field);

double L2normOnPatch(fvMesh& mesh, volScalarField& field,
                     word patch)
{
    double L2 = 0;
    //Access the mesh information for the boundary
    label patchID = mesh.boundaryMesh().findPatchID(patch);
    const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
    //List of cells close to a boundary
    const labelUList& faceCells = cPatch.faceCells();
    forAll(cPatch, faceI)
    {
        //id of the owner cell having the face
        label faceOwner = faceCells[faceI] ;
        scalar faceArea = mesh.magSf().boundaryField()[patchID][faceI];
        L2 += faceArea * field[faceOwner] * field[faceOwner];
    }

    return Foam::sqrt(L2);
}

double L2productOnPatch(fvMesh& mesh, List<scalar>& field1,
                        List<scalar>& field2, word patch)
{
    M_Assert(field1.size() == field2.size(),
             "The two fields do not have the same size, code will abort");
    double L2 = 0;
    //Access the mesh information for the boundary
    label patchID = mesh.boundaryMesh().findPatchID(patch);
    const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
    M_Assert(field1.size() == cPatch.size(),
             "The filds must have the same size of the patch, code will abort");
    forAll(cPatch, faceI)
    {
        scalar faceArea = mesh.magSf().boundaryField()[patchID][faceI];
        L2 += faceArea * field1[faceI] * field2[faceI];
    }

    return L2;
}

double LinfNormOnPatch(fvMesh& mesh, volScalarField& field,
                       word patch)
{
    double Linf = 0;
    //Access the mesh information for the boundary
    label patchID = mesh.boundaryMesh().findPatchID(patch);
    const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
    //List of cells close to a boundary
    const labelUList& faceCells = cPatch.faceCells();
    forAll(cPatch, faceI)
    {
        label faceOwner = faceCells[faceI] ;

        if (faceI == 0)
        {
            Linf = std::abs(field[faceOwner]);
        }
        else if (std::abs(field[faceOwner]) > Linf)
        {
            Linf = std::abs(field[faceOwner]);
        }
    }

    return Linf;
}

double integralOnPatch(fvMesh& mesh, volScalarField& field,
                       word patch)
{
    double integral = 0;
    //Access the mesh information for the boundary
    label patchID = mesh.boundaryMesh().findPatchID(patch);
    const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
    //List of cells close to a boundary
    const labelUList& faceCells = cPatch.faceCells();
    forAll(cPatch, faceI)
    {
        //id of the owner cell having the face
        label faceOwner = faceCells[faceI] ;
        scalar faceArea = mesh.magSf().boundaryField()[patchID][faceI];
        integral += faceArea * field[faceOwner];
    }

    return integral;
}

double integralOnPatch(fvMesh& mesh, List<scalar> field,
                       word patch)
{
    double integral = 0;
    //Access the mesh information for the boundary
    label patchID = mesh.boundaryMesh().findPatchID(patch);
    const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
    //List of cells close to a boundary
    const labelUList& faceCells = cPatch.faceCells();
    int patchSize = mesh.magSf().boundaryField()[patchID].size();
    std::string message = "The input list (size = " + std::to_string(field.size())
                          + ") must have the same size of the patch mesh (size = "
                          + std::to_string(patchSize) + ")";
    M_Assert( patchSize == field.size(), message.c_str());
    forAll(cPatch, faceI)
    {
        //id of the owner cell having the face
        label faceOwner = faceCells[faceI] ;
        scalar faceArea = mesh.magSf().boundaryField()[patchID][faceI];
        integral += faceArea * field[faceI];
    }

    return integral;
}


}
