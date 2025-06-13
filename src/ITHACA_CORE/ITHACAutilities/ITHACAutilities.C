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

#include "ITHACAutilities.H"
#include "ITHACAstream.H"
#include "ITHACAparameters.H"
#include "turbulentTransportModel.H"

/// \file
/// Source file of the ITHACAutilities namespace.

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace ITHACAutilities
{

Eigen::MatrixXd rand(label rows, label cols, double min,
                     double max)
{
    std::srand(static_cast<long unsigned int>
               (std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    Eigen::MatrixXd matr = Eigen::MatrixXd::Random(rows, cols);
    matr = (matr.array() + 1) / 2;
    matr = matr.array() * (max - min);
    matr = matr.array() + min;
    return matr;
}

Eigen::MatrixXd rand(label rows, Eigen::MatrixXd minMax)
{
    std::srand(static_cast<long unsigned int>
               (std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    label cols = minMax.rows();
    Eigen::MatrixXd matr = Eigen::MatrixXd::Random(rows, cols);
    matr = (matr.array() + 1) / 2;

    for (label i = 0; i < cols; i++)
    {
        matr.col(i) = matr.col(i).array() * (minMax(i, 1) - minMax(i, 0));
        matr.col(i) = matr.col(i).array() + (minMax(i, 0));
    }

    return matr;
}


bool isInteger(double ratio)
{
    bool checkResult = 0;

    if (abs(round(ratio) - ratio) < std::sqrt(SMALL))
    {
        checkResult = true;
    }
    else
    {
        checkResult = false;
    }

    return checkResult;
}

bool isTurbulent()
{
    bool checkTurb;
    ITHACAparameters* para = ITHACAparameters::getInstance();
    auto& tur =
        para->mesh.lookupObject<incompressible::turbulenceModel>("turbulenceProperties");

    if (tur.type() == "Stokes" || tur.type() == "Maxwell"
            || tur.type() == "laminarModel")
    {
        checkTurb = false;
    }
    else
    {
        checkTurb = true;
    }

    return checkTurb;
}

template<typename T>
List<T> combineList(List<List<T >> & doubleList)
{
    List<T> a = ListListOps::combine<List<T >> (doubleList,
                accessOp<List<T >> ());
#if OPENFOAM >= 1812
    inplaceUniqueSort(a);
#else
    labelList order;
    uniqueOrder(a, order);
    List<T> b(order.size());

    for (label i = 0; i < order.size(); ++i)
    {
        b[i] = a[order[i]];
    }
    a.resize(order.size());
    a = b;
#endif
    return a;
}

template List<label> combineList(List<List<label >> & doubleList);


// Using the Eigen library, using the SVD decomposition method to solve the
// matrix pseudo-inverse, the default error er is 0
Eigen::MatrixXd pinv_eigen_based(Eigen::MatrixXd& origin, const float er)
{
    // perform svd decomposition
    Eigen::JacobiSVD<Eigen::MatrixXd> svd_holder(origin,
            Eigen::ComputeThinU | Eigen::ComputeThinV);
    // Build SVD decomposition results
    Eigen::MatrixXd U = svd_holder.matrixU();
    Eigen::MatrixXd V = svd_holder.matrixV();
    Eigen::MatrixXd D = svd_holder.singularValues();
    // Build the S matrix
    Eigen::MatrixXd S(V.cols(), U.cols());
    S.setZero();

    for (unsigned int i = 0; i < D.size(); ++i)
    {
        if (D(i, 0) > er)
        {
            S(i, i) = 1 / D(i, 0);
        }
        else
        {
            S(i, i) = 0;
        }
    }

    return V * S * U.transpose();
}


Eigen::MatrixXd invertMatrix(Eigen::MatrixXd& matrixToInvert,
                             const word inversionMethod)
{
    Info << "Inversion method : " << inversionMethod << endl;

    if (inversionMethod == "pinv_eigen_based")
    {
        return pinv_eigen_based(matrixToInvert);
    }
    else if (inversionMethod == "direct")
    {
        return matrixToInvert.inverse();
    }
    else if (inversionMethod == "fullPivLu")
    {
        return matrixToInvert.fullPivLu().inverse();
    }
    else if (inversionMethod == "partialPivLu")
    {
        return matrixToInvert.partialPivLu().inverse();
    }
    else if (inversionMethod == "householderQr")
    {
        return matrixToInvert.householderQr().solve(Eigen::MatrixXd::Identity(
                    matrixToInvert.rows(), matrixToInvert.cols()));
    }
    else if (inversionMethod == "colPivHouseholderQr")
    {
        return matrixToInvert.colPivHouseholderQr().inverse();
    }
    else if (inversionMethod == "fullPivHouseholderQr")
    {
        return matrixToInvert.fullPivHouseholderQr().inverse();
    }
    else if (inversionMethod == "completeOrthogonalDecomposition")
    {
        return matrixToInvert.completeOrthogonalDecomposition().pseudoInverse();
    }
    else if (inversionMethod == "jacobiSvd")
    {
        return matrixToInvert.jacobiSvd(Eigen::ComputeThinU |
                                        Eigen::ComputeThinV).solve(Eigen::MatrixXd::Identity(matrixToInvert.rows(),
                                                matrixToInvert.cols()));
    }
    else if (inversionMethod == "llt")
    {
        return matrixToInvert.llt().solve(Eigen::MatrixXd::Identity(
                                              matrixToInvert.rows(), matrixToInvert.cols()));
    }
    else if (inversionMethod == "ldlt")
    {
        Eigen::LLT<Eigen::MatrixXd> lltOfA(matrixToInvert);
        return lltOfA.solve(Eigen::MatrixXd::Identity(matrixToInvert.rows(),
                            matrixToInvert.cols()));
    }
    else if (inversionMethod == "bdcSvd")
    {
        return matrixToInvert.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(
                   Eigen::MatrixXd::Identity(matrixToInvert.rows(), matrixToInvert.cols()));
    }
    else
    {
        Info << "Unkwown inversion method, solving with : completeOrthogonalDecomposition"
             << endl;
        return matrixToInvert.completeOrthogonalDecomposition().pseudoInverse();
    }
}

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

template<typename T>
void setToZero(T& f1)
{
  multField(f1, 0.0);
}
template void setToZero(volScalarField& f1);
template void setToZero(volVectorField& f1);
template void setToZero(volTensorField& f1);

volTensorField tensorFieldProduct(const volScalarField& coef, const volTensorField& S)
{
  return (coef * S);
}

volTensorField tensorFieldProduct(const volTensorField& coef, const volTensorField& S)
{
  return (coef & S);
}

//Compute the dimensions :

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
double norm_L2(const T v)
{
  return std::sqrt(dot_product_L2(v,v));
}

// Specialization
template double norm_L2(const volTensorField& v);
template double norm_L2(const volVectorField& v);
template double norm_L2(const volScalarField& v);
template double norm_L2(const volTensorField v);
template double norm_L2(const volVectorField v);
template double norm_L2(const volScalarField v);
template double norm_L2(const tmp<volTensorField> v);
template double norm_L2(const tmp<volVectorField> v);
template double norm_L2(const tmp<volScalarField> v);




}
