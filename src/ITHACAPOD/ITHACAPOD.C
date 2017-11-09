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
/// source file for the ITHACAPOD class

#include "ITHACAPOD.H"

void ITHACAPOD::getModes(PtrList<volVectorField>& snapshotsU, PtrList<volVectorField>& modes, bool podex, bool supex, bool sup, int nmodes)
{

  if (nmodes == 0)
  {
    nmodes = snapshotsU.size();
  }

  if (podex == 0)
  {
    PtrList<volVectorField> Bases;
    Bases.resize(nmodes);
    modes.resize(nmodes);
    Eigen::MatrixXd _corMatrix;
    scalarField eigenValues(snapshotsU.size());
    scalarField cumEigenValues(snapshotsU.size());
    List<scalarField> eigenVector(snapshotsU.size());

    for (label i = 0; i < snapshotsU.size(); i++)
    {
      eigenVector[i].setSize(snapshotsU.size());
    }
    _corMatrix = ITHACAPOD::corMatrix(snapshotsU);
    Info << "####### Performing the POD decomposition for " << snapshotsU[0].name() << " #######" << endl;
        Eigen::EigenSolver<Eigen::MatrixXd> es(_corMatrix);
    Info << "####### End of the POD decomposition for " << snapshotsU[0].name() << " #######" << endl;
    Eigen::MatrixXd eigenVectoreig = es.eigenvectors().real();

    cumEigenValues[0] = es.eigenvalues().real()(0) / es.eigenvalues().real().sum();
    eigenValues[0] = es.eigenvalues().real()(0) / es.eigenvalues().real().sum();
    for (label i = 1; i < snapshotsU.size(); i++)
    {
      cumEigenValues[i] = cumEigenValues[i - 1] + es.eigenvalues().real()(i) / es.eigenvalues().real().sum();
      eigenValues[i] = es.eigenvalues().real()(i) / es.eigenvalues().real().sum();
    }

    for (label i = 0; i < snapshotsU.size(); i++)
    {
      for (label k = 0; k < snapshotsU.size(); k++)
      {
        eigenVector[i][k] = eigenVectoreig(k, i);
      }
    }

    volVectorField tmb_bu
    (
      snapshotsU[0].name(),
      snapshotsU[0]
    );
    for (label i = 0; i < nmodes; i++)
    {
      volVectorField tmb_bu(eigenVector[i][0]*snapshotsU[0]);
      for (label k = 0; k < snapshotsU.size(); k++)
      {
        tmb_bu += eigenVector[i][k] * snapshotsU[k];
      }
      Bases.set(i, tmb_bu);
      Info << "creating the bases " << i << " for " << snapshotsU[0].name() << endl;
    }
    ITHACAPOD::normalizeBases(Bases);
    for (label i = 0; i < Bases.size(); i++)
    {
      volVectorField tmp
      (
        snapshotsU[0].name(),
        Bases[i]
      );
      modes.set(i, tmp);
    }
    Info << "####### Saving the POD bases for " << snapshotsU[0].name() << " #######" << endl;
    ITHACAPOD::exportBases(modes, snapshotsU, sup);
    ITHACAPOD::exportEigenvalues(eigenValues, snapshotsU[0].name());
    ITHACAPOD::exportcumEigenvalues(cumEigenValues, snapshotsU[0].name());
  }
  else
  {
    Info << "Reading the existing modes" << endl;
    if (sup == 1)
    {
      ITHACAstream::read_fields(modes, snapshotsU[0], "./ITHACAoutput/supremizer/");
    }
    else
    {
      ITHACAstream::read_fields (modes, snapshotsU[0], "./ITHACAoutput/POD/");
    }
  }
}

void ITHACAPOD::getModes(PtrList<volScalarField>& snapshotsP, PtrList<volScalarField>& modes, bool podex, bool supex, bool sup, int nmodes)
{
  if (nmodes == 0)
  {
    nmodes = snapshotsP.size();
  }
  if (podex == 0)
  {
    PtrList<volScalarField> Bases;
    Bases.resize(nmodes);
    modes.resize(nmodes);
    Eigen::MatrixXd _corMatrix;
    scalarField eigenValues(snapshotsP.size());
    scalarField cumEigenValues(snapshotsP.size());
    List<scalarField> eigenVector(snapshotsP.size());
    for (label i = 0; i < snapshotsP.size(); i++)
    {
      eigenVector[i].setSize(snapshotsP.size());
    }
    _corMatrix = ITHACAPOD::corMatrix(snapshotsP);
    Info << "####### Performing the POD decomposition for " << snapshotsP[0].name() << " #######" << endl;
    Eigen::EigenSolver<Eigen::MatrixXd> es(_corMatrix);
    Info << "####### End of the POD decomposition for " << snapshotsP[0].name() << " #######" << endl;
    Eigen::MatrixXd eigenVectoreig = es.eigenvectors().real();

    cumEigenValues[0] = es.eigenvalues().real()(0) / es.eigenvalues().real().sum();
    eigenValues[0] = es.eigenvalues().real()(0) / es.eigenvalues().real().sum();
    for (label i = 1; i < snapshotsP.size(); i++)
    {
      cumEigenValues[i] = cumEigenValues[i - 1] + es.eigenvalues().real()(i) / es.eigenvalues().real().sum();
      eigenValues[i] = es.eigenvalues().real()(i) / es.eigenvalues().real().sum();
    }

    for (label i = 0; i < snapshotsP.size(); i++)
    {
      for (label k = 0; k < snapshotsP.size(); k++)
      {
        eigenVector[i][k] = eigenVectoreig(k, i);
      }
    }

    volScalarField tmb_bu
    (
      snapshotsP[0].name(),
      snapshotsP[0]
    );
    for (label i = 0; i < nmodes; i++)
    {
      volScalarField tmb_bu(eigenVector[i][0]*snapshotsP[0]);
      for (label k = 0; k < snapshotsP.size(); k++)
      {
        tmb_bu += eigenVector[i][k] * snapshotsP[k];
      }
      Bases.set(i, tmb_bu);
      Info << "creating the bases " << i << " for " << snapshotsP[0].name() << endl;
    }
    ITHACAPOD::normalizeBases(Bases);
    for (label i = 0; i < Bases.size(); i++)
    {
      volScalarField tmp
      (
        snapshotsP[0].name(),
        Bases[i]
      );
      modes.set(i, tmp);
    }
    Info << "####### Saving the POD bases for " << snapshotsP[0].name() << " #######" << endl;
    ITHACAPOD::exportBases(modes, snapshotsP, sup);
    ITHACAPOD::exportEigenvalues(eigenValues, snapshotsP[0].name());
    ITHACAPOD::exportcumEigenvalues(cumEigenValues, snapshotsP[0].name());
  }
  else
  {
    Info << "Reading the existing modes" << endl;
    if (sup == 1)
    {
      ITHACAstream::read_fields(modes, "Usup", "./ITHACAoutput/supremizer/");
    }
    else
    {
      ITHACAstream::read_fields (modes, snapshotsP[0], "./ITHACAoutput/POD/");
    }
  }
}

/// Normalize the bases
void ITHACAPOD::normalizeBases(PtrList<volScalarField>& Bases)
{
  scalar magSumSquare;
  for (label j = 0; j < Bases.size(); j++)
  {
    magSumSquare = Foam::sqrt(fvc::domainIntegrate(Bases[j] * Bases[j]).value());
    if (magSumSquare > SMALL)
    {
      Bases[j] /= magSumSquare;
    }
    Bases[j].correctBoundaryConditions();
  }
}

void ITHACAPOD::normalizeBases(PtrList<volVectorField>& Bases)
{
  scalar magSumSquare;
  for (label j = 0; j < Bases.size(); j++)
  {
    magSumSquare = Foam::sqrt(fvc::domainIntegrate(Bases[j] & Bases[j]).value());
    if (magSumSquare > SMALL)
    {
      Bases[j] /= magSumSquare;
    }
    Bases[j].correctBoundaryConditions();
  }
}

/// Normalize the bases
void ITHACAPOD::normalizeBases(PtrList<volVectorField>& BasesU, PtrList<volScalarField>& BasesP)
{
  scalar magSumSquare;
  for (label j = 0; j < BasesU.size(); j++)
  {
    magSumSquare = Foam::sqrt(fvc::domainIntegrate(BasesU[j] & BasesU[j]).value());
    if (magSumSquare > SMALL)
    {
      BasesP[j] /= magSumSquare;
    }
    BasesP[j].correctBoundaryConditions();
  }
}


/// Construct the Correlation Matrix for Scalar Field
Eigen::MatrixXd ITHACAPOD::corMatrix(PtrList<volScalarField>& snapshots)
{
  Info << "########## Filling the correlation matrix for " << snapshots[0].name() << "##########" << endl;
  Eigen::MatrixXd matrix( snapshots.size(), snapshots.size());
  for (label i = 0; i < snapshots.size(); i++)
  {
    Info << "Filling row " << i << " of the " << snapshots[0].name() << " correlation matrix" << endl;
    for (label j = 0; j <= i; j++)
    {
      matrix(i, j) = fvc::domainIntegrate(snapshots[i] * snapshots[j]).value();
    }
  }
  for (label i = 1; i < snapshots.size(); i++)
  {
    for (label j = 0; j < i; j++)
    {
      matrix(j, i) = matrix(i, j);
    }
  }
  return matrix;
}


/// Construct the Correlation Matrix for Vector Field
Eigen::MatrixXd ITHACAPOD::corMatrix(PtrList<volVectorField>& snapshots)
{
  Info << "########## Filling the correlation matrix for " << snapshots[0].name() << "##########" << endl;
  Eigen::MatrixXd matrix( snapshots.size(), snapshots.size());
  for (label i = 0; i < snapshots.size(); i++)
  {
    Info << "Filling row " << i << " of the " << snapshots[0].name() << " correlation matrix" << endl;
    for (label j = 0; j <= i; j++)
    {
      matrix(i, j) = fvc::domainIntegrate(snapshots[i] & snapshots[j]).value();
    }
  }
  for (label i = 1; i < snapshots.size(); i++)
  {
    for (label j = 0; j < i; j++)
    {
      matrix(j, i) = matrix(i, j);
    }
  }
  return matrix;
}

/// Export the Bases
void ITHACAPOD::exportBases(PtrList<volVectorField>& s, PtrList<volVectorField>& _snapshots, bool sup)
{
  if (sup)
  {
    fileName fieldname;
    for (label i = 0; i < s.size(); i++)
    {
      mkDir("./ITHACAoutput/supremizer/" + name(i + 1));
      fieldname = "./ITHACAoutput/supremizer/" + name(i + 1) + "/" + _snapshots[i].name();
      OFstream os(fieldname);
      _snapshots[i].writeHeader(os);
      os << s[i] << endl;
    }

  }
  else
  {
    fileName fieldname;
    for (label i = 0; i < s.size(); i++)
    {
      mkDir("./ITHACAoutput/POD/" + name(i + 1));
      fieldname = "./ITHACAoutput/POD/" + name(i + 1) + "/" + _snapshots[i].name();
      OFstream os(fieldname);
      _snapshots[i].writeHeader(os);
      os << s[i] << endl;
    }
  }
}


/// Export the Bases
void ITHACAPOD::exportBases(PtrList<volScalarField>& s, PtrList<volScalarField>& _snapshots, bool sup)
{
  if (sup)
  {
    fileName fieldname;
    for (label i = 0; i < s.size(); i++)
    {
      mkDir("./ITHACAoutput/supremizer/" + name(i + 1));
      fieldname = "./ITHACAoutput/supremizer/" + name(i + 1) + "/" + _snapshots[i].name();
      OFstream os(fieldname);
      _snapshots[i].writeHeader(os);
      os << s[i] << endl;
    }

  }
  else
  {
    fileName fieldname;
    for (label i = 0; i < s.size(); i++)
    {
      mkDir("./ITHACAoutput/POD/" + name(i + 1));
      fieldname = "./ITHACAoutput/POD/" + name(i + 1) + "/" + _snapshots[i].name();
      OFstream os(fieldname);
      _snapshots[i].writeHeader(os);
      os << s[i] << endl;
    }
  }
}

void ITHACAPOD::exportEigenvalues(scalarField Eigenvalues, fileName name, bool sup)
{
  if (sup)
  {
    fileName fieldname;
    mkDir("./ITHACAoutput/supremizer/");
    fieldname = "./ITHACAoutput/supremizer/Eigenvalues_" + name;
    OFstream os(fieldname);
    for (label i = 0; i < Eigenvalues.size(); i++)
    {
      os << Eigenvalues[i] << endl;
    }
  }
  else {
    fileName fieldname;
    mkDir("./ITHACAoutput/POD/");
    fieldname = "./ITHACAoutput/POD/Eigenvalues_" + name;
    OFstream os(fieldname);
    for (label i = 0; i < Eigenvalues.size(); i++)
    {
      os << Eigenvalues[i] << endl;
    }
  }
}

void ITHACAPOD::exportcumEigenvalues(scalarField cumEigenvalues, fileName name, bool sup)
{
  if (sup)
  {
    fileName fieldname;
    mkDir("./ITHACAoutput/supremizer/");
    fieldname = "./ITHACAoutput/supremizer/cumEigenvalues_" + name;
    OFstream os(fieldname);
    for (label i = 0; i < cumEigenvalues.size(); i++)
    {
      os << cumEigenvalues[i] << endl;
    }
  }
  else {
    fileName fieldname;
    mkDir("./ITHACAoutput/POD/");
    fieldname = "./ITHACAoutput/POD/cumEigenvalues_" + name;
    OFstream os(fieldname);
    for (label i = 0; i < cumEigenvalues.size(); i++)
    {
      os << cumEigenvalues[i] << endl;
    }
  }
}

// * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * * * //







