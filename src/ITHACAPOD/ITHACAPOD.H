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
    ITHACAPOD

Description
    Implementation of a POD decomposition of a general field

SourceFiles
    ITHACAPOD.C

\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ITHACAPOD class.
/// \dir
/// Directory containing the header and source files for the ITHACAPOD class.

#ifndef ITHACAPOD_H
#define ITHACAPOD_H

#include "fvCFD.H"
#include "ITHACAutilities.H"
#include "ITHACAstream.H"
#include "../thirdparty/Eigen/Eigen/Eigen"

/*---------------------------------------------------------------------------*\
                        Class reductionProblem Declaration
\*---------------------------------------------------------------------------*/

/// Class of the computation of the POD, it exploits the method of snapshots.
class ITHACAPOD
{

public:
	// Public Members
		/// Compute the bases or read them for a vector field
        ///
        /// @param[in]  snapshotsU  a PtrList of volVectorField snapshots.
        /// @param[out] modesU      a PtrList of volVectorField where the modes are stored (it must be passed empty).
        /// @param[in]  podex       boolean variable 1 if the modes have been already computed and stored (in this case the function is reading them) 0 elsewhere.
        /// @param[in]  supex       boolean variable 1 if the supremizer modes have been already computed and stored (in this case the function is reading them) 0 elsewhere.
		/// @param[in]  sup         boolean variable 1 if you want to compute the supremizer modes 0 elsewhere.
        /// @param[in]  nmodes      int variable to set the number of modes to be stored, if set to 0 the maximum number of modes will computed.
        ///
        static void getModes(PtrList<volVectorField>& snapshotsU, PtrList<volVectorField>& modesU, bool podex, bool supex=0, bool sup=0, int nmodes=0);
		
        /// Compute the bases or read them for a scalar field
        /// 
        /// @param[in]  snapshotsP  a PtrList of volScalarField snapshots.
        /// @param[out] modesP      a PtrList of volScalarField where the modes are stored (it must be passed empty).
        /// @param[in]  podex       boolean variable 1 if the modes have been already computed and stored (in this case the function is reading them) 0 elsewhere.
        /// @param[in]  supex       boolean variable 1 if the supremizer modes have been already computed and stored (in this case the function is reading them) 0 elsewhere.
        /// @param[in]  sup         boolean variable 1 if you want to compute the supremizer modes 0 elsewhere.
        /// @param[in]  nmodes      int variable to set the number of modes to be stored, if set to 0 the maximum number of modes will computed.
        /// 
		static void getModes(PtrList<volScalarField>& snapshotsP, PtrList<volScalarField>& modesP, bool podex, bool supex=0, bool sup=0, int nmodes=0);
		
        /// Computes the correlation matrix given a vector field snapshot Matrix using the L2 norm
        /// 
        /// @param[in] snapshots    a PtrList of volVectorField snapshots.
        /// 
        /// @return    the Eigen::MatrixXd correlation matrix.
        /// 
		static Eigen::MatrixXd corMatrix(PtrList<volVectorField>& snapshots);
		
        /// Computes the correlation matrix given a scalar field snapshot Matrix using the L2 norm.
        /// 
        /// @param[in] snapshots    a PtrList of volScalarField snapshots.
        /// 
        /// @return    the Eigen::MatrixXd correlation matrix.
        /// 
		static Eigen::MatrixXd corMatrix(PtrList<volScalarField>& snapshots);
		
        /// Export the basis for a vector field into the ITHACAOutput/POD or ITHACAOutput/supremizer
        /// 
        /// @param[in] s        a PtrList of volVectorField where the snapshots associated with the bases are stored.
        /// @param[in] bases    a PtrList of volVectorField where the bases are stored.
        /// @param[in] sup      a boolean variable 1 if you want to export the supremizer bases (in ITHACAOutput/supremizer) 0 elsewhere (Default is 0).
        /// 
		static void exportBases(PtrList<volVectorField>& s, PtrList<volVectorField>& bases, bool sup = 0);
        
        /// Export the basis for a scalar field into the ITHACAOutput/POD or ITHACAOutput/supremizer
        /// 
        /// @param[in] s        a PtrList of volScalarField where the snapshots associated with the bases are stored.
        /// @param[in] bases    a PtrList of volScalarField where the bases are stored.
        /// @param[in] sup      a boolean variable 1 if you want to export the supremizer bases (in ITHACAOutput/supremizer) 0 elsewhere (Default is 0).
        /// 
		static void exportBases(PtrList<volScalarField>& s, PtrList<volScalarField>& bases, bool sup = 0);
		
        /// Normalize the Bases for a vector field
        /// 
        /// @param[in] s    a PtrList of volVectorField where the bases are stored.
        /// 
		static void normalizeBases(PtrList<volVectorField>& s);
		
        /// Normalize the Bases for a scalar field
        /// 
        /// @param[in] s    a PtrList of volScalarField where the bases are stored.
        /// 
		static void normalizeBases(PtrList<volScalarField>& s);
		
        /// Normalize a volScalarField field using a volVectorField
        /// 
        /// @param[in] s    a PtrList of volVectorField where the bases used for normalization are stored.
        /// @param[in] p    a PtrList of volScalarField  where the bases are stored.
        /// 
		static void normalizeBases(PtrList<volVectorField>& s, PtrList<volScalarField>& p);
		
        /// Export the eigenvalues as a txt file
        /// 
        /// @param[in] eigenValues  a scalarField of eigenvalues.
        /// @param[in] name         a fileName to indicate the name of the txt file.         
        /// @param[in] sup          a boolean variable 1 if you want to export the eigenValues (in ITHACAOutput/supremizer) 0 elsewhere (Default is 0).
        /// 
		static void exportEigenvalues(scalarField eigenValues, fileName name, bool sup = 0);
		
        /// Export the cumulative eigenvalues as a txt file
        /// 
        /// @param[in] cumEigenValues   a scalarField of cumEigenValues.
        /// @param[in] name             a fileName to indicate the name of the txt file.         
        /// @param[in] sup              a boolean variable 1 if you want to export the cumEigenValues (in ITHACAOutput/supremizer) 0 elsewhere (Default is 0).
        /// 
		static void exportcumEigenvalues(scalarField cumEigenValues,  fileName name, bool sup = 0);

protected:

};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif



