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
    Example of the reconstruction of a non-linear function using the DEIM
SourceFiles
    08DEIM.C
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "simpleControl.H"
#include "fvMeshSubset.H"
#include "ITHACAutilities.H"
#include "Foam2Eigen.H"
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include "hyperReduction.templates.H"
#include <chrono>
class DEIM_function : public HyperReduction<PtrList<volScalarField> & >
{
    public:
        using HyperReduction::HyperReduction;
        static volScalarField evaluate_expression(volScalarField& S, Eigen::MatrixXd mu)
        {
            volScalarField yPos = S.mesh().C().component(vector::Y).ref();
            volScalarField xPos = S.mesh().C().component(vector::X).ref();

            for (auto i = 0; i < S.size(); i++)
            {
                S[i] = std::exp(- 2 * std::pow(xPos[i] - mu(0) - 1,
                                               2) - 2 * std::pow(yPos[i] - mu(1) - 0.5, 2));
            }

            return S;
        }

        Eigen::VectorXd onlineCoeffs(Eigen::MatrixXd mu)
        {
            theta.resize(nodePoints().size());
            auto f = evaluate_expression(subField(), mu);

            for (int i = 0; i < nodePoints().size(); i++)
            {
                // double on_coeff = f[localMagicPoints[i]];
                theta(i) = f[localNodePoints[i]];
            }

            return theta;
        }

        Eigen::VectorXd theta;
        PtrList<volScalarField> fields;
        autoPtr<volScalarField> subField;
};

int main(int argc, char* argv[])
{
    // Read parameters from ITHACAdict file
#include "setRootCase.H"
    Foam::Time runTime(Foam::Time::controlDictName, args);
    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );
    ITHACAparameters* para = ITHACAparameters::getInstance(mesh, runTime);
    int NDEIM = para->ITHACAdict->lookupOrDefault<int>("NDEIM", 15);
    simpleControl simple(mesh);
#include "createFields.H"
    // List of volScalarField where the snapshots are stored
    PtrList<volScalarField> Sp;
    // Non linear function
    volScalarField S
    (
        IOobject
        (
            "S",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, -1, 1, 0, 0, 0), 0)
    );
    // Parameters used to train the non-linear function
    Eigen::MatrixXd pars = ITHACAutilities::rand(100, 2, -0.5, 0.5);

    // To compare with the same parameters, save data and then load it
    // cnpy::load(pars, "./random100.npy");
    // cnpy::save(pars, "./random100.npy");

    // Perform the offline phase
    for (int i = 0; i < 100; i++)
    {
        DEIM_function::evaluate_expression(S, pars.row(i));
        Sp.append((S).clone());
        ITHACAstream::exportSolution(S, name(i + 1), "./ITHACAoutput/Offline/");
    }

    Eigen::MatrixXd snapshotsModes;
    // To use the DEIMmodes method from ITHACAPOD :
    Eigen::VectorXd normalizingWeights = ITHACAutilities::getMassMatrixFV(
            Sp[0]).array();
    PtrList<volScalarField> modes = ITHACAPOD::DEIMmodes(Sp, NDEIM, "GappyDEIM",
                                    S.name());
    snapshotsModes = Foam2Eigen::PtrList2Eigen(modes);
    // To use the SVD modes from the HyperReduction class :
    // Eigen::VectorXd normalizingWeights;
    // c.getModesSVD(c.snapshotsListTuple, snapshotsModes, normalizingWeights);
    // Create DEIM object with given number of basis functions
    Eigen::VectorXi initSeeds(0);
    DEIM_function c(NDEIM, NDEIM, initSeeds, "DEIM", Sp);
    // Compute the DEIM
    c.offlineGappyDEIM(snapshotsModes, normalizingWeights);
    // To access the same folder hierarchy as before :
    // c.problemName = "Gaussian_function";
    // c.offlineGappyDEIM(snapshotsModes, normalizingWeights, "ITHACAoutput/DEIM/"+c.problemName);
    // Generate the submeshes with the depth of the layer
    c.generateSubmesh(2, mesh);
    c.subField = c.interpolateField<volScalarField>(Sp[0]);
    // Define a new online parameter
    Eigen::MatrixXd par_new(2, 1);
    par_new(0, 0) = 0;
    par_new(1, 0) = 0;
    // Online evaluation of the non linear function
    Eigen::VectorXd aprfield = c.MatrixOnline * c.onlineCoeffs(par_new);
    // Transform to an OpenFOAM field and export
    volScalarField S2("S_online", Foam2Eigen::Eigen2field(S, aprfield));
    ITHACAstream::exportSolution(S2, name(1), "./ITHACAoutput/Online/");
    // Evaluate the full order function and export it
    DEIM_function::evaluate_expression(S, par_new);
    ITHACAstream::exportSolution(S, name(1), "./ITHACAoutput/Online/");
    // Compute the L2 error and print it
    Info << ITHACAutilities::errorL2Rel(S2, S) << endl;
    return 0;
}