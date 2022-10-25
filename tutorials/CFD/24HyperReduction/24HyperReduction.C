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
    Example of the reconstruction of a non-linear function using the HyperReduction
SourceFiles
    24HyperReduction.C
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

class HyperReduction_function : public HyperReduction<PtrList<volScalarField>&>
{
    public:
    using HyperReduction::HyperReduction;
        static volScalarField evaluate_expression(volScalarField& S, Eigen::MatrixXd mu)
        {
            volScalarField yPos = S.mesh().C().component(vector::Y).ref();
            volScalarField xPos = S.mesh().C().component(vector::X).ref();

            for (auto i = 0; i < S.size(); i++)
            {
                S[i] = std::exp( - 2 * std::pow(xPos[i] - mu(0) - 1,
                                                2) - 2 * std::pow(yPos[i] - mu(1) - 0.5, 2));
            }

            return S;
        }

        static Eigen::VectorXd evaluate_expression(volScalarField& S, Eigen::MatrixXd mu, List<label> nodesList)
        {
            Eigen::VectorXd ret(nodesList.size());
            volScalarField yPos = S.mesh().C().component(vector::Y).ref();
            volScalarField xPos = S.mesh().C().component(vector::X).ref();

            for (auto i = 0; i < nodesList.size(); i++)
            {
                ret[i] = std::exp( - 2 * std::pow(xPos[nodesList[i]] - mu(0) - 1,
                                                2) - 2 * std::pow(yPos[nodesList[i]] - mu(1) - 0.5, 2));
            }

            return ret;
        }
        
        Eigen::VectorXd onlineCoeffs(volScalarField& S, Eigen::MatrixXd mu, volScalarField& SS)
        {
            Eigen::VectorXd f = evaluate_expression(S, mu, localNodePoints);
            return pinvPU*f;
        }

};

class HyperReduction_vectorFunction : public HyperReduction<PtrList<volVectorField>&>
{
    public:
    using HyperReduction::HyperReduction;
        static volVectorField evaluate_expression(volVectorField& S, Eigen::MatrixXd mu)
        {
            volScalarField yPos = S.mesh().C().component(vector::Y).ref();
            volScalarField xPos = S.mesh().C().component(vector::X).ref();

            for (auto i = 0; i < S.size(); i++)
            {
                S[i][0] = std::exp( - 2 * std::pow(xPos[i] - mu(0) - 1,
                                                2) - 2 * std::pow(yPos[i] - mu(1) - 0.5, 2));
                S[i][1] = std::exp( - 2 * std::pow(xPos[i] - mu(0) - 0.5,
                                                2) - 2 * std::pow(yPos[i] - mu(1) - 0.5, 2));
                S[i][2] = std::exp( - 2 * std::pow(xPos[i] - mu(0) - 1.5,
                                                2) - 2 * std::pow(yPos[i] - mu(1) - 0., 2));
            }

            return S;
        }

        static Eigen::VectorXd evaluate_expression(volVectorField& S, Eigen::MatrixXd mu, List<label> nodesList)
        {
            Eigen::VectorXd ret(nodesList.size()*3);
            volScalarField yPos = S.mesh().C().component(vector::Y).ref();
            volScalarField xPos = S.mesh().C().component(vector::X).ref();

            for (auto i = 0; i < nodesList.size(); i++)
            {
                ret(i*3) = std::exp( - 2 * std::pow(xPos[nodesList[i]] - mu(0) - 1,
                                                2) - 2 * std::pow(yPos[nodesList[i]] - mu(1) - 0.5, 2));
                ret(i*3 + 1) = std::exp( - 2 * std::pow(xPos[nodesList[i]] - mu(0) - 0.5,
                                                2) - 2 * std::pow(yPos[nodesList[i]] - mu(1) - 0.5, 2));
                ret(i*3 + 2) = std::exp( - 2 * std::pow(xPos[nodesList[i]] - mu(0) - 1.5,
                                                2) - 2 * std::pow(yPos[nodesList[i]] - mu(1) - 0., 2));
            }

            return ret;
        }
        
        Eigen::VectorXd onlineCoeffs(volVectorField& S, Eigen::MatrixXd mu)
        {
            Eigen::VectorXd f = evaluate_expression(S, mu, localNodePoints);
            return pinvPU*f;
        }

};

class HyperReduction_vectorScalarFunction : public HyperReduction<PtrList<volVectorField>&, PtrList<volScalarField>&>
{
    public:
    using HyperReduction::HyperReduction;
        static std::pair<volVectorField, volScalarField> evaluate_expression(volVectorField& V, volScalarField& S, Eigen::MatrixXd mu)
        {
            volScalarField yPos = S.mesh().C().component(vector::Y).ref();
            volScalarField xPos = S.mesh().C().component(vector::X).ref();

            for (auto i = 0; i < S.size(); i++)
            {
                V[i][0] = std::exp( - 2 * std::pow(xPos[i] - mu(0) - 1,
                                                2) - 2 * std::pow(yPos[i] - mu(1) - 0.5, 2));
                V[i][1] = std::exp( - 2 * std::pow(xPos[i] - mu(0) - 0.5,
                                                2) - 2 * std::pow(yPos[i] - mu(1) - 0.5, 2));
                V[i][2] = std::exp( - 2 * std::pow(xPos[i] - mu(0) - 1.5,
                                                2) - 2 * std::pow(yPos[i] - mu(1) - 0., 2));
                S[i] = std::exp( - 2 * std::pow(xPos[i] - mu(0) - 1,
                                                2) - 2 * std::pow(yPos[i] - mu(1) - 0.5, 2));
            }

            return std::make_pair(V, S);
        }

        static Eigen::VectorXd evaluate_expression(volScalarField& S, Eigen::MatrixXd mu, List<label> nodesList)
        {
            Eigen::VectorXd ret(nodesList.size()*4);
            volScalarField yPos = S.mesh().C().component(vector::Y).ref();
            volScalarField xPos = S.mesh().C().component(vector::X).ref();

            for (auto i = 0; i < nodesList.size(); i++)
            {
                ret(i*4) = std::exp( - 2 * std::pow(xPos[nodesList[i]] - mu(0) - 1,
                                                2) - 2 * std::pow(yPos[nodesList[i]] - mu(1) - 0.5, 2));
                ret(i*4 + 1) = std::exp( - 2 * std::pow(xPos[nodesList[i]] - mu(0) - 0.5,
                                                2) - 2 * std::pow(yPos[nodesList[i]] - mu(1) - 0.5, 2));
                ret(i*4 + 2) = std::exp( - 2 * std::pow(xPos[nodesList[i]] - mu(0) - 1.5,
                                                2) - 2 * std::pow(yPos[nodesList[i]] - mu(1) - 0., 2));
                ret(i*4 + 3) = std::exp( - 2 * std::pow(xPos[nodesList[i]] - mu(0) - 1,
                                                2) - 2 * std::pow(yPos[nodesList[i]] - mu(1) - 0.5, 2));
            }

            return ret;
        }
        
        Eigen::VectorXd onlineCoeffs(volScalarField& S, Eigen::MatrixXd mu)
        {
            Eigen::VectorXd f = evaluate_expression(S, mu, localNodePoints);
            return pinvPU*f;
        }

};

void test_scalar(ITHACAparameters* para, Foam::fvMesh& mesh, Foam::Time& runTime){
    int n_modes = para->ITHACAdict->lookupOrDefault<int>("Modes", 15);
    int n_nodes = para->ITHACAdict->lookupOrDefault<int>("Nodes", 15);
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
    Eigen::MatrixXd pars;
    cnpy::load(pars, "trainingPars.npy");

    // Perform the offline phase
    for (int i = 0; i < 100; i++)
    {
        HyperReduction_function::evaluate_expression(S, pars.row(i));
        Sp.append((S).clone());
        ITHACAstream::exportSolution(S, name(i + 1), "./ITHACAoutput/Offline/");
    }

    // Create HyperReduction object with given number of basis functions
    Eigen::VectorXi initSeeds(0);

    HyperReduction_function c(n_modes, n_nodes, initSeeds, "Gaussian_function", Sp);
    c.offlineStage();

    // Generate the submeshes with the depth of the layer
    c.generateSubmesh(2, mesh);
    auto sfield = c.interpolateField<volScalarField>(Sp[0]);
    // Define a new online parameter
    Eigen::MatrixXd par_new(2, 1);
    par_new(0, 0) = 0;
    par_new(1, 0) = 0;
    // Online evaluation of the non linear function
    Eigen::VectorXd aprfield = c.renormalizedBasisMatrix*c.onlineCoeffs(sfield(), par_new, Sp[0]);
    // Transform to an OpenFOAM field and export
    volScalarField S2("S_online", Foam2Eigen::Eigen2field(S, aprfield));
    ITHACAstream::exportSolution(S2, name(1), "./ITHACAoutput/Online/");
    // Evaluate the full order function and export it
    HyperReduction_function::evaluate_expression(S, par_new);
    ITHACAstream::exportSolution(S, name(1), "./ITHACAoutput/Online/");
    // Compute the L2 error and print it
    Info << ITHACAutilities::errorL2Rel(S2, S) << endl;
}

void test_vector(ITHACAparameters* para, Foam::fvMesh& mesh, Foam::Time& runTime){
    int n_modes = para->ITHACAdict->lookupOrDefault<int>("Modes", 15);
    int n_nodes = para->ITHACAdict->lookupOrDefault<int>("Nodes", 15);
    simpleControl simple(mesh);

#include "createFields.H"
    // List of volVectorField where the snapshots are stored
    PtrList<volVectorField> Sp;
    // Non linear function
    volVectorField S
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
        dimensionedVector("zero", dimensionSet(0, 0, -1, 1, 0, 0, 0), vector(0, 0, 0))
    );
    // Parameters used to train the non-linear function
    Eigen::MatrixXd pars;
    cnpy::load(pars, "trainingPars.npy");

    // Perform the offline phase
    for (int i = 0; i < 100; i++)
    {
        HyperReduction_vectorFunction::evaluate_expression(S, pars.row(i));
        Sp.append((S).clone());
        ITHACAstream::exportSolution(S, name(i + 1), "./ITHACAoutput/Offline/");
    }

    // Create HyperReduction object with given number of basis functions
    Eigen::VectorXi initSeeds(0);

    HyperReduction_vectorFunction c(n_modes, n_nodes, initSeeds, "Gaussian_function", Sp);
    c.offlineStage();

    // Generate the submeshes with the depth of the layer
    c.generateSubmesh(2, mesh);
    auto sfield = c.interpolateField<volVectorField>(Sp[0]);
    // Define a new online parameter
    Eigen::MatrixXd par_new(2, 1);
    par_new(0, 0) = 0;
    par_new(1, 0) = 0;
    // Online evaluation of the non linear function
    Eigen::VectorXd aprfield = c.renormalizedBasisMatrix*c.onlineCoeffs(sfield(), par_new);
    // Transform to an OpenFOAM field and export
    volVectorField S2("S_online", Foam2Eigen::Eigen2field(S, aprfield));
    ITHACAstream::exportSolution(S2, name(1), "./ITHACAoutput/Online/");
    // Evaluate the full order function and export it
    HyperReduction_vectorFunction::evaluate_expression(S, par_new);
    ITHACAstream::exportSolution(S, name(1), "./ITHACAoutput/Online/");
    // Compute the L2 error and print it
    Info << ITHACAutilities::errorL2Rel(S2, S) << endl;
}

void test_vector_scalar(ITHACAparameters* para, Foam::fvMesh& mesh, Foam::Time& runTime){
    int n_modes = para->ITHACAdict->lookupOrDefault<int>("Modes", 15);
    int n_nodes = para->ITHACAdict->lookupOrDefault<int>("Nodes", 15);
    simpleControl simple(mesh);

#include "createFields.H"
    // List of volVectorField where the snapshots are stored
    PtrList<volVectorField> Vp;
    PtrList<volScalarField> Sp;
    // Non linear function
    volVectorField V
    (
        IOobject
        (
            "V",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T.mesh(),
        dimensionedVector("zero", dimensionSet(0, 0, -1, 1, 0, 0, 0), vector(0, 0, 0))
    );
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
    Eigen::MatrixXd pars;
    cnpy::load(pars, "trainingPars.npy");

    // Perform the offline phase
    for (int i = 0; i < 100; i++)
    {
        HyperReduction_vectorScalarFunction::evaluate_expression(V, S, pars.row(i));
        Vp.append((V).clone());
        Sp.append((S).clone());
        ITHACAstream::exportSolution(S, name(i + 1), "./ITHACAoutput/Offline/");
        ITHACAstream::exportSolution(V, name(i + 1), "./ITHACAoutput/Offline/");
    }

    // Create HyperReduction object with given number of basis functions
    Eigen::VectorXi initSeeds(0);

    HyperReduction_vectorScalarFunction c(n_modes, n_nodes, initSeeds, "Gaussian_function", Vp, Sp);
    c.offlineStage();

    // Generate the submeshes with the depth of the layer
    c.generateSubmesh(2, mesh);
    auto vfield = c.interpolateField<volVectorField>(Vp[0]);
    auto sfield = c.interpolateField<volScalarField>(Sp[0]);
    // Define a new online parameter
    Eigen::MatrixXd par_new(2, 1);
    par_new(0, 0) = 0;
    par_new(1, 0) = 0;
    // Online evaluation of the non linear function
    auto theta = c.onlineCoeffs(sfield(), par_new);
    Eigen::VectorXd aprfield = c.renormalizedBasisMatrix*theta;
    // Transform to an OpenFOAM field and export
    // c.eigen2fields(aprfield, V, S);
    // volVectorField V2("V_online", V);
    // volScalarField S2("S_online", S);
    Eigen::VectorXd recvec = aprfield.head(3*S.size());
    Eigen::VectorXd recsca = aprfield.tail(S.size());
    volVectorField V2("V_online", Foam2Eigen::Eigen2field(V, recvec));
    volScalarField S2("S_online", Foam2Eigen::Eigen2field(S, recsca));

    ITHACAstream::exportSolution(V2, name(1), "./ITHACAoutput/Online/");
    ITHACAstream::exportSolution(S2, name(1), "./ITHACAoutput/Online/");
    // Evaluate the full order function and export it
    HyperReduction_vectorScalarFunction::evaluate_expression(V, S, par_new);
    ITHACAstream::exportSolution(V, name(1), "./ITHACAoutput/Online/");
    ITHACAstream::exportSolution(S, name(1), "./ITHACAoutput/Online/");
    // Compute the L2 error and print it
    Info << "Scalar L2 error: " << ITHACAutilities::errorL2Rel(S2, S) << endl;
    Info << "Vector L2 error: " << ITHACAutilities::errorL2Rel(V2, V) << endl;
}

int main(int argc, char* argv[])
{
    // load stage to perform
    argList::addOption("test", "scalar", "Perform scalar test");
    argList::addOption("test", "vector", "Perform vector test");
    argList::addOption("test", "vs", "Perform (vector, scalar) test");
    // add options for parallel run
    HashTable<string> validParOptions;
    validParOptions.set
    (
        "test",
        "scalar"
    );
    validParOptions.set
    (
        "test",
        "vector"
    );
    validParOptions.set
    (
        "test",
        "vs"
    );

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

    if (args.get("test").match("scalar"))
    {   
        Info << "Init scalar testcase\n";
        test_scalar(para, mesh, runTime);
    }
    else if (args.get("test").match("vector"))
    {
        Info << "Init vector testcase\n";
        test_vector(para, mesh, runTime);
    }
    else if (args.get("test").match("vs"))
    {
        Info << "Init vector scalar testcase\n";
        test_vector_scalar(para, mesh, runTime);
    }
    else
    {
        Info << "Pass '-test scalar', '-test vector', '-test vs'" << endl;
    }
    
    return 0;
}