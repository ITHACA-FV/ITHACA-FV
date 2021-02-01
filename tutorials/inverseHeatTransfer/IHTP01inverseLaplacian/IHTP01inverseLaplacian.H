#ifndef IHTP01inverseLaplacian_H
#define IHTP01inverseLaplacian_H

/// \brief Class where the first inverse heat transfer problem
/// tutorial is solved using the Parameterization of the heat flux.
/// \details It is a child of the inverseLaplacianProblem_paramBC
/// class and some of its functions are overridden to be adapted
/// to the specific case.
class IHTP01inverseLaplacian_paramBC: public inverseLaplacianProblem_paramBC
{
    public:
        explicit IHTP01inverseLaplacian_paramBC(int argc, char* argv[])
            :
            inverseLaplacianProblem_paramBC(argc, argv),
            T(_T()),
            mesh(_mesh()),
            runTime(_runTime())
        {
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            thermocouplesRead = 0;
        }

        /// Temperature field
        volScalarField& T;

        /// Mesh
        fvMesh& mesh;

        /// RunTime
        Time& runTime;

        /// Constant to define boundary conditions
        double a;

        /// Constant to define boundary conditions
        double b;

        /// Constant to define boundary conditions
        double c;

        /// Constant to define boundary conditions
        double d;

        /// Boundary heat flux at the gammaEx1 patch
        List<scalar> heatFlux_gammaEx1;

        /// Boundary heat flux at the gammaEx2 patch
        List<scalar> heatFlux_gammaEx2;

        /// Boundary heat flux at the gammaEx3 patch
        List<scalar> heatFlux_gammaEx3;

        /// Boundary heat flux at the gammaEx4 patch
        List<scalar> heatFlux_gammaEx4;

        /// gammaEx1 boundary patches label
        label gammaEx1_ind;

        /// gammaEx2 boundary patches label
        label gammaEx2_ind;

        /// gammaEx3 boundary patches label
        label gammaEx3_ind;

        /// gammaEx4 boundary patches label
        label gammaEx4_ind;

        /// List of estimated heat fluxes used for post processing
        PtrList<volScalarField> gField;


        //--------------------------------------------------------------------------
        /// Solve the additional problem
        ///
        void solveAdditional()
        {
            restart();
            fvMesh& mesh = _mesh();
            simpleControl& simple = _simple();
            volScalarField Tad(_T);
            Foam::Time& runTime = _runTime();
            set_valueFraction();
            List<scalar> RobinBC = - Tf;
            forAll(mesh.boundaryMesh(), patchI)
            {
                if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
                {
                    ITHACAutilities::assignMixedBC(Tad, patchI, RobinBC, refGrad, valueFraction);
                }
                else if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
                {
                    ITHACAutilities::assignBC(Tad, patchI, homogeneousBC);
                }
                else if (patchI == mesh.boundaryMesh().findPatchID("gammaEx1"))
                {
                    ITHACAutilities::assignBC(Tad, patchI, -heatFlux_gammaEx1);
                }
                else if (patchI == mesh.boundaryMesh().findPatchID("gammaEx2"))
                {
                    ITHACAutilities::assignBC(Tad, patchI, -heatFlux_gammaEx2);
                }
                else if (patchI == mesh.boundaryMesh().findPatchID("gammaEx3"))
                {
                    ITHACAutilities::assignBC(Tad, patchI, -heatFlux_gammaEx3);
                }
                else if (patchI == mesh.boundaryMesh().findPatchID("gammaEx4"))
                {
                    ITHACAutilities::assignBC(Tad, patchI, -heatFlux_gammaEx4);
                }
            }
#if OFVER == 6

            while (simple.loop(runTime))
#else
            while (simple.loop())
#endif
            {
                while (simple.correctNonOrthogonal())
                {
                    fvScalarMatrix TEqn
                    (
                        fvm::laplacian(DT, Tad)
                    );
                    TEqn.solve();
                }
            }

            Tad_base.append(Tad.clone());
            //Reinitializing runTime
            instantList Times = runTime.times();
            runTime.setTime(Times[1], 1);
            addSol = fieldValueAtThermocouples(Tad);
            ITHACAstream::exportSolution(Tad, "1",
                                         folderOffline,
                                         "Tad");
        }

        //--------------------------------------------------------------------------
        /// Post process all the estimated heat fluxes heatFluxFieldName inside folder
        /// computing the relative error field with respect to gTrue with the related
        /// L2 and Linfty norms
        ///
        /// @param[in]  folder              Folder where the estimatied heat fluxes are read and the post processed results saved
        /// @param[in]  heatFluxFieldName   Name of the heat flux to post process
        /// @param[in]  innerField          Value of the output field in the interior of the domain (0 by default)
        ///
        void postProcess(word folder, word heatFluxFieldName, scalar innerField = 0.0)
        {
#include"postProcess.H"
        }


        //--------------------------------------------------------------------------
        /// Assign the correct boundary condition to direct problem and assign -g/k to the hotSide patch
        ///
        void assignDirectBC()
        {
#include"directBC.H"
        }

        //--------------------------------------------------------------------------
        /// Perform a solution of the direct problem with the correct boundary conditions
        ///
        void solveTrue()
        {
            word outputFolder = "./ITHACAoutput/true/";
#include"solveTrue.H"
        }

        //--------------------------------------------------------------------------
        /// Compute the best approximation of the true heat flux gTrue in the parameterized space
        ///
        /// @return     Best approximation of the true heat flux
        ///
        Eigen::VectorXd bestApproximator()
        {
            fvMesh& mesh = _mesh();
            int Nbasis = gBaseFunctions.size();
            Eigen::MatrixXd mass(Nbasis, Nbasis);
            Eigen::VectorXd source(Nbasis);
            forAll(gBaseFunctions, baseI)
            {
                forAll(gBaseFunctions, baseJ)
                {
                    mass(baseI, baseJ) = ITHACAutilities::L2productOnPatch(mesh,
                                         gBaseFunctions[baseI], gBaseFunctions[baseJ], "hotSide");
                }
                source(baseI) = ITHACAutilities::L2productOnPatch(mesh, gBaseFunctions[baseI],
                                gTrue, "hotSide");
            }
            return mass.fullPivLu().solve(source); //Vector of coefficients
        }


        //--------------------------------------------------------------------------
        /// Compute the best interpolation of the true heat flux gTrue in the parameterized space
        ///
        /// @return     Best interpolation of the true heat flux
        ///
        Eigen::VectorXd bestInterpolation()
        {
            int Nbasis = gBaseFunctions.size();
            Eigen::VectorXd coeff(Nbasis);
            Eigen::MatrixXd Phi = Eigen::MatrixXd::Zero(Nbasis, Nbasis);
            forAll(gBaseFunctions, baseI)
            {
                scalar centerX = thermocouplesPos[baseI].x();
                scalar centerZ = thermocouplesPos[baseI].z();
                forAll(thermocouplesPos, TCi)
                {
                    scalar tcX = thermocouplesPos[TCi].x();
                    scalar tcZ = thermocouplesPos[TCi].z();
                    double radius = Foam::sqrt((tcX - centerX) * (tcX - centerX) +
                                               (tcZ - centerZ) * (tcZ - centerZ));
                    Phi(TCi, baseI) = Foam::exp(- (shapeParameter *
                                                   shapeParameter
                                                   * radius * radius));
                }
                coeff(baseI) = k * (b * thermocouplesPos[baseI].x() + c);
            }
            return Phi.fullPivLu().solve(coeff);
        }

};



/// \brief Class where the first inverse heat transfer problem
/// tutorial is solved using Alifanov's regularization.
/// \details It is a child of the inverseLaplacianProblem_CG
/// class and some of its functions are overridden to be adapted
/// to the specific case.
class IHTP01inverseLaplacian_CG: public inverseLaplacianProblem_CG
{
    public:
        explicit IHTP01inverseLaplacian_CG(int argc, char* argv[])
            :
            inverseLaplacianProblem_CG(argc, argv),
            T(_T()),
            lambda(_lambda()),
            deltaT(_deltaT()),
            mesh(_mesh()),
            runTime(_runTime())
        {
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            interpolationPlaneDefined = 0;
            cgIter = 0;
            thermocouplesRead = 0;
        }

        /// Temperature field
        volScalarField& T;

        /// Adjoint field
        volScalarField& lambda;

        /// Sensitivity field
        volScalarField& deltaT;

        /// Mesh
        fvMesh& mesh;

        /// RunTime
        Time& runTime;

        /// Constant to define boundary conditions
        double a;

        /// Constant to define boundary conditions
        double b;

        /// Constant to define boundary conditions
        double c;

        /// Constant to define boundary conditions
        double d;

        /// Boundary heat flux at the gammaEx1 patch
        List<scalar> heatFlux_gammaEx1;

        /// Boundary heat flux at the gammaEx2 patch
        List<scalar> heatFlux_gammaEx2;

        /// Boundary heat flux at the gammaEx3 patch
        List<scalar> heatFlux_gammaEx3;

        /// Boundary heat flux at the gammaEx4 patch
        List<scalar> heatFlux_gammaEx4;

        /// gammaEx1 boundary patches label
        label gammaEx1_ind;

        /// gammaEx2 boundary patches label
        label gammaEx2_ind;

        /// gammaEx3 boundary patches label
        label gammaEx3_ind;

        /// gammaEx4 boundary patches label
        label gammaEx4_ind;

        /// List of estimated heat fluxes used for post processing
        PtrList<volScalarField> gField;

        //--------------------------------------------------------------------------
        /// Post process all the estimated heat fluxes heatFluxFieldName inside folder
        /// computing the relative error field with respect to gTrue with the related
        /// L2 and Linfty norms
        ///
        /// @param[in]  folder              Folder where the estimatied heat fluxes are read and the post processed results saved
        /// @param[in]  heatFluxFieldName   Name of the heat flux to post process
        /// @param[in]  innerField          Value of the output field in the interior of the domain (0 by default)
        ///
        void postProcess(word folder, word heatFluxFieldName, scalar innerField = 0.0)
        {
#include"postProcess.H"
        }


        //--------------------------------------------------------------------------
        /// Assign the correct boundary condition to direct problem and assign -g/k to the hotSide patch
        ///
        void assignDirectBC()
        {
#include"directBC.H"
        }

        //--------------------------------------------------------------------------
        /// Perform a solution of the direct problem with the correct boundary conditions
        ///
        void solveTrue()
        {
            word outputFolder = "./ITHACAoutput/true/";
#include"solveTrue.H"
        }
};
#endif
