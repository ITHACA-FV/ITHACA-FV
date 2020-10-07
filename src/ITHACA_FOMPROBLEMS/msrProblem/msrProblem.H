#ifndef msrProblem_H
#define msrProblem_H
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "reductionProblem.H"
#include "ITHACAstream.H"
#include "ITHACAforces.H"
#include "volFields.H"
#include <iostream>
#include "IOmanip.H"
#include "IFstream.H"
#include "primitiveFields.H"
#include "FieldFields.H"
#include "scalarMatrices.H"
#include "SortableList.H"
#include "volFieldsFwd.H"
#include "forces.H"
#include "forceCoeffs.H"
#include "volFields.H"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdio.h>
#include "ITHACAPOD.H"
#include <math.h>



/// Class to implement Molten Salt Reactor multiphysics problem
class msrProblem: public reductionProblem
{

    public:
        //Constructors
        msrProblem();
        ///Construct with argc e argv
        msrProblem(int argc, char* argv[]);
        ~msrProblem() {};

        // Dummy variables to transform msrFoam into a class, createTime.H

        autoPtr<Time> _runTime;

        // Dummy variables to transform msrFoam into a class, createMesh.H

        mutable autoPtr<fvMesh> _mesh;

        // Dummy variables to transform msrFoam into a class, createFields.H

        autoPtr<volVectorField> _U;
        autoPtr<surfaceScalarField> _phi;
        autoPtr<singlePhaseTransportModel> _laminarTransport;
        autoPtr<incompressible::turbulenceModel> turbulence;
        autoPtr<IOMRFZoneList> _MRF;
        autoPtr<volScalarField> _p;
        label pRefCell;
        scalar pRefValue;
        autoPtr<dimensionedScalar> _nu;
        autoPtr<dimensionedScalar> _Pr;
        autoPtr<dimensionedScalar> _Prt;
        autoPtr<dimensionedScalar> _rhoRef;
        autoPtr<volScalarField> _v;
        autoPtr<dimensionedScalar> _CpRef;
        autoPtr<dimensionedScalar> _Tref;
        autoPtr<dimensionedScalar> _TrefXS;
        autoPtr<volScalarField> _TXS;
        /// Initial fields (for restart purposes)
        autoPtr<volScalarField> _p0;
        autoPtr<volVectorField> _U0;
        autoPtr<surfaceScalarField> _phi0;
        autoPtr<volScalarField> _flux0;
        autoPtr<volScalarField> _prec10;
        autoPtr<volScalarField> _prec20;
        autoPtr<volScalarField> _prec30;
        autoPtr<volScalarField> _prec40;
        autoPtr<volScalarField> _prec50;
        autoPtr<volScalarField> _prec60;
        autoPtr<volScalarField> _prec70;
        autoPtr<volScalarField> _prec80;
        autoPtr<volScalarField> _T0;
        autoPtr<volScalarField> _dec10;
        autoPtr<volScalarField> _dec20;
        autoPtr<volScalarField> _dec30;
        autoPtr<dimensionedScalar> _K0;
        autoPtr<volScalarField> _v0;
        autoPtr<volScalarField> _NSF0;
        autoPtr<volScalarField> _A0;
        autoPtr<volScalarField> _D0;
        autoPtr<volScalarField> _SP0;
        autoPtr<volScalarField> _TXS0;


        // Dummy variables to transform msrFoam into a class, createFields_neutronics.H

        autoPtr<dimensionedScalar> _Keff;
        autoPtr<dimensionedScalar> _IV1;
        autoPtr<dimensionedScalar> _D1_0;
        autoPtr<volScalarField> _D;
        autoPtr<dimensionedScalar> _alfa_D1;
        autoPtr<dimensionedScalar> _A1_0;
        autoPtr<volScalarField> _A;
        autoPtr<dimensionedScalar> _alfa_A1;
        autoPtr<dimensionedScalar> _NSF1_0;
        autoPtr<volScalarField> _NSF;
        autoPtr<dimensionedScalar> _alfa_NSF1;
        autoPtr<dimensionedScalar> _SP1_0;
        autoPtr<volScalarField> _SP;
        autoPtr<dimensionedScalar> _alfa_SP1;
        autoPtr<dimensionedScalar> _lam1;
        autoPtr<dimensionedScalar> _lam2;
        autoPtr<dimensionedScalar> _lam3;
        autoPtr<dimensionedScalar> _lam4;
        autoPtr<dimensionedScalar> _lam5;
        autoPtr<dimensionedScalar> _lam6;
        autoPtr<dimensionedScalar> _lam7;
        autoPtr<dimensionedScalar> _lam8;
        autoPtr<dimensionedScalar> _beta1;
        autoPtr<dimensionedScalar> _beta2;
        autoPtr<dimensionedScalar> _beta3;
        autoPtr<dimensionedScalar> _beta4;
        autoPtr<dimensionedScalar> _beta5;
        autoPtr<dimensionedScalar> _beta6;
        autoPtr<dimensionedScalar> _beta7;
        autoPtr<dimensionedScalar> _beta8;
        autoPtr<dimensionedScalar> _betaTot;
        autoPtr<volScalarField> _flux;
        autoPtr<dimensionedScalar> _decLam1;
        autoPtr<dimensionedScalar> _decLam2;
        autoPtr<dimensionedScalar> _decLam3;
        autoPtr<dimensionedScalar> _decBeta1;
        autoPtr<dimensionedScalar> _decBeta2;
        autoPtr<dimensionedScalar> _decBeta3;
        autoPtr<dimensionedScalar> _decbetaTot;
        autoPtr<dimensionedScalar> _betaTE;
        autoPtr<volScalarField> _prec1;
        autoPtr<volScalarField> _prec2;
        autoPtr<volScalarField> _prec3;
        autoPtr<volScalarField> _prec4;
        autoPtr<volScalarField> _prec5;
        autoPtr<volScalarField> _prec6;
        autoPtr<volScalarField> _prec7;
        autoPtr<volScalarField> _prec8;
        autoPtr<dimensionedScalar> _Sc;
        autoPtr<dimensionedScalar> _Sct;
        autoPtr<volScalarField> _alphat;
        autoPtr<volScalarField> _difft;

        //recirculation time
        double tau;


        // Dummy variables to transform msrFoam into a class, createFields_T.H

        autoPtr<volScalarField> _T;
        autoPtr<volScalarField> _logT;
        autoPtr<volScalarField> _dec1;
        autoPtr<volScalarField> _dec2;
        autoPtr<volScalarField> _dec3;

        // Dummy variables to transform msrFoam into a class, createFvOptions.H

        /// fvOptions
        autoPtr<fv::options> _fvOptions;

        // Dummy variables to transform msrFoam into a class, remaining in .C

        /// simpleControl
        autoPtr<simpleControl> _simple;
        /// Tolerance for the residual of the stationary problems, there is the same tolerance for velocity and pressure
        scalar tolerance;
        /// Number of maximum iterations to be done for the computation of the truth solution
        scalar maxIter;
        /// continuity error
        scalar cumulativeContErr = 0;
        label counter = 1;

        /// boolean variable to decide if apply prec inlet BC
        bool precInBool = false;
        /// indexes of inlet and outlet to adopt for precursors boundary conditions
        int precinIndex;
        int precoutIndex;
        /// matrix to store the values of precursors BC inlet conditions
        Eigen::MatrixXd bc_prec;




        // Fields containing the offline solution i.e. the snapshot matrix

        /// List of pointers used to form the pressure snapshots matrix
        PtrList<volScalarField> Pfield;

        /// List of pointers used to form the velocity snapshots matrix
        PtrList<volVectorField> Ufield;

        /// List of pointers used to form the flux snapshots matrix
        PtrList<volScalarField> Fluxfield;

        /// List of pointers used to form the prec1 snapshots matrix
        PtrList<volScalarField> Prec1field;

        /// List of pointers used to form the prec2 snapshots matrix
        PtrList<volScalarField> Prec2field;

        /// List of pointers used to form the prec3 snapshots matrix
        PtrList<volScalarField> Prec3field;

        /// List of pointers used to form the prec4 snapshots matrix
        PtrList<volScalarField> Prec4field;

        /// List of pointers used to form the prec5 snapshots matrix
        PtrList<volScalarField> Prec5field;

        /// List of pointers used to form the prec6 snapshots matrix
        PtrList<volScalarField> Prec6field;

        /// List of pointers used to form the prec7 snapshots matrix
        PtrList<volScalarField> Prec7field;

        /// List of pointers used to form the prec8 snapshots matrix
        PtrList<volScalarField> Prec8field;

        /// List of pointers used to form the temperature snapshots matrix
        PtrList<volScalarField> Tfield;

        /// List of pointers used to form the dec1 snapshots matrix
        PtrList<volScalarField> Dec1field;

        /// List of pointers used to form the dec2 snapshots matrix
        PtrList<volScalarField> Dec2field;

        /// List of pointers used to form the dec3 snapshots matrix
        PtrList<volScalarField> Dec3field;

        /// List of pointers used to form the powerDens snapshots matrix
        PtrList<volScalarField> PowerDensfield;

        /// List of pointers used to form the v snapshosts matrix
        PtrList<volScalarField> vFields;

        /// List of pointers used to form the D snapshosts matrix
        PtrList<volScalarField> DFields;

        /// List of pointers used to form the NSF snapshosts matrix
        PtrList<volScalarField> NSFFields;

        /// List of pointers used to form the A snapshosts matrix
        PtrList<volScalarField> AFields;

        /// List of pointers used to form the SP snapshosts matrix
        PtrList<volScalarField> SPFields;

        /// List of pointers used to form the SP snapshosts matrix
        PtrList<volScalarField> TXSFields;

        /// List of pointers used to form the pressure modes
        PtrList<volScalarField> Pmodes;

        /// List of pointers used to form the velocity modes
        PtrList<volVectorField> Umodes;

        /// List of pointers used to form the flux modes
        PtrList<volScalarField> Fluxmodes;

        /// List of pointers used to form the prec1 modes
        PtrList<volScalarField> Prec1modes;

        /// List of pointers used to form the prec1 modes
        PtrList<volScalarField> Prec2modes;

        /// List of pointers used to form the prec1 modes
        PtrList<volScalarField> Prec3modes;

        /// List of pointers used to form the prec1 modes
        PtrList<volScalarField> Prec4modes;

        /// List of pointers used to form the prec1 modes
        PtrList<volScalarField> Prec5modes;

        /// List of pointers used to form the prec1 modes
        PtrList<volScalarField> Prec6modes;

        /// List of pointers used to form the prec1 modes
        PtrList<volScalarField> Prec7modes;

        /// List of pointers used to form the prec1 modes
        PtrList<volScalarField> Prec8modes;

        /// List of pointers used to form the temperature modes
        PtrList<volScalarField> Tmodes;

        /// List of pointers used to form the dec1 modes
        PtrList<volScalarField> Dec1modes;

        /// List of pointers used to form the dec2 modes
        PtrList<volScalarField> Dec2modes;

        /// List of pointers used to form the dec3 modes
        PtrList<volScalarField> Dec3modes;

        /// List of pointers used to form the v modes
        PtrList<volScalarField> vmodes;

        /// List of pointers used to form the D modes
        PtrList<volScalarField> Dmodes;

        /// List of pointers used to form the NSF snapshosts matrix
        PtrList<volScalarField> NSFmodes;

        /// List of pointers used to form the A snapshosts matrix
        PtrList<volScalarField> Amodes;

        /// List of pointers used to form the SP snapshosts matrix
        PtrList<volScalarField> SPmodes;

        /// List of pointers used to form the SP snapshosts matrix
        PtrList<volScalarField> TXSmodes;

        //----------------------------------------------------------------------//
        // Coefficient matrices:

        //Fluid-dynamics matrices:

        /// Diffusion term
        Eigen::MatrixXd B_matrix;

        /// Mass Matrix
        Eigen::MatrixXd M_matrix;

        /// Gradient of pressure matrix
        Eigen::MatrixXd K_matrix;

        /// Non linear term
        List <Eigen::MatrixXd> C_matrix;

        /// Div of velocity
        Eigen::MatrixXd P_matrix;

        /// Laplacian term PPE
        Eigen::MatrixXd D_matrix;

        /// Divergence of momentum PPE
        List <Eigen::MatrixXd> G_matrix;

        /// PPE BC1
        Eigen::MatrixXd BC1_matrix;

        /// PPE BC2
        List <Eigen::MatrixXd> BC2_matrix;

        /// PPE BC3
        Eigen::MatrixXd BC3_matrix;


        //Neutronics matrices:

        /// laplacian_flux
        List<Eigen::MatrixXd> LF_matrix;
        /// mass flux
        Eigen::MatrixXd MF_matrix;
        /// production flux
        List<Eigen::MatrixXd> PF_matrix;
        /// absorption flux
        List<Eigen::MatrixXd> AF_matrix;

        ITHACAparameters* para;

        /// prec_source 1
        Eigen::MatrixXd PS1_matrix;
        /// prec_source 2
        Eigen::MatrixXd PS2_matrix;
        /// prec_source 3
        Eigen::MatrixXd PS3_matrix;
        /// prec_source 4
        Eigen::MatrixXd PS4_matrix;
        /// prec_source 5
        Eigen::MatrixXd PS5_matrix;
        /// prec_source 6
        Eigen::MatrixXd PS6_matrix;
        /// prec_source 7
        Eigen::MatrixXd PS7_matrix;
        /// prec_source 8
        Eigen::MatrixXd PS8_matrix;

        /// precursor stream term-1
        List<Eigen::MatrixXd> ST1_matrix;
        /// precursor stream term-2
        List<Eigen::MatrixXd> ST2_matrix;
        /// precursor stream term-3
        List<Eigen::MatrixXd> ST3_matrix;
        /// precursor stream term-4
        List<Eigen::MatrixXd> ST4_matrix;
        /// precursor stream term-5
        List<Eigen::MatrixXd> ST5_matrix;
        /// precursor stream term-6
        List<Eigen::MatrixXd> ST6_matrix;
        /// precursor stream term-7
        List<Eigen::MatrixXd> ST7_matrix;
        /// precursor stream term-8
        List<Eigen::MatrixXd> ST8_matrix;
        /// precursor mass term-1
        Eigen::MatrixXd MP1_matrix;
        /// precursor mass term-2
        Eigen::MatrixXd MP2_matrix;
        /// precursor mass term-3
        Eigen::MatrixXd MP3_matrix;
        /// precursor mass term-4
        Eigen::MatrixXd MP4_matrix;
        /// precursor mass term-5
        Eigen::MatrixXd MP5_matrix;
        /// precursor mass term-6
        Eigen::MatrixXd MP6_matrix;
        /// precursor mass term-7
        Eigen::MatrixXd MP7_matrix;
        /// precursor mass term-8
        Eigen::MatrixXd MP8_matrix;
        /// precursor laplacian term-1
        Eigen::MatrixXd LP1_matrix;
        /// precursor laplacian term-2
        Eigen::MatrixXd LP2_matrix;
        /// precursor laplacian term-3
        Eigen::MatrixXd LP3_matrix;
        /// precursor laplacian term-4
        Eigen::MatrixXd LP4_matrix;
        /// precursor laplacian term-5
        Eigen::MatrixXd LP5_matrix;
        /// precursor laplacian term-6
        Eigen::MatrixXd LP6_matrix;
        /// precursor laplacian term-7
        Eigen::MatrixXd LP7_matrix;
        /// precursor laplacian term-8
        Eigen::MatrixXd LP8_matrix;
        /// precursor flux source term-1
        List<Eigen::MatrixXd> FS1_matrix;
        /// precursor flux source term-2
        List<Eigen::MatrixXd> FS2_matrix;
        /// precursor flux source term-3
        List<Eigen::MatrixXd> FS3_matrix;
        /// precursor flux source term-4
        List<Eigen::MatrixXd> FS4_matrix;
        /// precursor flux source term-5
        List<Eigen::MatrixXd> FS5_matrix;
        /// precursor flux source term-6
        List<Eigen::MatrixXd> FS6_matrix;
        /// precursor flux source term-7
        List<Eigen::MatrixXd> FS7_matrix;
        /// precursor flux source term-8
        List<Eigen::MatrixXd> FS8_matrix;


        //Energy matrices:

        /// decay heat stream term-1
        List<Eigen::MatrixXd> SD1_matrix;
        /// decay heat stream term-2
        List<Eigen::MatrixXd> SD2_matrix;
        /// decay heat stream term-3
        List<Eigen::MatrixXd> SD3_matrix;
        /// decay heat mass term-1
        Eigen::MatrixXd MD1_matrix;
        /// decay heat mass term-2
        Eigen::MatrixXd MD2_matrix;
        /// decay heat mass term-3
        Eigen::MatrixXd MD3_matrix;
        /// decay heat laplacian term-1
        Eigen::MatrixXd LD1_matrix;
        /// decay heat laplacian term-2
        Eigen::MatrixXd LD2_matrix;
        /// decay heat laplacian term-3
        Eigen::MatrixXd LD3_matrix;
        /// decay heat flux source term-1
        List<Eigen::MatrixXd> DFS1_matrix;
        /// decay heat flux source term-2
        List<Eigen::MatrixXd> DFS2_matrix;
        /// decay heat flux source term-3
        List<Eigen::MatrixXd> DFS3_matrix;

        ///temperature mass term
        Eigen::MatrixXd TM_matrix;
        /// temperature stream term
        List<Eigen::MatrixXd> TS_matrix;
        /// temperature laplacian term
        Eigen::MatrixXd LT_matrix;
        /// temperature flux source term TXS
        List<Eigen::MatrixXd> TXS_matrix;
        /// temperature decay heat source term-1
        List<Eigen::MatrixXd> THS1_matrix;
        /// temperature decay heat source term-2
        List<Eigen::MatrixXd> THS2_matrix;
        /// temperature decay heat source term-3
        List<Eigen::MatrixXd> THS3_matrix;
        //--------------------------------------------------------------------------------//

        /// Number of modes adopted during Galerkin projection

        label NUmodes;
        label NPmodes;
        label NFluxmodes;
        Eigen::VectorXi NPrecmodes;
        Eigen::VectorXi NDecmodes;
        label NTmodes;
        label NCmodes;


        //--------------------------------------------------------------------------------//
        //Lift-functions and homogenization related terms:

        /// List of pointers used to form the list of lifting functions
        PtrList<volVectorField> liftfield;
        /// List of pointers used to form the list of lifting functions
        PtrList<volScalarField> liftfieldT;
        /// List of pointers used to form the homogeneous velocity snapshots
        PtrList<volVectorField> Uomfield;
        /// List of pointers used to form the homogeneous velocity snapshots
        PtrList<volScalarField> Tomfield;
        /// boolean variables to check if the homogenization of U and T is
        /// performed (true) or not (false)
        bool homboolU = false;
        bool homboolT = false;

        //-----------------------------------------------------------------------------//
        //RBF interpolation procedure quantities, i.e. constants changing with
        //temperature

        // Create a SAMPLES for interpolation
        std::vector<SPLINTER::DataTable*> SAMPLES_v;
        // Create a SAMPLES for interpolation
        std::vector<SPLINTER::DataTable*> SAMPLES_D;
        // Create a SAMPLES for interpolation
        std::vector<SPLINTER::DataTable*> SAMPLES_NSF;
        // Create a SAMPLES for interpolation
        std::vector<SPLINTER::DataTable*> SAMPLES_A;
        // Create a SAMPLES for interpolation
        std::vector<SPLINTER::DataTable*> SAMPLES_SP;
        // Create a SAMPLES for interpolation
        std::vector<SPLINTER::DataTable*> SAMPLES_TXS;
        // Create a RBF splines for interpolation
        std::vector<SPLINTER::RBFSpline*> rbfsplines_v;
        // Create a RBF splines for interpolation
        std::vector<SPLINTER::RBFSpline*> rbfsplines_D;
        // Create a RBF splines for interpolation
        std::vector<SPLINTER::RBFSpline*> rbfsplines_NSF;
        // Create a RBF splines for interpolation
        std::vector<SPLINTER::RBFSpline*> rbfsplines_A;
        // Create a RBF splines for interpolation
        std::vector<SPLINTER::RBFSpline*> rbfsplines_SP;
        // Create a RBF splines for interpolation
        std::vector<SPLINTER::RBFSpline*> rbfsplines_TXS;

        //-----------------------------------------------------------------------------------------//
        // Methods:

        /// Perform the offline solution of MSR steady state problem
        void truthSolve(List<scalar> mu_now);
        //----------------------------------------------------------------------------------

        /// method to apply RBF interpolation procedure
        /// NC is the number of modes to adopt for construncting the basis
        void msrcoeff(label& NC);
        //----------------------------------------------------------------------------------

        /// Perform a lift solve for the velocity field
        void liftSolve();
        //----------------------------------------------------------------------------------

        /// Perform a lift solve for the temperature
        void liftSolveT();
        //----------------------------------------------------------------------------------

        /// Method to read all the fieds of the MSR problem in the offline folder,
        /// it also reads mu_samples matrix. First method reads from standard "Offline" folder
        /// second one reads from folder specified by the user, solutions for mu.col(0), mu.col(1), mu.col(2), etc...
        ///are sought in dir/0, dir/1, dir/2, etc...
        void readMSRfields();
        void readMSRfields(std::string& dir);
        //-----------------------------------------------------------------------------------

        /// Method to compute the homogenized velocity field, it also sets homboolU=true
        void homogenizeU();
        //-----------------------------------------------------------------------------------

        /// Method to compute the homogenized temperature field, it also sets homboolT=true
        void homogenizeT();
        //--------------------------------------------------------------------------

        /// Method to compute the modes for all the fields in the MSR problem,
        /// if hombool==false the velocity modes are computed form Uomfield snapshost, Ufields otherwise.
        /// Same for temperature field as regards homboolT.
        /// Number of modes are specified in ITHACAdict file, if not set, 10 is used.
        /// Singular Value Decomposition method
        void msrgetModesSVD();
        //--------------------------------------------------------------------------

        /// Method to compute the modes for all the fields in the MSR problem,
        /// if hombool==false the velocity modes are computed form Uomfield snapshost, Ufields otherwise.
        /// Same for temperature field as regards homboolT.
        /// Number of modes are specified in ITHACAdict file, if not set, 10 is used.
        /// EigenValue Decomposition method
        void msrgetModesEVD();
        //--------------------------------------------------------------------------

        /// Project using the Poisson Equation for pressure
        ///
        /// @param[in]  folder      The folder used to save the reduced matrices.
        /// @param[in]  NUmodes     The number of velocity modes.
        /// @param[in]  NPmodes     The number of pressure modes.
        /// @param[in]  NFluxmodes  The number of flux modes.
        /// @param[in]  Nprecmodes  The number of precursors modes. Nprecmodes.size() must be 8
        /// @param[in]  NFluxmodes  The number of temperature modes.
        /// @param[in]  Nprecmodes  The number of decay heat modes. Ndecmodes.size() must be 3
        /// @param[in]  NFluxmodes  The number of temperature changing constants modes. (see file updateConstants.H)
        void projectPPE(fileName folder, label NUmodes, label NPmodes, label NFluxmodes,
                        Eigen::VectorXi Nprecmodes, label NTmodes, Eigen::VectorXi Ndecmodes,
                        label NCmodes);
        //--------------------------------------------------------------------------

        /// sub-functions needed by projectPPE

        ///momentum eq. methods:
        Eigen::MatrixXd diffusive_term(label NUmodes, label NPmodes);
        Eigen::MatrixXd pressure_gradient_term(label NUmodes, label NPmodes);
        List < Eigen::MatrixXd > convective_term(label NUmodes, label NPmodes);
        Eigen::MatrixXd mass_term(label NUmodes, label NPmodes);

        /// continuity eq. methods:
        Eigen::MatrixXd divergence_term(label NUmodes, label NPmodes);
        List < Eigen::MatrixXd > div_momentum(label NUmodes, label NPmodes);
        Eigen::MatrixXd laplacian_pressure(label NPmodes);
        Eigen::MatrixXd pressure_BC1(label NUmodes, label NPmodes);
        List < Eigen::MatrixXd > pressure_BC2(label NUmodes, label NPmodes);
        Eigen::MatrixXd pressure_BC3(label NUmodes, label NPmodes);

        /// diffusion eq. methods:
        List<Eigen::MatrixXd> laplacian_flux(label NFluxmodes, label NCmodes);
        Eigen::MatrixXd mass_flux(label NFluxmodes);
        List<Eigen::MatrixXd> prod_flux(label NFluxmodes, label NCmodes);
        List<Eigen::MatrixXd> abs_flux(label NFluxmodes, label NCmodes);
        Eigen::MatrixXd prec_source(label NFluxmodes, label NPrecmodes, label family);

        /// precursor eq. methods:
        List<Eigen::MatrixXd> stream_term(label NUmodes, label NPrecmodes,
                                          label family);
        Eigen::MatrixXd prec_mass(label NPrecmodes, label family);
        Eigen::MatrixXd laplacian_prec(label NPrecmodes, label family);
        List<Eigen::MatrixXd> flux_source(label NFluxmodes, label NPrecmodes,
                                          label NCmodes, label family);

        /// decay heat eq. methods:
        List<Eigen::MatrixXd> stream_dec(label NUmodes, label NDecmodes,
                                         label decgroup);
        Eigen::MatrixXd dec_mass(label NDecmodes, label decgroup);
        Eigen::MatrixXd laplacian_dec(label NDecmodes, label decgroup);
        List<Eigen::MatrixXd> dec_fluxsource(label NFluxmodes, label NDecmodes,
                                             label NCmodes, label decgroup);

        /// temperature eq. methods
        Eigen::MatrixXd mass_temp(label NTmodes);
        List<Eigen::MatrixXd> temp_stream(label NUmodes, label NTmodes);
        Eigen::MatrixXd laplacian_temp(label NTmodes);
        List<Eigen::MatrixXd> temp_heatsource(label NTmodes, label NDecmodes,
                                              label NCmodes, label decgroup);
        List<Eigen::MatrixXd> temp_XSfluxsource(label NTmodes, label NFluxmodes,
                                                label NCmodes);

        //--------------------------------------------------------------------------
        /// method to change the viscosity in UEqn
        void change_viscosity(double mu);

        //--------------------------------------------------------------------------
        /// method to set all fields back to values in 0 folder
        void restart();



    protected:

        //--------------------------------------------------------------------------
        /// method to choose one field among precs & decs
        /// field can be "prec" or "dec" only
        /// if field==prec then ith can range from 1 to 8 included
        /// else ith can range from 1 to 3

        PtrList<volScalarField> choose_group(string field, int ith);

        //--------------------------------------------------------------------------
        /// method to save matrices for precs and decs M can be an Eigen::MatrixXd or
        /// List<Eigen::MatrixXd>
        /// nome is the name of the matrix, n can vary from 1 to 8 and it appended to matrix's
        /// name to store the matrix belonging to group n
        /// folder is the folder in which the matrix is stored

        template<typename M>
        void savegroupMatrix(string nome, int n, word folder, M matrice);

        ///note: these methods are protected since they are very peculiar and problem
        ///dependent, in this way the user cannot use it in the main of the programm causing
        ///errors.


};

#endif











