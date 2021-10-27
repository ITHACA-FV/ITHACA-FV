# Introduction to tutorial 2

The problems consists of a parametrized POD-Galerkin ROM for a steady heat transfer phenomenon. The parametrization is on the diffusivity constant. The OpenFOAM full order problem is based on **laplacianFoam**.
The problem equation are
<center>
<img src="https://render.githubusercontent.com/render/math?math=\nabla \cdot (k \nabla T) = S"> 
</center>
where <img src="https://render.githubusercontent.com/render/math?math=k"> is the diffusivity, <img src="https://render.githubusercontent.com/render/math?math=T"> is the temperature and <img src="https://render.githubusercontent.com/render/math?math=k"> is the diffusivity, <img src="https://render.githubusercontent.com/render/math?math=S"> is the source term. The problem discretised and formalized in matrix equation reads:
<center>
<img src="https://render.githubusercontent.com/render/math?math=AT = S"></center>
where <img src="https://render.githubusercontent.com/render/math?math=A"> is the matrix of interpolation coefficients, <img src="https://render.githubusercontent.com/render/math?math=T"> is the vector of unknowns and <img src="https://render.githubusercontent.com/render/math?math=S"> is the vector representing the source term. The domain is subdivided in 9 different parts and each part has parametrized diffusivity. See the image below for a clarification.
<center>
![drawing](../../../docs/images/drawing.png)
</center>
Both the full order and the reduced order problem are solved exploiting the parametric affine decomposition of the differential operators:
<center>
<img src="https://render.githubusercontent.com/render/math?math=A = \sum_{i=1}^N \theta_i(\mu) A_i"> 
</center>

For the operations performed by each command check the comments in the source [**02thermalBlock.C**](https://mathlab.github.io/ITHACA-FV/02thermalBlock_8C.html) file.

## The necessary header files
First of all let's have a look to the header files that needs to be included and what they are responsible for:

The standard C++ header for input/output stream objects:

    #include <iostream>
The OpenFOAM header files:

    #include "fvCFD.H"
    #include "IOmanip.H"
    #include "Time.H"

The header file of ITHACA-FV necessary for this tutorial

    #include "ITHACAPOD.H"
    #include "ITHACAutilities.H"

The [**Eigen**](https://mathlab.github.io/ITHACA-FV/namespaceEigen.html) library for matrix manipulation and linear and non-linear algebra operations:

    #include <Eigen/Dense>

## Implementation of the tutorial02 class

Then we can define the [**tutorial02**](https://mathlab.github.io/ITHACA-FV/classtutorial02.html) class as a child of the [**laplacianProblem**](https://mathlab.github.io/ITHACA-FV/classlaplacianProblem.html) class

    class tutorial02: public laplacianProblem
    {
        public:
            explicit tutorial02(int argc, char* argv[])
                :
                laplacianProblem(argc, argv),
                T(_T()),
                nu(_nu()),
                S(_S())
            {}

The members of the class are the fields that need to be manipulated during the resolution of the problem

Inside the class it is defined the offline solve method according to the specific parametrized problem that needs to be solved.

    
    void offlineSolve(word folder = "./ITHACAoutput/Offline/")
        {

If the offline solve has already been performed than read the existing snapshots


    if (offline)
    {
        ITHACAstream::read_fields(Tfield, "T", folder);
        mu_samples =
            ITHACAstream::readMatrix(folder + "/mu_samples_mat.txt");
    }

else perform the offline solve where a loop over all the parameters is performed:

    for (label i = 0; i < mu.rows(); i++)
    {
        for (label j = 0; j < mu.cols() ; j++)
        {
            // mu_now[i] =
            theta[j] = mu(i, j);
        }

a 0 internal constant value is assigned before each solve command with the lines

    assignIF(T, IF);

and the solve operation is performed, see also the [**laplacianProblem**](https://mathlab.github.io/ITHACA-FV/classlaplacianProblem.html) class for the definition of the methods

    truthSolve(mu_now, folder);

We need also to implement a method to set/define the source term that may be problem dependent. In this case the source term is defined with a hat function:

    void SetSource()
    {
        volScalarField yPos = T.mesh().C().component(vector::Y);
        volScalarField xPos = T.mesh().C().component(vector::X);
        forAll(S, counter)
        {
            S[counter] = Foam::sin(xPos[counter] / 0.9 * M_PI) + Foam::sin(
                                yPos[counter] / 0.9 * M_PI);
        }
    }
Define by:
<center>
<img src="https://render.githubusercontent.com/render/math?math=S = \sin(\frac{\pi}{L}\cdot x) + \sin(\frac{\pi}{L}\cdot y)"> 
</center>

where <img src="https://render.githubusercontent.com/render/math?math=L"> is the dimension of the thermal block which is equal to 0.9.

<center>
![hat](../../../docs/images/hat.png)
</center>

With the following is defined a method to set compute the parameter of the affine expansion:

    void compute_nu()
    {
The list of parameters is resized according to number of parametrized regions

    nu_list.resize(9);

The nine different volScalarFields to identify the viscosity in each domain are initialized:

    volScalarField nu1(nu);
    volScalarField nu2(nu);
    volScalarField nu3(nu);
    volScalarField nu4(nu);
    volScalarField nu5(nu);
    volScalarField nu6(nu);
    volScalarField nu7(nu);
    volScalarField nu8(nu);
    volScalarField nu9(nu);

and the 9 different boxes are defined:

    Eigen::MatrixXd Box1(2, 3);
    Box1 << 0, 0, 0, 0.3, 0.3, 0.1;
    Eigen::MatrixXd Box2(2, 3);
    Box2 << 0.3, 0, 0, 0.6, 0.3, 0.1;
    Eigen::MatrixXd Box3(2, 3);
    Box3 << 0.6, 0, 0, 0.91, 0.3, 0.1;
    Eigen::MatrixXd Box4(2, 3);
    Box4 << 0, 0.3, 0, 0.3, 0.6, 0.1;
    Eigen::MatrixXd Box5(2, 3);
    Box5 << 0.3, 0.3, 0, 0.6, 0.6, 0.1;
    Eigen::MatrixXd Box6(2, 3);
    Box6 << 0.6, 0.3, 0, 0.91, 0.6, 0.1;
    Eigen::MatrixXd Box7(2, 3);
    Box7 << 0, 0.6, 0, 0.3, 0.91, 0.1;
    Eigen::MatrixXd Box8(2, 3);
    Box8 << 0.3, 0.61, 0, 0.6, 0.91, 0.1;
    Eigen::MatrixXd Box9(2, 3);

and for each of the defined boxes the relative diffusivity field is set to 1 inside the box and remain 0 elsewhere:

    ITHACAutilities::setBoxToValue(nu1, Box1, 1.0);
    ITHACAutilities::setBoxToValue(nu2, Box2, 1.0);
    ITHACAutilities::setBoxToValue(nu3, Box3, 1.0);
    ITHACAutilities::setBoxToValue(nu4, Box4, 1.0);
    ITHACAutilities::setBoxToValue(nu5, Box5, 1.0);
    ITHACAutilities::setBoxToValue(nu6, Box6, 1.0);
    ITHACAutilities::setBoxToValue(nu7, Box7, 1.0);
    ITHACAutilities::setBoxToValue(nu8, Box8, 1.0);
    ITHACAutilities::setBoxToValue(nu9, Box9, 1.0);

See also the [**ITHACAutilities::setBoxToValue**](https://mathlab.github.io/ITHACA-FV/namespaceITHACAutilities.html#afe130b8e3e00ec1ab84026ae867f92cb) for more details.

The list of diffusivity fields is set with:

            nu_list.set(0, (nu1).clone());
            nu_list.set(1, (nu2).clone());
            nu_list.set(2, (nu3).clone());
            nu_list.set(3, (nu4).clone());
            nu_list.set(4, (nu5).clone());
            nu_list.set(5, (nu6).clone());
            nu_list.set(6, (nu7).clone());
            nu_list.set(7, (nu8).clone());
            nu_list.set(8, (nu9).clone());
        }

## Definition of the main function

Once the [**tutorial02**](https://mathlab.github.io/ITHACA-FV/classtutorial02.html) class is defined the main function is defined, an example of type tutorial02 is constructed:

    tutorial02 example(argc, argv);

the number of parameter is set:

    example.Pnumber = 9;
    example.setParameters();

the range of the parameters is defined:

    example.mu_range.col(0) = Eigen::MatrixXd::Ones(9, 1) * 0.001;
    example.mu_range.col(1) = Eigen::MatrixXd::Ones(9, 1) * 0.1;

and 500 random combinations of the parameters are generated:

    example.genRandPar(500);
the size of the list of values that are multiplying the affine forms is set:

    example.theta.resize(9);

the source term is defined, the compute_nu and assemble_operator functions are called

    example.SetSource();
    example.compute_nu();
    example.assemble_operator();

then the Offline full order Solve is performed:

    example.offlineSolve();

Once the Offline solve is performed the modes ar obtained using the ITHACAPOD::getModes function:

    ITHACAPOD::getModes(example.Tfield, example.Tmodes, example._T().name(),

and the projection is performed onto the POD modes using 10 modes

    example.project(NmodesTproj);

Once the projection is performed we can construct a reduced object:

    reducedLaplacian reduced(example);

and solve the reduced problem for some values of the parameters:

    for (int i = 0; i < FOM_test.mu.rows(); i++)
    {
        reduced.solveOnline(FOM_test.mu.row(i));
    }

Finally, once the online solve has been performed we can reconstruct the solution:

reduced.reconstruct("./ITHACAoutput/Reconstruction");