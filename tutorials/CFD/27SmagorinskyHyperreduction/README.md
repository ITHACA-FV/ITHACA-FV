# Tutorial 27
## Introduction to tutorial 27

In this tutorial we perform the hyperreduction of the Smagorinsky term from LES simulations inside a reduced order model. The test case is the 2D wake flow past a cylinder at Reynolds 100.

The following image depicts the magnitude of the Smagorinsky term for one time step.


![image info](https://github.com/ITHACA-FV/ITHACA-FV/blob/master/docs/images/SmagorinskyCylinder.png)


First, the LES simulation is ran for 6 vortex shedding cycles, saving 20 snapshots per cycle. The initial conditions correspond to fields value during the steady state (the transitional period has been precomputed). The first 5 vortex shedding cycle are then used for learning, and the final one serves as validation data. Hence, we are trying to extrapolate in time but the simulation is not parametric.

The nonpolynomial term that we are considering is $F(u) = 2 \nabla \cdot ({\nu}_t (u) {\nabla}^s u)$ with ${\nu}_t(u) = C_S {\delta}^2 \|{\nabla}^s u\|$. We denote by $F$ the Smagorinsky term (or `fullStressFunction` in the code) and ${\nu}_t$ by `nut`.
Given the reduced basis $\Phi$ and the associated temporal modes $a$ for the velocity $u$, we are trying to approximate at a low computational cost $\Phi^T F(\Phi a)$.

Two methods are possible:
- DEIM (on $F$ or on $\nu_t$), see [Chaturantabut et al., 2010]
- Empirical Quadrature (keeping or not a linear dependency in $a$), see [Hernández et al, 2017]


# A detailed look into the code
The code details for tutorial 27 are here presented.

## Offline stage

This part is done through the class `tutorial27_offline` in the files `27Offline.H` and `27Offline.C`.
The constructor sets up the offline hyper-reduction problem by reading all the key settings from `StoredParameters` (parameters, problem, the turbulent quantity to be hyper-reduced, the magic points):
```cpp
tutorial27_offline::tutorial27_offline(int argc, char* argv[]):
    m_parameters(new StoredParameters(argc, argv)),
    m_UnsteadyNSTurb(new UnsteadyNSTurb(m_parameters)),
    nModesU(m_parameters->get_nModes()["U"]),
    interpolatedField(m_parameters->get_DEIMInterpolatedField()),
    nModesHR(m_parameters->get_nModes()[interpolatedField]),
    nMagicPoints(m_parameters->get_nMagicPoints()),
    magicPoints(m_parameters->get_magicPoints())
{
    if (!ITHACAutilities::containsSubstring(interpolatedField, "reduced"))
    {
        Info << "Error: Hyperreduction InterpolatedField not valid: " << interpolatedField << endl;
        Info << "Expects reducedFullStressFunction or reducedNut" << endl;
        abort();
    }
}
```
The following part computes the velocity reduced basis using POD, and extracts the velocity spatial modes, the time coefficients, the covariance matrix.
```cpp
void tutorial27_offline::decompose()
{
    ITHACAPOD::PODTemplate<volVectorField> podU(m_parameters, "U");
    
    // compute spatial modes, temporal modes and cov matrix for U
    podU.getModes(spatialModesU, temporalModesU, temporalModesUSimulation, covMatrixU);
    
    // get mean of U field computed when getModes was used
    meanU = new volVectorField(podU.get_mean());
    m_parameters->set_meanU(*meanU);

    volVectorField template_HRInterpField(interpolatedField, m_UnsteadyNSTurb->initSmagFunction());
    m_parameters->set_template_field_fullStressFunction(template_HRInterpField);
```

If the chosen hyper-reduced field is the full stress function, the code builds the hyper-reduction machinery for that quantity. Depending on the method (DEIM or ECP), it constructs a different object, the `TurbDiffusionHyperreduction`. This object has both interpolation and snapshot fields as vector fields (DEIM) or the interpolation vector-valued field and the snapshot scalar-values field (ECP).
```cpp
    if (ITHACAutilities::containsSubstring(interpolatedField, "fullStressFunction"))
    {
        bool mgPoints_0_Neighborhoods_1 = true;

        if (m_parameters->get_HRMethod() == "DEIM")
        {
            TurbDiffusionHyperreduction<volVectorField,volVectorField,volVectorField>*
                turbDiffHR_fullStressFunction = new TurbDiffusionHyperreduction<volVectorField,volVectorField,volVectorField>
                (m_parameters, mgPoints_0_Neighborhoods_1, template_HRInterpField, template_HRInterpField, meanU, spatialModesU);

            turbDiffHR_fullStressFunction->computeTurbDiffusionHyperreduction();
        }
        else if (m_parameters->get_HRMethod() == "ECP")
        {
            volScalarField template_HRSnapshotsField(m_parameters->get_HRSnapshotsField(),
                                                        m_UnsteadyNSTurb->initProjSmagFunction());
            TurbDiffusionHyperreduction<volVectorField,volVectorField,volScalarField>*
                turbDiffHR_projFullStressFunction = new TurbDiffusionHyperreduction<volVectorField,volVectorField,volScalarField>
                (m_parameters, mgPoints_0_Neighborhoods_1, template_HRInterpField, template_HRSnapshotsField, meanU, spatialModesU);

            turbDiffHR_projFullStressFunction->computeTurbDiffusionHyperreduction();
        }
    }
```
In both cases, it calls:

```cpp
computeTurbDiffusionHyperreduction();
```
This is the actual offline hyper-reduction decomposition for the turbulent stress term.
The boolean
```cpp
mgPoints_0_Neighborhoods_1 = true;
```
indicates that this case uses neighborhood-based support around magic points rather than isolated points.

Then, we have the block regarding a different chosen field (the `nut` interpolated field):
```cpp
    else if (ITHACAutilities::containsSubstring(interpolatedField, "nut"))
    {
        bool mgPoints_0_Neighborhoods_1 = false;
        volScalarField template_HRInterpField(interpolatedField, m_UnsteadyNSTurb->initNutFunction());
        volScalarField template_HRSnapshotsField(m_UnsteadyNSTurb->initProjSmagFromNutFunction());      
        TurbDiffusionHyperreduction<volVectorField,volScalarField,volScalarField>*
            turbDiffHR_nut = new TurbDiffusionHyperreduction<volVectorField,volScalarField,volScalarField>
            (m_parameters, mgPoints_0_Neighborhoods_1, template_HRInterpField, template_HRSnapshotsField, meanU, spatialModesU);
        
        turbDiffHR_nut->computeTurbDiffusionHyperreduction();
    }
    else
    {
        Info << "Error: Hyperreduction InterpolatedField not valid: " << interpolatedField << endl;
        Info << "Expects reducedFullStressFunction or reducedNut" << endl;
        abort();
    }
}
```
If the chosen field is nut (eddy viscosity), the code follows a different branch. The interpolated field is scalar; the hyper-reduction snapshots field is also scalar; the velocity basis is still the reduced basis used for projection.

Then it computes the turbulent diffusion hyper-reduction for `nut`.
The flag `mgPoints_0_Neighborhoods_1 = false;` means this branch uses magic points directly, not neighborhoods.

The following block builds the projected hyper-reduction operator for the `fullStressFunction` case. It:

- allocates `K_proj` (theprojected DEIM/ECP matrix);
- loads `K_DEIM`, the previously computed hyper-reduction matrix;
- computes the projected matrix (if not saved)

```cpp
void tutorial27_offline::project()
{
    word folder_DEIM = m_parameters->get_folder_DEIM();

    int d = 3;
    Eigen::MatrixXd K_proj(nModesU, d*(nMagicPoints+1));
    Eigen::MatrixXd K_DEIM = m_parameters->get_K_DEIM();

    if (ITHACAutilities::containsSubstring(interpolatedField, "fullStressFunction"))
    {
        
        word file_projK_DEIM = folder_DEIM + "projected_K_DEIM_" + std::to_string(nModesU) + ".npy";
        bool exist_projK_DEIM = ITHACAutilities::check_file(file_projK_DEIM);

        if (!exist_projK_DEIM)
        {
            
            const volScalarField weight = sqrt(m_parameters->get_volume());
            volVectorField meanSmag(m_parameters->get_meanVectorDEIM());

            for (int i = 0 ; i < nModesU ; i++)
            {
                Eigen::VectorXd K_proj_line_i(d*nMagicPoints);
                if (m_parameters->get_HRMethod() == "DEIM")
                {
                    volVectorField weightedSpatialModesU_OF = weight * spatialModesU[i];
                    Eigen::VectorXd weightedSpatialModesU = Foam2Eigen::field2Eigen(weightedSpatialModesU_OF);
                    K_proj_line_i = K_DEIM.transpose() * weightedSpatialModesU;
                }
                else if (m_parameters->get_HRMethod() == "ECP")
                {
                    K_proj_line_i = K_DEIM.row(i);
                }

                K_proj.row(i).head(d*nMagicPoints) = K_proj_line_i;
                K_proj(i, d*nMagicPoints) = ITHACAutilities::dot_product_L2(spatialModesU[i], meanSmag);
                K_proj(i, d*nMagicPoints + 1) = 0.0;
                K_proj(i, d*nMagicPoints + 2) = 0.0;
            }

            cnpy::save(K_proj, file_projK_DEIM);
            m_parameters->set_projected_K_DEIM(K_proj);
        }
        else
        {
            cnpy::load(K_proj, file_projK_DEIM);
            m_parameters->set_projected_K_DEIM(K_proj);
        }
    }
```

This is instead the offline projection stage for the `nut` hyper-reduction case:
```cpp
    else if (ITHACAutilities::containsSubstring(interpolatedField, "nut"))
    {
        Eigen::Tensor<double, 3> xi_onMagicPts(nModesU + 1, nMagicPoints + 1, nModesU);
        word xi_onMagicPts_file("xi_onMagicPts.npy");
        bool exist_xi_onMagicPts = ITHACAutilities::check_file(folder_DEIM + xi_onMagicPts_file);

        if (!exist_xi_onMagicPts)
        {
            volScalarField meanNut = m_parameters->get_meanScalarDEIM();

            if (m_parameters->get_HRMethod() == "DEIM")
            {
                K_DEIM.array().colwise() /= Foam2Eigen::field2Eigen(m_parameters->get_volume()).array().sqrt();
                for (label i = 0; i < nModesU; i++)
                {
                    for (label q = 0; q < nModesU+1; q++)
                    {
                        for (label j = 0; j < nMagicPoints; j++)
                        {
                            Eigen::VectorXd K_DEIM_col_j = K_DEIM.col(j);
                            if (q == 0)
                            {
                                xi_onMagicPts(q,j,i) = ITHACAutilities::dot_product_L2(spatialModesU[i],
                                    m_UnsteadyNSTurb->diffusion(Foam2Eigen::Eigen2field(meanNut,K_DEIM_col_j),*meanU));
                            }
                            else
                            {
                                xi_onMagicPts(q,j,i) = ITHACAutilities::dot_product_L2(spatialModesU[i],
                                    m_UnsteadyNSTurb->diffusion(Foam2Eigen::Eigen2field(meanNut,K_DEIM_col_j),spatialModesU[q-1]));
                            }
                        }
                        if (q == 0)
                        {
                            xi_onMagicPts(q,nMagicPoints,i) = ITHACAutilities::dot_product_L2(spatialModesU[i],
                                m_UnsteadyNSTurb->diffusion(meanNut,*meanU));
                        }
                        else
                        {
                            xi_onMagicPts(q,nMagicPoints,i) = ITHACAutilities::dot_product_L2(spatialModesU[i],
                                m_UnsteadyNSTurb->diffusion(meanNut,spatialModesU[q-1]));
                        }
                    }
                }
            }
            else if (m_parameters->get_HRMethod() == "ECP")
            {
                PtrList<volTensorField> defTensorModesOnMgPts = m_parameters->get_deformationTensorOfModesOnMagicPoints();
                labelList localMagicPoints(m_parameters->get_localMagicPoints());
                
                for (label i = 0; i < nModesU; i++)
                {
                    for (label q = 0; q < nModesU+1; q++)
                    {
                        for (label j = 0; j < nMagicPoints; j++)
                        {
                            label row_id = ((pow(max(i,q-1),2)+3*max(i,q-1))/2 + min(i+1,q)) * (m_parameters->get_ECPAlgo()!="Global");
                            xi_onMagicPts(q,j,i) = -2*K_DEIM(row_id,j)*(dev(defTensorModesOnMgPts[i+1][localMagicPoints[j]]) 
                                                    && dev(defTensorModesOnMgPts[q][localMagicPoints[j]]));
                        }

                        if (q == 0)
                        {
                            xi_onMagicPts(q,nMagicPoints,i) = -fvc::domainIntegrate(m_UnsteadyNSTurb->
                                    projDiffusionIBP(meanNut, *meanU, spatialModesU[i])).value();
                        }
                        else
                        {
                            xi_onMagicPts(q,nMagicPoints,i) = -fvc::domainIntegrate(m_UnsteadyNSTurb->
                                    projDiffusionIBP(meanNut, spatialModesU[q-1], spatialModesU[i])).value();
                        }
                    }
                }
            }
            cnpy::save(xi_onMagicPts, folder_DEIM + xi_onMagicPts_file);
            m_parameters->set_xi_onMagicPts(xi_onMagicPts);
        }
        else
        {
            cnpy::load(xi_onMagicPts, folder_DEIM + xi_onMagicPts_file);
            m_parameters->set_xi_onMagicPts(xi_onMagicPts);
        }
    }
}
```
Note that instead of a 2D matrix, it builds a 3D tensor:
```cpp
xi_onMagicPts(nModesU + 1, nMagicPoints + 1, nModesU)
```
This tensor stores how the turbulent diffusion operator couples:

- trial functions ($q$)
- magic points ($j$)
- test modes ($i$)

In the DEIM case: for each mode combination and each magic point, the code computes an L2 projection of the diffusion operator built from one DEIM basis column; either the mean velocity or one POD mode.
The extra slot `nMagicPoints` stores the contribution from the mean `nut`.

In the ECP case: the code uses: deformation tensors evaluated on magic points, cubature weights from `K_DEIM`, integrated-by-parts projected diffusion terms for the mean contribution.
Then the tensor is saved or reloaded from disk.

So the overall purpose is to prepare all the reduced and hyper-reduced ingredients needed for a later online ROM with efficient turbulent diffusion treatment.


## Online stage

This part is done through the class `tutorial27_online` in the files `27Online.H` and `27Online.C`.

At first, the online constructor loads the offline info (parameters, modes, coefficients, magic points, Smagorinsky constants) and initialize the online Hyper-reduction object:
```cpp
tutorial27_online::tutorial27_online(StoredParameters* parameters):
    m_parameters(parameters),
    m_UnsteadyNSTurb(new UnsteadyNSTurb(m_parameters)),
    runTime2(Foam::Time::controlDictName, ".", m_parameters->get_casenameData()),
    meanU(new volVectorField(m_parameters->get_meanU())),
    nModesU(m_parameters->get_nModes()["U"]),
    interpolatedField(m_parameters->get_DEIMInterpolatedField()),
    nModesHR(m_parameters->get_nModes()[interpolatedField]),
    nMagicPoints(m_parameters->get_nMagicPoints()),
    magicPoints(m_parameters->get_magicPoints()),
    localMagicPoints(m_parameters->get_localMagicPoints()),
    Ck(m_parameters->get_Ck()),
    Ce(m_parameters->get_Ce())
    ITHACAstream::read_fields(spatialModesU, "U", "./ITHACAoutput/spatialModes_" + std::to_string(nModesU) + "modes/");
    massMatrixInv = ITHACAutilities::getMassMatrix(spatialModesU).completeOrthogonalDecomposition().pseudoInverse();
    cnpy::load(temporalModesUSimulation, "./ITHACAoutput/temporalModesSimulation_" + std::to_string(nModesU) + "modes/U.npy");

    weightAtMg = Eigen::VectorXd(nMagicPoints);
    for (int i = 0; i < nMagicPoints; i++)
    {
        weightAtMg(i) = std::pow(m_parameters->get_submesh().V()[localMagicPoints[i]], 0.5);
    }
```
This part specializes the online setup depending on which quantity is being hyper-reduced:
```cpp
    if (ITHACAutilities::containsSubstring(interpolatedField, "fullStressFunction"))
    {
        d = 3;
        meanSmagOnMagicNeighborhoods = new volVectorField(m_parameters->get_meanVectorDEIMMagic());
        for (auto& defTensorOfMode:m_parameters->get_deformationTensorOfModesOnMagicNeighborhoods())
        {
            defTensorOfModesAtMg.push_back(defTensorOfMode);
        }
    }
    else if (ITHACAutilities::containsSubstring(interpolatedField, "nut"))
    {
        d = 1;
        meanNutOnMagicPoints = new volScalarField(m_parameters->get_meanScalarDEIMMagic());
        for (auto& defTensorOfMode:m_parameters->get_deformationTensorOfModesOnMagicPoints())
        {
            defTensorOfModesAtMg.push_back(defTensorOfMode);
        }
    }

    nut0 = new volScalarField("nut", volScalarField(
        IOobject("nut0", m_parameters->get_casenameData() + "0", m_parameters->get_submesh(), IOobject::MUST_READ),
        m_parameters->get_submesh()));

    delta = new volScalarField(m_parameters->get_magicDelta());
    aaa = new volScalarField(Ce/(*delta));
    inv2aaa = new volScalarField(1/(2*(*aaa)));

    fullDefField = new volTensorField(defTensorOfModesAtMg[0]);
    bbb = new volScalarField((2.0/3.0)*tr(*fullDefField));
    ccc = new volScalarField(2*Ck*(*delta)*(dev(*fullDefField) && (*fullDefField)));
    sqrtk = new volScalarField((-(*bbb) + sqrt(sqr(*bbb) + 4*(*aaa)*(*ccc)))*(*inv2aaa));
    nut = new volScalarField(Ck*(*delta)*(*sqrtk));
}
```
If the interpolated field is `fullStressFunction`:

- the online output dimension is set to d = 3
- the mean Smagorinsky contribution on the magic neighborhoods is loaded
- the deformation tensors of the velocity modes on the magic neighborhoods are stored

If the interpolated field is `nut`:

- the output dimension is set to d = 1
- the mean eddy viscosity on the magic points is loaded
- the deformation tensors on the magic points are stored.

This other function reconstructs approximate `nut` from reduced coefficients:
```cpp
void tutorial27_online::evaluateApproxNut(const Eigen::VectorXd& reducedCoeffs)
{
    *fullDefField = defTensorOfModesAtMg[0];
    for (int c = 0 ; c < reducedCoeffs.size() ; c++)
    {
        *fullDefField += reducedCoeffs(c)*defTensorOfModesAtMg[c+1];
    }
    *bbb = (2.0/3.0) * tr(*fullDefField);
    *ccc = 2*Ck * (*delta) * ((dev(*fullDefField) && (*fullDefField)));
    *sqrtk = mag(-(*bbb) + sqrt(sqr(*bbb) + 4*(*aaa) * (*ccc))) * (*inv2aaa);
    *nut = Ck * (*delta) * (*sqrtk);
    for (int i = 0; i < nut0->size(); i++)
    {
        (*nut0)[i] = (*nut)[i];
    }
    nut0->correctBoundaryConditions();
}
```
It rebuilds the deformation tensor, compure the scalar quantities entering the Smagorinsky closure (`bbb`, `ccc`, `sqrtk`), evalues the approximate eddy viscosity, adn copy the result into `nut0`.


Then, there are two other functions which build the online hyper-reduced vector evaluated at the selected points: `computeApproxSmagMg` is used for the `fullStressFunction` case, while `computeApproxNutMg` is used for the `nut` case.
```cpp
Eigen::VectorXd tutorial27_online::computeApproxSmagMg(const Eigen::VectorXd& reducedCoeffs)
{
    evaluateApproxNut(reducedCoeffs);
    volVectorField stressField = fvc::div(2*(*nut0) * dev(*fullDefField)) - *meanSmagOnMagicNeighborhoods;

    Eigen::VectorXd output = Eigen::VectorXd::Zero(d * (nMagicPoints + 1));
    for (int i = 0; i < nMagicPoints; i++)
    {
        vector stressVector = weightAtMg(i) * (stressField)[localMagicPoints[i]];

        output(i * d) = stressVector.x();
        output(i * d + 1) = stressVector.y();
        output(i * d + 2) = stressVector.z();
    }
    output(d * nMagicPoints) = 1.0;
    return output;
}
```

```cpp
Eigen::VectorXd tutorial27_online::computeApproxNutMg(const Eigen::VectorXd& reducedCoeffs)
{
    evaluateApproxNut(reducedCoeffs);
    *nut0 -= *meanNutOnMagicPoints;

    Eigen::VectorXd output(nMagicPoints + 1);
    for (int i = 0; i < nMagicPoints; i++)
    {
        output(i) = weightAtMg(i) * (*nut0)[localMagicPoints[i]];
    }
    output(nMagicPoints) = 1.0;
    return output;
}
```

Then, the ROM coefficients of the Smagorinsky term are predicted, again in the two cases:
```cpp
Eigen::VectorXd tutorial27_online::predictSmagROMCoeffs(const Eigen::VectorXd& reducedCoeffs)
{
    if (ITHACAutilities::containsSubstring(interpolatedField, "fullStressFunction"))
    {
        Eigen::MatrixXd projected_K_DEIM = m_parameters->get_projected_K_DEIM();
        Eigen::VectorXd smagVect(computeApproxSmagMg(reducedCoeffs));

        return massMatrixInv*projected_K_DEIM*smagVect;
    }

    else if (ITHACAutilities::containsSubstring(interpolatedField, "nut"))
    {
        Eigen::VectorXd nutVect(computeApproxNutMg(reducedCoeffs));
        Eigen::Tensor<double, 3 > xi_complete_onMagicPts(m_parameters->get_xi_onMagicPts());
        Eigen::VectorXd b_extended(reducedCoeffs.size() + 1);
        b_extended << 1.0,reducedCoeffs;

        Eigen::VectorXd output(nModesU);
        for (int i = 0; i < nModesU; i++)
        {
            output(i) = 0;
            for (int j = 0 ; j < nMagicPoints +1; j++)
            {
                for (int q = 0 ; q < nModesU + 1 ; q++)
                {
                    output(i) += xi_complete_onMagicPts(q,j,i) * b_extended(q) * nutVect(j);
                }
            }
        }
        return massMatrixInv*output;
    }
}
```
In the `fullStressFunction case`: it computes the sampled stress vector `smagVect`, multiplies it by the projected DEIM matrix, applies the inverse reduced mass matrix to obtain the reduced coefficients.
In the `nut` case: it computes the sampled eddy-viscosity vector `nutVect`, loads the precomputed tensor `xi_onMagicPts`, extends the reduced state with a leading 1.0 to include the mean contribution, assembles the reduced turbulent term coefficient-by-coefficient through a triple sum, and apply the inverse reduced mass matrix.

Finally, there are some functions to define the reconstruction of the solution, and the comparison with the true solution.
```cpp
volVectorField tutorial27_online::computeROMproj_fromCoeffs(const Eigen::VectorXd& reducedCoeffs, bool stressUnit) {...}
Eigen::VectorXd tutorial27_online::computeROMcoeffs_fromFullDim(volVectorField& f_full){...}
volVectorField tutorial27_online::computeSmagTerm_fromChronos(const Eigen::VectorXd& reducedCoeffs){...}
void tutorial27_online::prediction(){...}
```
In particular, the helper functions:

- `computeROMproj_fromCoeffs`: reconstruct a vector field from reduced coefficients
- `computeROMcoeffs_fromFullDim`: project a full-order vector field onto the velocity basis
- `computeSmagTerm_fromChronos`: compute the exact Smagorinsky term starting from the reduced velocity reconstruction,
`prediction` is instead the online validation loop, reconstructing the solutions, saving them, and computing the L2 errors.

## Changing the hyperreduction method and parameters

All these information are declared in the `ITHACAdict` file inside the system folder. Here are the tunable entrances:

- `HyperReduction`: can be set to `DEIM` or `ECP` to choose the hyperreduction method.
- `DEIMInterpolatedField`: can be set to `reducedFullStressFunction` or `reducedNut` to select the hyperreduced field. In the case of the quadrature method, the former option preserves a linear dependency with the temporal modes while the latter does not.
- `nmodes`: integer specifying the number of modes kept after the POD. Can be changed both for the velocity and for the hyperreduced fields. The number of magic points automatically equals the number of modes for the hyperreduced field.

## The main file
The main file is quite simple now, because it is just performing the offline, online, and the validation steps:
```cpp
int main(int argc, char* argv[])
{
    tutorial27_offline offlinePart(argc, argv);

    // Compute POD on velocity and on the hyperreduced term. Then select the magic points
    offlinePart.decompose();

    // Compute the matrix for fast online prediction of the nonpolynomial term 
    offlinePart.project();

    tutorial27_online onlinePart(offlinePart.m_parameters);

    // Perfom the prediction on the test time steps and compare the results to the reference
    onlinePart.prediction();
}
```


## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/27SmagorinskyHyperreduction/27SmagorinskyHyperreduction.C).