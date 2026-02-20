# Introduction to tutorial 27

In this tutorial we perform the hyperreduction of the Smagorinsky term from LES simulations inside a reduced order model. The test case is the 2D wake flow past a cylinder at Reynolds 100.

The following image depicts the magnitude of the Smagorinsky term for one time step.


![image info](https://github.com/ITHACA-FV/ITHACA-FV/blob/master/docs/images/SmagorinskyCylinder.png)


First, the LES simulation is ran for 6 vortex shedding cycles, saving 20 snapshots per cycle. The initial conditions correspond to fields value during the steady state (the transitional period has been precomputed). The first 5 vortex shedding cycle are then used for learning, and the final one serves as validation data. Hence, we are trying to extrapolate in time but the simulation is not parametric.

The nonpolynomial term that we are considering is $F(u) = 2\nabla \cdot (\nu_{t}(u) \nabla^s u)$
with $\nu_{t}(u) = C_S \delta^2 ||\nabla^s u||$. We denote $F$ by Smagorinsky term (or fullStressFunction in the code) and $\nu_{t}$ by nut.
Given the reduced basis $\Phi$ and the associated temporal modes $a$ for the velocity $u$, we are trying to approximate at a low computational cost $\Phi^T F(\Phi a)$.

Two methods are possible:
- DEIM (on $F$ or on $\nu_t$), see [Chaturantabut et al., 2010]
- Empirical Quadrature (keeping or not a linear dependency in $a$), see [Hernández et al, 2017]


## A look into the code


### Offline stage

This part is done through the class tutorial27_offline in the files 27Offline.H and 27Offline.C.
We first compute the POD reduced basis $\Phi$ for the velocity and its associated temporal modes.
Then, we compute a POD reduced basis for the hyperreduced field $F(\Phi a)$ or $\nu_t(\Phi a)$.
Magic points are then selected, depending on the chosen hyperreduction method.
Finally, the method project() precomputes the matrix such that online the prediction of $\Phi^T F(\Phi a)$ is fast.


### Online stage

This part is done through the class tutorial27_online in the files 27Online.H and 27Online.C.
Two operations must be performed at each time step to obtain the final prediction.
First, the nonpolynomial term is evaluated on the submesh of magic points.
Next, given these evaluations and the precomputed matrix, the approximation is computed.
These computations are performed for the 20 time steps of the testing interval in the projection() method.


### Changing the hyperreduction method and parameters

All these information are declared in the ITHACAdict file inside the system folder. Here are the tunable entrances:
- HyperReduction: can be set to "DEIM" or "ECP" to choose the hyperreduction method.
- DEIMInterpolatedField: can be set to "reducedFullStressFunction" or "reducedNut" to select the hyperreduced field. In the case of the quadrature method, the former option preserves a linear dependency with the temporal modes while the latter does not.
- nmodes: integer specifying the number of modes kept after the POD. Can be changed both for the velocity U and for the hyperreduced fields. The number of magic points automatically equals the number of modes for the hyperreduced field.


### Results

Results are stored in the folder ITHACAoutput/Online. Additional subfolders are created depending on the chosen parameters (see above section).
Saved results include:
- a file 'exactCoeffs.npy' containing the (vectorial) values of $\Phi^T F(\Phi a)$ at each time step of the test interval
- a file 'predictionCoeffs.npy' containing the approximated values $\widehat{\Phi^T F(\Phi a)}$ at each time step of the test interval
- a file 'relativeSmagProjError.npy' containing the L2 relative error between the reference $\Phi\Phi^T F(\Phi a)$ and the prediction
- subfolders containg the OpenFOAM field for the reference $\Phi\Phi^T F(\Phi a)$ and the prediction, for visualization