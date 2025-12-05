# Tutorial 09
Implementation of tutorial 9 which presents DEIM for a Heat Conduction Problem.

## Introduction
In this tutorial we implement a test where we use the Discrete Empirical Interpolation Method for a case where we have a non-linear dependency with respect to the input parameters.
The following image illustrates the computational domain which is the same as the previous example:

![domain](../../../docs/images/domain_deim.png)

The physical problem is given by a heat transfer problem which is described by the Poisson equation:
$$ \nabla \cdot (\nu \nabla T) = S $$

The parametric diffusivity is described by a parametric Gaussian function:
$$ \nu(\mathbf{x},\mathbf{\mu}) = e^{-2(x-\mu_x-1)^2 - 2(y-\mu_y-0.5)^2}, $$

The problem is then discretized as:
$$ A(\mu)T = b.$$

In this case, even if the problem is linear, due to non-linearity with respect to the input parameter of the conductivity constant it is not possible to have an affine decomposition of the discretized differential operator.

We seek therefore an approximate affine expansion of the differential operator of this type:
$$ A(\mu) = \sum_{i = 1}^{N_D} \theta_i(\mu) A_i  $$
using the Discrete Empirical Interpolation Method

## The plain program
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/09DEIM_ROM/09DEIM_ROM.C).
