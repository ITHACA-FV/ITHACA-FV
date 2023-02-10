<p align="center">
  <a href="http://ithaca-fv.github.io/ITHACA-FV/" target="_blank" >
    <img alt="ITHACA-FV" src="./docs/logo/ithaca-fv-small.png" width="200" />
  </a>
</p>

## ITHACA-FV - In real Time Highly Advanced Computational Applications for Finite Volumes - ROMs for OpenFOAM ##

<p align="center">
<a href="https://www.gnu.org/licenses/lgpl-3.0" target="_blank">
<img alt="Software License" src="https://img.shields.io/badge/License-LGPL%20v3-blue.svg">
</a>
<a href="https://www.codacy.com/gh/ithaca-fv/ITHACA-FV/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ithaca-fv/ITHACA-FV&amp;utm_campaign=Badge_Grade">
<img alt="Codacy Badge" src="https://app.codacy.com/project/badge/Grade/f51e1de0fd7c42a6a0be13889cbd4b0e">
</a>
</p>
<p align="center">
<a href="https://github.com/ithaca-fv/ITHACA-FV/actions?query=workflow%3AOF1812"><img alt="OF1812" src="https://github.com/ithaca-fv/ITHACA-FV/workflows/OF1812/badge.svg"></a>
<a href="https://github.com/ithaca-fv/ITHACA-FV/actions?query=workflow%3AOF1906"><img alt="OF1906" src="https://github.com/ithaca-fv/ITHACA-FV/workflows/OF1906/badge.svg"></a>
<a href="https://github.com/ithaca-fv/ITHACA-FV/actions?query=workflow%3AOF1912"><img alt="OF1912" src="https://github.com/ithaca-fv/ITHACA-FV/workflows/OF1912/badge.svg"></a>
<a href="https://github.com/ithaca-fv/ITHACA-FV/actions?query=workflow%3AOF2006"><img alt="OF2006" src="https://github.com/ithaca-fv/ITHACA-FV/workflows/OF2006/badge.svg"></a>
</p>

### 0. Introduction
**ITHACA-FV** is an implementation in **OpenFOAM** of several reduced order modelling techniques. **ITHACA-FV** is designed for [**OpenFOAM v2106**](https://www.openfoam.com/releases/openfoam-v2106/) and older openfoam.com versions but it can be easily adapted also to other versions of OpenFOAM.

**ITHACA-FV** can also be used as a basis for more advanced projects that would like to assess the capability of reduced order models in their existing **OpenFOAM**-based software, thanks to the availability of several reduced order methods and algorithms.

Linear and non-linear algebra operations which are not already implemented in OpenFOAM are performed with the external library [**Eigen**](http://eigen.tuxfamily.org/index.php?title=Main_Page). The source code of Eigen 3.3.7 is provided together with ITHACA-FV and is located in the [src/thirdyparty/Eigen](./src/thirdparty/Eigen) folder.  For the EigenValue decomposition it is also possible to rely on the [**Spectra-0.7.0**](https://spectralib.org/) library and the source code is provided in the [src/thirdyparty/spectra](./src//thirdparty/spectra) folder. Numerical optimization can be performed using the external library [**OptimLib**](https://www.kthohr.com/optimlib.html) and the header based source code is provided in the [src/thirdyparty/OptimLib](./src/thirdparty/OptimLib) folder.

**ITHACA-FV** has been tested on ubuntu 16.04, CentOS 7, ArchLinux but can be easily compiled on any linux distribution with a compiled version of OpenFOAM 1812, OpenFOAM 1906, OpenFOAM 1912, OpenFOAM 2006, OpenFOAM 2012 and OpenFOAM 2106.

### 1. Prerequisites
**ITHACA-FV** requires
*   [**OpenFOAM 1812**](https://www.openfoam.com/releases/openfoam-v1812/) or
*   [**OpenFOAM 1906**](https://www.openfoam.com/releases/openfoam-v1906/) or
*   [**OpenFOAM 1912**](https://www.openfoam.com/releases/openfoam-v1912/) or
*   [**OpenFOAM 2006**](https://www.openfoam.com/releases/openfoam-v2006/) or
*   [**OpenFOAM 2012**](https://www.openfoam.com/releases/openfoam-v2012/) or
*   [**OpenFOAM 2106**](https://www.openfoam.com/news/main-news/openfoam-v2106)


### 2. Installation and usage
First of all you need to source the bashrc file of your installation **OpenFOAM**. This is of course depending on the location of your OpenFOAM installation and of your particular version of OpenFOAM
```
source $HOME/OpenFOAM/OpenFOAM-v2106/etc/bashrc
```
Then navigate to the folder where you want to install ITHACA-FV such as, for example, your HOME folder
```
cd ~
```
Now you can clone the **ITHACA-FV** repository inside the selected folder
```
git clone --depth 1 https://github.com/ITHACA-FV/ITHACA-FV
```
and you can compile **ITHACA-FV** by navigating inside the src folder, sourcing the bashrc file of ITHACA-FV and compiling using wmake:
```
cd ITHACA-FV
git submodule update --init
source etc/bashrc
./Allwmake
```
By default the `Allwmake` routine will only compile the ITHACA-FV library with standard packages and without applications, unitTests and tutorials. For more informations concerning the possible options one can execute the command with the `-h` flag.

For example in order to compile the library with also all the applications, tutorials and unit tests the following commands needs to be executed:
```
cd ITHACA-FV
source etc/bashrc
./Allwmake -tau
```

To enable parallel compilation the `-j` flag followed by the number of processors that you want to use can be added. For example, to compile the library with `4` processors the following command can be issued:
```
cd ITHACA-FV
source etc/bashrc
./Allwmake -tau -j 4
```

In the near future the ITHACA-FV will also be linked with the pytorch package for machine learning. Some basic functions are already available. In order to compile these additional functionalities one will need to have torch installed and compile the library with the `-m` options. Moreover one will need to install a version of libtorch with ABI enabled. The one available at the following link for example has it:
```
    curl https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-1.4.0%2Bcpu.zip > libtorch.zip && \
    unzip libtorch.zip -d opt/ && \
```

For a brief description of the classes and methods, you can check the official ITHACA-FV doxygen [documentation](https://ithaca-fv.github.io/ITHACA-FV).


### 3. Docker Image

Docker images for linux/arm64 and linux/amd64  are available on [Docker Hub](https://hub.docker.com/u/ithacafv).
These images are based **OpenFOAM-v2106**, and provided an isolated environment, where you can find a compiled version of the master branch of **ITHACA-FV**.
In order to pull the image, run the following command:

```
docker pull ithacafv/ithacafv:manifest-latest
```

Once the image is downloaded, you can start the container and mount the $HOME directory by executing:

```
docker run -ti --rm -v "${HOME}:/home/ithacafv/${USER}" ithacafv/ithacafv:manifest-latest
```

### 4. Singularity


Create a local directory (which you own), followed by using the variables, `SINGULARITY_TMPDIR` and `SINGULARITY_CACHEDIR` for temporary and cache directory during built proccess. This step is highly recommended if you do NOT have `root` access, as singularity uses a temporary directory to build the squashfs filesystem, and this temp space needs to be large enough to hold the entire resulting Singularity image.

```
export SINGULARITY_TMPDIR=$HOME/mycontainter
```
and 

```
export SINGULARITY_CACHEDIR=$HOME/mycontainter
```

Buidling singularity image file `.sif` from the docker image, which is build and stored in `$HOME/mycontainter`, The below image used is based **OpenFOAM-v2106**, and where you can find a compiled version of the master branch of **ITHACA-FV**.

```
singularity build ithacafv.sif docker://ithacafv/ithacafv:manifest-latest
```

To view / list all the images/cache, 

```
singularity cache list -v
```

To run the singulairty image interactively, use `shell` from your working directory. This mounts your working directory to the container. Activate the the OpenFOAM environment after`shell`, using `source /usr/lib/openfoam/openfoam2106/etc/bashrc` (or, `source /etc/bash.bashrc`)

```
singularity shell $HOME/mycontainter/ithacafv.sif
```

To build singularity images from scratch use defination files, `singularity-reciepe.def` provided in singularity directory. This require `sudo` prilvildges. Recommened to build in detached mode, by passing the flag `-d`.

```
sudo singularity build -d ithacafv.sif singularity/singularity-reciepe.def
```

Additionally, `--fakeroot` can be passed if you do NOT have root access for the build.

```
singularity build --fakeroot ithacafv.sif singularity/singularity-reciepe.def
```

Running singualrity in batch mode, you can add the following in your the batch script, 

```
singularity exec <image>.sif <command>
```

For example, running the any tutorial, add the following in the batch script,

```
singularity exec ithacafv.sif /bin/bash -c Of.sh
``` 

the shell script `Of.sh` is provided within the singularity directory. A sample batch scipt is also proived. Please note, over slurm machine, mpi binds without passing any flags.



### 5. [Tutorials](https://ithaca-fv.github.io/ITHACA-FV//examples.html)
Several tutorials are provided in the [**tutorials** subfolder](./tutorials).
*   [**CFD/Tutorial 1**](https://github.com/ithaca-fv/ITHACA-FV/tree/master/tutorials/CFD/01POD): In this tutorial it is shown how to perform POD on an already run standard **OpenFOAM** case.
*   [**CFD/Tutorial 2**](https://github.com/ithaca-fv/ITHACA-FV/tree/master/tutorials/CFD/02thermalBlock): In this tutorial the development of a parametrized POD-Galerkin ROM for a steady heat transfer problem is implemented. The parametrization is on the diffusivity constant. The OpenFOAM full order problem is based on **laplacianFoam**.
*   [**CFD/Tutorial 3**](https://github.com/ithaca-fv/ITHACA-FV/tree/master/tutorials/CFD/03steadyNS): In this tutorial, the development of a parametrized POD-Galerkin ROM for a steady NS-problem. The parametrization is on the viscosity is implemented. The OpenFOAM full order problem is based on **simpleFoam**.
*   [**CFD/Tutorial 4**](https://ithaca-fv.github.io/ITHACA-FV/04unsteadyNS_8C-example.html): In this tutorial, the development of a parametrized POD-Galerkin method for an unsteady Navier-Stokes problem is implemented. The parametrization is on the viscosity. The OpenFOAM full order problem is based on **pimpleFoam**.

*   [**CFD/Tutorial 6**](https://ithaca-fv.github.io/ITHACA-FV/06POD_RBF_8C-example.html): This tutorial presents the application of data driven POD-Galerkin model to RANS equation in steady state setting. The components of the velocity at the inlet are paramertized. The OpenFOAM full order problem is based on **simpleFoam**.

*   [**CFD/Tutorial 8**](https://ithaca-fv.github.io/ITHACA-FV/08DEIM_8C-example.html): In this tutorial we propose an example concerning the usage of the Discrete Empirical Interpolation Methods for the reconstruction of a non-linear function. In this case we do not perform model reduction but we test just the motodology a non-linear function.

*   [**CFD/Tutorial 9**](https://ithaca-fv.github.io/ITHACA-FV/09DEIM_ROM_8C-example.html): In this tutorial we propose an example concerning the usage of the Discrete Empirical Interpolation Methods for model reduction purposes. The non-linearity is in the forcing term of a heat transfer problem. The OpenFOAM full order problem is based on **laplacianFoam**.


### 6. Authors and contributors
**ITHACA-FV** is currently developed and maintained at the [University of Urbino Carlo Bo](https://www.uniurb.it/) by [Dr. Giovanni Stabile](www.giovannistabile.com) and at [SISSA mathLab](http://mathlab.sissa.it/) in collaboration with [Prof. Gianluigi Rozza](https://people.sissa.it/~grozza/)'s group.

Contact us by email for further information or questions about **ITHACA-FV**, or open an ''Issue'' on this website. **ITHACA-FV** is at an early development stage, so contributions improving either the code or the documentation are welcome, both as patches or merge requests on this website. If you are willing to contribute please follow the [developer instructions](https://github.com/ithaca-fv/ITHACA-FV/tree/master/src).

### 7. How to cite
Most of the theoretical aspects behind **ITHACA-FV** are deeply explained in [<b> Stabile2017CAIM </b>](https://arxiv.org/pdf/1701.03424.pdf) and [<b> Stabile2017CAF </b>](https://arxiv.org/pdf/1710.11580.pdf).
For this reason, if you use this software, please consider citing the mentioned works, reported in the following bibtex entries:
```
@Article{Stabile2017CAIM,
Title                    = {{POD-Galerkin reduced order methods for CFD using Finite Volume Discretisation: vortex shedding around a circular cylinder}},
Author                   = {Stabile, Giovanni and Hijazi, Saddam and Mola, Andrea and Lorenzi, Stefano and Rozza, Gianluigi},
Journal                  = {Communications in Applied and Industrial Mathematics},
Year                     = {(2017)},
Volume                   = {8},
Number                   = {1},
pages                    = {210-236},
Doi                      = {10.1515/caim-2017-0011}}
```

```
@Article{Stabile2017CAF,
Title                    = {{Finite volume POD-Galerkin stabilised reduced order methods for the parametrised incompressible Navier-Stokes equations}},
Author                   = {Stabile, Giovanni and Rozza, Gianluigi},
Journal                  = {Computers & Fluids},
Year                     = {2018},
Doi                      = {10.1016/j.compfluid.2018.01.035}}
```


and cite the [ITHACA-FV website](https://ithaca-fv.github.io/ITHACA-FV).


### 8. License
**ITHACA-FV** is freely available under the GNU LGPL, version 3.
