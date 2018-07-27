<p align="center">
  <a href="http://mathlab.github.io/ITHACA-FV/" target="_blank" >
    <img alt="ITHACA-FV" src="./docs/logo/ithaca-fv-small.png" width="200" />
  </a>
</p>

## ITHACA-FV - In real Time Highly Advanced Computational Applications for Finite Volumes - ROMs for OpenFOAM ##

<p align="center">
    <a href="https://www.gnu.org/licenses/lgpl-3.0" target="_blank">
        <img alt="Software License" src="https://img.shields.io/badge/License-LGPL%20v3-blue.svg">
    </a>
    <a href="https://travis-ci.org/mathLab/ITHACA-FV" target="_blank">
        <img alt="Build Status" src="https://travis-matrix-badges.herokuapp.com/repos/giovastabile/ITHACA-FV/branches/master/1">
    </a>
      <a href="https://travis-ci.org/mathLab/ITHACA-FV" target="_blank">
        OF6<img alt="Build Status" src="https://travis-matrix-badges.herokuapp.com/repos/giovastabile/ITHACA-FV/branches/master/2">
    </a>
    <a href="https://www.codacy.com/project/giovastabile/ITHACA-FV/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=mathLab/ITHACA-FV&amp;utm_campaign=Badge_Grade_Dashboard">
        <img alt="Codacy Badge" src="https://api.codacy.com/project/badge/Grade/d7ff770dfb954819a0e691ea03de281b">
    </a>
</p>

### 0. Introduction
**ITHACA-FV** is an implementation in **OpenFOAM** of several reduced order modelling techniques. **ITHACA-FV** is designed for [**OpenFOAM 6.0**](https://openfoam.org/version/6) and [**OpenFOAM 5.0**](https://openfoam.org/version/5-0) but it can be easily adapted also to other versions of OpenFOAM. 

**ITHACA-FV** can also be used as a basis for more advanced projects that would like to assess the capability of reduced order models in their existing **OpenFOAM**-based software, thanks to the availability of several reduced order methods and algorithms.

Linear and non-linear algebra operations which are not already implemented in OpenFOAM are performed with the external library [**Eigen**](http://eigen.tuxfamily.org/index.php?title=Main_Page). The source code of Eigen 3.3.4 is provided together with ITHACA-FV and is located in the [src/thirdyparty/Eigen](./src/thirdparty/Eigen) folder.  For the EigenValue decomposition it is also possible to rely on the [**Spectra-0.6.1**](https://spectralib.org/) library and the source code is provided in the [src/thirdyparty/spectra-0.6.1](./src//thirdparty/spectra-0.6.1) folder.

**ITHACA-FV** has been tested on ubuntu 16.04 but can be easily compiled on any linux distribution with a compiled version of OpenFOAM 6.0 or OpenFOAM 5.0. 

### 1. Prerequisites
**ITHACA-FV** requires
* [**OpenFOAM 6.0**](https://openfoam.org/version/6) or
* [**OpenFOAM 5.0**](https://openfoam.org/version/5-0) 


### 2. Installation and usage
First of all you need to source the bashrc file of your installation of **OpenFOAM 6.0** or **OpenFOAM 5.0**. This is of course depending on the location of your OpenFOAM installation and of your particular version of OpenFOAM
```
source $HOME/OpenFOAM/OpenFOAM-6/etc/bashrc
``` 
Then navigate to the folder where you want to install ITHACA-FV such as, for example, the utilities folder of OpenFOAM
```
cd ${FOAM_APP}/utilities
``` 
Now you can clone the **ITHACA-FV** repository inside the selected folder
```
git clone https://github.com/mathLab/ITHACA-FV
```
and you can compile **ITHACA-FV** by navigating inside the src folder and compiling using wmake
```
cd ITHACA-FV
./Allwmake 
```
For a brief description of the classes and methods, you can check the official ITHACA-FV doxygen [documentation](https://mathlab.github.io/ITHACA-FV/).


### 3. [Tutorials](https://mathlab.github.io/ITHACA-FV//examples.html)
Several tutorials are provided the [**tutorials** subfolder](./tutorials).
* [**Tutorial 1**](https://mathlab.github.io/ITHACA-FV/01POD_8C-example.html): In this tutorial it is shown how to perform POD on an already run standard **OpenFOAM** case. 
* [**Tutorial 2**](https://mathlab.github.io/ITHACA-FV/02thermalBlock_8C-example.html): In this tutorial is implemented the development of a parametrized POD-Galerkin ROM for a steady heat transfer problem. The parametrization is on the diffusivity constant. The OpenFOAM full order problem is based on **laplacianFoam**. 
* [**Tutorial 3**](https://mathlab.github.io/ITHACA-FV/03steadyNS_8C-example.html) In this tutorial is implemented the development of a parametrized POD-Galerkin ROM for a steady NS-problem. The parametrization is on the viscosity. The OpenFOAM full order problem is based on **simpleFoam**.
* **Tutorial 4** In this tutorial is implemented the development of a parametrized POD-Galerkin method for an unsteady Navier-Stokes problem. The parametrization is on the viscosity. The OpenFOAM full order problem is based on **pimpleFoam**.


### 4. Authors and contributors
**ITHACA-FV** is currently developed and mantained at [SISSA mathLab](http://mathlab.sissa.it/) by [Dr. Giovanni Stabile](mailto:gstabile@sissa.it) under the supervision of [Prof. Gianluigi Rozza](mailto:gianluigi.rozza@sissa.it)

Contact us by email for further information or questions about **ITHACA-FV**, or open an ''Issue'' on this website. **ITHACA-FV** is at an early development stage, so contributions improving either the code or the documentation are welcome, both as patches or merge requests on this website.

### 5. How to cite
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


and cite the [ITHACA-FV website](http://mathlab.sissa.it/ITHACA-FV).


### 6. License
**ITHACA-FV** is freely available under the GNU LGPL, version 3.

![Google Analytics](https://ga-beacon.appspot.com/UA-66224794-1/rbnics/readme?pixel)
