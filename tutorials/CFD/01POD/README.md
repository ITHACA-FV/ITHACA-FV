# Tutorial 01
## Introduction

In this tutorial the POD is applied on a simple lid-driven cavity example.

The case is solved using icoFoam and then the modes are extracted using
the `perform_POD` utility. Run the code using the standard icoFoam solver from OpenFOAM
```cpp
    icoFoam
```
Once the solution is performed you can the ITHACA-FV utility to extract the modes running:
```cpp
    perform_POD
```

The modes, together with the eigenvalues and the cumulative eigenvalues are stored inside
the `ITHACAoutput/POD` folder.

You can eventually use also the Allrun file to perform both operation together. For more
details see also the `perform_POD.C` file and check the `system/ITHACAPODdict` file in the
system folder to understand how this file should be prepared.

## The ITHACAPODdict file

The ITHACAPODdict file located in the system folder is used to define the characteristics
of the fields on which you want to perform the POD and how many modes you want to extract.
Let's have a detailed look to it.

With the following lines you decide on which fields you want to perform the POD, put the
name of the field followed by the `_pod` string:
```cpp
    fields
    (
    U_pod
    p_pod
    );
```
Then for each field we have to specify the exact file name of the field, how many modes
we want to extract and the type of the field (vector or scalar). In this case we have 10
modes for both velocity and pressure field:
```cpp
    U_pod
    {
    field_name U;
    nmodes 2;
    field_type vector;
    }
    p_pod
    {
    field_name p;
    nmodes 2;
    field_type scalar;
    }
```
Then you have to select the `InitialTime` from which you want to start acquiring the snapshots and the `FinalTime`
```cpp
    InitialTime 0;
    FinalTime 50;
```

Eventually, instead of the finalTime you could define the number of snapshots:
```cpp
    Nsnapshots 40
```