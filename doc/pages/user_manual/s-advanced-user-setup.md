title: Advanced User Setup
---


Advanced users can have additional control over a simulation without having to modify the source code and recompile the code. To do that, the user can provide a set of routines that are called in different stages of the simulation via the Problem file (*ProblemFile.f90*). A description of the routines of the Problem File can be found in the following section.

## Routines of the Problem File: *ProblemFile.f90* 


- UserDefinedStartup: Called before any other routines

- UserDefinedFinalSetup: Called after the mesh is read in to allow mesh related initializations or memory allocations.

- UserDefinedInitialCondition: called to set the initial condition for the flow. By default it sets an uniform initial condition, but the user can change it.

- UserDefinedState1, UserDefinedNeumann: Used to define an user-defined boundary condition.

- UserDefinedPeriodicOperation: Called before every time-step to allow periodic operations to be performed.

- UserDefinedSourceTermNS: Called to apply source terms to the equation.

- UserDefinedFinalize: Called after the solution computed to allow, for example error tests to be performed.

- UserDefinedTermination: Called at the the end of the main driver after everything else is done.


## Compiling the Problem File

The Problem File file must be compiled using a specific Makefile that links it with the libraries of the code. If you are using the `horses/dev` environment module, you can get templates of the *Problemfile.f90* and *Makefile* with the following commands:

```bash
	$ horses-get-makefile
	$ horses-get-problemfile
```

Otherwise, search the test cases for examples.\\

To run a simulation using user-defined operations, create a folder called SETUP on the path were the simulation is going to be run. Then, store the modified *ProblemFile.f90* and the *Makefile* in SETUP, and compile using:

```bash
	$ make <<Options>>
```

where again the options are (bold are default):

- MODE=DEBUG/HPC/**RELEASE**
- COMPILER=ifort/**gfortran**
- COMM=PARALLEL/**SEQUENTIAL**
- ENABLE\_THREADS=NO/**YES**

