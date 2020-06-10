# LIFE
### Overview
A lattice Boltzmann-immersed boundary-finite element (LIFE) solver for fluid-structure interaction problems involving slender structures.

### Key Features
 - BGK collision operator
 - central moments collision operator (De Rosis et al. 2019)
 - velocity and pressure boundary conditions set via the regularised method (Latt et al. 2008)
 - convective outlet condition (Yang 2013)
 - implicit immersed boundary method (IBM) for rigid and flexible structures (Li et al. 2016, Pinelli et al. 2010)
 - flexible motion solved via a corotational finite element method (FEM) for slender structures
 - implicit Newmark time integration scheme for FEM
 - force/displacement mapping between LBM and FEM grids
 - strongly-coupled block Gauss-Seidel implicit fluid-structure coupling scheme (Degroote 2013)
 - dynamic under-relaxation with Aitken's method (Irons and Tuck 1969)
 - shared memory parallelisation with OpenMP

# Instruction for Users
### Project Structure
The LIFE project is organised into five folders:

 - **examples**: contains the input files for some example cases
 - **inc**: contains the project header files
 - **input**: contains the geometry.config file that specifies the immersed boundary bodies
 - **obj**: contains the intermediate object files that are generated during the build
 - **src**: contains the project source files

### Case Definition
LIFE uses two files to define the case setup:

 - **params.h**: sets up the domain and flow conditions
 - **geometry.config**: specifies the immersed boundary bodies (rigid/flexible)

The **params.h** file is a header file and therefore the project must be rebuilt every time this file is modified. The **geometry.config** is read by LIFE during initialisation and therefore does not require a rebuild after every modification. To set up a case with no immersed boundary bodies then **geometry.config** can be left blank or removed entirely. If present, it must exist within the **input** folder within the working directory. The instructions for the case setup are given in the comments within each of these files.

There are three steps to adding a new type of geometry for LIFE:

 1. Create a new keyword configuration describing the geometry in the **geometry.config** file.
 2. Add a new condition in the **ObjectsClass::geometryReadIn()** routine to read this configuration.
 3. Implement a new constructor for **IBMBodyClass** to build the new geometry type.

### Dependencies
LIFE requires LAPACK and the Boost library. To install these on Ubuntu type the following into a terminal:

	sudo apt-get install libblas-dev liblapack-dev libboost-all-dev

### Building and running LIFE
A makefile is provided to build LIFE. To use the makefile just navigate to the top level LIFE directory and type:

	make

An executable named **LIFE** will be generated. To run this just type:

	./LIFE

### Post-Processing
LIFE creates a directory called **Results**. This folder contains VTK and log files describing the results of the simulation.

### Restarting
LIFE has support for restarting a simulation from where it left off. To do this just re-run the executable. There is also support for changing the **geometry.config** between restarts to add/remove bodies mid-simulation.

# Contributing
To contribute to LIFE please follow these steps:

 1. Fork the repository and then clone the fork to your local machine.
 2. On the current master branch run the **store-ref-data** script to generate the test data which will be used for testing your new changes.
 3. Implement and test your new features/functionality.
 4. If implementing a substantial new feature then add a new example to the **examples** directory to showcase this.
 5. Run the **run-tests** script to make sure LIFE still produces the same results for the existing example cases.
 6. If the test script passes all of the tests then push the changes and open a pull request to incorporate the changes into this repository.

If your changes require adding/removing options from the **params.h** or **geometry.config** files then make sure to update these files in each of the example case directories otherwise the example/testing scripts will fail.
