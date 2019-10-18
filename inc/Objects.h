/*
    LIFE: Lattice boltzmann-Immersed boundary-Finite Element
    Copyright (C) 2019 Joseph O'Connor

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OBJECTS_H // OBJECTS_H
#define OBJECTS_H

// Includes
#include "params.h"
#include "IBMBody.h"
#include "IBMNode.h"

// Forward declarations
class GridClass;

// Objects class
class ObjectsClass {

	// Friend classes
	friend class GridClass;
	friend class IBMBodyClass;

	// Default constructor and destructor
public:
	ObjectsClass() {gPtr = NULL; nFlex = 0; simDOFs = 0; hasIBM = false; hasFlex = false; subIt = 0; relax = relaxMax; subRes = 0.0; subNum = 0.0; subDen = 0.0;};
	~ObjectsClass() {};

	// Custom constructor for reading in geometries
	ObjectsClass(GridClass &grid);

	// Public members
public:

	// Pointer to grid
	GridClass *gPtr;

	// Flag to say if there are IBM bodies
	bool hasIBM;

	// Private members
private:

	// Set helper members
	bool hasFlex;					// Flag to say if there are flexible bodies in case
	int nFlex;						// Number of flexible bodies
	int simDOFs;					// Number of FEM DOFs in full system

	// Vectors of IBM bodies and nodes
	vector<IBMBodyClass> iBody;
	vector<IBMNodeClass> iNode;

	// Subiteration loop parameters
	int subIt;			// Number of iterations
	double relax;		// Aitken relaxation parameter
	double subRes;		// Residual
	double subNum;		// Numerator
	double subDen;		// Denominator

	// Public methods
public:

	// I/O
	void writeInfo();						// Write info every tinfo time steps
	void writeLog();						// Write out log info at start
	void writeVTK(bool writeIBM = true);	// Write IBM VTK
	void writeTotalForces();				// Write out forces on bodies
	void writeTips();						// Write out tip positions
	void writeRestart();					// Write restart file
	void readRestart();						// Read restart file

	// Object routines
	void objectKernel();					// Main kernel for objects

	// Private methods
private:

	// IBM routines
	void ibmKernelInterp();					// Interpolate and force calc
	void ibmKernelSpread();					// Force spread and update macro
	void femKernel();						// Do FEM and update IBM positions and velocities
	void recomputeObjectVals();				// Recompute objects support, ds, and epsilon
	void computeEpsilon();					// Compute epsilon

	// Initialisation
	void removeOverlapMarkers();			// Remove overlapping markers
	void initialiseObjects();				// Initialise objects
	void initialDeflect();					// Give FEM body an initial deflection

	// I/O
	void geometryReadIn();					// Geometry read in
	int getBodyIdxFromID(int id);			// Get body index from the ID
	void restartCleanup();					// Clean up files from previous simulations before restart
	void fillEmptyVTK();					// Fill in empty VTK files
	void writeEmptyVTK(string &fname);		// Write empty VTK file
	void deleteTipsAndForces();				// Delete tip and force data that was written after last restart
	void deleteFileOutput(string fname);	// Delete data from output files that was written after last restart
};

#endif // OBJECTS_H
