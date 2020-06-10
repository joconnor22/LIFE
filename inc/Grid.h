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

#ifndef GRID_H // GRID_H
#define GRID_H

// Includes
#include "defs.h"

// Forward declarations
class ObjectsClass;

// Grid class
class GridClass {

	// Friend classes
	friend class ObjectsClass;
	friend class IBMNodeClass;

	// Default constructor and destructor
public:
	GridClass();
	~GridClass() {};

	// Public members
public:

	// Grid values
	int t;									// Time step
	int tOffset;							// Offset time (for restarts)

	// Object pointer
	ObjectsClass *oPtr;						// Pointer to objects class

	// Scaling parameters
	double Dx;								// Length scaling
	double Dt;								// Time scaling
	double Dm;								// Mass scaling
	double Drho;							// Density scaling

	// Restart flag
	bool restartFlag;						// Restart flag

	// Big or little endian
	bool bigEndian;							// For binary output

	// Private members
private:

	// Lattice parameters
	double c_s;								// Speed of sound
	array<double, nVels> w;					// Weightings
	array<int, nVels * dims> c;				// Direction vectors
	array<int, nVels> opposite;				// Opposite vectors

	// Grid values
	double tau;								// Relaxation time
	double nu;								// Lattice viscosity

	// Flattened kernel arrays
	vector<double> u;						// Velocity
	vector<double> u_n;						// Velocity (start of timestep)
	vector<double> rho;						// Density
	vector<double> rho_n;					// Density (start of timestep)
	vector<double> force_xy;				// Cartesian force (pressure and gravity)
	vector<double> force_ibm;				// Cartesian force (IBM)
	vector<eLatType> type;					// Lattice type matrix
	vector<double> f;						// Populations
	vector<double> f_n;						// Populations (start of timestep)

	// Boundary conditions
	vector<int> BCVec;						// Vector of site IDs to apply boundary conditions
	vector<double> delU;					// Convective speed through boundary

	// Helper arrays
	vector<double> u_in;					// Inlet velocity profile
	vector<double> rho_in;					// Inlet density profile

	// Timings
	double startTime;							// Start time when clock is called
	double loopTime;							// Average loop time

	// Public methods
public:

	// LBM methods
	void solver();							// Main solver

	// I/O
	void writeInfo();						// Write info every tinfo time steps
	void writeLog();						// Write log data at start
	void writeVTK();						// Write grid VTK
	void writeRestart();					// Write restart file
	void readRestart();						// Read restart file
	void startClock();						// Start the clock for getting MLUPS

	// Private methods
private:

	// LBM methods
	void lbmKernel();															// Main LBM kernel
	void streamCollide(int i, int j, int id);									// Stream and collide in one go (push algorithm)
	double equilibrium(int id, int v);											// Equilibrium function
	double latticeForce(int id, int v);											// Discretise lattice force (BGK only)
	void macroscopic(int id);													// Compute macroscopic quantities
	void applyBCs(int i, int j, int id);										// Apply boundary conditions
	void convectiveSpeed();														// Get convective speed through right boundary
	void convectiveBC(int j, int id);											// Convective BC
	void regularisedBC(int i, int j, int id,
			array<int, dims> &normalVector, eDirectionType normalDirection);	// Regularised BC


	// Initialisation
	void initialiseGrid();									// Initialise grid values

	// Helper routines
	array<int, dims> getNormalVector(int i, int j, eDirectionType &normalDirection);	// Get normal vector for boundary site
	double getRampCoefficient();														// Get inlet ramp coefficient
};

#endif // GRID_H
