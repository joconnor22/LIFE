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

#ifndef IBMNODE_H // IBMNODE_H
#define IBMNODE_H

// Includes
#include "defs.h"
#include "IBMSupport.h"

// Forward declarations
class IBMBodyClass;

// IBM node class
class IBMNodeClass {

	// Friend classes
	friend class ObjectsClass;
	friend class FEMElementClass;
	friend class FEMBodyClass;
	friend class IBMBodyClass;

	// Default constructor and destructor
public:
	IBMNodeClass() {iPtr = NULL; ID = 0; ds = 0.0; epsilon = 0.0; interpRho = 0.0; suppCount = 0;};
	~IBMNodeClass() {};

	// Custom constructor for building node
	IBMNodeClass(IBMBodyClass *iPtr, int nodeID, const array<double, dims> &pos);

	// Private members
private:

	// Pointer to body
	IBMBodyClass *iPtr;						// Pointer to owning IBM body

	// ID
	int ID;									// Global node ID

	// Positions, velocities and forces
	array<double, dims> pos;				// Position of node
	array<double, dims> vel;				// Velocity at node
	array<double, dims> force;				// Force on IBM marker

	// Spacing and force multiplier
	double ds; 								// Spacing between points (lattice units)
	double epsilon;							// Multiplier for IBM

	// Interpolated values
	array<double, dims> interpMom;			// Interpolated momentum
	double interpRho;						// Interpolated density

	// Vector of support points
	size_t suppCount;						// Number of support points for this IBM marker
	array<IBMSupportClass, suppSize> supp;	// Vector of support points

	// Private methods
private:

	// IBM methods
	void findSupport();						// Find support
	void computeDs();						// Compute spacing between IBM nodes
	void interpolate();						// Interpolation
	void forceCalc();						// Force calculation
	void spread();							// Spread force back
	void updateMacroscopic();				// Update macroscopic values at support points
};

#endif // IBMNODE_H
