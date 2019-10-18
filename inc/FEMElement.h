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

#ifndef FEMELEMENT_H // FEMELEMENT_H
#define FEMELEMENT_H

// Includes
#include "defs.h"

// Forward declarations
class FEMBodyClass;
class FEMNodeClass;

// FEM element class
class FEMElementClass {

	// Friend classes
	friend class FEMBodyClass;
	friend class ObjectsClass;

	// Nested class for mapping forces from IBM
	class forceMapClass {

		// Friend classes
		friend class FEMElementClass;

		// Default constructor and destructor
	public:
		forceMapClass() {nodeID = 0; zeta1 = 0.0; zeta2 = 0.0;};
		~forceMapClass() {};

		// Custom constructor for initialising members
		forceMapClass(int node, double zetaA, double zetaB) {nodeID = node; zeta1 = zetaA; zeta2 = zetaB;};

		// Members
	private:
		int nodeID;			// IBM node ID
		double zeta1;		// Local coordinate of start of IBM force
		double zeta2;		// Local coordinate of end of IBM force
	};

	// Default constructor and destructor
public:
	FEMElementClass() {fPtr = NULL; L0 = 0.0; L = 0.0; angle = 0.0; A = 0.0; I = 0.0; E = 0.0; rho = 0.0;};
	~FEMElementClass() {};

	// Custom constructor for building elements
	FEMElementClass(FEMBodyClass *fBody, int i, const array<double, dims> &geom, double angle, double L, double den, double youngMod);

	// Private members
private:

	// Pointer to fBody
	FEMBodyClass *fPtr;

	// Geometry properties
	double L0;						// Initial length of element
	double L;						// Current length of element
	double angle;					// Orientation of element at start of timestep
	double A;						// Cross-sectional area of element
	double I;						// Second moment areas

	// Structural properties
	double E;						// Youngs modulus
	double rho;						// Material density

	// Transformation matrix
	array<array<double, elementDOFs>, elementDOFs> T;		// Local transformation matrix

	// Local mass and stiffness matrices (which are constant throughout simulation)
	array<array<double, elementDOFs>, elementDOFs> M;		// Mass matrix
	array<array<double, elementDOFs>, elementDOFs> K_L;		// Stiffness matrix (linear)
	array<array<double, elementDOFs>, elementDOFs> K_NL;	// Stiffness matrix (non-linear)
	array<double, elementDOFs> R;							// Load vector
	array<double, elementDOFs> F;							// Vector of internal forces

	// Global DOFs
	array<int, elementDOFs> DOFs;							// Global DOFs for this element

	// Pointers to the nodes belonging to this element
	array<FEMNodeClass*, elementNodes> node;

	// Vector of force mappings from IBM to element
	vector<forceMapClass> forceMap;

	// Private methods
private:

	// FEM routines
	void loadVector();																			// Construct load vector
	void massMatrix();																			// Construct mass matrix
	void stiffMatrix();																			// Construct stiffness matrix
	void forceVector();																			// Construct internal force vector
	array<double, dims> shapeFuns(const array<double, elementDOFs> &vec, double zeta);			// Multiply by shape functions

	// Helper routines
	void setLocalMatrices();																							// Set local mass and stiffness matrices
	void setElementTransform();																							// Set the transformation matrix for the element
	void assembleGlobalMat(const array<double, elementDOFs> &localVec, vector<double> &globalVec);						// Assemble into global vector
	void assembleGlobalMat(const array<array<double, elementDOFs>, elementDOFs> &localVec, vector<double> &globalVec);	// Assemble into global matrix
	array<double, elementDOFs> disassembleGlobalMat(const vector<double> &globalVec);									// Disassemble global vector
};

#endif // FEMELEMENT_H
