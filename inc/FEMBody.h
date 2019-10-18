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

#ifndef FEMBODY_H // FEMBODY_H
#define FEMBODY_H

// Includes
#include "defs.h"
#include "FEMNode.h"
#include "FEMElement.h"

// Forward declarations
class IBMBodyClass;

// FEM body class
class FEMBodyClass {

	// Friend classes
	friend class FEMElementClass;
	friend class ObjectsClass;

	// Nested class for mapping to IBM positions
	class posMapClass {

		// Set friend class
		friend class FEMBodyClass;

		// Default constructor and destructor
	public:
		posMapClass() {elID = 0; zeta = 0.0;};
		~posMapClass() {};

		// Members
	private:
		int elID;		// Element ID
		double zeta;	// Node position in local coordinates
	};

	// Default constructor and destructor
public:
	FEMBodyClass() {iPtr = NULL; L0 = 0.0; angle0 = 0.0; bodyDOFs = 0; bcDOFs = 0; itNR = 0; resNR = 0.0; tCrit = 0.0; subRes = 0.0; subNum = 0.0; subDen = 0.0;};
	~FEMBodyClass() {};

	// Custom constructor for building corotational FEM body
	FEMBodyClass(IBMBodyClass *iBodyPtr, const array<double, dims> &pos, const array<double, dims> &geom, double angle, string nElementsStr, string BC, double rho, double E);

	// Private members
private:

	// Pointer to IBM body
	IBMBodyClass *iPtr;

	// Geometry
	double L0;					// Initial length
	double angle0;				// Initial angle

	// System values
	int bodyDOFs;				// DOFs for whole body
	int bcDOFs;					// Number of DOFs removed when applying BCs
	int itNR;					// Number of iterations for Newton-Raphson solver
	double resNR;				// Residual Newton-Raphson solver reached
	double tCrit;				// Critical time step

	// System matrices
	vector<double> M;						// Mass matrix
	vector<double> K;						// Stiffness matrix
	vector<double> R;						// Load vector
	vector<double> F;						// Vector of internal forces
	vector<double> U;						// Vector of displacements
	vector<double> delU;					// Vector of incremental displacements
	vector<double> Udot;					// Vector velocities
	vector<double> Udotdot;					// Vector of accelerations
	vector<double> U_n;						// Vector of displacements at start of timestep
	vector<double> U_nm1;					// Vector of displacements at last time step
	vector<double> U_nm2;					// Vector of displacements at two time steps ago
	vector<double> Udot_n;					// Vector velocities at start of timestep
	vector<double> Udotdot_n;				// Vector of accelerations at start of timestep
	vector<double> U_km1;					// Vector of displacements at last iteration

	// Vector of nodes and elements
	vector<FEMNodeClass> node;				// Vector of nodes
	vector<FEMElementClass> element;		// Vector of elements

	// Vector of parent elements for each IBM node
	vector<posMapClass> posMap;

	// Parameters for Aitken relaxation
	double subRes;				// Residual
	double subNum;				// Numerator
	double subDen;				// Denominator
	vector<double> R_k;			// Vector of nodal residuals
	vector<double> R_km1;		// Vector of nodal residuals at last iteration

	// Private methods
private:

	// FEM methods
	void dynamicFEM();				// Dynamic FEM routine
	void newtonRaphsonDynamic();	// Newton Raphson iterator
	void constructRVector();		// Get load vector
	void buildGlobalMatrices();		// Build global matrices
	void setNewmark();				// Set Newmark scheme
	void finishNewmark();			// Finish Newmark
	void updateIBMValues();			// Update IBM values
	void updateFEMValues();			// Update FEM parameters
	double checkNRConvergence();	// Check convergence of Newton Raphson iterator
	void setInitialDeflection();	// Set initial deflection
	void staticFEM();				// Static FEM routine
	void newtonRaphsonStatic();		// Static Newton Raphson iterator
	void predictor();				// Structural predictor at start of time step
	void subResidual();				// Subiteration residual for this body

	// Helper routines
	void computeNodeMapping(int nIBMNodes, int nFEMNodes);	// Get node mappings between FEM and IBM grids
	void resetValues();										// Reset the start of time step values
};

#endif // FEMBODY_H
