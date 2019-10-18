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

// Includes
#include "../inc/FEMElement.h"
#include "../inc/FEMBody.h"
#include "../inc/Grid.h"
#include "../inc/Objects.h"
#include "../inc/Utils.h"

// Construct load vector
void FEMElementClass::loadVector() {

	// Declare vectors
	array<array<double, dims>, dims> Tsub;
	array<double, dims> F;
	array<double, elementDOFs> RGlobal;

	// Get weight vector
	array<double, dims> weight = {rho * A * gravityX, rho * A * gravityY};

	// Get force scaling parameter
	double forceScale = fPtr->iPtr->oPtr->gPtr->Dm / SQ(fPtr->iPtr->oPtr->gPtr->Dt);

	// Loop through all force mapping nodes
	for (size_t n = 0; n < forceMap.size(); n++) {

		// Get node and integration ranges
		IBMNodeClass *node = fPtr->iPtr->node[forceMap[n].nodeID];
		double a = forceMap[n].zeta1;
		double b = forceMap[n].zeta2;

		// Get subset of transpose matrix
		Tsub = {{{T[0][0], T[0][1]},
				 {T[1][0], T[1][1]}}};

		// Convert force to local coordinates
		F = Tsub * (((-node->epsilon * 1.0 * forceScale) * node->force) + weight);

		// Get the nodal values by integrating over range of IB point
		R[0] = F[0] * 0.5 * L * (0.5 * b - 0.5 * a + 0.25 * SQ(a) - 0.25 * SQ(b));
		R[1] = F[1] * 0.5 * L * (0.5 * b - 0.5 * a - SQ(a) * SQ(a) / 16.0 + SQ(b) * SQ(b) / 16.0 + 3.0 * SQ(a) / 8.0 - 3.0 * SQ(b) / 8.0);
		R[2] = F[1] * 0.5 * L * (L * (-SQ(a) * SQ(a) + SQ(b) * SQ(b)) / 32.0 - L * (-TH(a) + TH(b)) / 24.0 - L * (-SQ(a) + SQ(b)) / 16.0 + L * (b - a) / 8.0);
		R[3] = F[0] * 0.5 * L * (-0.25 * SQ(a) + 0.25 * SQ(b) + 0.5 * b - 0.5 * a);
		R[4] = F[1] * 0.5 * L * (0.5 * b - 0.5 * a + SQ(a) * SQ(a) / 16.0 - SQ(b) * SQ(b) / 16.0 - 3.0 * SQ(a) / 8.0 + 3.0 * SQ(b) / 8.0);
		R[5] = F[1] * 0.5 * L * (L * (-SQ(a) * SQ(a) + SQ(b) * SQ(b)) / 32.0 + L * (-TH(a) + TH(b)) / 24.0 - L * (-SQ(a) + SQ(b)) / 16.0 - L * (b - a) / 8.0);

		// Get element load vector
		RGlobal = Utils::Transpose(T) * R;

		// Assemble into global vector
		assembleGlobalMat(RGlobal, fPtr->R);
	}
}

// Construct internal force vector
void FEMElementClass::forceVector() {

	// Calculate the local displacements
	double u = (SQ(L) - SQ(L0)) / (L + L0);
	double theta1 = Utils::shiftAngle(node[0]->angle - angle);
	double theta2 = Utils::shiftAngle(node[1]->angle - angle);

	// Calculate internal forces for beam
	double F0 = (E * A / L0) * u;
	double M1 = (2 * E * I / L0) * (2.0 * theta1 + theta2);
	double M2 = (2 * E * I / L0) * (theta1 + 2.0 * theta2);

	// Set the internal local nodal forces
	F[0] = -F0;
	F[1] = (1.0 / L0) * (M1 + M2);
	F[2] = M1;
	F[3] = F0;
	F[4] = -(1.0 / L0) * (M1 + M2);
	F[5] = M2;

	// Get element internal forces
	array<double, elementDOFs> FGlobal = Utils::Transpose(T) * F;

	// Assemble into global vector
	assembleGlobalMat(FGlobal, fPtr->F);
}

// Construct mass matrix
void FEMElementClass::massMatrix() {

	// Get linear stiffness matrix in global coordinates
	array<array<double, elementDOFs>, elementDOFs> MGlobal = Utils::Transpose(T) * M * T;

	// Assemble into global matrix
	assembleGlobalMat(MGlobal, fPtr->M);
}

// Construct stiffness matrix
void FEMElementClass::stiffMatrix() {

	// Get linear stiffness matrix in global coordinates
	array<array<double, elementDOFs>, elementDOFs> KGlobal = Utils::Transpose(T) * K_L * T;

	// Internal forces
	double F0 = -F[0];
	double V0 = F[4];

	// Construct upper half of local stiffness matrix for single element
	K_NL[0][1] = -V0 / L0;
	K_NL[0][4] = V0 / L0;
	K_NL[1][0] = -V0 / L0;
	K_NL[1][1] = F0 / L0;
	K_NL[1][3] = V0 / L0;
	K_NL[1][4] = -F0 / L0;
	K_NL[3][1] = V0 / L0;
	K_NL[3][4] = -V0 / L0;
	K_NL[4][0] = V0 / L0;
	K_NL[4][1] = -F0 / L0;
	K_NL[4][3] = -V0 / L0;
	K_NL[4][4] = F0 / L0;

	// Multiply by transformation matrices to get global matrix for single element
	KGlobal = KGlobal + Utils::Transpose(T) * K_NL * T;

	// Assemble into global matrix
	assembleGlobalMat(KGlobal, fPtr->K);
}

// Multiply by shape functions
array<double, dims> FEMElementClass::shapeFuns(const array<double, elementDOFs> &vec, double zeta) {

	// Results vector
	array<double, dims> resVec;

	// Use shape functions to calculate values
	double N0 = 1.0 - (zeta + 1.0) / 2.0;
	double N1 = 1.0 - 3.0 * SQ((zeta + 1.0) / 2.0) + 2.0 * TH((zeta + 1.0) / 2.0);
	double N2 = ((zeta + 1.0) / 2.0 - 2.0 * SQ((zeta + 1.0) / 2.0) + TH((zeta + 1.0) / 2.0)) * L;
	double N3 = (zeta + 1.0) / 2.0;
	double N4 = 3.0 * SQ((zeta + 1.0) / 2.0) - 2.0 * TH((zeta + 1.0) / 2.0);
	double N5 = (-SQ((zeta + 1.0) / 2.0) + TH((zeta + 1.0) / 2.0)) * L;

	// Calculate values using shape functions
	resVec[eX] = vec[0] * N0 + vec[3] * N3;
	resVec[eY] = vec[1] * N1 + vec[2] * N2 + vec[4] * N4 + vec[5] * N5;

	// Return
	return resVec;
}

// Assemble into global vector
void FEMElementClass::assembleGlobalMat(const array<double, elementDOFs> &localVec, vector<double> &globalVec) {

	// Loop through and set
	for (size_t i = 0; i < localVec.size(); i++)
		globalVec[DOFs[i]] += localVec[i];
}

// Assemble into global matrix
void FEMElementClass::assembleGlobalMat(const array<array<double, elementDOFs>, elementDOFs> &localVec, vector<double> &globalVec) {

	// Get dimensions
	int dim = fPtr->bodyDOFs;

	// Get rows and cols
	size_t rows = localVec.size();
	size_t cols = localVec[0].size();

	// Now loop through and set
	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < cols; j++) {
			globalVec[DOFs[i] * dim + DOFs[j]] += localVec[i][j];
		}
	}
}

// Assemble into global vector
array<double, elementDOFs> FEMElementClass::disassembleGlobalMat(const vector<double> &globalVec) {

	// Declare vector
	array<double, elementDOFs> localVec;

	// Loop through and set
	for (size_t i = 0; i < localVec.size(); i++)
		localVec[i] = globalVec[DOFs[i]];

	// Return
	return localVec;
}

// Set local mass and stiffness matrices
void FEMElementClass::setLocalMatrices() {

	// Size the local matrices
	M.fill({0.0});
	K_L.fill({0.0});
	K_NL.fill({0.0});
	R.fill(0.0);
	F.fill(0.0);

	// Coefficients
	double C1 = rho * A * L0 / 420.0;

	// Set mass matrix first
	M[0][0] = C1 * 140.0;
	M[0][3] = C1 * 70.0;
	M[1][1] = C1 * 156.0;
	M[1][2] = C1 * 22.0 * L0;
	M[1][4] = C1 * 54;
	M[1][5] = C1 * (-13.0 * L0);
	M[2][2] = C1 * 4.0 * SQ(L0);
	M[2][4] = C1 * 13.0 * L0;
	M[2][5] = C1 * (-3.0 * SQ(L0));
	M[3][3] = C1 * 140.0;
	M[4][4] = C1 * 156.0;
	M[4][5] = C1 * (-22.0 * L0);
	M[5][5] = C1 * 4.0 * SQ(L0);

	// Now set stiffness matrix
	K_L[0][0] = E * A / L0;
	K_L[0][3] = -E * A / L0;
	K_L[1][1] = 12.0 * E * I / TH(L0);
	K_L[1][2] = 6.0 * E * I / SQ(L0);
	K_L[1][4] = -12.0 * E * I / TH(L0);
	K_L[1][5] = 6.0 * E * I / SQ(L0);
	K_L[2][2] = 4.0 * E * I / L0;
	K_L[2][4] = -6.0 * E * I / SQ(L0);
	K_L[2][5] = 2.0 * E * I / L0;
	K_L[3][3] = E * A / L0;
	K_L[4][4] = 12.0 * E * I / TH(L0);
	K_L[4][5] = -6.0 * E * I / SQ(L0);
	K_L[5][5] = 4.0 * E * I / L0;

	// Copy to the lower half (symmetrical matrix)
	for (size_t i = 1; i < M.size(); i++) {
		for (size_t j = 0; j < i; j++) {
			M[i][j] = M[j][i];
			K_L[i][j] = K_L[j][i];
		}
	}
}

// Set transformation matrix for element
void FEMElementClass::setElementTransform() {

	// Set to correct values
	T[0][0] = T[1][1] =  T[3][3] = T[4][4] = cos(angle);
	T[0][1] = T[3][4] = sin(angle);
	T[1][0] = T[4][3] = -sin(angle);
	T[2][2] = T[5][5] =  1.0;
}

// Custom constructor for building elements
FEMElementClass::FEMElementClass(FEMBodyClass *fBody, int i, const array<double, dims> &geom, double angleRad, double length, double den, double youngMod) {

	// Set pointer to fBody
	fPtr = fBody;

	// Set values
	L0 = length;
	L = length;
	angle = angleRad;
	A = fPtr->iPtr->oPtr->gPtr->Dx * geom[1];
	I = fPtr->iPtr->oPtr->gPtr->Dx * TH(geom[1]) / 12.0;
	E = youngMod;
	rho = den;

	// Vector of node indices
	array<int, elementNodes> nodeIdx;

	// Get nodal indices
	nodeIdx[0] = i;
	nodeIdx[1] = i + 1;

	// Set pointers and global DOFs
	for (int i = 0; i < elementNodes; i++) {

		// Insert pointer to node
		node[i] = &(fPtr->node[nodeIdx[i]]);

		// Get DOFs
		for (int j = 0; j < nodeDOFs; j++)
			DOFs[i * nodeDOFs + j] = nodeIdx[i] * nodeDOFs + j;
	}

	// Set transformation matrix
	T.fill({0.0});
	setElementTransform();

	// Set local mass and stiffness matrices
	setLocalMatrices();
}
