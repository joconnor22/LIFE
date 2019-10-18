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
#include "../inc/FEMBody.h"
#include "../inc/Grid.h"
#include "../inc/Objects.h"
#include "../inc/Utils.h"

// Dynamic FEM routine
void FEMBodyClass::dynamicFEM() {

	// Reset back to start of time step
	U = U_n;
	Udot = Udot_n;
	Udotdot = Udotdot_n;

	// Update FEM elements
	updateFEMValues();

	// Construct load vector as invariant during Newton-Raphson iterations
	constructRVector();

	// While loop parameters
	double TOL = 1e-10;
	double MAXIT = 20;

	// Set while counter to zero
	itNR = 0;

	// While loop for FEM solver
	do {

		// Solve and iterate over the system
		newtonRaphsonDynamic();

		// Check residual
		resNR = checkNRConvergence();

		// Increment counter
		itNR++;

	} while (resNR > TOL && itNR < MAXIT);

	// Compute new velocities and accelerations
	finishNewmark();

	// Update IBM nodes
	updateIBMValues();

	// Get subiteration residual for this body (will be summed later on)
	subResidual();
}

// Newton raphson iterator
void FEMBodyClass::newtonRaphsonDynamic() {

	// Build global matrices
	buildGlobalMatrices();

	// Apply Newmark scheme (using Newmark coefficients)
	setNewmark();

	// Solve linear system using LAPACK library
	delU = Utils::solveLAPACK(K, F, bcDOFs);

	// Add deltaU to U
	U = U + delU;

	// Update FEM positions
	updateFEMValues();
}

// Build global matrices
void FEMBodyClass::buildGlobalMatrices() {

	// Set matrices to zero
	fill(M.begin(), M.end(), 0.0);
	fill(K.begin(), K.end(), 0.0);
	fill(F.begin(), F.end(), 0.0);

	// Loop through and build global matrices
	for (size_t el = 0; el < element.size(); el++) {

		// Build force vector
		element[el].forceVector();

		// Build mass matrix
		element[el].massMatrix();

		// Build stiffness matrix
		element[el].stiffMatrix();
	}
}

// Set Newmark
void FEMBodyClass::setNewmark() {

	// Newmark-beta method for time integration
	double Dt = iPtr->oPtr->gPtr->Dt;
	double a0, a2, a3;
	a0 = 1.0 / (alpha * SQ(Dt));
	a2 = 1.0 / (alpha * Dt);
	a3 = 1.0 / (2.0 * alpha) - 1.0;

	// Calculate effective load vector
	F = R - F + Utils::MatMultiply(M, a0 * (U_n - U) + a2 * Udot + a3 * Udotdot);

	// Calculate effective stiffness matrix
	K = K + a0 * M;
}

// Finish Newmark
void FEMBodyClass::finishNewmark() {

	// Get timestep
	double Dt = iPtr->oPtr->gPtr->Dt;

	// Newmark coefficients
	double a6 = 1.0 / (alpha * SQ(Dt));
	double a7 = -1.0 / (alpha * Dt);
	double a8 = -(1.0 / (2.0 * alpha) - 1.0);
	double a9 = Dt * (1.0 - delta);
	double a10 = delta * Dt;

	// Update velocities and accelerations
	Udotdot = a6 * (U - U_n) + a7 * Udot_n + a8 * Udotdot_n;
	Udot = Udot_n + a9 * Udotdot_n + a10 * Udotdot;
}

// Update IBM markers
void FEMBodyClass::updateIBMValues() {

	// Declare local vectors
	array<double, elementDOFs> dashU;
	array<double, elementDOFs> dashUdot;
	array<double, dims> dashUShape;
	array<double, dims> dashUdotShape;
	array<array<double, dims>, dims> Tsub;

	// Loop through posMap points
	for (size_t i = 0; i < posMap.size(); i++) {

		// Get pointer to element and zeta value
		FEMElementClass *el = &(element[posMap[i].elID]);
		double zeta = posMap[i].zeta;

		// Disassemble positions, velocities and accelerations
		dashU = el->disassembleGlobalMat(U);
		dashUdot = el->disassembleGlobalMat(Udot);

		// Convert the displacement into local coordinates (shift angle as well)
		for (size_t n = 0; n < el->node.size(); n++) {
			dashU[n*nodeDOFs+eX] += el->node[n]->pos0[eX] - el->node[0]->pos[eX];
			dashU[n*nodeDOFs+eY] += el->node[n]->pos0[eY] - el->node[0]->pos[eY];
			dashU[n*nodeDOFs+(nodeDOFs-1)] += el->node[n]->angle0 - el->angle;
			dashU[n*nodeDOFs+(nodeDOFs-1)] = Utils::shiftAngle(dashU[n*nodeDOFs+(nodeDOFs-1)]);
		}

		// Get element values in local coordinates
		dashU = el->T * dashU;
		dashUdot = el->T * dashUdot;

		// Multiply by shape functions
		dashUShape = el->shapeFuns(dashU, zeta);
		dashUdotShape = el->shapeFuns(dashUdot, zeta);

		// Get subset of transformation matrix
		Tsub = {{{el->T[0][0], el->T[0][1]},
				 {el->T[1][0], el->T[1][1]}}};

		// Shift back to global coordinates
		dashUShape = Utils::Transpose(Tsub) * dashUShape;
		dashUdotShape = Utils::Transpose(Tsub) * dashUdotShape;

		// Set the IBM nodes
		iPtr->node[i]->pos = el->node[0]->pos + dashUShape;
		iPtr->node[i]->vel = dashUdotShape;
	}
}

// Update FEM values
void FEMBodyClass::updateFEMValues() {

	// Set the new positions
	for (size_t n = 0; n < node.size(); n++) {

		// Positions
		for (int d = 0; d < dims; d++)
			node[n].pos[d] = node[n].pos0[d] + U[node[n].DOFs[d]];

		// Set angle
		node[n].angle = node[n].angle0 + U[node[n].DOFs[dims]];
	}

	// Set the new angles and lengths of the elements
	array<double, dims> elVector;
	for (size_t el = 0; el < element.size(); el++) {

		// Set the new angles and lengths of the elements
		elVector = element[el].node[1]->pos - element[el].node[0]->pos;
		element[el].angle = atan2(elVector[eY], elVector[eX]);
		element[el].L = sqrt(elVector * elVector);

		// Set new transformation matrix
		element[el].setElementTransform();
	}
}

// Get load vector
void FEMBodyClass::constructRVector() {

	// Set R to zero
	fill(R.begin(), R.end(), 0.0);

	// Loop through all FEM elements
	for (size_t el = 0; el < element.size(); el++)
		element[el].loadVector();
}

// Check convergence of Newton Raphson iterator
inline double FEMBodyClass::checkNRConvergence () {

	// Get the norm of delU
	return sqrt(delU * delU) / (ref_L * sqrt(static_cast<double>(delU.size())));
}

// Do a sum reduction to get the subiteration residual
void FEMBodyClass::subResidual() {

	// Get the residual from this time step and reassign old time step value
	R_km1 = R_k;
	R_k = U - U_km1;

	// Get residual for this body
	subRes = R_k * R_k;

	// Get numerator and denominator for calculating relaxation factor
	subNum = R_km1 * (R_k - R_km1);
	subDen = (R_k - R_km1) * (R_k - R_km1);
}

// predictor
void FEMBodyClass::predictor() {

	// Extrapolate the values
	if (iPtr->oPtr->gPtr->t > 2) {

		// Do 2rd order extrapolation
		U = 2.5 * U_n - 2.0 * U_nm1 + 0.5 * U_nm2;
	}
	else if (iPtr->oPtr->gPtr->t == 2) {

		// Do 1st order extrapolation
		U = 2.0 * U_n - U_nm1;
	}
	else if (iPtr->oPtr->gPtr->t == 1) {

		// Do zeroth order extrapolation
		U = U_n;
	}

	// Update FEM elements
	updateFEMValues();

	// Update the velocity
	finishNewmark();

	// Update IBM values
	updateIBMValues();

	// Set U_km1
	U_km1.swap(U);
}

// Get node mappings between FEM and IBM grids
void FEMBodyClass::computeNodeMapping(int nIBMNodes, int nFEMNodes) {

	// Size the map from FEM to IBM (1 for each IBM node)
	posMap.resize(nIBMNodes);

	// Set first node first
	posMap[0].elID = 0;
	posMap[0].zeta = -1.0;

	// Now loop through and set values
	for (int i = 1; i < nIBMNodes; i++) {

		// Check if remainder is zero
		if (fmod(i * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)), 1.0) == 0.0) {
			posMap[i].elID = static_cast<int>(std::floor(i * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)))) - 1;
			posMap[i].zeta = 1.0;
		}
		else {
			posMap[i].elID = static_cast<int>(std::floor(i * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0))));
			posMap[i].zeta = fmod(i * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)), 1.0) * 2.0 - 1.0;
		}
	}

	// Loop through elements
	double node1, node2;
	for (size_t el = 0; el < element.size(); el++) {

		// Loop through all nodes and scale range to local coordinates for element
		for (int node = 0; node < nIBMNodes; node++) {
			node1 = -1 + 2.0 * ((static_cast<double>(node) - 0.5) * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)) - static_cast<double>(el)) / 1.0;
			node2 = -1 + 2.0 * ((static_cast<double>(node) + 0.5) * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)) - static_cast<double>(el)) / 1.0;

			// Check if any points lie within element coordinate range
			if ((node1 > -1.0 && node1 < 1.0) || (node2 > -1.0 && node2 < 1.0)) {

				// Sort out end nodes where one point lie outside element
				if (node1 < -1.0)
					node1 = -1.0;
				if (node2 > 1.0)
					node2 = 1.0;

				// Call constructor for chile IB point
				element[el].forceMap.emplace_back(node, node1, node2);
			}
		}
	}
}

// Reset the start of time step values
void FEMBodyClass::resetValues() {

	// Reset start of time step values
	U_nm2.swap(U_nm1);
	U_nm1.swap(U_n);
	U_n.swap(U);
	Udot_n.swap(Udot);
	Udotdot_n.swap(Udotdot);
}

// Set initial deflection using static FEM
void FEMBodyClass::setInitialDeflection() {

	// Write out
	cout << endl << "Starting static FEM for body " << iPtr->ID << "...";

	// Initial deflection
	double deflect = 0.0;

	// Set initial deflection
#ifdef INITIAL_DEFLECT
	deflect = -INITIAL_DEFLECT;
#endif

	// Get beam properties
	double E = element[0].E;
	double I = element[0].I;

	// Calculate linear force required to give specified deflection
	double linearForce = 3.0 * E * I * deflect / SQ(L0);

	// Set load steps and get max deltaForce
	int nSteps = 100;
	double deltaForceMax = linearForce / static_cast<double>(nSteps);

	// Get rotation matrix
	array<array<double, dims>, dims> T = Utils::getRotationMatrix(angle0);

	// Get initial force at tip
	array<double, dims> localForce = {0.0, deltaForceMax};
	array<double, dims> globalForce = T * localForce;

	// Set it to load vector
	R[bodyDOFs-nodeDOFs] = globalForce[0];
	R[bodyDOFs-nodeDOFs+1] = globalForce[1];

	// Declare feedback loop parameters
	array<double, dims> deflectVector;
	double gain = 1e-3;
	double TOL = 1e-10;
	int MAXIT = 10000;
	int it = 0;
	double error;

	// Start feedback loop
	do {

		// Static FEM calculation
		staticFEM();

		// Get deflection vector
		deflectVector =  {U[bodyDOFs-nodeDOFs], U[bodyDOFs-nodeDOFs+1]};
		deflectVector = Utils::Transpose(T) * deflectVector;

		// Get error
		error = (deflect - (deflectVector[1] / L0));

		// Force correction
		double forceCorrect = error * gain;

		// Check force correction is smaller than force interval
		if (fabs(forceCorrect) > fabs(deltaForceMax))
			forceCorrect = Utils::sgn(forceCorrect) * fabs(deltaForceMax);

		// Add to R vector
		localForce[1] += forceCorrect;
		globalForce = T * localForce;

		// Add to R vector
		R[bodyDOFs-nodeDOFs] = globalForce[0];
		R[bodyDOFs-nodeDOFs+1] = globalForce[1];

		// Increment counter
		it++;

		// Check if done max iterations
		if (it >= MAXIT)
			ERROR("max iterations hit! Try changing the gain...exiting");

	} while (fabs(error) > TOL);

	// Write out
	cout << "finished in " << it << " iterations";

	// Reset values to zero
	itNR = 0;
	resNR = 0.0;
	fill(R.begin(), R.end(), 0.0);
	fill(delU.begin(), delU.end(), 0.0);
}

// Newton raphson iterator
void FEMBodyClass::staticFEM() {

	// While loop parameters
	double TOL = 1e-10;
	double MAXIT = 20;

	// Set while counter to zero
	itNR = 0;

	// While loop for FEM solver
	do {

		// Solve and iterate over the system
		newtonRaphsonStatic();

		// Check residual
		resNR = checkNRConvergence();

		// Increment counter
		itNR++;

	} while (resNR > TOL && itNR < MAXIT);

	// Update IBM nodes
	updateIBMValues();
}

// Newton raphson iterator
void FEMBodyClass::newtonRaphsonStatic() {

	// Build global matrices
	buildGlobalMatrices();

	// Solve linear system using LAPACK library
	delU = Utils::solveLAPACK(K, R - F, bcDOFs);

	// Add deltaU to U
	U = U + delU;

	// Update FEM positions
	updateFEMValues();
}

// Custom constructor for building corotational FEM body
FEMBodyClass::FEMBodyClass(IBMBodyClass *iBodyPtr, const array<double, dims> &pos, const array<double, dims> &geom, double angle, string nElementsStr, string BC, double rho, double E) {

	// Set pointer
	iPtr = iBodyPtr;

	// Set to initial value
	itNR = 0;
	resNR = 0.0;

	// Set sub residaul, numerator and denominator to initial value
	subRes = 0.0;
	subNum = 0.0;
	subDen = 0.0;

	// Get number of DOFs required for BC
	if (BC == "CLAMPED")
		bcDOFs = 3;
	else if (BC == "SUPPORTED")
		bcDOFs = 2;

	// Unpack geometry
	angle0 = angle;
	L0 = geom[0];

	// Get number of elements
	int numElements;
	if (nElementsStr == "CONFORMING")
		numElements = static_cast<int>(floor(L0 / iPtr->oPtr->gPtr->Dx));
	else
		numElements = static_cast<int>(stod(nElementsStr));

	// Get number of nodes
	int numNodes = numElements + 1;

	// Get number of DOFs in body
	bodyDOFs = numNodes * nodeDOFs;

	// Get critical time step for this body
	tCrit = (L0 / numElements) / sqrt(E / rho);

	// Get rotation matrix
	array<array<double, dims>, dims> T = Utils::getRotationMatrix(angle);

	// Position vector for marker
	array<double, dims> position;

	// Loop through and build nodes
	for (int i = 0; i < numNodes; i++) {

		// Get position of this marker
		position[eX] = i * L0 / numElements;
		position[eY] = 0.0;

		// Rotate
		position = T * position;

		// Add the start point
		position = pos + position;

		// Call node constructor
		node.emplace_back(i, position, angle);
	}

	// Loop through and build elements
	for (int i = 0; i < numElements; i++)
		element.emplace_back(this, i, geom, angle, L0 / numElements, rho, E);

	// Get number of IBM and FEM nodes
	int nIBMNodes = static_cast<int>(iPtr->node.size());
	int nFEMNodes = static_cast<int>(node.size());

	// Compute IBM-FEM conforming parameters
	computeNodeMapping(nIBMNodes, nFEMNodes);

	// Size the matrices
	M.resize(bodyDOFs * bodyDOFs, 0.0);
	K.resize(bodyDOFs * bodyDOFs, 0.0);
	R.resize(bodyDOFs, 0.0);
	F.resize(bodyDOFs, 0.0);
	U.resize(bodyDOFs, 0.0);
	delU.resize(bodyDOFs, 0.0);
	Udot.resize(bodyDOFs, 0.0);
	Udotdot.resize(bodyDOFs, 0.0);
	U_n.resize(bodyDOFs, 0.0);
	U_nm1.resize(bodyDOFs, 0.0);
	U_nm2.resize(bodyDOFs, 0.0);
	Udot_n.resize(bodyDOFs, 0.0);
	Udotdot_n.resize(bodyDOFs, 0.0);
	U_km1.resize(bodyDOFs, 0.0);
	R_k.resize(bodyDOFs, 0.0);
	R_km1.resize(bodyDOFs, 0.0);
}
