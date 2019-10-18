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
#include "../inc/IBMNode.h"
#include "../inc/Grid.h"
#include "../inc/Objects.h"
#include "../inc/Utils.h"

// Interpolate
void IBMNodeClass::interpolate() {

	// Get pointer to grid
	GridClass *gPtr = iPtr->oPtr->gPtr;

	// Set to zero
	interpRho = 0.0;
	fill(interpMom.begin(), interpMom.end(), 0.0);

	// Loop through support points and interpolate
	for (size_t s = 0; s < suppCount; s++) {

		// Get ID
		int id = supp[s].idx * Ny + supp[s].jdx;

		// Interpolate density
		interpRho += gPtr->rho[id] * supp[s].diracVal * 1.0 * 1.0;

		// Interpolate momentum
		for (int d = 0; d < dims; d++)
			interpMom[d] += gPtr->rho[id] * gPtr->u[id * dims + d] * supp[s].diracVal * 1.0 * 1.0;
	}
}

// Force calculation
void IBMNodeClass::forceCalc() {

	// Double get vel scaling
	double velScale = iPtr->oPtr->gPtr->Dt / iPtr->oPtr->gPtr->Dx;

	// Loop through dimensions (divided by timestep = 1)
	force = 2.0 * (velScale * interpRho * vel  - interpMom);
}

// Force spread
void IBMNodeClass::spread() {

	// Get pointer to grid
	GridClass *gPtr = iPtr->oPtr->gPtr;

	// Loop through support
	for (size_t s = 0; s < suppCount; s++) {

		// Get ID
		int id = supp[s].idx * Ny + supp[s].jdx;

		// Get forces
		double Fx = force[eX] * epsilon * ds * 1.0 * supp[s].diracVal;
		double Fy = force[eY] * epsilon * ds * 1.0 * supp[s].diracVal;

		// Do either atomic (for performance) or ordered (for repeatability) reduction
#ifdef ORDERED
#pragma omp ordered
		{
			// Ordered write so is repeatable for different runs
			gPtr->force_ibm[id * dims + eX] += Fx;
			gPtr->force_ibm[id * dims + eY] += Fy;
		}
#else
		// Atomic operation for concurrent writes
#pragma omp atomic update
		gPtr->force_ibm[id * dims + eX] += Fx;

		// Atomic operation for concurrent writes
#pragma omp atomic update
		gPtr->force_ibm[id * dims + eY] += Fy;
#endif
	}
}

// Update macroscopic
void IBMNodeClass::updateMacroscopic() {

	// Get pointer to grid
	GridClass *gPtr = iPtr->oPtr->gPtr;

	// Loop through support points
	for (size_t s = 0; s < suppCount; s++) {

		// Reset temp values
		double rhoTmp = 0.0;
		double uTmp = 0.0;
		double vTmp = 0.0;

		// Get ID
		int id = supp[s].idx * Ny + supp[s].jdx;

		// Sum to find rho and momentum
		for (int v = 0; v < nVels; v++) {
			rhoTmp += gPtr->f[id * nVels + v];
			uTmp += gPtr->c[v * dims + eX] * gPtr->f[id * nVels + v];
			vTmp += gPtr->c[v * dims + eY] * gPtr->f[id * nVels + v];
		}

		// Add forces and divide by rho
		uTmp = (uTmp + 0.5 * (gPtr->force_xy[id * dims + eX] + gPtr->force_ibm[id * dims + eX])) / rhoTmp;
		vTmp = (vTmp + 0.5 * (gPtr->force_xy[id * dims + eY] + gPtr->force_ibm[id * dims + eY])) / rhoTmp;

		// Atomic operation for concurrent writes
#pragma omp atomic write
		gPtr->rho[id] = rhoTmp;

		// Atomic operation for concurrent writes
#pragma omp atomic write
		gPtr->u[id * dims + eX] = uTmp;

		// Atomic operation for concurrent writes
#pragma omp atomic write
		gPtr->u[id * dims + eY] = vTmp;
	}
}

// Find support
void IBMNodeClass::findSupport() {

	// Clear the current support points
	supp.fill(IBMSupportClass(0, 0, 0.0));

	// Stencil width
	double stencilWidth = 1.5;

	// Get lattice spacing
	double Dx = iPtr->oPtr->gPtr->Dx;

	// Get closest lattice sites to marker
	int inear = static_cast<int>(round(pos[eX] / Dx));
	int jnear = static_cast<int>(round(pos[eY] / Dx));

	// Loop through x
	suppCount = 0;
	for (int i = inear - 2; i <= inear + 2; i++) {

		// Get distance in x
		double distX = fabs(pos[eX] / Dx - i);

		// Loop through y
		for (int j = jnear - 2; j <= jnear + 2; j++) {

			// Get distance in y
			double distY = fabs(pos[eY] / Dx - j);

			// Check distance and if it is within grid
			if (distX < stencilWidth && distY < stencilWidth && i >= 0 && i <= Nx-1 && j >= 0 && j <= Ny-1) {

				// Check if buffer size is big enough
				if (suppCount == suppSize)
					ERROR("Support buffer size is not big enough for number of support points...exiting");

				// Add to support
				supp[suppCount++] = IBMSupportClass(i, j, Utils::diracDelta(distX) * Utils::diracDelta(distY));
			}
		}
	}
}

// Compute spacing between IBM nodes
void IBMNodeClass::computeDs() {

	// Set ds to large value
	double currentDs = 10.0;

	// Loop through all markers within own body
	for (size_t n = 0; n < iPtr->node.size(); n++) {

		// Don't check self
		if (ID != iPtr->node[n]->ID) {

			// Get distance between nodes
			double mag = sqrt((pos - iPtr->node[n]->pos) * (pos - iPtr->node[n]->pos)) / iPtr->oPtr->gPtr->Dx;

			// Check if min found so far
			if (mag < currentDs)
				currentDs = mag;
		}
	}

	// Set ds value
	ds = currentDs;
}

// Custom constructor for building node
IBMNodeClass::IBMNodeClass(IBMBodyClass *iPointer, int nodeID, const array<double, dims> &position) {

	// Set pointer
	iPtr = iPointer;

	// Set ID
	ID = nodeID;

	// Set positions
	pos = position;

	// Set velocities and forces to zero
	vel.fill(0.0);
	force.fill(0.0);

	// These values will be set later
	ds = 0.0;
	epsilon = 0.0;

	// Set interpolated values
	interpMom.fill(0.0);
	interpRho = 0.0;

	// Set number of support points
	suppCount = 0;
}
