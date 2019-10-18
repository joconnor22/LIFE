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
#include "../inc/IBMBody.h"
#include "../inc/FEMBody.h"
#include "../inc/Grid.h"
#include "../inc/Objects.h"
#include "../inc/Utils.h"

// Custom constructor for creating one object from vector of all objects
IBMBodyClass::IBMBodyClass(vector<IBMNodeClass> &iNode) {

	// Set values (some of them don't matter)
	oPtr = iNode[0].iPtr->oPtr;
	ID = 0;
	flex = eFlexible;
	bodyType = eCircle;
	sBody = NULL;

	// Reserve space
	node.reserve(iNode.size());

	// Loop through all markers and copy into this vector
	for (size_t n = 0; n < iNode.size(); n++)
		node.push_back(&iNode[n]);
}

// Custom constructor for building circle
IBMBodyClass::IBMBodyClass(ObjectsClass *objects, int bodyID,  const array<double, dims> &pos, double radius) {

	// Set body values
	oPtr = objects;
	ID = bodyID;
	flex = eRigid;
	bodyType = eCircle;
	sBody = NULL;

	// Get number of markers
	int numNodes = static_cast<int>(floor(2.0 * M_PI * radius / oPtr->gPtr->Dx));

	// Position vector for marker
	array<double, dims> position;

	// Reserve space
	node.reserve(numNodes);

	// Loop through and build body
	for (int i = 0; i < numNodes; i++) {

		// Get position of this marker
		position[eX] = pos[eX] + radius * cos(i * 2.0 * M_PI / numNodes);
		position[eY] = pos[eY] + radius * sin(i * 2.0 * M_PI / numNodes);

		// Call node constructor
		oPtr->iNode.emplace_back(this, oPtr->iNode.size(), position);

		// Insert pointer into node vector
		node.push_back(&(oPtr->iNode.back()));
	}
}

// Custom constructor for building filament
IBMBodyClass::IBMBodyClass(ObjectsClass *objects, int bodyID, const array<double, dims> &pos, const array<double, dims> &geom,
		double angle, string flexStr, string nElements, string BC, double rho, double E) {

	// Set body values
	oPtr = objects;
	ID = bodyID;
	bodyType = eFilament;
	sBody = NULL;

	// Set flex type
	if (flexStr == "FLEXIBLE")
		flex = eFlexible;
	else if (flexStr == "RIGID")
		flex = eRigid;
	else
		ERROR("Not a proper flexible type set...exiting");

	// Unpack geometry and convert angle to radians
	double length = geom[0];
	angle = angle * M_PI / 180.0;

	// Get rotation matrix
	array<array<double, dims>, dims> T = Utils::getRotationMatrix(angle);

	// Build the IBM body
	int numNodes = static_cast<int>(floor(length / oPtr->gPtr->Dx)) + 1;

	// Position vector for marker
	array<double, dims> position;

	// Reserve space
	node.reserve(numNodes);

	// Loop through and build body
	for (int i = 0; i < numNodes; i++) {

		// Get position of this marker
		position[eX] = i * length / (numNodes - 1);
		position[eY] = 0.0;

		// Rotate
		position = T * position;

		// Add the start point
		position = pos + position;

		// Call node constructor
		oPtr->iNode.emplace_back(this, oPtr->iNode.size(), position);

		// Insert pointer into node vector
		node.push_back(&(oPtr->iNode.back()));
	}

	// If flexible then build body
	if (flex == eFlexible)
		sBody = new FEMBodyClass(this, pos, geom, angle, nElements, BC, rho, E);
}
