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

#ifndef IBMBODY_H // IBMBODY_H
#define IBMBODY_H

// Includes
#include "defs.h"
#include "IBMNode.h"

// Forward declarations
class ObjectsClass;
class FEMBodyClass;

// IBM body class
class IBMBodyClass {

	// Friend classes
	friend class ObjectsClass;
	friend class IBMNodeClass;
	friend class FEMBodyClass;
	friend class FEMElementClass;

	// Default constructor and destructor
public:
	IBMBodyClass() {oPtr = NULL; ID = 0; flex = eRigid; bodyType = eCircle; sBody = NULL;};
	~IBMBodyClass() {};

	// Custom constructor for creating one object from vector of all objects
	IBMBodyClass(vector<IBMNodeClass> &iNode);

	// Custom constructor for building circle
	IBMBodyClass(ObjectsClass *objects, int bodyID, const array<double, dims> &pos, double radius);

	// Custom constructor for building filament
	IBMBodyClass(ObjectsClass *objects, int bodyID, const array<double, dims> &pos, const array<double, dims> &geom,
			double angle, string flex, string nElements, string BC, double rho, double E);

	// Private members
private:

	// Pointer to objects
	ObjectsClass *oPtr;

	// Body parameters
	int ID;								// Body ID
	eFlexibleType flex;					// Flexible or rigid
	eBodyType bodyType;					// Type of body

	// Vector of IBM nodes
	vector<IBMNodeClass*> node;			// Vector of IBM nodes

	// Pointer to FEM body
	FEMBodyClass *sBody;				// Pointer to structural solver
};

#endif // IBMBODY_H
