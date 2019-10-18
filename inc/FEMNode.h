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

#ifndef FEMNODE_H // FEMNODE_H
#define FEMNODE_H

// Includes
#include "defs.h"

// FEM node class
class FEMNodeClass {

	// Friend classes
	friend class FEMElementClass;
	friend class FEMBodyClass;
	friend class ObjectsClass;

	// Default constructor and destructor
public:
	FEMNodeClass() {angle0 = 0; angle = 0;};
	~FEMNodeClass() {};

	// Custom constructor for building nodes
	FEMNodeClass(int ID, const array<double, dims> &position, double angleRad);

	// Private members
private:

	// Values
	array<double, dims> pos0;			// Initial position of FEM node
	array<double, dims> pos;			// Current position of FEM node
	double angle0;						// Initial angles of FEM node
	double angle;						// Current angles of FEM node
	array<int, nodeDOFs> DOFs;			// Global DOFs for this node
};

#endif // FEMNODE_H
