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
#include "../inc/FEMNode.h"

// Custom constructor for building nodes
FEMNodeClass::FEMNodeClass(int ID, const array<double, dims> &position, double angleRad) {

	// Set values
	pos0 = position;
	pos = position;
	angle0 = angleRad;
	angle = angleRad;

	// Set the global DOFs
	for (int i = 0; i < nodeDOFs; i++)
		DOFs[i] = ID * nodeDOFs + i;
}
