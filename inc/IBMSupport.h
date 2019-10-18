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

#ifndef IBMSUPPORT_H // IBMSUPPORT_H
#define IBMSUPPORT_H

// IBM support class
class IBMSupportClass {

	// Friend classes
	friend class ObjectsClass;
	friend class IBMNodeClass;

	// Default constructor and destructor
public:
	IBMSupportClass() {idx = 0; jdx = 0; diracVal = 0.0;};
	~IBMSupportClass() {};

	// Custom constructor for building node
	IBMSupportClass(int i, int j, double deltaVal);

	// Private members
private:

	// Values
	int idx;					// i-index
	int jdx;					// j-index
	double diracVal;			// Dirac delta for support site
};

#endif // IBMSUPPORT_H
