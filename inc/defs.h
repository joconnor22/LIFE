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

#ifndef DEFS_H	// DEFS_H
#define DEFS_H

// Include headers and set namespace
#include <iostream>
#include <iomanip>
#include <cmath>
#include <omp.h>
#include <boost/filesystem.hpp>
using namespace std;

// Version and release date
const string version = "v1.0.3";
const string date = "10th June 2020";

// Number of dimensions and lattice velocities
const int dims = 2;
const int nVels = 9;

// IBM support buffer size
const int suppSize = 9;

// Set default values for 2D beam
const int nodeDOFs = 3;
const int elementNodes = 2;
const int elementDOFs = elementNodes * nodeDOFs;

// Enumerations
enum eDirectionType {eX, eY};
enum eFlexibleType {eFlexible, eRigid};
enum eBCType {eClamped, eSupported};
enum eBodyType {eCircle, eFilament};
enum eLatType {eFluid, eWall, eVelocity, eFreeSlip, ePressure, eConvective};
enum eProfileType {eParabolic, eShear, eBoundaryLayer};

// Macros
#define SQ(x) ((x) * (x))
#define TH(x) ((x) * (x) * (x))
#define QU(x) ((x) * (x) * (x) * (x))
#define ERROR Utils::errorExit
#define WARN Utils::warning

#endif	// DEFS_H
