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
#include "../inc/Utils.h"
#include "../inc/Grid.h"
#include "../inc/Objects.h"

// Write header at start of run
void Utils::writeHeader() {

	// Write out version number and release date
	cout << endl << "*** LIFE " << version << " (" << date << ") ***" << endl << endl;

	// Write out license header
	cout << "LIFE: Lattice boltzmann-Immersed boundary-Finite Element" << endl;
	cout << "Copyright (C) 2019 Joseph O'Connor" << endl << endl;

	cout << "This program is free software: you can redistribute it and/or modify" << endl;
	cout << "it under the terms of the GNU General Public License as published by" << endl;
	cout << "the Free Software Foundation, either version 3 of the License, or" << endl;
	cout << "(at your option) any later version." << endl << endl;

	cout << "This program is distributed in the hope that it will be useful," << endl;
	cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl;
	cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << endl;
	cout << "GNU General Public License for more details." << endl << endl;

	cout << "You should have received a copy of the GNU General Public License" << endl;
	cout << "along with this program.  If not, see <http://www.gnu.org/licenses/>." << endl;
}

// Create directories at start
bool Utils::createDirectories() {

	// First check if restart file exists
	if (boost::filesystem::exists("Results/Restart/Fluid.restart"))
		return true;

	// Remove directories
	if (boost::filesystem::exists("Results")) {
		if (!boost::filesystem::remove_all("Results"))
			ERROR("Problem removing results directory...exiting");
	}

	// Create results directory
	if (!boost::filesystem::create_directory("Results"))
		ERROR("Problem creating results directory...exiting");

	// Create VTK directory
#ifdef VTK
	if (!boost::filesystem::create_directory("Results/VTK"))
		ERROR("Problem creating vtk directory...exiting");
#endif

	// Only create restart directory if we need to
	if (tRestart > 0) {
		if (!boost::filesystem::create_directory("Results/Restart"))
			ERROR("Problem creating restart directory...exiting");
	}

	// If we make it here then this is not a restart case
	return false;
}

// Read in restart file
void Utils::readRestart(GridClass &grid) {

	// Starting main algorithm
	cout << endl << endl << endl << "*** READ IN RESTART DATA ***";

	// Reading in restart files
	cout << endl << endl << "Reading in restart files...";

	// Read grid restart data
	grid.readRestart();

	// Read bodies restart data
	grid.oPtr->readRestart();

	// Delete VTK files that were written after last restart written
	Utils::deleteVTKs(grid);

	// Write out header
	cout << "finished";
}

// Write out restart file
void Utils::writeRestart(GridClass &grid) {

	// Wrting restart files
	cout << endl << "Writing restart data...";

	// Write grid info
	grid.writeRestart();

	// Write IBM info
	grid.oPtr->writeRestart();

	// Write header
	cout << "finished";
}

// Write info
void Utils::writeInfo(GridClass &grid) {

	// Write grid info
	grid.writeInfo();

	// Write IBM info
	grid.oPtr->writeInfo();

	// Write out forces on bodies
#ifdef FORCES
	grid.oPtr->writeTotalForces();
#endif

	// Write out forces on bodies
#ifdef TIPS
	grid.oPtr->writeTips();
#endif
}

// Write log
void Utils::writeLog(GridClass &grid) {

	// Write grid info
	grid.writeLog();

	// Write IBM info
	grid.oPtr->writeLog();
}

// Write VTK
void Utils::writeVTK(GridClass &grid) {

	// Write VTK
#ifdef VTK

	// Write header
	cout << endl << "Writing VTK data...";

	// Write grid VTK
	grid.writeVTK();

	// Write IBM VTK
	grid.oPtr->writeVTK();

	// Write header
	cout << "finished";
#endif
}

// Delete future VTKs
void Utils::deleteVTKs(GridClass &grid) {

	// Results path
	string path = "Results/VTK";

	// Check if it exists
	if (boost::filesystem::exists(path)) {

		// Loop through three different file types
		for (int i = 0; i < 3; i++) {

			// File strings
			string fileStr, extStr;

			// Get types
			if (i == 0) {
				fileStr = "Fluid";
				extStr = ".vti";
			}
			else if (i == 1) {
				fileStr = "IBM";
				extStr = ".vtp";
			}
			else if (i == 2) {
				fileStr = "FEM";
				extStr = ".vtp";
			}

			// Get directory iterator
			boost::filesystem::directory_iterator endit;

			// Loop through all fluid vtk files
			for (boost::filesystem::directory_iterator it(path); it != endit; it++) {

				// Get string
				string fStr = it->path().stem().string();

				// Check to make sure it is the files we want
				if (fStr.find(fileStr) != string::npos && it->path().extension() == extStr) {

					// Get time step
					string tStepStr = fStr.substr(fileStr.size() + 1, fStr.size());

					// Convert to int
					int tStep = stoi(tStepStr);

					// If tStep is bigger than current time then delete
					if (tStep >= grid.t) {
						if (!boost::filesystem::remove(path + "/" + fStr + extStr))
							ERROR("Problem removing future VTK files...exiting");
					}
				}
			}
		}
	}
}

// Get number of omp threads (as built in doesn't work on GCC)
int Utils::omp_thread_count() {

	// Thread counters
	int n = 0;

	// Do parallel reduction on n
#pragma omp parallel reduction(+:n)
	n += 1;

	// Return total thread count
	return n;
}

// Convert seconds to hours:minutes:seconds
array<int, 3> Utils::secs2hms(double seconds) {

	// Round to nearest second
	seconds = round(seconds);

	// Get number of hours and number of minutes
	double hours = seconds / (60.0 * 60.0);
	double minutes = seconds / (60.0);

	// Declare vector
	array<int, 3> hms;

	// Get hours:minutes:seconds
	hms[0] = static_cast<int>(floor(hours));
	hms[1] = static_cast<int>(floor(minutes - (static_cast<double>(hms[0]) * 60.0)));
	hms[2] = static_cast<int>(floor(static_cast<double>(seconds) - static_cast<double>(hms[1]) * 60.0 - static_cast<double>(hms[0]) * 60.0 * 60.0));

	// Return
	return hms;
}

// Get string for boundary condition
string Utils::getBoundaryString(eLatType BCType) {

	// String
	string str;

	// Check against possible options
	if (BCType == eFluid)
		str = "Periodic BC";
	else if (BCType == eVelocity)
		str = "Velocity BC";
	else if (BCType == ePressure)
		str = "Pressure BC";
	else if (BCType == eWall)
		str = "Wall BC";
	else if (BCType == eFreeSlip)
		str = "Free Slip BC";
	else if (BCType == eConvective)
		str = "Convective BC";

	// Return
	return str;
}

// Solve linear system using LAPACK routines
vector<double> Utils::solveLAPACK(vector<double> A, vector<double> b, int BC) {

	// Set up the correct values
	char trans = 'T';
	int dim = static_cast<int>(sqrt(static_cast<int>(A.size())));
	int row = dim - BC;
	int col = dim - BC;
	int offset = BC * dim + BC;
    int nrhs = 1;
    int LDA = dim;
    int LDB = dim;
    int info;
	std::vector<int> ipiv(row, 0);

    // Factorise and solve
	dgetrf_(&row, &col, A.data() + offset, &LDA, ipiv.data(), &info);
	dgetrs_(&trans, &row, &nrhs, A.data() + offset, &LDA, ipiv.data(), b.data() + BC, &LDB, &info);

	// Set return values not included to zero
	fill(b.begin(), b.begin() + BC, 0.0);

	// Return RHS
	return b;
}
