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
#include "../inc/Grid.h"
#include "../inc/Objects.h"
#include "../inc/Utils.h"

// ***** Main function ***** //
int main() {

	// Write out header
	Utils::writeHeader();

	// Start initialising
	cout << endl << endl << "*** INITIALISING LIFE ***";

	// Time the code
	double start = omp_get_wtime();

	// Set number of threads
#ifdef THREADS
	omp_set_num_threads(THREADS);
#endif

	// Create grid and initialise values
	GridClass grid;

	// If IBM then setup the body object
	ObjectsClass objects(grid);

	// Read in restart file if required
	if (grid.restartFlag == true) {

		// Read in restart file
		Utils::readRestart(grid);
	}

	// Write log
	Utils::writeLog(grid);

	// Starting main algorithm
	cout << endl << endl << endl << "*** STARTING LIFE ***";

	// Write out info
	Utils::writeInfo(grid);

	// Write VTK
	Utils::writeVTK(grid);

	// Start the clock
	grid.startClock();

	// ***** MAIN LBM ALGORITHM ***** //
	for (grid.t = grid.tOffset + 1; grid.t <= grid.tOffset + nSteps; grid.t++) {

		// Run the main LBM solver
		grid.solver();

		// Write out info
		if (grid.t % tinfo == 0)
			Utils::writeInfo(grid);

		// Write VTK
		if (grid.t % tVTK == 0)
			Utils::writeVTK(grid);

		// Write restart data
		if (tRestart > 0 && grid.t % tRestart == 0)
			Utils::writeRestart(grid);
	}

	// Starting main algorithm
	cout << endl << endl << endl << "*** FINISHED LIFE ***";

	// Finish timing
	double end = omp_get_wtime();
	cout << fixed << setprecision(2) << endl << endl << "Simulation took " << end - start << " seconds" << endl << endl << endl;
}
