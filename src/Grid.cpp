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

// Main solver
void GridClass::solver() {

	// Do grid kernel (LBM)
	lbmKernel();

	// Do object kernel (IBM + FEM)
	if (oPtr->hasIBM == true)
		oPtr->objectKernel();
}

// Main LBM kernel
void GridClass::lbmKernel() {

	// Calculate convective speed
	if (WALL_RIGHT == eConvective)
		convectiveSpeed();

	// Swap to start of timestep
	f_n.swap(f);
	u_n.swap(u);
	rho_n.swap(rho);

	// Start parallel section
#pragma omp parallel
	{

		// Set forcing if Womersley pressure gradient is used
#ifdef WOMERSLEY

		// Loop through all points
#pragma omp for schedule(guided) collapse(2)
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {

				// ID
				int id = i * Ny + j;

				// Calculate forcing
				force_xy[id * dims + eX] = (rho_n[id] * Drho * gravityX + dpdx * cos(2.0 * M_PI * t * Dt / ((SQ(height_p) * M_PI) / (2.0 * SQ(WOMERSLEY) * nu_p)))) * SQ(Dx * Dt) / Dm;
				force_xy[id * dims + eY] = (rho_n[id] * Drho * gravityY + dpdy * cos(2.0 * M_PI * t * Dt / ((SQ(height_p) * M_PI) / (2.0 * SQ(WOMERSLEY) * nu_p)))) * SQ(Dx * Dt) / Dm;
			}
		}
#endif

		// Loop through all points
#pragma omp for schedule(guided) collapse(2)
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {

				// ID
				int id = i * Ny + j;

				// Stream and collide in one go
				streamCollide(i, j, id);

				// Update fluid site macroscopic values
				if (type[id] == eFluid)
					macroscopic(id);
			}
		}

		// Loop through all BC sites and apply BCs before getting macroscopic
#pragma omp for schedule(guided)
		for (size_t bc = 0; bc < BCVec.size(); bc++) {

			// ID
			int id = BCVec[bc];
			int i = id / Ny;
			int j = id - (i * Ny);

			// Apply boundary condition and update macroscopic
			applyBCs(i, j, id);
			macroscopic(id);
		}
	}
}

// Stream and collide in one go (pull algorithm)
inline void GridClass::streamCollide(int i, int j, int id) {

	// Loop through all velocity directions
	for (int v = 0; v < nVels; v++) {

		// Calculate newx and newy then stream and collide
		int src_id = ((i - c[v * dims + eX] + Nx) % Nx) * Ny + ((j - c[v * dims + eY] + Ny) % Ny);

		// Only update if within the box
		f[id * nVels + v] = f_n[src_id * nVels + v] + omega * (equilibrium(src_id, v) - f_n[src_id * nVels + v]) + latticeForce(src_id, v);
	}
}

// Equilibrium function
inline double GridClass::equilibrium(int id, int v) {

	// Get constants
	double C1 = SQ(u_n[id * dims + eX]) + SQ(u_n[id * dims + eY]);
	double C2 = u_n[id * dims + eX] * c[v * dims + eX] + u_n[id * dims + eY] * c[v * dims + eY];

	// Return value
	return rho_n[id] * w[v] * (1.0 + 3.0 * C2 + 4.5 * SQ(C2) - 1.5 * C1);
}

// Discretise cartesian force onto lattice
inline double GridClass::latticeForce(int id, int v) {

	// Return value
	return (1.0 - 0.5 * omega) * (w[v] / SQ(c_s)) *
				((force_xy[id * dims + eX] + force_ibm[id * dims + eX]) * (c[v * dims + eX] - u_n[id * dims + eX] + c[v * dims + eX] * (c[v * dims + eX] * u_n[id * dims + eX] + c[v * dims + eY] * u_n[id * dims + eY]) / SQ(c_s)) +
						((force_xy[id * dims + eY] + force_ibm[id * dims + eY]) * (c[v * dims + eY] - u_n[id * dims + eY] + c[v * dims + eY] * (c[v * dims + eX] * u_n[id * dims + eX] + c[v * dims + eY] * u_n[id * dims + eY]) / SQ(c_s))));
}

// Compute macroscopic quantities
inline void GridClass::macroscopic(int id) {

	// Reset
	rho[id] = 0.0;
	u[id * dims + eX] = 0.0;
	u[id * dims + eY] = 0.0;

	// Sum to find rho and momentum
	for (int v = 0; v < nVels; v++) {
		rho[id] += f[id * nVels + v];
		u[id * dims + eX] += c[v * dims + eX] * f[id * nVels + v];
		u[id * dims + eY] += c[v * dims + eY] * f[id * nVels + v];
	}

	// Divide by rho to get velocity
	u[id * dims + eX] = (u[id * dims + eX] + 0.5 * force_xy[id * dims + eX]) / rho[id];
	u[id * dims + eY] = (u[id * dims + eY] + 0.5 * force_xy[id * dims + eY]) / rho[id];
}

// Apply boundary conditions
void GridClass::applyBCs(int i, int j, int id) {

	// Get ramp coefficient
	double rampCoefficient = getRampCoefficient();

	// Get normal direction
	eDirectionType normalDirection;
	array<int, dims> normalVector = getNormalVector(i, j, normalDirection);

	// Boundary type
	switch (type[id]) {

		// Wall BC
		case eWall:

			// Set velocity to zero
			u_n[id * dims + eX] = 0.0;
			u_n[id * dims + eY] = 0.0;

			// Do regularised BC
			regularisedBC(i, j, id, normalVector, normalDirection);
			break;


		// Velocity BC
		case eVelocity:

			// Set velocity to boundary values
			u_n[id * dims + eX] = u_in[j * dims + eX] * rampCoefficient;
			u_n[id * dims + eY] = u_in[j * dims + eY] * rampCoefficient;

			// Do regularised BC
			regularisedBC(i, j, id, normalVector, normalDirection);
			break;


		// Free slip BC
		case eFreeSlip:

			// Normal velocity is zero
			u_n[id * dims + normalDirection] = 0.0;

			// Extrapolate tangential velocities
			for (int d = 0; d < dims; d++) {
				if (d != normalDirection)
					u_n[id * dims + d] = Utils::zeroGradient(u, normalVector, 2, i, j, d, dims);
			}

			// Do regularised BC
			regularisedBC(i, j, id, normalVector, normalDirection);
			break;


		// Uniform pressure BC
		case ePressure:

			// Set density to boundary values
			rho_n[id] = rho_in[j];

			// Extrapolate tangential velocities
			for (int d = 0; d < dims; d++) {
				if (d != normalDirection)
					u_n[id * dims + d] = Utils::zeroGradient(u, normalVector, 2, i, j, d, dims);
			}

			// Do regularised BC
			regularisedBC(i, j, id, normalVector, normalDirection);
			break;


		// Convective BC
		case eConvective:

			// Do regularised BC
			convectiveBC(j, id);
			break;


		// Fluid (don't do anything)
		case eFluid:
			break;
	}
}

// Regularised BC
void GridClass::regularisedBC(int i, int j, int id, array<int, dims> &normalVector, eDirectionType normalDirection) {

	// If it is a corner then we need to extrapolate
	if (normalVector[eX] != 0 && normalVector[eY] != 0) {

		// If velocity BC then extrapolate density
		if (type[id] == eVelocity || type[id] == eWall || type[id] == eFreeSlip)
			rho_n[id] = Utils::extrapolate(rho, normalVector, 1, i, j);
	}

	// Otherwise normal edge
	else {

		// Declare fplus and fzero
		double fplus = 0.0, fzero = 0.0;

		// Loop through velocities
		for (int v = 0; v < nVels; v++) {

			// If normal is opposite then add to fplus
			if (c[v * dims + normalDirection] == -normalVector[normalDirection])
				fplus += f[id * nVels + v];

			// If it is perpendicular to wall then add to fzero
			else if (c[v * dims + normalDirection] == 0)
				fzero += f[id * nVels + v];
		}

		// Velocity condition
		if (type[id] == eVelocity || type[id] == eWall || type[id] == eFreeSlip) {
			rho_n[id] = (2.0 * fplus + fzero) / (1.0 - normalVector[normalDirection] * u_n[id * dims + normalDirection]);
		}

		// Pressure condition
		else if (type[id] == ePressure) {
			u_n[id * dims + normalDirection] = normalVector[normalDirection] * (1.0 - (2.0 * fplus + fzero) / rho_n[id]);
		}
	}

	// Declare stresses
	double Sxx = 0.0, Syy = 0.0, Sxy = 0.0;

	// Update f values
	for (int v = 0; v < nVels; v++) {

		// Get feq
		double feq = equilibrium(id, v);

		// If corner then unknowns share at least one of normal components
		if (normalVector[eX] != 0 && normalVector[eY] != 0) {
			if (c[v * dims + eX] == normalVector[eX] || c[v * dims + eY] == normalVector[eY]) {

				// If buried link just set to feq
				if (normalVector[eX] * c[v * dims + eX] + normalVector[eY] * c[v * dims + eY] == 0)
					f[id * nVels + v] = feq;
				else
					f[id * nVels + v] = feq + (f[id * nVels + opposite[v]] - equilibrium(id, opposite[v]));
			}
		}

		// If other then unknowns share the normal vector component
		else {
			if (c[v * dims + normalDirection] == normalVector[normalDirection])
				f[id * nVels + v] = feq + (f[id * nVels + opposite[v]] - equilibrium(id, opposite[v]));
		}

		// Store off-equilibrium
		double fneq = f[id * nVels + v] - feq;

		// Compute off-equilbrium stress components
		Sxx += c[v * dims + eX] * c[v * dims + eX] * fneq;
		Syy += c[v * dims + eY] * c[v * dims + eY] * fneq;
		Sxy += c[v * dims + eX] * c[v * dims + eY] * fneq;
	}

	// Compute regularised non-equilibrium components and add to feq to get new populations
	for (int v = 0; v < nVels; v++)
		f[id * nVels + v] = equilibrium(id, v) + (w[v] / (2.0 * QU(c_s))) * (((SQ(c[v * dims + eX]) - SQ(c_s)) * Sxx) + ((SQ(c[v * dims + eY]) - SQ(c_s)) * Syy) + (2.0 * c[v * dims + eX] * c[v * dims + eY] * Sxy));
}

// Convective BC
void GridClass::convectiveBC(int j, int id) {

	// Set the values
	f[id * nVels + 2] = f_n[id * nVels + 2] + 3.0 * w[2] * (delU[j * dims + eX] * c[2 * dims + eX] + delU[j * dims + eY] * c[2 * dims + eY]);
	f[id * nVels + 6] = f_n[id * nVels + 6] + 3.0 * w[6] * (delU[j * dims + eX] * c[6 * dims + eX] + delU[j * dims + eY] * c[6 * dims + eY]);
	f[id * nVels + 8] = f_n[id * nVels + 8] + 3.0 * w[8] * (delU[j * dims + eX] * c[8 * dims + eX] + delU[j * dims + eY] * c[8 * dims + eY]);
}

// Calculate convective speed
void GridClass::convectiveSpeed() {

	// Set uOut to zero
	double uOut = 0.0;

	// Loop through outlet
	for (int j = 0; j < Ny; j++)
		uOut += u[((Nx - 1) * Ny + j) * dims + eX];

	// Get average
	uOut /= static_cast<double>(Ny);

	// Now loop through again and get delU
#pragma omp parallel for schedule(guided)
	for (int j = 0; j < Ny; j++) {
		delU[j * dims + eX] = (-uOut / 2.0) * (3.0 * u[((Nx - 1) * Ny + j) * dims + eX] - 4.0 * u[((Nx - 2) * Ny + j) * dims + eX] + u[((Nx - 3) * Ny + j) * dims + eX]);
		delU[j * dims + eY] = (-uOut / 2.0) * (3.0 * u[((Nx - 1) * Ny + j) * dims + eY] - 4.0 * u[((Nx - 2) * Ny + j) * dims + eY] + u[((Nx - 3) * Ny + j) * dims + eY]);
	}
}

// Get normal vector for bounday site
inline array<int, dims> GridClass::getNormalVector(int i, int j, eDirectionType &normalDirection) {

	// Declare vector
	array<int, dims> normalVector = {0};

	// Get normal direction in x
	if (i == 0) {
		normalVector[eX] = 1;
		normalDirection = eX;
	}
	else if (i == Nx - 1) {
		normalVector[eX] = -1;
		normalDirection = eX;
	}

	// Get normal direction in x
	if (j == 0) {
		normalVector[eY] = 1;
		normalDirection = eY;
	}
	else if (j == Ny - 1) {
		normalVector[eY] = -1;
		normalDirection = eY;
	}

	// Do some extra checking for corners
	if (normalVector[eX] != 0 && normalVector[eY] != 0) {

		// If surrounded by periodic cells then something has gone wrong
		if (type[(i + normalVector[eX]) * Ny + j] == eFluid && type[i * Ny + j + normalVector[eY]] == eFluid)
			ERROR("Corner node is surrounded by fluid lattice sites...exiting");

		// Check in x-direction
		else if (type[(i + normalVector[eX]) * Ny + j] == eFluid) {
			normalDirection = eX;
			normalVector[eY] = 0;
		}

		// Check in y-direction
		else if (type[i * Ny + j + normalVector[eY]] == eFluid) {
			normalDirection = eY;
			normalVector[eX] = 0;
		}
	}

	// Return normal direction
	return normalVector;
}

// Get inlet ramp coefficient
double GridClass::getRampCoefficient() {

	// Get ramp coefficient
#ifdef INLET_RAMP
	if (Dt * t <= INLET_RAMP)
		return (1.0 - cos(M_PI * Dt * t / INLET_RAMP)) / 2.0;
#endif
	return 1.0;
}

// Write info at tInfo frequency
void GridClass::writeInfo() {

	// Calculate max velocity
	double maxVel = 0.0;
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {

			// Get id
			int id = i * Ny + j;

			// Get velocity
			double vel = sqrt(SQ(u[id * dims + eX]) + SQ(u[id * dims + eY]));

			// Check if isnan then break
			if (std::isnan(vel) == true) {

				// Write VTK data
#ifdef VTK
				Utils::writeVTK(*this);
#endif
				// Tell user where it blew up
				ERROR("Simulation blew up (t = " + to_string(t) + ") at i = " + to_string(i) + ", j = " + to_string(j) + "...exiting");
			}

			// If bigger then set
			if (vel > maxVel) {
				maxVel = vel;
			}
		}
	}

	// Calculate values
	double maxRe = (maxVel * Dx / Dt) * ref_L / ref_nu;

	// Calculate new average time per time step
	if (t == tOffset) {
		loopTime = 0.0;
	}
	else {
		loopTime = omp_get_wtime() - startTime;
		loopTime /= (t - tOffset);
	}

	// Get estimated time left
	array<int, 3> hms = Utils::secs2hms(loopTime * (tOffset + nSteps - t));

	// Write out performance information
	cout << endl << endl;
	cout << "Time step " << t << " of " << tOffset + nSteps << endl;
	cout << setprecision(4) << "Simulation has done " << t * Dt << " of " << (tOffset + nSteps) * Dt << " seconds" << endl;
	cout << "Time to finish = " << hms[0] << " [h] " << hms[1] << " [m] " << hms[2] << " [s]" << endl;
	cout << setprecision(4) << "MLUPS = " << Nx * Ny / (1000000.0 * loopTime) << endl;

	// Write out max velocity information
	cout << setprecision(5) << "Max Velocity = " << maxVel << endl;
	cout << setprecision(5) << "Max Velocity (m/s) = " << maxVel * Dx / Dt << endl;
	cout << setprecision(5) << "Max Reynolds number = " << maxRe;
}

// Write log data at start
void GridClass::writeLog() {

	// Add a triple line break if this is a restart
	if (restartFlag == true) {
		ofstream output("Results/Log.out", Utils::io_ofstream::app);
		output << endl << endl << endl;
		output.close();
	}

	// Create stream for output to stdout and file
	Utils::io_ofstream output("Results/Log.out", Utils::io_ofstream::app);
	int ss = static_cast<int>(output.precision());

	// Set precision
	output << setprecision(PRECISION);

	// HEADER
	cout << endl << endl << endl;
	output << "*** LOG INFO ***\n";

	// Get values to write to header
	time_t now = time(0);
	struct tm *localnow = localtime(&now);
	char buffer[80];
	strftime (buffer, 80, "%r on %e/%m/%Y", localnow);

	// TIME
	output << "\nSimulation started at " << buffer;

	// Number of threads
	output << "\n\nRunning with " << Utils::omp_thread_count() << " threads\n";

	// OPTIONS
	output << "\nOPTIONS:\n";

	// Restart
	if (restartFlag == true)
		output << "RESTART = ON\n";
	else
		output << "RESTART = OFF\n";

	// Inlet ramp
#ifdef INLET_RAMP
	output << "Inlet Ramp = " << INLET_RAMP << " s\n";
#else
	output << "Inlet Ramp = OFF\n";
#endif

	// Universal epsilon calculation
#ifdef UNI_EPSILON
	output << "Universal Epsilon Calculation = ON\n";
#else
	output << "Universal Epsilon Calculation = OFF\n";
#endif

	// Ordered reductions
#ifdef ORDERED
	output << "Ordered Reductions = ON\n";
#else
	output << "Ordered Reductions = OFF\n";
#endif

	// Initial deflect
#ifdef INITIAL_DEFLECT
	output << "Initial Deflection = " << 100*INITIAL_DEFLECT << " %L\n";
#else
	output << "Initial Deflection = OFF\n";
#endif

	// Womersley flow
#ifdef WOMERSLEY
	output << "Womersley Number = " << WOMERSLEY << "\n";
#else
	output << "Womersley Number = OFF\n";
#endif

	// OUTPUT OPTIONS
	output << "\nOUTPUT OPTIONS:\n";

	// VTK option
#ifdef VTK
	output << "VTK Output = ON\n";
#else
	output << "VTK Output = OFF\n";
#endif

	// FEM VTK option
#ifdef VTK_FEM
	output << "FEM VTK Output = ON\n";
#else
	output << "FEM VTK Output = OFF\n";
#endif

	// IBM forces
#ifdef FORCES
	output << "Write IBM Forces = ON\n";
#else
	output << "Write IBM Forces = OFF\n";
#endif

	// Tip
#ifdef TIPS
	output << "Write Tip Positions = ON\n";
#else
	output << "Write Tip Positions = OFF\n";
#endif

	// LATTICE VALUES
	output << "\nLATTICE VALUES:\n";
	output << "Nx = "  << Nx << "\n";
	output << "Ny = "  << Ny << "\n";
	output << "omega = "  << omega << "\n";
	output << "tau = "  << tau << "\n";
	output << "nu = "  << nu << "\n";
	output << "Dx = "  << Dx << "\n";
	output << "Dt = "  << Dt << "\n";
	output << "Dm = "  << Dm << "\n";
	output << "Drho = "  << Drho << "\n";

	// PHYSICAL VALUES
	output << "\nPHYSICAL VALUES:\n";
	output << "Length (m) = " << Dx * (Nx - 1) << "\n";
	output << "Height (m) = " << Dx * (Ny - 1) << "\n";
	output << "Density (kg/m^3) = " << rho_p << "\n";
	output << "Viscosity (m^2/s) = " << nu_p << "\n";
	output << "Initial Velocity (m/s) = { " << ux0_p << ", " << uy0_p << " }\n";
	output << "Inlet Velocity (m/s) = { " << uxInlet_p << ", " << uyInlet_p << " }\n";
	output << "Gravity vector (m/s^2) = { " << gravityX << ", " << gravityY << " }\n";
	output << "Pressure gradient vector (Pa/m) = { " << dpdx << ", " << dpdy << " }\n";

	// BOUNDARY CONDITIONS
	output << "\nBOUNDARY CONDITIONS\n";
	output << "Left Wall = " << Utils::getBoundaryString(WALL_LEFT) << "\n";
	output << "Right Wall = " << Utils::getBoundaryString(WALL_RIGHT) << "\n";
	output << "Bottom Wall = " << Utils::getBoundaryString(WALL_BOTTOM) << "\n";
	output << "Top Wall = " << Utils::getBoundaryString(WALL_TOP) << "\n";

	// TIME STEP
	output << "\nTIME STEP:\n";
	output << "Time = " << Dt * nSteps << " s\n";
	output << "nSteps = " << nSteps << "\n";
	output << "tInfo = " << tinfo << "\n";
#ifdef VTK
	output << "tVTK = " << tVTK << "\n";
#endif
	output << "tRestart = " << tRestart << "\n";

	// If doing a restart then write out how many time steps for this run
	if (restartFlag == true) {
		output << "Total Time = " << Dt * (tOffset + nSteps) << " s\n";
		output << "Total nSteps = " << (tOffset + nSteps) << "\n";
	}

	// REFERENCE VALUES
	output << "\nREFERENCE VALUES:\n";
	output << "Reference Viscosity (m^2/s) = " << ref_nu << "\n";
	output << "Reference Density (kg/m^3) = " << ref_rho << "\n";
	output << "Reference Pressure (Pa) = " << ref_P << "\n";
	output << "Reference Length (m) = " << ref_L << "\n";
	output << "Reference Velocity = " << ref_U * Dt / Dx << "\n";
	output << "Reference Velocity (m/s) = " << ref_U << "\n";
	output << "Reynolds Number = " << ref_U * ref_L / ref_nu;

	// Set precision
	output << setprecision(ss);

	// Close
	output.close();
}

// Write fluid VTK
void GridClass::writeVTK() {

	// Get the endianness
	string endianStr = (bigEndian ? "BigEndian" : "LittleEndian");

	// Create file
	ofstream output;
	output.open("Results/VTK/Fluid." + to_string(t) + ".vti", ios::binary);

	// Handle failure to open
	if (!output.is_open())
		ERROR("Error opening fluid VTK file...exiting");

	// Write XML header
	output << "<?xml version=\"1.0\"?>\n";

	// Begin VTK file
	int level = 0;
	output << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"" << endianStr << "\" header_type=\"UInt64\">\n";

	// New level -> ImageData
	level = 1;
	output << string(level, '\t') << "<ImageData "
								  << "WholeExtent=\"" << 0.0 << " " << Nx - 1 << " " << 0.0 << " " << Ny - 1 << " " << 0.0 << " " << 0.0 << "\" "
								  << "Origin=\"" << 0.0 << " " << 0.0 << " " << 0.0 << "\" "
								  << "Spacing=\"" << Dx << " " << Dx << " " << Dx << "\">\n";

	// New level -> Piece
	level = 2;
	output << string(level, '\t') << "<Piece Extent=\"" << 0.0 << " " << Nx - 1 << " " << 0.0 << " " << Ny - 1 << " " << 0.0 << " " << 0.0 << "\">\n";

	// New level -> PointData
	level = 3;
	output << string(level, '\t') << "<PointData>\n";

	// New level -> DataArray
	level = 4;

	// Density
	output << string(level, '\t') << "<DataArray type=\"Float64\" Name=\"Density\" format=\"appended\" offset=\"" << 0 << "\"/>\n";

	// Pressure
	output << string(level, '\t') << "<DataArray type=\"Float64\" Name=\"Pressure\" format=\"appended\" offset=\"" << 1*(Nx*Ny*sizeof(double) + sizeof(unsigned long long)) << "\"/>\n";

	// Velocity
	output << string(level, '\t') << "<DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << 2*(Nx*Ny*sizeof(double) + sizeof(unsigned long long)) << "\"/>\n";

	// Back level -> PointData
	level = 3;
	output << string(level, '\t') << "</PointData>\n";

	// Back level -> Piece
	level = 2;
	output << string(level, '\t') << "</Piece>\n";

	// Back level -> ImageData
	level = 1;
	output << string(level, '\t') << "</ImageData>\n";

	// New level -> AppendedData
	level = 1;
	output << string(level, '\t') << "<AppendedData encoding=\"raw\">\n";

	// New level -> Raw data
	level = 2;
	output << string(level, '\t') << "_";

	// Density
	unsigned long long size = Nx * Ny * sizeof(double);
	output.write((char*)&size, sizeof(unsigned long long));
	for (int j = 0; j < Ny; j++) {
		for (int i = 0; i < Nx; i++) {
			double density = rho[i * Ny + j] * Drho;
			output.write((char*)&density, sizeof(double));
		}
	}

	// Pressure
	size = Nx * Ny * sizeof(double);
	output.write((char*)&size, sizeof(unsigned long long));
	for (int j = 0; j < Ny; j++) {
		for (int i = 0; i < Nx; i++) {
			double pressure = ref_P + (rho[i * Ny + j] - rho_p / Drho) * SQ(c_s) * Dm / (Dx * SQ(Dt));
			output.write((char*)&pressure, sizeof(double));
		}
	}

	// Velocity
	size = 3 * Nx * Ny * sizeof(double);
	output.write((char*)&size, sizeof(unsigned long long));
	for (int j = 0; j < Ny; j++) {
		for (int i = 0; i < Nx; i++) {
			double ux = u[(i * Ny + j) * dims + eX] * (Dx / Dt);
			double uy = u[(i * Ny + j) * dims + eY] * (Dx / Dt);
			double uz = 0.0;
			output.write((char*)&ux, sizeof(double));
			output.write((char*)&uy, sizeof(double));
			output.write((char*)&uz, sizeof(double));
		}
	}

	// Back level -> AppendedData
	level = 1;
	output << "\n" << string(level, '\t') << "</AppendedData>\n";

	// Back level -> VTKFile
	level = 0;
	output << string(level, '\t') << "</VTKFile>\n";

	// Close file
	output.close();

	// Check if we should write some blank VTK files for the body
	if (!oPtr->hasIBM && boost::filesystem::exists("Results/VTK/IBM.0.vtp")) {
		string fname = "Results/VTK/IBM." + to_string(t) + string(".vtp");
		oPtr->writeEmptyVTK(fname);
	}
#ifdef VTK_FEM
	if (!oPtr->hasFlex && boost::filesystem::exists("Results/VTK/FEM.0.vtp")) {
		string fname = "Results/VTK/FEM." + to_string(t) + string(".vtp");
		oPtr->writeEmptyVTK(fname);
	}
#endif
}

// Initialise grid values
void GridClass::initialiseGrid() {

	// First check to make we don't have convective BC anywhere other than RHS
	if (WALL_LEFT == eConvective || WALL_BOTTOM == eConvective || WALL_TOP == eConvective)
		ERROR("Currently convective BC only supported for right boundary...exiting");
	else if (WALL_RIGHT == eConvective)
		delU.resize(Ny * dims, 0.0);

	// Set type matrix
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {

			// Get id
			int id = i * Ny + j;

			// Left wall
			if (i == 0)
				type[id] = WALL_LEFT;

			// Right wall
			else if (i == Nx - 1)
				type[id] = WALL_RIGHT;

			// Bottom wall
			if (j == 0)
				type[id] = WALL_BOTTOM;

			// Top wall
			else if (j == Ny - 1)
				type[id] = WALL_TOP;

			// Add to the BC vector
			if (type[id] != eFluid)
				BCVec.push_back(id);
		}
	}

	// Loop through and set inlet profile
	for (int j = 0; j < Ny; j++) {

#ifndef PROFILE

		// Set initial velocity to be uniform
		u_in[j * dims + eX] = uxInlet_p * Dt / Dx;
		u_in[j * dims + eY] = uyInlet_p * Dt / Dx;


#else
		// Generate velocity profile
		if (PROFILE == eParabolic) {

			// Get parameters for working out parabolic profile
			double R = height_p / 2.0;
			double YPos = j * Dx - R;

			// Set velocity
			u_in[j * dims + eX] = 1.5 * (uxInlet_p * Dt / Dx) * (1.0 - SQ(YPos / R));
			u_in[j * dims + eY] = 1.5 * (uyInlet_p * Dt / Dx) * (1.0 - SQ(YPos / R));
		}
		else if (PROFILE == eShear) {

			// Get parameters for working out shear profile
			double H = height_p;
			double YPos = j * Dx;

			// Set velocity
			u_in[j * dims + eX] = (uxInlet_p * Dt / Dx) * (YPos / H);
			u_in[j * dims + eY] = (uyInlet_p * Dt / Dx) * (YPos / H);
		}
		else if (PROFILE == eBoundaryLayer) {

			// Get parameters for working out shear profile
			double H = height_p;
			double YPos = j * Dx;

			// Set velocity
			u_in[j * dims + eX] = ((1.5 * uxInlet_p * Dt / Dx) / SQ(H)) * YPos * (2.0 * H - YPos);
			u_in[j * dims + eY] = ((1.5 * uyInlet_p * Dt / Dx) / SQ(H)) * YPos * (2.0 * H - YPos);
		}
#endif
	}

	// Set initial velocity
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {

			// Get id
			int id = i * Ny + j;

#ifdef INLET_RAMP
			u[id * dims + eX] = 0.0;
			u[id * dims + eY] = 0.0;
#else
#ifdef PROFILE

			// Set velocity
			u[id * dims + eX] = u_in[j * dims + eX];
			u[id * dims + eY] = u_in[j * dims + eY];
#else

			// Set initial velocity to be uniform
			u[id * dims + eX] = ux0_p * Dt / Dx;
			u[id * dims + eY] = uy0_p * Dt / Dx;
#endif
#endif

			// If a wall then set to zero
			if (type[id] == eWall) {
				u[id * dims + eX] = 0.0;
				u[id * dims + eY] = 0.0;
			}
		}
	}

	// Set start of time step values
	u_n = u;
	rho_n = rho;

	// Set the xy forces
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {

			// ID
			int id = i * Ny + j;

			// Set the xyz cartesian forces
			force_xy[id * dims + eX] = (rho[id] * Drho * gravityX + dpdx) * SQ(Dx * Dt) / Dm;
			force_xy[id * dims + eY] = (rho[id] * Drho * gravityY + dpdy) * SQ(Dx * Dt) / Dm;
		}
	}

	// Set f values to equilibrium
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {

			// ID
			int id = i * Ny + j;

			// Loop though vels
			for (int v = 0; v < nVels; v++)
				f[id * nVels + v] = equilibrium(id, v);
		}
	}

	// Set start of time step values
	f_n = f;
}

// Start the clock for getting MLUPS
void GridClass::startClock() {

	// Start clock
	startTime = omp_get_wtime();
}

// Read in restart file
void GridClass::readRestart() {

	// Open file
	ifstream file;
	file.open("Results/Restart/Fluid.restart", ios::binary);

	// Handle failure to open
	if (!file.is_open())
		ERROR("Error opening Fluid.restart file...exiting");

	// Declare values
	int tRead, nxRead, nyRead;
	double omegaRead, dxRead, dtRead, dmRead;

	// Read in global info
	file.read((char*)&tRead, sizeof(int));
	file.read((char*)&nxRead, sizeof(int));
	file.read((char*)&nyRead, sizeof(int));
	file.read((char*)&omegaRead, sizeof(double));
	file.read((char*)&dxRead, sizeof(double));
	file.read((char*)&dtRead, sizeof(double));
	file.read((char*)&dmRead, sizeof(double));

	// Swap byte order if bigEndian
	int nxSwap = (bigEndian ? Utils::swapEnd(nxRead) : nxRead);
	int nySwap = (bigEndian ? Utils::swapEnd(nyRead) : nyRead);
	double omegaSwap = (bigEndian ? Utils::swapEnd(omegaRead) : omegaRead);
	double dxSwap = (bigEndian ? Utils::swapEnd(dxRead) : dxRead);
	double dtSwap = (bigEndian ? Utils::swapEnd(dtRead) : dtRead);
	double dmSwap = (bigEndian ? Utils::swapEnd(dmRead) : dmRead);

	// Check they match up
	if (Nx != nxSwap || Ny != nySwap || omega != omegaSwap || Dx != dxSwap || Dt != dtSwap || Dm != dmSwap)
		ERROR("Grid size/scaling has changed between runs...this is not supported");

	// Get the time to start from
	tOffset = (bigEndian ? Utils::swapEnd(tRead) : tRead);
	t = tOffset;

	// Now loop through each lattice site and read necessary data
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {

			// ID
			int id = i * Ny + j;

			// Declare values
			int iRead, jRead;
			double rhoRead, uxRead, uyRead, fxRead, fyRead;

			// Read in data
			file.read((char*)&iRead, sizeof(int));
			file.read((char*)&jRead, sizeof(int));
			file.read((char*)&rhoRead, sizeof(double));
			file.read((char*)&uxRead, sizeof(double));
			file.read((char*)&uyRead, sizeof(double));
			file.read((char*)&fxRead, sizeof(double));
			file.read((char*)&fyRead, sizeof(double));

			// Swap byte order if bigEndian
			int iSwap = (bigEndian ? Utils::swapEnd(iRead) : iRead);
			int jSwap = (bigEndian ? Utils::swapEnd(jRead) : jRead);
			double rhoSwap = (bigEndian ? Utils::swapEnd(rhoRead) : rhoRead);
			double uxSwap = (bigEndian ? Utils::swapEnd(uxRead) : uxRead);
			double uySwap = (bigEndian ? Utils::swapEnd(uyRead) : uyRead);
			double fxSwap = (bigEndian ? Utils::swapEnd(fxRead) : fxRead);
			double fySwap = (bigEndian ? Utils::swapEnd(fyRead) : fyRead);

			// Check they match up
			if (i != iSwap || j != jSwap)
				ERROR("Grid indices do not match Fluid.restart file...exiting");

			// Read into grid data
			rho[id] = rhoSwap;
			u[id * dims + eX] = uxSwap;
			u[id * dims + eY] = uySwap;
			force_ibm[id * dims + eX] = fxSwap;
			force_ibm[id * dims + eY] = fySwap;

			// Read in f values
			for (int v = 0; v < nVels; v++) {
				double fRead;
				file.read((char*)&fRead, sizeof(double));
				double fSwap = (bigEndian ? Utils::swapEnd(fRead) : fRead);
				f[id * nVels + v] = fSwap;
			}
		}
	}
}

// Read in restart file
void GridClass::writeRestart() {

	// Create file
	ofstream output;
	output.open("Results/Restart/Fluid.restart", ios::binary);

	// Handle failure to open
	if (!output.is_open())
		ERROR("Error opening Fluid.restart file...exiting");

	// Swap byte order if bigEndian
	int tWrite = (bigEndian ? Utils::swapEnd(t) : t);
	int nxWrite = (bigEndian ? Utils::swapEnd(Nx) : Nx);
	int nyWrite = (bigEndian ? Utils::swapEnd(Ny) : Ny);
	double omegaWrite = (bigEndian ? Utils::swapEnd(omega) : omega);
	double dxWrite = (bigEndian ? Utils::swapEnd(Dx) : Dx);
	double dtWrite = (bigEndian ? Utils::swapEnd(Dt) : Dt);
	double dmWrite = (bigEndian ? Utils::swapEnd(Dm) : Dm);

	// Write out global information
	output.write((char*)&tWrite, sizeof(int));
	output.write((char*)&nxWrite, sizeof(int));
	output.write((char*)&nyWrite, sizeof(int));
	output.write((char*)&omegaWrite, sizeof(double));
	output.write((char*)&dxWrite, sizeof(double));
	output.write((char*)&dtWrite, sizeof(double));
	output.write((char*)&dmWrite, sizeof(double));

	// Now loop through each lattice site and write necessary data
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {

			// ID
			int id = i * Ny + j;

			// Swap byte order if bigEndian
			int iWrite = (bigEndian ? Utils::swapEnd(i) : i);
			int jWrite = (bigEndian ? Utils::swapEnd(j) : j);
			double rhoWrite = (bigEndian ? Utils::swapEnd(rho[id]) : rho[id]);
			double uxWrite = (bigEndian ? Utils::swapEnd(u[id * dims + eX]) : u[id * dims + eX]);
			double uyWrite = (bigEndian ? Utils::swapEnd(u[id * dims + eY]) : u[id * dims + eY]);
			double fxWrite = (bigEndian ? Utils::swapEnd(force_ibm[id * dims + eX]) : force_ibm[id * dims + eX]);
			double fyWrite = (bigEndian ? Utils::swapEnd(force_ibm[id * dims + eY]) : force_ibm[id * dims + eY]);

			// Write out global information
			output.write((char*)&iWrite, sizeof(int));
			output.write((char*)&jWrite, sizeof(int));
			output.write((char*)&rhoWrite, sizeof(double));
			output.write((char*)&uxWrite, sizeof(double));
			output.write((char*)&uyWrite, sizeof(double));
			output.write((char*)&fxWrite, sizeof(double));
			output.write((char*)&fyWrite, sizeof(double));

			// Write out f values
			for (int v = 0; v < nVels; v++) {
				double fWrite = (bigEndian ? Utils::swapEnd(f[id * nVels + v]) : f[id * nVels + v]);
				output.write((char*)&fWrite, sizeof(double));
			}
		}
	}

	// Close file
	output.close();
}

// Constructor
GridClass::GridClass() {

	// Create directories
	restartFlag = Utils::createDirectories();

	// Write out
	cout << endl << endl << "Initialising grid...";

	// Initially set objects to NULL
	oPtr = NULL;

	// Some default parameters
	double rho0 = 1.0;

	// D2Q9 parameters (can be hard-coded)
	c_s = 1.0 / sqrt(3.0);
	w = {4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
	c = {0, 0, 1, 0, -1, 0, 0, 1, 0, -1, 1, 1, -1, -1, 1, -1, -1, 1};
	opposite = {0, 2, 1, 4, 3, 6, 5, 8, 7};

	// Initialise parameters
	t = 0;
	tOffset = 0;
	tau = 1.0 / omega;
	nu = (tau - 0.5) * SQ(c_s);
	Dx = height_p / (Ny - 1);
	Dt = SQ(Dx) * nu / nu_p;
	Dm = (rho_p / rho0) * TH(Dx);
	Drho = (rho_p / rho0);

	// Call clock
	startTime = omp_get_wtime();
	loopTime = 0.0;

	// Set the sizes and initialise arrays
	u.resize(Nx * Ny * dims, 0.0);
	u_n.resize(Nx * Ny * dims, 0.0);
	rho.resize(Nx * Ny, rho0);
	rho_n.resize(Nx * Ny, rho0);
	force_xy.resize(Nx * Ny * dims, 0.0);
	force_ibm.resize(Nx * Ny * dims, 0.0);
	type.resize(Nx * Ny, eFluid);
	f.resize(Nx * Ny * nVels, 0.0);
	f_n.resize(Nx * Ny * nVels, 0.0);

	// Set sizes of helper arrays
	u_in.resize(Ny * dims, 0.0);
	rho_in.resize(Ny, rho0);

	// Initialise grid
	initialiseGrid();

	// Check for big or little endian
	bigEndian = Utils::isBigEndian();

	// Write out header
	cout << "finished";
}
