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
#include "../inc/Objects.h"
#include "../inc/FEMBody.h"
#include "../inc/Grid.h"
#include "../inc/Utils.h"

// Main kernel for objects
void ObjectsClass::objectKernel() {

	// While loop parameters
	subIt = 0;
	int MAXIT = 20;

	// Subiteration loop
	do {

		// Do predictor, recompute support, ds, and epsilon
		if (hasFlex == true)
			recomputeObjectVals();

		// Do IBM interpolation step
		ibmKernelInterp();

		// If only rigid bodies then we can break
		if (hasFlex == false)
			break;

		// Do FEM step
		femKernel();

		// Increase subIt
		subIt++;

	} while (subIt < MAXIT && subRes > subTol);

	// Do IBM spreading step
	ibmKernelSpread();

	// If it reached max iterations then exit
	if (subIt == MAXIT)
		ERROR("Subiteration scheme hit " + to_string(subIt) + " iterations...exiting");
}

// Do FEM and update IBM positions and velocities
void ObjectsClass::femKernel() {

	// Declare residual parameters
	double res = 0.0;
	double num = 0.0;
	double den = 0.0;

	// Loop through all bodies, do FEM, and then sum to get residual values
#ifdef ORDERED
	#pragma omp parallel for ordered schedule(dynamic,1)
#else
	#pragma omp parallel for schedule(guided) reduction(+:res, num, den)
#endif
	for (size_t ib = 0; ib < iBody.size(); ib++) {
		if (iBody[ib].flex == eFlexible) {

			// Do the dynamic FEM routine
			iBody[ib].sBody->dynamicFEM();

			// Either do ordered or non-deterministic sum (it can affect results)
#ifdef ORDERED
		#pragma omp ordered
#endif
			{
				// Sum to get global values
				res += iBody[ib].sBody->subRes;
				num += iBody[ib].sBody->subNum;
				den += iBody[ib].sBody->subDen;
			}
		}
	}

	// Set global values
	subRes = sqrt(res) / (ref_L * sqrt(static_cast<double>(simDOFs)));
	subNum = num;
	subDen = den;
}

// Interpolate and force calc
void ObjectsClass::ibmKernelInterp() {

	// Reset IBM forces
	fill(gPtr->force_ibm.begin(), gPtr->force_ibm.end(), 0.0);

	// Loop through all bodies and nodes
#pragma omp parallel for schedule(guided)
	for (size_t i = 0; i < iNode.size(); i++) {

		// Interpolate
		iNode[i].interpolate();

		// Force calculation
		iNode[i].forceCalc();
	}
}

// Force spread and update macro
void ObjectsClass::ibmKernelSpread() {

	// Reset IBM forces
	fill(gPtr->force_ibm.begin(), gPtr->force_ibm.end(), 0.0);

	// Start parallel section
#pragma omp parallel
	{

		// Loop through all bodies and nodes
#ifdef ORDERED
	#pragma omp for ordered schedule(dynamic,1)
#else
	#pragma omp for schedule(guided)
#endif
		for (size_t i = 0; i < iNode.size(); i++) {

			// Force spread
			iNode[i].spread();
		}

		// Loop through all bodies and nodes
#pragma omp for schedule(guided)
		for (size_t i = 0; i < iNode.size(); i++) {

			// Update macroscopic
			iNode[i].updateMacroscopic();
		}
	}
}

// Recompute support, ds, and epsilon and subiteration values
void ObjectsClass::recomputeObjectVals() {

	// Start parallel section
#pragma omp parallel
	{

		// Do predictor step if first iteration
		if (subIt == 0) {

			// Loop through bodies, set start of time step, and do predictor (if on)
#pragma omp for schedule(guided)
			for (size_t ib = 0; ib < iBody.size(); ib++) {
				if (iBody[ib].flex == eFlexible) {

					// Set the start of timestep values
					iBody[ib].sBody->resetValues();

					// Predictor step
					iBody[ib].sBody->predictor();
				}
			}
		}

		// Else calculate new relaxation parameter and do the relaxation
		else {

			// Get relaxation value (only one thread)
#pragma omp single
			{
				if (subIt == 1) {

					// Use max relax if need be on first iteration
					relax = static_cast<double>(Utils::sgn(relax) * min(fabs(relax), relaxMax));
				}
				else {

					// Use global residual values to get next relaxation factor
					relax = -relax * subNum / subDen;
				}
			}

			// Relax displacements and update IBM
#pragma omp for schedule(guided)
			for (size_t ib = 0; ib < iBody.size(); ib++) {
				if (iBody[ib].flex == eFlexible) {

					// Apply relaxation
					iBody[ib].sBody->U = iBody[ib].sBody->U_km1 + relax * (iBody[ib].sBody->U - iBody[ib].sBody->U_km1);

					// Update FEM elements
					iBody[ib].sBody->updateFEMValues();

					// Update the velocity
					iBody[ib].sBody->finishNewmark();

					// Update IBM values
					iBody[ib].sBody->updateIBMValues();

					// Set previous iteration value
					iBody[ib].sBody->U_km1.swap(iBody[ib].sBody->U);
				}
			}
		}

		// Loop through all bodies and nodes
#pragma omp for schedule(guided)
		for (size_t i = 0; i < iNode.size(); i++) {
			if (iNode[i].iPtr->flex == eFlexible) {

				// Find support
				iNode[i].findSupport();

				// Compute ds
				iNode[i].computeDs();
			}
		}
	}

	// Compute epsilon
	computeEpsilon();
}

// Compute epsilon
void ObjectsClass::computeEpsilon() {

	// If universal calculation then need to move all markers into tmp iBody holder
#ifdef UNI_EPSILON

	// Temporary iBody for storing all markers in whole simulation
	vector<IBMBodyClass> iBodyTmp;

	// Call constructor with vector of all iBodies
	iBodyTmp.emplace_back(iNode);

	// Set pointer to stop copying of data
	vector<IBMBodyClass> *iBodyPtr = &iBodyTmp;
#else

	// If separated epsilon then just set pointer to point to existing vector
	vector<IBMBodyClass> *iBodyPtr = &iBody;
#endif

	// Get lattice spacing
	double Dx = gPtr->Dx;

	// Loop through all bodies and get epsilon
#pragma omp parallel for schedule(guided)
	for (size_t ib = 0; ib < (*iBodyPtr).size(); ib++) {

		// Do if first time step; if not first time step then only do if flexible
		if(gPtr->t == 0 || (*iBodyPtr)[ib].flex == eFlexible) {

			// Get size of A matrix
			size_t dim = (*iBodyPtr)[ib].node.size();

			// Set A matrix
			vector<double> A(dim * dim, 0.0);

			// Loop through all nodes
			for (size_t i = 0; i < dim; i++) {

				// Set node i
				IBMNodeClass *nodei = (*iBodyPtr)[ib].node[i];

				// Loop though all nodes again
				for (size_t j = 0; j < dim; j++) {

					// Set node j
					IBMNodeClass *nodej = (*iBodyPtr)[ib].node[j];

					// Now loop through all support markers for node i
					for (size_t s = 0; s < nodei->suppCount; s++) {

						// Dirac delta value of support marker for node i and support s
						double diracVal_i = nodei->supp[s].diracVal;

						// Delta value between node i support s and node j
						double distX = fabs(nodej->pos[eX] / Dx - nodei->supp[s].idx);
						double distY = fabs(nodej->pos[eY] / Dx - nodei->supp[s].jdx);
						double diracVal_j = Utils::diracDelta(distX) * Utils::diracDelta(distY);

						// Add to A matrix
						A[i * dim + j] += diracVal_i * diracVal_j;
					}

					// Mulitply by volume
					A[i * dim + j] *= 1.0 * 1.0 * nodej->ds;
				}
			}

			// Set RHS
			vector<double> b(dim, 1.0);

			// Solve system
			vector<double> epsilon = Utils::solveLAPACK(A, b);

			// Set to node values
			for (size_t i = 0; i < dim; i++)
				(*iBodyPtr)[ib].node[i]->epsilon = epsilon[i];
		}
	}

	// If universal calculation then need to feed them back to iBody vector
#ifdef UNI_EPSILON

	// Loop through all bodies and nodes
	for (size_t n = 0; n < iNode.size(); n++)
		iNode[n].epsilon = (*iBodyPtr)[0].node[n]->epsilon;
#endif
}

// Read in geometry file
void ObjectsClass::geometryReadIn() {

	// First check if restart file exists
	if (!boost::filesystem::exists("input/geometry.config"))
		return;

	// Open config file
	ifstream file;
	file.open("input/geometry.config");

	// Handle failure to open
	if (!file.is_open())
		ERROR("Error opening geometry configuration file...exiting");

	// Skip comment lines in config file
	streamoff fileOffset;
	string line;
	file.seekg(ios::beg);
	do {

		// Get the current position within the file
		fileOffset = file.tellg();

		// Get the whole line
		getline(file, line);

	} while (line[0] == '#' && !file.eof());

	// Reset file position to the start of the last read line
	file.seekg(fileOffset, ios::beg);

	// Type of case (the first entry on each line is the keyword describing the body case)
	string bodyCase;

	// Do a quick scan to see how many bodies there are for reserving memory
	int bodyCount = 0, nodeCount = 0;
	while (file >> bodyCase) {

		// Get number of bodies
		int nBodies; file >> nBodies;

		// Fast forward to dimensions of body
		string dummy;
		for (int i = 0; i < 5; i++)
			file >> dummy;

		// Read in specific values
		double dim; file >> dim;

		// ** CIRCLE ** //
		if (bodyCase == "CIRCLE") {
			bodyCount += nBodies;
			nodeCount += nBodies * static_cast<int>(floor(2.0 * M_PI * dim / gPtr->Dx));
		}
		// ** FILAMENT ** //
		else if (bodyCase == "FILAMENT") {
			bodyCount += nBodies;
			nodeCount += nBodies * (static_cast<int>(floor(dim / gPtr->Dx)) + 1);
		}

		// Skip to end of line
		file.ignore(numeric_limits<streamsize>::max(), '\n');
	}

	// Reserve the space in iBody
	iBody.reserve(bodyCount);
	iNode.reserve(nodeCount);

	// Clear the end of file error within the ifstream
	file.clear();

	// Reset to start of geometry configurations
	file.seekg(fileOffset, ios::beg);

	// Start reading in config file
	while(file) {

		// Get type of body
		file >> bodyCase;

		// Read in general values first
		int number; file >> number;
		int ID; file >> ID;
		array<double, dims> start; file >> start[eX]; file >> start[eY];
		array<double, dims> space; file >> space[eX]; file >> space[eY];

		// ** CIRCLE ** //
		if (bodyCase == "CIRCLE") {

			// Read in specific values
			double radius; file >> radius;

			// Loop through and build
			for (int i = 0; i < number; i++) {

				// Get position vector
				array<double, dims> pos = {start[eX] + i * space[eX], start[eY] + i * space[eY]};

				// Call constructor to build it
				iBody.emplace_back(this, ID, pos, radius);

				// Increment body ID
				ID++;
			}
		}

		// ** FILAMENT ** //
		else if (bodyCase == "FILAMENT") {

			// Read in specific values
			array<double, dims> geom; file >> geom[eX]; file >> geom[eY];
			double angle; file >> angle;
			string flex; file >> flex;
			string nElements; file >> nElements;
			string BC; file >> BC;
			double rho; file >> rho;
			double E; file >> E;

			// Loop through and build
			for (int i = 0; i < number; i++) {

				// Get position vector
				array<double, dims> pos = {start[eX] + i * space[eX], start[eY] + i * space[eY]};

				// Call constructor to build it
				iBody.emplace_back(this, ID, pos, geom, angle, flex, nElements, BC, rho, E);

				// Increment body ID
				ID++;
			}
		}

		// Set case to none
		bodyCase = "NONE";
	}

	// Check there is no duplicate IDs
	for (size_t i = 0; i < iBody.size(); i++) {
		for (size_t j = i + 1; j < iBody.size(); j++) {
			if (iBody[i].ID == iBody[j].ID)
				ERROR("Duplicate body IDs in geometry.config...exiting");
		}
	}

	// Remove overlapping markres
	removeOverlapMarkers();

	// Set IBM body flag
	if (iBody.size() > 0)
		hasIBM = true;

	// Loop through all bodies
	for (size_t ib = 0; ib < iBody.size(); ib++) {

		// Set flag
		if (iBody[ib].flex == eFlexible)
			hasFlex = true;
	}

	// Print number of bodies
	if (hasIBM == true)
		cout << "found " << iBody.size() << " object(s)";
	else if (hasIBM == false)
		cout << "no bodies found";
}

// Write out forces on bodies
void ObjectsClass::writeTotalForces() {

	// Only write if there are IBM bodies
	if (hasIBM == true) {

		// File name
		string fname = "Results/TotalForces.out";

		// Check if file already exists
		bool existing = false;
		if (boost::filesystem::exists(fname))
			existing = true;

		// Open the file
		ofstream output;
		output.open(fname.c_str(), ios::app);
		output.precision(PRECISION);

		// Handle failure to open
		if (!output.is_open())
			ERROR("Error opening forces file...exiting");

		// Write out header
		if (existing == false)
			output << "Timestep\tTime (s)\tFx (N)\tFy (N)\tCx\tCy";

		// Get force scaling
		double forceScale = gPtr->Dm * gPtr->Dx / SQ(gPtr->Dt) * 1.0 / gPtr->Dx;

		// Force vector
		array<double, dims> force = {0.0};

		// Get total force
		for (size_t n = 0; n < iNode.size(); n++) {
			for (int d = 0; d < dims; d++) {
				force[d] -= iNode[n].force[d] * 1.0 * iNode[n].epsilon * iNode[n].ds * forceScale;
			}
		}

		// Get ND force
		double ND = 0.5 * ref_rho * SQ(ref_U) * ref_L;

		// Write out
		output << endl << gPtr->t << "\t" << gPtr->Dt * gPtr->t << "\t" << force[eX] << "\t" << force[eY] << "\t" << force[eX] / ND << "\t" << force[eY] / ND;

		// Close file
		output.close();
	}
}

// Write out tip positions
void ObjectsClass::writeTips() {

	// Only write if there are FEM bodies
	if (hasFlex == true) {

		// File names
		string fnamePos = "Results/TipPositions.out";
		string fnameVel = "Results/TipVelocities.out";

		// Check if files already exist
		bool existingPos = false;
		bool existingVel = false;
		if (boost::filesystem::exists(fnamePos))
			existingPos = true;
		if (boost::filesystem::exists(fnameVel))
			existingVel = true;

		// Open the file
		ofstream outputPos, outputVel;
		outputPos.open(fnamePos.c_str(), ios::app);
		outputVel.open(fnameVel.c_str(), ios::app);
		outputPos.precision(PRECISION);
		outputVel.precision(PRECISION);

		// Handle failure to open
		if (!outputPos.is_open() || !outputVel.is_open())
			ERROR("Error opening tip positions/velocities files...exiting");

		// Write out header
		if (existingPos == false)
			outputPos << "Timestep\tTime (s)\tTipX\tTipY";
		if (existingVel == false)
			outputVel << "Timestep\tTime (s)\tTipX\tTipY";

		// Write out
		outputPos << endl << gPtr->t << "\t" << gPtr->Dt * gPtr->t;
		outputVel << endl << gPtr->t << "\t" << gPtr->Dt * gPtr->t;

		// Now loop through all flexible bodies
		for (size_t ib = 0; ib < iBody.size(); ib++) {

			// Only do if flexible
			if (iBody[ib].flex == eFlexible) {
				outputPos << "\t" << iBody[ib].node[iBody[ib].node.size()-1]->pos[eX] << "\t" << iBody[ib].node[iBody[ib].node.size()-1]->pos[eY];
				outputVel << "\t" << iBody[ib].node[iBody[ib].node.size()-1]->vel[eX] << "\t" << iBody[ib].node[iBody[ib].node.size()-1]->vel[eY];
			}
		}

		// Close file
		outputPos.close();
		outputVel.close();
	}
}

// Write info data at tInfo frequency
void ObjectsClass::writeInfo() {

	// Only write if there are FEM bodies
	if (hasFlex == true) {

		// Average values
		int count = 0;
		int aveIt = 0;
		double aveRes = 0.0;
		double aveDisp = 0.0;
		array<array<double, dims>, dims> T;
		array<double, dims> deflectVector;

		// Loop through all FEM bodies and get average iterations, residual and deflection
		for (size_t ib = 0; ib < iBody.size(); ib++) {

			// Only do if body is flexible
			if (iBody[ib].flex == eFlexible) {

				// Increment count
				count++;

				// Sum to get average FEM loop parameters
				aveIt += iBody[ib].sBody->itNR;
				aveRes += iBody[ib].sBody->resNR;

				// Get rotation matrix
				T = Utils::getRotationMatrix(iBody[ib].sBody->angle0);

				// Get deflection vector and rotate into local coordinates
				int idx = iBody[ib].sBody->bodyDOFs - nodeDOFs;
				deflectVector = {iBody[ib].sBody->U[idx], iBody[ib].sBody->U[idx+1]};
				deflectVector = Utils::Transpose(T) * deflectVector;

				// Sum to get average displacement
				aveDisp -= deflectVector[eY] / iBody[ib].sBody->L0;
			}
		}

		// Get average
		aveIt /= count;
		aveRes /= static_cast<double>(count);
		aveDisp /= static_cast<double>(count);

		// Write out sub iterations
		cout << endl;
		cout << "Current relaxation factor = " << relax << endl;
		cout << "Subiterations taking " << subIt << " iterations to reach a residual of " << subRes << endl;

		// Write out average values
		cout << "FEM taking " << aveIt << " iterations to reach a residual of " << aveRes << " on average" << endl;
		cout << setprecision(2) << "Average tip deflection = " << 100.0 * aveDisp << " (% of L)";
	}
}

// Write info data at tInfo frequency
void ObjectsClass::writeLog() {

	// Only write if we have IBM bodies
	if (hasIBM == true) {

		// Test to see if over tCrit is met
		bool tCritTest = false;

		// Create stream for output to stdout and file
		Utils::io_ofstream output("Results/Log.out", Utils::io_ofstream::app);

		// Loop through all bodies
		for (size_t ib = 0; ib < iBody.size(); ib++) {

			// BODY HEADER
			output << "\n\nBODY " << iBody[ib].ID;

			// Write out body type
			if (iBody[ib].bodyType == eCircle)
				output << " (CIRCLE):\n";
			else if (iBody[ib].bodyType == eFilament)
				output << " (FILAMENT):\n";

			// Write out number of nodes
			output << "IBM Nodes = " << iBody[ib].node.size() << "\n";

			// Write out flexible type
			if (iBody[ib].flex == eFlexible)
				output << "FlexType = FLEXIBLE";
			if (iBody[ib].flex == eRigid)
				output << "FlexType = RIGID";

			// If flexible then write out some more information
			if (iBody[ib].flex == eFlexible) {

				// Write out number of elements
				output << "\n";
				output << "Number of elements = " << iBody[ib].sBody->element.size() << "\n";
				output << "Density (kg/m^3) = " << iBody[ib].sBody->element[0].rho << "\n";
				output << "Young Modulus (Pa) = " << iBody[ib].sBody->element[0].E << "\n";
				output << "tCrit (s) = " << iBody[ib].sBody->tCrit;

				// Do a check to see if tCrit has been met
				if (iBody[ib].sBody->tCrit < gPtr->Dt)
					tCritTest = true;
			}
		}

		// If tCrit has been met then write warning
		if (tCritTest == true)
			WARN("SOME TCRITS ARE SMALLER THAN TIMESTEP");

		// Close
		output.close();
	}
}

// Write VTK (IBM or FEM)
void ObjectsClass::writeVTK(bool writeIBM) {

	// Only write if we have IBM/FEM bodies to write
	if ((writeIBM && hasIBM == true) || (!writeIBM && hasFlex == true)) {

		// Get the endianness
		string endianStr = (gPtr->bigEndian ? "BigEndian" : "LittleEndian");

		// Create file
		string fStr = (writeIBM ? "IBM" : "FEM");
		ofstream output;
		output.open("Results/VTK/" + fStr + "." + to_string(gPtr->t) + ".vtp", ios::binary);

		// Handle failure to open
		if (!output.is_open())
			ERROR("Error opening body VTK file...exiting");

		// Write XML header
		output << "<?xml version=\"1.0\"?>\n";

		// Begin VTK file
		int level = 0;
		output << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"" << endianStr << "\" header_type=\"UInt64\">\n";

		// New level -> PolyData
		level = 1;
		output << string(level, '\t') << "<PolyData>\n";

		// Now loop through each IBM body
		unsigned long long offset = 0;
		for (size_t ib = 0; ib < iBody.size(); ib++) {

			// If FEM then only write flexible bodies
			if (writeIBM || (!writeIBM && iBody[ib].flex == eFlexible)) {

				// Get number of nodes and lines
				size_t nNodes = (writeIBM ? iBody[ib].node.size() : iBody[ib].sBody->node.size());
				size_t nLines = (iBody[ib].bodyType == eCircle ? nNodes : nNodes - 1);

				// New level -> Piece
				level = 2;
				output << string(level, '\t') << "<Piece "
											  << "NumberOfPoints=\"" << nNodes << "\" "
											  << "NumberOfVerts=\"" << 0 << "\" "
											  << "NumberOfLines=\"" << nLines << "\" "
											  << "NumberOfStrips=\"" << 0 << "\" "
											  << "NumberOfPolys=\"" << 0 << "\">\n";

				// New level -> Points
				level = 3;
				output << string(level, '\t') << "<Points>\n";

				// New level -> Points data (increment offset)
				level = 4;
				output << string(level, '\t') << "<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\"/>\n";
				offset += 3 * nNodes * sizeof(double) + sizeof(unsigned long long);

				// Back level -> Points
				level = 3;
				output << string(level, '\t') << "</Points>\n";

				// New level -> Lines
				level = 3;
				output << string(level, '\t') << "<Lines>\n";

				// New level -> connectivity and offsets
				level = 4;
				output << string(level, '\t') << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\" offset=\"" << offset << "\"/>\n";
				offset += 2 * nLines * sizeof(long long) + sizeof(unsigned long long);
				output << string(level, '\t') << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\"" << offset << "\"/>\n";
				offset += nLines * sizeof(long long) + sizeof(unsigned long long);

				// Back level -> Lines
				level = 3;
				output << string(level, '\t') << "</Lines>\n";

				// Back level -> Piece
				level = 2;
				output << string(level, '\t') << "</Piece>\n";
			}
		}

		// Back level -> PolyData
		level = 1;
		output << string(level, '\t') << "</PolyData>\n";

		// New level -> AppendedData
		level = 1;
		output << string(level, '\t') << "<AppendedData encoding=\"raw\">\n";

		// New level -> Raw data
		level = 2;
		output << string(level, '\t') << "_";

		// Now loop through each IBM body
		for (size_t ib = 0; ib < iBody.size(); ib++) {

			// If FEM then only write flexible bodies
			if (writeIBM || (!writeIBM && iBody[ib].flex == eFlexible)) {

				// Get number of nodes and lines
				size_t nNodes = (writeIBM ? iBody[ib].node.size() : iBody[ib].sBody->node.size());
				size_t nLines = (iBody[ib].bodyType == eCircle ? nNodes : nNodes - 1);

				// Positions
				unsigned long long size = 3 * nNodes * sizeof(double);
				output.write((char*)&size, sizeof(unsigned long long));
				for (size_t n = 0; n < nNodes; n++) {
					double posX = (writeIBM ? iBody[ib].node[n]->pos[eX] : iBody[ib].sBody->node[n].pos[eX]);
					double posY = (writeIBM ? iBody[ib].node[n]->pos[eY] : iBody[ib].sBody->node[n].pos[eY]);
					double posZ = 0.0;
					output.write((char*)&posX, sizeof(double));
					output.write((char*)&posY, sizeof(double));
					output.write((char*)&posZ, sizeof(double));
				}

				// Connectivity
				size = 2 * nLines * sizeof(long long);
				output.write((char*)&size, sizeof(unsigned long long));
				for (size_t n = 0; n < nLines; n++) {
					long long node1 = n;
					long long node2 = (n + 1) % nNodes;
					output.write((char*)&node1, sizeof(long long));
					output.write((char*)&node2, sizeof(long long));
				}

				// Offsets
				size = nLines * sizeof(long long);
				output.write((char*)&size, sizeof(unsigned long long));
				for (size_t n = 0; n < nLines; n++) {
					long long offsets = 2 * (n + 1);
					output.write((char*)&offsets, sizeof(long long));
				}
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

		// If FEM writing on then write
#ifdef VTK_FEM
		if(writeIBM)
			writeVTK(false);
#endif
	}
}

// Remove markers that are too close to each other
void ObjectsClass::removeOverlapMarkers() {

	// IDs of nodes to erase
	vector<int> eraseID;

	// Loop through all IBM bodies and nodes
	for (size_t n1 = 0; n1 < iNode.size(); n1++) {

		// Loop through all nodes ahead in the array
		for (size_t n2 = n1 + 1; n2 < iNode.size(); n2++) {

			// Get distance between nodes
			double mag = sqrt((iNode[n1].pos - iNode[n2].pos) * (iNode[n1].pos - iNode[n2].pos)) / gPtr->Dx;

			// If distance is smaller than 0.5Dx then delete
			if (mag < 0.5) {

				// If both bodies are flexible then error
				if ((iNode[n1].iPtr->flex == eFlexible) && (iNode[n2].iPtr->flex == eFlexible))
					ERROR("Cannot place two flexible bodies so close to each other...exiting");

				// If body1 is flex and body2 is rigid then remove from body2
				if ((iNode[n1].iPtr->flex == eFlexible) && (iNode[n2].iPtr->flex == eRigid))
					eraseID.push_back(iNode[n2].ID);

				// If body1 is rigid and body2 is flex then remove from body1
				if ((iNode[n1].iPtr->flex == eRigid) && (iNode[n2].iPtr->flex == eFlexible))
					eraseID.push_back(iNode[n1].ID);

				// If both bodies are rigid then remove from body1 (there is probably a more elegant way to choose)
				if ((iNode[n1].iPtr->flex == eRigid) && (iNode[n2].iPtr->flex == eRigid))
					eraseID.push_back(iNode[n1].ID);
			}
		}
	}

	// Now loop through and delete
	for (size_t i = 0; i < eraseID.size(); i++) {

		// First delete pointer stored in iBody
		bool foundMatch = false;
		for (size_t ib = 0; ib < iBody.size(); ib++) {
			for (size_t n = 0; n < iBody[ib].node.size(); n++) {

				// Check ID
				if (iBody[ib].node[n]->ID == eraseID[i]) {
					iBody[ib].node.erase(iBody[ib].node.begin() + n);
					foundMatch = true;
				}

				// Move rest of pointers back
				if (foundMatch)
					iBody[ib].node[n]--;
			}
		}

		// Now delete the actual node
		for (size_t n = 0; n < iNode.size(); n++) {

			// Check ID
			if (iNode[n].ID == eraseID[i]) {
				iNode.erase(iNode.begin() + n);
				break;
			}
		}
	}

	// Now reset the IDs of all nodes
	int ID = 0;
	for (size_t n = 0; n < iNode.size(); n++)
		iNode[n].ID = ID++;
}

// Initialise objects
void ObjectsClass::initialiseObjects() {

	// Call static FEM to give it an initial deflection
#ifdef INITIAL_DEFLECT
	if (hasFlex == true)
		initialDeflect();
#endif

	// Loop through all nodes
	for (size_t n = 0; n < iNode.size(); n++) {

		// Find support
		iNode[n].findSupport();

		// Compute ds
		iNode[n].computeDs();
	}

	// Compute epsilon
	computeEpsilon();
}

// Call static FEM to give initial deflection
void ObjectsClass::initialDeflect() {

	// Starting main algorithm
	cout << endl << endl << endl << "*** INITIALISING DEFLECTIONS ***" << endl;

	// Loop through all bodies
	for (size_t ib = 0; ib < iBody.size(); ib++) {

		// Only do if flexible
		if (iBody[ib].flex == eFlexible)
			iBody[ib].sBody->setInitialDeflection();
	}
}

// Read in restart file
void ObjectsClass::readRestart() {

	// ** IBM FILE ** //
	if (boost::filesystem::exists("Results/Restart/IBM.restart")) {

		// Get flag for endianess
		bool bigEndian = gPtr->bigEndian;

		// Open file
		ifstream file;
		file.open("Results/Restart/IBM.restart", ios::binary);

		// Handle failure to open
		if (!file.is_open())
			ERROR("Error opening IBM.restart file...exiting");

		// Read in number of IBM bodies
		size_t nIBMRead;
		file.read((char*)&nIBMRead, sizeof(size_t));
		size_t nIBMSwap = (bigEndian ? Utils::swapEnd(nIBMRead) : nIBMRead);

		// Loop through lines in IBM.restart file
		for (size_t i = 0; i < nIBMSwap; i++) {

			// Read in body ID and number of nodes
			int idRead;
			size_t nodeSizeRead;
			file.read((char*)&idRead, sizeof(int));
			file.read((char*)&nodeSizeRead, sizeof(size_t));

			// Swap byte order if bigEndian
			int idSwap = (bigEndian ? Utils::swapEnd(idRead) : idRead);
			size_t nodeSizeSwap = (bigEndian ? Utils::swapEnd(nodeSizeRead) : nodeSizeRead);

			// Check to make sure it exists in this run
			int ib = getBodyIdxFromID(idSwap);

			// If it doesn't exist then skip the rest of the data for this body
			if (ib < 0) {
				size_t shift = (6 * sizeof(double) * nodeSizeSwap);
				file.seekg (shift, ios::cur);
				continue;
			}

			// Check to make sure number of nodes match up
			if (nodeSizeSwap != iBody[ib].node.size())
				ERROR("Number of nodes in IBM.restart do not match existing number of nodes...exiting");

			// Read in rest of the line
			for (size_t n = 0; n < nodeSizeSwap; n++) {

				// Read in positions and velocities
				double posxRead, posyRead, velxRead, velyRead, forcexRead, forceyRead;
				file.read((char*)&posxRead, sizeof(double));
				file.read((char*)&posyRead, sizeof(double));
				file.read((char*)&velxRead, sizeof(double));
				file.read((char*)&velyRead, sizeof(double));
				file.read((char*)&forcexRead, sizeof(double));
				file.read((char*)&forceyRead, sizeof(double));

				// Swap byte order if bigEndian
				iBody[ib].node[n]->pos[eX] = (bigEndian ? Utils::swapEnd(posxRead) : posxRead);
				iBody[ib].node[n]->pos[eY] = (bigEndian ? Utils::swapEnd(posyRead) : posyRead);
				iBody[ib].node[n]->vel[eX] = (bigEndian ? Utils::swapEnd(velxRead) : velxRead);
				iBody[ib].node[n]->vel[eY] = (bigEndian ? Utils::swapEnd(velyRead) : velyRead);
				iBody[ib].node[n]->force[eX] = (bigEndian ? Utils::swapEnd(forcexRead) : forcexRead);
				iBody[ib].node[n]->force[eY] = (bigEndian ? Utils::swapEnd(forceyRead) : forceyRead);
			}
		}
	}

	// ** FEM FILE ** //
	if (boost::filesystem::exists("Results/Restart/FEM.restart")) {

		// Get flag for endianess
		bool bigEndian = gPtr->bigEndian;

		// Open file
		ifstream file;
		file.open("Results/Restart/FEM.restart", ios::binary);

		// Handle failure to open
		if (!file.is_open())
			ERROR("Error opening FEM.restart file...exiting");

		// Read in body ID and number of nodes
		int nFEMRead;
		double relaxRead;
		file.read((char*)&nFEMRead, sizeof(int));
		file.read((char*)&relaxRead, sizeof(double));

		// Swap byte order if bigEndian
		int nFEMSwap = (bigEndian ? Utils::swapEnd(nFEMRead) : nFEMRead);
		relax = (bigEndian ? Utils::swapEnd(relaxRead) : relaxRead);

		// Loop through lines in FEM.restart file
		for (int i = 0; i < nFEMSwap; i++) {

			// Read in body ID and number of nodes
			int idRead;
			size_t elSizeRead, nodeSizeRead;
			file.read((char*)&idRead, sizeof(int));
			file.read((char*)&elSizeRead, sizeof(size_t));
			file.read((char*)&nodeSizeRead, sizeof(size_t));

			// Swap byte order if bigEndian
			int idSwap = (bigEndian ? Utils::swapEnd(idRead) : idRead);
			size_t elSizeSwap = (bigEndian ? Utils::swapEnd(elSizeRead) : elSizeRead);
			size_t nodeSizeSwap = (bigEndian ? Utils::swapEnd(nodeSizeRead) : nodeSizeRead);

			// Check to make sure it exists in this run
			int ib = getBodyIdxFromID(idSwap);

			// If it doesn't exist or if it is now rigid then skip this line
			if (ib < 0 || iBody[ib].flex == eRigid) {
				size_t shift = (6 * sizeof(double) * nodeSizeSwap * nodeDOFs + 2 * sizeof(double) * elSizeSwap + 6 * sizeof(double) * nodeSizeSwap);
				file.seekg (shift, ios::cur);
				continue;
			}

			// Check to make sure number of nodes match up
			if (elSizeSwap != iBody[ib].sBody->element.size() || nodeSizeSwap != iBody[ib].sBody->node.size())
				ERROR("Number of nodes/elements in FEM.restart do not match existing number of nodes/elements...exiting");

			// Read in U
			for (int i = 0; i < iBody[ib].sBody->bodyDOFs; i++) {
				double valRead;
				file.read((char*)&valRead, sizeof(double));
				iBody[ib].sBody->U[i] = (bigEndian ? Utils::swapEnd(valRead) : valRead);
			}

			// Read in Udot
			for (int i = 0; i < iBody[ib].sBody->bodyDOFs; i++) {
				double valRead;
				file.read((char*)&valRead, sizeof(double));
				iBody[ib].sBody->Udot[i] = (bigEndian ? Utils::swapEnd(valRead) : valRead);
			}

			// Read in Udotdot
			for (int i = 0; i < iBody[ib].sBody->bodyDOFs; i++) {
				double valRead;
				file.read((char*)&valRead, sizeof(double));
				iBody[ib].sBody->Udotdot[i] = (bigEndian ? Utils::swapEnd(valRead) : valRead);
			}

			// Read in U_n
			for (int i = 0; i < iBody[ib].sBody->bodyDOFs; i++) {
				double valRead;
				file.read((char*)&valRead, sizeof(double));
				iBody[ib].sBody->U_n[i] = (bigEndian ? Utils::swapEnd(valRead) : valRead);
			}

			// Read in U_nm1
			for (int i = 0; i < iBody[ib].sBody->bodyDOFs; i++) {
				double valRead;
				file.read((char*)&valRead, sizeof(double));
				iBody[ib].sBody->U_nm1[i] = (bigEndian ? Utils::swapEnd(valRead) : valRead);
			}

			// Read in R_k
			for (int i = 0; i < iBody[ib].sBody->bodyDOFs; i++) {
				double valRead;
				file.read((char*)&valRead, sizeof(double));
				iBody[ib].sBody->R_k[i] = (bigEndian ? Utils::swapEnd(valRead) : valRead);
			}

			// Read in elements
			for (size_t el = 0; el < elSizeSwap; el++) {
				double lengthRead, angleRead;
				file.read((char*)&lengthRead, sizeof(double));
				file.read((char*)&angleRead, sizeof(double));
				iBody[ib].sBody->element[el].L = (bigEndian ? Utils::swapEnd(lengthRead) : lengthRead);
				iBody[ib].sBody->element[el].angle = (bigEndian ? Utils::swapEnd(angleRead) : angleRead);
			}

			// Read in nodes
			for (size_t n = 0; n < nodeSizeSwap; n++) {
				double pos0xRead, pos0yRead, posxRead, posyRead, angle0Read, angleRead;
				file.read((char*)&pos0xRead, sizeof(double));
				file.read((char*)&pos0yRead, sizeof(double));
				file.read((char*)&posxRead, sizeof(double));
				file.read((char*)&posyRead, sizeof(double));
				file.read((char*)&angle0Read, sizeof(double));
				file.read((char*)&angleRead, sizeof(double));
				iBody[ib].sBody->node[n].pos0[eX] = (bigEndian ? Utils::swapEnd(pos0xRead) : pos0xRead);
				iBody[ib].sBody->node[n].pos0[eY] = (bigEndian ? Utils::swapEnd(pos0yRead) : pos0yRead);
				iBody[ib].sBody->node[n].pos[eX] = (bigEndian ? Utils::swapEnd(posxRead) : posxRead);
				iBody[ib].sBody->node[n].pos[eY] = (bigEndian ? Utils::swapEnd(posyRead) : posyRead);
				iBody[ib].sBody->node[n].angle0 = (bigEndian ? Utils::swapEnd(angle0Read) : angle0Read);
				iBody[ib].sBody->node[n].angle = (bigEndian ? Utils::swapEnd(angleRead) : angleRead);
			}
		}
	}

	// Loop through all nodes
	for (size_t n = 0; n < iNode.size(); n++) {

		// Find support
		iNode[n].findSupport();

		// Compute ds
		iNode[n].computeDs();
	}

	// Compute epsilon
	if (hasIBM)
		computeEpsilon();

	// If flexible then recalculate FEM values too
	for (size_t ib = 0; ib < iBody.size(); ib++) {
		if (iBody[ib].flex == eFlexible) {

			// Loop through elements
			for (size_t el = 0; el < iBody[ib].sBody->element.size(); el++) {

				// Recalculate values
				iBody[ib].sBody->element[el].setElementTransform();
				iBody[ib].sBody->element[el].setLocalMatrices();
			}
		}
	}

	// Clean up files from previous runs
	restartCleanup();
}

// Write out restart file
void ObjectsClass::writeRestart() {

	// ** IBM FILE ** //
	if (hasIBM == true) {

		// Get flag for endianess
		bool bigEndian = gPtr->bigEndian;

		// Open file
		ofstream outputIBM;
		outputIBM.open("Results/Restart/IBM.restart", ios::binary);

		// Handle failure to open
		if (!outputIBM.is_open())
			ERROR("Error opening IBM.restart file...exiting");

		// Swap byte order if bigEndian
		size_t ibSizeWrite = (bigEndian ? Utils::swapEnd(iBody.size()) : iBody.size());

		// Write out number of bodies
		outputIBM.write((char*)&ibSizeWrite, sizeof(size_t));

		// Loop through IBM bodies
		for (size_t ib = 0; ib < iBody.size(); ib++) {

			// Swap byte order if bigEndian
			int idWrite = (bigEndian ? Utils::swapEnd(iBody[ib].ID) : iBody[ib].ID);
			size_t nodeSizeWrite = (bigEndian ? Utils::swapEnd(iBody[ib].node.size()) : iBody[ib].node.size());

			// Write out body ID and node
			outputIBM.write((char*)&idWrite, sizeof(int));
			outputIBM.write((char*)&nodeSizeWrite, sizeof(size_t));

			// Loop through nodes and write out ID, pos, pos0, vel
			for (size_t n = 0; n < iBody[ib].node.size(); n++) {

				// Swap byte order if bigEndian
				double posxWrite = (bigEndian ? Utils::swapEnd(iBody[ib].node[n]->pos[eX]) : iBody[ib].node[n]->pos[eX]);
				double posyWrite = (bigEndian ? Utils::swapEnd(iBody[ib].node[n]->pos[eY]) : iBody[ib].node[n]->pos[eY]);
				double velxWrite = (bigEndian ? Utils::swapEnd(iBody[ib].node[n]->vel[eX]) : iBody[ib].node[n]->vel[eX]);
				double velyWrite = (bigEndian ? Utils::swapEnd(iBody[ib].node[n]->vel[eY]) : iBody[ib].node[n]->vel[eY]);
				double forcexWrite = (bigEndian ? Utils::swapEnd(iBody[ib].node[n]->force[eX]) : iBody[ib].node[n]->force[eX]);
				double forceyWrite = (bigEndian ? Utils::swapEnd(iBody[ib].node[n]->force[eY]) : iBody[ib].node[n]->force[eY]);

				// Write out number of bodies
				outputIBM.write((char*)&posxWrite, sizeof(double));
				outputIBM.write((char*)&posyWrite, sizeof(double));
				outputIBM.write((char*)&velxWrite, sizeof(double));
				outputIBM.write((char*)&velyWrite, sizeof(double));
				outputIBM.write((char*)&forcexWrite, sizeof(double));
				outputIBM.write((char*)&forceyWrite, sizeof(double));
			}
		}

		// Close file
		outputIBM.close();
	}

	// ** FEM FILE ** //
	if (hasFlex == true) {

		// Get flag for endianess
		bool bigEndian = gPtr->bigEndian;

		// Open file
		ofstream outputFEM;
		outputFEM.open("Results/Restart/FEM.restart", ios::binary);

		// Handle failure to open
		if (!outputFEM.is_open())
			ERROR("Error opening FEM.restart file...exiting");

		// Swap byte order if bigEndian
		int flexWrite = (bigEndian ? Utils::swapEnd(nFlex) : nFlex);
		double relaxWrite = (bigEndian ? Utils::swapEnd(relax) : relax);

		// Write out global information
		outputFEM.write((char*)&flexWrite, sizeof(int));
		outputFEM.write((char*)&relaxWrite, sizeof(double));

		// Loop through IBM bodies
		for (size_t ib = 0; ib < iBody.size(); ib++) {
			if (iBody[ib].flex == eFlexible) {

				// Swap byte order if bigEndian
				int idWrite = (bigEndian ? Utils::swapEnd(iBody[ib].ID) : iBody[ib].ID);
				size_t elSizeWrite = (bigEndian ? Utils::swapEnd(iBody[ib].sBody->element.size()) : iBody[ib].sBody->element.size());
				size_t nodeSizeWrite = (bigEndian ? Utils::swapEnd(iBody[ib].sBody->node.size()) : iBody[ib].sBody->node.size());

				// Write out body ID and node
				outputFEM.write((char*)&idWrite, sizeof(int));
				outputFEM.write((char*)&elSizeWrite, sizeof(size_t));
				outputFEM.write((char*)&nodeSizeWrite, sizeof(size_t));

				// Write out U
				for (int i = 0; i < iBody[ib].sBody->bodyDOFs; i++) {
					double valWrite = (bigEndian ? Utils::swapEnd(iBody[ib].sBody->U[i]) : iBody[ib].sBody->U[i]);
					outputFEM.write((char*)&valWrite, sizeof(double));
				}

				// Write out Udot
				for (int i = 0; i < iBody[ib].sBody->bodyDOFs; i++) {
					double valWrite = (bigEndian ? Utils::swapEnd(iBody[ib].sBody->Udot[i]) : iBody[ib].sBody->Udot[i]);
					outputFEM.write((char*)&valWrite, sizeof(double));
				}

				// Write out Udotdot
				for (int i = 0; i < iBody[ib].sBody->bodyDOFs; i++) {
					double valWrite = (bigEndian ? Utils::swapEnd(iBody[ib].sBody->Udotdot[i]) : iBody[ib].sBody->Udotdot[i]);
					outputFEM.write((char*)&valWrite, sizeof(double));
				}

				// Write out U_n
				for (int i = 0; i < iBody[ib].sBody->bodyDOFs; i++) {
					double valWrite = (bigEndian ? Utils::swapEnd(iBody[ib].sBody->U_n[i]) : iBody[ib].sBody->U_n[i]);
					outputFEM.write((char*)&valWrite, sizeof(double));
				}

				// Write out U_nm1
				for (int i = 0; i < iBody[ib].sBody->bodyDOFs; i++) {
					double valWrite = (bigEndian ? Utils::swapEnd(iBody[ib].sBody->U_nm1[i]) : iBody[ib].sBody->U_nm1[i]);
					outputFEM.write((char*)&valWrite, sizeof(double));
				}

				// Write out R_k
				for (int i = 0; i < iBody[ib].sBody->bodyDOFs; i++) {
					double valWrite = (bigEndian ? Utils::swapEnd(iBody[ib].sBody->R_k[i]) : iBody[ib].sBody->R_k[i]);
					outputFEM.write((char*)&valWrite, sizeof(double));
				}

				// Write out elements
				for (size_t el = 0; el < iBody[ib].sBody->element.size(); el++) {
					double lengthWrite = (bigEndian ? Utils::swapEnd(iBody[ib].sBody->element[el].L) : iBody[ib].sBody->element[el].L);
					double angleWrite = (bigEndian ? Utils::swapEnd(iBody[ib].sBody->element[el].angle) : iBody[ib].sBody->element[el].angle);
					outputFEM.write((char*)&lengthWrite, sizeof(double));
					outputFEM.write((char*)&angleWrite, sizeof(double));
				}

				// Write out nodes
				for (size_t n = 0; n < iBody[ib].sBody->node.size(); n++) {
					double pos0xWrite = (bigEndian ? Utils::swapEnd(iBody[ib].sBody->node[n].pos0[eX]) : iBody[ib].sBody->node[n].pos0[eX]);
					double pos0yWrite = (bigEndian ? Utils::swapEnd(iBody[ib].sBody->node[n].pos0[eY]) : iBody[ib].sBody->node[n].pos0[eY]);
					double posxWrite = (bigEndian ? Utils::swapEnd(iBody[ib].sBody->node[n].pos[eX]) : iBody[ib].sBody->node[n].pos[eX]);
					double posyWrite = (bigEndian ? Utils::swapEnd(iBody[ib].sBody->node[n].pos[eY]) : iBody[ib].sBody->node[n].pos[eY]);
					double angle0Write = (bigEndian ? Utils::swapEnd(iBody[ib].sBody->node[n].angle0) : iBody[ib].sBody->node[n].angle0);
					double angleWrite = (bigEndian ? Utils::swapEnd(iBody[ib].sBody->node[n].angle) : iBody[ib].sBody->node[n].angle);
					outputFEM.write((char*)&pos0xWrite, sizeof(double));
					outputFEM.write((char*)&pos0yWrite, sizeof(double));
					outputFEM.write((char*)&posxWrite, sizeof(double));
					outputFEM.write((char*)&posyWrite, sizeof(double));
					outputFEM.write((char*)&angle0Write, sizeof(double));
					outputFEM.write((char*)&angleWrite, sizeof(double));
				}
			}
		}

		// Close file
		outputFEM.close();
	}
}

// Get body index from the ID
int ObjectsClass::getBodyIdxFromID(int id) {

	// Loop through all bodies
	for (size_t i = 0; i < iBody.size(); i++) {
		if (iBody[i].ID == id)
			return static_cast<int>(i);
	}

	// If we get here then it doesn't exist
	return -1;
}

// Clean up files from previous simulations before restart
void ObjectsClass::restartCleanup() {

	// If there are fluid VTKs with no corresponding body VTK then write some empty ones
	fillEmptyVTK();

	// Sort out tips/forces
	deleteTipsAndForces();

	// If this part of simulation has no IBM or FEM then get rid of their old restart files
	if (!hasIBM && boost::filesystem::exists("Results/Restart/IBM.restart")) {
		if(!boost::filesystem::remove("Results/Restart/IBM.restart"))
			ERROR("Problem removing IBM.restart...exiting");;
	}
	if (!hasFlex && boost::filesystem::exists("Results/Restart/FEM.restart")) {
		if(!boost::filesystem::remove("Results/Restart/FEM.restart"))
			ERROR("Problem removing FEM.restart...exiting");;
	}
}

// Fill in empty VTK files
void ObjectsClass::fillEmptyVTK() {

	// Only write if we have IBM bodies
	if (hasIBM == true) {

		// Results path
		string path = "Results/VTK";

		// Check if it exists
		if (boost::filesystem::exists(path)) {

			// Extract string
			string rmStr = "Fluid.";
			string exStrGrid = ".vti";
			string exStrBody = ".vtp";

			// Get directory iterator
			boost::filesystem::directory_iterator endit;

			// Loop through all fluid vtk files
			for (boost::filesystem::directory_iterator it(path); it != endit; it++) {

				// Get string
				string fStr = it->path().stem().string();

				// Check to make sure it is a fluid VTK file
				if (fStr.find(rmStr) != string::npos && it->path().extension() == exStrGrid) {

					// Get time step
					string tStep = fStr.substr(rmStr.size(), fStr.size());

					// Create fname string
					string fname = path + "/IBM." + tStep + exStrBody;

					// Check if file exists
					if (!boost::filesystem::exists(fname))
						writeEmptyVTK(fname);

					// If writing FEM we need to do the same for them as well
#ifdef VTK_FEM

					// Create fname string
					fname = path + "/FEM." + tStep + exStrBody;

					// Check if file exists
					if (!boost::filesystem::exists(fname) && hasFlex == true)
						writeEmptyVTK(fname);
#endif
				}
			}
		}
	}
}

// Write empty VTK polydata
void ObjectsClass::writeEmptyVTK(string &fname) {

	// Get the endianness
	string endianStr = (gPtr->bigEndian ? "BigEndian" : "LittleEndian");

	// Create file
	ofstream output;
	output.open(fname, ios::binary);

	// Handle failure to open
	if (!output.is_open())
		ERROR("Error opening body VTK file...exiting");

	// Write XML header
	output << "<?xml version=\"1.0\"?>\n";

	// Begin VTK file
	int level = 0;
	output << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"" << endianStr << "\" header_type=\"UInt64\">\n";

	// New level -> PolyData
	level = 1;
	output << string(level, '\t') << "<PolyData>\n";

	// New level -> Piece
	level = 2;
	output << string(level, '\t') << "<Piece "
								  << "NumberOfPoints=\"" << 0 << "\" "
								  << "NumberOfVerts=\"" << 0 << "\" "
								  << "NumberOfLines=\"" << 0 << "\" "
								  << "NumberOfStrips=\"" << 0 << "\" "
								  << "NumberOfPolys=\"" << 0 << "\">\n";

	// Back level -> Piece
	level = 2;
	output << string(level, '\t') << "</Piece>\n";

	// Back level -> PolyData
	level = 1;
	output << string(level, '\t') << "</PolyData>\n";

	// Back level -> VTKFile
	level = 0;
	output << string(level, '\t') << "</VTKFile>\n";

	// Close file
	output.close();
}

// Delete the tip and force data that was written after last restart
void ObjectsClass::deleteTipsAndForces() {

	// Delete the tip positions that were written after last restart
	deleteFileOutput("Results/TipPositions.out");

	// Delete the tip velocities that were written after last restart
	deleteFileOutput("Results/TipVelocities.out");

	// Delete the forces that were written after last restart
	deleteFileOutput("Results/TotalForces.out");
}

// Delete data from output files that were written after last restart
void ObjectsClass::deleteFileOutput(string fname) {

	// Check file exists
	if (boost::filesystem::exists(fname)) {

		// Set vector
		vector<string> fileStr;

		// Open file
		ifstream file;
		file.open(fname);

		// Handle failure to open
		if (!file.is_open())
			ERROR("Error opening file for deleting future tips/forces...exiting");

		// String for holding the line
		string lineStr;

		// Get the header
		getline(file, lineStr);
		fileStr.push_back(lineStr);

		// Loop through each line in file
		while (getline(file, lineStr)) {

			// Get the timestep
			int tStep = stoi(lineStr);

			// If tStep is bigger than gPtr->t then break
			if (tStep >= gPtr->t)
				break;

			// Else we add the string to the vector
			fileStr.push_back(lineStr);
		}

		// File name
		string fout = "Results/Output.temp";

		// Open the file
		ofstream output;
		output.open(fout.c_str());

		// Handle failure to open
		if (!output.is_open())
			ERROR("Error opening temporary output file...exiting");

		// Write out header
		output << fileStr[0];

		// Loop through and write out
		for (size_t i = 1; i < fileStr.size(); i++)
			output << endl << fileStr[i];

		// Close file
		output.close();

		// Now rename the temp file
		boost::filesystem::rename(fout, fname);
	}
}

// Custom constructor
ObjectsClass::ObjectsClass(GridClass &g) {

	// Write out header
	cout << endl << endl << "Initialising objects...";

	// Set pointers
	gPtr = &g;
	gPtr->oPtr = this;

	// Initial sub iteration members
	subIt = 0;
	subRes = 0.0;
	relax = relaxMax;
	subNum = 0.0;
	subDen = 0.0;

	// Set flag to false initially
	hasIBM = false;
	hasFlex = false;

	// Read in geometry file
	geometryReadIn();

	// Initialise objects
	if (hasIBM)
		initialiseObjects();

	// Get numer of FEM DOFs in whole simulation
	nFlex = simDOFs = 0;
	for (size_t ib = 0; ib < iBody.size(); ib++) {
		if (iBody[ib].flex == eFlexible) {
			nFlex++;
			simDOFs += iBody[ib].sBody->bodyDOFs;
		}
	}
}
