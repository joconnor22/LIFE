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

#ifndef UTILS_H // UTILS_H
#define UTILS_H

// Includes
#include "params.h"

// Forward declarations
class GridClass;

// Utility namespace for helper functions
namespace Utils {

// DECLARATIONS

// Write header at start of run
void writeHeader();

// Create directories at start
bool createDirectories();

// Read in restart file
void readRestart(GridClass &grid);

// Read in restart file
void writeRestart(GridClass &grid);

// Write out info
void writeInfo(GridClass &grid);

// Write out log
void writeLog(GridClass &grid);

// Write out VTK
void writeVTK(GridClass &grid);

// Delete future VTKs
void deleteVTKs(GridClass &grid);

// Get number of omp threads (as built in doesn't work on GCC)
int omp_thread_count();

// Convert seconds to hours:minutes:seconds
array<int, 3> secs2hms(double seconds);

// Get string for boundary condition
string getBoundaryString(eLatType BCType);

// Solve linear system using LAPACK routines
vector<double> solveLAPACK(vector<double> A, vector<double> b, int BC = 0);


// DEFINITIONS
// Error and exit
inline void errorExit(const string &msg) {

	// Write message then exit
	cout << endl << endl << msg << endl << endl;
	exit(99);
}

// Warning
inline void warning(const string &msg) {

	// Write message then exit
	cout << endl << endl << endl << "** WARNING: " << msg << " **";
}

// Derived struct from ofstream for writing to stdout and file at same time
struct io_ofstream : ofstream {

	// Constructor
	io_ofstream(const string& fileNameIn, openmode mode = out) : ofstream(fileNameIn, mode) , fileName(fileNameIn) {};

	// Name of file
    const string fileName;
};

// Overloaded << operator for writing to file and stdout
template <typename T>
io_ofstream& operator<<(Utils::io_ofstream& strm, const T &var) {

	// Write to stout and file
    std::cout << var;
    static_cast<ofstream&>(strm) << var;

    // Return
    return strm;
}

// Check for endianness
inline bool isBigEndian()
{
	// Create union and assign bytes
    union {uint32_t i; char c[4];} bint = {0x01020304};

    // Return fist byte
    return bint.c[0] == 1;
}

// Swap byte orders between big and little endian
template <typename T>
inline T swapEnd(const T &val) {

	// Results
	T res = val;

	// Create pointer to char arrays
	char* resArray = reinterpret_cast<char*>(&res);

	// Loop through and reorder
	for(size_t i = 0; i < sizeof(res) / 2; i++)
		swap(resArray[i], resArray[sizeof(res) - 1 - i]);

	// Return
	return res;
}

// Extrapolate value
inline double extrapolate(const vector<double> &vec, const array<int, dims> &normal, int order, int i, int j, int d = 0, int arrayDims = 1) {

	// 0th order extrapolation
	if (order == 0) {

		// Get indices for extrapolation
		int i1 = i + normal[eX];
		int j1 = j + normal[eY];

		// Return value
		return vec[(i1 * Ny + j1) * arrayDims + d];
	}

	// 1st order extrapolation
	else if (order == 1) {

		// Get indices for extrapolation
		int i1 = i + normal[eX];
		int i2 = i + 2 * normal[eX];
		int j1 = j + normal[eY];
		int j2 = j + 2 * normal[eY];

		// Return value
		return 2.0 * vec[(i1 * Ny + j1) * arrayDims + d] - vec[(i2 * Ny + j2) * arrayDims + d];
	}

	// 2st order extrapolation
	else if (order == 2) {

		// Get indices for extrapolation
		int i1 = i + normal[eX];
		int i2 = i + 2 * normal[eX];
		int i3 = i + 3 * normal[eX];
		int j1 = j + normal[eY];
		int j2 = j + 2 * normal[eY];
		int j3 = j + 3 * normal[eY];

		// Return value
		return 2.5 * vec[(i1 * Ny + j1) * arrayDims + d] - 2.0 * vec[(i2 * Ny + j2) * arrayDims + d] + 0.5 * vec[(i3 * Ny + j3) * arrayDims + d];
	}

	// Can't do other extrapolations
	else {
        ERROR("Invalid order of extrapolation...exiting");
        return 99;
	}
}

// Get value by applying a zero gradient
inline double zeroGradient(const vector<double> &vec, const array<int, dims> &normal, int order, int i, int j, int d = 0, int arrayDims = 1) {

	// 1st order extrapolation
	if (order == 1) {

		// Get indices for extrapolation
		int i1 = i + normal[eX];
		int j1 = j + normal[eY];

		// Return value
		return vec[(i1 * Ny + j1) * arrayDims + d];
	}

	// 2nd order extrapolation
	else if (order == 2) {

		// Get indices for extrapolation
		int i1 = i + normal[eX];
		int i2 = i + 2 * normal[eX];
		int j1 = j + normal[eY];
		int j2 = j + 2 * normal[eY];

		// Return value
		return (4.0 / 3.0) * vec[(i1 * Ny + j1) * arrayDims + d] - (1.0 / 3.0) * vec[(i2 * Ny + j2) * arrayDims + d];
	}

	// Can't do other extrapolations
	else {
		ERROR("Invalid order of zero gradient...exiting");
        return 99;
	}
}

// Get discretised dirac delta
inline double diracDelta(double dist) {

	// Get absolute value
    double absDist = fabs(dist);

    // Get discretised 3-point dirac delta
    if (absDist > 1.5)
    	return 0.0;
    else if (absDist > 0.5)
    	return (5.0 - 3.0 * absDist - sqrt(-3.0 * SQ(1.0 - absDist) + 1.0)) / 6.0;
    else
    	return (1.0 + sqrt(1.0 - 3.0 * SQ(absDist))) / 3.0;
}

// Shift angle into -pi -> pi range
inline double shiftAngle(double angle) {

	// Get angle
	angle = fmod(angle + M_PI, 2.0 * M_PI);

	// Check if it is below zero
	if (angle < 0.0)
		angle += 2.0 * M_PI;

	// Return new angle
	return angle - M_PI;
}

// Get sign of number
template <typename T>
inline int sgn(T val) {

	// Return sign of number
	return (T(0) < val) - (val < T(0));
}

// Get rotation matrix
template <typename T>
inline array<array<T, dims>, dims> getRotationMatrix(T angle) {

	// Declare rotation matrix
	array<array<T, dims>, dims> R = {{{cos(angle), -sin(angle)}, {sin(angle),  cos(angle)}}};

	// Return
	return R;
}

// Get transpose of matrix
template <typename T, size_t N1, size_t N2>
inline array<array<T, N1>, N2> Transpose(const array<array<T, N2>, N1> &Mat) {

	// Declare results matrix
	array<array<T, N1>, N2> ResMat;

	// Loop through and switch i's and j's
	for (size_t i = 0; i < N2; i++) {
		for (size_t j = 0; j < N1; j++)
			ResMat[i][j] = Mat[j][i];
	}

	// Return
	return ResMat;
}

// Matrix multiply (flat matrix x vector)
template <typename T>
inline vector<T> MatMultiply(const vector<T> &LMat, const vector<T> &RVec) {

	// Check they are the right size
	if (LMat.size() % RVec.size() != 0)
		ERROR("Matrix and vector not right size for multiplying...exiting");

	// Get dimensions
	size_t cols = RVec.size();
	size_t rows = LMat.size() / cols;

	// Declare res vector
	vector<T> ResVec(rows, 0.0);

	// Loop through
	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < cols; j++) {
			ResVec[i] += LMat[i * cols + j] * RVec[j];
		}
	}

	// Return
	return ResVec;
}
}


// OVERLOADED OPERATORS FOR ARRAYS
// Overloaded + operator for vectors
template <typename T, size_t N>
inline array<T, N> operator+(const array<T, N> &LVec, const array<T, N> &RVec) {

	// Declare res vector
	array<T, N> ResVec;

	// Add the left and right vectors
	for (size_t i = 0; i < N; i++)
		ResVec[i] = LVec[i] + RVec[i];

	// Return value
	return ResVec;
}

// Overloaded - operator for vectors
template <typename T, size_t N>
inline array<T, N> operator-(const array<T, N> &LVec, const array<T, N> &RVec) {

	// Declare res vector
	array<T, N> ResVec;

	// Add the left and right vectors
	for (size_t i = 0; i < N; i++)
		ResVec[i] = LVec[i] - RVec[i];

	// Return value
	return ResVec;
}

// Overloaded * operator for scalar * vector
template <typename T, size_t N>
inline array<T, N> operator*(double LScalar, const array<T, N> &RVec) {

	// Declare res vector
	array<T, N> ResVec;

	// Add the left and right vectors
	for (size_t i = 0; i < N; i++)
		ResVec[i] = LScalar * RVec[i];

	// Return value
	return ResVec;
}

// Overloaded * operator for vector * vector
template <size_t N>
inline double operator*(const array<double, N> &LVec, const array<double, N> &RVec) {

	// Declare res vector
	double res = 0.0;

	// Multiply the left and right vectors
	for (size_t i = 0; i < N; i++)
		res += LVec[i] * RVec[i];

	// Return value
	return res;
}

// Overloaded * operator for matrix * vector
template <typename T, size_t N1, size_t N2>
inline array<T, N2> operator*(const array<array<T, N1>, N2> &LMat, const array<T, N1> &RVec) {

	// Declare res vector
	array<T, N2> ResVec = {0.0};

	// Multiply the matrix and vector
	for (size_t i = 0; i < N2; i++) {
		for (size_t j = 0; j < N1; j++) {
			ResVec[i] += LMat[i][j] * RVec[j];
		}
	}

	// Return value
	return ResVec;
}

// Overloaded * operator for matrix * matrix
template <typename T, size_t N1, size_t N2, size_t N3>
inline array<array<T, N3>, N1> operator*(const array<array<T, N2>, N1> &LMat, const array<array<T, N3>, N2> &RMat) {

	// Declare results matrix
	array<array<T, N3>, N1> ResMat = {0.0};

	// Multiply the left and right matrices
	for (size_t i = 0; i < N1; i++) {
		for (size_t j = 0; j < N3; j++) {
			for (size_t k = 0; k < N2; k++) {
				ResMat[i][j] += LMat[i][k] * RMat[k][j];
			}
		}
	}

	// Return
	return ResMat;
}


// OVERLOADED OPERATORS FOR VECTORS
// Overloaded + operator for vectors
template <typename T>
inline vector<T> operator+(const vector<T> &LVec, const vector<T> &RVec) {

	// Check they are the same size
	if (LVec.size() != RVec.size())
		ERROR("Vectors not right size for adding...exiting");

	// Get rows
	size_t rows = LVec.size();

	// Declare res vector
	vector<T> ResVec;
	ResVec.reserve(rows);

	// Add the left and right vectors
	for (size_t i = 0; i < rows; i++)
		ResVec.push_back(LVec[i] + RVec[i]);

	// Return value
	return ResVec;
}

// Overloaded - operator for vectors
template <typename T>
inline vector<T> operator-(const vector<T> &LVec, const vector<T> &RVec) {

	// Check they are the same size
	if (LVec.size() != RVec.size())
		ERROR("Vectors not right size for subtracting...exiting");

	// Get rows
	size_t rows = LVec.size();

	// Declare res vector
	vector<T> ResVec;
	ResVec.reserve(rows);

	// Add the left and right vectors
	for (size_t i = 0; i < rows; i++)
		ResVec.push_back(LVec[i] - RVec[i]);

	// Return value
	return ResVec;
}

// Overloaded * operator for scalar * vector
template <typename T>
inline vector<T> operator*(const T LScalar, const vector<T> &RVec) {

	// Get rows
	size_t rows = RVec.size();

	// Declare res vector
	vector<T> ResVec;
	ResVec.reserve(rows);

	// Multiply the left and right matrices (right matrix is a column vector)
	for (size_t i = 0; i < rows; i++)
		ResVec.push_back(LScalar * RVec[i]);

	// Return value
	return ResVec;
}

// Overloaded * operator for vector * vector
template <typename T>
inline T operator*(const vector<T> &LVec, const vector<T> &RVec) {

	// Check they are the same size
	if (LVec.size() != RVec.size())
		ERROR("Vectors not right size for multiplying...exiting");

	// Get cols
	size_t cols = LVec.size();

	// Declare res scalar
	T res = 0.0;

	// Multiply vectors
	for (size_t i = 0; i < cols; i++)
		res += LVec[i] * RVec[i];

	// Return value
	return res;
}

// LAPACK interfaces
extern "C" void dgetrf_(int* dim1, int* dim2, double* a, int* lda, int* ipiv, int* info);
extern "C" void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO );

#endif // UTILS_H
