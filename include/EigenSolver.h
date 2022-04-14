/*===========================================================================
This file is part of AC4DC.

    AC4DC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AC4DC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AC4DC.  If not, see <https://www.gnu.org/licenses/>.
===========================================================================*/
#pragma once

#include <vector>
#include <eigen3/Eigen/Dense>

using namespace std;

class EigenSolver
{
public:
	// Class for general purpose linear algebra routines. Interface between Eigen an1d standard library.
	EigenSolver();
	~EigenSolver() {}

	// Solves inhomogenous system of equation of the form Ay = b.
	// The answer is written to b, so that b = y after solution is obtained.
	void SolveSystem(vector<vector<double>> &A, vector<double> &B, int Length);

	// Generalized eigenvalue problem. Solves FC = eSC, where
	// F is the matrix of interest (Fock matrix in Roothan-Hartree-Fock). Matrix has to be self-adjoint.
	// S is the overlap matrix.
	// e is the diagonal eigenvalue matrix.
	// C is the eigenvector matrix (eigenvectors form columns of this matrix).
	void SolveGenEig(vector<std::vector<double>> &F, vector<vector<double>> &S);
	// Same as above, but S is diagonal. It's elements are conmonents of input vector S.
	void SolveGenEig(vector<std::vector<double>> &F, std::vector<double> &S);
	vector<double> EigenVals();
	vector<vector<double>> EigenVecs();
private:
	Eigen::MatrixXd EigenVectors;
	Eigen::VectorXd EigenValues;
};
