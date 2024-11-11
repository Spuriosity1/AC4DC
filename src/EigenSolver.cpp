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
#include "EigenSolver.h"
#include <vector>
#include "stdafx.h"

EigenSolver::EigenSolver()
{
	//this class would test sever different pakages for matrix diagonalization. It uses just "Eigen" for now. Reminder: "Armadillo" in future
}

void EigenSolver::SolveSystem(std::vector<std::vector<double>> &A, std::vector<double> &B, int Length)
{
	Eigen::MatrixXd M(Length, Length);
	Eigen::VectorXd Y(Length);
	Eigen::VectorXd X(Length);

	for (int i = 0; i < Length; i++)
	{
		Y(i) = B[i];
		if (std::isnan(B[i])){throw std::runtime_error("SolveSystem encountered nan!");}
		for (int j = 0; j < Length; j++)
		{
			M(i, j) = A[i][j];
			if (std::isnan(A[i][j])){throw std::runtime_error("SolveSystem encountered nan!");}
		}
	}

	X = M.fullPivLu().solve(Y);

	for (int i = 0; i < Length; i++)
	{
		B[i] = X(i);
		if (std::isnan(B[i])){
			typedef Eigen::MatrixXd MatrixXd;
			cout << "Here is the matrix M:" << endl << M << endl;
			Eigen::FullPivLU<MatrixXd> lu(M);			

			cout << "And here is the RHS, the vector Y:" << endl;
			Eigen::IOFormat fmt(4, 0, ", ", "\n", "", "");
			cout << Y.transpose().format(fmt) << endl;
			
			throw std::runtime_error("SolveSystem encountered nan!");
		}
		if (std::isinf(Y[0])){
			throw std::runtime_error("SolveSystem encountered inf!");
		}
	}
}

void EigenSolver::SolveGenEig(std::vector<std::vector<double>> &F, std::vector<std::vector<double>> &S)
{
	int Length = F.size();
	Eigen::MatrixXd M(Length, Length);
	Eigen::MatrixXd Q(Length, Length);

	for (int i = 0; i < Length; i++) {
		for (int j = 0; j < Length; j++) {
			M(i, j) = F[i][j];
			Q(i, j) = S[i][j];
		}
	}

	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> Magic(M, Q);

	EigenValues = Magic.eigenvalues();
	EigenVectors = Magic.eigenvectors();
}

vector<double> EigenSolver::EigenVals()
{
	vector<double> Result(EigenValues.size(), 0);

	for (size_t i = 0; i < Result.size(); i++) {
		Result[i] = EigenValues(i);
	}

	return Result;
}

vector<vector<double>> EigenSolver::EigenVecs()
{
	vector<vector<double>> Result(EigenValues.size(), vector<double>(EigenValues.size(), 0));
	for (size_t i = 0; i < Result.size(); i++) {
		for (size_t j = 0; j < Result.size(); j++) {
			Result[i][j] = EigenVectors(j, i);
		}
	}

	return Result;
}
