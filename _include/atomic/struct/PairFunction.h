/**
 * @file PairFunction.h
 * @brief
 */
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

class PairFunction
{
public:
	PairFunction(int size = 0)
	{
		clear();
		if (size > 0) {
			resize(size);
		}
	}

	std::vector<double> F, G;  // Large and small components of radial wavefunction and it's derivatives

	void clear()
	{
		F.clear();
		G.clear();
	}

	void resize(int size)
	{
		F.resize(size);
		G.resize(size);
	}

	int size() { return static_cast<int>(F.size()); }

	void scale(double norm)
	{
		for (int i = 0; i < size(); i++) {
			F[i] *= norm;
			G[i] *= norm;
		}
	}

	const PairFunction& operator=(const PairFunction& Psi)
	{
		F = Psi.F;
		G = Psi.G;

		return *this;
	}

	~PairFunction() {}
};
