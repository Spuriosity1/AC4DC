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
#include "RadialWF.h"

RadialWF::RadialWF(const RadialWF& other) : PairFunction(other), l(other.l), Energy(other.Energy), n(other.n), infinity(other.infinity), turn(other.turn), occup_number(other.occup_number)
{}

int RadialWF::check_nodes()
{
	int Result = 0;
	for (int i = 1; i < turn; i++) {
		if ((F[i] * F[i - 1]) < 0) {
			Result++;
		}
	}

	return Result;
}

void RadialWF::set_L(int X)
{
	l = X;
	occup_number = 4 * l + 2;
}
