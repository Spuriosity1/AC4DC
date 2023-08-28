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
#include <vector>
#include <algorithm>
#include "Wigner3j.h"
#include <cmath>


double LogFactorialFraction(double Num, double Denom);

namespace Constant
{
	double Wigner3j(double j1, double j2, double j3, double m1, double m2, double m3)
	{
		double Result = 0.;
		double fracpart, intpart, J = j1 + j2 + j3;

		//some particular cases can be calculated in simpler and faster way.
		if (m1 + m2 + m3 != 0) return 0.;
		if (j1 < 0. || j2 < 0. || j3 < 0.) return 0.;
		fracpart = std::modf(J, &intpart);// J should be integer
		if (fracpart != 0) return 0;
		else if (m1 == 0 && m2 == 0)
		{
			// /j1 j2 j3\ - Edmonds, "Angular momentum in qantum mechanics", p 50
			// \0  0  0 /
			fracpart = std::modf(0.5*J, &intpart);
			if (fracpart != 0) return 0;
			else
			{
				if (j1 > j2 + j3 || j2 > j1 + j3 || j3 > j1 + j2 ||
					j1 < fabs(j2 - j3) || j2 < fabs(j1 - j3) || j3 < fabs(j2 - j1)) return 0;
				std::vector<double> factor(3);
				factor[0] = J / 2. - j1;
				factor[1] = J / 2. - j2;
				factor[2] = J / 2. - j3;
				std::sort(factor.begin(), factor.end());//smallest first, largest last

				for (unsigned int i = 0; i < factor.size(); i++)
				{
					fracpart = std::modf(factor[i], &intpart);// J should be integer
					if (fracpart != 0) return 0;
				}

				Result += LogFactorialFraction(2 * factor[0], factor[2]);//taking the largest factor under the square root in denominator (factor[2]!)*(factor[2]!)
				Result += LogFactorialFraction(2 * factor[1], factor[2]);
				Result += LogFactorialFraction(2 * factor[2], J + 1);
				Result *= 0.5;// end of square root

				Result += LogFactorialFraction(J / 2, factor[1]);
				Result += LogFactorialFraction(1, factor[0]);//smallest factor goes unpaired with numerator

				fracpart = std::modf(0.25*J, &intpart);//fracpart = 0 if J/2 is even
				if (fracpart == 0) { Result = exp(Result); }
				else { Result = -exp(Result); }

				return Result;
			}
		}
		else if ((m1 == 0 || m2 == 0 || m3 == 0) && (m1 == 0.5 || m2 == 0.5 || m3 == 0.5))
		{
			// /ja   jb   J\ - Brink and Satcher, "Angular momentum", p 138
			// \0.5 -0.5  0 /
			double ja, jb;
			int columns_permutation = 1;//account for sign flip
			if (j1 > j2 + j3 || j2 > j1 + j3 || j3 > j1 + j2 ||
				j1 < fabs(j2 - j3) || j2 < fabs(j1 - j3) || j3 < fabs(j2 - j1)) return 0;

			if (m1 == 0)
			{
				J = j1;
				if (m2 == 0.5)
				{
					ja = j2;
					jb = j3;
				}
				else
				{
					ja = j3;
					jb = j2;
					columns_permutation = -1;
				}
			}
			else if (m1 == 0.5)
			{
				ja = j1;
				if (m2 == 0)
				{
					J = j2;
					jb = j3;
					columns_permutation = -1;
				}
				else
				{
					J = j3;
					jb = j2;
				}
			}
			else
			{
				jb = j1;
				if (m2 == 0)
				{
					J = j2;
					ja = j3;
				}
				else
				{
					J = j3;
					ja = j2;
					columns_permutation = -1;
				}
			}

			double K;
			fracpart = std::modf(0.5*(ja + jb + J), &intpart);
			if (fracpart == 0) { K = J; }
			else  { K = J + 1; }

			std::vector<double> factor_qsrt(3);//sort factorials under square root
			factor_qsrt[0] = ja + jb - J;
			factor_qsrt[1] = ja + J - jb;
			factor_qsrt[2] = jb + J - ja;
			std::sort(factor_qsrt.begin(), factor_qsrt.end());

			std::vector<double> factor(3);//sort factorials under square root
			factor[0] = (ja + jb - K) / 2;
			factor[1] = (ja + K - jb - 1) / 2;
			factor[2] = (jb + J - ja) / 2;
			std::sort(factor.begin(), factor.end());

			Result += LogFactorialFraction(factor_qsrt[0], factor[2]);
			Result += LogFactorialFraction(factor_qsrt[1], factor[2]);
			Result += LogFactorialFraction(factor_qsrt[2], (J + ja + jb + 1));
			Result *= 0.5;

			Result += LogFactorialFraction(1, factor[0]);
			Result += LogFactorialFraction(0.5*(K + ja + jb), factor[1]);

			fracpart = std::modf(0.25*(K + ja + jb) - 0.5, &intpart);
			if (fracpart != 0) columns_permutation *= -1;
			Result = columns_permutation*exp(Result)*2. / sqrt((2 * ja + 1)*(2 * jb + 1));

			return Result;
		}
		else
		{
			std::vector<double> Regge(9);
			std::vector<double> ReggeInt(9);
			Regge[0] = j2 + j3 - j1;//R(11)
			Regge[1] = j1 - j2 + j3;//R(12)
			Regge[2] = j1 + j2 - j3;//R(13)
			Regge[3] = j1 + m1;//R(21)
			Regge[4] = j2 + m2;//R(22)
			Regge[5] = j3 + m3;//R(23)
			Regge[6] = j1 - m1;//R(31)
			Regge[7] = j2 - m2;//R(32)
			Regge[8] = j3 - m3;//R(33)

			// if any Regge symbol is <0 or not integer, 3j symbol = 0
			for (unsigned int i = 0; i < Regge.size(); i++)
			{
				if (Regge[i] < 0.) return 0.;
				fracpart = std::modf(Regge[i], &intpart);
				if (Regge[i] != intpart) return 0.;
				ReggeInt[i] = (int)intpart;
			}

			//using 8.3.29 of Varshalovich. Outer part is everything under the square root==================NEEDS TO BE FINISHED=======================

			return 0.;
		}

	}
}

//this function evaluates log(Num!/Denom!). It is used to calculate factorial summations in Wigner3j funstion. Allows to avoid getting huge factorial products
double LogFactorialFraction(double Num, double Denom)
{
	double Result = 0.;
	if (Num > Denom)
	{
		for (int i = Denom + 1; i <= Num; i++)
		{
			Result += log((double)i);
		}
	}
	else if (Num < Denom)
	{
		for (int i = Num + 1; i <= Denom; i++)
		{
			Result -= log((double)i);
		}
	}
	return Result;
}


