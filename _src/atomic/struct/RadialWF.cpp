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
#include <numeric>
#include <algorithm>
#include <iostream>

RadialWF::RadialWF(const RadialWF& other) : PairFunction(other), l(other.l), Energy(other.Energy), n(other.n), infinity(other.infinity), turn(other.turn), occup_number(other.occup_number),shell_flag(other.shell_flag) // Necessary due to use of push_back. Cursed. - S.P. 
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

void RadialWF::set_L(int X, bool reset_occupancy)
{
    //std::cout<<"SETTING L"<<l<<std::endl;
	l = X;
    if (reset_occupancy){
	    occup_number = 4 * l + 2;
    }
}

void RadialWF::flag_shell(bool force_flag)
{
    if (l == -10 || force_flag) shell_flag = true;   

}

void RadialWF::shell_flag_check()
{
    if (shell_flag)set_L(-10,false);
    
}

/**
 * @brief Recover{s, p, d, f} orbital counts, assuming that atom in ground state. TODO I made this really complicated because I was scared of adding variables -S.P.
 * @return a vector with each element corresponding to the number of electrons with angular momentum equal to the index. e.g. occupancies for 2p 2 -> {0,2,0,0},  3N 5 -> {2,3,0,0}
 */
std::vector<int> RadialWF::get_subshell_occupancies()
{
    shell_flag_check();
    // Orbital case
    if (l != -10){  //TODO replace with shell_flag - S.P.
        std::vector<int> shell_orbital_occupancy = {0,0,0,0};
        if (l > 3){
            std::cerr<<"ERROR: invalid orbital angular momentum encountered! Ensure that all [atom].inp subshells are s, p, d, f, or N.";
        }
        shell_orbital_occupancy[l] = occup_number;
        // std::cout<<"l = "<<l<<" occup_number = "<<occup_number<<std::endl;
        // std::cout<<"Returning {";
        // for(auto i: shell_orbital_occupancy) std::cout << i << ' ';
        // std::cout<<"}"<<std::endl;
        return shell_orbital_occupancy;
    }
    // Shell (averaged orbitals) case
    std::vector<int> shell_orbital_occupancy = {2,6,10,14};
    for (int i = shell_orbital_occupancy.size() - 1; i >= 0; i--)
    {
        // If higher subshell has electrons, lower subshells are filled.
        if(i != shell_orbital_occupancy.size() - 1 && shell_orbital_occupancy[i+1] > 0){
            continue;
        }
        // Higher subshells empty
        else{
            // Determine how many electrons are free to add to this subshell after accounting for lower orbitals.
            shell_orbital_occupancy[i] = 0;
            for(int j = i - 1; j >= 0;j--){
                shell_orbital_occupancy[i] -= shell_orbital_occupancy[j];
            }					
            shell_orbital_occupancy[i] = std::max(0, occup_number + shell_orbital_occupancy[i]);
        }
        
    }
    if (shell_orbital_occupancy.back() > 14)
    {
        std::cerr<<"WARNING: g subshell not supported (> 32 electrons in one shell detected in [atom].inp file).";
    }
    // std::cout<<"l = "<<l<<" occup_number = "<<occup_number<<std::endl;
    // std::cout<<"Returning {";
    // for(auto i: shell_orbital_occupancy) std::cout << i << ' ';
    // std::cout<<"}"<<std::endl;

    return shell_orbital_occupancy;
}
