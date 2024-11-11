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

#include "Pulse.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include "Constant.h"


void Pulse::save(const std::vector<double>& Tvec, const std::string& fname) {
    std::ofstream f;
    std::cout << "[ Flux ] Saving to file "<<fname<<"..."<<endl;
    f.open(fname);
    f << "# Time (fs) | Intensity (pht/cm2/fs)" <<endl;
    for (auto& t : Tvec) {
        f << t*Constant::fs_per_au << " ";
        double intensity = (*this)(t);
        intensity /= Constant::fs_per_au;
        intensity /= Constant::cm_per_au*Constant::cm_per_au;
        f << intensity <<std::endl;
    }
    f.close();
}

void Pulse::set_pulse(double fluence, double fwhm_param) {
    // The photon flux model
    std::cout<<"[ Flux ] fluence="<<fluence<<", fwhm="<<fwhm_param<<endl;
    this->I0 = fluence/fwhm_param;
    this->fwhm = fwhm_param;

}
void Pulse::set_pulse(double fluence, double fwhm_param, double probe_delay) {
    set_pulse(fluence,fwhm_param);
    this->probe_delay = probe_delay;

}

/**
 * @brief Returns pulse's flux at time t
 * @details If pulse shape is gaussian, t = 0 is defined as the time of peak flux. If it's square, it's defined as the end of the pulse.
 * 
 * @param t time, in same units as the time-integrated flux (fluence)
 * @return double 
 */
double Pulse::operator()(double t) {
    const double norm = sqrt(Constant::Pi/4/log(2));
    switch (shape)
    {
    case PulseShape::gaussian:
        return this->I0/norm*pow(2,-t*t*4/this->fwhm/this->fwhm);
        break;
    case PulseShape::square:
        if (t < -this->fwhm || t >0) {
            return 0;
        } else {
            return I0;
        }
        break;
    case PulseShape::pumpProbeGaussians:{
        double t_offset = t - this->probe_delay;
        return this->I0/norm*pow(2,-t*t*4/this->fwhm/this->fwhm) + this->I0/norm*pow(2,-t_offset*t_offset*4/this->fwhm/this->fwhm);
        break;
    }
    case PulseShape::pumpProbeSquares:{
        double t_offset = t - this->probe_delay;
        double total = 0; 
        if (t < -this->fwhm || t >0) {
        } else {
            total+=I0;
        }
        if (t_offset < -this->fwhm || t_offset >0) {
            return total;
        } else {
            return total + I0;
        }
        break;
    }
    default:
        throw runtime_error("Pulse shape has not been set.");
        break;
    }   
}