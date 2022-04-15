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
#include <iostream>
#include <cassert>
#ifndef LOSS_GEOMETRY_CXX_H
#define LOSS_GEOMETRY_CXX_H

// Overkill, to be sure

struct LossGeometry {
    const static int sphere = 1;
    const static int cylinder = 2;
    const static int plane = 3;
    const static int none = 0;
    int mode = 1;
    double L0 = 1;
    double factor() const {
        assert(mode >=0 && mode <=2);
        return 1.*constants[mode]/L0;
    }
    private:
    const double constants[4] = {0, 2, 3, 4};
};

namespace{
    [[maybe_unused]] std::ostream& operator<<(std::ostream& os, const LossGeometry& g) {
        switch (g.mode)
        {
        case LossGeometry::none:
            os << "closed (lossless)";
            break;
        case LossGeometry::sphere:
            os << "Sphere";
            break;
        case LossGeometry::cylinder:
            os << "Thin Cylinder";
            break;
        case LossGeometry::plane:
            os << "Thin Plane";
            break;
        default:
            os << "Unknown shape";
            break;
        }
        return os;
    }

    [[maybe_unused]] std::istream& operator>>(std::istream& is, LossGeometry& lg) {
        std::string tmp;
        is >> tmp;
        if (tmp.length() == 0) {
            std::cerr<<"No geometry specifier provided, defaulting to spherical..."<<std::endl;
            lg.mode = LossGeometry::sphere;
            return is;
        }
        switch ((char) tmp[0])
        {
        case 'n':
            lg.mode = LossGeometry::none;
            break;
        case 's':
            lg.mode = LossGeometry::sphere;
            break;
        case 'c':
            lg.mode = LossGeometry::cylinder;
            break;
        case 'p':
            lg.mode = LossGeometry::plane;
            break;
        default:
            std::cerr<<"Unrecognised grid type \""<<tmp<<"\""<<std::endl;
            lg.mode = LossGeometry::sphere;
            break;
        }
        return is;
    }
}

#endif