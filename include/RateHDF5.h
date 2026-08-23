/**
 * @file RateHDF5.h
 * @brief Serialization of per-atom Hartree-Fock rate data to/from HDF5 files.
 * @details This is the file-based seam between the two stages of the suite:
 *   `atomic_rate_data` (Stage 1, Hartree-Fock) writes a RateData::Atom + form
 *   factors for a given (element, photon energy), and `ac4dc` (Stage 2, plasma
 *   dynamics) reads them back. Uses the portable HDF5 C API.
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

#include <string>
#include <vector>
#include "Constant.h"

namespace RateHDF5 {
    /// Default directory in which atomic rate files live. Matches where the
    /// atomic_rate_data binary writes them when pointed at input/atoms/<el>.inp
    /// with no -o override, so the two stages line up out of the box.
    const std::string default_dir = "input/atoms";

    /// Returns the canonical file path for a given (element, photon energy),
    /// e.g. filename("C", 9000) -> "input/atoms/C_9000eV.h5" (using default_dir).
    /// The photon energy is rounded to the nearest eV for the (searchable) name;
    /// the exact value is stored inside the file as an attribute.
    std::string filename(const std::string& element, double omega_eV,
                         const std::string& dir = default_dir);

    /// Writes the full atomic rate data (everything Stage 2 consumes, plus form
    /// factors for the scattering code) to an HDF5 file, overwriting any existing
    /// file. atom.nAtoms / atom.name are molecule-specific and are NOT stored.
    /// Throws std::runtime_error on any HDF5 failure.
    void write(const std::string& fname, const RateData::Atom& atom,
               const std::vector<CustomDataType::ffactor>& ff,
               int Z, double omega_eV);

    /// Reads atomic rate data written by write(). Populates every RateData::Atom
    /// field except nAtoms and name (which the caller overlays from the .mol).
    /// Returns false if the file does not exist. Prints a warning (non-fatal) if
    /// the stored photon energy differs from expected_omega_eV by more than ~1 eV.
    /// Throws std::runtime_error on a corrupt/unreadable file.
    bool read(const std::string& fname, RateData::Atom& atom,
              std::vector<CustomDataType::ffactor>& ff, double expected_omega_eV);
}
