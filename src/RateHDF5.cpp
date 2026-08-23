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

#include "RateHDF5.h"
#include <hdf5.h>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>

using CustomDataType::energy_config;
using CustomDataType::ffactor;

namespace {

// --- error handling -------------------------------------------------------
void fail(const std::string& fname, const std::string& what) {
    throw std::runtime_error("[ RateHDF5 ] " + what + " (file: " + fname + ")");
}

bool file_exists(const std::string& path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0;
}

// --- attribute writers ----------------------------------------------------
void write_scalar_attr(hid_t loc, const char* name, hid_t type, const void* value) {
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate2(loc, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, type, value);
    H5Aclose(attr);
    H5Sclose(space);
}

void write_str_attr(hid_t loc, const char* name, const std::string& value) {
    hid_t type = H5Tcopy(H5T_C_S1);
    H5Tset_size(type, value.size() + 1); // include NUL
    H5Tset_strpad(type, H5T_STR_NULLTERM);
    write_scalar_attr(loc, name, type, value.c_str());
    H5Tclose(type);
}

// --- simple 1D dataset writer (native contiguous type) --------------------
void write_1d(hid_t loc, const char* name, hid_t type, size_t n, const void* data) {
    hsize_t dim = n;
    hid_t space = H5Screate_simple(1, &dim, nullptr);
    hid_t dset = H5Dcreate2(loc, name, type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (n > 0) H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(dset);
    H5Sclose(space);
}

// --- compound type for RateData::Rate -------------------------------------
hid_t make_rate_type() {
    hid_t t = H5Tcreate(H5T_COMPOUND, sizeof(RateData::Rate));
    H5Tinsert(t, "val",    HOFFSET(RateData::Rate, val),    H5T_NATIVE_DOUBLE);
    H5Tinsert(t, "from",   HOFFSET(RateData::Rate, from),   H5T_NATIVE_LONG);
    H5Tinsert(t, "to",     HOFFSET(RateData::Rate, to),     H5T_NATIVE_LONG);
    H5Tinsert(t, "energy", HOFFSET(RateData::Rate, energy), H5T_NATIVE_DOUBLE);
    return t;
}

// --- compound type for CustomDataType::energy_config ----------------------
hid_t make_energy_config_type() {
    hid_t t = H5Tcreate(H5T_COMPOUND, sizeof(energy_config));
    H5Tinsert(t, "index",          HOFFSET(energy_config, index),          H5T_NATIVE_INT);
    H5Tinsert(t, "valence_energy", HOFFSET(energy_config, valence_energy), H5T_NATIVE_DOUBLE);
    H5Tinsert(t, "receiver_index", HOFFSET(energy_config, receiver_index), H5T_NATIVE_INT);
    H5Tinsert(t, "donator_index",  HOFFSET(energy_config, donator_index),  H5T_NATIVE_INT);
    return t;
}

void write_rate_dataset(hid_t loc, const char* name, const std::vector<RateData::Rate>& v) {
    hid_t type = make_rate_type();
    write_1d(loc, name, type, v.size(), v.data());
    H5Tclose(type);
}

// --- variable-length string dataset (index_names) -------------------------
void write_vlen_strings(hid_t loc, const char* name, const std::vector<std::string>& strs) {
    hid_t type = H5Tcopy(H5T_C_S1);
    H5Tset_size(type, H5T_VARIABLE);
    std::vector<const char*> ptrs(strs.size());
    for (size_t i = 0; i < strs.size(); ++i) ptrs[i] = strs[i].c_str();
    write_1d(loc, name, type, ptrs.size(), ptrs.data());
    H5Tclose(type);
}

// --- 2D form-factor dataset -----------------------------------------------
void write_form_factors(hid_t loc, const std::vector<ffactor>& ff) {
    // index dataset
    std::vector<int> idx(ff.size());
    for (size_t i = 0; i < ff.size(); ++i) idx[i] = ff[i].index;
    write_1d(loc, "form_factor_index", H5T_NATIVE_INT, idx.size(), idx.data());

    // 2D form_factors dataset (rows = configs, cols = q-grid); {0,0} when ff is empty.
    const size_t rows = ff.size();
    const size_t q = ff.empty() ? 0 : ff[0].val.size();
    std::vector<double> flat(rows * q);
    for (size_t i = 0; i < rows; ++i) {
        if (ff[i].val.size() != q)
            throw std::runtime_error("[ RateHDF5 ] non-uniform form-factor length; "
                                     "2D storage assumption violated");
        for (size_t j = 0; j < q; ++j) flat[i * q + j] = ff[i].val[j];
    }

    hsize_t dims[2] = {rows, q};
    hid_t space = H5Screate_simple(2, dims, nullptr);
    hid_t dset = H5Dcreate2(loc, "form_factors", H5T_NATIVE_DOUBLE, space,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (rows * q > 0)
        H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat.data());
    H5Dclose(dset);
    H5Sclose(space);
}

// --- flattened EII group --------------------------------------------------
void write_eii(hid_t file, const std::vector<RateData::EIIdata>& eii) {
    hid_t grp = H5Gcreate2(file, "EII", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    const size_t M = eii.size();
    std::vector<int> init(M), counts(M);
    std::vector<int> fin_flat, occ_flat;
    std::vector<float> ionB_flat, kin_flat;
    size_t total = 0;
    for (const auto& e : eii) total += e.fin.size();
    fin_flat.reserve(total);
    occ_flat.reserve(total);
    ionB_flat.reserve(total);
    kin_flat.reserve(total);
    for (size_t i = 0; i < M; ++i) {
        const RateData::EIIdata& e = eii[i];
        init[i] = e.init;
        counts[i] = static_cast<int>(e.fin.size());
        fin_flat.insert(fin_flat.end(), e.fin.begin(), e.fin.end());
        occ_flat.insert(occ_flat.end(), e.occ.begin(), e.occ.end());
        ionB_flat.insert(ionB_flat.end(), e.ionB.begin(), e.ionB.end());
        kin_flat.insert(kin_flat.end(), e.kin.begin(), e.kin.end());
    }
    write_1d(grp, "init",   H5T_NATIVE_INT,   init.size(),      init.data());
    write_1d(grp, "counts", H5T_NATIVE_INT,   counts.size(),    counts.data());
    write_1d(grp, "fin",    H5T_NATIVE_INT,   fin_flat.size(),  fin_flat.data());
    write_1d(grp, "occ",    H5T_NATIVE_INT,   occ_flat.size(),  occ_flat.data());
    write_1d(grp, "ionB",   H5T_NATIVE_FLOAT, ionB_flat.size(), ionB_flat.data());
    write_1d(grp, "kin",    H5T_NATIVE_FLOAT, kin_flat.size(),  kin_flat.data());
    H5Gclose(grp);
}

// ==== readers =============================================================
double read_double_attr(hid_t loc, const char* name) {
    hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
    double v = 0;
    H5Aread(attr, H5T_NATIVE_DOUBLE, &v);
    H5Aclose(attr);
    return v;
}

int read_int_attr(hid_t loc, const char* name) {
    hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
    int v = 0;
    H5Aread(attr, H5T_NATIVE_INT, &v);
    H5Aclose(attr);
    return v;
}

// number of elements in a 1D dataset
size_t read_len(hid_t loc, const char* name) {
    hid_t dset = H5Dopen2(loc, name, H5P_DEFAULT);
    hid_t space = H5Dget_space(dset);
    hsize_t dim = 0;
    H5Sget_simple_extent_dims(space, &dim, nullptr);
    H5Sclose(space);
    H5Dclose(dset);
    return dim;
}

template <typename T>
std::vector<T> read_1d(hid_t loc, const char* name, hid_t type) {
    size_t n = read_len(loc, name);
    std::vector<T> out(n);
    if (n > 0) {
        hid_t dset = H5Dopen2(loc, name, H5P_DEFAULT);
        H5Dread(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
        H5Dclose(dset);
    }
    return out;
}

std::vector<RateData::Rate> read_rate_dataset(hid_t loc, const char* name) {
    size_t n = read_len(loc, name);
    std::vector<RateData::Rate> out(n);
    if (n > 0) {
        hid_t type = make_rate_type();
        hid_t dset = H5Dopen2(loc, name, H5P_DEFAULT);
        H5Dread(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
        H5Dclose(dset);
        H5Tclose(type);
    }
    return out;
}

std::vector<std::string> read_vlen_strings(hid_t loc, const char* name) {
    size_t n = read_len(loc, name);
    std::vector<std::string> out;
    if (n == 0) return out;
    hid_t type = H5Tcopy(H5T_C_S1);
    H5Tset_size(type, H5T_VARIABLE);
    std::vector<char*> ptrs(n, nullptr);
    hid_t dset = H5Dopen2(loc, name, H5P_DEFAULT);
    H5Dread(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptrs.data());
    for (size_t i = 0; i < n; ++i) out.emplace_back(ptrs[i] ? ptrs[i] : "");
    // reclaim vlen buffers allocated by the HDF5 library
    hid_t space = H5Dget_space(dset);
    H5Dvlen_reclaim(type, space, H5P_DEFAULT, ptrs.data());
    H5Sclose(space);
    H5Dclose(dset);
    H5Tclose(type);
    return out;
}

std::vector<ffactor> read_form_factors(hid_t loc) {
    std::vector<int> idx = read_1d<int>(loc, "form_factor_index", H5T_NATIVE_INT);
    std::vector<ffactor> out;
    hid_t dset = H5Dopen2(loc, "form_factors", H5P_DEFAULT);
    hid_t space = H5Dget_space(dset);
    hsize_t dims[2] = {0, 0};
    H5Sget_simple_extent_dims(space, dims, nullptr);
    H5Sclose(space);
    const size_t rows = dims[0], q = dims[1];
    std::vector<double> flat(rows * q);
    if (rows * q > 0)
        H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat.data());
    H5Dclose(dset);
    out.reserve(rows);
    for (size_t i = 0; i < rows; ++i) {
        ffactor f;
        f.index = (i < idx.size()) ? idx[i] : static_cast<int>(i);
        f.val.assign(flat.begin() + i * q, flat.begin() + (i + 1) * q);
        out.push_back(std::move(f));
    }
    return out;
}

std::vector<RateData::EIIdata> read_eii(hid_t file) {
    hid_t grp = H5Gopen2(file, "EII", H5P_DEFAULT);
    std::vector<int> init   = read_1d<int>(grp, "init",   H5T_NATIVE_INT);
    std::vector<int> counts = read_1d<int>(grp, "counts", H5T_NATIVE_INT);
    std::vector<int> fin    = read_1d<int>(grp, "fin",    H5T_NATIVE_INT);
    std::vector<int> occ    = read_1d<int>(grp, "occ",    H5T_NATIVE_INT);
    std::vector<float> ionB = read_1d<float>(grp, "ionB", H5T_NATIVE_FLOAT);
    std::vector<float> kin  = read_1d<float>(grp, "kin",  H5T_NATIVE_FLOAT);
    H5Gclose(grp);

    std::vector<RateData::EIIdata> out(init.size());
    size_t off = 0;
    for (size_t i = 0; i < init.size(); ++i) {
        const size_t c = counts[i];
        out[i].init = init[i];
        out[i].fin.assign(fin.begin() + off, fin.begin() + off + c);
        out[i].occ.assign(occ.begin() + off, occ.begin() + off + c);
        out[i].ionB.assign(ionB.begin() + off, ionB.begin() + off + c);
        out[i].kin.assign(kin.begin() + off, kin.begin() + off + c);
        off += c;
    }
    return out;
}

} // namespace

namespace RateHDF5 {

std::string filename(const std::string& element, double omega_eV, const std::string& dir) {
    std::ostringstream ss;
    ss << dir << "/" << element << "_" << (long)std::llround(omega_eV) << "eV.h5";
    return ss.str();
}

void write(const std::string& fname, const RateData::Atom& atom,
           const std::vector<ffactor>& ff, int Z, double omega_eV) {
    hid_t file = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) fail(fname, "could not create file");

    // --- root attributes ---
    write_str_attr(file, "element", atom.name);
    write_scalar_attr(file, "Z", H5T_NATIVE_INT, &Z);
    write_scalar_attr(file, "photon_energy_eV", H5T_NATIVE_DOUBLE, &omega_eV);
    unsigned int num_conf = atom.num_conf;
    write_scalar_attr(file, "num_conf", H5T_NATIVE_UINT, &num_conf);
    write_str_attr(file, "code_version", "AC4DC 0.2.0");
    write_str_attr(file, "units", "atomic units");

    // --- rate datasets ---
    write_rate_dataset(file, "Photo", atom.Photo);
    write_rate_dataset(file, "Fluor", atom.Fluor);
    write_rate_dataset(file, "Auger", atom.Auger);

    // --- configuration metadata ---
    write_vlen_strings(file, "index_names", atom.index_names);
    hid_t ec_type = make_energy_config_type();
    write_1d(file, "EnergyConfig", ec_type, atom.EnergyConfig.size(), atom.EnergyConfig.data());
    H5Tclose(ec_type);

    // --- EII (flattened group) & form factors ---
    write_eii(file, atom.EIIparams);
    write_form_factors(file, ff);

    H5Fclose(file);
}

bool read(const std::string& fname, RateData::Atom& atom,
          std::vector<ffactor>& ff, double expected_omega_eV) {
    if (!file_exists(fname)) return false;

    hid_t file = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) fail(fname, "exists but could not be opened");

    double stored_omega = read_double_attr(file, "photon_energy_eV");
    if (std::fabs(stored_omega - expected_omega_eV) > 1.0) {
        std::cerr << "\033[1;35m[ RateHDF5 ] Warning:\033[0m " << fname
                  << " was computed for photon energy " << stored_omega
                  << " eV but " << expected_omega_eV << " eV was requested." << std::endl;
    }

    atom.num_conf = static_cast<unsigned int>(read_int_attr(file, "num_conf"));
    atom.Photo = read_rate_dataset(file, "Photo");
    atom.Fluor = read_rate_dataset(file, "Fluor");
    atom.Auger = read_rate_dataset(file, "Auger");
    atom.index_names = read_vlen_strings(file, "index_names");

    hid_t ec_type = make_energy_config_type();
    atom.EnergyConfig = read_1d<energy_config>(file, "EnergyConfig", ec_type);
    H5Tclose(ec_type);

    atom.EIIparams = read_eii(file);
    ff = read_form_factors(file);

    H5Fclose(file);
    return true;
}

} // namespace RateHDF5
