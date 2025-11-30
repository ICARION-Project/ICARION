// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "fieldArrayLoader.h"
#include <iostream>
#include <stdexcept>

/** Helper: linear index from i,j,k - use same convention as Grid3D */
static inline size_t idx_3d(size_t i, size_t j, size_t k, size_t nx, size_t ny) {
    return i + nx * (j + ny * k);  // Same as Grid3D convention
}

FieldArray load_field_array(const std::string& path) {
    FieldArray fld;

    // Open file
    hid_t file_id = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        std::cerr << "ERROR: cannot open HDF5 file " << path << "\n";
        return fld;
    }

    // Helper lambda for loading a 1D array
    auto load_1d = [&](const char* dset_name, std::vector<double>& vec_out) {
        hid_t dset = H5Dopen2(file_id, dset_name, H5P_DEFAULT);
        if (dset < 0) {
            throw std::runtime_error(std::string("Dataset ") + dset_name + " not found");
        }
        hid_t space = H5Dget_space(dset);
        int ndims = H5Sget_simple_extent_ndims(space);
        if (ndims != 1) {
            throw std::runtime_error(std::string("Dataset ") + dset_name + " is not 1D");
        }
        hsize_t dims[1];
        H5Sget_simple_extent_dims(space, dims, nullptr);
        vec_out.resize(dims[0]);
        herr_t status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec_out.data());
        H5Sclose(space);
        H5Dclose(dset);
        if (status < 0) {
            throw std::runtime_error(std::string("Failed to read dataset ") + dset_name);
        }
    };

    // Helper lambda for loading a 3D array (into std::vector<double>)
    auto load_3d = [&](const char* dset_name, std::vector<double>& vec_out,
                       size_t &out_nx, size_t &out_ny, size_t &out_nz) {
        hid_t dset = H5Dopen2(file_id, dset_name, H5P_DEFAULT);
        // Check if dataset exists
        if (dset < 0) {
            throw std::runtime_error(std::string("Dataset ") + dset_name + " not found");
        }
        hid_t space = H5Dget_space(dset);
        int ndims = H5Sget_simple_extent_ndims(space);
        if (ndims != 3) {
            throw std::runtime_error(std::string("Dataset ") + dset_name + " is not 3D");
        }
        hsize_t dims[3];
        H5Sget_simple_extent_dims(space, dims, nullptr);
        out_nx = dims[0];
        out_ny = dims[1];
        out_nz = dims[2];
        size_t total = out_nx * out_ny * out_nz;
        vec_out.resize(total);
        herr_t status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec_out.data());
        H5Sclose(space);
        H5Dclose(dset);
        if (status < 0) {
            throw std::runtime_error(std::string("Failed to read dataset ") + dset_name);
        }
    };

    try {
        load_1d("x", fld.xs);
        load_1d("y", fld.ys);
        load_1d("z", fld.zs);
        load_3d("Ex", fld.Ex, fld.nx, fld.ny, fld.nz);
        size_t nx2, ny2, nz2;
        load_3d("Ey", fld.Ey, nx2, ny2, nz2);
        load_3d("Ez", fld.Ez, nx2, ny2, nz2);
        // Expectation: nx2 == fld.nx, etc.
        
        // Load phi if present (optional)
        size_t nx3, ny3, nz3;
        try {
            load_3d("phi", fld.phi, nx3, ny3, nz3);
        } catch (...) {
            // phi is optional, ignore if not found
        }

    } catch (const std::exception& ex) {
        std::cerr << "ERROR loading field array: " << ex.what() << "\n";
        H5Fclose(file_id);
        // Return an invalid field
        FieldArray empty;
        return empty;
    }

    H5Fclose(file_id);
    return fld;
}

Vec3 interpolate_field(const FieldArray& fld, const Vec3& pos) {
    if (!fld.is_valid()) {
        return Vec3{0.0, 0.0, 0.0};
    }

    // Locate the surrounding grid points: i0 ≤ pos.x ≤ i1, etc.
    // FFor simple implementation: linear search + clipping — you can optimize (bin search) later.

    auto locate = [](const std::vector<double>& arr, double v, size_t &i0, size_t &i1, double &alpha) {
        // arr monoton increasing
        size_t n = arr.size();
        if (v <= arr.front()) {
            i0 = 0; i1 = 1; alpha = 0.0;
        } else if (v >= arr.back()) {
            i0 = n-2; i1 = n-1; alpha = 1.0;
        } else {
            // linear search (or bisection optimal)
            for (size_t i = 0; i < n-1; ++i) {
                if (v >= arr[i] && v <= arr[i+1]) {
                    i0 = i;
                    i1 = i+1;
                    double dv = arr[i+1] - arr[i];
                    alpha = (dv > 0.0) ? (v - arr[i]) / dv : 0.0;
                    return;
                }
            }
            // fallback
            i0 = n-2; i1 = n-1; alpha = 1.0;
        }
    };

    size_t ix0, ix1, iy0, iy1, iz0, iz1;
    double ax, ay, az;
    locate(fld.xs, pos.x, ix0, ix1, ax);
    locate(fld.ys, pos.y, iy0, iy1, ay);
    locate(fld.zs, pos.z, iz0, iz1, az);
    

    // interpolation in 3D: Sum over 8 corners
    Vec3 result{0.0, 0.0, 0.0};
    for (size_t dx = 0; dx <= 1; ++dx) {
        double wx = (dx == 0 ? (1.0 - ax) : ax);
        size_t ix = (dx == 0 ? ix0 : ix1);
        for (size_t dy = 0; dy <= 1; ++dy) {
            double wy = (dy == 0 ? (1.0 - ay) : ay);
            size_t iy = (dy == 0 ? iy0 : iy1);
            for (size_t dz = 0; dz <= 1; ++dz) {
                double wz = (dz == 0 ? (1.0 - az) : az);
                size_t iz = (dz == 0 ? iz0 : iz1);
                double w = wx * wy * wz;
                size_t idx = idx_3d(ix, iy, iz, fld.nx, fld.ny);

                result.x += w * fld.Ex[idx];
                result.y += w * fld.Ey[idx];
                result.z += w * fld.Ez[idx];
            }
        }
    }
    return result;
}

/** 
 * @brief Convert a Grid3D to a FieldArray representation
 * @details Important to make Space-Charge solver compatible
 * with field interpolation routines.
 * @param grid Input Grid3D
 * @return FieldArray Converted field array
 */
FieldArray grid_to_fieldarray(const Grid3D& grid)
{
    FieldArray fld;
    fld.nx = grid.Nx;
    fld.ny = grid.Ny;
    fld.nz = grid.Nz;

    // coordinate axes (assuming uniform grid)
    fld.xs.resize(grid.Nx);
    fld.ys.resize(grid.Ny);
    fld.zs.resize(grid.Nz);
    for (int i = 0; i < grid.Nx; ++i) fld.xs[i] = grid.origin_m.x + i * grid.dx;
    for (int j = 0; j < grid.Ny; ++j) fld.ys[j] = grid.origin_m.y + j * grid.dy;
    for (int k = 0; k < grid.Nz; ++k) fld.zs[k] = grid.origin_m.z + k * grid.dz;

    // flatten E-field components with correct indexing
    size_t N = grid.size();
    fld.Ex.resize(N);
    fld.Ey.resize(N);
    fld.Ez.resize(N);
    fld.phi.resize(N);

    // Convert from Grid3D indexing to FieldArray indexing
    for (int i = 0; i < grid.Nx; ++i) {
        for (int j = 0; j < grid.Ny; ++j) {
            for (int k = 0; k < grid.Nz; ++k) {
                size_t grid_idx = grid.index(i, j, k);          // Grid3D indexing
                size_t field_idx = idx_3d(i, j, k, grid.Nx, grid.Ny); // Same indexing as Grid3D
                
                fld.Ex[field_idx] = grid.E[grid_idx].x;
                fld.Ey[field_idx] = grid.E[grid_idx].y;
                fld.Ez[field_idx] = grid.E[grid_idx].z;
                fld.phi[field_idx] = grid.phi[grid_idx];
                

            }
        }
    }

    return fld;
}