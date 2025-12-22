// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include <string>
#include <vector>
#include "core/utils/mathUtils.h"  
#include "H5Cpp.h"
#include "core/types/Grid3D.h"
/**
 * @brief Container for a static 3D electric field on a regular grid
 * 
 * Stores precomputed electric field data from BEM/FEM solvers or analytical solutions.
 * Field components (Ex, Ey, Ez) and potential (phi) are stored in row-major order.
 * Grid axes (xs, ys, zs) define the spatial sampling points.
 * 
 * Typical usage:
 * 1. Load from HDF5 file using load_field_array()
 * 2. Validate with is_valid()
 * 3. Interpolate at ion positions with interpolate_field()
 */
struct FieldArray {
    std::vector<double> xs, ys, zs;  ///< Grid axes [m] (1D coordinate arrays)
    std::vector<double> Ex, Ey, Ez;  ///< Field components [V/m] (flattened 3D arrays)
    std::vector<double> phi;         ///< Electric potential [V] at grid points
    size_t nx = 0, ny = 0, nz = 0;   ///< Number of grid points per axis

    /** @brief Check if field array is properly initialized and consistent */
    bool is_valid() const {
        return (nx > 0 && ny > 0 && nz > 0
                && xs.size() == nx
                && ys.size() == ny
                && zs.size() == nz
                && Ex.size() == nx * ny * nz
                && Ey.size() == nx * ny * nz
                && Ez.size() == nx * ny * nz);
                // Note: phi is optional, not checked here
    }
};

/**
 * @brief Load 3D field array from HDF5 file
 * 
 * @param path Path to HDF5 file containing field data
 * 
 * @return Populated FieldArray (nx=0 on error)
 * 
 * Expected HDF5 structure:
 * - Datasets: "x", "y", "z" (1D coordinate arrays)
 * - Datasets: "Ex", "Ey", "Ez" (3D field component arrays)
 * - Dataset: "phi" (3D potential array, optional)
 * 
 * Field data typically comes from:
 * - BEM solver (Boundary Element Method)
 * - FEM solver (Finite Element Method)
 * - Analytical field calculations exported to HDF5
 * 
 * Units must be SI: coordinates in meters, fields in V/m, potential in V.
 */
FieldArray load_field_array(const std::string& path);

/**
 * @brief Interpolate electric field at arbitrary position using trilinear interpolation
 * 
 * @param field Loaded field array
 * @param pos Position (x,y,z) [m] in same coordinate system as field grid
 * 
 * @return Interpolated electric field vector [V/m]
 * 
 * Uses trilinear interpolation on the 3D grid. For positions outside the grid,
 * returns zero field (no extrapolation).
 * 
 * Interpolation weights are computed from fractional grid indices, then applied
 * to the 8 surrounding grid points (corners of interpolation cube).
 */
Vec3 interpolate_field(const FieldArray& field, const Vec3& pos);

/**
 * @brief Convert Grid3D structure to FieldArray format
 * 
 * @param grid Grid3D object (alternative grid representation)
 * 
 * @return FieldArray with copied data
 * 
 * Utility function for compatibility between different grid representations
 * used in fieldsolver and trajectory simulation modules.
 */
FieldArray grid_to_fieldarray(const Grid3D& grid);
