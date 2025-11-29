// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   -------------------------------------
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        grid3D.h
 *   @brief       Defines a 3D grid structure for potentials and fields.
 *
 * @details
 * Provides a Grid3D struct encapsulating a regular 3D grid with
 * uniform spacing in each dimension. The grid stores the electric
 * potential (phi) and electric field (E) at each grid point.
 *
 *   @date        2025-10-18
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */
#pragma once
#include "core/utils/mathUtils.h"
#include <vector>

struct Grid3D {
    int Nx, Ny, Nz;
    double dx, dy, dz;
    Vec3 origin_m;

    std::vector<double> phi;  // potential [V]
    std::vector<Vec3>   E;    // electric field [V/m]

    Grid3D(int nx, int ny, int nz, double dx_, double dy_, double dz_)
        : Nx(nx), Ny(ny), Nz(nz), dx(dx_), dy(dy_), dz(dz_)
    {
        phi.resize(Nx*Ny*Nz, 0.0);
        E.resize(Nx*Ny*Nz, Vec3{0.0,0.0,0.0});
    }

    inline int index(int i, int j, int k) const {
        return i + Nx * (j + Ny * k);
    }

    inline size_t size() const { return phi.size(); }
};
