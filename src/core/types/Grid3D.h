// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
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
