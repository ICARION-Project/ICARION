// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#pragma once
#include <vector>

/**
 * @brief Placeholder for boundary condition definitions.
 * Extend later for Dirichlet / Neumann / periodic boundaries.
 */
struct BoundaryConditions {
    enum class Type { Dirichlet, Neumann };
    Type type = Type::Dirichlet;
    double value = 0.0;

    static BoundaryConditions Dirichlet(double val) {
        BoundaryConditions bc;
        bc.type = Type::Dirichlet;
        bc.value = val;
        return bc;
    }

    static BoundaryConditions Neumann(double val) {
        BoundaryConditions bc;
        bc.type = Type::Neumann;
        bc.value = val;
        return bc;
    }
};

