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
 *   @file        BoundaryConditions.h
 *   @brief       Defines boundary conditions for a 3D grid.
 *
 * @details
 * Provides a BoundaryConditions struct to specify Dirichlet and Neumann
 * boundary conditions for the simulation grid.
 *
 *   @date        2025-10-18
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

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

