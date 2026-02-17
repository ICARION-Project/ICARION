// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
/**
 * @file IFieldProvider.h
 * @brief Interface for interchangeable electric field providers
 * 
 * Defines abstract base class for field evaluation at arbitrary positions.
 * Implementations include grid-based interpolation (GridFieldProvider) and
 * BEM surface charge integration (BEMFieldProvider).
 */
#pragma once

#include "core/types/Vec3.h"

/**
 * @class IFieldProvider
 * @brief Abstract interface for electric field evaluation
 * 
 * Provides unified API for different field sources:
 * - Precomputed field arrays (HDF5 from BEM/FEM solvers)
 * - BEM surface charge distributions (direct integration)
 * - Analytical field expressions
 * - Space-charge fields from Poisson solver
 * 
 * Enables runtime polymorphism for field evaluation in integrator.
 */
class IFieldProvider {
public:
    virtual ~IFieldProvider() = default;
    
    /**
     * @brief Evaluate electric field at given position
     * @param pos Position [m] in simulation domain
     * @return Electric field E [V/m] at position
     * 
     * Pure virtual function - must be implemented by derived classes.
     * May use interpolation (grid-based) or direct calculation (BEM).
     */
    virtual Vec3 get_E(const Vec3& pos) const = 0;
    
    /**
     * @brief Evaluate time-dependent electric field at given position
     * @param pos Position [m] in simulation domain
     * @param t Current simulation time [s]
     * @return Electric field E [V/m] at position and time
     * 
     * Optional method for time-varying fields (RF modulation, time-varying voltages).
     * Default implementation ignores time and calls get_E(pos).
     * Override in derived classes for true time-dependent fields.
     */
    virtual Vec3 get_E(const Vec3& pos, double t) const {
        (void)t;
        return get_E(pos);
    }
    
    /**
     * @brief Evaluate electric potential at given position
     * @param pos Position [m] in simulation domain
     * @return Electric potential φ [V] at position
     * 
     * Optional method - default returns 0.0.
     * Useful for diagnostics and energy calculations.
     * May not be available for all field providers (e.g., if only E-field is stored).
     */
    virtual double get_phi(const Vec3& pos) const { (void)pos; return 0.0; }
};
