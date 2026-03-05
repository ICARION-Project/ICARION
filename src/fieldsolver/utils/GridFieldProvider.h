// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
/**
 * @file GridFieldProvider.h
 * @brief Grid-based field provider using precomputed field arrays
 */
#pragma once

#include "IFieldProvider.h"
#include "core/physics/fields/fieldSampling.h"
#include "core/io/fieldArrayLoader.h"

// Convert between Vec3 and Vec3d
inline ICARION::fields::Vec3d to_Vec3d(const Vec3& v) {
    return ICARION::fields::Vec3d(v.x, v.y, v.z);
}

inline Vec3 from_Vec3d(const ICARION::fields::Vec3d& v) {
    return Vec3{v.x, v.y, v.z};
}

/**
 * @class GridFieldProvider
 * @brief Field provider using trilinear interpolation on regular 3D grid
 * 
 * Wraps precomputed electric field data from:
 * - FieldArray: Loaded from HDF5 files (BEM/FEM solver output)
 * - FieldServer snapshots: Live field updates (for adaptive meshes - Full build only)
 * 
 * Provides fast field evaluation via trilinear interpolation.
 * Suitable for static fields or slowly-varying time-dependent fields.
 */
class GridFieldProvider : public IFieldProvider {
public:
    /**
     * @brief Construct from FieldServer snapshot (for adaptive field updates - Full build only)
     * @param snapshot Pointer to FieldServer snapshot (must remain valid during usage)
     */
    GridFieldProvider(const ICARION::fields::FieldSnapshot* snapshot) 
        : snapshot_(snapshot), fld_(nullptr) {}

    /**
     * @brief Construct from preloaded FieldArray (most common usage)
     * @param fld Pointer to FieldArray with Ex, Ey, Ez data (must remain valid)
     */
    GridFieldProvider(const FieldArray* fld) : snapshot_(nullptr), fld_(fld) {}

    /**
     * @brief Evaluate electric field at position using trilinear interpolation
     * @param pos Position [m] in simulation domain
     * @return Electric field E [V/m] interpolated from grid data
     * 
     * Interpolation steps:
     * 1. Find enclosing grid cell from position
     * 2. Compute fractional coordinates within cell
     * 3. Interpolate E-field from 8 surrounding grid points
     * 4. Return zero field if position is outside grid bounds
     */
    Vec3 get_E(const Vec3& pos) const override {
        if (snapshot_) {
            auto ed = ICARION::fields::sample_field(*snapshot_, to_Vec3d(pos));
            return from_Vec3d(ed);
        }
        if (fld_) {
            return interpolate_field(*fld_, pos);
        }
        return Vec3{0.0,0.0,0.0};
    }

    /**
     * @brief Get potential at position
     * @param pos Position [m]
     * @return 0.0 (potential not stored in current FieldArray format)
     * 
     * FieldArray stores E-field components but not potential.
     * Potential reconstruction would require line integration of E-field.
     */
    double get_phi(const Vec3& pos) const override {
        (void)pos;
        return 0.0;
    }

    /**
     * @brief Get underlying FieldArray (for GPU upload)
     * @return Pointer to FieldArray or nullptr if using FieldSnapshot
     * 
     * Used by GPU integration to extract field data for texture upload.
     * Returns nullptr if provider was constructed from FieldSnapshot.
     */
    const FieldArray* get_field_array() const {
        return fld_;
    }

private:
    const ICARION::fields::FieldSnapshot* snapshot_ = nullptr;  ///< Live field server snapshot (optional)
    const FieldArray* fld_ = nullptr;          ///< Preloaded field array (alternative source)
};
