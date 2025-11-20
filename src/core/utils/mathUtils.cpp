/**
 * =====================================================================
 *
 *   ICARION: A Modular Framework for Ion Collision and Reaction Integration
 * =====================================================================
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        mathUtils.cpp
 *   @brief       Utility functions for vector and matrix calculations in ICARION.
 *
 * @details
 * Provides:
 * - 3D vector (`Vec3`) arithmetic operations (addition, subtraction, scaling).
 *   - Dot and cross products, norm, normalization.
 *   - Multiplication of vectors by 3×3 matrices.
 *   - Element-wise operations on `std::vector<double>` (add, subtract, scale).
 *   - L2-norm calculation.
 *   - Runge–Kutta intermediate state updates for ODE integration.
 *
 *   @date        2025-10-06
 *   @version     1.0
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

#include "core/utils/mathUtils.h"

#include <cmath>
#include <random>
#include <stdexcept>

// =============================
// std::vector<double> operations
// =============================

std::vector<double> add_vectors(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) throw std::invalid_argument("Vector sizes must match");
    std::vector<double> result(a.size());
    for (size_t i=0; i<a.size(); i++) result[i] = a[i] + b[i];
    return result;
}

std::vector<double> sub_vectors(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) throw std::invalid_argument("Vector sizes must match");
    std::vector<double> result(a.size());
    for (size_t i=0; i<a.size(); i++) result[i] = a[i] - b[i];
    return result;
}

std::vector<double> scale_vector(const std::vector<double>& v, double scalar) {
    std::vector<double> result(v.size());
    for (size_t i=0; i<v.size(); i++) result[i] = v[i] * scalar;
    return result;
}

double get_vector_magnitude(const std::vector<double>& y) {
    double sum = 0.0;
    for (double val : y) sum += val*val;
    return std::sqrt(sum);
}