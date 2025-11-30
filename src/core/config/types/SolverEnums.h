// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#ifndef ICARION_CONFIG_SOLVER_ENUMS_H
#define ICARION_CONFIG_SOLVER_ENUMS_H

namespace ICARION::config {

/**
 * @brief Numerical integrator types
 * 
 * Available ODE solvers for trajectory integration.
 */
enum class SolverType {
    RK4,                ///< 4th-order Runge-Kutta (fixed step)
    RK45,               ///< Runge-Kutta-Fehlberg (adaptive step)
    Boris               ///< Boris algorithm (for strong magnetic fields)
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_SOLVER_ENUMS_H
