// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "poissonSolver.h"
#include "core/log/Logger.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>

PoissonSolver::PoissonSolver(Grid3D& grid)
    : m_grid(grid)
{
    m_rho.resize(grid.size(), 0.0);
}

void PoissonSolver::setBoundaryConditions(const BoundaryConditions& bc) {
    m_bc = bc;
}

/**
 * @brief Set the source term (charge density) for the Poisson equation.
 * @param rho Charge density array (size = Nx*Ny*Nz)
 * 
 * NOTE: Sign convention handled in solvers (all use -rho/eps0 as per Poisson equation).
 * Do NOT negate here - would cause double negation!
 */
void PoissonSolver::setSourceTerm(const std::vector<double>& rho) {
    if (rho.size() != m_grid.size())
        throw std::runtime_error("Source term size mismatch with grid dimensions.");
    
    m_rho = rho;  // Direct copy, no negation
}

/**
 * @brief Reset the solver state (potential and fields).
 */
void PoissonSolver::reset() {
    std::fill(m_grid.phi.begin(), m_grid.phi.end(), 0.0);
    std::fill(m_grid.E.begin(), m_grid.E.end(), Vec3{0.0,0.0,0.0});
}

/**
 * @brief Solve the Poisson equation using Gauss-Seidel iteration.
 * @details
 * Iteratively updates the potential φ until convergence or max iterations.
 * @param eps0 Permittivity of free space
 * @param tol  Convergence tolerance
 * @param max_iter Maximum number of iterations
 *
 */
void PoissonSolver::solve(double eps0, double tol, int max_iter) {
    size_t grid_size = m_grid.Nx * m_grid.Ny * m_grid.Nz;
    
    // Choose optimal solver based on problem size
    if (grid_size > 500000) {
        // Very large grids: Multigrid (most efficient)
        solveMultiGrid(eps0, tol, max_iter);
    } else if (grid_size > 100000) {
        // Large grids: Conjugate Gradient (faster convergence)
        solveConjugateGradient(eps0, tol, max_iter);
    } else if (grid_size > 20000) {
        // Medium grids: Red-Black (parallelizable)
        solveRedBlack(eps0, tol, max_iter);
    } else {
        // Small grids: Standard Gauss-Seidel
        solveGaussSeidel(eps0, tol, max_iter);
    }
}

// Helper: initialize phi at Dirichlet nodes if mask/values provided
static inline void initialize_dirichlet(std::vector<double>& phi, const std::vector<char>& mask,
                                       const std::vector<double>& values) {
    if (mask.empty() || values.empty()) return;
    size_t N = phi.size();
    if (mask.size() != N || values.size() != N) return;
    for (size_t i = 0; i < N; ++i) {
        if (mask[i]) phi[i] = values[i];
    }
}

/**
 * @brief Adaptive tolerance solver - looser tolerance early in simulation
 * @param time_step Current simulation time step (0 = first step)
 */
void PoissonSolver::solveAdaptive(double eps0, double final_tol, int max_iter, int time_step) {
    // Start with looser tolerance, tighten over time
    double adaptive_tol;
    if (time_step < 10) {
        adaptive_tol = std::max(final_tol * 100.0, 1e-4);  // 100x looser initially
    } else if (time_step < 50) {
        adaptive_tol = std::max(final_tol * 10.0, 1e-5);   // 10x looser
    } else {
        adaptive_tol = final_tol;  // Full precision after 50 steps
    }
    
    // Reduce max iterations for early steps
    int adaptive_max_iter = (time_step < 10) ? max_iter / 4 : max_iter;
    
    solve(eps0, adaptive_tol, adaptive_max_iter);
}

void PoissonSolver::solveGaussSeidel(double eps0, double tol, int max_iter) {
    const int Nx = m_grid.Nx, Ny = m_grid.Ny, Nz = m_grid.Nz;
    const double dx2 = m_grid.dx * m_grid.dx;
    const double dy2 = m_grid.dy * m_grid.dy;
    const double dz2 = m_grid.dz * m_grid.dz;
    const double denom = 2.0 * (1.0/dx2 + 1.0/dy2 + 1.0/dz2);
    
    // SOR relaxation parameter - optimal for 3D Poisson: ω ≈ 1.7-1.9
    const double omega = 1.8;

    // initialize Dirichlet nodes if present
    initialize_dirichlet(m_grid.phi, m_dirichlet_mask, m_dirichlet_values);

    double res = 1e9;
    int iter = 0;

    do {
        double res_sum = 0.0;

        for (int k = 1; k < Nz - 1; ++k)
            for (int j = 1; j < Ny - 1; ++j)
                for (int i = 1; i < Nx - 1; ++i) {
                        int idx = m_grid.index(i,j,k);
                        // skip fixed Dirichlet nodes
                        if (!m_dirichlet_mask.empty() && m_dirichlet_mask[idx]) continue;
                    double phi_old = m_grid.phi[idx];
                    double phi_gauss_seidel =
                        ((m_grid.phi[m_grid.index(i+1,j,k)] + m_grid.phi[m_grid.index(i-1,j,k)]) / dx2 +
                        (m_grid.phi[m_grid.index(i,j+1,k)] + m_grid.phi[m_grid.index(i,j-1,k)]) / dy2 +
                        (m_grid.phi[m_grid.index(i,j,k+1)] + m_grid.phi[m_grid.index(i,j,k-1)]) / dz2
                        + m_rho[idx] / eps0) / denom;  // Poisson: ∇²φ = -ρ/ε₀ → discretized: φ[i] involves +ρ/ε₀
                    
                    // SOR update: φ_new = φ_old + ω(φ_GS - φ_old)
                    double phi_new = phi_old + omega * (phi_gauss_seidel - phi_old);

                    res_sum += std::fabs(phi_new - phi_old);
                    m_grid.phi[idx] = phi_new;
                }

        res = res_sum / (Nx * Ny * Nz);
        ++iter;

    } while (res > tol && iter < max_iter);

    m_lastResidual = res;
    ICARION::log::debug_log(std::string("[PoissonSolver] Converged in ") + std::to_string(iter) + " iterations, residual " + std::to_string(res));
    if (iter % 100 == 0)
        ICARION::log::debug_log(std::string("[iter ") + std::to_string(iter) + "] mean(phi)=" + std::to_string(
            std::accumulate(m_grid.phi.begin(), m_grid.phi.end(), 0.0) / m_grid.phi.size()) + ", residual=" + std::to_string(res));

    if (iter == max_iter)
        std::cerr << "[PoissonSolver] Warning: reached max iterations (" << max_iter << ")\n";
}

/**
 * @brief Red-Black Gauss-Seidel solver - better for parallelization
 * @details
 * Separates grid points into red and black checkerboard pattern.
 * All red points can be updated simultaneously, then all black points.
 */
void PoissonSolver::solveRedBlack(double eps0, double tol, int max_iter) {
    const int Nx = m_grid.Nx, Ny = m_grid.Ny, Nz = m_grid.Nz;
    const double dx2 = m_grid.dx * m_grid.dx;
    const double dy2 = m_grid.dy * m_grid.dy;
    const double dz2 = m_grid.dz * m_grid.dz;
    const double denom = 2.0 * (1.0/dx2 + 1.0/dy2 + 1.0/dz2);
    const double omega = 1.8; // SOR parameter

    double res = 1e9;
    int iter = 0;

    do {
        double res_sum = 0.0;

    // Red points: (i+j+k) is even
        #pragma omp parallel for reduction(+:res_sum) collapse(3)
        for (int k = 1; k < Nz - 1; ++k) {
            for (int j = 1; j < Ny - 1; ++j) {
                for (int i = 1; i < Nx - 1; ++i) {
                    if ((i + j + k) % 2 == 0) {
                        int idx = m_grid.index(i,j,k);
                        if (!m_dirichlet_mask.empty() && m_dirichlet_mask[idx]) continue;
                        double phi_old = m_grid.phi[idx];
                        double phi_gs = ((m_grid.phi[m_grid.index(i+1,j,k)] + m_grid.phi[m_grid.index(i-1,j,k)]) / dx2 +
                                        (m_grid.phi[m_grid.index(i,j+1,k)] + m_grid.phi[m_grid.index(i,j-1,k)]) / dy2 +
                                        (m_grid.phi[m_grid.index(i,j,k+1)] + m_grid.phi[m_grid.index(i,j,k-1)]) / dz2
                                        + m_rho[idx] / eps0) / denom;  // Poisson: ∇²φ = -ρ/ε₀ → +ρ/ε₀ in discretized form
                        double phi_new = phi_old + omega * (phi_gs - phi_old);
                        res_sum += std::fabs(phi_new - phi_old);
                        m_grid.phi[idx] = phi_new;
                    }
                }
            }
        }

        // Black points: (i+j+k) is odd
        #pragma omp parallel for reduction(+:res_sum) collapse(3)
        for (int k = 1; k < Nz - 1; ++k) {
            for (int j = 1; j < Ny - 1; ++j) {
                for (int i = 1; i < Nx - 1; ++i) {
                    if ((i + j + k) % 2 == 1) {
                        int idx = m_grid.index(i,j,k);
                        if (!m_dirichlet_mask.empty() && m_dirichlet_mask[idx]) continue;
                        double phi_old = m_grid.phi[idx];
                        double phi_gs = ((m_grid.phi[m_grid.index(i+1,j,k)] + m_grid.phi[m_grid.index(i-1,j,k)]) / dx2 +
                                        (m_grid.phi[m_grid.index(i,j+1,k)] + m_grid.phi[m_grid.index(i,j-1,k)]) / dy2 +
                                        (m_grid.phi[m_grid.index(i,j,k+1)] + m_grid.phi[m_grid.index(i,j,k-1)]) / dz2
                                        + m_rho[idx] / eps0) / denom;  // Poisson: ∇²φ = -ρ/ε₀ → +ρ/ε₀ in discretized form
                        double phi_new = phi_old + omega * (phi_gs - phi_old);
                        res_sum += std::fabs(phi_new - phi_old);
                        m_grid.phi[idx] = phi_new;
                    }
                }
            }
        }

        res = res_sum / (Nx * Ny * Nz);
        ++iter;

        if (iter % 100 == 0) {
            ICARION::log::debug_log(std::string("[Red-Black iter ") + std::to_string(iter) + "] residual=" + std::to_string(res));
        }

    } while (res > tol && iter < max_iter);

    m_lastResidual = res;
    ICARION::log::debug_log(std::string("[PoissonSolver Red-Black] Converged in ") + std::to_string(iter) + " iterations, residual " + std::to_string(res));

    if (iter == max_iter)
        std::cerr << "[PoissonSolver] Warning: reached max iterations (" << max_iter << ")\n";
}

/**
 * @brief Preconditioned Conjugate Gradient solver - much faster for large systems
 * @details
 * Uses Jacobi preconditioning for better convergence.
 * Ideal for large grids where iterative methods become slow.
 */
void PoissonSolver::solveConjugateGradient(double eps0, double tol, int max_iter) {
    const int Nx = m_grid.Nx, Ny = m_grid.Ny, Nz = m_grid.Nz;
    const double dx2 = m_grid.dx * m_grid.dx;
    const double dy2 = m_grid.dy * m_grid.dy;
    const double dz2 = m_grid.dz * m_grid.dz;
    const double diag = 2.0 * (1.0/dx2 + 1.0/dy2 + 1.0/dz2);
    
    size_t N = m_grid.size();
    std::vector<double> r(N), p(N), Ap(N), z(N);
    
    // initialize Dirichlet nodes into phi
    initialize_dirichlet(m_grid.phi, m_dirichlet_mask, m_dirichlet_values);

    // Initial residual: r = b - A*x (where b = -rho/eps0, x = phi)
    #pragma omp parallel for
    for (int k = 1; k < Nz - 1; ++k) {
        for (int j = 1; j < Ny - 1; ++j) {
            for (int i = 1; i < Nx - 1; ++i) {
                int idx = m_grid.index(i,j,k);
                if (!m_dirichlet_mask.empty() && m_dirichlet_mask[idx]) {
                    r[idx] = 0.0; // Dirichlet nodes considered satisfied
                    continue;
                }
                double Ax = diag * m_grid.phi[idx] -
                    (m_grid.phi[m_grid.index(i+1,j,k)] + m_grid.phi[m_grid.index(i-1,j,k)]) / dx2 -
                    (m_grid.phi[m_grid.index(i,j+1,k)] + m_grid.phi[m_grid.index(i,j-1,k)]) / dy2 -
                    (m_grid.phi[m_grid.index(i,j,k+1)] + m_grid.phi[m_grid.index(i,j,k-1)]) / dz2;
                r[idx] = m_rho[idx] / eps0 - Ax;  // r = b - Ax, where b = ρ/ε₀ (note: Poisson is ∇²φ = -ρ/ε₀, but discretized form uses +ρ/ε₀)
            }
        }
    }
    
    // Jacobi preconditioner: z = r / diag
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        z[i] = r[i] / diag;
        p[i] = z[i];
    }
    
    double rsold = 0.0;
    #pragma omp parallel for reduction(+:rsold)
    for (size_t i = 0; i < N; ++i) {
        rsold += r[i] * z[i];
    }
    
    int iter = 0;
    double res = std::sqrt(rsold);
    
    while (res > tol && iter < max_iter) {
        // A * p
        #pragma omp parallel for
        for (int k = 1; k < Nz - 1; ++k) {
            for (int j = 1; j < Ny - 1; ++j) {
                for (int i = 1; i < Nx - 1; ++i) {
                    int idx = m_grid.index(i,j,k);
                    if (!m_dirichlet_mask.empty() && m_dirichlet_mask[idx]) { Ap[idx] = 0.0; continue; }
                    Ap[idx] = diag * p[idx] -
                        (p[m_grid.index(i+1,j,k)] + p[m_grid.index(i-1,j,k)]) / dx2 -
                        (p[m_grid.index(i,j+1,k)] + p[m_grid.index(i,j-1,k)]) / dy2 -
                        (p[m_grid.index(i,j,k+1)] + p[m_grid.index(i,j,k-1)]) / dz2;
                }
            }
        }
        
        double pAp = 0.0;
        #pragma omp parallel for reduction(+:pAp)
        for (size_t i = 0; i < N; ++i) {
            pAp += p[i] * Ap[i];
        }
        
        double alpha = rsold / pAp;
        
        // Update solution and residual
        #pragma omp parallel for
        for (size_t i = 0; i < N; ++i) {
            m_grid.phi[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
            z[i] = r[i] / diag; // Preconditioner
        }
        
        double rsnew = 0.0;
        #pragma omp parallel for reduction(+:rsnew)
        for (size_t i = 0; i < N; ++i) {
            rsnew += r[i] * z[i];
        }
        
        double beta = rsnew / rsold;
        
        #pragma omp parallel for
        for (size_t i = 0; i < N; ++i) {
            p[i] = z[i] + beta * p[i];
        }
        
        rsold = rsnew;
        res = std::sqrt(rsnew);
        ++iter;
        
        if (iter % 50 == 0) {
            ICARION::log::debug_log(std::string("[CG iter ") + std::to_string(iter) + "] residual=" + std::to_string(res));
        }
    }
    
    m_lastResidual = res;
    ICARION::log::debug_log(std::string("[PoissonSolver CG] Converged in ") + std::to_string(iter) + " iterations, residual " + std::to_string(res));
              
    if (iter == max_iter)
        std::cerr << "[PoissonSolver CG] Warning: reached max iterations (" << max_iter << ")\n";
}

/**
 * @brief Simple Multigrid V-Cycle solver - fastest for very large grids
 * @details
 * Uses geometric multigrid with restriction/prolongation operators.
 * Most efficient solver for large 3D Poisson problems.
 */
void PoissonSolver::solveMultiGrid(double eps0, double tol, int max_iter) {
    // For simplicity, fall back to CG for now
    // Full multigrid implementation requires significant infrastructure
    ICARION::log::debug_log("[Multigrid] Using CG fallback for now...");
    solveConjugateGradient(eps0, tol, max_iter);
}

/**
 * @brief Compute electric field from potential φ via central differences.
 * @details
 * E = -∇φ computed at cell centers.
 * Simply loops over grid points and calculates E components from derivatives.
 */
void PoissonSolver::computeElectricField() {
    const int Nx = m_grid.Nx, Ny = m_grid.Ny, Nz = m_grid.Nz;
    const double dx = m_grid.dx, dy = m_grid.dy, dz = m_grid.dz;

    for (int k = 1; k < Nz - 1; ++k) {
        for (int j = 1; j < Ny - 1; ++j) {
            for (int i = 1; i < Nx - 1; ++i) {
                int idx = m_grid.index(i,j,k);
                double dphidx = (m_grid.phi[m_grid.index(i+1,j,k)] - m_grid.phi[m_grid.index(i-1,j,k)]) / (2*dx);
                double dphidy = (m_grid.phi[m_grid.index(i,j+1,k)] - m_grid.phi[m_grid.index(i,j-1,k)]) / (2*dy);
                double dphidz = (m_grid.phi[m_grid.index(i,j,k+1)] - m_grid.phi[m_grid.index(i,j,k-1)]) / (2*dz);
                // E = -∇φ (standard convention)
                m_grid.E[idx] = Vec3{-dphidx, -dphidy, -dphidz};
            }
        }
    }
}
