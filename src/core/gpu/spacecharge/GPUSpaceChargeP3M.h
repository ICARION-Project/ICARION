// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file GPUSpaceChargeP3M.h
 * @brief GPU-accelerated space charge field calculation using P³M algorithm
 * 
 * Experimental prototype for GPU P³M. Not wired into SimulationEngine; not
 * validated against CPU results; geometry handling and boundary conditions are
 * incomplete. Avoid using for production until integration and validation are
 * completed.
 * 
 * **Algorithm:**
 * 1. Particle-to-Grid (P2G): Scatter ion charges to 3D grid (CIC or NGP)
 * 2. Poisson Solver: Compute potential via FFT (cuFFT)
 * 3. Grid-to-Particle (G2P): Interpolate E-field to ion positions
 * 
 * Previous speed claims were speculative; no validated benchmarks exist.
 */

#pragma once

#ifdef ICARION_USE_GPU

#include "core/types/Vec3.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"
#include "core/gpu/core/GPUContext.h"
#include <vector>
#include <memory>
#include <cufft.h>

namespace icarion {
namespace gpu {

/**
 * @brief GPU-accelerated P³M space charge solver
 * 
 * **Status:** Experimental helper only; not used anywhere by default.
 * **When to use:** Prefer CPU direct until this path is validated.
 * 
 * **Limitations:**
 * - Rectangular domain only; no cylindrical/Orbitrap geometry support.
 * - Boundary conditions and masking incomplete; periodic/Dirichlet not verified.
 * - P2G/G2P mapping and charge assignment not tested against CPU.
 * - Fixed grid resolution (no adaptive mesh refinement).
 * 
 * **Future extensions:**
 * - FMM algorithm for non-uniform ion distributions
 * - PME (Particle-Mesh-Ewald) for periodic boundaries
 * - Multi-GPU via domain decomposition
 */
class GPUSpaceChargeP3M {
public:
    /**
     * @brief P²G (Particle-to-Grid) interpolation method
     */
    enum class InterpolationMethod {
        NGP,  ///< Nearest Grid Point (0th order, fast but noisy)
        CIC,  ///< Cloud-in-Cell (1st order, smooth)
        TSC   ///< Triangular-Shaped Cloud (2nd order, very smooth but slower)
    };
    
    /**
     * @brief Configuration for P³M solver
     */
    struct Config {
        int grid_nx = 64;        ///< Grid points in x-direction
        int grid_ny = 64;        ///< Grid points in y-direction
        int grid_nz = 64;        ///< Grid points in z-direction
        
        Vec3 domain_min{0, 0, 0};    ///< Domain lower bounds [m]
        Vec3 domain_max{0.1, 0.1, 0.1}; ///< Domain upper bounds [m]
        
        InterpolationMethod p2g_method = InterpolationMethod::CIC;  ///< P²G interpolation
        InterpolationMethod g2p_method = InterpolationMethod::CIC;  ///< G²P interpolation
        
        double epsilon_0 = 8.854187817e-12;  ///< Vacuum permittivity [F/m]
        
        bool use_direct_cutoff = false;  ///< Add direct summation for nearby particles
        double direct_cutoff_radius = 0.0;  ///< Cutoff radius for direct sum [m]
        
        /**
         * @brief Validate configuration
         * @return Error message if invalid, empty string if OK
         */
        std::string validate() const;
    };
    
    /**
     * @brief Create P³M solver
     * @param context GPU context
     * @param config Solver configuration
     * @return Unique pointer to solver, nullptr on failure
     */
    static std::unique_ptr<GPUSpaceChargeP3M> create(
        const GPUContext& context,
        const Config& config
    );
    
    /**
     * @brief Destructor - frees GPU resources
     */
    ~GPUSpaceChargeP3M();
    
    // No copy/move (manages GPU resources)
    GPUSpaceChargeP3M(const GPUSpaceChargeP3M&) = delete;
    GPUSpaceChargeP3M& operator=(const GPUSpaceChargeP3M&) = delete;
    GPUSpaceChargeP3M(GPUSpaceChargeP3M&&) = delete;
    GPUSpaceChargeP3M& operator=(GPUSpaceChargeP3M&&) = delete;
    
    /**
     * @brief Compute space charge field for ion ensemble
     * 
     * @param ions Ion ensemble (CPU memory)
     * @param E_field_out Output E-field per ion [V/m] (CPU memory, must be sized)
     * @return true on success, false on error
     * 
     * **Workflow:**
     * 1. Upload ion positions/charges to GPU
     * 2. P²G: Scatter charges to grid
     * 3. FFT: Transform charge density to Fourier space
     * 4. Solve Poisson equation in Fourier space
     * 5. IFFT: Transform potential back to real space
     * 6. Compute E-field: E = -∇φ (finite differences)
     * 7. G²P: Interpolate E-field to ion positions
     * 8. Download E-fields to CPU
     * 
     * **Error handling:**
     * - Returns false on CUDA/cuFFT errors
     * - Prints error message to stderr
     * - Safe to retry after error
     */
    bool compute_space_charge_field(
        const std::vector<IonState>& ions,
        std::vector<Vec3>& E_field_out
    );

    /**
     * @brief SoA wrapper for space charge field computation
     */
    bool compute_space_charge_field(
        const ICARION::core::IonEnsemble& ions,
        std::vector<Vec3>& E_field_out
    );

private:
    bool compute_space_charge_field_raw(
        const std::vector<Vec3>& positions,
        const std::vector<double>& charges,
        std::vector<Vec3>& E_field_out
    );

public:
    /**
     * @brief Performance statistics
     */
    struct Stats {
        size_t total_calls = 0;              ///< Total field computations
        double total_time_ms = 0.0;          ///< Total GPU time [ms]
        double avg_time_ms = 0.0;            ///< Average time per call [ms]
        double last_p2g_time_ms = 0.0;       ///< Last P²G time [ms]
        double last_fft_time_ms = 0.0;       ///< Last FFT time [ms]
        double last_poisson_time_ms = 0.0;   ///< Last Poisson solve time [ms]
        double last_g2p_time_ms = 0.0;       ///< Last G²P time [ms]
        size_t last_n_ions = 0;              ///< Ion count in last call
        
        /**
         * @brief Estimate speedup vs CPU direct summation
         * @return Speedup factor (e.g., 100× means 100 times faster)
         */
        double speedup_vs_direct_cpu() const;
    };
    
    /**
     * @brief Get performance statistics
     */
    const Stats& get_stats() const { return stats_; }
    
private:
    /**
     * @brief Reset performance counters
     */
    void reset_stats() { stats_ = Stats{}; }
    
private:
    /**
     * @brief Private constructor (use create() factory)
     */
    GPUSpaceChargeP3M(const GPUContext& context, const Config& config);
    
    /**
     * @brief Initialize GPU resources (grids, FFT plan)
     * @return true on success
     */
    bool initialize();
    
    /**
     * @brief Free GPU resources
     */
    void cleanup();
    
    // Configuration
    const GPUContext& context_;
    Config config_;
    
    // GPU buffers
    Vec3* d_ion_positions_ = nullptr;     ///< Device: ion positions [N]
    double* d_ion_charges_ = nullptr;     ///< Device: ion charges [N]
    Vec3* d_ion_E_fields_ = nullptr;      ///< Device: output E-fields [N]
    
    double* d_charge_density_ = nullptr;  ///< Device: charge grid [nx×ny×nz]
    cufftDoubleComplex* d_potential_fft_ = nullptr;  ///< Device: FFT of potential [nx×ny×(nz/2+1)]
    double* d_potential_ = nullptr;       ///< Device: potential grid [nx×ny×nz]
    double* d_Ex_ = nullptr;              ///< Device: E_x grid [nx×ny×nz]
    double* d_Ey_ = nullptr;              ///< Device: E_y grid [nx×ny×nz]
    double* d_Ez_ = nullptr;              ///< Device: E_z grid [nx×ny×nz]
    
    // cuFFT plan handles
    cufftHandle fft_plan_forward_ = 0;    ///< R2C FFT plan
    cufftHandle fft_plan_inverse_ = 0;    ///< C2R FFT plan
    
    // Grid geometry (device constant memory)
    float* d_inv_spacing_ = nullptr;      ///< 1 / cell_size [3 floats]
    
    // Statistics
    mutable Stats stats_;
    
    // Initialization flag
    bool initialized_ = false;
    
    // Kernel launch wrappers (implemented in .cu file)
    void launch_p2g_cic_kernel(
        const Vec3* d_positions,
        const double* d_charges,
        size_t n_ions,
        double* d_charge_density
    );
    
    void launch_poisson_solve_fourier_kernel(
        cufftDoubleComplex* d_rho_hat,
        cufftDoubleComplex* d_phi_hat,
        double epsilon_0
    );
    
    void launch_compute_E_field_kernel(
        const double* d_potential,
        double* d_Ex,
        double* d_Ey,
        double* d_Ez
    );
    
    void launch_g2p_cic_kernel(
        const Vec3* d_positions,
        size_t n_ions,
        const double* d_Ex,
        const double* d_Ey,
        const double* d_Ez,
        Vec3* d_E_fields_out
    );
    
    void launch_scale_potential_kernel(
        double* d_potential,
        double scale_factor
    );
};

} // namespace gpu
} // namespace icarion

#endif // ICARION_USE_GPU
