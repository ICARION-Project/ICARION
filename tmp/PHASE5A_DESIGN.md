# Phase 5A: SimulationEngine Design Document

**Date:** 2025-11-22  
**Branch:** `refactor/simulation-engine`  
**Status:** 🔵 In Progress

---

## 🎯 Goals

1. **Clean Orchestrator:** Replace monolithic `integrate_trajectory()` with modular `SimulationEngine`
2. **Dependency Injection:** All physics modules injected (ForceRegistry, CollisionHandler, ReactionHandler, IntegrationStrategy)
3. **SSOT Compliance:** Use `FullConfig` (no legacy `GlobalParams` dependencies)
4. **Separation of Concerns:** Separate physics, I/O, and domain management
5. **Testability:** Each component unit-testable in isolation

---

## 🏗️ Architecture Overview

```
SimulationEngine (orchestrator)
    │
    ├── ForceRegistry (compute total force)
    ├── IntegrationStrategy (RK4/RK45/Boris)
    ├── CollisionHandler (EHSS/HSS/Langevin/Friction)
    ├── ReactionHandler (ion-molecule reactions)
    ├── OutputManager (HDF5 trajectory writer)
    └── DomainManager (boundary checks, transitions)
```

---

## 📦 Component Specifications

### **1. SimulationEngine**

**Responsibilities:**
- Orchestrate simulation loop (birth, integration, output)
- Manage global time progression (`t_global`)
- Handle ion birth events
- Coordinate parallel execution (OpenMP)
- Exit detection (all ions inactive)

**API Design:**

```cpp
namespace ICARION::integrator {

class SimulationEngine {
public:
    /**
     * @brief Construct simulation engine with full configuration
     * @param config Complete simulation settings (SSOT)
     */
    explicit SimulationEngine(const config::FullConfig& config);
    
    /**
     * @brief Run complete simulation
     * @param ions Initial ion ensemble (modified in-place)
     * @return Final ion states after simulation
     */
    std::vector<IonState> run(std::vector<IonState>& ions);
    
    /**
     * @brief Get simulation statistics
     */
    SimulationStats get_stats() const;
    
private:
    const config::FullConfig& config_;
    
    // Physics modules (injected)
    std::unique_ptr<physics::ForceRegistry> force_registry_;
    std::unique_ptr<strategies::IIntegrationStrategy> integration_strategy_;
    std::unique_ptr<physics::ICollisionHandler> collision_handler_;
    std::unique_ptr<physics::IReactionHandler> reaction_handler_;
    
    // I/O and domain management
    std::unique_ptr<OutputManager> output_manager_;
    std::unique_ptr<DomainManager> domain_manager_;
    
    // RNG pool (one per ion)
    std::vector<EhssRng> rng_pool_;
    
    // Simulation state
    double t_start_;
    double t_end_;
    double dt_;
    
    // Helper methods
    void initialize_physics_modules();
    void handle_ion_birth(std::vector<IonState>& ions, double t_global);
    void process_ion(IonState& ion, double t_global, EhssRng& rng, 
                     const std::vector<IonState>& all_ions);
    bool all_ions_exited(const std::vector<IonState>& ions) const;
};

} // namespace ICARION::integrator
```

**Key Design Decisions:**
- **Constructor:** Takes `FullConfig` (SSOT), initializes all subsystems
- **run():** Main simulation loop (analogous to old `integrate_trajectory()`)
- **Private helpers:** `process_ion()` handles single ion integration step
- **RNG pool:** One per ion (avoids race conditions in OpenMP)

---

### **2. OutputManager**

**Responsibilities:**
- Wrap HDF5Writer v2 for trajectory output
- Manage write buffers (time-based or size-based flushing)
- Handle metadata (species, parameters)
- Periodic output during simulation

**API Design:**

```cpp
namespace ICARION::integrator {

class OutputManager {
public:
    /**
     * @brief Construct output manager
     * @param hdf5_filename HDF5 trajectory output file path
     * @param log_filename Text log file path (optional, nullptr = no text log)
     * @param write_interval_dt Time interval between HDF5 writes [s]
     * @param buffer_max Max timesteps in RAM before flush
     */
    OutputManager(const std::string& hdf5_filename,
                  const std::string* log_filename,
                  double write_interval_dt, 
                  size_t buffer_max = 50);
    
    /**
     * @brief Initialize HDF5 file and write metadata
     * @param config Simulation configuration (for parameter export)
     * @param ions Initial ion ensemble (for species metadata)
     */
    void initialize(const config::FullConfig& config, 
                    const std::vector<IonState>& ions);
    
    /**
     * @brief Log progress message to text log (optional)
     * @param message Progress message (e.g., "50% completed")
     */
    void log_progress(const std::string& message);
    
    /**
     * @brief Log timestep snapshot (buffers in RAM)
     * @param t Current time [s]
     * @param ions Current ion states
     */
    void log_step(double t, const std::vector<IonState>& ions);
    
    /**
     * @brief Log progress message (to text log)
     * @param message Progress message (e.g., "50% completed")
     */
    void log_progress(const std::string& message);
    
    /**
     * @brief Check if write is needed (based on time or buffer size)
     * @param t_current Current simulation time [s]
     * @return true if flush should be triggered
     */
    bool should_write(double t_current) const;
    
    /**
     * @brief Flush buffers to HDF5 file
     */
    void flush();
    
    /**
     * @brief Finalize output (write remaining data, close file)
     */
    void finalize(double t_final, const std::vector<IonState>& final_ions);
    
private:
    std::string hdf5_filename_;
    double write_interval_dt_;
    double next_write_time_;
    size_t buffer_max_;
    
    // HDF5 trajectory buffers
    std::vector<double> times_buffer_;
    std::vector<std::vector<IonState>> trajectory_buffer_;
    
    // Text log output (optional)
    std::unique_ptr<io::RunLogger> text_logger_;
};

} // namespace ICARION::integrator
```

**Key Design Decisions:**
- **Unified output:** Manages both HDF5 trajectory data AND text logging
- **Buffered writes:** Avoid HDF5 overhead on every timestep
- **Time-based + size-based:** Flush when interval reached OR buffer full
- **Metadata on init:** Species and params written once at start
- **Finalize on exit:** Ensures last timestep always written + completion summary
- **Optional text log:** Can disable RunLogger for performance (pass nullptr)

---

### **3. DomainManager**

**Responsibilities:**
- Find which domain an ion is in (spatial lookup)
- Handle coordinate transforms (global ↔ local)
- Check aperture crossings (domain transitions)
- Apply boundary conditions

**API Design:**

```cpp
namespace ICARION::integrator {

class DomainManager {
public:
    /**
     * @brief Construct domain manager
     * @param domains Vector of instrument domains
     */
    explicit DomainManager(const std::vector<InstrumentDomain>& domains);
    
    /**
     * @brief Find which domain contains the ion position
     * @param pos Global position [m]
     * @return Domain index (-1 if outside all domains)
     */
    int find_domain_index(const Vec3& pos) const;
    
    /**
     * @brief Get domain by index
     * @param idx Domain index
     * @return Reference to domain (SSOT: config::DomainConfig)
     */
    const config::DomainConfig& get_domain(int idx) const;
    
    /**
     * @brief Transform position from global to local coordinates
     * @param pos Global position [m]
     * @param domain_idx Domain index
     * @return Local position [m]
     */
    Vec3 global_to_local_pos(const Vec3& pos, int domain_idx) const;
    
    /**
     * @brief Transform velocity from global to local coordinates
     * @param vel Global velocity [m/s]
     * @param domain_idx Domain index
     * @return Local velocity [m/s]
     */
    Vec3 global_to_local_vel(const Vec3& vel, int domain_idx) const;
    
    /**
     * @brief Transform from local to global coordinates
     */
    Vec3 local_to_global_pos(const Vec3& pos_local, int domain_idx) const;
    Vec3 local_to_global_vel(const Vec3& vel_local, int domain_idx) const;
    
    /**
     * @brief Check if ion crossed aperture when transitioning domains
     * @param ion Ion state (modified if crossing detected)
     * @param pos_before Position before integration step [m]
     * @param pos_after Position after integration step [m]
     */
    void check_aperture_crossing(IonState& ion, const Vec3& pos_before, 
                                  const Vec3& pos_after);
    
    /**
     * @brief Update ion domain-specific properties (temperature, pressure, etc.)
     * @param ion Ion state (modified)
     * @param domain_idx New domain index
     */
    void update_domain_properties(IonState& ion, int domain_idx);
    
private:
    const std::vector<InstrumentDomain>& domains_;
    
    // Helper: check if position is inside domain bounding box
    bool is_inside_domain(const Vec3& pos, const InstrumentDomain& domain) const;
};

} // namespace ICARION::integrator
```

**Key Design Decisions:**
- **Const reference to domains:** No ownership (owned by FullConfig)
- **Coordinate transforms:** Explicit API for readability
- **Aperture crossing:** Deactivates ions that miss transition aperture
- **Domain properties:** Updates temperature, pressure, gas velocity on transition

---

## 🔄 Integration with Existing Systems

### **ForceRegistry Integration**

**Current State (Phase 1):**
- `ForceRegistry` computes `F_total` via `compute_total_force()`
- Returns `Vec3` force in Newtons

**Integration Strategy:**
```cpp
// Inside SimulationEngine::process_ion()
auto compute_accel = [&](const IonState& ion, double t) -> Vec3 {
    // Build force context (position, velocity, time, domain)
    physics::ForceContext ctx;
    ctx.domain_config = &domain_manager_->get_domain(ion.current_domain_index).dom;
    ctx.t = t;
    // Add more context as needed (space charge grid, etc.)
    
    Vec3 F_total = force_registry_->compute_total_force(ion, ctx);
    return F_total / ion.mass_kg;  // Convert to acceleration
};

// Pass to integration strategy
integration_strategy_->step(ion, t, dt, compute_accel, &domain_config);
```

### **IntegrationStrategy Selection**

**Runtime Selection:**
```cpp
// In SimulationEngine::initialize_physics_modules()
const auto& sim_config = config_.simulation;
std::string solver_type = sim_config.integration_method;  // "RK4", "RK45", "Boris"

integration_strategy_ = strategies::IntegrationStrategyFactory::create(
    solver_type, &sim_config.domain);
```

### **CollisionHandler Integration**

**Phase 2 API:**
```cpp
// Inside SimulationEngine::run() - Main loop with progress logging
size_t step_count = 0;
while (t_global < t_end_) {
    // Ion processing...
    
    // HDF5 output
    if (output_manager_->should_write(t_global)) {
        output_manager_->log_step(t_global, ions);
        output_manager_->flush();
    }
    
    // Progress logging (every 1000 steps)
    if (step_count % 1000 == 0) {
        double progress = (t_global - t_start_) / (t_end_ - t_start_) * 100.0;
        std::ostringstream msg;
        msg << std::fixed << std::setprecision(1)
            << progress << "% completed (t = " << t_global * 1000.0 << " ms)";
        output_manager_->log_progress(msg.str());
    }
    
    t_global += dt_;
    step_count++;
}

// Inside SimulationEngine::process_ion() - Collision handling
if (collision_handler_) {
    config::EnvironmentConfig env;
    env.temperature_K = ion.domain_temperature_K;
    env.pressure_Pa = ion.domain_particle_density_m3 * BOLTZMANN_CONSTANT * ion.domain_temperature_K;
    env.gas_species = domain_manager_->get_domain(ion.current_domain_index).environment.gas_species;
    env.gas_velocity_m_s = ion.domain_gas_velocity_m_s;
    env.compute_derived_properties();
    
    collision_handler_->handle_collision(ion, dt, rng, env);
}
```

### **ReactionHandler Integration**

**Phase 3 API:**
```cpp
// In SimulationEngine::process_ion()
if (reaction_handler_) {
    physics::ReactionContext react_ctx;
    react_ctx.dt = dt;
    react_ctx.temperature_K = ion.domain_temperature_K;
    react_ctx.pressure_Pa = ion.domain_particle_density_m3 * BOLTZMANN_CONSTANT * ion.domain_temperature_K;
    
    reaction_handler_->handle_reaction(ion, rng, react_ctx);
}
```

---

## 🧪 Testing Strategy

### **Unit Tests**

**test_output_manager.cpp:**
- ✅ Buffer management (fill, flush, overflow)
- ✅ Time-based write triggering
- ✅ Metadata writing (species, params)
- ✅ Finalize writes remaining data

**test_domain_manager.cpp:**
- ✅ Domain lookup (inside/outside domains)
- ✅ Coordinate transforms (global ↔ local)
- ✅ Aperture crossing detection
- ✅ Domain property updates

**test_simulation_engine.cpp:**
- ✅ Ion birth timing
- ✅ Single-ion integration (with mocked modules)
- ✅ Parallel execution (OpenMP correctness)
- ✅ Exit detection (all ions inactive)
- ✅ Full simulation (simplified config)

### **Integration Tests**

**Test Scenario:**
```cpp
// tests/integrator/test_simulation_engine_integration.cpp
TEST(SimulationEngineIntegration, CompareWithLegacy) {
    // Load test config
    auto config = load_test_config("examples/ims_basic.json");
    
    // Run with SimulationEngine
    SimulationEngine engine(config);
    auto ions_new = config.generate_ions(rng);
    auto final_new = engine.run(ions_new);
    
    // Run with legacy integrate_trajectory()
    auto ions_old = config.generate_ions(rng);
    auto final_old = integrate_trajectory(ions_old, ...);
    
    // Compare results (positions, velocities, arrival times)
    EXPECT_NEAR(final_new[0].pos.z, final_old[0].pos.z, 1e-6);
    // ... more comparisons
}
```

### **Performance Benchmarks**

**Acceptance Criteria:**
- ✅ CPU time within 5% of legacy `integrate_trajectory()`
- ✅ Memory overhead <10%
- ✅ HDF5 output identical to legacy

---

## 📅 Implementation Plan

### **Step 1: DomainManager** ✅ COMPLETE (Nov 22, 2025)
- ✅ Extract domain lookup, coordinate transforms from legacy code
- ✅ Implement aperture crossing logic
- ✅ Write unit tests (test_domain_manager.cpp - 8 tests, 100% passing)
- ✅ SSOT migration: Uses config::DomainConfig (not InstrumentDomain)

### **Step 2: OutputManager** (5 hours)
- Implement buffered HDF5 output with time-based triggering
- Wire up HDF5Writer v2 API
- Integrate RunLogger for text logging (progress, summaries)
- Write unit tests (test_output_manager.cpp)

### **Step 3: SimulationEngine Core** (8 hours)
- Implement constructor, module initialization
- Implement `process_ion()` with physics integration
- Implement main `run()` loop (ion birth, exit detection)
- Wire up ForceRegistry, IntegrationStrategy, CollisionHandler, ReactionHandler

### **Step 4: Testing & Validation** (6 hours)
- Write SimulationEngine unit tests
- Integration test against legacy code
- Performance benchmarks
- Fix any issues, optimize hotspots

**Total Estimated Time:** ~22 hours (~3 days)

---

## 🚧 Migration Notes

### **Legacy Code Removal**

**Files to be deprecated (Phase 6):**
- `src/core/integrator/integrator.cpp` (old `integrate_trajectory()`)
- `src/core/integrator/integrator_helpers.cpp` (old `integrate_one_step()`)

**Temporary Coexistence:**
- Both `SimulationEngine` and legacy `integrate_trajectory()` will exist during Phase 5
- `main.cpp` still uses legacy (Phase 6 will switch)

### **FullConfig vs GlobalParams**

**SSOT Enforcement:**
- SimulationEngine uses **only** `FullConfig`
- Legacy code still uses `GlobalParams` (temporary)
- Phase 6 migration will remove `GlobalParams` entirely

---

## 📊 Success Metrics

- ✅ All unit tests pass (90%+ coverage)
- ✅ Integration test matches legacy output (bit-exact for deterministic cases)
- ✅ Performance within 5% of legacy
- ✅ Code is modular (can swap physics modules via DI)
- ✅ No global state (all state in SimulationEngine)

---

**Next Steps:**
1. Start with DomainManager implementation
2. Then OutputManager
3. Finally SimulationEngine core
4. Test, benchmark, iterate

**Status:** Design complete, ready to implement 🚀
