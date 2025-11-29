# ICARION main.cpp Refactoring Plan

**Date:** 2024-01-XX  
**Status:** 🟡 PLANNING (NOT YET IMPLEMENTED)  
**Target:** Reduce main.cpp from 557 lines to <100 lines by extracting setup helpers

---

## 1. Executive Summary

### Current State
- **File:** `src/main/main.cpp`
- **Line Count:** 557 lines in single function
- **Problem:** Monolithic main() function mixing:
  - CLI parsing & validation (100+ lines)
  - Configuration loading & overrides (50+ lines)
  - Logging setup (30+ lines)
  - Ion generation (80+ lines)
  - Physics module creation (150+ lines)
  - Special modes (--validate-config, --dry-run, --dump-* flags) (100+ lines)

### Proposed Solution
**Extract setup helpers into `src/main/setup/` directory:**
```
src/main/setup/
├── CliHandler.h/cpp          - Command-line parsing & special modes
├── LoggingSetup.h/cpp        - Logging initialization
├── ConfigSetup.h/cpp         - Config loading, validation, overrides
├── IonSetup.h/cpp            - Ion generation from config
├── PhysicsSetup.h/cpp        - Force registries, collision/reaction handlers
└── EngineSetup.h/cpp         - SimulationEngine initialization
```

### Benefits
- ✅ **Maintainability:** Each setup phase in separate file (easier to modify)
- ✅ **Testability:** Unit tests for ion generation, physics setup logic
- ✅ **Readability:** main() becomes high-level orchestration (<100 lines)
- ✅ **Reusability:** Setup helpers usable by tests, benchmarks, GUI apps

### Risk Assessment
- **Risk Level:** 🟢 LOW
- **Complexity:** Medium (pure refactoring, no logic changes)
- **Testing:** Existing CTest suite validates correctness (48/48 must pass)
- **Estimated Time:** ~2-3 hours

---

## 2. Current main.cpp Structure Analysis

### Section Breakdown (557 lines)
```cpp
Lines    | Section                              | Complexity | Extract To
---------|--------------------------------------|------------|------------------
1-20     | Includes & namespace                 | Simple     | (keep in main.cpp)
21-60    | CLI Options struct                   | Simple     | setup/CliHandler.h
61-150   | CLI parsing (--help, --version, etc.)| High       | setup/CliHandler.cpp
151-180  | Profiler & logging setup             | Medium     | setup/LoggingSetup.cpp
181-230  | Special modes (--validate-config, etc.)| High     | setup/CliHandler.cpp
231-280  | Config loading & overrides           | Medium     | setup/ConfigSetup.cpp
281-340  | Ion generation                       | High       | setup/IonSetup.cpp
341-480  | Physics module setup (forces, etc.)  | High       | setup/PhysicsSetup.cpp
481-520  | Engine creation & simulation run     | Medium     | setup/EngineSetup.cpp
521-557  | Results reporting & cleanup          | Simple     | (keep in main.cpp)
```

### Key Observations
1. **Special Modes Complexity:** Lines 181-230 handle 7+ special flags (--validate-config, --dump-build-info, --dump-hdf5-schema, --list-collision-models, --dry-run, --help, --version)
   - Each mode exits early with return 0
   - Makes control flow hard to follow

2. **Physics Setup Complexity:** Lines 341-480 (~140 lines)
   - Creates ForceRegistry per domain
   - Auto-selects space charge method (Direct vs Grid based on N)
   - Loads molecular geometries for EHSS collision model
   - Creates integration strategy, collision handler, reaction handler
   - **This is the most complex section** (needs careful extraction)

3. **Ion Generation Complexity:** Lines 281-340 (~60 lines)
   - Handles ion_clouds vs inline ion definitions
   - Loads HDF5 ion clouds from files
   - Generates ions from species specs
   - Validates ion positions vs domain bounds

---

## 3. Proposed Helper Classes

### 3.1 CliHandler (CLI Parsing & Special Modes)
**File:** `src/main/setup/CliHandler.h/cpp`  
**Purpose:** Parse command-line arguments, handle special modes (--help, --validate-config, etc.)

```cpp
// CliHandler.h
#pragma once
#include <string>
#include <optional>
#include <vector>

namespace icarion::setup {

struct CliOptions {
    std::string config_file;
    bool benchmark = false;
    bool profile = false;
    std::optional<std::string> profile_output;
    bool verbose = false;
    bool quiet = false;
    
    // Special modes
    bool help = false;
    bool version = false;
    bool validate_config = false;
    bool dump_build_info = false;
    bool dump_hdf5_schema = false;
    bool list_collision_models = false;
    bool dry_run = false;
    
    // Config overrides
    std::optional<size_t> override_timesteps;
    std::optional<double> override_dt;
    std::optional<std::string> override_output_file;
};

class CliHandler {
public:
    // Parse CLI args -> CliOptions
    static CliOptions parse(int argc, char* argv[]);
    
    // Handle special modes (--help, --version, --validate-config, etc.)
    // Returns true if program should exit (special mode handled)
    static bool handle_special_modes(const CliOptions& opts);
    
private:
    static void print_help();
    static void print_version();
    static void dump_build_info();
    static void dump_hdf5_schema();
    static void list_collision_models();
};

} // namespace icarion::setup
```

**Key Methods:**
- `parse()`: Parse argc/argv → CliOptions struct
- `handle_special_modes()`: Execute --help, --version, --validate-config, etc. (returns true if early exit)
- Helper methods for each special mode (print_help, dump_build_info, etc.)

---

### 3.2 LoggingSetup (Logging Initialization)
**File:** `src/main/setup/LoggingSetup.h/cpp`  
**Purpose:** Initialize logging system based on CLI options

```cpp
// LoggingSetup.h
#pragma once
#include "CliHandler.h"

namespace icarion::setup {

class LoggingSetup {
public:
    // Initialize logging (sinks, verbosity) based on CLI options
    static void initialize(const CliOptions& opts);
    
    // Shutdown logging (call at program exit)
    static void shutdown();
};

} // namespace icarion::setup
```

**Key Methods:**
- `initialize()`: Setup console/file sinks, set log level (debug/info/warn based on verbose/quiet flags)
- `shutdown()`: Flush and close log files

---

### 3.3 ConfigSetup (Config Loading & Overrides)
**File:** `src/main/setup/ConfigSetup.h/cpp`  
**Purpose:** Load JSON config, apply CLI overrides, validate

```cpp
// ConfigSetup.h
#pragma once
#include "core/config/Config.h"
#include "CliHandler.h"

namespace icarion::setup {

class ConfigSetup {
public:
    // Load config from file, apply CLI overrides
    static config::Config load(const CliOptions& opts);
    
    // Validate config (schema validation, physics checks)
    // Throws std::runtime_error if invalid
    static void validate(const config::Config& cfg);
    
private:
    static void apply_overrides(config::Config& cfg, const CliOptions& opts);
};

} // namespace icarion::setup
```

**Key Methods:**
- `load()`: Load JSON config, apply CLI overrides (--timesteps, --dt, --output-file)
- `validate()`: Schema validation, physics checks (throws on error)
- `apply_overrides()`: Apply CLI overrides to config struct

---

### 3.4 IonSetup (Ion Generation)
**File:** `src/main/setup/IonSetup.h/cpp`  
**Purpose:** Generate initial ion states from config (ion_clouds or inline species)

```cpp
// IonSetup.h
#pragma once
#include "core/config/Config.h"
#include "core/integrator/IonState.h"
#include <vector>

namespace icarion::setup {

class IonSetup {
public:
    // Generate initial ion states from config
    // Handles both ion_clouds (HDF5 files) and inline species definitions
    static std::vector<IonState> generate_ions(const config::Config& cfg);
    
private:
    // Load ions from HDF5 ion cloud file
    static std::vector<IonState> load_ion_cloud(
        const std::string& filename,
        const std::optional<size_t>& max_ions
    );
    
    // Generate ions from inline species definition (uniform distribution)
    static std::vector<IonState> generate_from_species(
        const config::IonSpecies& species,
        const config::Config& cfg
    );
    
    // Validate ion positions vs domain bounds
    static void validate_ion_positions(
        const std::vector<IonState>& ions,
        const config::Config& cfg
    );
};

} // namespace icarion::setup
```

**Key Methods:**
- `generate_ions()`: Main entry point (dispatches to ion_clouds or inline species)
- `load_ion_cloud()`: Load from HDF5 file (uses io::IonCloudReader)
- `generate_from_species()`: Generate uniform distribution in domain
- `validate_ion_positions()`: Check all ions within domain bounds (warn if outside)

---

### 3.5 PhysicsSetup (Physics Module Creation)
**File:** `src/main/setup/PhysicsSetup.h/cpp`  
**Purpose:** Create force registries, collision handler, reaction handler

```cpp
// PhysicsSetup.h
#pragma once
#include "core/config/Config.h"
#include "core/physics/ForceRegistry.h"
#include "core/physics/ICollisionHandler.h"
#include "core/physics/IReactionHandler.h"
#include "core/integrator/IIntegrationStrategy.h"
#include <vector>
#include <memory>

namespace icarion::setup {

struct PhysicsModules {
    std::vector<std::shared_ptr<physics::ForceRegistry>> force_registries;
    std::shared_ptr<integrator::IIntegrationStrategy> integration_strategy;
    std::shared_ptr<physics::ICollisionHandler> collision_handler;
    std::shared_ptr<physics::IReactionHandler> reaction_handler;
};

class PhysicsSetup {
public:
    // Create all physics modules (forces, integration, collisions, reactions)
    static PhysicsModules create(
        const config::Config& cfg,
        const std::vector<IonState>& ions  // Needed for space charge auto-selection
    );
    
private:
    // Create ForceRegistry per domain (electric, magnetic, damping, space charge)
    static std::vector<std::shared_ptr<physics::ForceRegistry>> 
    create_force_registries(const config::Config& cfg, size_t n_ions);
    
    // Auto-select space charge method (Direct vs Grid based on N)
    static void add_space_charge_forces(
        std::vector<std::shared_ptr<physics::ForceRegistry>>& registries,
        const config::Config& cfg,
        const std::vector<IonState>& ions
    );
    
    // Create integration strategy (RK4, RK45, Boris)
    static std::shared_ptr<integrator::IIntegrationStrategy>
    create_integrator(const config::Config& cfg);
    
    // Create collision handler (HSS, EHSS, Friction, etc.)
    static std::shared_ptr<physics::ICollisionHandler>
    create_collision_handler(const config::Config& cfg);
    
    // Create reaction handler
    static std::shared_ptr<physics::IReactionHandler>
    create_reaction_handler(const config::Config& cfg);
    
    // Load molecular geometries for EHSS collision model
    static std::unique_ptr<physics::GeometryMap>
    load_geometry_map(const config::Config& cfg);
};

} // namespace icarion::setup
```

**Key Methods:**
- `create()`: Main entry point (creates all physics modules, returns struct)
- `create_force_registries()`: Create ForceRegistry per domain with electric, magnetic, damping forces
- `add_space_charge_forces()`: Auto-select Direct (N<1000) vs Grid (N≥1000) space charge
- `create_integrator()`: RK4/RK45/Boris strategy from config
- `create_collision_handler()`: HSS/EHSS/Friction handler (load geometries for EHSS)
- `create_reaction_handler()`: Enable/disable reactions

**Notes:**
- This is the **most complex helper** (~150 lines of logic)
- Auto-selection logic for space charge (Direct vs Grid) needs careful extraction
- EHSS geometry loading (fallback to HSS on error) needs error handling

---

### 3.6 EngineSetup (SimulationEngine Initialization)
**File:** `src/main/setup/EngineSetup.h/cpp`  
**Purpose:** Create and configure SimulationEngine

```cpp
// EngineSetup.h
#pragma once
#include "core/integrator/SimulationEngine.h"
#include "PhysicsSetup.h"
#include "core/config/Config.h"

namespace icarion::setup {

class EngineSetup {
public:
    // Create SimulationEngine with physics modules
    static integrator::SimulationEngine create(
        const config::Config& cfg,
        const PhysicsModules& physics
    );
};

} // namespace icarion::setup
```

**Key Methods:**
- `create()`: Instantiate SimulationEngine with config + physics modules

---

## 4. Refactored main.cpp (Target: <100 lines)

### After Refactoring
```cpp
// src/main/main.cpp
#include "setup/CliHandler.h"
#include "setup/LoggingSetup.h"
#include "setup/ConfigSetup.h"
#include "setup/IonSetup.h"
#include "setup/PhysicsSetup.h"
#include "setup/EngineSetup.h"
#include "utils/profiling/Profiler.h"
#include <chrono>

using namespace icarion;

int main(int argc, char* argv[]) {
    try {
        // 1. Parse CLI arguments
        auto opts = setup::CliHandler::parse(argc, argv);
        
        // 2. Handle special modes (--help, --version, --validate-config, etc.)
        if (setup::CliHandler::handle_special_modes(opts)) {
            return 0;  // Early exit for special modes
        }
        
        // 3. Initialize logging & profiling
        setup::LoggingSetup::initialize(opts);
        if (opts.benchmark || opts.profile) {
            profiling::Profiler::getInstance().enable();
        }
        
        // 4. Load & validate config
        auto config = setup::ConfigSetup::load(opts);
        setup::ConfigSetup::validate(config);
        
        log::Logger::main()->info("Config loaded: {}", opts.config_file);
        log::Logger::main()->info("Total time: {:.2e} s, dt: {:.2e} s",
                                  config.simulation.total_time_s,
                                  config.simulation.dt_s);
        
        // 5. Generate initial ions
        auto ions = setup::IonSetup::generate_ions(config);
        log::Logger::main()->info("Generated {} ions", ions.size());
        
        // 6. Create physics modules (forces, integration, collisions, reactions)
        auto physics = setup::PhysicsSetup::create(config, ions);
        log::Logger::main()->info("Physics modules initialized");
        
        // 7. Create simulation engine
        auto engine = setup::EngineSetup::create(config, physics);
        log::Logger::main()->info("SimulationEngine initialized");
        
        // 8. Run simulation
        log::Logger::main()->info("Starting simulation (t_max = {:.2e} s)",
                                  config.simulation.total_time_s);
        
        auto start = std::chrono::high_resolution_clock::now();
        auto final_ions = engine.run(ions);
        auto end = std::chrono::high_resolution_clock::now();
        
        double elapsed_s = std::chrono::duration<double>(end - start).count();
        
        // 9. Report results
        size_t active_count = std::count_if(
            final_ions.begin(), final_ions.end(),
            [](const auto& ion) { return ion.active; }
        );
        
        log::Logger::main()->info("");
        log::Logger::main()->info("=== Simulation Complete ===");
        log::Logger::main()->info("CPU time:     {:.3f} s", elapsed_s);
        log::Logger::main()->info("Active ions:  {}/{}", active_count, final_ions.size());
        log::Logger::main()->info("Output file:  {}", config.output.trajectory_file);
        log::Logger::main()->info("===========================");
        
        // 10. Export profiling data (if enabled)
        if ((opts.benchmark || opts.profile) && opts.profile_output.has_value()) {
            std::string filename = opts.profile_output.value();
            if (filename.find(".csv") != std::string::npos) {
                profiling::Profiler::getInstance().exportCSV(filename);
            } else {
                if (filename.find(".json") == std::string::npos) {
                    filename += ".json";
                }
                profiling::Profiler::getInstance().exportJSON(filename);
            }
            log::Logger::main()->info("Profile data written to: {}", filename);
        }
        
    } catch (const std::exception& e) {
        log::Logger::main()->error("Fatal error: {}", e.what());
        setup::LoggingSetup::shutdown();
        return 1;
    }
    
    setup::LoggingSetup::shutdown();
    return 0;
}
```

**Line Count:** ~85 lines (vs 557 original) ✅  
**Readability:** High-level orchestration, easy to follow  
**Testability:** Each helper independently testable

---

## 5. Implementation Strategy

### Phase 1: Create Setup Directory Structure
**Estimated Time:** 5 minutes

1. Create directory: `src/main/setup/`
2. Add to CMakeLists.txt:
   ```cmake
   add_subdirectory(setup)
   ```
3. Create `src/main/setup/CMakeLists.txt`:
   ```cmake
   add_library(icarion_setup STATIC
       CliHandler.cpp
       LoggingSetup.cpp
       ConfigSetup.cpp
       IonSetup.cpp
       PhysicsSetup.cpp
       EngineSetup.cpp
   )
   target_link_libraries(icarion_setup PUBLIC
       icarion_core
       icarion_io
       icarion_utils
   )
   ```

### Phase 2: Extract CliHandler (CLI Parsing)
**Estimated Time:** 30 minutes

1. Create `CliHandler.h/cpp`
2. Copy CLI parsing logic from main.cpp (lines 61-150)
3. Extract `CliOptions` struct
4. Implement `parse()` method
5. Implement `handle_special_modes()` and helper methods
6. Test: Run with --help, --version, --validate-config flags

### Phase 3: Extract LoggingSetup
**Estimated Time:** 15 minutes

1. Create `LoggingSetup.h/cpp`
2. Copy logging initialization from main.cpp (lines 151-180)
3. Implement `initialize()` and `shutdown()` methods
4. Test: Verify log files created correctly

### Phase 4: Extract ConfigSetup
**Estimated Time:** 20 minutes

1. Create `ConfigSetup.h/cpp`
2. Copy config loading logic (lines 231-280)
3. Implement `load()`, `validate()`, `apply_overrides()` methods
4. Test: Load sample config, apply overrides (--timesteps, --dt)

### Phase 5: Extract IonSetup
**Estimated Time:** 30 minutes

1. Create `IonSetup.h/cpp`
2. Copy ion generation logic (lines 281-340)
3. Implement `generate_ions()`, `load_ion_cloud()`, `generate_from_species()` methods
4. Test: Generate ions from ion_clouds and inline species

### Phase 6: Extract PhysicsSetup (Most Complex)
**Estimated Time:** 45 minutes

1. Create `PhysicsSetup.h/cpp`
2. Copy physics module creation (lines 341-480)
3. Implement `create()` method (main entry point)
4. Implement `create_force_registries()` (electric, magnetic, damping forces)
5. Implement `add_space_charge_forces()` (auto-selection logic)
6. Implement `create_integrator()`, `create_collision_handler()`, `create_reaction_handler()`
7. Implement `load_geometry_map()` (EHSS geometry loading)
8. Test: Verify force registries, collision handlers created correctly

### Phase 7: Extract EngineSetup
**Estimated Time:** 10 minutes

1. Create `EngineSetup.h/cpp`
2. Copy engine creation logic (lines 481-520)
3. Implement `create()` method
4. Test: Verify SimulationEngine initializes correctly

### Phase 8: Refactor main.cpp
**Estimated Time:** 20 minutes

1. Replace monolithic code with helper calls (see Section 4)
2. Remove extracted code (keep only orchestration)
3. Update includes
4. Test: Run full simulation, verify identical behavior

### Phase 9: Testing & Validation
**Estimated Time:** 20 minutes

1. Run CTest suite: `ctest --output-on-failure`
   - **Success Criteria:** 48/48 tests pass ✅
2. Run sample simulations: `examples/ims_basic.json`
   - **Success Criteria:** Identical output to baseline
3. Test special modes: `--help`, `--version`, `--validate-config`, `--dry-run`
   - **Success Criteria:** All modes work correctly
4. Test CLI overrides: `--timesteps 2000 --dt 1e-9`
   - **Success Criteria:** Overrides applied correctly

---

## 6. Success Criteria

### ✅ Functional Requirements
- [ ] 48/48 CTest suite passes (no regressions)
- [ ] Sample simulations produce identical output (bit-exact results)
- [ ] All CLI flags work correctly (--help, --version, --validate-config, etc.)
- [ ] CLI overrides apply correctly (--timesteps, --dt, --output-file)
- [ ] Profiling/benchmarking work correctly (--profile, --benchmark)

### ✅ Code Quality Requirements
- [ ] main.cpp reduced to <100 lines (target: ~85 lines)
- [ ] Each setup helper <200 lines (single responsibility)
- [ ] No code duplication (DRY principle)
- [ ] Clear separation of concerns (CLI, config, ions, physics, engine)

### ✅ Documentation Requirements
- [ ] Each helper class documented (purpose, methods, usage)
- [ ] Update DEVELOPERS_GUIDE.md (new setup/ directory structure)
- [ ] Add unit test examples for setup helpers

---

## 7. Risk Mitigation

### Potential Issues
1. **Dependency Ordering:** Setup helpers call each other (e.g., PhysicsSetup needs IonSetup)
   - **Mitigation:** Define clear interfaces (PhysicsModules struct passed between helpers)

2. **Error Handling:** Current code uses exceptions, helpers must preserve this
   - **Mitigation:** Use `try-catch` in main.cpp, helpers throw std::runtime_error

3. **Special Modes Complexity:** --validate-config, --dry-run interact with config loading
   - **Mitigation:** handle_special_modes() after config loading (where needed)

4. **EHSS Geometry Loading:** Can fail, current code falls back to HSS
   - **Mitigation:** PhysicsSetup::load_geometry_map() handles fallback logic

5. **Space Charge Auto-Selection:** Depends on ion count, needs access to ions vector
   - **Mitigation:** PhysicsSetup::create() takes ions vector as parameter

---

## 8. Testing Checklist

### Unit Tests (Optional but Recommended)
- [ ] CliHandler::parse() with valid/invalid args
- [ ] ConfigSetup::apply_overrides() with CLI overrides
- [ ] IonSetup::generate_ions() with ion_clouds and inline species
- [ ] PhysicsSetup::create_force_registries() with different domain configs
- [ ] PhysicsSetup::add_space_charge_forces() with N<1000 and N≥1000

### Integration Tests (CTest Suite)
- [ ] Run full CTest suite: `ctest --output-on-failure`
- [ ] Verify 48/48 tests pass ✅

### Manual Tests
- [ ] Run `icarion --help` (verify help text)
- [ ] Run `icarion --version` (verify version info)
- [ ] Run `icarion --validate-config examples/ims_basic.json` (verify validation)
- [ ] Run `icarion --dry-run examples/ims_basic.json` (verify dry-run)
- [ ] Run `icarion examples/ims_basic.json` (full simulation)
- [ ] Run `icarion examples/ims_basic.json --timesteps 5000 --dt 1e-9` (CLI overrides)
- [ ] Run `icarion examples/ims_basic.json --profile --profile-output profile.json` (profiling)

---

## 9. Next Steps After Refactoring

### Follow-up Tasks
1. **Unit Tests for Setup Helpers:** Write tests for ion generation, physics setup logic
2. **ConfigLoader/OutputManager Review:** Evaluate if these need refactoring (lower priority)
3. **GUI Integration:** Setup helpers reusable for future GUI application
4. **Performance Optimization:** After refactoring, ready for CPU optimizations
5. **GPU Implementation:** After refactoring, ready for modular GPU components

---

## 10. Notes & Assumptions

### Assumptions
- Current code is correct (557-line main.cpp works, but hard to maintain)
- CTest suite validates correctness (48/48 tests must pass after refactoring)
- Pure refactoring (no logic changes, only code reorganization)

### Design Decisions
- **Why separate helpers?** Single responsibility, testability, reusability
- **Why inline in headers?** Some helpers are simple (1-2 methods), can be header-only
- **Why PhysicsModules struct?** Clean interface for passing physics components to EngineSetup
- **Why not use singletons?** Setup helpers are stateless, static methods are cleaner

### Edge Cases
- **Empty ion_clouds:** IonSetup::generate_ions() should throw error
- **Invalid integrator name:** PhysicsSetup::create_integrator() defaults to RK45 (with warning)
- **EHSS geometry load failure:** PhysicsSetup::create_collision_handler() falls back to HSS
- **N=0 ions:** Space charge auto-selection should handle gracefully (no space charge forces)

---

## 11. Conclusion

**Status:** 🟢 READY TO IMPLEMENT (low risk, ~2-3 hours)

This refactoring will:
- ✅ Reduce main.cpp from 557 → ~85 lines (84% reduction)
- ✅ Improve maintainability (each setup phase in separate file)
- ✅ Enable unit testing for setup logic
- ✅ Prepare for GUI integration (setup helpers reusable)
- ✅ No performance impact (pure refactoring, no logic changes)

**Recommendation:** Proceed with implementation after approval.

---

**Author:** GitHub Copilot (Claude Sonnet 4.5)  
**Review Required:** Yes (user approval before implementation)
