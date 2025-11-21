# HDF5 Output System - Refactoring Plan

**Status:** HDF5Writer v2 implementiert und getestet ✅  
**Datum:** 2025-11-21  
**Branch:** `io-refactor`

---

## 🎯 Strategie: Parallele Systeme bis Integrator-Refactor

### **Entscheidung:**
Legacy HDF5Writer **parallel laufen lassen** statt Hybrid-Ansatz oder Reverse-Adapter.

**Begründung:**
- ✅ Saubere Trennung (keine Vermischung von Alt/Neu)
- ✅ Legacy bleibt stable während Entwicklung
- ✅ Weniger Fehleranfälligkeit
- ✅ Klarer Migrationspfad

---

## 📋 3-Phasen Plan

### **Phase 1: HDF5Writer v2 fertigstellen** ✅ ERLEDIGT

**Implementiert:**
```
src/core/io/
├── hdf5Writer_v2.h         ✅ Modern API (3 public methods)
├── hdf5Writer_v2.cpp       ✅ Complete implementation (772 lines)
└── [Legacy bleibt]
    ├── hdf5Writer.h        ⚠️  Deprecated (für alte Tests)
    └── hdf5Writer.cpp      ⚠️  Deprecated (vom Integrator genutzt)
```

**Features:**
- ✅ `create_file()` - Schreibt komplette Metadaten (config, reproducibility, system, species, domains, ions)
- ✅ `append_trajectory()` - 3D Arrays [T × N × 3] mit Chunking + GZIP
- ✅ `finalize()` - Completion metadata mit Timestamp
- ✅ 10 Tests mit 48 Assertions (alle bestehen)
- ✅ Dokumentation: `docs/HDF5_OUTPUT_STRUCTURE.md` (447 lines)

**Status:** Production-ready! Wartet auf Integrator-Refactor.

---

### **Phase 2: Integrator refactorieren** ⏳ TODO (später)

**Ziel:** Clean SimulationEngine mit FullConfig

#### **Neue Architektur:**

```cpp
// src/core/integrator/SimulationEngine.h (NEW)

namespace ICARION::sim {

class SimulationEngine {
public:
    /**
     * @brief Construct simulation engine
     * @param config Complete simulation configuration (replaces GlobalParams)
     */
    SimulationEngine(const config::FullConfig& config);
    
    /**
     * @brief Run complete simulation with HDF5 output
     * @param ions Initial ion ensemble
     * @return Final ion states
     */
    std::vector<IonState> run(std::vector<IonState>& ions);
    
    /**
     * @brief Run single timestep (for custom control)
     */
    void step(std::vector<IonState>& ions, double t, double dt);
    
private:
    const config::FullConfig& config_;
    std::string hdf5_file_;
    
    void setup_output(const std::vector<IonState>& ions);
    void write_snapshot(double t, const std::vector<IonState>& ions);
    void finalize(const std::vector<IonState>& ions);
};

} // namespace ICARION::sim
```

#### **Implementation Highlights:**

```cpp
// src/core/integrator/SimulationEngine.cpp

#include "core/io/hdf5Writer_v2.h"  // ✅ Modern writer only

void SimulationEngine::setup_output(const std::vector<IonState>& ions) {
    // ✅ Clean HDF5 v2 API
    io::HDF5Writer::create_file(
        hdf5_file_,
        config_,      // FullConfig (not GlobalParams!)
        ions,
        GIT_HASH,     // From CMake
        __VERSION__   // Compiler info
    );
}

void SimulationEngine::write_snapshot(double t, const std::vector<IonState>& ions) {
    // ✅ No buffer complexity - direct append
    io::HDF5Writer::append_trajectory(hdf5_file_, t, ions);
}

void SimulationEngine::finalize(const std::vector<IonState>& ions) {
    size_t active = std::count_if(ions.begin(), ions.end(), 
                                   [](auto& ion) { return ion.active; });
    
    // ✅ Completion metadata
    io::HDF5Writer::finalize(hdf5_file_, true, config_.simulation.total_time_s, active);
}
```

#### **Vorteile gegenüber Legacy:**

| **Aspekt** | **Legacy (integrator.cpp)** | **Modern (SimulationEngine)** |
|------------|----------------------------|-------------------------------|
| Config | GlobalParams (flat struct) | FullConfig (hierarchisch) |
| HDF5 Init | `write_params_to_HDF5()` | `HDF5Writer::create_file()` |
| Trajectory | Buffer + `append_to_HDF5()` | Direct `append_trajectory()` |
| Finalize | `write_simulation_metadata()` | `HDF5Writer::finalize()` |
| Metadata | Minimal | Vollständig (reproducibility, system, species, domains) |
| Format | Legacy (undokumentiert) | HDF5 v2.0 (pandas-kompatibel) |

---

### **Phase 3: Legacy löschen** ⏳ TODO (nach Phase 2)

**Cleanup Steps:**

```bash
# 1. Delete legacy HDF5 writer
rm src/core/io/hdf5Writer.h
rm src/core/io/hdf5Writer.cpp

# 2. Delete legacy integrator
rm src/core/integrator/integrator.h
rm src/core/integrator/integrator.cpp

# 3. Rename v2 to primary
mv src/core/io/hdf5Writer_v2.h src/core/io/hdf5Writer.h
mv src/core/io/hdf5Writer_v2.cpp src/core/io/hdf5Writer.cpp

# 4. Update includes
find src -name "*.cpp" -o -name "*.h" | xargs sed -i 's/hdf5Writer_v2.h/hdf5Writer.h/g'

# 5. Update CMakeLists.txt
# Remove: src/core/integrator/integrator.cpp
```

**Ergebnis:**

```
src/core/io/
├── hdf5Writer.h         ✅ Modern API (renamed from v2)
└── hdf5Writer.cpp       ✅ Modern implementation

src/core/integrator/
├── SimulationEngine.h   ✅ Clean API with FullConfig
└── SimulationEngine.cpp ✅ Uses modern HDF5Writer
```

---

## 📊 Vergleich: Alternativen (NICHT gewählt)

### **❌ Alternative A: Reverse-Adapter**
```cpp
FullConfig FullConfigBuilder::from_global_params(GlobalParams& gParams);
```
**Problem:** Informationsverlust, doppelte Konversion (FullConfig → GlobalParams → FullConfig)

### **❌ Alternative B: Hybrid-Ansatz**
```cpp
void integrate_trajectory(..., const FullConfig* full_config = nullptr) {
    if (full_config != nullptr) {
        HDF5Writer::create_file(...);  // Modern
    } else {
        write_params_to_HDF5(...);     // Legacy
    }
}
```
**Problem:** Zwei Code-Pfade in einer Funktion, schwer zu testen

### **✅ Alternative C: Parallel (GEWÄHLT)**
- Legacy bleibt komplett unverändert
- Modern system komplett getrennt
- Cleanup nach Integrator-Refactor

---

## 🚀 Timeline

| **Phase** | **Task** | **Status** | **Zeit** |
|-----------|----------|------------|---------|
| **Phase 1** | HDF5Writer v2 Implementation | ✅ DONE | ~5h |
| | - API Design & Header | ✅ | 1h |
| | - Metadata Writers | ✅ | 2h |
| | - Trajectory Append | ✅ | 1h |
| | - Tests (10 cases) | ✅ | 1h |
| **Phase 2** | Integrator Refactoring | ⏳ TODO | ~2 weeks |
| | - SimulationEngine Design | ⏳ | 2h |
| | - Physics Integration | ⏳ | 3 days |
| | - Domain Management | ⏳ | 2 days |
| | - Testing & Validation | ⏳ | 1 week |
| **Phase 3** | Legacy Cleanup | ⏳ TODO | ~1h |
| | - Delete old files | ⏳ | 10 min |
| | - Rename v2 → primary | ⏳ | 10 min |
| | - Update includes | ⏳ | 20 min |
| | - Test suite update | ⏳ | 20 min |

---

## 📦 Commits (Phase 1)

```
800a339  feat(io): Add modern HDF5Writer v2 for FullConfig
90d9111  feat(io): Fix HDF5Writer v2 field name alignment
6e4613c  test(io): Add comprehensive HDF5Writer v2 integration tests
d3db845  feat(io): Document HDF5Writer v2 integration path in main.cpp
```

**Files Changed:**
- `src/core/io/hdf5Writer_v2.h` (169 lines)
- `src/core/io/hdf5Writer_v2.cpp` (772 lines)
- `tests/io/test_hdf5_writer_v2.cpp` (460 lines)
- `docs/HDF5_OUTPUT_STRUCTURE.md` (447 lines)
- `src/main/main.cpp` (integration docs)

**Total:** ~1,850 lines of new code

---

## 🔍 HDF5 v2.0 Format Specification

### **File Structure:**
```
simulation.h5
├── /metadata/
│   ├── /config/              # Simulation parameters (dt, integrator, collision_model)
│   ├── /reproducibility/     # Git hash, RNG seed, compiler, build info
│   ├── /system/              # Hostname, CPU, GPU, memory, timestamp
│   ├── /species/             # Tabular: names[], mass_kg[], charge_C[], mobility[], CCS[]
│   ├── /reactions/           # Tabular: reactants[], products[], rate_constants[]
│   └── /completion/          # Success flag, final_time, active_ions, timestamp
├── /trajectory/
│   ├── time [T]              # Time array [s]
│   ├── positions [T×N×3]     # Ion positions [m]
│   ├── velocities [T×N×3]    # Ion velocities [m/s]
│   ├── domain_indices [T×N]  # Domain location per ion
│   └── species_ids [T×N]     # Species ID per ion (variable-length strings)
├── /ions/
│   ├── initial_species_id [N]  # Initial species per ion
│   ├── initial_pos_x/y/z [N]   # Initial positions [m]
│   ├── initial_vel_x/y/z [N]   # Initial velocities [m/s]
│   ├── birth_time_s [N]        # Creation time [s]
│   └── charge_C [N]            # Charge [C]
└── /domains/
    ├── /domain_0/
    │   ├── name, instrument, solver
    │   ├── /geometry/        # length_m, radius_m, origin_m
    │   ├── /environment/     # pressure_Pa, temperature_K, gas_species
    │   └── /fields/          # DC, RF, AC field configurations
    └── /domain_1/...
```

### **Key Features:**
- ✅ **Pandas-compatible:** Tabular species/reactions data
- ✅ **Compressed:** GZIP level 6, chunked storage [100, N, 3]
- ✅ **Extensible:** H5S_UNLIMITED dimension for trajectory
- ✅ **Reproducible:** Git hash, RNG seed, compiler version
- ✅ **Complete:** All metadata for simulation reconstruction

---

## 📚 Documentation

### **Files:**
1. `docs/HDF5_OUTPUT_STRUCTURE.md` - Complete format specification
2. `tests/io/test_hdf5_writer_v2.cpp` - API usage examples
3. `src/main/main.cpp` - Integration path documentation

### **Python Analysis Example:**

```python
import h5py
import pandas as pd

# Load HDF5 file
with h5py.File('simulation.h5', 'r') as f:
    # Read species metadata as DataFrame
    species_df = pd.DataFrame({
        'name': f['/metadata/species/names'][:].astype(str),
        'mass_kg': f['/metadata/species/mass_kg'][:],
        'charge_C': f['/metadata/species/charge_C'][:],
        'mobility_m2Vs': f['/metadata/species/mobility_m2Vs'][:],
        'CCS_m2': f['/metadata/species/CCS_m2'][:]
    })
    
    # Read trajectory
    time = f['/trajectory/time'][:]
    positions = f['/trajectory/positions'][:]  # [T × N × 3]
    
    # Get config
    dt = f['/metadata/config/dt_s'][()]
    integrator = f['/metadata/config/integrator'][()].decode()
    
    print(f"Timestep: {dt} s")
    print(f"Integrator: {integrator}")
    print(f"Timesteps: {len(time)}")
    print(f"Ions: {positions.shape[1]}")
    print("\nSpecies:")
    print(species_df)
```

---

## ✅ Acceptance Criteria (Phase 1)

- [x] HDF5Writer v2 API implementiert
- [x] Alle Metadaten-Gruppen geschrieben
- [x] Trajectory append mit Compression
- [x] Finalize mit Timestamp
- [x] 10 Tests (48 Assertions) bestehen
- [x] Dokumentation vollständig
- [x] Build erfolgreich
- [x] Integration-Doku in main.cpp

**Status:** ✅ PRODUCTION READY

---

## 🔄 Next Steps

1. **Jetzt:** Branch `io-refactor` mergen nach `main` (oder `develop`)
2. **Später:** Integrator-Refactoring starten (Phase 2)
3. **Cleanup:** Legacy HDF5Writer löschen (Phase 3)

---

## 📞 Contact

**Fragen zu HDF5Writer v2?**
- Siehe: `docs/HDF5_OUTPUT_STRUCTURE.md`
- Tests: `tests/io/test_hdf5_writer_v2.cpp`
- API: `src/core/io/hdf5Writer_v2.h`
