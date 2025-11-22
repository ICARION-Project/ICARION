# 🚧 VERBLEIBENDE SSOT MIGRATION ARBEITEN

**Branch:** `refactor/force-system-ssot`  
**Stand:** 2025-11-22 nach Step 4 (Tests komplett)  
**Status:** ✅ Force-Klassen fertig | ⏳ Integrator & main.cpp ausstehend

---

## ✅ BEREITS ERLEDIGT (Steps 1-4)

### **Step 1: MagneticFieldForce** ✅ (Commit 9e89b36)
- Constructor: `MagneticFieldForce(const config::MagneticFieldConfig&)`
- Speichert `const config::MagneticFieldConfig& magnetic_`
- Liest `magnetic_.field_strength_T` (Vec3)
- Deleted: `MagneticFieldParams` struct (~80 Zeilen)

### **Step 2: ElectricFieldForce** ✅ (Commit 90b02aa)
- Constructor: `ElectricFieldForce(const config::DomainConfig&)`
- Speichert `const config::DomainConfig* domain_`
- Liest `domain_->fields.dc.axial_V`, `domain_->fields.rf.*`, etc.
- Deleted: `AnalyticalFieldParams` struct (~100 Zeilen)
- Fix: Added `NoFixedInstrument` support

### **Step 3: DampingForce** ✅ (Commit 23327fd)
- Constructor: `DampingForce(const config::EnvironmentConfig&, DampingModel)`
- Speichert `const config::EnvironmentConfig& env_`
- Liest `env_.pressure_Pa`, `env_.temperature_K`, `env_.particle_density_m_3`
- Deleted: `DampingParams` struct (~40 Zeilen)

### **Step 4: Test-Updates** ✅ (Commits 7d555c4, fda11d2, 89dfc03)
- Alle 6 Force-Test-Dateien aktualisiert
- Test-Erwartungen korrigiert (TOF, AC field, Friction, Orbitrap)
- Alle Force-Tests bestehen (100%)
- **Gesamt gelöscht:** ~200 Zeilen Parameter-Duplikation

---

## ⏳ NOCH ZU ERLEDIGEN (Steps 5-10)

### **Step 5: Update compute_accelerations()** ⚠️ KRITISCH

**Dateien:**
- `src/core/physics/computeAccelerations.h`
- `src/core/physics/computeAccelerations.cpp`

**Änderungen:**

1. **Signatur ändern:**
   ```cpp
   // ALT (Legacy):
   IonState compute_accelerations(
       const IonState& ion,
       double t,
       const GlobalParams& gp,
       const InstrumentDomain& domain,
       const std::vector<IonState>& all_ions,
       bool enable_space_charge,
       const std::string& solver_type,
       IFieldProvider* field_provider = nullptr
   );
   
   // NEU (SSOT):
   IonState compute_accelerations(
       const IonState& ion,
       double t,
       const config::DomainConfig& domain,
       const std::vector<IonState>& all_ions,
       bool enable_space_charge,
       const std::string& solver_type,
       IFieldProvider* field_provider = nullptr
   );
   ```

2. **Feld-Berechnungen aktualisieren:**
   ```cpp
   // ALT:
   Vec3 E = ElectricFieldForce::compute_field(ion, t, gp, domain);
   Vec3 B = domain.uniform_field_T;  // von InstrumentDomain
   
   // NEU:
   ElectricFieldForce ef(domain);  // nimmt DomainConfig
   Vec3 E = ef.compute(ion, t, ctx);
   
   MagneticFieldForce mf(domain.fields.magnetic);
   Vec3 B = mf.compute(ion, t, ctx);
   ```

3. **Damping aktualisieren:**
   ```cpp
   // ALT:
   DampingParams dp;
   dp.pressure_Pa = gp.pressure_Pa;
   dp.temperature_K = gp.temperature_K;
   DampingForce df(dp, model);
   
   // NEU:
   DampingForce df(domain.environment, model);
   ```

4. **Space Charge aktualisieren:**
   - Keine Änderung nötig (verwendet bereits `all_ions` direkt)

**Status:** ❌ Nicht begonnen  
**Aufwand:** 45 min  
**Blocker:** Keine

---

### **Step 6: Update integrator_helpers.cpp** ⚠️ KRITISCH

**Dateien:**
- `src/core/integrator/integrator_helpers.h`
- `src/core/integrator/integrator_helpers.cpp`

**Änderungen:**

1. **`integrate_one_step()` Signatur:**
   ```cpp
   // ALT:
   IonState integrate_one_step(
       const IonState& ion_state,
       const GlobalParams& gp,
       const InstrumentDomain& domain,
       ...
   );
   
   // NEU:
   IonState integrate_one_step(
       const IonState& ion_state,
       const config::DomainConfig& domain,
       ...
   );
   ```

2. **`integrate_trajectory()` Signatur:**
   ```cpp
   // ALT:
   std::vector<IonState> integrate_trajectory(
       const IonState& ion_init,
       const GlobalParams& gp,
       const std::vector<InstrumentDomain>& domains,
       ...
   );
   
   // NEU:
   std::vector<IonState> integrate_trajectory(
       const IonState& ion_init,
       const config::FullConfig& config,  // enthält domains
       ...
   );
   ```

3. **Multi-Domain Handling:**
   ```cpp
   // NEU:
   for (size_t domain_idx = 0; domain_idx < config.domains.size(); ++domain_idx) {
       const auto& domain = config.domains[domain_idx];
       // Integration mit domain
   }
   ```

4. **Collision Handler Factory:**
   ```cpp
   // ALT:
   auto handler = CollisionHandlerFactory::create(gp.collision_model, ...);
   
   // NEU:
   auto handler = CollisionHandlerFactory::create(
       domain.environment.collision_model,
       domain.environment
   );
   ```

**Status:** ❌ Nicht begonnen  
**Aufwand:** 60 min  
**Blocker:** Wartet auf Step 5

---

### **Step 7: Update integrator.cpp** ⚠️ KRITISCH

**Dateien:**
- `src/core/integrator/integrator.cpp`

**Änderungen:**

1. **RK4 Integrator:**
   ```cpp
   // ALT:
   k1 = compute_accelerations(ion, t, gp, domain, ...);
   k2 = compute_accelerations(ion_temp, t + dt/2, gp, domain, ...);
   
   // NEU:
   k1 = compute_accelerations(ion, t, domain, ...);
   k2 = compute_accelerations(ion_temp, t + dt/2, domain, ...);
   ```

2. **DOPRI5 (RK45) Integrator:**
   - Gleiche Änderungen wie RK4
   - Alle 7 k-Stufen aktualisieren

3. **Boris Integrator:**
   ```cpp
   // ALT:
   Vec3 B = domain.uniform_field_T;
   
   // NEU:
   MagneticFieldForce mf(domain.fields.magnetic);
   Vec3 B = mf.compute(ion, t, ctx);
   ```

**Status:** ❌ Nicht begonnen  
**Aufwand:** 30 min  
**Blocker:** Wartet auf Step 5

---

### **Step 8: Update main.cpp** ⚠️ KRITISCH

**Dateien:**
- `src/main.cpp`

**Änderungen:**

1. **LegacyAdapter entfernen:**
   ```cpp
   // ALT (zu löschen):
   #include "core/config/adapter/LegacyAdapter.h"
   
   auto gp = LegacyAdapter::to_global_params(full_config);
   auto domains = LegacyAdapter::to_instrument_domains(full_config);
   
   // Simulation
   auto results = integrate_trajectory(ion, gp, domains, ...);
   
   // NEU:
   auto results = integrate_trajectory(ion, full_config, ...);
   ```

2. **Direkter FullConfig-Zugriff:**
   ```cpp
   // Kein Adapter mehr - direkte Verwendung
   const auto& domain = full_config.domains[0];
   
   // Force-Erstellung
   ElectricFieldForce ef(domain);
   MagneticFieldForce mf(domain.fields.magnetic);
   DampingForce df(domain.environment, model);
   ```

3. **Output-Pfad (falls nötig):**
   ```cpp
   // ALT:
   std::string output_path = gp.output_path;
   
   // NEU:
   std::string output_path = full_config.simulation.output_path;
   ```

**Status:** ❌ Nicht begonnen  
**Aufwand:** 30 min  
**Blocker:** Wartet auf Steps 5, 6, 7

---

### **Step 9: Legacy-Code löschen** 🗑️

**Zu löschende Dateien:**

1. **LegacyAdapter:**
   - `src/core/config/adapter/LegacyAdapter.h`
   - `src/core/config/adapter/LegacyAdapter.cpp`
   - `tests/config/test_legacy_adapter.cpp`

2. **GlobalParams/InstrumentDomain:**
   - `src/core/param/paramUtils.h`
   - `src/core/param/paramUtils.cpp`
   - Alle Includes von `paramUtils.h` entfernen

3. **Parameter-Structs (bereits gelöscht):**
   - ✅ `MagneticFieldParams` (gelöscht in Step 1)
   - ✅ `AnalyticalFieldParams` (gelöscht in Step 2)
   - ✅ `DampingParams` (gelöscht in Step 3)

**Verifikation:**
```bash
# Sicherstellen, dass keine Legacy-Referenzen bleiben:
grep -r "GlobalParams" src/
grep -r "InstrumentDomain" src/
grep -r "LegacyAdapter" src/
grep -r "MagneticFieldParams" src/
grep -r "AnalyticalFieldParams" src/
grep -r "DampingParams" src/
```

**Status:** ❌ Nicht begonnen  
**Aufwand:** 15 min  
**Blocker:** Wartet auf Steps 5-8

---

### **Step 10: Dokumentation** 📝

**Dateien zu aktualisieren:**

1. **`docs/ARCHITECTURE.md`:**
   - Section "Force System" aktualisieren
   - SSOT-Prinzip dokumentieren
   - Neue Constructor-Signaturen zeigen

2. **`docs/DEVELOPERS_GUIDE.md`:**
   - "Adding New Forces" aktualisieren
   - Config-Zugriff dokumentieren
   - Alte Beispiele entfernen

3. **`docs/PUBLIC_CPP_API_v1.0.md`:**
   - API-Änderungen dokumentieren
   - Breaking Changes auflisten

4. **`RELEASE_NOTES_v1.1.md` (neu):**
   - Migration Guide schreiben
   - Breaking Changes auflisten
   - Code-Beispiele (alt vs. neu)

**Status:** ❌ Nicht begonnen  
**Aufwand:** 30 min  
**Blocker:** Wartet auf Steps 5-9

---

## 📊 FORTSCHRITT

| Step | Beschreibung | Status | Zeit |
|------|--------------|--------|------|
| 1 | MagneticFieldForce | ✅ Fertig | 30 min |
| 2 | ElectricFieldForce | ✅ Fertig | 45 min |
| 3 | DampingForce | ✅ Fertig | 30 min |
| 4 | Tests | ✅ Fertig | 90 min |
| **5** | **compute_accelerations()** | ❌ Offen | **45 min** |
| **6** | **integrator_helpers** | ❌ Offen | **60 min** |
| **7** | **integrator.cpp** | ❌ Offen | **30 min** |
| **8** | **main.cpp** | ❌ Offen | **30 min** |
| **9** | **Legacy löschen** | ❌ Offen | **15 min** |
| **10** | **Dokumentation** | ❌ Offen | **30 min** |
| | **GESAMT** | **40% fertig** | **6h 15min** |
| | **Verbleibend** | | **3h 30min** |

---

## 🎯 NÄCHSTE SCHRITTE (Reihenfolge wichtig!)

### **1. Step 5: compute_accelerations() (45 min)** 🔴 HÖCHSTE PRIORITÄT
- Zentrale Funktion, blockiert alle anderen
- Signatur ändern: `GlobalParams` → `DomainConfig`
- Force-Konstruktoren aktualisieren
- Tests schreiben/aktualisieren

### **2. Step 6: integrator_helpers (60 min)**
- `integrate_one_step()` aktualisieren
- `integrate_trajectory()` aktualisieren
- Multi-Domain Handling modernisieren

### **3. Step 7: integrator.cpp (30 min)**
- RK4, DOPRI5, Boris aktualisieren
- Alle `compute_accelerations()` Aufrufe ändern

### **4. Step 8: main.cpp (30 min)**
- LegacyAdapter entfernen
- Direkter FullConfig-Zugriff
- Smoke-Test durchführen

### **5. Step 9: Legacy löschen (15 min)**
- LegacyAdapter löschen
- paramUtils löschen
- Grep-Checks durchführen

### **6. Step 10: Dokumentation (30 min)**
- Architecture Guide aktualisieren
- Developer Guide aktualisieren
- Release Notes schreiben

---

## ⚠️ KRITISCHE ABHÄNGIGKEITEN

```
Step 5 (compute_accelerations)
    ↓
Step 6 (integrator_helpers)
    ↓
Step 7 (integrator.cpp)
    ↓
Step 8 (main.cpp)
    ↓
Step 9 (Legacy löschen)
    ↓
Step 10 (Dokumentation)
```

**WICHTIG:** Steps 5-8 müssen sequentiell erfolgen!

---

## ✅ ERFOLGS-KRITERIEN

### **Code:**
- [ ] Keine `GlobalParams` Referenzen in `src/`
- [ ] Keine `InstrumentDomain` Referenzen in `src/`
- [ ] Keine `LegacyAdapter` Dateien
- [ ] Alle Forces verwenden `const DomainConfig&` oder Sub-Configs
- [ ] Kein Parameter-Duplikation (SSOT!)

### **Tests:**
- [ ] Alle Unit-Tests bestehen (24/24)
- [ ] Integration-Test mit vollem Workflow
- [ ] Performance-Regression < 5%

### **Kompilierung:**
- [ ] Keine Compiler-Warnungen
- [ ] Keine Link-Fehler
- [ ] Clean build von scratch

### **Dokumentation:**
- [ ] Architecture Guide aktualisiert
- [ ] Developer Guide aktualisiert
- [ ] API-Dokumentation aktualisiert
- [ ] Release Notes geschrieben

---

## 🚀 GESCHÄTZTER ZEITAUFWAND

**Verbleibend:** 3h 30min  
**Optimal:** 1 Session (halber Tag)  
**Realistisch:** 2 Sessions (über 2 Tage)

---

## 📝 COMMIT-STRATEGIE

```bash
# Step 5
git commit -m "refactor(physics): Update compute_accelerations to use DomainConfig (Step 5)"

# Step 6
git commit -m "refactor(integrator): Update integrator_helpers to use FullConfig (Step 6)"

# Step 7
git commit -m "refactor(integrator): Update integrator.cpp to use DomainConfig (Step 7)"

# Step 8
git commit -m "refactor(main): Remove LegacyAdapter, use FullConfig directly (Step 8)"

# Step 9
git commit -m "refactor: Delete legacy code (LegacyAdapter, paramUtils) (Step 9)"

# Step 10
git commit -m "docs: Update documentation for SSOT migration (Step 10)"
```

---

**Letzte Aktualisierung:** 2025-11-22, 16:45  
**Nächster Schritt:** Step 5 - compute_accelerations() aktualisieren
