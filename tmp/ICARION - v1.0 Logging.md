# **4) High-Quality Console Output**

Einfach, klar, nicht überladen — im Stil moderner HPC Codes:
```yaml
=========================================================
 ICARION v1.0 — Simulation Started
---------------------------------------------------------
 Title: Orbitrap Example
 Instrument: Orbitrap
 GPU: Enabled
 Timestep: 1.00e-10 s
 Total time: 1.00e-4 s (1,000,000 steps)
 Stochastic collisions: Enabled (EHSS)
 Reactions: Enabled
 Output file: orbitrap.h5
=========================================================
```
Zwischenstände:
```yaml
Step 10000 / 1000000 (1.0%)  - avg dt: 1.2 µs   ions active: 998
Step 20000 / ...             - collisions: 1.2M  reactions: 200
```
Abschluss:
```yaml
Simulation completed in 12.4 seconds
HDF5 file written successfully: orbitrap.h5
```

# ✅ **5) High-Quality Logging System**

Du bietest:

### **logging/debug.log**

- RK4 step statistics
    
- Collision counters
    
- Energy drift
    
- CPU/GPU divergence (falls aktiviert)
    
- Warning system
    

### **logging/errors.log**

- invalid configs
    
- missing species
    
- divergence / instability
    
- reaction errors
    

### **logging/performance.log**

- GPU timing
    
- Kernel occupancy
    
- Integrator throughput
    
- Memory usage
    

Der Nutzer kann Logging-Level einstellen:
```pgsql
"logging": {
  "level": "info",  // info | debug | warning | error
  "collision_debug": false,
  "write_performance": true
}
```