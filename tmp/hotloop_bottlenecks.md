## Hotloop Bottlenecks (SimulationEngine & Hybrid)

### Main Findings
- **Doppelte Domain-Suche:** `find_domain_index` vor Integration und erneut für Boundary → zwei Geometrie-Lookups pro Ion/Step.
- **Viele kleine OpenMP-Loops:** Birth/Domain-Find/Integration/Boundary/Time-Update in separaten Schleifen → hoher OMP-Overhead/Synchronisation.
- **Memory-Bound SoA:** Pro Ion viele Arrays (pos/vel/active/born/domain/temp/density/mass/charge), keine Tiling/Blocking → Bandbreite limitiert, OMP skaliert flach.
- **Integrator Stage-Overhead:** RK4/RK45 bauen pro Stage Scratch-Ensembles (`resize(1)` + Kopie) und setzen `ForceContext` 6–7x neu → viel Kleinarbeit, wenig Rechenarbeit.
- **Collisions/Reactions:** Bei fehlendem Batch per-Ion Aufrufe mit hohem Funktions-/RNG-Overhead; `per_domain` Vektoren werden pro Step neu allokiert/gefüllt.
- **Space-Charge-Update pro Step:** Teure Modelle (Direct O(N²), Grid-Update) laufen jedes Mal, kein Throttle.
- **I/O/Locking:** Output/HDF5-Locks können blockieren (s. Benchmark); Progress/Safety-Logging erhöht Overhead.
- **Hybrid -30%:** Host↔Device Transfers + Kernel-Launch/Synchronisation übersteigen GPU-Gewinn, CPU erledigt weiter Domain/Boundary/Output → Wartezeiten/PCIe-Bottleneck.
- **OMP-Scaling flach:** Kleine Workpackages, Bandbreite-limitiert, viele Syncs → mehr Threads kaum schneller.

### TODO / Maßnahmen (Vorschläge)
- Profiling (perf/VTune, nvprof) auf `process_timestep`, `perform_collisions`, `perform_reactions`, Integrator-Stages.
- Domain-Suche zusammenführen/zwischenspeichern; Boundary-Check ohne zweiten vollen Lookup.
- Loop-Fusion/Chunking; größere `schedule(static, …)`-Chunks, weniger OMP-Loops.
- Integratoren: Scratch-Ensemble voralloziieren, Context reuse; Stages direkt auf SoA-Daten statt Kopien.
- Collisions/Reactions: Batch-Pfade erzwingen/ausbauen; `per_domain` Vektoren recyceln.
- Space-Charge: Update-Frequenz drosseln oder adaptiv; Direct nur für kleine N.
- I/O: HDF5-Locking vermeiden (unique tmp), `write_interval` hoch; Progress/Safety aus für Benchmarks.
- Hybrid: Transfers minimieren (pinned buffers, kein per-Step Feld-Transfer), größere Kernel-Batches; sonst GPU deaktivieren für diese Workloads.
