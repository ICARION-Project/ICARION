# Domain Refactor Plan (Geometrie/Feld-Strategien)

## Ziel
- Domain-Handling modularisieren: Geometrie/Boundaries und Feldmodelle als Strategien
- DomainManager delegiert an Strategien statt hardcoded Instrument-Typen
- Bessere Erweiterbarkeit für neue Geometrien/Felder, klare Separation of Concerns

## Soll-Struktur
- `IGeometry`/`IBoundary`:
  - Methoden: `contains(pos)`, `global_to_local(pos/vel)`, `local_to_global`, `handle_boundary(pos, vel)` (Absorption/Reflection/…)
  - Konkrete Klassen: `CylindricalGeometry`, `OrbitrapGeometry` (aktuell unterstützt), ggf. `BoxGeometry` (für Space Charge) später
- `IFieldModel`:
  - Methoden: `E(pos, t)`, ggf. `B(pos, t)`; optional `has_field()`
  - Konkrete Klassen: `AnalyticalUniformField`, `AnalyticalQuadrupoleField`, `OrbitrapField`, `MapField` (HDF5 CPU; GPU Textur optional)
- `Domain`:
  - Hält `std::unique_ptr<IGeometry>`, `std::unique_ptr<IFieldModel>`, `Environment`
  - Optional BoundaryAction/Reflection params
- `DomainFactory`:
  - Baut pro DomainConfig die passenden Geometrie- und Feldstrategien (InstrumentType → Geom+Field Kombi)
- `DomainManager`:
  - Hält `std::vector<Domain>`
  - Delegiert `contains/transform/boundary` an `IGeometry`
  - Feldzugriff via `IFieldModel`
  - Domain-Wechsel/Index-Management

## Schritte
1) Interfaces definieren (`IGeometry`, `IFieldModel`), minimale Methoden
2) Implementiere aktuelle Typen:
   - Geometrie: `CylindricalGeometry`, `OrbitrapGeometry`
   - Feld: `Analytical...` (Uniform/RF/Orbitrap), `MapField` (CPU HDF5); GPU Map bleibt optional
3) Domain-Klasse einführen (oder DomainConfig+Strategien kombinieren) mit Owning-Pointern
4) DomainFactory: map InstrumentType → Geometrie/Feld Strategien; Environment bleibt aus Config
5) DomainManager refactor: remove hardcoded geom checks; delegiere an `IGeometry`/`IFieldModel`; BoundaryAction weiter nutzen
6) Tests: Domain contains/transform/boundaries; Feld-Sampling (Analytisch/Map)
7) Doku: update Domain limitations (Geom/Feld), Beschreibung der Strategie-Struktur

## Notizen
- Environment bleibt Datencontainer (T/P/Dichte/Flow)
- GPU Feldmodelle weiterhin als experimentell kennzeichnen; CPU first
- Space Charge Geometrie (Box) optional ergänzen, falls sinnvoll
