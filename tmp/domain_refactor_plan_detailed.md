# Domain Refactor Plan – Detailed (IGeometry/IFieldModel)

## Ziel
Modularisierung des Domain-Handling durch Strategien für Geometrie/Boundaries und Feldmodelle. DomainManager delegiert an Strategien; neue Geometrien/Felder lassen sich hinzufügen, ohne den Manager/Integrator anzupassen.

## Soll-Struktur
- `IGeometry`
  - Methoden: `contains(const Vec3&)`, `global_to_local(pos/vel)`, `local_to_global(pos/vel)`, `handle_boundary(pos, vel)` (liefert ggf. BoundaryActionResult oder Flag)
  - Daten: Domain-Index, Origin, Rotation (oder von außen übergeben)
  - Konkrete Klassen: `CylindricalGeometry`, `OrbitrapGeometry` (aktuell unterstützt), optional `BoxGeometry` (für Space Charge)
- `IFieldModel`
  - Methoden: `E(const Vec3&, double t)`, optional `B(const Vec3&, double t)`, `has_field()`
  - Konkrete Klassen: `AnalyticalUniformField`, `AnalyticalQuadrupoleField`, `OrbitrapField`, `MapField` (HDF5 CPU; GPU Textur optional/experimentell)
- `Domain`
  - Hält `std::unique_ptr<IGeometry>`, `std::unique_ptr<IFieldModel>`, Environment
  - Domain-Index/Name, evtl. BoundaryAction-Konfiguration
- `DomainFactory`
  - Baut pro `DomainConfig` die passenden Strategien (InstrumentType → Geometrie/Feld Kombi) + Environment
- `DomainManager`
  - Hält `std::vector<Domain>`
  - Delegiert `contains/transform/boundary` an `IGeometry`, Feldzugriff an `IFieldModel`
  - Domain-Wechsel/Index-Verwaltung bleibt

## Umsetzungsschritte
1) **Interfaces anlegen**
   - `IGeometry`-Interface mit obigen Methoden, Basisklasse mit Domain-Index/Origin/Rotation.
   - `IFieldModel`-Interface mit `E/B` und optionalem `has_B` Flag.

2) **Konkrete Implementierungen (Minimal)**
   - `CylindricalGeometry`: enthält Radius/Länge/Epsilon; implementiert contains, transforms, boundary (Absorption/Reflection via BoundaryAction).
   - `OrbitrapGeometry`: implementiert contains/transforms; Boundaries weiterhin auf CPU (keine GPU).
   - Felder: `AnalyticalUniformField`, `AnalyticalQuadrupoleField`, `OrbitrapField` (basierend auf vorhandenen analytischen Formeln), `MapField` (HDF5 CPU) – GPU-Map bleibt optional.

3) **Domain-Klasse einführen**
   - Neue Struct/Classe `Domain` (nicht zu verwechseln mit DomainConfig): Owning-Pointer auf Geometrie/Feld, Environment (kopiert aus Config).
   - Optional: BoundaryAction pro Domain (wie bisher in DomainManager vorbereitet).

4) **DomainFactory**
   - Neue Factory, die aus `DomainConfig` eine `Domain` baut: wählt Geometrie/Feld-Strategie anhand InstrumentType/Config (z.B. Cylindrical + QuadrupoleField, Orbitrap + OrbitrapField, etc.).
   - Fehlerfälle: Unsupported InstrumentType → Exception.

5) **DomainManager refactor**
   - Statt `domains_` als `std::vector<DomainConfig>` -> `std::vector<Domain>`.
   - `find_domain_index`, `global_to_local`, `boundary` delegieren an `Domain`/`IGeometry`.
   - Field Access in Forces ggf. über Domain/FieldModel (ForceContext aktualisieren, falls nötig).

6) **ForceContext/ForceRegistry Abgleich**
   - Sicherstellen, dass ForceContext `domain` noch verwenden kann; ggf. Zeiger auf Domain/FieldModel statt DomainConfig bereitstellen.
   - Analytische Felder in ForceRegistry → FieldModel nutzen, MapField für HDF5.

7) **Tests**
   - Unit-Tests für `CylindricalGeometry.contains/transform/boundary` und Orbitrap.
   - Tests für FieldModel (Analytical vs. bekannte Werte; MapField Sampling).
   - DomainManager-Tests: Domain Lookup/Transitions/Boundary Actions.
   - Integrations-Smoke-Test (kleiner Lauf) zur Parität mit altem Verhalten.

8) **Doku**
   - Update Docs (Architecture/README) mit neuer Domain-Strategie-Struktur und Limits (Geom: Cylinder/Orbitrap; MapField CPU; GPU optional/experimentell).

## Risiken/Abhängigkeiten
- ForceContext/AnalyticalField-Kopplung muss angepasst werden, damit Feldmodelle genutzt werden können.
- GPU-Feldpfade sind optional; CPU zuerst stabilisieren.
- BoundaryAction-Integration in neue Geometrie-Strategien sauber einbauen.

## Optional (später)
- BoxGeometry für Space-Charge-Boxen.
- GPU MapField (Texturen) als eigenes FieldModel.
