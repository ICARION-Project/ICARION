# Field Model Extraction Plan (from ElectricFieldForce)

## Ziel
- Repliziere die bestehende analytische Feldlogik aus `ElectricFieldForce` als `IFieldModel`-Strategien.
- Sicherstellen, dass Physik/Skalierung identisch bleibt.
- Tests: alter Pfad vs. neuer FieldModel-Ausgabe für jede Feldart.

## Feldtypen aus ElectricFieldForce
- Uniform/axial DC (IMS/Drift): Feld aus Spannung/Länge (params in DomainConfig.fields.dc_field_strength_Vm oder abgeleitet)
- Quadrupole RF (LQIT/QuadrupoleRF): radiales Quadrupolfeld, zeitabhängig (RF-Spannung, Frequenz, Phase, Radius) + ggf. axial DC/Radial DC
- Orbitrap: quadrologarithmisches Potential, existierende Formel in ElectricFieldForce
- AC/LQIT: AC-Feldkomponenten (Spannung, Frequenz), ggf. Sweep
- Map-Field (HDF5): existierendes `interpolate_field`

## Vorgehen
1) Inventar ElectricFieldForce
   - Extrahiere die analytischen Feldberechnungen pro InstrumentType/FieldConfig.
   - Notiere benötigte Parameter aus DomainConfig/FieldParams.

2) FieldModel-Implementierungen
   - `UniformFieldModel` (axial/gerichtet): E = (V/length) entlang z oder spezifizierter Richtung.
   - `QuadrupoleFieldModel` (RF/DC): E(r,t) entsprechend bestehender Formel (radial ∝ V/r0^2, inkl. Phase/Frequenz), optional axial DC.
   - `OrbitrapFieldModel`: nutzt vorhandene Orbitrap-Feldformel (quadrologarithmisch) aus ElectricFieldForce.
   - `ACFieldModel`/Erweiterung für LQIT: AC-Komponenten aus DomainConfig (Frequenz, Amplitude, Phase/Sweep).
   - `MapFieldModel`: Wrap um `FieldArray` + `interpolate_field` (CPU).

3) DomainFactory wiring
   - DomainConfig → wählt passende FieldModel(s) basierend auf InstrumentType/FieldConfig.
   - Kombinieren falls nötig (z.B. DC + RF): evtl. ein CompositeFieldModel, das mehrere Komponenten summiert.

4) Tests
   - Für jeden Feldtyp einen kleinen Testfall aufsetzen:
     - Uniform: DomainConfig mit Spannung/Länge → E vergleichen (alt ElectricFieldForce vs. FieldModel::E)
     - Quadrupole RF: fester t, Position (x,y) → Vergleich alt vs. neu (inkl. Phase)
     - Orbitrap: Position (x,y,z) → Vergleich alt vs. neu
     - AC/LQIT: Position, t → Vergleich alt vs. neu
     - MapField: HDF5/FieldArray Sampling identisch zu `interpolate_field`
   - Toleranzbasierter Vergleich (kein Bitgleichheit nötig).

5) Doku
   - Beschreibe Feldmodelle und Abdeckung; CPU-first, GPU-Map optional/experimentell.

## Hinweise
- Physik unverändert lassen: Formeln 1:1 aus ElectricFieldForce übernehmen.
- Composite/Summation von Feldkomponenten ggf. als eigener FieldModel-Typ.
- Tests erst aktivieren, wenn FieldModels verdrahtet sind (DomainManager/ForceRegistry).
