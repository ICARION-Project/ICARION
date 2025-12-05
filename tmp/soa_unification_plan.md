# SoA Unification Plan (SimulationEngine)

## Ziel
- Einen einzigen Laufpfad in der SimulationEngine: intern SoA (`IonEnsemble`).
- AoS (`std::vector<IonState>`) nur als Eingabe/Ausgabe-Wrapper (oder deprecaten).
- Paritäts-/Toleranztests gegen bisherigen AoS-Pfad.

## Schritte
1) Eintrittspfad vereinheitlichen
   - In `main`: nach `config.generate_ions` AoS → `IonEnsemble` konvertieren.
   - `SimulationEngine::run` nur noch SoA nehmen (oder AoS-Wrapper ruft SoA intern).

2) Engine bereinigen
   - `process_timestep_soa` als einzige Implementierung behalten.
   - `run(std::vector<IonState>&)` entweder entfernen oder als Thin-Wrapper: AoS→SoA→AoS zurück.
   - Alle Kraft-/Kollisions-/Reaktions-Aufrufe im SoA-Pfad nutzen (kein doppelter Code).

3) GPU-Hilfen ausrichten
   - `IonStateGPU`/Upload/Download am SoA-Pfad ausrichten; nur eine Definition behalten.
   - `ForceRegistry`/CollisionHandler/Integrator bleiben unverändert, arbeiten auf SoA.

4) Tests
   - Paritätstests AoS vs. SoA für kleine Fälle (Position/Velocity innerhalb Toleranz).
   - Boundaries/Domain-Wechsel, Kollisions-/Reaktionspfade mit SoA prüfen.

5) Doku
   - README/Docs: SoA als interner Pfad dokumentieren; AoS nur für Kompatibilität (optional).

## Hinweise
- CPU bleibt Hauptpfad; GPU bleibt optional/experimentell.
- AoS→SoA-Konvertierung ggf. in einem Helper kapseln.
- Legacy AoS-Pfade deprecaten, nicht duplizieren.
